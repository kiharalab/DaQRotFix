"""Automated probability-ordered rotamer trials with local ISOLDE MDFF."""

from dataclasses import dataclass, field
from typing import List, Optional

import numpy

from .correction import (
    BACKBONE_ATOM_NAMES,
    apply_rotamer_target,
    mean_sidechain_density,
    residue_label,
)


DEFAULT_LOCAL_RADIUS = 5.0
DEFAULT_COORDINATE_UPDATES = 3
DEFAULT_MAX_ITERATIONS_PER_UPDATE = 100
DEFAULT_OUTLIER_THRESHOLD = 0.0005
# ISOLDE/OpenMM stores position-restraint k in kJ mol-1 nm-2. This is
# equivalent to 500 kJ mol-1 A-2 and keeps the backbone nearly fixed.
DEFAULT_BACKBONE_SPRING_CONSTANT = 50000.0


@dataclass(frozen=True)
class LocalMDFFCandidate:
    target_index: int
    target_name: str
    probability: float
    local_outlier_count: Optional[int] = None
    target_is_outlier: Optional[bool] = None
    density_score: Optional[float] = None
    outlier_labels: tuple = ()
    rejection_reason: Optional[str] = None

    @property
    def completed(self):
        return self.rejection_reason is None and self.local_outlier_count is not None


@dataclass
class LocalMDFFResidueResult:
    residue: object
    candidates: List[LocalMDFFCandidate] = field(default_factory=list)
    selected: Optional[LocalMDFFCandidate] = None

    @property
    def corrected(self):
        return self.selected is not None


@dataclass
class LocalMDFFReport:
    original_outlier_count: int
    residues: List[LocalMDFFResidueResult]
    outlier_threshold: float

    @property
    def corrected_count(self):
        return sum(result.corrected for result in self.residues)


def target_probability(target):
    """Return an ISOLDE target's expected frequency as a numeric value."""
    for key in ("Frequency", "Probability", "frequency", "probability"):
        if key in target:
            return float(target[key])
    return 0.0


def probability_ordered_targets(rotamer):
    """Return ``(original index, target)`` pairs in decreasing probability."""
    targets = [(index, rotamer.get_target(index)) for index in range(rotamer.num_targets)]
    return sorted(
        targets,
        key=lambda item: (-target_probability(item[1]), item[0]),
    )


def _local_protein_residues(residue, radius):
    from chimerax.atomic import Residues
    from chimerax.isolde.selections import get_shell_of_residues

    center = Residues([residue])
    shell = get_shell_of_residues(center, float(radius))
    candidates = list(center) + list(shell)
    # A local OpenMM selection containing waters/ligands is more likely to
    # fail parameterisation and does not help rotamer validation. Keep amino
    # acid residues only, while ISOLDE supplies its own fixed context shell.
    protein_type = getattr(residue, "PT_AMINO", None)
    if protein_type is not None:
        candidates = [
            r for r in candidates if getattr(r, "polymer_type", None) == protein_type
        ]
    unique = []
    seen = set()
    for candidate in candidates:
        marker = id(candidate)
        if marker not in seen:
            seen.add(marker)
            unique.append(candidate)
    return Residues(unique)


def _outliers_in_residues(rotamer_mgr, residues, outlier_threshold=None):
    rotamers = rotamer_mgr.get_rotamers(residues)
    if outlier_threshold is None:
        outlier_threshold = float(rotamer_mgr.cutoffs[1])
    scores = rotamer_mgr.validate_rotamers(rotamers)
    return [
        rotamer
        for rotamer, score in zip(rotamers, scores)
        if 0.0 <= float(score) < float(outlier_threshold)
    ]


class IsoldeLocalMDFFRunner:
    """Run one short, explicitly requested local ISOLDE minimisation."""

    def __init__(
        self,
        session,
        model,
        volume,
        *,
        coordinate_updates=DEFAULT_COORDINATE_UPDATES,
        max_iterations_per_update=DEFAULT_MAX_ITERATIONS_PER_UPDATE,
        backbone_spring_constant=DEFAULT_BACKBONE_SPRING_CONSTANT,
        timeout_ms=180000,
    ):
        self.session = session
        self.model = model
        self.volume = volume
        self.coordinate_updates = max(1, int(coordinate_updates))
        self.max_iterations_per_update = max(1, int(max_iterations_per_update))
        self.backbone_spring_constant = float(backbone_spring_constant)
        self.timeout_ms = int(timeout_ms)

    @staticmethod
    def _backbone_atoms(residues):
        atoms = residues.atoms
        return atoms[numpy.isin(atoms.names, BACKBONE_ATOM_NAMES)]

    @staticmethod
    def _save_restraints(manager, atoms):
        saved = {}
        for atom in atoms:
            restraint = manager.get_restraint(atom)
            if restraint is not None:
                saved[id(atom)] = (
                    bool(restraint.enabled),
                    float(restraint.spring_constant),
                    numpy.array(restraint.target, copy=True),
                )
        return saved

    @staticmethod
    def _restore_restraints(manager, atoms, saved):
        for atom in atoms:
            restraint = manager.get_restraint(atom)
            if restraint is None:
                continue
            prior = saved.get(id(atom))
            if prior is None:
                restraint.enabled = False
            else:
                restraint.enabled = prior[0]
                restraint.spring_constant = prior[1]
                restraint.target = prior[2]

    def run(self, local_residues):
        from Qt.QtCore import QEventLoop, QTimer
        from chimerax.core.errors import UserError
        from chimerax.core.triggerset import DEREGISTER
        from chimerax.isolde import session_extensions

        isolde = getattr(self.session, "isolde", None)
        if isolde is None:
            raise UserError(
                "Open ISOLDE and select the model before starting automated local MDFF."
            )
        if isolde.simulation_running:
            raise UserError("Stop the current ISOLDE simulation before automation.")

        mdff_mgr = session_extensions.get_mdff_mgr(
            self.model, self.volume, create=False
        )
        if mdff_mgr is None or not mdff_mgr.enabled:
            raise UserError(
                "The selected map must already have an enabled ISOLDE MDFF potential."
            )

        previous_model = isolde.selected_model
        from chimerax.atomic import AtomicStructure

        selection_snapshots = []
        for structure in self.session.models.list(type=AtomicStructure):
            selection_snapshots.append(
                (
                    structure,
                    numpy.array(structure.atoms.selected, copy=True),
                    numpy.array(structure.bonds.selected, copy=True),
                )
            )
        old_padding = isolde.params.num_selection_padding_residues
        old_soft_shell = isolde.params.soft_shell_cutoff_distance
        old_max_iterations = isolde.sim_params.minimization_max_iterations
        backbone = self._backbone_atoms(local_residues)
        restraint_mgr = session_extensions.get_position_restraint_mgr(self.model)
        saved_restraints = self._save_restraints(restraint_mgr, backbone)
        stopped = False
        timed_out = False
        handler = None

        def clean_failed_start():
            """Undo ISOLDE's partially constructed state after start failure."""
            sim_manager = isolde.sim_manager
            if sim_manager is None:
                return
            try:
                sim_manager._sim_end_cb(None, "startup failure")
            finally:
                isolde._sim_end_cb(None, "startup failure")

        try:
            if isolde.selected_model is not self.model:
                isolde.change_selected_model(self.model)

            self.session.selection.clear()
            local_residues.atoms.selected = True
            # Show the exact mobile region while its candidate is minimized.
            # The camera intentionally remains there so progress is visible.
            from chimerax.core.commands import run

            run(self.session, "view sel", log=False)
            # The local residues have already been chosen explicitly. Prevent
            # ISOLDE from expanding the mobile set a second time; its hard
            # shell remains available as fixed physical context.
            isolde.params.num_selection_padding_residues = 0
            isolde.params.soft_shell_cutoff_distance = 0.0
            # Bound every OpenMM minimisation batch. ISOLDE's default of 1000
            # iterations is too long for an automated candidate trial.
            isolde.sim_params.minimization_max_iterations = (
                self.max_iterations_per_update
            )

            restraints = restraint_mgr.add_restraints(backbone)
            restraints.targets = restraints.atoms.coords
            restraints.spring_constants = self.backbone_spring_constant
            restraints.enableds = True

            loop = QEventLoop()
            timer = QTimer()
            timer.setSingleShot(True)
            updates = {"count": 0}

            def simulation_terminated(*_):
                nonlocal stopped
                stopped = True
                loop.quit()
                return DEREGISTER

            def coordinate_updated(*_):
                updates["count"] += 1
                if updates["count"] >= self.coordinate_updates:
                    isolde.commit_sim()
                    return DEREGISTER

            def timeout():
                nonlocal timed_out
                timed_out = True
                if isolde.simulation_running:
                    isolde.discard_sim(revert_to="start", warn=False)
                else:
                    loop.quit()

            try:
                isolde.start_sim()
            except Exception:
                clean_failed_start()
                raise
            if not isolde.simulation_running:
                raise RuntimeError("ISOLDE did not start the requested local simulation")
            sim_handler = isolde.sim_handler
            # Keep the short run in minimisation mode. Setting the documented
            # SimHandler flag avoids Isolde.minimize() in ISOLDE 1.11.1, whose
            # implementation calls a GUI callback that is not present on the
            # Isolde object. MDFF remains an ISOLDE-configured active force.
            sim_handler.minimize = True
            handler = sim_handler.triggers.add_handler(
                "coord update", coordinate_updated
            )
            isolde.triggers.add_handler(
                isolde.SIMULATION_TERMINATED, simulation_terminated
            )
            timer.timeout.connect(timeout)
            timer.start(self.timeout_ms)
            loop.exec()
            timer.stop()

            if timed_out:
                raise RuntimeError("The local ISOLDE MDFF run timed out")
            if not stopped:
                raise RuntimeError("The local ISOLDE MDFF run ended unexpectedly")
        finally:
            # A failed start can leave a partial simulation alive.
            if isolde.simulation_running:
                isolde.discard_sim(revert_to="start", warn=False)
            if handler is not None:
                try:
                    handler.remove()
                except Exception:
                    pass
            self._restore_restraints(restraint_mgr, backbone, saved_restraints)
            isolde.params.num_selection_padding_residues = old_padding
            isolde.params.soft_shell_cutoff_distance = old_soft_shell
            isolde.sim_params.minimization_max_iterations = old_max_iterations
            if isolde.selected_model is not previous_model:
                isolde.change_selected_model(previous_model)
            for structure, atom_selection, bond_selection in selection_snapshots:
                if not structure.deleted:
                    structure.atoms.selected = atom_selection
                    structure.bonds.selected = bond_selection


def optimize_sidechains_with_local_mdff(
    session,
    model,
    volume,
    *,
    local_radius=DEFAULT_LOCAL_RADIUS,
    coordinate_updates=DEFAULT_COORDINATE_UPDATES,
    backbone_spring_constant=DEFAULT_BACKBONE_SPRING_CONSTANT,
    outlier_threshold=None,
    progress_callback=None,
    runner=None,
):
    """Try canonical rotamers by probability and assess them after local MDFF.

    The original outlier list is snapshotted once. For each residue, every
    candidate starts from the coordinates accepted for all preceding residues.
    The first candidate producing zero local rotamer outliers is accepted. If
    none reaches zero, the completed candidate with the fewest local outliers
    is retained; probability order is the deterministic tie-breaker.
    """
    from chimerax.core.errors import UserError
    from chimerax.isolde import session_extensions

    isolde = getattr(session, "isolde", None)
    if isolde is None:
        raise UserError("Open ISOLDE before starting automated local MDFF.")
    if isolde.simulation_running:
        raise UserError("Stop the current ISOLDE simulation before automation.")
    mdff_mgr = session_extensions.get_mdff_mgr(model, volume, create=False)
    if mdff_mgr is None or not mdff_mgr.enabled:
        raise UserError(
            "The selected map must already have an enabled ISOLDE MDFF manager."
        )

    rotamer_mgr = session_extensions.get_rotamer_mgr(session)
    if outlier_threshold is None:
        outlier_threshold = float(rotamer_mgr.cutoffs[1])
    outlier_threshold = float(outlier_threshold)
    if not 0.0 < outlier_threshold <= 1.0:
        raise UserError(
            "The rotamer outlier threshold must be greater than 0 and at most 1."
        )
    original_outliers = _outliers_in_residues(
        rotamer_mgr, model.residues, outlier_threshold
    )
    original_residues = [rotamer.residue for rotamer in original_outliers]
    total = len(original_residues)
    if runner is None:
        runner = IsoldeLocalMDFFRunner(
            session,
            model,
            volume,
            coordinate_updates=coordinate_updates,
            backbone_spring_constant=backbone_spring_constant,
        )
    if progress_callback is not None:
        progress_callback("starting", 0, total, None, None, None)

    results = []
    for residue_index, residue in enumerate(original_residues, start=1):
        result = LocalMDFFResidueResult(residue)
        results.append(result)
        if getattr(residue, "deleted", False):
            continue
        rotamer = rotamer_mgr.get_rotamer(residue)
        if rotamer is None:
            continue

        local_residues = _local_protein_residues(residue, local_radius)
        baseline_coords = numpy.array(model.atoms.coords, copy=True)
        best_coords = None
        best_candidate = None

        try:
            targets = probability_ordered_targets(rotamer)
            for trial_index, (target_index, target) in enumerate(targets, start=1):
                model.atoms.coords = baseline_coords
                name = str(target.get("Name", f"target-{target_index}"))
                probability = target_probability(target)
                if progress_callback is not None:
                    progress_callback(
                        "candidate",
                        residue_index - 1,
                        total,
                        residue,
                        trial_index,
                        len(targets),
                    )
                try:
                    apply_rotamer_target(rotamer, target)
                    density = mean_sidechain_density(model, residue, volume)
                    runner.run(local_residues)
                    local_outliers = _outliers_in_residues(
                        rotamer_mgr, local_residues, outlier_threshold
                    )
                    labels = tuple(
                        residue_label(outlier.residue) for outlier in local_outliers
                    )
                    target_is_outlier = any(
                        outlier.residue is residue for outlier in local_outliers
                    )
                    candidate = LocalMDFFCandidate(
                        target_index,
                        name,
                        probability,
                        local_outlier_count=len(local_outliers),
                        target_is_outlier=target_is_outlier,
                        density_score=density,
                        outlier_labels=labels,
                    )
                    result.candidates.append(candidate)

                    if (
                        best_candidate is None
                        or candidate.local_outlier_count
                        < best_candidate.local_outlier_count
                    ):
                        best_candidate = candidate
                        best_coords = numpy.array(model.atoms.coords, copy=True)
                    if candidate.local_outlier_count == 0:
                        break
                except Exception as error:
                    result.candidates.append(
                        LocalMDFFCandidate(
                            target_index,
                            name,
                            probability,
                            rejection_reason=str(error),
                        )
                    )
        finally:
            model.atoms.coords = baseline_coords

        if best_candidate is not None and best_coords is not None:
            model.atoms.coords = best_coords
            result.selected = best_candidate
        if progress_callback is not None:
            progress_callback(
                "residue", residue_index, total, residue, result, None
            )

    report = LocalMDFFReport(total, results, outlier_threshold)
    if progress_callback is not None:
        progress_callback("complete", total, total, None, report, None)
    return report
