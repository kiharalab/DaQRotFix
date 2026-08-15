"""Density-only canonical rotamer selection using ISOLDE validation data."""

from dataclasses import dataclass, field
from math import degrees, isfinite
from typing import List, Optional

import numpy


BACKBONE_ATOM_NAMES = ("N", "CA", "C", "O")
SIDECHAIN_SCORE_EXCLUDED_NAMES = BACKBONE_ATOM_NAMES + ("OXT",)


@dataclass(frozen=True)
class CandidateEvaluation:
    target_index: int
    target_name: str
    density_score: Optional[float] = None
    rejection_reason: Optional[str] = None

    @property
    def valid(self):
        return self.rejection_reason is None and self.density_score is not None


@dataclass
class ResidueCorrection:
    residue: object
    candidates: List[CandidateEvaluation] = field(default_factory=list)
    selected: Optional[CandidateEvaluation] = None

    @property
    def corrected(self):
        return self.selected is not None


@dataclass
class CorrectionReport:
    original_outlier_count: int
    residues: List[ResidueCorrection]

    @property
    def corrected_count(self):
        return sum(result.corrected for result in self.residues)


def residue_label(residue):
    """Return a stable, human-readable residue identifier."""
    chain = getattr(residue, "chain_id", "") or "?"
    number = getattr(residue, "number", "?")
    insertion = getattr(residue, "insertion_code", "") or ""
    name = getattr(residue, "name", "UNK")
    return f"{chain}:{name}{number}{insertion}"


def sidechain_heavy_atoms(residue):
    """Return non-hydrogen side-chain atoms, excluding the protein backbone."""
    atoms = residue.atoms
    mask = numpy.logical_and(
        atoms.element_names != "H",
        numpy.logical_not(
            numpy.isin(atoms.names, SIDECHAIN_SCORE_EXCLUDED_NAMES)
        ),
    )
    return atoms[mask]


def mean_sidechain_density(model, residue, volume):
    """Score a residue with the unweighted mean map value at side-chain atoms.

    Atom coordinates are local to ``model``. Passing its scene transform lets
    ``Volume.interpolated_values`` correctly handle independently positioned
    model and map instances.
    """
    atoms = sidechain_heavy_atoms(residue)
    if len(atoms) == 0:
        raise ValueError("residue has no side-chain heavy atoms")

    values, outside = volume.interpolated_values(
        atoms.coords,
        point_xform=model.scene_position,
        out_of_bounds_list=True,
    )
    if len(outside):
        raise ValueError(
            f"{len(outside)} side-chain atom(s) are outside the selected map"
        )
    score = float(numpy.mean(values))
    if not isfinite(score):
        raise ValueError("map interpolation produced a non-finite density score")
    return score


def apply_rotamer_target(rotamer, target):
    """Set chi angles directly without minimization or backbone movement.

    This mirrors the coordinate construction used by ISOLDE's rotamer preview:
    each ISOLDE-provided moving-atom set is rotated around its chi bond. The
    caller is responsible for restoring original coordinates first.
    """
    target_angles = numpy.asarray(target["Angles"], dtype=float)
    if len(target_angles) != rotamer.num_chi_dihedrals:
        raise ValueError("rotamer target has the wrong number of chi angles")

    from chimerax.geometry import matrix, rotation

    for chi_index, target_angle in enumerate(target_angles):
        chi = rotamer.chi_dihedrals[chi_index]
        moving_atoms = rotamer.moving_atoms(chi_index)
        chi_coords = chi.atoms.coords
        axis = chi_coords[2] - chi_coords[1]
        center = matrix.project_to_axis(chi_coords[3], axis, chi_coords[1])
        transform = rotation(
            axis,
            degrees(float(target_angle) - float(chi.angle)),
            center,
        )
        moving_atoms.coords = transform.transform_points(moving_atoms.coords)


def _rotamer_validation(rotamer_mgr, residue):
    from chimerax.atomic import Residues

    score = float(rotamer_mgr.validate_residues(Residues([residue]))[0])
    outlier_cutoff = float(rotamer_mgr.cutoffs[1])
    return score, score < outlier_cutoff


def _target_name(target, target_index):
    return str(target.get("Name", f"target-{target_index}"))


def evaluate_outlier(model, volume, rotamer_mgr, rotamer):
    """Evaluate every canonical target and permanently apply only the best."""
    residue = rotamer.residue
    residue_atoms = residue.atoms
    original_coords = numpy.array(residue_atoms.coords, copy=True)
    evaluations = []
    valid = []

    targets = [rotamer.get_target(i) for i in range(rotamer.num_targets)]
    try:
        for target_index, target in enumerate(targets):
            # This restoration is deliberately inside the loop. No candidate
            # is allowed to inherit coordinates from the preceding candidate.
            residue_atoms.coords = original_coords
            name = _target_name(target, target_index)
            try:
                apply_rotamer_target(rotamer, target)
                validation_score, is_outlier = _rotamer_validation(
                    rotamer_mgr, residue
                )
                if is_outlier:
                    evaluations.append(
                        CandidateEvaluation(
                            target_index,
                            name,
                            rejection_reason=(
                                "remains a rotamer outlier "
                                f"(validation score {validation_score:.6g})"
                            ),
                        )
                    )
                    continue

                density_score = mean_sidechain_density(model, residue, volume)
                evaluation = CandidateEvaluation(
                    target_index, name, density_score=density_score
                )
                evaluations.append(evaluation)
                valid.append((evaluation, target))
            except (ArithmeticError, IndexError, KeyError, TypeError, ValueError) as error:
                evaluations.append(
                    CandidateEvaluation(
                        target_index, name, rejection_reason=str(error)
                    )
                )
    finally:
        residue_atoms.coords = original_coords

    result = ResidueCorrection(residue=residue, candidates=evaluations)
    if not valid:
        return result

    # ISOLDE stores targets in decreasing probability order. Probability must
    # not break density ties, so target name (and then its angles) provides a
    # deterministic, probability-independent tie break.
    selected_eval, selected_target = max(
        valid,
        key=lambda item: (
            item[0].density_score,
            item[0].target_name,
            tuple(float(angle) for angle in item[1]["Angles"]),
        ),
    )

    residue_atoms.coords = original_coords
    try:
        apply_rotamer_target(rotamer, selected_target)
    except Exception:
        residue_atoms.coords = original_coords
        raise
    result.selected = selected_eval
    return result


def correct_rotamer_outliers(
    session,
    model,
    volume,
    *,
    progress_callback=None,
):
    """Correct original rotamer outliers without starting ISOLDE simulation."""
    from chimerax.core.errors import UserError
    from chimerax.isolde import session_extensions

    isolde = getattr(session, "isolde", None)
    if isolde is not None and isolde.simulation_running:
        raise UserError("Stop the current ISOLDE simulation before correcting rotamers.")

    # Required preflight. In particular, create=False is intentional and must
    # not be changed: scoring alone must never create an MDFF manager.
    mdff_mgr = session_extensions.get_mdff_mgr(model, volume, create=False)
    if mdff_mgr is None:
        raise UserError(
            "The selected map has no ISOLDE MDFF manager for this structure. "
            "Add/enable the map in ISOLDE first; DaQRotFix will not "
            "create a manager merely for scoring."
        )

    rotamer_mgr = session_extensions.get_rotamer_mgr(session)
    all_rotamers = rotamer_mgr.get_rotamers(model.residues)
    original_outliers, _scores = rotamer_mgr.outliers(all_rotamers)
    # Snapshot the original objects before making any coordinate changes.
    original_outliers = list(original_outliers)
    total_outliers = len(original_outliers)
    if progress_callback is not None:
        progress_callback("enumerating", 0, total_outliers, None, None)

    results = []
    for index, rotamer in enumerate(original_outliers, start=1):
        result = evaluate_outlier(model, volume, rotamer_mgr, rotamer)
        results.append(result)
        if progress_callback is not None:
            progress_callback(
                "enumerating", index, total_outliers, rotamer.residue, result
            )

    report = CorrectionReport(total_outliers, results)
    if progress_callback is not None:
        progress_callback("complete", total_outliers, total_outliers, None, None)
    return report
