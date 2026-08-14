"""Command registration for CorrectRotamer."""

from chimerax.atomic import AtomicStructureArg
from chimerax.core.commands import CmdDesc, FloatArg, IntArg
from chimerax.map import MapArg

from .correction import correct_rotamer_outliers, residue_label
from .local_mdff import (
    DEFAULT_ANNEALING_UPDATES,
    DEFAULT_FINAL_MINIMIZATION_UPDATES,
    DEFAULT_FINAL_TEMPERATURE,
    DEFAULT_INITIAL_TEMPERATURE,
    optimize_sidechains_with_local_mdff,
)


def _log_report(session, report):
    logger = session.logger
    logger.info(
        "Correct Rotamers: found "
        f"{report.original_outlier_count} original rotamer outlier(s)."
    )
    for result in report.residues:
        label = residue_label(result.residue)
        if result.selected is None:
            logger.warning(f"{label}: no valid in-map canonical rotamer; unchanged")
            continue
        logger.info(
            f"{label}: selected {result.selected.target_name}, "
            f"mean side-chain density {result.selected.density_score:.6g}"
        )
    logger.info(
        f"Correct Rotamers: corrected {report.corrected_count}/"
        f"{report.original_outlier_count} original outlier(s). "
        "No ISOLDE simulation or minimization was started."
    )


def correct_rotamers_cmd(
    session,
    structure,
    volume,
):
    report = correct_rotamer_outliers(
        session,
        structure,
        volume,
    )
    _log_report(session, report)
    return report


correct_rotamers_desc = CmdDesc(
    required=[("structure", AtomicStructureArg), ("volume", MapArg)],
    synopsis="correct ISOLDE rotamer outliers using side-chain map density",
)


def _log_local_mdff_report(session, report):
    logger = session.logger
    logger.info(
        f"Rotamer Optimizer: found {report.original_outlier_count} "
        f"original rotamer outlier(s) at threshold "
        f"{report.outlier_threshold:.6g}."
    )
    for result in report.residues:
        selected = result.selected
        label = residue_label(result.residue)
        if result.resolved_before_trial:
            logger.info(f"{label}: already resolved by an earlier local optimization")
        elif selected is None:
            logger.warning(f"{label}: no safe local improvement; unchanged")
        else:
            logger.info(
                f"{label}: selected {selected.target_name} "
                f"({selected.probability:.3g}%), "
                f"local outliers {result.baseline_local_outlier_count} -> "
                f"{selected.local_outlier_count}"
            )
    logger.info(
        f"Rotamer Optimizer: optimized {report.corrected_count}/"
        f"{report.original_outlier_count} original outlier(s)."
    )


def correct_rotamers_mdff_cmd(
    session,
    structure,
    volume,
    *,
    radius=5.0,
    initial_temperature=DEFAULT_INITIAL_TEMPERATURE,
    annealing_updates=DEFAULT_ANNEALING_UPDATES,
    final_temperature=DEFAULT_FINAL_TEMPERATURE,
    final_minimization_updates=DEFAULT_FINAL_MINIMIZATION_UPDATES,
    coordinate_updates=None,
    backbone_spring_constant=500.0,
    outlier_threshold=None,
):
    report = optimize_sidechains_with_local_mdff(
        session,
        structure,
        volume,
        local_radius=radius,
        initial_temperature=initial_temperature,
        annealing_updates=annealing_updates,
        final_temperature=final_temperature,
        final_minimization_updates=(
            coordinate_updates
            if coordinate_updates is not None
            else final_minimization_updates
        ),
        # Command value is in kJ mol-1 A-2; ISOLDE stores kJ mol-1 nm-2.
        backbone_spring_constant=backbone_spring_constant * 100.0,
        outlier_threshold=outlier_threshold,
    )
    _log_local_mdff_report(session, report)
    return report


correct_rotamers_mdff_desc = CmdDesc(
    required=[("structure", AtomicStructureArg), ("volume", MapArg)],
    keyword=[
        ("radius", FloatArg),
        ("initial_temperature", FloatArg),
        ("annealing_updates", IntArg),
        ("final_temperature", FloatArg),
        ("final_minimization_updates", IntArg),
        # Legacy v0.5.x alias for final_minimization_updates.
        ("coordinate_updates", IntArg),
        ("backbone_spring_constant", FloatArg),
        ("outlier_threshold", FloatArg),
    ],
    synopsis=(
        "optimize rotamer outliers with local ISOLDE MDFF annealing and minimization"
    ),
)
