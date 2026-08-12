"""Command registration for CorrectRotamer."""

from chimerax.atomic import AtomicStructureArg
from chimerax.core.commands import CmdDesc
from chimerax.map import MapArg

from .correction import correct_rotamer_outliers, residue_label


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
