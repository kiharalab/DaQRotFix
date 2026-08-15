"""Manual ChimeraX integration smoke test (run with the ChimeraX executable)."""

from types import SimpleNamespace
from pathlib import Path

from chimerax.atomic import AtomicStructure
from chimerax.core.commands import run
from chimerax.isolde import session_extensions, tests as isolde_tests
from chimerax.map import Volume


PDB_PATH = Path(isolde_tests.__file__).with_name("1pmx_1.pdb")

run(session, f'open "{PDB_PATH}"')
model = next(m for m in session.models.list() if isinstance(m, AtomicStructure))
run(session, f"molmap #{model.id_string} 3 gridSpacing 1")
volume = next(m for m in session.models.list() if isinstance(m, Volume))

# Test setup is allowed to create this manager. The bundle under test must find
# it with create=False and never creates one itself.
mdff_mgr = session_extensions.get_mdff_mgr(model, volume, create=True)
mdff_enabled_before = mdff_mgr.enabled

session.isolde = SimpleNamespace(
    simulation_running=False,
)

from chimerax.daq_rot_fix.correction import (
    correct_rotamer_outliers,
)

report = correct_rotamer_outliers(
    session, model, volume
)
assert report.original_outlier_count >= 0
assert all(
    len(result.candidates) > 0 for result in report.residues
), "Every original outlier should enumerate canonical targets"
assert mdff_mgr.enabled == mdff_enabled_before

print(
    "DAQROTFIX_INTEGRATION_OK",
    report.original_outlier_count,
    report.corrected_count,
)
