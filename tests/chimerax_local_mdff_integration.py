"""Manual GUI-mode integration smoke test for one short local MDFF trial."""

from pathlib import Path

from chimerax.atomic import AtomicStructure
from chimerax.core.commands import run
from chimerax.isolde import session_extensions, tests as isolde_tests
from chimerax.isolde.cmd.cmd import isolde_start
from chimerax.map import Volume


isolde = isolde_start(session, show_splash=False)
# This smoke test must be portable across machines. Do not depend on the
# user's previously saved GPU platform setting.
isolde.sim_params.platform = "OpenCL"
# Startup scripts run before ChimeraX starts its normal frame timer, whereas
# the GUI workflow runs after it is active. Start it explicitly for this test
# because ISOLDE schedules OpenMM work from the ``new frame`` trigger.
session.update_loop.start_redraw_timer()
PDB_PATH = Path(isolde_tests.__file__).with_name("1pmx_1.pdb")
run(session, f'open "{PDB_PATH}"')
model = next(m for m in session.models.list() if isinstance(m, AtomicStructure))
isolde.change_selected_model(model)
run(session, f"molmap #{model.id_string} 3 gridSpacing 1")
volume = next(m for m in session.models.list() if isinstance(m, Volume))
mdff_mgr = session_extensions.get_mdff_mgr(model, volume, create=True)
mdff_mgr.enabled = True

rotamer_mgr = session_extensions.get_rotamer_mgr(session)
outliers, _scores = rotamer_mgr.outliers(rotamer_mgr.get_rotamers(model.residues))
from chimerax.daq_rot_fix.local_mdff import (
    IsoldeLocalMDFFRunner,
    _local_protein_residues,
)

residue = outliers[0].residue if len(outliers) else model.residues[0]
local_residues = _local_protein_residues(residue, 3.0)
original_temperature = isolde.sim_params.temperature
original_smoothing = isolde.sim_params.trajectory_smoothing
runner = IsoldeLocalMDFFRunner(
    session,
    model,
    volume,
    initial_temperature=20.0,
    annealing_updates=2,
    final_temperature=5.0,
    final_minimization_updates=1,
    max_iterations_per_update=25,
    timeout_ms=180000,
)
runner.run(local_residues)

assert not isolde.simulation_running
assert mdff_mgr.enabled
assert isolde.sim_params.temperature == original_temperature
assert isolde.sim_params.trajectory_smoothing == original_smoothing
session.logger.info(f"DAQROTFIX_LOCAL_MDFF_OK {len(outliers)}")
