"""Manual ChimeraX smoke test for camera-only view restoration."""

from pathlib import Path

import numpy

from chimerax.core.commands import run
from chimerax.isolde import tests as isolde_tests

from chimerax.daq_rot_fix.tool import _capture_view, _restore_view


pdb_path = Path(isolde_tests.__file__).with_name("1pmx_1.pdb")
run(session, f'open "{pdb_path}"')
run(session, "view")

state = _capture_view(session)
position_before = session.main_view.camera.position.matrix
run(session, "turn y 37")
run(session, "zoom 2")
assert not numpy.allclose(session.main_view.camera.position.matrix, position_before)

_restore_view(session, state)
assert numpy.allclose(session.main_view.camera.position.matrix, position_before)
session.logger.info("DAQROTFIX_VIEW_RESTORE_OK")
