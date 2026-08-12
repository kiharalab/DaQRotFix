import importlib.util
import pathlib
import sys
import types
import unittest

import numpy


MODULE_PATH = pathlib.Path(__file__).parents[1] / "src" / "correction.py"
SPEC = importlib.util.spec_from_file_location("correct_rotamer_correction", MODULE_PATH)
correction = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = correction
SPEC.loader.exec_module(correction)


class FakeAtoms:
    def __init__(self, names, elements, coords):
        self.names = numpy.asarray(names)
        self.element_names = numpy.asarray(elements)
        self._coords = numpy.asarray(coords, dtype=float)

    def __len__(self):
        return len(self.names)

    def __getitem__(self, item):
        return FakeAtoms(
            self.names[item], self.element_names[item], self._coords[item]
        )

    @property
    def coords(self):
        return numpy.array(self._coords, copy=True)

    @coords.setter
    def coords(self, values):
        self._coords[...] = values


class FakeResidue:
    name = "LYS"
    chain_id = "A"
    number = 10
    insertion_code = ""

    def __init__(self):
        self.atoms = FakeAtoms(
            ["N", "CA", "C", "O", "OXT", "CB", "CG", "HB2"],
            ["N", "C", "C", "O", "O", "C", "C", "H"],
            numpy.zeros((8, 3)),
        )


class DensityVolume:
    def __init__(self):
        self.last_coords = None
        self.last_xform = None

    def interpolated_values(
        self, coords, point_xform=None, out_of_bounds_list=False
    ):
        self.last_coords = numpy.array(coords, copy=True)
        self.last_xform = point_xform
        marker = int(round(float(coords[0, 0])))
        if marker == 2:
            return numpy.zeros(len(coords)), numpy.array([0])
        return numpy.full(len(coords), marker * 2.0), numpy.array([], dtype=int)


class FakeModel:
    scene_position = object()


class FakeRotamer:
    num_chi_dihedrals = 1

    def __init__(self, residue):
        self.residue = residue
        self.targets = [
            {"Name": "outlier", "Angles": [1.0]},
            {"Name": "outside", "Angles": [2.0]},
            {"Name": "lower", "Angles": [3.0]},
            {"Name": "winner", "Angles": [4.0]},
        ]
        self.num_targets = len(self.targets)

    def get_target(self, index):
        return self.targets[index]


class CorrectionTests(unittest.TestCase):
    def test_density_uses_only_sidechain_heavy_atoms_and_model_transform(self):
        residue = FakeResidue()
        residue.atoms._coords[5, 0] = 3.0
        residue.atoms._coords[6, 0] = 3.0
        volume = DensityVolume()
        model = FakeModel()

        score = correction.mean_sidechain_density(model, residue, volume)

        self.assertEqual(score, 6.0)
        self.assertEqual(volume.last_coords.shape, (2, 3))
        self.assertIs(volume.last_xform, model.scene_position)

    def test_each_target_starts_from_original_and_only_winner_is_permanent(self):
        residue = FakeResidue()
        rotamer = FakeRotamer(residue)
        model = FakeModel()
        volume = DensityVolume()
        starts = []

        old_apply = correction.apply_rotamer_target
        old_validate = correction._rotamer_validation

        def fake_apply(_rotamer, target):
            starts.append(numpy.array(residue.atoms._coords, copy=True))
            marker = float(target["Angles"][0])
            residue.atoms._coords[5:7, 0] = marker

        def fake_validate(_manager, _residue):
            marker = int(round(float(residue.atoms._coords[5, 0])))
            return (0.0, marker == 1)

        correction.apply_rotamer_target = fake_apply
        correction._rotamer_validation = fake_validate
        try:
            result = correction.evaluate_outlier(
                model, volume, object(), rotamer
            )
        finally:
            correction.apply_rotamer_target = old_apply
            correction._rotamer_validation = old_validate

        self.assertEqual(len(starts), 5)  # four trials plus final application
        for start in starts:
            numpy.testing.assert_allclose(start, 0.0)
        self.assertEqual(result.selected.target_name, "winner")
        self.assertEqual(result.selected.density_score, 8.0)
        numpy.testing.assert_allclose(residue.atoms._coords[5:7, 0], 4.0)
        self.assertIn("rotamer outlier", result.candidates[0].rejection_reason)
        self.assertIn("outside", result.candidates[1].rejection_reason)

    def test_missing_mdff_manager_is_checked_without_creation(self):
        class UserError(Exception):
            pass

        calls = []
        extensions = types.ModuleType("chimerax.isolde.session_extensions")

        def get_mdff_mgr(model, volume, create=False):
            calls.append((model, volume, create))
            return None

        extensions.get_mdff_mgr = get_mdff_mgr
        chimerax = types.ModuleType("chimerax")
        core = types.ModuleType("chimerax.core")
        errors = types.ModuleType("chimerax.core.errors")
        errors.UserError = UserError
        isolde_module = types.ModuleType("chimerax.isolde")
        isolde_module.session_extensions = extensions

        injected = {
            "chimerax": chimerax,
            "chimerax.core": core,
            "chimerax.core.errors": errors,
            "chimerax.isolde": isolde_module,
            "chimerax.isolde.session_extensions": extensions,
        }
        old_modules = {name: sys.modules.get(name) for name in injected}
        sys.modules.update(injected)

        session = types.SimpleNamespace(
            isolde=types.SimpleNamespace(simulation_running=False)
        )
        model = object()
        volume = object()
        try:
            with self.assertRaises(UserError):
                correction.correct_rotamer_outliers(
                    session, model, volume
                )
        finally:
            for name, previous in old_modules.items():
                if previous is None:
                    sys.modules.pop(name, None)
                else:
                    sys.modules[name] = previous

        self.assertEqual(calls, [(model, volume, False)])


if __name__ == "__main__":
    unittest.main()
