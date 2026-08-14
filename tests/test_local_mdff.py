import importlib.util
import pathlib
import sys
import types
import unittest

import numpy


ROOT = pathlib.Path(__file__).parents[1]
PACKAGE = "correctrotamer_test_package"
package = types.ModuleType(PACKAGE)
package.__path__ = [str(ROOT / "src")]
sys.modules.setdefault(PACKAGE, package)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


correction = _load(f"{PACKAGE}.correction", ROOT / "src" / "correction.py")
local_mdff = _load(f"{PACKAGE}.local_mdff", ROOT / "src" / "local_mdff.py")


class FakeAtoms:
    def __init__(self):
        self._coords = numpy.zeros((2, 3), dtype=float)

    @property
    def coords(self):
        return numpy.array(self._coords, copy=True)

    @coords.setter
    def coords(self, values):
        self._coords[...] = values


class FakeResidue:
    def __init__(self, name="LYS", number=7):
        self.name = name
        self.chain_id = "A"
        self.number = number
        self.insertion_code = ""
        self.deleted = False


class FakeRotamer:
    def __init__(self, residue):
        self.residue = residue
        self.targets = [
            {"Name": "low", "Frequency": 20.0, "Angles": [1.0]},
            {"Name": "high", "Frequency": 80.0, "Angles": [2.0]},
            {"Name": "middle", "Frequency": 50.0, "Angles": [3.0]},
        ]
        self.num_targets = len(self.targets)

    def get_target(self, index):
        return self.targets[index]


class FakeRotamerMgr:
    def __init__(self, rotamer):
        self.rotamer = rotamer
        self.cutoffs = (0.02, 0.0005)

    def get_rotamers(self, residues):
        return [self.rotamer]

    def outliers(self, rotamers):
        return [self.rotamer], numpy.array([0.0])

    def validate_rotamers(self, rotamers):
        return numpy.array([0.0 for _rotamer in rotamers])

    def get_rotamer(self, residue):
        return self.rotamer


class FakeRunner:
    def __init__(self, model):
        self.model = model
        self.starts = []

    def run(self, local_residues):
        marker = self.model.atoms._coords[0, 0]
        self.starts.append(numpy.array(self.model.atoms._coords, copy=True))
        self.model.atoms._coords[1, 0] = marker * 10.0


class LocalMDFFAlgorithmTests(unittest.TestCase):
    def test_default_two_stage_schedule(self):
        self.assertEqual(local_mdff.DEFAULT_INITIAL_TEMPERATURE, 20.0)
        self.assertEqual(local_mdff.DEFAULT_ANNEALING_UPDATES, 20)
        self.assertEqual(local_mdff.DEFAULT_FINAL_TEMPERATURE, 5.0)
        self.assertEqual(local_mdff.DEFAULT_FINAL_MINIMIZATION_UPDATES, 5)

    def test_annealing_temperature_includes_both_endpoints(self):
        temperatures = [
            local_mdff.annealing_temperature(20.0, 5.0, index, 4)
            for index in range(4)
        ]
        self.assertEqual(temperatures, [20.0, 15.0, 10.0, 5.0])

    def test_runner_rejects_heating_schedule(self):
        with self.assertRaisesRegex(ValueError, "Final temperature"):
            local_mdff.IsoldeLocalMDFFRunner(
                object(), object(), object(),
                initial_temperature=20.0,
                final_temperature=25.0,
            )

    def setUp(self):
        self.residue = FakeResidue()
        self.rotamer = FakeRotamer(self.residue)
        self.rotamer_mgr = FakeRotamerMgr(self.rotamer)
        self.model = types.SimpleNamespace(atoms=FakeAtoms(), residues=[self.residue])
        self.volume = object()
        self.session = types.SimpleNamespace(
            isolde=types.SimpleNamespace(simulation_running=False)
        )
        self.runner = FakeRunner(self.model)

        class UserError(Exception):
            pass

        extensions = types.ModuleType("chimerax.isolde.session_extensions")
        extensions.get_mdff_mgr = lambda model, volume, create=False: types.SimpleNamespace(
            enabled=True
        )
        extensions.get_rotamer_mgr = lambda session: self.rotamer_mgr
        chimerax = types.ModuleType("chimerax")
        core = types.ModuleType("chimerax.core")
        errors = types.ModuleType("chimerax.core.errors")
        errors.UserError = UserError
        isolde_module = types.ModuleType("chimerax.isolde")
        isolde_module.session_extensions = extensions
        self.injected = {
            "chimerax": chimerax,
            "chimerax.core": core,
            "chimerax.core.errors": errors,
            "chimerax.isolde": isolde_module,
            "chimerax.isolde.session_extensions": extensions,
        }
        self.old_modules = {name: sys.modules.get(name) for name in self.injected}
        sys.modules.update(self.injected)

        self.old_local = local_mdff._local_protein_residues
        self.old_apply = local_mdff.apply_rotamer_target
        self.old_density = local_mdff.mean_sidechain_density
        self.old_outliers = local_mdff._outliers_in_residues
        local_mdff._local_protein_residues = lambda residue, radius: [residue]

        def apply_target(rotamer, target):
            self.model.atoms._coords[0, 0] = target["Angles"][0]

        local_mdff.apply_rotamer_target = apply_target
        local_mdff.mean_sidechain_density = (
            lambda model, residue, volume: float(model.atoms._coords[0, 0])
        )

    def tearDown(self):
        local_mdff._local_protein_residues = self.old_local
        local_mdff.apply_rotamer_target = self.old_apply
        local_mdff.mean_sidechain_density = self.old_density
        local_mdff._outliers_in_residues = self.old_outliers
        for name, previous in self.old_modules.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous

    def _outlier_objects(self, count):
        return [types.SimpleNamespace(residue=self.residue) for _ in range(count)]

    @staticmethod
    def _outliers_for(*residues):
        return [types.SimpleNamespace(residue=residue) for residue in residues]

    def test_probability_order_and_first_zero_outlier_candidate_wins(self):
        counts = {2: 2, 3: 0, 1: 1}
        local_mdff._outliers_in_residues = (
            lambda manager, residues, threshold=None: self._outlier_objects(
                1
                if int(self.model.atoms._coords[0, 0]) == 0
                else counts[int(self.model.atoms._coords[0, 0])]
            )
        )

        report = local_mdff.optimize_sidechains_with_local_mdff(
            self.session, self.model, self.volume, runner=self.runner
        )

        result = report.residues[0]
        self.assertEqual([c.target_name for c in result.candidates], ["high", "middle"])
        self.assertEqual(result.selected.target_name, "middle")
        self.assertEqual(result.selected.local_outlier_count, 0)
        self.assertEqual(self.model.atoms._coords[0, 0], 3.0)
        self.assertEqual(self.model.atoms._coords[1, 0], 30.0)
        # Neighbor coordinates from the prior MDFF trial must be reset.
        self.assertEqual(self.runner.starts[1][1, 0], 0.0)

    def test_rejects_new_neighbor_outlier_then_accepts_safe_improvement(self):
        neighbor_one = FakeResidue("LEU", 8)
        neighbor_two = FakeResidue("ARG", 9)
        new_neighbor = FakeResidue("PHE", 10)
        calls = {"count": 0}

        def local_outliers(manager, residues, threshold=None):
            calls["count"] += 1
            if calls["count"] == 1:
                return self._outliers_for(self.residue)
            marker = int(self.model.atoms._coords[0, 0])
            if marker == 0:
                return self._outliers_for(
                    self.residue, neighbor_one, neighbor_two
                )
            if marker == 2:
                # Count improves, but a previously sound neighbor became bad.
                return self._outliers_for(neighbor_one, new_neighbor)
            if marker == 3:
                return self._outliers_for(neighbor_one, neighbor_two)
            return self._outliers_for(self.residue)

        local_mdff._outliers_in_residues = local_outliers

        report = local_mdff.optimize_sidechains_with_local_mdff(
            self.session, self.model, self.volume, runner=self.runner
        )

        result = report.residues[0]
        self.assertEqual([c.target_name for c in result.candidates], ["high", "middle"])
        self.assertEqual(result.selected.target_name, "middle")
        self.assertEqual(result.baseline_local_outlier_count, 3)
        self.assertEqual(result.selected.local_outlier_count, 2)
        self.assertEqual(result.candidates[0].new_outlier_labels, ("A:PHE10",))
        self.assertEqual(self.model.atoms._coords[0, 0], 3.0)
        self.assertEqual(self.model.atoms._coords[1, 0], 30.0)

    def test_no_safe_improvement_leaves_baseline_coordinates_unchanged(self):
        neighbor = FakeResidue("LEU", 8)
        new_neighbor = FakeResidue("PHE", 10)
        calls = {"count": 0}

        def local_outliers(manager, residues, threshold=None):
            calls["count"] += 1
            if calls["count"] == 1:
                return self._outliers_for(self.residue)
            marker = int(self.model.atoms._coords[0, 0])
            if marker == 0:
                return self._outliers_for(self.residue, neighbor)
            if marker == 2:
                return self._outliers_for(self.residue)
            if marker == 3:
                return self._outliers_for(neighbor, new_neighbor)
            return self._outliers_for(self.residue, neighbor)

        local_mdff._outliers_in_residues = local_outliers
        report = local_mdff.optimize_sidechains_with_local_mdff(
            self.session, self.model, self.volume, runner=self.runner
        )

        result = report.residues[0]
        self.assertIsNone(result.selected)
        self.assertEqual(len(result.candidates), 3)
        self.assertEqual(self.model.atoms._coords[0, 0], 0.0)
        self.assertEqual(self.model.atoms._coords[1, 0], 0.0)

    def test_custom_threshold_is_used_for_initial_and_local_validation(self):
        thresholds = []
        calls = {"count": 0}

        def threshold_outliers(manager, residues, threshold=None):
            thresholds.append(threshold)
            calls["count"] += 1
            if calls["count"] <= 2:
                return self._outlier_objects(1)
            return []

        local_mdff._outliers_in_residues = threshold_outliers
        report = local_mdff.optimize_sidechains_with_local_mdff(
            self.session,
            self.model,
            self.volume,
            outlier_threshold=0.01,
            runner=self.runner,
        )

        self.assertEqual(report.outlier_threshold, 0.01)
        self.assertEqual(thresholds, [0.01, 0.01, 0.01])

    def test_residue_resolved_by_previous_trial_is_skipped(self):
        calls = {"count": 0}

        def local_outliers(manager, residues, threshold=None):
            calls["count"] += 1
            return self._outlier_objects(1) if calls["count"] == 1 else []

        local_mdff._outliers_in_residues = local_outliers
        report = local_mdff.optimize_sidechains_with_local_mdff(
            self.session, self.model, self.volume, runner=self.runner
        )

        result = report.residues[0]
        self.assertTrue(result.resolved_before_trial)
        self.assertEqual(result.candidates, [])
        self.assertEqual(self.runner.starts, [])


if __name__ == "__main__":
    unittest.main()
