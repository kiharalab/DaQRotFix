<p align="center">
  <img src="img/daqrotfix-logo.png" alt="DaQRotFix logo" width="280">
</p>

# DaQRotFix

**Density-assisted Quick Rotamer Fix**

DaQRotFix is a plugin for [UCSF ChimeraX](https://www.rbvi.ucsf.edu/chimerax/)
that corrects protein rotamer outliers using an experimental density map and
[ISOLDE](https://tristanic.github.io/isolde/). It provides a graphical local
MDFF optimizer and a separate, simulation-free density-ranking command.

> DaQRotFix is currently alpha software. Save your model before running an
> automated optimization.

## Requirements

- UCSF ChimeraX
- ISOLDE installed in ChimeraX
- An atomic structure and a density map already associated through ISOLDE MDFF

DaQRotFix never creates an MDFF manager. Configure the structure/map pair in
ISOLDE before running either workflow. Stop any active ISOLDE simulation first.

## Installation

Download or clone this repository, then enter the following command in the
**ChimeraX command line** (not a system terminal):

```text
devel install /absolute/path/to/DaQRotFix
```

Restart ChimeraX if the tool does not appear immediately after its first
installation.

## Graphical workflow

1. Open the structure and density map in ChimeraX.
2. Open ISOLDE, select the structure, and configure and enable the map's MDFF
   potential.
3. Open **Tools → Validation → DaQRotFix**.
4. Select the structure and map, adjust options in **PARAMETERS**, and click
   **Optimize Side Chains with Local MDFF**.

The tool asks for confirmation before starting any simulations. During each
candidate trial it focuses the local region; after the run it restores the
original camera view and displays a residue-by-residue results table.

Default GUI settings are:

| Setting | Default |
| --- | ---: |
| Local radius | 5 Å |
| Initial temperature | 20 K |
| Annealing updates | 20 |
| Final temperature | 5 K |
| Final minimization updates | 5 |
| Backbone restraint | 500 kJ mol⁻¹ Å⁻² |
| Rotamer-outlier threshold | ISOLDE cutoff, normally 0.05% |

The update controls accept values from 1 to 100. The GUI expresses the outlier
threshold as a percentage.

## Commands

### Density-only correction

```text
correctrotamers #1 #2
```

Here `#1` is an atomic structure and `#2` is its ISOLDE-associated density map.
This command snapshots the current ISOLDE rotamer outliers and tries every
canonical target for each residue. It rejects targets that remain outliers or
place a side-chain heavy atom outside the map, then applies the valid target
with the highest unweighted mean density at its side-chain heavy-atom centers.
Rotamer probability is not used for ranking.

This workflow does not start or minimize a simulation, change the MDFF state,
or create or modify position restraints. Its MDFF manager must exist, but does
not need to be enabled and the ISOLDE tool does not need to be open.

### Automated local MDFF

```text
correctrotamers mdff #1 #2
```

ISOLDE must be open, and the selected structure/map pair must have an existing,
enabled MDFF manager. The command accepts these optional keywords:

```text
correctrotamers mdff #1 #2 radius 5 initialTemperature 20 annealingUpdates 20 finalTemperature 5 finalMinimizationUpdates 5 backboneSpringConstant 500 outlierThreshold 0.0005
```

| Keyword | Meaning | Default |
| --- | --- | ---: |
| `radius` | Local protein-region radius in Å | `5` |
| `initialTemperature` | Starting annealing temperature in K | `20` |
| `annealingUpdates` | Dynamics coordinate updates during cooling | `20` |
| `finalTemperature` | Ending annealing temperature in K | `5` |
| `finalMinimizationUpdates` | Minimization coordinate updates after cooling | `5` |
| `backboneSpringConstant` | N/CA/C/O restraint in kJ mol⁻¹ Å⁻² | `500` |
| `outlierThreshold` | Rotamer P-value on a 0–1 scale | current ISOLDE cutoff |

`coordinateUpdates` remains available as a legacy alias for
`finalMinimizationUpdates`.

## How automated local MDFF selects a rotamer

DaQRotFix snapshots the original outlier list, then processes those residues in
turn. For each residue it:

1. Builds a local protein region and records its current rotamer outliers.
2. Tries canonical targets in decreasing ISOLDE expected-frequency order.
3. Restores the full pre-trial model coordinates before each target.
4. Applies the target, checks that its scored side-chain atoms lie in the map,
   and runs local MDFF annealing followed by final minimization.
5. Accepts the first target that fixes the target residue, reduces the local
   outlier count, and creates no new neighboring outlier.

If no target meets all three acceptance conditions, the residue remains at its
pre-trial coordinates. A residue fixed indirectly by an earlier local trial is
recognized and skipped.

During a trial, only the local protein region is mobile and ISOLDE supplies its
fixed physical context shell. Backbone N, CA, C, and O atoms receive temporary
position restraints; side chains do not. Rejected coordinates cannot leak into
later trials. DaQRotFix also restores the user's prior backbone-restraint
settings, selection, selected ISOLDE model, temperature, trajectory smoothing,
and relevant simulation parameters after each trial.

## Development

Run the unit tests from a Python environment with NumPy installed:

```text
python3 -m unittest discover -s tests -p 'test_*.py'
```

Build the bundle from the ChimeraX command line:

```text
devel build /absolute/path/to/DaQRotFix
```

The `tests/chimerax_*_integration.py` files are ChimeraX smoke tests and must be
run with ChimeraX rather than a standard Python interpreter.

## License

[MIT](LICENSE)
