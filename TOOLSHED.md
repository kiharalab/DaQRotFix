# DaQRotFix

**Density-assisted Quick Rotamer Fix**

DaQRotFix is a UCSF ChimeraX plugin for correcting protein rotamer outliers
with an experimental density map and ISOLDE. It provides an automated local
MDFF workflow through a graphical interface and a separate, simulation-free
command that ranks canonical rotamers by map density.

## Features

- Detects rotamer outliers using ISOLDE validation.
- Tries canonical rotamers in decreasing expected-frequency order.
- Runs short, local ISOLDE MDFF annealing and minimization trials.
- Temporarily restrains backbone N, CA, C, and O atoms while side chains remain
  mobile.
- Accepts a candidate only when it fixes the target residue, reduces the local
  outlier count, and creates no new neighboring outlier.
- Restores rejected coordinates and temporary ISOLDE settings after each
  trial.
- Provides a residue-by-residue results table with candidate details.
- Includes a density-only command that does not start an ISOLDE simulation.

## Requirements

- UCSF ChimeraX
- ISOLDE
- An atomic structure and density map already associated through ISOLDE MDFF

Before running DaQRotFix, open ISOLDE, select the structure, configure the
structure/map MDFF association, and stop any active simulation. Automated local
MDFF also requires the selected MDFF potential to be enabled.

## Graphical Interface

Open:

**Tools → Structure Editing → DaQRotFix**

Select the structure and map, adjust the local-MDFF settings under
**PARAMETERS**, and click **Start Optimization**. DaQRotFix
asks for confirmation before starting simulations and restores the original
camera view when the complete run finishes.

Default settings:

| Setting | Default |
| --- | ---: |
| Local radius | 5 Å |
| Initial temperature | 20 K |
| Annealing coordinate updates | 20 |
| Final temperature | 5 K |
| Final minimization coordinate updates | 5 |
| Backbone restraint | 500 kJ mol⁻¹ Å⁻² |
| Rotamer-outlier threshold | Current ISOLDE cutoff, normally 0.05% |

## Commands

Density-only canonical-rotamer selection:

```text
correctrotamers #1 #2
```

Automated local MDFF:

```text
correctrotamers mdff #1 #2
```

In both commands, `#1` is the atomic structure and `#2` is its
ISOLDE-associated density map.

The local-MDFF command accepts the following optional keywords:

```text
radius 5 initialTemperature 20 annealingUpdates 20 finalTemperature 5 finalMinimizationUpdates 5 backboneSpringConstant 500 outlierThreshold 0.0005
```

The command-line outlier threshold is a P-value on a 0–1 scale; the graphical
interface expresses the same threshold as a percentage.

## Density-Only Selection

The `correctrotamers` command snapshots the current outlier list and tries every
canonical target from the same original residue coordinates. Targets are
rejected if they remain rotamer outliers or place a scored side-chain heavy
atom outside the selected map. The valid target with the highest unweighted
mean density at its side-chain heavy-atom centers is applied. Rotamer
probability is not used for ranking.

This workflow does not start or minimize an ISOLDE simulation and does not
change MDFF state or position restraints.

## Automated Local MDFF Selection

For each original outlier, DaQRotFix restores the full pre-trial model
coordinates before trying a candidate. Each candidate receives local MDFF
annealing followed by final minimization. The first probability-ordered
candidate that fixes the target, reduces the local outlier count, and creates
no new neighboring outlier is retained. If no candidate qualifies, the residue
is left unchanged.

## Author

- Genki Terashi
- Purdue University
- <gterashi@purdue.edu>

## License

DaQRotFix is distributed under the GNU General Public License v3.0 only
(`GPL-3.0-only`).

## Source Code

<https://github.com/gterashi/DaQRotFix>
