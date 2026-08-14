# CorrectRotamer for ChimeraX

CorrectRotamer is a ChimeraX bundle for automated local-MDFF optimization of
ISOLDE rotamer outliers.

It provides two deliberately separate workflows.

### Density-only correction

For each outlier present at the start of a run, the bundle restores the
original residue coordinates before trying every canonical ISOLDE target. It
discards targets that remain validation outliers or place any side-chain heavy
atom outside the map, then keeps the valid target with the highest unweighted
mean interpolated density. Rotamer probability is never used for ranking.

After enumeration, this workflow stops. It does not start, pause, or minimize
an ISOLDE simulation and does not create or modify position restraints.

### Automated local MDFF

This explicit workflow snapshots the original outlier list, then handles each
residue in turn. Canonical rotamers are tried in decreasing ISOLDE `Frequency`
order. Each candidate starts from the coordinates accepted before that
residue, receives local ISOLDE MDFF annealing followed by final minimization,
and is scored by the number of rotamer outliers in the local protein region.

A candidate is accepted only if the target residue is no longer an outlier,
the total local outlier count decreases, and no previously sound neighboring
residue becomes a new outlier. The first candidate meeting all three conditions
is retained in probability order. If none qualifies, the residue remains at
its pre-trial coordinates. Rejected local-MDFF coordinates are restored.
Protein backbone N, CA, C and O atoms are strongly and temporarily
position-restrained while side-chain atoms remain mobile.

## Install

ISOLDE must already be installed. From the ChimeraX command line:

```text
devel install /path/to/CorrectRotamer
```

Restart ChimeraX after the first installation if the tool is not immediately
listed. Open **Tools → Validation → Rotamer Optimizer**, or run:

```text
correctrotamers #1 #2
```

The local-MDFF workflow is also available as:

```text
correctrotamers mdff #1 #2 radius 5 initialTemperature 20 annealingUpdates 20 finalTemperature 5 finalMinimizationUpdates 20 backboneSpringConstant 500 outlierThreshold 0.0005
```

The graphical tool uses **Optimize Side Chains with Local MDFF** as its standard
workflow. Its separate **PARAMETERS** tab contains controls for local radius,
initial temperature, annealing coordinate updates, final temperature, final
minimization coordinate updates, backbone restraint strength and the
rotamer-outlier threshold. The default cooling schedule is 20 K to 5 K over
20 coordinate updates, followed by 20 minimization coordinate updates. Both
update controls accept values from 1 to 100. The threshold is shown as a
percentage (ISOLDE's default is 0.05%).
The current local region is focused automatically during every MDFF trial, and
the original camera view is restored after the complete run. The tool always
asks for confirmation before starting ISOLDE simulations.

Here `#1` is an atomic structure and `#2` is the density map already associated
with that structure in ISOLDE.

Both workflows refuse to run if an ISOLDE simulation is active or if
`get_mdff_mgr(structure, volume, create=False)` returns `None`.
Automated local MDFF additionally requires that manager to be enabled and that
ISOLDE is already open.

## Algorithm guarantees

- No OpenMM minimization occurs while candidates are being tested.
- Only ISOLDE chi moving-atom sets are rotated; backbone coordinates are not
  moved during enumeration.
- Original residue coordinates are restored before every candidate and again
  before applying the winner.
- Only side-chain non-hydrogen atoms contribute to the density score.
- Any candidate with a scored atom outside the selected map is rejected.
- The MDFF manager is verified with `create=False`; scoring never creates one.
- ISOLDE simulation, minimization, MDFF state and restraints are not changed.

## Automated local-MDFF guarantees

- Candidates are tried in decreasing ISOLDE expected-frequency order.
- The full model coordinates accepted before a residue are restored before
  every candidate, so rejected MDFF motion cannot leak into the next trial.
- Only the local protein region is mobile; ISOLDE still creates its fixed
  physical context shell.
- Backbone N, CA, C and O restraints are temporary and their prior settings
  are restored after every trial.
- Trajectory smoothing is disabled during automated annealing. The user's
  original temperature and smoothing settings are restored after every trial.
- Candidate acceptance is assessed only after cooling and final minimization.
- Side chains are not position-restrained by this bundle.
- The selected map must have an existing, enabled MDFF manager. The workflow
  does not create or enable one.

## Development check

Build the bundle with the installed ChimeraX executable:

```text
/Applications/ChimeraX.app/Contents/bin/ChimeraX --nogui --exit \
  --cmd "devel build /path/to/CorrectRotamer"
```
