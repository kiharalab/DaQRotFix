# CorrectRotamer for ChimeraX

CorrectRotamer is a ChimeraX bundle for automatic density-guided and local-MDFF
optimization of ISOLDE rotamer outliers.

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
residue, receives a short local ISOLDE MDFF minimization, and is scored by the
number of rotamer outliers in the local protein region.

The first candidate with zero local outliers is accepted. If all candidates
leave outliers, the candidate with the fewest is retained; probability order
breaks ties. Rejected local-MDFF coordinates are restored. Protein backbone
N, CA, C and O atoms are strongly and temporarily position-restrained while
side-chain atoms remain mobile.

## Install

ISOLDE must already be installed. From the ChimeraX command line:

```text
devel install /path/to/CorrectRotamer
```

Restart ChimeraX after the first installation if the tool is not immediately
listed. Open **Tools → Validation → Correct Rotamers**, or run:

```text
correctrotamers #1 #2
```

The local-MDFF workflow is also available as:

```text
correctrotamers mdff #1 #2 radius 5 coordinateUpdates 3 backboneSpringConstant 500
```

The graphical tool provides separate buttons for density-only correction and
automated local MDFF. The latter has controls for local radius, MDFF coordinate
updates and backbone restraint strength, and always asks for confirmation
before starting ISOLDE simulations.

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
- Side chains are not position-restrained by this bundle.
- The selected map must have an existing, enabled MDFF manager. The workflow
  does not create or enable one.

## Development check

Build the bundle with the installed ChimeraX executable:

```text
/Applications/ChimeraX.app/Contents/bin/ChimeraX --nogui --exit \
  --cmd "devel build /path/to/CorrectRotamer"
```
