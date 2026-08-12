# CorrectRotamer for ChimeraX

CorrectRotamer is a ChimeraX bundle for automatic, density-guided correction
of ISOLDE rotamer outliers.

For each outlier present at the start of a run, the bundle restores the
original residue coordinates before trying every canonical ISOLDE target. It
discards targets that remain validation outliers or place any side-chain heavy
atom outside the map, then keeps the valid target with the highest unweighted
mean interpolated density. Rotamer probability is never used for ranking.

After enumeration, the bundle stops. It does not start, pause, or minimize an
ISOLDE simulation and does not create or modify position restraints. Final
refinement remains under the user's control in ISOLDE.

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

The graphical tool provides model and map selectors, live validation of the
existing ISOLDE MDFF association, run progress, and a residue-level results
table. Selecting a result row focuses that residue in ChimeraX and shows the
density score or rejection reason for every canonical target that was tested.

Here `#1` is an atomic structure and `#2` is the density map already associated
with that structure in ISOLDE.

The command refuses to run if an ISOLDE simulation is active or if
`get_mdff_mgr(structure, volume, create=False)` returns `None`.

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

## Development check

Build the bundle with the installed ChimeraX executable:

```text
/Applications/ChimeraX.app/Contents/bin/ChimeraX --nogui --exit \
  --cmd "devel build /path/to/CorrectRotamer"
```
