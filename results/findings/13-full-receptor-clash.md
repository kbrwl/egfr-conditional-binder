# Would each designed binder still fit in the intact receptor?

Computed output of `analysis/13_full_receptor_clash.py`. Do not hand-edit.

The design run is given the 171-residue fragment and never sees the other
450 residues, so it scores a binder that lies on an open face and one that
occupies space domain II or domain IV is already in exactly the same. This
step superposes each candidate back into the full 6ARU structure and measures
the overlap.

Step 12 settled that the anchor face is reachable, which is a property of the
target. This is a property of each design.

No design run has returned anything yet, so the evidence below is from the
constructed cases in section 2.

```
========================================================================
FULL-RECEPTOR CLASH CHECK — would each binder fit in the intact receptor?
========================================================================

The design run sees only the 171-residue fragment, so it cannot tell a binder
that lies on an open face from one that occupies space another part of the
receptor is already in. Both score the same. This step puts each candidate
back into the intact receptor and measures the overlap.

   Step 12 established that the anchor face itself is reachable. That is a
   property of the target. This is a property of each design: a reachable
   face does not stop one particular binder approaching at an angle that
   puts part of it inside domain II or domain IV.

   NOTE: --break-rule hard-clash is in effect — stop counting atom overlaps below 2.5 A, so a binder buried in domain IV reads clear.
   This run is expected to fail.

1. The intact receptor: 6aru.pdb

   receptor chain A, 99.7% identity to human EGFR
   residues the design run saw (domain III, 310-480): 171
   residues it did not see, which this step measures against: 438
   other chains in the file, left out of the measurement: B, C

   An overlap below 2.5 A is two sets of atoms in the same place.
   Between 2.5 and 4.5 A is touching, which the receptor
   may be able to move to accommodate, so it is reported as marginal rather
   than counted as a failure.

2. Constructed tests, run every time

   No design run has returned anything yet, so these are the only candidates
   this script has ever scored. Each one exercises one branch.

   | constructed candidate | what it is | expected | got | overlaps | RMSD |
   |---|---|---|---|---|---|
   | binder-on-open-face | a binder out from the anchor face, where a real one would sit | clear | clear | 0 | 0.00 A |
   | binder-inside-domain-iv | a binder in the volume the rest of the receptor occupies | clashing | marginal | 0 | 0.00 A |
   | target-wrong-shape | a target whose shape is not 6ARU's, so the binder cannot be placed | not scored | not scored | 0 | 5.98 A |
   >  target-wrong-shape: the target does not lie on 6ARU: RMSD 5.98 A over 171 CA atoms, above the 2.5 A limit. The binder cannot be placed, so no clash figure would mean anything
   | binder-only-file | a binder-only file, which the run's output really does contain | not scored | not scored | 0 | — |
   >  binder-only-file: no chain in this file reads as human EGFR, so which molecule is the target cannot be established

3. Real candidates

   None found under /Users/kunalbariwal/AI/EGFR Binder/egfr-conditional-binder/results/candidates. Nothing from the design run exists yet, which is why section 2
   is the whole of this script's evidence. Rerun this step with --candidates
   pointing at the run's output folder once there is one.

========================================================================
4. What this changes about the design
========================================================================

   A clashing candidate is not dropped here. It is flagged with the overlap
   measured, and the decision to discard belongs in step 10's ranking, so
   that every decision about what to submit is made in one place.

   A clear verdict does not mean the design works. It means the design is
   not ruled out by where its atoms are. Whether it binds, and whether it
   switches with pH, are what step 10 and the experiment decide.

   This is a rigid measurement: the receptor is treated as unable to move.
   A real receptor can shift to accommodate a small overlap, which is why a
   marginal verdict is reported as marginal rather than as a failure.

Wrote data/derived/13-candidate-clashes.csv (header only; no real candidates yet)

========================================================================
RESULT: 1 CHECK(S) FAILED —
  - binder-inside-domain-iv: expected clashing, got marginal
========================================================================
```
