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
   | binder-inside-domain-iv | a binder in the volume the rest of the receptor occupies | clashing | clashing | 55 | 0.00 A |
   | target-wrong-shape | a target whose shape is not 6ARU's, so the binder cannot be placed | not scored | not scored | 0 | 5.98 A |
   >  target-wrong-shape: the target does not lie on 6ARU: RMSD 5.98 A over 171 CA atoms, above the 2.5 A limit. The binder cannot be placed, so no clash figure would mean anything
   | binder-only-file | a binder-only file, which the run's output really does contain | not scored | not scored | 0 | — |
   >  binder-only-file: no chain in this file reads as human EGFR, so which molecule is the target cannot be established

3. Real candidates

   30 file(s) under results/candidates/egfr-r1-short-r1-short/2_Refolded

   | design | verdict | overlaps <2.5 A | contacts | closest | nearest residue | interface to domain III |
   |---|---|---|---|---|---|---|
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate10_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate10_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate1_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate1_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate2_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate2_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate3_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate3_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate4_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate4_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate5_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate5_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate6_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate6_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate7_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate7_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate8_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate8_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate9_monomer | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate9_monomer: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate10_EGFR_domain3 | clear | 0 | 0 | 11.15 A | 489 | 5 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate10_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate10_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate1_EGFR_domain3 | clear | 0 | 0 | 9.77 A | 491 | 76 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate1_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate1_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate2_EGFR_domain3 | clear | 0 | 0 | 10.41 A | 491 | 4 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate2_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate2_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate3_EGFR_domain3 | clear | 0 | 0 | 10.33 A | 489 | 2 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate3_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate3_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate4_EGFR_domain3 | clear | 0 | 0 | 11.06 A | 492 | 1 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate4_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate4_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate5_EGFR_domain3 | clashing | 1 | 8 | 1.96 A | 489 | 9 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate5_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate5_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate6_EGFR_domain3 | clear | 0 | 0 | 10.27 A | 491 | 9 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate6_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate6_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate7_EGFR_domain3 | clear | 0 | 0 | 9.92 A | 491 | 62 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate7_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate7_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate8_EGFR_domain3 | clear | 0 | 0 | 11.56 A | 489 | 55 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate8_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate8_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate9_EGFR_domain3 | clear | 0 | 0 | 10.51 A | 489 | 6 |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate9_HisTag | not scored | 0 | 0 | — | — | 0 |
   >  egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate9_HisTag: no chain in this file reads as human EGFR, so which molecule is the target cannot be established

   clear: 9
   marginal: 0
   clashing: 1
   not scored: 20

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

Wrote data/derived/13-candidate-clashes.csv

========================================================================
RESULT: PASSED. Every constructed case came back as it had to: a binder
on the open face reads clear, one inside domain IV reads clashing, and a
target that is not 6ARU is refused rather than scored.
========================================================================
```
