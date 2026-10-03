# The His-tag counter-screen

Computed output of `analysis/15_histag_counterscreen.py`. Do not hand-edit.

Both targets carry a C-terminal His tag the organisers expect to leave on, and a
binder that grips it looks pH-selective and binds anything with a His tag. Written
before any candidate exists. No candidate has been screened: the evidence below is
the rule's own tests and a check of the campaign files.

```
========================================================================
HIS-TAG COUNTER-SCREEN — does any candidate grip the tag on the target?
========================================================================

Both targets carry a C-terminal His tag the organisers expect to leave on.
A binder that grips it looks pH-selective and binds anything with a His
tag. Our pairing rule builds acidic pockets to grip the target's own
histidines, which a tag's histidines fit as readily.

1. What the screen rests on

   Predictor: BindCraft2's own AlphaFold2, through its detargeting
   objective, with the tag given as a sequence target at a negative
   weight. Not `bindcraft score`, which reads existing structures and
   cannot predict a complex.
   Accepts it: i_pTM_detarget of at least 0.4, BindCraft2's own default ceiling
   for rejecting an off-target. That is the tool's default and not a
   threshold calibrated for this tag. Calibrating one needs known tag
   binders and known non-binders, and neither exists yet.
   The tag: GSHHHHHHGS, weight -0.5. The serine-glycine flank is there
   because BindCraft2 samples a window of at least 10 residues from a
   sequence target by default, and six residues is shorter than that.

2. The rule, tested every time

   The rule, on constructed i_pTM_detarget readings (ceiling 0.4):

     [PASS] '0.12'  -> no tag binding predicted   (well below the ceiling)
     [PASS] '0.39'  -> no tag binding predicted   (just below the ceiling)
     [PASS] '0.40'  -> accepts the tag   (exactly at the ceiling)
     [PASS] '0.71'  -> accepts the tag   (well above the ceiling)
     [PASS] ''      -> not recorded   (an empty cell)
     [PASS] none    -> not recorded   (no column at all)
     [PASS] 'nan'   -> not recorded   (not a number)

   Reading a constructed campaign folder, end to end:

     [PASS] design_a1b2c3_seq1: no tag binding predicted
     [PASS] design_d4e5f6_seq1: accepts the tag
     [PASS] design_g7h8i9_seq1: not recorded
     [PASS] design_j0k1l2_seq1: accepts the tag

3. The campaign files this depends on

   [PASS] egfr-domain3-h370.json: termini_accessible is on
   [PASS] egfr-domain3-h370-notag.json: termini_accessible is on
   [PASS] egfr-domain3-h370-notag.json: one target named HisTag
   [PASS] the tag target has a negative weight (-0.5), so it is avoided
   [PASS] the tag sequence file exists (his-tag-offtarget.fasta)
   [PASS] the tag sequence GSHHHHHHGS contains six consecutive histidines
   [PASS] egfr-domain3-h370.json: no off-targets, so it is the plain campaign
   [PASS] both files give the EGFR target and hotspots identically

4. Real candidates

   Metrics table: results/candidates/egfr-r1-short-r1-short/2_Refolded/!_Refolded.csv

   | design | i_pTM_detarget | verdict |
   |---|---|---|
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate1 | 0.08 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate10 | 0.11 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate2 | 0.10 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate3 | 0.10 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate4 | 0.08 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate5 | 0.08 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate6 | 0.09 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate7 | 0.10 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate8 | 0.09 | no tag binding predicted |
   | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate9 | 0.09 | no tag binding predicted |

   accepts the tag: 0
   no tag binding predicted: 10
   not recorded: 0

========================================================================
5. What this changes about the design
========================================================================

   Use `egfr-domain3-h370-notag.json` for the campaign, not the plain file, if
   the tag-avoiding variant starts. Designs generated without the tag
   cannot be screened afterwards without a prediction run that does not
   yet exist, so the choice has to be made before the campaign and not
   after it.

   Not run: detargeting against a sequence target has not been tried on
   real hardware, and the weight and the ceiling are starting values.
   A clear verdict is weak evidence, since a prediction that finds no
   interface does not show there is none. It can show a risk is present
   and cannot show it is absent.

Wrote data/derived/15-tag-screen.csv

========================================================================
RESULT: PASSED. 10 candidates screened.
========================================================================
```
