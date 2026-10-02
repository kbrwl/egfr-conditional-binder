# Is the H370 anchor face reachable in the intact receptor?

Computed output of `analysis/12_fragment_context_check.py`. Do not hand-edit.

Step 09 cuts domain III out and gives the design run the fragment alone, so
the design run cannot see whether the rest of the receptor folds across the
face we are aiming at. Step 06 asked this of the old 415-466 epitope and found
two anchors in a groove against domain IV. Step 08 then moved the face to the
cluster around H370, and six of those eight anchors had never been checked.

```
========================================================================
FRAGMENT CONTEXT CHECK — is the H370 face reachable in intact EGFR?
========================================================================

Step 09 hands the design run a 171-residue fragment, residues 310-480,
so the design run never sees the other 450 residues of the extracellular
region. This step asks whether those 450 residues fold across the face we
are aiming at. If they do, a design can score well and still be unable to
reach its target, and nothing in the design run would say so.

   anchors under test (from step 08's committed output): E344, H358, D368, H370, E391, E400, E421, E424
   fragment the design run receives: 310-480

Two measurements per anchor, because they answer different questions:
  contacts at 4.5 A   — is this residue covered over
  residues within 8 A and 12 A — can a body the size of a binder get here

   NOTE: --break-numbering 1 is in effect. This run is
   expected to fail at the identity check below.

1. 6ARU — extended (open) form

   receptor chain A, 609 residues mapped, 99.7% identity to human EGFR
   numbering offset: ours = file +24 (100.0% of residues) — derived, not assumed
   other protein chains present, excluded from the measurement: B, C

   [FAIL] identity check against 6ARU: E344: file has G at that position; H358: file has K at that position; D368: file has G at that position; H370: file has L at that position; E391: file has Q at that position; E400: file has K at that position; E421: file has F at that position; E424: file has L at that position

1. 1NQL — tethered (closed) form

   receptor chain A, 612 residues mapped, 99.8% identity to human EGFR
   numbering offset: ours = file +24 (100.0% of residues) — derived, not assumed
   other protein chains present, excluded from the measurement: B

   [FAIL] identity check against 1NQL: E344: file has G at that position; H358: file has K at that position; D368: file has G at that position; H370: file has L at that position; E391: file has Q at that position; E400: file has K at that position; E421: file has F at that position; E424: file has L at that position

========================================================================
2. What this changes about the design
========================================================================

Wrote data/derived/12-anchor-clearance.csv

========================================================================
RESULT: 2 CHECK(S) FAILED —
  - identity check failed against 6ARU
  - identity check failed against 1NQL
========================================================================
```
