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

1. 6ARU — extended (open) form

   receptor chain A, 609 residues mapped, 99.7% identity to human EGFR
   numbering offset: ours = file +24 (100.0% of residues) — derived, not assumed
   other protein chains present, excluded from the measurement: B, C

   [PASS] identity check against 6ARU: all 8 positions hold the amino acid our numbering says they should

   The current anchors, against the rest of the receptor:

   | anchor | touched at <4.5 A by | nearest outside domain III | within 8 A | within 12 A | reading |
   |---|---|---|---|---|---|
   | E344 | nothing | 24.1 A to 309 | 0 | 0 | open face, room for a binder |
   | H358 | nothing | 23.8 A to 309 | 0 | 0 | open face, room for a binder |
   | D368 | nothing | 13.0 A to 486 | 0 | 0 | open face, room for a binder |
   | H370 | nothing | 17.0 A to 486 | 0 | 0 | open face, room for a binder |
   | E391 | nothing | 17.4 A to 521 | 0 | 0 | open face, room for a binder |
   | E400 | nothing | 13.8 A to 516 | 0 | 0 | open face, room for a binder |
   | E421 | nothing | 4.6 A to 521 | 2 | 9 | clear of contact, some crowding |
   | E424 | nothing | 7.7 A to 533 | 1 | 7 | clear of contact, some crowding |

   None of the 8 anchors is touched from outside domain III.

   Control — the old 415-466 anchors through this same code.
   Step 06 measured these and found two of them in a groove against
   domain IV. If this script cannot reproduce that, its answer above
   is not to be believed either.

     D416: touched by nothing; nearest 5.5 A
     H418: touched by nothing; nearest 8.0 A
     E421: touched by nothing; nearest 4.6 A
     E424: touched by nothing; nearest 7.7 A
     H433: touched by nothing; nearest 11.3 A
     E455: touched by nothing; nearest 5.3 A
     D458: touched by 486 at 2.98 A; nearest 3.0 A
     D460: touched by 487 at 2.77 A; nearest 2.8 A

   [PASS] cross-check anchors in contact from outside domain III in 6ARU: 2 residues, identical to 06-intra-chain-occlusion.csv

1. 1NQL — tethered (closed) form

   receptor chain A, 612 residues mapped, 99.8% identity to human EGFR
   numbering offset: ours = file +24 (100.0% of residues) — derived, not assumed
   other protein chains present, excluded from the measurement: B

   [PASS] identity check against 1NQL: all 8 positions hold the amino acid our numbering says they should

   The current anchors, against the rest of the receptor:

   | anchor | touched at <4.5 A by | nearest outside domain III | within 8 A | within 12 A | reading |
   |---|---|---|---|---|---|
   | E344 | nothing | 23.6 A to 39 | 0 | 0 | open face, room for a binder |
   | H358 | nothing | 19.1 A to 39 | 0 | 0 | open face, room for a binder |
   | D368 | nothing | 14.8 A to 483 | 0 | 0 | open face, room for a binder |
   | H370 | nothing | 17.4 A to 489 | 0 | 0 | open face, room for a binder |
   | E391 | nothing | 19.3 A to 522 | 0 | 0 | open face, room for a binder |
   | E400 | nothing | 14.0 A to 516 | 0 | 0 | open face, room for a binder |
   | E421 | nothing | 7.1 A to 522 | 3 | 12 | clear of contact, some crowding |
   | E424 | nothing | 8.5 A to 516 | 0 | 5 | clear of contact, some crowding |

   None of the 8 anchors is touched from outside domain III.

   Control — the old 415-466 anchors through this same code.
   Step 06 measured these and found two of them in a groove against
   domain IV. If this script cannot reproduce that, its answer above
   is not to be believed either.

     D416: touched by nothing; nearest 5.1 A
     H418: touched by nothing; nearest 4.8 A
     E421: touched by nothing; nearest 7.1 A
     E424: touched by nothing; nearest 8.5 A
     H433: touched by nothing; nearest 10.8 A
     E455: touched by nothing; nearest 5.7 A
     D458: touched by 483 at 2.98 A; nearest 3.0 A
     D460: touched by 487 at 2.94 A; nearest 2.9 A


========================================================================
2. What this changes about the design
========================================================================

   Anchors touched from outside domain III, in either form: 0 — none
   Anchors with any neighbour inside 8 A, in either form: 2 — E421, E424
   Anchors on open face in both forms: 6 — E344, H358, D368, H370, E391, E400

   The face the design run is aimed at is not covered by the rest of
   the receptor, in either measured conformation. Cutting domain III out
   does not hide the anchors, so the occlusion risk that cost the old
   415-466 epitope two of its anchors does not apply to this face.

   What it does not license: this says the anchors are reachable, not
   that the fragment folds the way the intact protein does. That stays
   open, and only the design run's own predictor can close it.

   Worth carrying forward: E421, E424 have neighbours close by.
   Step 09 separately found E424 sits 6.6 A from the cut at 480. An
   anchor in both lists is the weakest of the eight on two
   independent grounds, and step 10 already ranks designs that lean on
   it below equivalent designs that do not.

Wrote data/derived/12-anchor-clearance.csv

========================================================================
RESULT: PASSED. The H370 anchor face is not occluded by the rest of the
receptor in either measured conformation. Whether the cut fragment holds
its shape is a separate question and is still open.
========================================================================
```
