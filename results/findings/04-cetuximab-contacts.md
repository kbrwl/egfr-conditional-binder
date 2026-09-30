# Cetuximab contact set (computed)

Computed output of `analysis/04_cetuximab_contacts.py`. Do not hand-edit.

Settles a claim that had been carried as UNVERIFIED, meaning stated from
memory rather than computed: that the differences Q390R, E412D, R414W and
K467R explain why cetuximab does not bind mouse EGFR.

A residue counts as in contact when any of its heavy atoms, meaning any
atom except hydrogen, comes within
4.5 angstroms of any heavy atom of either chain of the antibody
fragment, in structure 6ARU. Hydrogen is left out because it is too light
to appear in most measured structures.

```
========================================================================
CETUXIMAB CONTACT SET — resolving an UNVERIFIED claim
========================================================================

Structure: 6ARU, receptor chain A, Fab chains B, C
Definition: an EGFR residue is in contact if any heavy atom is within
            4.5 A of any heavy atom of either Fab chain. A Fab is
            the gripping arm of an antibody, on its own; a heavy atom is
            any atom except hydrogen, and hydrogens are absent from this
            structure file.
Method: Bio.PDB.NeighborSearch builds a spatial index, so we avoid
        comparing all ~25 million atom pairs one by one. Same answer.
Numbering: UniProt position = PDB position + 24, measured by step 02.

EGFR residues in contact with the Fab: 24

1. Full contact list

   One row per EGFR residue cetuximab touches. 'aa' is the amino acid in
   one-letter code, and 'min dist' is the closest heavy-atom approach
   between that residue and the Fab.

   | UniProt | aa | PDB# | min dist A | Fab chain | EGFR atom | Fab atom |
   |---|---|---|---|---|---|---|
   | 373 | P | 349 | 3.42 | C | CB | O |
   | 374 | V | 350 | 4.47 | C | N | O |
   | 377 | R | 353 | 3.61 | C | NH2 | O |
   | 406 | L | 382 | 3.55 | C | CD1 | OH |
   | 408 | Q | 384 | 3.43 | C | O | ND2 |
   | 432 | Q | 408 | 2.60 | C | NE2 | OH |
   | 433 | H | 409 | 3.47 | C | ND1 | OH |
   | 435 | Q | 411 | 3.88 | C | OE1 | CD2 |
   | 436 | F | 412 | 3.85 | C | CE2 | CE2 |
   | 439 | A | 415 | 4.13 | C | CB | CE2 |
   | 441 | V | 417 | 3.59 | C | CG2 | CD1 |
   | 442 | S | 418 | 3.26 | C | OG | ND2 |
   | 462 | I | 438 | 3.45 | C | CD1 | CB |
   | 464 | S | 440 | 2.74 | C | OG | OH |
   | 465 | G | 441 | 3.30 | C | N | OH |
   | 467 | K | 443 | 3.39 | C | CE | OD1 |
   | 489 | K | 465 | 2.77 | C | NZ | OD2 |
   | 490 | I | 466 | 3.92 | B | O | ND2 |
   | 491 | I | 467 | 3.05 | B | CB | O |
   | 492 | S | 468 | 2.82 | B | N | O |
   | 493 | N | 469 | 2.70 | B | O | N |
   | 495 | G | 471 | 3.44 | B | CA | ND2 |
   | 496 | E | 472 | 4.24 | B | N | OD1 |
   | 497 | N | 473 | 3.08 | B | CB | OE1 |

   Footprint spans UniProt 373-497.
   Residues by nearest Fab chain: chain B: 7, chain C: 17

2. Question (a): do human/mouse differences fall in the contact set?

   Of the 16 domain III differences, 3 are in contact:
     R377K  —  3.61 A from Fab chain C
     S442G  —  3.26 A from Fab chain C
     K467R  —  3.39 A from Fab chain C

   The other 13 are not: A313P, S315Y, M318V, V323I, E330D, S348T, N361Y, S364A, H383R, Q390R, D393E, E412D, R414W

   Testing the four differences the earlier assumption named:
     Q390R: not in contact (nearest contact residue is 377)
     E412D: not in contact (nearest contact residue is 408)
     R414W: not in contact (nearest contact residue is 408)
     K467R: in contact — 3.39 A from chain C

   Verdict on the earlier assumption:
   Partly confirmed: 1 of 4 are in
   contact. The assumption pointed in the right direction but was
   inaccurate as written, so the decisions log now carries the
   measured set in its place.

   Confirmed in contact: K467R
   Not in contact: Q390R, E412D, R414W

   What this establishes, and what it does not. A differing residue
   sitting inside the footprint is consistent with it causing the
   species failure, and is far better evidence than an assumption.
   It is not proof: whether a given substitution actually abolishes
   binding depends on how much that particular contact contributes
   to the grip, which this calculation does not measure. So we can
   cite it as the likely explanation, but not as a settled one.

3. Question (b): overlap between cetuximab's footprint and 415-466

   Cetuximab touches 10 of the 52 residues in
   our block (19.2%).
   Of cetuximab's 24 total contacts, 10 (41.7%) are in our block.

   Which ones:
     Q432: 2.60 A from chain C
     H433: 3.47 A from chain C  <-- ANCHOR
     Q435: 3.88 A from chain C
     F436: 3.85 A from chain C
     A439: 4.13 A from chain C
     V441: 3.59 A from chain C
     S442: 3.26 A from chain C
     I462: 3.45 A from chain C
     S464: 2.74 A from chain C
     G465: 3.30 A from chain C

   Anchors inside cetuximab's footprint: 1 of 8 — H433

   What this means for the design:
   Partial overlap. Our block sits alongside the druggable surface and
   shares part of it, so a binder here plausibly blocks the same
   function, while most of our contact positions are ones cetuximab
   never touches. That makes the novelty claim defensible on residue
   identity alone, and we can keep the block as it stands.

4. Question (c): residues in 415-466 cetuximab leaves untouched

   42 of 52 residues are cetuximab-free.
   These are the positions our design can build on to look unlike
   cetuximab while still sitting on the same surface.

   415 416 417 418 419 420 421 422 423 424 425 426 427
   428 429 430 431 434 437 438 440 443 444 445 446 447
   448 449 450 451 452 453 454 455 456 457 458 459 460
   461 463 466

   Anchors cetuximab leaves untouched: 7 of 8 — D416, H418, E421, E424, E455, D458, D460

   These free anchors are the most directly useful output of this step. An
   anchor that is both reachable (step 03) and outside cetuximab's
   footprint lets us build the pH switch on a surface no approved drug
   occupies, so those anchors get first preference in the design.

Wrote data/derived/04-cetuximab-contacts.csv
Wrote data/derived/04-epitope-overlap.csv

========================================================================
RESULT: 24 EGFR residues contact the Fab. 3 of 16 species differences are among them.
Prior speculation: 1 of 4 confirmed.
Overlap with our block: 10/52 residues.
========================================================================
```
