# Conformation check: is the epitope reachable when EGFR is closed?

Computed output of `analysis/06_tethered_occlusion.py`. Do not hand-edit.

Compares epitope accessibility between 6ARU and 1NQL, keeping
two sources of occlusion apart: other parts of the receptor, which is
intrinsic to the molecule and would be present in the assay, and a bound
ligand, which is not.

The `1NQL` structure choice was UNVERIFIED and is now verified
against the RCSB entry record: it is an EGFR extracellular-region structure
with EGF (epidermal growth factor) bound, described by the depositors as
inactive, solved at low pH.

```
========================================================================
CONFORMATION CHECK — is 415-466 reachable in the closed form?
========================================================================

Compares per-residue accessibility of the epitope between an open-form
structure and a closed-form structure, using the receptor protein alone in
both, so the difference comes from the conformation and not from whichever
partner happens to be bound in each file.

  open   / reference: 6ARU
  closed / tethered:  1NQL

  6ARU: already present
  1NQL: already present

1. What 1NQL actually contains (was UNVERIFIED)

   HEADER    HORMONE/GROWTH FACTOR RECEPTOR          21-JAN-03   1NQL
   TITLE     STRUCTURE OF THE EXTRACELLULAR DOMAIN OF HUMAN EPIDERMAL GROWTH FACTOR
   TITLE    2 (EGF) RECEPTOR IN AN INACTIVE (LOW PH) COMPLEX WITH EGF.
   COMPND   2 MOLECULE: EPIDERMAL GROWTH FACTOR RECEPTOR;
   COMPND   3 CHAIN: A;

   6ARU: receptor chain A — 609 observed residues, 99.7% identity to human EGFR
        numbering offset: UniProt = PDB +24 (100.0% of residues) — derived, not assumed
        other protein chains present: B, C
   1NQL: receptor chain A — 612 observed residues, 99.8% identity to human EGFR
        numbering offset: UniProt = PDB +24 (100.0% of residues) — derived, not assumed
        other protein chains present: B

   Epitope coverage in 1NQL: 52 of 52 residues resolved

2. Is domain III the same fold in both structures?

   An accessibility comparison would be measuring two things at once if
   domain III itself were folded differently in the two files. Laying the
   two copies of the domain on top of each other and measuring RMSD — the
   average distance left between matched atoms, in angstroms — tells us
   whether that is the case.

   Superposed 171 matched CA atoms (the alpha carbon, the one atom every residue has) across domain III (310-480).
   RMSD = 1.08 A
   Domain III has essentially the same fold in both structures: the
   two copies sit on top of each other to within a couple of
   angstroms on average, which for a protein means the same shape.
   Any accessibility difference below is therefore caused by what
   surrounds the domain, and not by the domain rearranging.

2b. Do these two structures actually represent different conformations?

   The rest of this step only means something if the two structures are in
   different conformations. If they both happen to be in the same one, then
   comparing them tells us nothing about the tethered state, and even a
   reassuring result would carry no information. Hence this check.

   Method: superpose on domain III only (done above), then measure how far
   the rest of the molecule sits from its counterpart. If the global
   arrangement is the same, those displacements are small.

   | region | matched CA | RMSD after domain III superposition |
   |---|---|---|
   | domains I-II (25-309) | 282 | 23.44 A |
   | domain III (310-480) | 171 | 1.08 A |
   | domain IV (481-620) | 140 | 2.12 A |

   The two structures do differ in global conformation: with domain
   III superposed, the other domains sit a long way from their
   counterparts (domains I-II (25-309) at 23.4 A).
   So this is a comparison between two different arrangements, and the
   accessibility comparison below is meaningful.

3. Accessibility of the epitope, receptor protein only

   RSA, relative solvent accessibility, is how much of a residue's surface
   water can reach, divided by the most that residue type could ever
   expose. Near 0 means covered over; near 1 means out in the open. A
   binder can only grip what water can reach.

   Comparable residues: 52
   Mean RSA, open   (6ARU):   0.173
   Mean RSA, closed (1NQL): 0.189
   Change: +0.016 (+9.5% relative)

   The eight anchors specifically:

   | anchor | RSA open | RSA closed | change | reading |
   |---|---|---|---|---|
   | D416 | 0.117 | 0.171 | +0.055 | little change |
   | H418 | 0.032 | 0.200 | +0.168 | more accessible |
   | E421 | 0.210 | 0.145 | -0.064 | little change |
   | E424 | 0.429 | 0.448 | +0.019 | little change |
   | H433 | 0.648 | 0.759 | +0.111 | more accessible |
   | E455 | 0.348 | 0.423 | +0.075 | little change |
   | D458 | 0.248 | 0.446 | +0.198 | more accessible |
   | D460 | 0.186 | 0.120 | -0.066 | little change |

4. What covers the epitope in the closed form?

   Source 1 — other parts of the receptor chain. This is the tether itself,
   and it would be present in the assay. Measured as: residues outside
   310-480, same chain, within 4.5 A of a residue in 415-466.
   Putting it that way avoids naming the domains that do the covering, so
   the answer does not rest on domain boundaries we are unsure of. The only
   boundary it needs is domain III's own range.

   6ARU: 19 epitope residue(s) contacted from outside domain III
        446 <- 504 at 4.45 A
        447 <- 516 at 3.22 A
        448 <- 516 at 3.57 A
        449 <- 516 at 3.25 A
        450 <- 516 at 2.98 A
        451 <- 522 at 2.64 A
        453 <- 516 at 4.02 A
        454 <- 524 at 3.60 A
        456 <- 486 at 3.02 A
        457 <- 483 at 3.07 A
        458 <- 486 at 2.98 A  <-- ANCHOR
        459 <- 486 at 2.54 A
        460 <- 487 at 2.77 A  <-- ANCHOR
        461 <- 489 at 2.86 A
        462 <- 489 at 3.38 A
        463 <- 491 at 2.88 A
        464 <- 491 at 3.15 A
        465 <- 492 at 3.29 A
        466 <- 493 at 2.84 A
   1NQL: 18 epitope residue(s) contacted from outside domain III
        447 <- 516 at 3.46 A
        448 <- 516 at 2.36 A
        449 <- 516 at 4.14 A
        450 <- 516 at 3.24 A
        451 <- 522 at 3.63 A
        453 <- 516 at 4.10 A
        454 <- 524 at 3.80 A
        456 <- 481 at 3.33 A
        457 <- 483 at 2.78 A
        458 <- 483 at 2.98 A  <-- ANCHOR
        459 <- 481 at 3.63 A
        460 <- 487 at 2.94 A  <-- ANCHOR
        461 <- 487 at 3.37 A
        462 <- 489 at 3.82 A
        463 <- 489 at 2.93 A
        464 <- 491 at 3.42 A
        465 <- 492 at 3.14 A
        466 <- 493 at 2.86 A

   Contacted from outside domain III in both structures: 18 residues, spanning 447-466.
   Contacted only in the closed structure: none

   The two structures give almost the same list, so this packing is a
   standing feature of how the protein folds rather than something the
   closed shape introduces. The partner residues are in the 481-524 range,
   which is domain IV, the domain that follows ours.

   Anchors sitting against domain IV in both structures: D458, D460

   What this changes for the design. These anchors are still reachable by
   water: step 03 measured accessibility on the whole receptor chain with
   domain IV already present, so its effect is already inside those
   numbers. What this adds is that they sit in a groove between two
   domains rather than on an open face. A binder reaching them has to fit
   into that groove, which is a harder shape to design against than a
   flat surface, and it makes those contacts more sensitive to any shift
   in how the two domains sit against each other.

   Worth weighing when choosing between the candidate anchor clusters in
   step 05, which does not have this information: it runs before this step
   and reads only the exposure and antibody-overlap tables.

   Source 2 — the bound ligand, meaning EGF, the epidermal growth factor
   that EGFR normally responds to. The assay presents the receptor without
   EGF, so occlusion by the ligand must not be counted against us: doing so
   would make the tethered form look far worse than it is.

   6ARU (cetuximab Fab, chains B, C): 10 epitope residue(s) contacted
        432 at 2.60 A from chain C
        433 at 3.47 A from chain C  <-- ANCHOR
        435 at 3.88 A from chain C
        436 at 3.85 A from chain C
        439 at 4.13 A from chain C
        441 at 3.59 A from chain C
        442 at 3.26 A from chain C
        462 at 3.45 A from chain C
        464 at 2.74 A from chain C
        465 at 3.30 A from chain C

   [PASS] cross-check cetuximab contacts inside the epitope, against step 04: 10 residues, identical to 04-epitope-overlap.csv
   1NQL (EGF, chains B): 0 epitope residue(s) contacted

5. VERDICT

   Anchors compared: 8
   Still accessible in the closed form (RSA > 0.05): 8 — D416, H418, E421, E424, H433, E455, D458, D460
   Become buried: 0 — none
   Substantially less accessible: 0 — none

   Anchors the two structures disagree about:
     H418: buried in 6ARU (0.032) but partial in 1NQL (0.200)
     D458: partial in 6ARU (0.248) but exposed in 1NQL (0.446)

   These are ambiguous rather than resolved. Two measured structures give
   different answers, and this comparison cannot say which one reflects
   the molecule in our assay. Do not round either reading into a
   conclusion.

   Worth noting specifically: H418 read as
   buried in 6ARU but accessible in 1NQL. 6ARU is an
   antibody complex, so burial there may be an artefact of the
   antibody holding a side chain in place rather than a property of
   the receptor on its own. Step 03 excluded H418 on the
   6ARU reading alone, and that exclusion should now be treated
   as uncertain rather than settled. It matters because H418 is a
   target histidine, and the method-novelty claim rests on it.

   The epitope survives the conformation check.
   Accessibility of the block is similar in both structures and at least
   three anchors remain reachable in the closed form, so a binder aimed
   here does not depend on the receptor being open.

   The block is slightly more accessible in the closed structure than in
   the open one. The worry that prompted this step was that domain II
   folds across our face of domain III when the receptor closes, and the
   numbers do not show that happening.

   What they do show is that 18 residues in the second half
   of our block sit against domain IV, in both structures and to within a
   few tenths of an angstrom of the same distances. That is a standing
   feature of the fold rather than something closing introduces, and it is
   already reflected in the accessibility numbers. Two anchors, D458 and
   D460, are in that group, so they sit in a groove between two domains
   rather than on an open face, and a binder aimed at them has to fit that
   groove.

   Limits of this comparison:

   - Two crystal structures are two snapshots. Neither tells us the
     proportion of open to closed in the assay buffer, which is the number
     that actually matters and the one we do not have.
   - 1NQL was solved at low pH and with EGF bound. Both are
     circumstances of getting the crystal to form, and neither describes our
     assay.
   - The two structures differ in construct, resolution and crystallisation
     conditions as well as in conformation. Some of the difference measured
     above is attributable to those.
   - Accessibility computed on the bare protein ignores glycans, the sugar
     chains attached to the protein surface, and step 03 showed one is
     attached inside this very block at N444.

   This step reduces the tethering risk without eliminating it. The honest
   statement for the write-up is that the epitope is accessible in both
   published conformations we could test, and that the balance of
   conformations in the assay remains unknown.

Wrote data/derived/06-conformation-comparison.csv
Wrote data/derived/06-intra-chain-occlusion.csv

========================================================================
RESULT: 8 of 8 anchors remain accessible in 1NQL. Mean epitope RSA 0.173 -> 0.189.
========================================================================
```
