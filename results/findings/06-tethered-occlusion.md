# Conformation check: is the epitope reachable when EGFR is closed?

Computed output of `analysis/06_tethered_occlusion.py`. Do not hand-edit.

Compares epitope accessibility between 6ARU and 1NQL,
separating occlusion by other parts of the receptor (intrinsic, would be
present in the assay) from occlusion by a bound ligand (not intrinsic).

`1NQL` structure choice was UNVERIFIED and is now verified against
the RCSB entry record: it is an EGFR extracellular-region structure with
EGF bound, described by the depositors as inactive, solved at low pH.

```
========================================================================
CONFORMATION CHECK — is 415-466 reachable in the CLOSED form?
========================================================================

Compares per-residue accessibility of the epitope between an OPEN-form
structure and a CLOSED-form structure, receptor protein only in both,
so the difference isolates conformation rather than bound partners.

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

   If domain III itself were folded differently, an accessibility
   comparison would be measuring two things at once. Superposing the
   domain and reporting RMSD (the average leftover distance between
   matched atoms, in angstroms) separates those.

   Superposed 171 matched CA atoms across domain III (310-480).
   RMSD = 1.08 A
   Domain III has essentially the same fold in both structures.
   Any accessibility difference below is therefore caused by what
   surrounds the domain, not by the domain rearranging.

2b. Do these two structures actually represent DIFFERENT conformations?

   This check is necessary, not decorative. If both structures happen to
   be in the same conformation, then comparing them tells us NOTHING
   about the tethered state, and a reassuring result would be vacuous.

   Method: superpose on domain III only (done above), then measure how
   far the REST of the molecule sits from its counterpart. If the global
   arrangement is the same, those displacements are small.

   | region | matched CA | RMSD after domain III superposition |
   |---|---|---|
   | domains I-II (25-309) | 282 | 23.44 A |
   | domain III (310-480) | 171 | 1.08 A |
   | domain IV (481-620) | 140 | 2.12 A |

   The two structures DO differ in global conformation: with
   domain III superposed, other domains sit far from their
   counterparts (domains I-II (25-309) at 23.4 A).
   So this IS a comparison between two different arrangements, and
   the accessibility comparison below is meaningful.

3. Accessibility of the epitope, receptor protein only

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

   Source 1 — OTHER PARTS OF THE RECEPTOR CHAIN. This is the tether
   itself and would be present in the assay. Measured as: residues
   outside 310-480, same chain, within 4.5 A of
   a residue in 415-466. No domain-boundary definition
   needed beyond domain III's own range.

   6ARU: 1 epitope residue(s) contacted from outside domain III
        461 <- 487 at 3.10 A
   1NQL: 1 epitope residue(s) contacted from outside domain III
        461 <- 487 at 3.61 A

   Source 2 — THE BOUND LIGAND. Not intrinsic: the assay presents the
   receptor without EGF, so this must NOT be counted against us.

   6ARU (cetuximab Fab, chains B, C): 11 epitope residue(s) contacted
        416 at 2.74 A from chain C  <-- ANCHOR
        417 at 3.30 A from chain C
        419 at 3.39 A from chain C
        441 at 2.77 A from chain C
        442 at 3.92 A from chain B
        443 at 3.05 A from chain B
        444 at 2.82 A from chain B
        445 at 2.70 A from chain B
        447 at 3.44 A from chain B
        448 at 4.24 A from chain B
        449 at 3.08 A from chain B
   1NQL (EGF, chains B): 0 epitope residue(s) contacted

5. VERDICT

   Anchors compared: 8
   Still accessible in the closed form (RSA > 0.05): 8 — D416, H418, E421, E424, H433, E455, D458, D460
   Become buried: 0 — none
   Substantially less accessible: 0 — none

   ANCHORS THE TWO STRUCTURES DISAGREE ABOUT:
     H418: buried in 6ARU (0.032) but partial in 1NQL (0.200)
     D458: partial in 6ARU (0.248) but exposed in 1NQL (0.446)

   These are AMBIGUOUS, not resolved. Two experimental structures
   give different answers, and this comparison cannot say which
   reflects the molecule in our assay. Do not round either reading
   into a conclusion.

   Worth noting specifically: H418 read as
   buried in 6ARU but accessible in 1NQL. Since
   6ARU is an antibody complex, burial there may be an
   artefact of that antibody holding a side chain in place
   rather than an intrinsic property. Step 03 excluded H418 on
   the 6ARU reading alone; that exclusion should now be
   treated as UNCERTAIN rather than settled. It matters, because
   H418 is a target histidine and therefore carries the
   method-novelty claim.

   THE EPITOPE SURVIVES THE CONFORMATION CHECK.
   Accessibility of the block is similar in both structures and at
   least three anchors remain reachable in the closed form, so a
   binder aimed here is not dependent on the receptor being open.

   Note the direction of the result: the block is slightly MORE
   accessible in the closed structure, not less. The specific fear
   that motivated this step -- domain II folding across our face of
   domain III -- is not borne out. Only one epitope residue (461) is
   contacted from outside domain III, and that is true in BOTH
   structures, so it is a fixed feature of the fold rather than
   something the tether introduces.

   LIMITS OF THIS COMPARISON, stated plainly:

   - Two crystal structures are two snapshots. Neither tells us the
     PROPORTION of open to closed in the assay buffer, which is the
     number that actually matters and which we do not have.
   - 1NQL was solved at low pH and with EGF bound. Both are
     crystallisation circumstances, not descriptions of our assay.
   - The two structures differ in construct, resolution and
     crystallisation conditions, not only conformation. Some of the
     difference measured above is attributable to those.
   - Accessibility computed on the bare protein ignores glycans, which
     step 03 showed are attached inside this very block at N444.

   This step reduces a risk; it does not eliminate it. The honest
   statement for the write-up is that the epitope is accessible in both
   published conformations we could test, and that the conformational
   equilibrium in the assay remains unknown.

Wrote data/derived/06-conformation-comparison.csv

========================================================================
RESULT: 8 of 8 anchors remain accessible in 1NQL. Mean epitope RSA 0.173 -> 0.189.
========================================================================
```
