# Is the fragment cut from the tethered structure different?

Computed output of `analysis/16_tethered_fragment.py`. Do not hand-edit.

The organisers say the screen uses the tethered form of EGFR. Our fragment is
cut from 6ARU, the extended form. This cuts the same residues, 310 to 480, from
1NQL, the tethered form, lays it on the 6ARU cut, and decides on the measured
movement of each anchor's charged tip.

```
========================================================================
TETHERED FRAGMENT — is the cut from 1NQL different where the anchors are?
========================================================================

The organisers say the screen uses the tethered form. Our fragment is cut
from 6ARU, the extended form. Step 12 found no anchor covered in either,
and step 06 found domain III superposes at 1.08 A. Neither says whether the
anchors' charged groups sit in the same places, which is what the pairing
rule is built on. This measures it, and decides on the measurement.

1. The deciding rule, tested every time

   The deciding rule, on constructed displacements of anchor tips (A):

     [PASS] no movement: same
     [PASS] the largest at 1.9 A, just under the tolerance: same
     [PASS] the largest at exactly the tolerance: differs
     [PASS] one anchor 5 A out, the rest still: differs

2. The two structures

   6ARU: receptor chain A, 609 residues mapped, 99.7% identity; ours = file +24 (100.0%), derived
   1NQL: receptor chain A, 612 residues mapped, 99.8% identity; ours = file +24 (100.0%), derived

   [PASS] identity check, anchors, 6ARU: all 8 hold the amino acid our numbering says
   [PASS] identity check, anchors, 1NQL: all 8 hold the amino acid our numbering says

3. The fragment, residues 310-480 (171 positions)

   6ARU: 171 residues present
   1NQL: 171 residues present
   shared by both and used for the fit: 171

4. The whole fragment: 1NQL laid on 6ARU over 171 CA atoms

   RMSD 1.08 A
   [PASS] control against step 06: 1.08 A over 171 residues here, 1.08 A over 171 there
   The six residues that deviate most (CA, after the fit): G322 4.5 A, V323 4.1 A, A310 3.5 A, R324 3.3 A, E319 3.2 A, M318 3.0 A

5. The anchor face

   79 residues have an atom within 8 A of an anchor's charged tip in 6ARU.
   Their CA RMSD under the whole-fragment fit: 0.69 A (whole fragment 1.08 A)

   Fitting on the 8 anchors' CA atoms alone gives 0.55 A. A fit on 8 points is flattering by construction,
   so this is shown and not relied on; the whole-fragment fit above is the one the
   decision uses.

6. Each anchor, in the common frame

   | anchor | CA moved | atom compared | that atom moved | basis |
   |---|---|---|---|---|
   | E344 | 1.79 A | CB | 2.14 A | the tip is not resolved in 1NQL, so the nearer atom CB is compared in both |
   | H358 | 0.37 A | imidazole centroid | 0.32 A | charged tip, resolved in both |
   | D368 | 0.46 A | CG | 0.70 A | charged tip, resolved in both |
   | H370 | 0.26 A | imidazole centroid | 0.67 A | charged tip, resolved in both |
   | E391 | 0.30 A | CD | 1.19 A | charged tip, resolved in both |
   | E400 | 0.64 A | CD | 0.97 A | charged tip, resolved in both |
   | E421 | 0.37 A | CD | 0.80 A | charged tip, resolved in both |
   | E424 | 0.60 A | CD | 0.42 A | charged tip, resolved in both |

   Not measurable at the charged tip: E344. The atom compared is the same in both structures, so the figure is a real movement
   of that atom. The tip sits further out along the side chain and usually moves at least
   as much, though a side chain can rotate so that it moves less, and this does not
   establish which. These anchors are excluded from the deciding rule, which is about
   tips, and listed as a gap.

   Charged tip measured for 7 of 8 anchors. Largest move among them: 1.19 A (E391); mean 0.72 A; tolerance 2.0 A

========================================================================
7. Decision, and what it changes about the design
========================================================================

   The two cuts are the same where the anchors' charged tips can be
   compared, by the rule fixed above: no tip moves by 2.0 A or more. The fragment stays
   cut from 6ARU. The reason is the measurement and not the argument that domain III
   is a rigid blob.

   The verdict covers 7 of 8 anchors. Not covered: E344 (CA moved 1.79 A, CB moved 2.14 A).
   The charged tip of each is unresolved in one structure, so the rule could not be
   applied to it. This is a stated gap in the verdict and not a finding of no movement.

   The tethered cut is not needed for the design run. It is on the shelf at
   data/structures/1nql_domain3.pdb in case the design run's own predictor
   disagrees.

   H418 is a separate matter and changes nothing for the H370 face: it reads as
   buried in 6ARU (RSA 0.032) and partly exposed in 1NQL (0.200), and 6ARU being an
   antibody complex is the suspected cause. If the tethered form is what is tested,
   the 1NQL reading is the relevant one. That belongs to the 415-466 fallback.

   What this does not settle: it compares two crystal structures of a molecule that
   moves, at 3.2 A and 2.8 A resolution. It says nothing about whether the cut
   fragment keeps its shape once separated, which is the open question from step 11.

Wrote data/derived/16-anchor-displacement.csv and data/structures/1nql_domain3.pdb

========================================================================
RESULT: PASSED. verdict same for 7 of 8 anchors; largest charged-tip move 1.19 A; fragment RMSD 1.08 A
========================================================================
```
