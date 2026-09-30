# Is the antibody in 6ARU mutated at the interface?

Computed output of `analysis/07_fab_mutant_check.py`. Do not hand-edit.

6ARU is titled a cetuximab Fab mutant, and step 04 measured the
antibody footprint using it. This checks whether any modified residue is
in the interface, which decides whether that footprint describes cetuximab
or only this modified version of it.

Reference for comparison: 1YY9, Li S. et al. (2005), "Structural
basis for inhibition of the epidermal growth factor receptor by
cetuximab", Cancer Cell 7:301-311, doi:10.1016/j.ccr.2005.03.003.

```
========================================================================
IS THE ANTIBODY IN 6ARU MUTATED AT THE INTERFACE?
========================================================================

Compares the two antibody chains of 6ARU, which is titled a
cetuximab Fab mutant, against 1YY9, the reference structure of
cetuximab (Li et al. 2005, Cancer Cell 7:301-311). Fab means the gripping
arm of an antibody, separated from the rest of it.

The mutations are not recorded anywhere we could read them: the file's own
SEQADV records, which list differences from a reference sequence, exist
only for the receptor chain; the Protein Data Bank entry names the mutant
but lists no substitutions; and the citation is 'To Be Published'.

   6ARU: already present
   1YY9: already present

1. Which chain is which, decided by measured identity to human EGFR

   6ARU:
     chain A: 609 residues,  99.7% identity to human EGFR  -> receptor
     chain B: 210 residues,   3.8% identity to human EGFR  -> antibody
     chain C: 212 residues,   3.8% identity to human EGFR  -> antibody
   1YY9:
     chain A: 613 residues,  99.5% identity to human EGFR  -> receptor
     chain D: 220 residues,   5.0% identity to human EGFR  -> antibody
     chain C: 211 residues,   3.8% identity to human EGFR  -> antibody

2. Pairing the antibody chains between the two files

   Matched by sequence similarity rather than by chain letter, because
   the same letter need not mean the same molecule in two files.

   6ARU chain B <-> 1YY9 chain C: 99.0% identical
   6ARU chain C <-> 1YY9 chain D: 98.6% identical

3. Antibody residues that touch the receptor in 6ARU

   Any antibody heavy atom within 4.5 angstroms of any
   receptor heavy atom, the same rule step 04 uses.
   Interface residues found: 20

4. Differences between the two antibodies

   chain B against 1YY9 chain C: 2 difference(s) over 210 aligned positions
     B 52: S in 1YY9 -> D in 6ARU
     B 56: S in 1YY9 -> D in 6ARU

   chain C against 1YY9 chain D: 3 difference(s) over 212 aligned positions
     C 28: S in 1YY9 -> D in 6ARU
     C 31: N in 1YY9 -> D in 6ARU  AT INTERFACE: 3.93 A from EGFR residue 433
     C 216: R in 1YY9 -> K in 6ARU

5. Verdict

   1 of 5 difference(s) sit within
   4.5 angstroms of the receptor:
     chain C 31: N -> D, 3.93 A from EGFR residue 433

   What this means for the project: step 04's contact set describes
   the modified antibody rather than cetuximab itself. Three results
   that rest on it need qualifying in the decisions log: the disproof
   of the earlier species-failure claim, the rule not to contact
   position 442, and the overlap figures in step 05's cluster
   comparison. Check in particular whether the EGFR residues listed
   above fall inside 415-466.

   Mutated antibody residues contact these residues inside our
   epitope: [433]

5b. Does the unmodified antibody touch the same residues?

   Recomputes the footprint on 1YY9, which has no engineered
   substitution at the interface, and compares the two inside our
   epitope 415-466. Both use the same rule: any receptor
   heavy atom within 4.5 angstroms of any antibody heavy
   atom.

   1YY9 receptor chain A: numbering offset
   UniProt = PDB +24, covering 100.0% of matched residues, derived by alignment.

   Inside 415-466:
     touched in both: [432, 433, 435, 436, 439, 441, 442, 462, 464, 465]
     touched only in 6ARU (the mutant): []
     touched only in 1YY9 (the reference): []

   Anchors touched in 6ARU: ['H433']
   Anchors touched in 1YY9: ['H433']

     H433: 6ARU 3.47 A; 1YY9 3.36 A

   Both structures touch the same anchors, so the overlap figures in
   step 05's cluster comparison hold for cetuximab and not only for
   this variant.

6. Limits

   The two structures were solved separately, at different resolutions
   and in different crystals, so a residue right at the cutoff can fall
   on either side for reasons unconnected to the mutation. Treat a
   difference of a few tenths of an angstrom as noise and a residue
   appearing in one list but not the other as worth checking rather than
   settled.

   1YY9's entry record does not state that its Fab is
   unmodified. It is used as the reference because it is the structure
   published with the account of how cetuximab works. A difference found
   here is a difference between two deposited structures, which is good
   evidence of a deliberate change where it is a single clean
   substitution and weaker where the files merely disagree.

   The comparison covers only residues both files resolve. A mutation in
   a stretch missing from either structure would not appear here.

   The interface test uses the positions in this one crystal structure.
   It answers whether a mutated residue touches the receptor in this
   structure, which is the question that matters for step 04, and does
   not measure how much any contact contributes to binding.

Wrote data/derived/07-fab-differences.csv

========================================================================
RESULT: 5 antibody difference(s), 1 at the interface.
========================================================================
```
