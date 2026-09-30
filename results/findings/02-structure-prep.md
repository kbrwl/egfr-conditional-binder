# Structure preparation: PDB 6ARU

Computed output of `analysis/02_structure_prep.py`. Do not hand-edit.

Resolves the PDB-to-UniProt numbering offset empirically, identifies the
chains from the file's own annotation and from measured identity to human
EGFR, reports which residues are actually resolved, and writes a
receptor-only structure for the solvent-accessibility step.

```
========================================================================
STRUCTURE PREPARATION — PDB 6ARU
========================================================================

0. Acquire structure
   https://files.rcsb.org/download/6ARU.pdb
   already present (1,372,383 bytes) -> data/structures/6aru.pdb
   already present (1,548,867 bytes) -> data/structures/6aru.cif  (mmCIF backup)

1. What the file says about itself
   HEADER    TRANSFERASE/IMMUNE SYSTEM               23-AUG-17   6ARU
   TITLE     STRUCTURE OF CETUXIMAB FAB MUTANT IN COMPLEX WITH EGFR EXTRACELLULAR
   TITLE    2 DOMAIN

   Chain annotation from the file's own COMPND records
   (the file's self-description — not our inference):
     chain A: EPIDERMAL GROWTH FACTOR RECEPTOR
     chain B: CETUXIMAB MUTANT LIGHT CHAIN,UNCHARACTERIZED PROTEIN
     chain C: CETUXIMAB MUTANT HEAVY CHAIN FAB FRAGMENT,IMMUNOGLOBULIN GAMMA-1 HEAVY CHAIN

2. Chain inventory, measured from the coordinates

   'Observed' means residues with actual coordinates. Residues listed in
   the file's sequence but too mobile to locate are absent here.

   chain A:  609 observed residues, numbered 4-612,  99.7% identity to human EGFR
              annotation: EPIDERMAL GROWTH FACTOR RECEPTOR
   chain B:  210 observed residues, numbered 1-210,   3.8% identity to human EGFR
              annotation: CETUXIMAB MUTANT LIGHT CHAIN,UNCHARACTERIZED PROTEIN
   chain C:  212 observed residues, numbered 1-219,   3.8% identity to human EGFR
              annotation: CETUXIMAB MUTANT HEAVY CHAIN FAB FRAGMENT,IMMUNOGLOBULIN GAMMA-1 HEAVY CHAIN

3. Which chain is the receptor?

   Decided by measurement, not by trusting the competition page.
  [PASS] receptor identified by EGFR identity: chain A at 99.7% — next best B at 3.8%
  [PASS] competition page claim 'chain A is the receptor': computed receptor is chain A — claim confirmed

   Fab chains (everything not the receptor): B, C
     chain B: CETUXIMAB MUTANT LIGHT CHAIN,UNCHARACTERIZED PROTEIN — 210 residues
     chain C: CETUXIMAB MUTANT HEAVY CHAIN FAB FRAGMENT,IMMUNOGLOBULIN GAMMA-1 HEAVY CHAIN — 212 residues

   Heavy vs light assignment is taken from the file's COMPND annotation
   above. Steps 03-04 do not depend on which is which: contacts are
   computed against 'either Fab chain', and both are removed together
   for the receptor-only file.

4. THE NUMBERING OFFSET, derived empirically

   Method: align the receptor chain's observed sequence against the full
   human UniProt sequence, then for every observed residue compute
      offset = UniProt position - PDB residue number
   A single dominant value means one consistent convention.

   | offset | aligned residues | share |
   |---|---|---|
   | +24 | 609 | 100.0% |

  [PASS] offset is consistent across the chain: 100.0% of aligned residues share offset +24
   Interpretation: PDB numbers by the MATURE protein. UniProt = PDB + 24. This is the trap the brief warned about.

5. Verification: are the eight anchors the residues we expect?

   This is the check that catches an offset error. If the offset were
   wrong by 24, these would read as the wrong amino acids.

  [PASS] anchor D416: PDB A/ASP392 reads D
  [PASS] anchor H418: PDB A/HIS394 reads H
  [PASS] anchor E421: PDB A/GLU397 reads E
  [PASS] anchor E424: PDB A/GLU400 reads E
  [PASS] anchor H433: PDB A/HIS409 reads H
  [PASS] anchor E455: PDB A/GLU431 reads E
  [PASS] anchor D458: PDB A/ASP434 reads D
  [PASS] anchor D460: PDB A/ASP436 reads D

6. Which residues are actually resolved?

   Receptor chain A, in UniProt numbering:
   spans 28-636, 609 residues observed in 1 continuous segment(s)
     28-636  (609 residues)

   Inside our epitope 415-466:
  [PASS] epitope fully resolved: all 52 residues have coordinates

7. Writing receptor-only structure for step 03

   Step 03 must compute solvent accessibility on the receptor ALONE.
   Cetuximab is sitting on the very surface we care about, so computing
   on the complex would report our epitope as buried when it is merely
   covered by an antibody that will not be present in our assay.

  [PASS] receptor-only file written: 609 residues written to data/structures/6aru_receptor_only.pdb

   Wrote data/derived/02-chain-inventory.csv
   Wrote data/derived/02-numbering-offset.csv
   Wrote data/derived/02-anchor-verification.csv

========================================================================
RESULT: PASSED.
Receptor is chain A. UniProt position = PDB residue number +24.
All eight anchors verified as the expected amino acids.
========================================================================
```
