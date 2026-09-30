# Numbering convention check

Computed output of `analysis/00_numbering_check.py`. Do not hand-edit.

All residue numbers in this project are positions in the full record for
human EGFR in UniProt, the public archive of protein sequences, where that
record is P00533. The official challenge constructs are the mature
extracellular region, UniProt 25-645, so:

    UniProt position = challenge-construct position + 24

```
========================================================================
NUMBERING CONVENTION CHECK
========================================================================

Convention: every residue number in this project is a position in the
full UniProt record. UniProt position = challenge position + 24.

1. Record lengths
  [PASS] human UniProt P00533 full length: 1210 aa (expected 1210)
  [PASS] mouse UniProt Q01279 full length: 1210 aa (no fixed expectation asserted)
  [PASS] human challenge construct length: 621 aa (expected 621)
  [PASS] mouse challenge construct length: 623 aa (expected 623)

2. Does the official human construct equal UniProt 25-645?
  [PASS] slice length: 621 aa
  [PASS] character-for-character identity: identical

3. Is UniProt residue 415 a threonine, by three independent routes?
  [PASS] full UniProt record, index 414: T
  [PASS] UniProt 25-645 slice: T
  [PASS] challenge construct, index 390: T

   Caveat: the mouse challenge construct is 623 amino acids, two longer
   than mouse UniProt 25-645 (621 aa), because mouse carries a 2-residue
   insertion near the position matching human 638. The +24 offset holds
   only upstream of that insertion. Domain III (310-480) and our epitope
   (415-466) sit well upstream, so the offset is safe for the work here.
   Do not reuse it for mouse positions past about 638 without rechecking.

4. The eight pH anchors, read in both numbering systems
   (the residue must also be the same in human and mouse, because these
    are the anchors we claim are conserved across the two species)
  [PASS] anchor D416: UniProt=D challenge=D mouse=D | acidic (binder gets HIS)
  [PASS] anchor H418: UniProt=H challenge=H mouse=H | target HIS (binder gets D/E)
  [PASS] anchor E421: UniProt=E challenge=E mouse=E | acidic (binder gets HIS)
  [PASS] anchor E424: UniProt=E challenge=E mouse=E | acidic (binder gets HIS)
  [PASS] anchor H433: UniProt=H challenge=H mouse=H | target HIS (binder gets D/E)
  [PASS] anchor E455: UniProt=E challenge=E mouse=E | acidic (binder gets HIS)
  [PASS] anchor D458: UniProt=D challenge=D mouse=D | acidic (binder gets HIS)
  [PASS] anchor D460: UniProt=D challenge=D mouse=D | acidic (binder gets HIS)

Wrote data/derived/numbering-check.csv

========================================================================
RESULT: ALL CHECKS PASSED.
The numbering convention holds. UniProt pos = challenge pos + 24.
========================================================================
```
