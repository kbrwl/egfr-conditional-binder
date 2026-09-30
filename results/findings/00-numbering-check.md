# Numbering convention check

Computed output of `analysis/00_numbering_check.py`. Do not hand-edit.

All residue numbers in this project are positions in the full UniProt
record (human P00533). The official challenge constructs are the mature
extracellular region, UniProt 25-645, so:

    UniProt position = challenge-construct position + 24

```
========================================================================
NUMBERING CONVENTION CHECK
========================================================================

Convention: all residue numbers in this project are positions in the
FULL UniProt record. UniProt position = challenge position + 24.

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

4. The eight pH anchors, read in both numbering systems
   (identity must also match between human and mouse -- these anchors
    are the ones we claim are cross-species conserved)
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
