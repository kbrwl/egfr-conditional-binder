# Structural novelty of the submitted designs

Computed output of `analysis/24_structural_novelty.py`. Do not hand-edit.

```
========================================================================
STRUCTURAL NOVELTY — each submitted design against the PDB
========================================================================

343 hits across 7 designs, from:
  foldseek easy-search <binder monomers, PDB format> <PDB database> hits.m8 tmp --alignment-type 1 --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits,alntmscore,qtmscore,ttmscore -e 10 --max-seqs 2000

The structure searched is the binder alone, not the complex. A complex would
match EGFR trivially and tell us nothing about whether our own design
resembles something already known.

1. Best structural match per design, by TM-score over the query

   | design | closest PDB entry | query coverage | TM-score | seq identity | structural verdict |
   |---|---|---|---|---|---|
   | egfr-ph-h370-01 | 3vs8-assembly1_C | 1.00 | 0.747 | 0.9% | moderate |
   | egfr-ph-h370-02 | 7m56-assembly1_A | 0.96 | 0.701 | 1.3% | moderate |
   | egfr-ph-h370-03 | 3vs9-assembly2_B | 1.00 | 0.698 | 1.4% | moderate |
   | egfr-ph-h370-04 | 7oom-assembly1_B | 0.94 | 0.644 | 11.5% | moderate |
   | egfr-ph-h370-05 | 3vs9-assembly2_B | 1.00 | 0.767 (near the 0.8 line) | 0.9% | moderate |
   | egfr-ph-h370-06 | 3vs9-assembly2_B | 1.00 | 0.722 | 0.9% | moderate |
   | egfr-ph-h370-07 | 3vs9-assembly2_B | 1.00 | 0.783 (near the 0.8 line) | 0.9% | moderate |

2. Novelty level

   Sequence half read as at or below 30%: no significant Swiss-Prot homology
   (best e-value 0.30 over 2,239 hits, none below 0.001). The reasoning and
   what turns on it are in this script's header.

   | design | level | why |
   |---|---|---|
   | egfr-ph-h370-01 | **3** | exactly one condition met: moderate structural |
   | egfr-ph-h370-02 | **3** | exactly one condition met: moderate structural |
   | egfr-ph-h370-03 | **3** | exactly one condition met: moderate structural |
   | egfr-ph-h370-04 | **3** | exactly one condition met: moderate structural |
   | egfr-ph-h370-05 | **3** | exactly one condition met: moderate structural |
   | egfr-ph-h370-06 | **3** | exactly one condition met: moderate structural |
   | egfr-ph-h370-07 | **3** | exactly one condition met: moderate structural |

3. What this changes about the submission

   level 3: 7 design(s)

   Every design clears the level 3 gate. None is removed on novelty.

   Reported as near a classification boundary rather than as the label given:
     egfr-ph-h370-05 at TM 0.767. A small change in the alignment would
     move it across, and the platform's own run may land the other side.
     egfr-ph-h370-07 at TM 0.783. A small change in the alignment would
     move it across, and the platform's own run may land the other side.

   This is our reading of a published rule. The platform scores novelty
   automatically on upload and that score is the one that counts.

========================================================================
RESULT: 7 designs scored; 7 at level 3; ALL CLEAR THE GATE.
========================================================================
```
