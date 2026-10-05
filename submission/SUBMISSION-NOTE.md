# Submission note — Challenge 1, EGFR

**Date:** 2026-10-06  
**Replaces:** none — this is the first upload  
**Designs:** 7  
**Distinct backbones:** 1  
**Lengths:** all 48 aa  
**Total billed compute for the whole project:** USD 27.18

## Files uploaded

| file | what it is |
|---|---|
| `challenge1-egfr-submission.csv` | the ranked designs, best first |
| `methods.md` | the written methods, with Appendix A |

## The designs

| rank | name | correct pairs | contacts (4.5 Å) | novelty | flags |
|---|---|---|---|---|---|
| 1 | egfr-ph-h370-01 | 2 | 31 | level 3 | contacts 442 |
| 2 | egfr-ph-h370-02 | 2 | 39 | level 3 | contacts 442 |
| 3 | egfr-ph-h370-03 | 1 | 19 | level 3 | — |
| 4 | egfr-ph-h370-04 | 1 | 3 | level 3 | — |
| 5 | egfr-ph-h370-05 | 0 | 3 | level 3 | — |
| 6 | egfr-ph-h370-06 | 0 | 1 | level 3 | — |
| 7 | egfr-ph-h370-07 | 0 | 1 | level 3 | — |

## What a reader should know before weighing this

- **No design reaches the three correct charge pairs our own rule asks for.** The best reach 2. The rule was set as a floor for a switch that is partial rather than complete at pH 6.5, and this submission does not meet it.
- **All 7 descend from 1 backbone.** Pairwise sequence identity within the set is reported by `analysis/20_submission_check.py`. A set this alike is variants of one idea, not independent attempts.
- **The pH switch is predicted, never measured.** Every pair in Appendix A is an inference from which residues face each other in a predicted structure.
- **Novelty is our own reading of a published scale.** The platform scores it automatically on upload and that score is the one that counts.

## Before uploading

Run `analysis/20_submission_check.py`. It must pass.

