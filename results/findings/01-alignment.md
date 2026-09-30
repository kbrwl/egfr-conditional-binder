# Human vs mouse EGFR alignment (computed)

Computed output of `analysis/01_alignment.py`. Do not hand-edit.

This is a regression test: the expected answer was known before the
script was written, and the script asserts it. See
`alignment-findings.md` in this directory for the interpretation.

```
========================================================================
HUMAN vs MOUSE EGFR ALIGNMENT  (regression test)
========================================================================

human P00533: 1210 aa
mouse Q01279: 1210 aa
method: global pairwise, BLOSUM62, gap open -11, gap extend -1
numbering: full human UniProt positions throughout

1. Identity by region

   | Region | Positions | Identical | Identity |
   |---|---|---|---|
   | Extracellular region | 25-645 | 551/621 | 88.7% |
   | Domain III | 310-480 | 155/171 | 90.6% |
   | Candidate epitope | 415-466 | 51/52 | 98.1% |

   Columns where human has a gap (mouse insertions): 2

2. Domain III differences, human UniProt numbering
   Notation: Q390R means human has Q at 390, mouse has R.

   A313P  S315Y  M318V  V323I  E330D  S348T  N361Y  S364A
   R377K  H383R  Q390R  D393E  E412D  R414W  S442G  K467R

   16 differences total; 14 fall before position 415, 2 at or after.

3. Differences inside the candidate epitope 415-466
   ['S442G']

4. Identical runs of 10+ consecutive residues overlapping domain III

   | Positions | Length | Sequence |
   |---|---|---|
   | 256-312 | 57 | DEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACG |
   | 415-441 | 27 | TDLHAFENLEIIRGRTKQHGQFSLAVV |
   | 443-466 | 24 | LNITSLGLRSLKEISDGDVIISGN |
   | 394-411 | 18 | ILKTVKEITGFLLIQAWP |
   | 331-347 | 17 | GPCRKVCNGIGIGEFKD |
   | 468-483 | 16 | NLCYANTINWKKLFGT |
   | 349-360 | 12 | LSINATNIKHFK |
   | 365-376 | 12 | ISGDLHILPVAF |

5. Regression assertions
  [PASS] domain III difference count: 16 (expected 16)
  [PASS] domain III difference identities: match expected list
  [PASS] epitope 415-466 differences: ['S442G'] (expected ['S442G'])

6. What this means for the design

   The 415-466 block is 52 positions with a single
   difference, S442G. Serine and glycine are both among the smallest
   amino acids, so the local shape barely changes -- this is about as
   close to species-identical as a real surface patch gets.

   Fourteen of the sixteen domain III differences fall before 415.
   That is the whole argument for aiming here rather than elsewhere in
   domain III: a binder confined to this block satisfies mouse
   cross-reactivity by construction.

   LIMIT OF THIS RESULT, stated plainly: this is sequence analysis. It
   says what each residue IS, not which direction it POINTS. Domain III
   folds into a solenoid (spiral-staircase) shape in which residues
   adjacent in sequence can point opposite ways. Whether these anchors
   are reachable is decided by steps 02-06, not here.

Wrote data/derived/01-domain3-differences.csv
Wrote data/derived/01-identical-runs.csv

========================================================================
RESULT: REGRESSION TEST PASSED. Environment reproduces prior work.
========================================================================
```
