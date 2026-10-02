# Occupancy at the assay's top concentration

Computed output of `analysis/17_occupancy_table.py`. Do not hand-edit.

K_D, the dissociation constant, is the concentration of the other
molecule at which half the binding sites are occupied; it is a
measure of affinity, with a smaller number meaning a tighter grip.
nM is nanomolar, a billionth of a mole per litre; µM is
micromolar, a thousand times looser for the same number. The
fraction of binder occupied at any concentration follows from the
definition of K_D alone: [analyte] / ([analyte] + K_D). Nothing
here is measured; it is the arithmetic consequence of what K_D
means, computed once rather than typed by hand into
`docs/explainers/07-the-assay-and-what-it-changes.md`,
`docs/decisions-log.md` and `CLAUDE.md`, which all previously
quoted it.

```
========================================================================
OCCUPANCY — how much of the binder is occupied at the assay's top
concentration, for a given K_D
========================================================================

K_D, the dissociation constant, is the concentration of the other
molecule at which half the binding sites are occupied. The fraction
occupied at any concentration follows from that definition alone:
[analyte] / ([analyte] + K_D). Nothing here is measured; it is the
arithmetic consequence of what K_D means, made regenerable rather
than typed into three separate documents by hand.

   Checks on the formula itself, every run:

   [PASS] K_D = 0 (an infinitely tight binder): occupancy = 1.000000, expected 1.0
   [PASS] K_D equal to the analyte concentration (the definition of K_D): occupancy = 0.500000, expected 0.5
   [PASS] occupancy strictly decreases as K_D grows, over 1 to 100000 nM: 0.999 > 0.990 > 0.909 > 0.500 > 0.091 > 0.010

1. Occupancy at the assay's top concentration, 1000 nM

   | K_D | occupied fraction |
   |---|---|
   | 100 nM | 91% |
   | 1 µM | 50% |
   | 3 µM | 25% |
   | 10 µM | 9% |

   A design in the high nanomolar range for human binding at pH 6.5
   sits near full occupancy and gives a strong, easily fitted
   signal. A design deliberately weakened toward the far end of the
   instrument's reportable range gives a progressively smaller
   one, and far enough along it reads as no detectable binding at
   pH 6.5 as well as at pH 7.4, which fails human binding rather
   than demonstrating pH selectivity. Where exactly the instrument
   stops returning a fittable curve is not established here: that
   depends on the size of the design and the density on the chip
   as well as on occupancy, and the organisers have not said.

========================================================================
RESULT: PASSED. Formula checks hold; occupancy tabulated for 4 illustrative K_D values at 1000 nM.
========================================================================
```
