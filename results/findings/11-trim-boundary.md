# Where the C-terminal trim boundary goes

Computed output of `analysis/11_trim_boundary.py`. Do not hand-edit.

Predicts each candidate target fragment from its sequence alone and
compares it against the same residues in 6ARU, the measured structure,
to decide whether the trim should stop at 480 as step 09 cut it or
extend past the severed disulfide staple at 499.

The predictor is ESMFold via its public interface, which works from a
single sequence with no multiple sequence alignment and is therefore
less accurate than AlphaFold2 in absolute terms. **The absolute numbers
here are not evidence about how well the real fragment folds.** The
differences between the three candidates are the usable result, because
all three pass through the same instrument. This was the available
instrument because the project's graphics-card access is still blocked.

```

========================================================================
WHERE THE C-TERMINAL TRIM BOUNDARY GOES
========================================================================

Residue numbers are positions in the human record P00533 in UniProt,
the public protein sequence archive. Each candidate fragment is
predicted from its own sequence alone and compared against the same
residues in 6ARU, the measured structure from the Protein Data Bank.

RMSD, root-mean-square deviation, is the average distance between two
structures' matching atoms once one has been placed on top of the other
as well as it can be. Smaller means the same shape. It is measured here
between alpha carbons, one backbone atom per residue, so it reports the
fold rather than the side chains.

1. What is being compared, and what the instrument is worth

   The predictor is ESMFold, reached through its public interface. It
   predicts from a single sequence with no multiple sequence alignment,
   which makes it less accurate than AlphaFold2 in absolute terms. It is
   used because it needs no graphics card and no account, and this
   project's graphics-card access is still blocked, so the alternative
   is no measurement at all.

   So the absolute numbers below are NOT evidence about how well the real
   fragment folds. The differences between the three candidates are the
   usable result, because all three pass through the same instrument.

     310-480: our 310-480, 171 residues — as step 09 cut it, severed staple, E424 near the end
     310-480-C470S: our 310-480, 171 residues — the unpaired cysteine mutated to serine, one deliberate difference from the real sequence
     310-499: our 310-499, 190 residues — extended past the staple's far end, whole staple, more margin for E424

2. Does the fragment come back the right shape?

   Two numbers per candidate, because the fragments are different
   lengths and one number would not be a fair comparison.

   'whole fragment' fits each candidate on its own residues, so the
   171-residue candidates are judged on 171 residues and the
   190-residue one on 190. That is each fragment's own quality but it
   is not like-for-like: the longer fragment is being asked to get more
   residues right.

   'common core' fits every candidate on our 310-480 only, which all
   three contain. This is the like-for-like comparison and it is the
   number that decides the question.

   | candidate | whole fragment | common core 310-480 | residues in core |
   |---|---|---|---|
   | 310-480 | 5.00 A | 5.00 A | 171 |
   | 310-480-C470S | 4.97 A | 4.97 A | 171 |
   | 310-499 | 4.90 A | 4.97 A | 171 |

3. How well does it do across the eight anchors specifically?

   The anchors are the eight residues the binder is aimed at, so their
   geometry matters more than the fragment's overall shape. Two
   measurements, which answer different questions:

   'under the core fit' is how far the anchors sit from where 6ARU puts
   them once the whole core has been lined up. It mixes the anchor
   geometry together with any overall drift.

   'anchor face fit' lines the structures up on the eight anchors alone
   and reports how well they can be made to agree. That is the number a
   designed binder depends on: the binder is built against this face, so
   what matters is whether the face has the right shape, which can hold
   even when the rest of the fragment drifts.

   | candidate | anchors under the core fit | anchor face fit |
   |---|---|---|
   | 310-480 | 2.53 A | 0.64 A |
   | 310-480-C470S | 2.50 A | 0.64 A |
   | 310-499 | 2.57 A | 0.65 A |

   Per-anchor distance from 6ARU under the common-core fit:

   | candidate | 344 | 358 | 368 | 370 | 391 | 400 | 421 | 424 |
   |---   |---   |---   |---   |---   |---   |---   |---   |---|
   | 310-480 | 2.4 | 2.8 | 3.1 | 1.6 | 2.4 | 3.3 | 1.6 | 2.4 |
   | 310-480-C470S | 2.4 | 2.8 | 3.1 | 1.6 | 2.3 | 3.3 | 1.6 | 2.3 |
   | 310-499 | 2.2 | 2.7 | 3.1 | 1.7 | 2.3 | 3.5 | 1.8 | 2.7 |

4. What happens around E424 and the C-terminal cut

   E424 is the anchor the question is about: step 09 measured it 6.6 angstroms from
   the cut at 480, and every other anchor at 12.3 angstroms or more.

   Four measurements. E424's own distance from where 6ARU puts it. The
   shift in its charged group, the pair of carboxyl oxygens that a
   charge pair actually depends on, measured under the anchor-face fit
   because that is the frame a binder would be built in. The RMSD of
   E424's neighbourhood, the 39 residues within 10 angstroms of it in
   6ARU, which is the local surface a binder reaching for E424 touches.
   And the last ten residues of each fragment, which is the cut edge
   itself.

   | candidate | E424 | E424 charge group | E424 neighbourhood | own last 10 | shared 471-480 |
   |---|---|---|---|---|---|
   | 310-480 | 2.36 A | 0.89 A | 3.47 A | 4.35 A (471-480) | 4.35 A |
   | 310-480-C470S | 2.32 A | 0.90 A | 3.44 A | 4.45 A (471-480) | 4.45 A |
   | 310-499 | 2.71 A | 0.90 A | 3.74 A | 3.74 A (490-499) | 2.95 A |

   The 'shared 471-480' column is the same ten residues for every
   candidate, so it is the comparable one. For 310-480 those ten
   residues are the fragment's own end; for 310-499 they sit 19
   residues inside it, which is the whole point of extending.

5. The staple, and whether the prediction forms it

   A disulfide bond measures about 2.05 angstroms between the two sulfur
   atoms, written SG. Below 2.5 angstroms counts as formed here.

     310-480: C499 is outside this fragment, so the staple cannot form and C470 is unpaired
     310-480-C470S: C470 is mutated to serine, so there is no unpaired cysteine and no staple to form
     310-499: C470 SG to C499 SG is 2.22 angstroms, so the staple is formed

6. The predictor's own confidence

   pLDDT is ESMFold's confidence at each residue, on a 0 to 100 scale
   after rescaling from the 0 to 1 its interface reports. Above about 70
   is usually read as a reliable backbone and below about 50 as
   essentially no information; those are the field's conventions rather
   than anything measured here. It is the model's opinion of itself, so
   it is reported beside the comparison against 6ARU and not instead of
   it: a prediction can be confident and wrong.

   | candidate | whole fragment | anchors | E424 | own last 10 |
   |---|---|---|---|---|
   | 310-480 | 80.6 | 82.9 | 83.4 | 54.4 |
   | 310-480-C470S | 80.5 | 82.9 | 83.3 | 53.1 |
   | 310-499 | 84.4 | 84.7 | 87.4 | 78.2 |

7. Files written

   data/derived/11-boundary-candidates.csv
   data/derived/11-boundary-anchor-deviations.csv

========================================================================
RESULT: best common-core fit is 310-480-C470S at 4.97 A; best anchor-face fit is 310-480 at 0.64 A.
Whether those differences are large enough to decide the boundary is
judged in results/findings/11-trim-boundary.md, not here.
========================================================================
```
