# How far is the signal peptide from our anchors?

Computed output of `analysis/23_signal_peptide_distance.py`. Do not hand-edit.

```
========================================================================
SIGNAL PEPTIDE DISTANCE — how far is residue 25 from our anchors?
========================================================================

The organisers confirmed on 4 October 2026 that the test article starts at
methionine 1, so the signal peptide (residues 1-24) is present rather than
cleaved, and said designs aimed near residue 25 should model it. This step
measures where our own anchors sit relative to that end.

1. The structure: 1nql.pdb, the tethered conformation, chain A
   identified as EGFR by sequence identity, not by the chain letter: 99.8%

   first residue resolved      27 (UniProt numbering)
   mature region begins at     25
   residues 25-26 are not in the file

   The signal peptide itself is resolved in no EGFR structure. The distances
   below are to the first resolved residue, which is where the unmodelled
   leader attaches, and are therefore LOWER BOUNDS on the distance to any
   part of it.

2. Distance from the mature N-terminus to each tight4 anchor

   Bar: 84 A, a 24-residue leader at 3.5 A a residue, fully extended.

   | anchor | distance (A) | reading |
   |---|---|---|
   | E344 (CB (FALLBACK)) | 49.3 | within a fully extended leader's span |
   | H358 (imidazole centroid) | 47.2 | within a fully extended leader's span |
   | D368 (CG) | 58.9 | within a fully extended leader's span |
   | H370 (imidazole centroid) | 61.3 | within a fully extended leader's span |

3. What this changes about the design

   The nearest anchor is 47.2 A from the signal peptide's attachment
   point, inside the 84 A a fully extended 24-residue leader could
   span, so distance alone does not rule out contact with the anchor face.
   It is a weak result in the unfavourable direction: full extension is the
   worst case and a disordered chain rarely approaches it, and nothing here
   measures what the leader actually does. The design does not account for
   it and the write-up says so rather than treating the distance as clearance.

   This is a lower bound measured on one structure, and the leader's own
   extent is unknown. It bounds the question rather than closing it.

========================================================================
RESULT: nearest tight4 anchor is 47.2 A from the mature N-terminus (lower bound).
========================================================================
```
