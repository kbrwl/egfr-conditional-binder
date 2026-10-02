# Glycan sequons near the H370 anchors

Computed output of `analysis/14_glycan_sequons.py`. Do not hand-edit.

The target in the assay is glycosylated. Step 03 found N444 bonded to a
sugar inside the old 415-466 epitope, and the measurement was never
repeated after step 08 moved the face to the cluster around H370. This
finds every sequon, N-X-S/T, by reading the sequence, compares them with
mouse, and measures each anchor's distance to each attachment nitrogen in
the extended and the tethered structure.

```
========================================================================
GLYCAN SEQUONS — which anchors sit near a place a sugar chain attaches?
========================================================================

The target in the assay is made in human cells and is glycosylated. A
sugar chain is large and mobile, so a binder aimed beside one may find the
site covered in the real molecule while the model shows it open. Step 03
checked this for the old 415-466 epitope and found N444 inside it, bonded
to a sugar. It was never repeated for the current face around H370.

1. The sequon rule, tested every time

   Sequon rule, on constructed sequences:

     [PASS] KNATQ   -> [(2, 'NAT')]   (serine or threonine / a plain sequon)
     [PASS] KNAVQ   -> []   (serine or threonine: valine is not S or T)
     [PASS] KNPTQ   -> []   (proline blocks attachment in the middle position)
     [PASS] KNASQ   -> [(2, 'NAS')]   (serine or threonine / serine accepted)
     [PASS] KNNSTQ  -> [(2, 'NNS'), (3, 'NST')]   (overlapping sequons are both counted)
     [PASS] NAT     -> [(1, 'NAT')]   (a sequon at the very start)
     [PASS] QQNA    -> []   (an asparagine too near the end to form a sequon)

   Mouse comparison, on constructed sequences:

     [PASS] sequon lost in mouse: {4: 'human only', 9: 'both'}
     [PASS] sequon gained in mouse: {9: 'both', 4: 'mouse only'}
     [PASS] identical sequences: {4: 'both', 9: 'both'}

2. Sequons in the human extracellular sequence, residues 25-645

   11 found by reading the sequence. Compared against the mouse sequence
   after aligning the two:

   | position | human | mouse, aligned | status |
   |---|---|---|---|
   | N128 | NKT | NRT | both |
   | N175 | NMS | NMS | both |
   | N196 | NGS | NGS | both |
   | N352 | NAT | NAT | both |
   | N361 | NCT | YCT | human only |
   | N413 | NRT | NWT | both |
   | N444 | NIT | NIT | both |
   | N528 | NVS | NVS | both |
   | N568 | NIT | NIT | both |
   | N603 | NNT | NNT | both |
   | N623 | NCT | NCT | both |

   1 of 11 differ between the species: N361 (human only)
   A sequon present in one species and absent in the other puts a sugar
   chain on one target and not the other. Where one of these lies near
   an anchor it is a cross-reactivity risk at that anchor, and section 4
   says which.

3. The two structures

   6ARU (extended): receptor chain A, 609 residues mapped, 99.7% identity to human EGFR
        numbering offset: ours = file +24 (100.0% of residues), derived and not assumed
        sugar atoms in the file: 164; other protein chains left out of the measurement: B, C
   1NQL (tethered): receptor chain A, 612 residues mapped, 99.8% identity to human EGFR
        numbering offset: ours = file +24 (100.0% of residues), derived and not assumed
        sugar atoms in the file: 176; other protein chains left out of the measurement: B

   Identity checks, by amino acid and not by arithmetic:
   [PASS] identity check, every sequon asparagine, 6ARU: all 11 positions hold the amino acid our numbering says
   [PASS] identity check, every anchor, 6ARU: all 14 positions hold the amino acid our numbering says
   [PASS] identity check, every sequon asparagine, 1NQL: all 11 positions hold the amino acid our numbering says
   [PASS] identity check, every anchor, 1NQL: all 14 positions hold the amino acid our numbering says

   Is a sugar seen bonded to each sequon in the structures? A sugar within
   2.0 A of the attachment nitrogen is a bond. Its absence is weak
   evidence of nothing there, because a structure shows only the sugars that
   hold still.

   | sequon | 6ARU nitrogen to nearest sugar | 1NQL nitrogen to nearest sugar |
   |---|---|---|
   | N128 | 42.71 A | 30.51 A |
   | N175 | 38.11 A | 19.61 A |
   | N196 | 60.52 A | 40.25 A |
   | N352 | 1.44 A — bonded | 1.43 A — bonded |
   | N361 | 1.44 A — bonded | 1.36 A — bonded |
   | N413 | 1.44 A — bonded | 4.97 A |
   | N444 | 1.44 A — bonded | 1.45 A — bonded |
   | N528 | 30.63 A | 1.44 A — bonded |
   | N568 | 45.74 A | 1.44 A — bonded |
   | N603 | 43.53 A | 1.50 A — bonded |
   | N623 | 50.11 A | 13.84 A |

4. Distance from each anchor to every sequon's attachment nitrogen

   Measured from the charged tip of each anchor's side chain, the same
   point step 05 uses, to the ND2 of every sequon. Bands, from step 05:
     under 15 A   likely shadowed at least some of the time
     under 25 A   within reach of an extended chain
     beyond      probably clear

   current H370 face

   | anchor | nearest sequon, 6ARU | nearest sequon, 1NQL | closest overall | assessment |
   |---|---|---|---|---|
   | E344 | N352 at 12.4 A | N352 at 12.0 A | N352 at 12.0 A (1NQL) | likely shadowed |
   | H358 | N352 at 10.8 A | N352 at 10.4 A | N352 at 10.4 A (1NQL) | likely shadowed |
   | D368 | N352 at 15.3 A | N352 at 14.3 A | N352 at 14.3 A (1NQL) | likely shadowed |
   | H370 | N352 at 11.1 A | N352 at 10.6 A | N352 at 10.6 A (1NQL) | likely shadowed |
   | E391 | N352 at 12.0 A | N352 at 11.3 A | N352 at 11.3 A (1NQL) | likely shadowed |
   | E400 | N361 at 17.2 A | N361 at 16.5 A | N361 at 16.5 A (1NQL) | within reach of an extended chain |
   | E421 | N413 at 16.3 A | N361 at 16.2 A | N361 at 16.2 A (1NQL) | within reach of an extended chain |
   | E424 | N361 at 17.6 A | N361 at 16.9 A | N361 at 16.9 A (1NQL) | within reach of an extended chain |

   old 415-466 epitope

   | anchor | nearest sequon, 6ARU | nearest sequon, 1NQL | closest overall | assessment |
   |---|---|---|---|---|
   | D416 | N444 at 11.4 A | N444 at 12.7 A | N444 at 11.4 A (6ARU) | likely shadowed |
   | H418 | N413 at 10.0 A | N413 at 15.1 A | N413 at 10.0 A (6ARU) | likely shadowed |
   | E421 | N413 at 16.3 A | N361 at 16.2 A | N361 at 16.2 A (1NQL) | within reach of an extended chain |
   | E424 | N361 at 17.6 A | N361 at 16.9 A | N361 at 16.9 A (1NQL) | within reach of an extended chain |
   | H433 | N352 at 19.1 A | N352 at 18.6 A | N352 at 18.6 A (1NQL) | within reach of an extended chain |
   | E455 | N361 at 25.8 A | N361 at 25.0 A | N361 at 25.0 A (1NQL) | within reach of an extended chain |
   | D458 | N352 at 27.6 A | N352 at 27.3 A | N352 at 27.3 A (1NQL) | probably clear |
   | D460 | N444 at 27.0 A | N352 at 26.9 A | N352 at 26.9 A (1NQL) | probably clear |

   The current anchors, one line each:

     E344: likely shadowed. Nearest attachment point N352, 12.0 A away in 1NQL; sequon both across species; a sugar is seen bonded there in 6ARU, 1NQL
     H358: likely shadowed. Nearest attachment point N352, 10.4 A away in 1NQL; sequon both across species; a sugar is seen bonded there in 6ARU, 1NQL
     D368: likely shadowed. Nearest attachment point N352, 14.3 A away in 1NQL; sequon both across species; a sugar is seen bonded there in 6ARU, 1NQL
     H370: likely shadowed. Nearest attachment point N352, 10.6 A away in 1NQL; sequon both across species; a sugar is seen bonded there in 6ARU, 1NQL
     E391: likely shadowed. Nearest attachment point N352, 11.3 A away in 1NQL; sequon both across species; a sugar is seen bonded there in 6ARU, 1NQL
     E400: within reach of an extended chain. Nearest attachment point N361, 16.5 A away in 1NQL; sequon human only across species; a sugar is seen bonded there in 6ARU, 1NQL
     E421: within reach of an extended chain. Nearest attachment point N361, 16.2 A away in 1NQL; sequon human only across species; a sugar is seen bonded there in 6ARU, 1NQL
     E424: within reach of an extended chain. Nearest attachment point N361, 16.9 A away in 1NQL; sequon human only across species; a sugar is seen bonded there in 6ARU, 1NQL

5. Controls against steps 03 and 05

   The old 415-466 anchors go through this same code. If it cannot
   reproduce what those steps committed, the answer above is not to be
   believed either.

   N444 in the human sequence: NIT — a sequon
   [PASS] N444 to its nearest sugar atom in 6ARU: 1.44 A here, 1.44 A in step 03
   [PASS] anchors to the N444 CB against step 05's table: all 7 distances agree to 0.1 A

========================================================================
6. What this changes about the design
========================================================================

   Likely shadowed (under 15 A): E344, H358, D368, H370, E391
   Within reach of an extended chain (under 25 A): E400, E421, E424
   Probably clear: none

   Nothing is dropped. These are reported so that step 10's ranking can
   demote a design leaning on a flagged anchor, the way it already demotes
   one leaning on E424 near the cut edge. A crystal structure cannot
   settle whether a mobile chain actually covers a site, so 'within reach'
   means covered some of the time and not blocked.

Wrote data/derived/14-sequons.csv and data/derived/14-anchor-glycan-distance.csv

========================================================================
RESULT: PASSED. 11 sequons in the extracellular region; 5 anchor(s) likely shadowed, 3 within reach of an extended chain, 0 probably clear.
========================================================================
```
