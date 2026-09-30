# Anchor geometry: is the epitope a real surface?

Computed output of `analysis/05_anchor_geometry.py`. Do not hand-edit.

Tests whether the anchors that survived step 03 sit within 25 A
of one another, which is roughly the span a single small binder face can
cover. Distances are measured between side-chain functional groups, not
backbones. Verdict rule fixed in advance: fewer than 3 clustered
anchors means the epitope is dead.

```
========================================================================
ANCHOR GEOMETRY — is the epitope a real surface?
========================================================================

Reach cutoff: 25 A — roughly what one small binder face
              can span. Chosen in advance, not fitted to the answer.
Minimum viable cluster: 3 anchors.

Measuring from side-chain functional groups, not the backbone:
  D -> CG (carboxylate carbon)
  E -> CD (carboxylate carbon)
  H -> centroid of the imidazole ring (CG, ND1, CD2, CE1, NE2)
  unresolved side chain -> CB, reported as a FALLBACK

Anchors carried forward from step 03 (exposed or partial): 7 of 8
  D416, E421, E424, H433, E455, D458, D460
Excluded as buried/unresolved: H418 (buried)

1. Atom used for each anchor

   | anchor | atom measured |
   |---|---|
   | D416 | CG |
   | E421 | CD |
   | E424 | CD |
   | H433 | imidazole centroid |
   | E455 | CD |
   | D458 | CG |
   | D460 | CG |

   No fallbacks — every functional group is fully resolved, so every
   distance below is measured from the atoms that actually form the
   charge pair.

2. Full pairwise distance matrix (angstroms)

   | | D416 | E421 | E424 | H433 | E455 | D458 | D460 |
   |---|---|---|---|---|---|---|---|
   | **D416** | — | 7.7 | 17.7 | 26.1 | 22.6 | 26.7 | 24.6 |
   | **E421** | 7.7 | — | 12.4 | 28.4 | 20.2 | 26.3 | 26.6 |
   | **E424** | 17.7 | 12.4 | — | 27.7 | 10.2 | 19.2 | 22.8 |
   | **H433** | 26.1 | 28.4 | 27.7 | — | 23.6 | 15.6 | 12.2 |
   | **E455** | 22.6 | 20.2 | 10.2 | 23.6 | — | 10.8 | 15.6 |
   | **D458** | 26.7 | 26.3 | 19.2 | 15.6 | 10.8 | — | 8.3 |
   | **D460** | 24.6 | 26.6 | 22.8 | 12.2 | 15.6 | 8.3 | — |

   Closest pair:  D416–E421 at 7.7 A
   Furthest pair: E421–H433 at 28.4 A
   Pairs within 25 A: 15 of 21

3. Largest subset with every pairwise distance under 25 A

   Exhaustive search over all subsets (there are few enough that this is
   exact rather than approximate).

   LARGEST CLUSTER: 4 anchors
   D416, E421, E424, E455
   Maximum internal distance: 22.6 A

   Composition: 4 acidic (each takes a HISTIDINE on the
   binder), 0 target histidine (takes a D or E on the binder).

   NOTE: 4 different subsets of size 4 qualify.
   The design is not forced to one choice of contact set:
     D416, E421, E424, E455  (span 22.6 A)
     E424, E455, D458, D460  (span 22.8 A)
     H433, E455, D458, D460  (span 23.6 A)
     D416, E424, E455, D460  (span 24.6 A)

4. How tight can a cluster be? (relevant to molecule size choice)

   A microbinder (<40 aa) presents a smaller face than a minibinder
   (40-100 aa), so it needs a tighter anchor cluster. Largest cluster
   at several spans:

   For each limit: the largest cluster that fits, and the TIGHTEST
   example of that size (not an arbitrary one, which would hide the
   most compact option available).

   | span limit | largest cluster | tightest such set | its span |
   |---|---|---|---|
   | 12 A | 2 | D416, E421 | 7.7 A |
   | 15 A | 2 | D416, E421 | 7.7 A |
   | 18 A | 3 | H433, D458, D460 | 15.6 A |
   | 20 A | 3 | H433, D458, D460 | 15.6 A |
   | 25 A | 4 | D416, E421, E424, E455 | 22.6 A |

   Reading this for the molecule-size decision: a microbinder (<40 aa)
   presents a small face and needs a tight cluster; a minibinder
   (40-100 aa) can span more. The row where the cluster size drops below
   three is the point at which a binder becomes too small to carry the
   stacked switch at all.

5. Distance from the glycosylation site N444

   Step 03 found a sugar chain covalently attached at N444, inside our
   block. Crystal structures resolve only the innermost sugars, but a
   real glycan extends well beyond and moves. Anchors close to N444 are
   therefore at greater risk of being shadowed in the real molecule.

   A complex N-linked glycan is a branched chain that can sweep
   20-30 A from where it attaches, so the bands below are cautious:
     under 15 A  — likely shadowed at least some of the time
     under 25 A  — within reach of an extended chain
     beyond that — probably clear

   | anchor | distance to N444 CB | assessment |
   |---|---|---|
   | D416 | 11.3 A | LIKELY SHADOWED — read exposure as an upper bound |
   | E421 | 18.2 A | possibly reached by an extended chain |
   | E424 | 27.1 A | probably clear |
   | H433 | 23.0 A | possibly reached by an extended chain |
   | E455 | 29.5 A | probably clear |
   | D458 | 29.7 A | probably clear |
   | D460 | 25.1 A | probably clear |

   Anchors flagged: D416, E421, H433

   This is a flagged risk, not a resolved question, and it cuts both
   ways: a glycan is flexible, so 'within reach' means 'sometimes
   covered', not 'blocked'. It cannot be settled from a crystal
   structure, and we are not going to pretend otherwise. Its practical
   use is as a tie-breaker between otherwise equivalent clusters.

6. VERDICT

   The 415-466 epitope SURVIVES the structure check.

   4 anchors cluster within 25 A (max internal span 22.6 A), against a minimum of 3.
   A single binder face can reach them, so the stacked switch the
   design depends on is geometrically possible.

   4 distinct cluster(s) of size 4 qualify, so the
   design is not forced to one contact set. Comparing them on the
   three things that matter:

   | cluster | span A | target His? | cetuximab overlap | glycan-risk anchors |
   |---|---|---|---|---|
   | D416, E421, E424, E455 | 22.6 | none | none | D416, E421 |
   | E424, E455, D458, D460 | 22.8 | none | none | none |
   | H433, E455, D458, D460 | 23.6 | H433 | H433 | H433 |
   | D416, E424, E455, D460 | 24.6 | none | none | D416 |

   THE METHOD-NOVELTY CLAIM SURVIVES, but it is a CHOICE, not a
   free consequence of the epitope. Specifically:

     H433, E455, D458, D460 (span 23.6 A) contains H433

   Pairing an acidic binder residue against a target histidine is
   what requires reasoning about the TARGET's protonation rather
   than only the binder's, which is the part off-the-shelf
   pipelines do not do. Choosing a cluster without a target
   histidine would give up that claim.

   THE TRADE-OFF, stated plainly. The clusters containing a target
   histidine are the ones that overlap cetuximab's footprint,
   because H433 is the single anchor cetuximab touches. So the
   mechanistically distinctive choice is also the one most open to
   a 'this is cetuximab's epitope' objection. The counter-argument
   is that cetuximab contacts H433 with no pH dependence at all,
   so sharing one residue with it is not sharing a mechanism --
   but that is an argument to make in the write-up, not a
   computed result, and it should not be presented as one.

   For reference, the tightest cluster (D416, E421, E424, E455) has 4 of
   4 anchors outside cetuximab's footprint.

   WHAT THIS DOES NOT ESTABLISH. Geometric reachability is a
   necessary condition, not a sufficient one. It does not show that a
   foldable binder exists that presents the right partner residues in
   the right orientations, nor that the resulting switch is large
   enough to clear the assay's detection floor at pH 7.4. Those are
   design and prediction questions, still open.

Wrote data/derived/05-anchor-distance-matrix.csv
Wrote data/derived/05-anchor-clusters.csv
Wrote explorer/anchor-viewer.html  (open in a browser)
Wrote explorer/anchor-view.pml     (PyMOL: run this file)

========================================================================
RESULT: largest cluster = 4 anchors within 25 A. Epitope SURVIVES.
========================================================================
```
