# Anchor geometry: is the epitope a real surface?

Computed output of `analysis/05_anchor_geometry.py`. Do not hand-edit.

Tests whether the anchors that survived step 03 sit within 25
angstroms of one another, which is roughly how far across a single small
binder can reach with one face. Distances are measured between the tips of
the side chains, the parts that actually form a charge pair, rather than
between backbone atoms, which are the same chemistry in every residue.
The rule was fixed before measuring: if fewer than 3 anchors
cluster, this epitope fails and the fallback stretches are tried instead.

```
========================================================================
ANCHOR GEOMETRY — is the epitope a real surface?
========================================================================

Reach cutoff: 25 A — roughly what one small binder face
              can span. Fixed before the measurement, so it cannot have
              been tuned to suit the answer.
Minimum viable cluster: 3 anchors.

Measured from the charged tip of each side chain, because that is what
forms a charge pair. Backbone-to-backbone distances were the alternative
and were rejected: the backbone is identical chemistry in every residue
and sits several angstroms from the charge.
  D (aspartate) -> CG, the carboxylate carbon
  E (glutamate) -> CD, the carboxylate carbon, one atom further out
  H (histidine) -> centre of the imidazole ring (CG, ND1, CD2, CE1, NE2)
  unresolved side chain -> CB, the first side-chain carbon, flagged as a
                          fallback

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

   Every anchor against every other anchor. Each cell is the distance
   between those two functional groups, so the table shows which anchors
   could share one binder face and which are on opposite sides.

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

   Every possible subset of the anchors is checked, so this answer is
   exact. With this few anchors that is quick, which is why no
   approximate clustering method is used.

   Largest cluster: 4 anchors
   D416, E421, E424, E455
   Maximum internal distance: 22.6 A

   Composition: 4 acidic (each one calls for a histidine on
   the binder), 0 target histidine (calls for a D or E on the
   binder). That is the shopping list the binder has to present.

   4 different subsets of size 4 qualify, so the
   design has a choice of contact set:
     D416, E421, E424, E455  (span 22.6 A)
     E424, E455, D458, D460  (span 22.8 A)
     H433, E455, D458, D460  (span 23.6 A)
     D416, E424, E455, D460  (span 24.6 A)

4. How tight can a cluster be? (relevant to molecule size choice)

   A microbinder (under 40 amino acids, abbreviated aa) presents a smaller
   face than a minibinder (40-100 aa), so it needs a tighter anchor
   cluster. Largest cluster at several spans:

   For each limit: the largest cluster that fits, and the most compact
   example of that size. An arbitrary example of that size would hide how
   tight the best available option actually is.

   | span limit | largest cluster | tightest such set | its span |
   |---|---|---|---|
   | 12 A | 2 | D416, E421 | 7.7 A |
   | 15 A | 2 | D416, E421 | 7.7 A |
   | 18 A | 3 | H433, D458, D460 | 15.6 A |
   | 20 A | 3 | H433, D458, D460 | 15.6 A |
   | 25 A | 4 | D416, E421, E424, E455 | 22.6 A |

   What this decides: how big the molecule has to be. A microbinder
   presents a small face and needs a tight cluster; a minibinder (40-100
   aa) can span more. The row where the cluster size drops below three is
   the size at which a binder becomes too small to carry the stacked
   switch at all, so pick a size above that row.

5. Distance from the sugar attachment point N444

   Step 03 found a sugar chain attached at N444, inside our block.

   What is measured here, and how it differs from step 03. The two steps
   measure different things, and read together without this note they
   look like they disagree:

     step 03 measures the distance from each anchor to the sugar atoms
     that are actually present in the structure file, with a 5 A cutoff.
     It found no anchor within 5 A, so no anchor touches the sugars we
     can see.

     this step measures the distance from each anchor to the point where
     the chain is attached, N444, with the 15 A and 25 A bands below.
     The chain itself is longer than the part the structure shows.

   Both are correct. Step 03 answers 'does an anchor touch a sugar atom we
   have coordinates for', and the answer is no. This step answers 'could
   the full chain reach an anchor', and the answer is that three of them
   are close enough that it might. The second question is worth asking
   because a structure only shows the first few sugars of a chain that
   continues past them, so step 03's answer does not settle it.

   A complex N-linked glycan is a branched chain of sugars that can
   sweep 20-30 A from where it attaches, so the bands below are
   deliberately cautious:
     under 15 A  — likely shadowed at least some of the time
     under 25 A  — within reach of an extended chain
     beyond that — probably clear

   | anchor | distance to N444 CB | assessment |
   |---|---|---|
   | D416 | 11.3 A | likely shadowed — read exposure as an upper bound |
   | E421 | 18.2 A | possibly reached by an extended chain |
   | E424 | 27.1 A | probably clear |
   | H433 | 23.0 A | possibly reached by an extended chain |
   | E455 | 29.5 A | probably clear |
   | D458 | 29.7 A | probably clear |
   | D460 | 25.1 A | probably clear |

   'Flagged' below means flagged by the attachment-point measure used
   here, at 15 A or 25 A. It does not contradict step 03, which found
   no anchor within 5 A of a sugar atom present in the file.
   Anchors flagged: D416, E421, H433

   This is a flagged risk rather than a settled question, and it cuts
   both ways: a glycan is flexible, so 'within reach' means 'covered
   some of the time' rather than 'blocked'. A crystal structure cannot
   settle it either way. In practice we use it as a tie-breaker between
   clusters that are otherwise equally good.

6. Verdict

   The 415-466 epitope survives the structure check.

   4 anchors cluster within 25 A (max internal span 22.6 A), against a minimum of 3.
   A single binder face can reach them, so the stacked switch the
   design depends on is geometrically possible.

   4 distinct cluster(s) of size 4 qualify, so the
   design has a choice of contact set. Comparing them on the three
   things that matter:

   The last column uses the attachment-point measure from section 5,
   at 15 A and 25 A. It is not step 03's 5 A check against the sugar
   atoms present in the file, which found no anchor within 5 A.

   | cluster | span A | target His? | cetuximab overlap | within 25 A of N444 |
   |---|---|---|---|---|
   | D416, E421, E424, E455 | 22.6 | none | none | D416, E421 |
   | E424, E455, D458, D460 | 22.8 | none | none | none |
   | H433, E455, D458, D460 | 23.6 | H433 | H433 | H433 |
   | D416, E424, E455, D460 | 24.6 | none | none | D416 |

   Clusters that include a target histidine, which is the half of
   the pairing rule that needs an acidic residue on the binder:

     H433, E455, D458, D460 (span 23.6 A) contains H433

   This arrangement is published, so it is not ours to claim. Liu X
   et al. (2022), Molecular Therapy - Oncolytics 27:256-269,
   doi:10.1016/j.omto.2022.11.001, mapped EGFR's own H370 and H433
   as the determinants of pH-dependent antibody binding and then
   put an acidic residue against H433 deliberately, on this target
   and at pH 6.5 against 7.4. See the decisions log under Resolved.

   That paper does support the rejection criterion this project
   derived from first principles. Their histidine-on-the-binder
   variant facing H433 changed neither pH-dependency nor affinity,
   while the acidic variants improved pH-dependency substantially.

   On tooling, what is checkable is narrower than what this file
   used to assert. Verified 1 October 2026 from the repositories:
   AlphaFold2 and AlphaFold-Multimer take sequences and database
   paths only; RFdiffusion's inference configuration conditions on
   backbone, contigs, hotspot residues and a closed set of
   radius-of-gyration and contact potentials; ProteinMPNN's
   arguments contain no pH, pKa or protonation option. BindCraft2
   is the exception worth stating rather than hiding: its
   filters.py carries a side-chain pKa table and reports Binder_pI
   and Binder_Net_Charge at a hard-coded REPORTED_PH of 7.4. That
   is the binder's own sequence charge, the pH is not a user
   setting in its 235-setting catalogue, and both are reported
   readouts rather than design objectives.

   The trade-off on cetuximab overlap still stands. The clusters
   containing a target histidine are the ones that overlap
   cetuximab's footprint, because H433 is the single anchor
   cetuximab touches. Step 07 confirmed that overlap holds for
   unmodified cetuximab and not merely for the engineered variant
   in 6ARU. The counter-argument, that cetuximab contacts H433 with
   no pH dependence so sharing a residue is not sharing a
   mechanism, is something to argue in the write-up; nothing here
   computes it and it should not be presented as a result.

   For reference, the tightest cluster (D416, E421, E424, E455) has 4 of
   4 anchors outside cetuximab's footprint.

   What this does not establish. Geometric reachability is necessary
   but not sufficient: it does not show that a foldable binder exists
   which presents the right partner residues in the right
   orientations, nor that the resulting switch is large enough to clear
   the assay's detection floor at pH 7.4. Those are design and
   prediction questions, and both are still open.

7. Second scenario: the same question with H418 included

   Everything above excludes H418, following step 03, which measured it
   as buried in 6ARU. Step 06 measured it as partially exposed in 1NQL,
   and domain III has the same fold in both structures, so the two
   measurements disagree and the question is open. Because 6ARU has the
   antibody clamped on, the burial there may be the antibody holding that
   side chain rather than a property of the receptor alone.

   This section answers what the clusters would be if H418 is usable. It
   is conditional on that unresolved question and is labelled as such
   wherever the numbers are used.

   Anchors against domain IV, read from step 06's table: D458, D460

   Anchors considered: D416, H418, E421, E424, H433, E455, D458, D460

   Distances from H418 to the others:
     H418 to D416: 6.4 A
     H418 to E421: 7.5 A
     H418 to E424: 18.0 A
     H418 to H433: 25.2 A
     H418 to E455: 24.0 A
     H418 to D458: 27.4 A
     H418 to D460: 26.2 A

   H418 to H433 is 25.2 A, against a reach cutoff of
   25 A. The two target histidines cannot both be
   reached by one binder, so a design uses one or the other.

   Largest cluster with H418 available: 5 anchors
   D416, H418, E421, E424, E455
   Maximum internal distance 24.0 A

   Compared with 4 anchors at 22.6 A when H418
   is excluded.

   All clusters of that size, with the same columns as section 6:

   | cluster | span A | target His? | cetuximab overlap | within 25 A of N444 | domain IV groove |
   |---|---|---|---|---|---|
   | D416, H418, E421, E424, E455 | 24.0 | H418 | none | D416, E421 | none |

   Clusters here that carry a target histidine, avoid the
   antibody footprint and avoid the domain IV groove:
     D416, H418, E421, E424, E455 (span 24.0 A)

   No cluster in section 6, where H418 is excluded, manages all
   three at once. That is what makes settling H418 worth doing
   before choosing a contact set.

   These numbers hold only if H418 is usable. Until that is settled they
   describe an option, not a decision.

Wrote data/derived/05-anchor-distance-matrix-with-h418.csv
Wrote data/derived/05-anchor-clusters-with-h418.csv
Wrote data/derived/05-anchor-distance-matrix.csv
Wrote data/derived/05-anchor-clusters.csv
Wrote explorer/anchor-viewer.html  (open in a browser)
Wrote explorer/anchor-view.pml     (PyMOL: run this file)

========================================================================
RESULT: largest cluster = 4 anchors within 25 A. Epitope survives.
========================================================================
```
