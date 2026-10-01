# An epitope built around H370

Computed output of `analysis/08_h370_epitope.py`. Do not hand-edit.

H370 is one of two histidines that Liu et al. 2022 (Molecular Therapy
- Oncolytics 27:256-269, doi:10.1016/j.omto.2022.11.001) showed carry EGFR's
pH-dependent antibody binding, found by mutating each histidine to alanine.
The other is H433, which sits in our current epitope and which that paper
built its own design against.

The candidate set is defined from the structure rather than from a sequence
window, because residues adjacent in sequence can point in opposite
directions and a window is therefore not a patch.

```
========================================================================
AN EPITOPE BUILT AROUND H370, COMPARED WITH 415-466
========================================================================

H370 is one of the two histidines Liu et al. 2022 showed experimentally
to carry EGFR's pH-dependent antibody binding, by mutating each histidine
to alanine. The other is H433, which sits in our current epitope and which
that paper built its own design against.

1. Species conservation (question b)

   The conserved sequence run containing H370, as step 01 found it:
   365-376, 12 residues, ISGDLHILPVAF
   Human/mouse differences inside it: none

   The residues immediately flanking that run, which is where it ends and
   why:
     364: human S, mouse A  <- differs
     377: human R, mouse K  <- differs

2. Every residue within 25 A of H370

   Measured between side-chain tips, the parts that form a charge pair.
   'conserved' compares human against the official mouse construct.

   154 residues within reach. Of those, the ones that could carry
   a charge pair — acidic (D or E) or histidine:

   | pos | aa | mouse | conserved | dist from H370 | RSA | exposure | cetuximab |
   |---|---|---|---|---|---|---|---|
   | 370 | H | H | yes | 0.0 A | 0.134 | partial | - |
   | 368 | D | D | yes | 5.5 A | 0.122 | partial | - |
   | 433 | H | H | yes | 9.7 A | 0.648 | exposed | touched |
   | 347 | D | D | yes | 12.4 A | 0.706 | exposed | - |
   | 344 | E | E | yes | 14.2 A | 0.480 | exposed | - |
   | 379 | D | D | yes | 14.4 A | 0.183 | partial | - |
   | 358 | H | H | yes | 16.6 A | 0.500 | exposed | - |
   | 460 | D | D | yes | 16.9 A | 0.186 | partial | - |
   | 458 | D | D | yes | 16.9 A | 0.248 | partial | - |
   | 400 | E | E | yes | 17.3 A | 0.147 | partial | - |
   | 418 | H | H | yes | 18.7 A | 0.032 | buried | - |
   | 391 | E | E | yes | 18.9 A | 0.283 | exposed | - |
   | 393 | D | E | NO | 20.4 A | 0.337 | exposed | - |
   | 383 | H | R | NO | 20.6 A | 0.732 | exposed | - |
   | 412 | E | D | NO | 21.1 A | 0.623 | exposed | - |
   | 455 | E | E | yes | 21.2 A | 0.348 | exposed | - |
   | 388 | D | D | yes | 21.2 A | 0.426 | exposed | - |
   | 416 | D | D | yes | 21.5 A | 0.117 | partial | - |
   | 421 | E | E | yes | 21.9 A | 0.210 | partial | - |
   | 424 | E | E | yes | 22.0 A | 0.429 | exposed | - |

3. Candidate anchors: conserved, reachable, and able to carry a pair

     H370: 0.0 A from H370, RSA 0.134 (partial), target histidine, so the binder gets an acidic residue here
     D368: 5.5 A from H370, RSA 0.122 (partial), acidic, so the binder gets a histidine here
     H433: 9.7 A from H370, RSA 0.648 (exposed), target histidine, so the binder gets an acidic residue here
     D347: 12.4 A from H370, RSA 0.706 (exposed), acidic, so the binder gets a histidine here
     E344: 14.2 A from H370, RSA 0.480 (exposed), acidic, so the binder gets a histidine here
     D379: 14.4 A from H370, RSA 0.183 (partial), acidic, so the binder gets a histidine here
     H358: 16.6 A from H370, RSA 0.500 (exposed), target histidine, so the binder gets an acidic residue here
     D460: 16.9 A from H370, RSA 0.186 (partial), acidic, so the binder gets a histidine here
     D458: 16.9 A from H370, RSA 0.248 (partial), acidic, so the binder gets a histidine here
     E400: 17.3 A from H370, RSA 0.147 (partial), acidic, so the binder gets a histidine here
     E391: 18.9 A from H370, RSA 0.283 (exposed), acidic, so the binder gets a histidine here
     E455: 21.2 A from H370, RSA 0.348 (exposed), acidic, so the binder gets a histidine here
     D388: 21.2 A from H370, RSA 0.426 (exposed), acidic, so the binder gets a histidine here
     D416: 21.5 A from H370, RSA 0.117 (partial), acidic, so the binder gets a histidine here
     E421: 21.9 A from H370, RSA 0.210 (partial), acidic, so the binder gets a histidine here
     E424: 22.0 A from H370, RSA 0.429 (exposed), acidic, so the binder gets a histidine here

   13 acidic, 3 histidine.

4. Cetuximab overlap (question a)

   Candidate anchors inside cetuximab's measured contact set: 1 of 16
     H433 at 3.46 A from the antibody

   For context, the contacts cetuximab makes near this region: [373, 374, 377]
   H370 itself: not touched by cetuximab.

5. Do they cluster within 25 A of each other? (question c)

   Pairwise distances (angstroms):
   | | E344 | D347 | H358 | D368 | H370 | D379 | D388 | E391 | E400 | D416 | E421 | E424 | H433 | E455 | D458 | D460 |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
   | **E344** | — | 12.6 | 6.0 | 13.8 | 14.2 | 19.0 | 21.3 | 16.1 | 16.7 | 27.4 | 23.9 | 22.6 | 23.2 | 25.6 | 24.9 | 28.6 |
   | **D347** | 12.6 | — | 15.2 | 14.4 | 12.4 | 13.1 | 25.7 | 22.3 | 24.8 | 32.3 | 31.6 | 31.2 | 17.7 | 31.0 | 25.8 | 27.4 |
   | **H358** | 6.0 | 15.2 | — | 17.8 | 16.6 | 17.6 | 16.9 | 11.4 | 20.3 | 26.2 | 22.8 | 24.3 | 26.2 | 29.1 | 29.3 | 32.0 |
   | **D368** | 13.8 | 14.4 | 17.8 | — | 5.5 | 19.6 | 24.8 | 21.7 | 13.1 | 22.9 | 22.1 | 19.1 | 11.3 | 17.0 | 12.5 | 14.9 |
   | **H370** | 14.2 | 12.4 | 16.6 | 5.5 | — | 14.4 | 21.2 | 18.9 | 17.3 | 21.5 | 21.9 | 22.0 | 9.7 | 21.2 | 16.9 | 16.9 |
   | **D379** | 19.0 | 13.1 | 17.6 | 19.6 | 14.4 | — | 17.7 | 17.4 | 29.9 | 27.3 | 28.9 | 33.5 | 19.2 | 35.1 | 31.1 | 29.6 |
   | **D388** | 21.3 | 25.7 | 16.9 | 24.8 | 21.2 | 17.7 | — | 5.8 | 27.8 | 17.4 | 17.5 | 27.0 | 28.9 | 33.3 | 34.9 | 34.0 |
   | **E391** | 16.1 | 22.3 | 11.4 | 21.7 | 18.9 | 17.4 | 5.8 | — | 23.8 | 18.7 | 16.8 | 24.1 | 27.8 | 30.5 | 32.3 | 32.7 |
   | **E400** | 16.7 | 24.8 | 20.3 | 13.1 | 17.3 | 29.9 | 27.8 | 23.8 | — | 22.3 | 17.8 | 8.5 | 23.3 | 9.8 | 14.7 | 20.4 |
   | **D416** | 27.4 | 32.3 | 26.2 | 22.9 | 21.5 | 27.3 | 17.4 | 18.7 | 22.3 | — | 7.7 | 17.7 | 26.1 | 22.6 | 26.7 | 24.6 |
   | **E421** | 23.9 | 31.6 | 22.8 | 22.1 | 21.9 | 28.9 | 17.5 | 16.8 | 17.8 | 7.7 | — | 12.4 | 28.4 | 20.2 | 26.3 | 26.6 |
   | **E424** | 22.6 | 31.2 | 24.3 | 19.1 | 22.0 | 33.5 | 27.0 | 24.1 | 8.5 | 17.7 | 12.4 | — | 27.7 | 10.2 | 19.2 | 22.8 |
   | **H433** | 23.2 | 17.7 | 26.2 | 11.3 | 9.7 | 19.2 | 28.9 | 27.8 | 23.3 | 26.1 | 28.4 | 27.7 | — | 23.6 | 15.6 | 12.2 |
   | **E455** | 25.6 | 31.0 | 29.1 | 17.0 | 21.2 | 35.1 | 33.3 | 30.5 | 9.8 | 22.6 | 20.2 | 10.2 | 23.6 | — | 10.8 | 15.6 |
   | **D458** | 24.9 | 25.8 | 29.3 | 12.5 | 16.9 | 31.1 | 34.9 | 32.3 | 14.7 | 26.7 | 26.3 | 19.2 | 15.6 | 10.8 | — | 8.3 |
   | **D460** | 28.6 | 27.4 | 32.0 | 14.9 | 16.9 | 29.6 | 34.0 | 32.7 | 20.4 | 24.6 | 26.6 | 22.8 | 12.2 | 15.6 | 8.3 | — |

   Largest cluster: 8 anchors, E344, H358, D368, H370, E391, E400, E421, E424, maximum internal distance 24.3 A
   Of those, 6 acidic and 2 histidine.
   Includes H370: yes

5b. Are they on one face, or wrapped around the protein?

   The distance test asks whether the anchors fit inside a ball. A binder
   presents something closer to a flat face, so anchors on opposite sides
   of that ball pass the distance test and are still unreachable together.

   Method: take the direction each anchor points away from the protein's
   centre, then measure the angle between those directions. Anchors on one
   face point roughly the same way. A wide spread means the set wraps
   around the protein, however close together the distances look.

   Widest angle within the 8-anchor cluster: 58 degrees, between H358 and E424

   Largest subset that also fits one face, under 90 degrees of spread:
     8 anchors, E344, H358, D368, H370, E391, E400, E421, E424
     span 24.3 A, angular spread 58 degrees
     6 acidic, 2 histidine; includes H370: yes

   The 90 degree limit is a working convention chosen
   here, not a measured property of binders, and it is cruder than
   docking a real backbone against the surface. Treat it as a filter
   that removes clearly unreachable sets rather than as a guarantee
   about the ones that remain.

6. Verdict, and how it compares with 415-466

   The same face test applied to the 415-466 clusters, so both epitopes
   are judged on identical terms:

     D416 E421 E424 E455: span 22.6 A, spread 39 deg — one face
     E424 E455 D458 D460: span 22.8 A, spread 24 deg — one face
     H433 E455 D458 D460: span 23.6 A, spread 30 deg — one face
     D416 E424 E455 D460: span 24.6 A, spread 38 deg — one face
     D416 H418 E421 E424 E455 (conditional on H418): span 24.0 A, spread 44 deg — one face

   415-466, from steps 03 and 05:
     52 residues, one human/mouse difference (S442G)
     7 of 8 anchors reachable, largest cluster 4 anchors within 22.6 A
     contains H433, which cetuximab touches and which Liu et al. used
     two anchors, D458 and D460, sit in the groove against domain IV
     a sugar chain is attached at N444, inside the block

   An epitope around H370, from this script:
     conserved run 365-376, 0 human/mouse difference(s) inside it
     16 candidate anchors conserved and reachable
     largest cluster 8 anchors within 24.3 A
     contains H370: yes
     candidate anchors cetuximab touches: 1

   This epitope is viable: 8 anchors cluster within reach and
   H370 is among them, against a minimum of 3.

   The argument for this epitope, if it is viable: H433 is the histidine
   Liu et al. designed against, so a binder built on it sits on top of
   published work, while H370 is the other histidine their experiment
   implicated and no published binder uses it. Both are experimentally
   supported; only one is taken.

   What this does not settle. Cluster geometry is a necessary condition and
   not a sufficient one, exactly as in step 05. It says nothing about
   whether a foldable binder can present the right partners in the right
   orientations, and nothing about whether the switch clears the assay's
   detection floor.

Wrote data/derived/08-h370-neighbourhood.csv
Wrote data/derived/08-h370-clusters.csv

========================================================================
RESULT: 16 candidate anchors around H370; largest cluster 8 within 24.3 A; includes H370: True.
========================================================================
```
