# Trimming the target to domain III

Computed output of `analysis/09_trim_target.py`. Do not hand-edit.

Cuts the EGFR extracellular region down to domain III and writes the
files a binder-design run receives: the trimmed structure, and the
hotspot list -- the positions on the target the run is told to aim at --
expressed in the trimmed file's own numbering rather than ours.

EGFR is the epidermal growth factor receptor, the protein being designed
against. UniProt is the public sequence archive whose numbering this
project uses throughout. The Protein Data Bank, abbreviated PDB, is the
public archive of measured three-dimensional structures.

The cut is made to bring the design run's cost down far enough to afford
a few hundred attempts rather than a few. Whether a fragment this size
holds the shape its anchors sit on is **not** measured here and is
recorded as unverified in `docs/decisions-log.md`. What is measured here
is what the cut breaks: which disulfide bonds it severs, and which
anchors sit near a cut edge.

```
========================================================================
TRIMMING THE TARGET TO DOMAIN III (310-480)
========================================================================

Residue numbers below are positions in the human record P00533 in
UniProt, the public protein sequence archive, unless the line says
otherwise. The structure is 6ARU from the Protein Data
Bank, the public archive of measured three-dimensional structures,
with the antibody already removed by analysis/02.

1. The numbering this script is working in

   Measured by analysis/02 and read back from its committed table,
   rather than assumed here:
     our number = 6ARU residue number +24, holding for 100.0% of the chain
     chain A, 609 residues with coordinates, file numbers 4-612

2. The anchor set, read back from step 08

   Step 08 decided which residues around H370 cluster on one face. It
   is not recomputed here and not retyped here. What this script adds
   is a written expectation of what that set should be, so that if
   step 08's answer changes, this run stops.

   [PASS] anchor set matches step 08's largest H370 cluster: step 08 gives [344, 358, 368, 370, 391, 400, 421, 424]
   8 anchors: E344, H358, D368, H370, E391, E400, E421, E424

3. Do all eight anchors fall inside 310-480?

   If any does not, the cut would remove a residue the design run has
   to aim at. The script stops rather than widening the boundary,
   because a boundary that moves to fit the answer is not a boundary.

     E344: inside
     H358: inside
     D368: inside
     H370: inside
     E391: inside
     E400: inside
     E421: inside
     E424: inside
   [PASS] all anchors inside 310-480: all eight within the boundary

4. What the fragment contains

   171 residues with coordinates, in 1 continuous piece:
     our numbering 310-480  (6ARU numbering 286-456), 171 residues

   The boundary asks for 171 residues and the structure
   resolves every one of them, so the fragment has no gap the
   trimming introduced and none it inherited.

   The fragment is bare protein. The sugar chains and ions in the
   deposited file were already dropped by analysis/02, including the
   one attached at N444, so the design run sees a surface with no
   sugars on it. That is a known simplification rather than a
   measurement, and it applies to the untrimmed target equally.

   For scale: the extracellular region is 621 residues and this fragment is 171, which is 28% of it.

5. Is the mapping right? Checked by residue identity, not arithmetic

   An offset error does not raise an error; it returns a plausible
   residue 24 positions away. The only way to catch it is to ask what
   amino acid each position actually turns out to be, and compare
   that against what it has to be. Three independent routes:

     written    the expectation recorded at the top of this script
     sequence   read out of the human UniProt sequence
     structure  read out of the 6ARU coordinates, via the
                offset analysis/02 measured

   | our pos | written | sequence | structure | file pos | agree |
   |---|---|---|---|---|---|
   | 344 | E | E | E | 320 | yes |
   | 358 | H | H | H | 334 | yes |
   | 368 | D | D | D | 344 | yes |
   | 370 | H | H | H | 346 | yes |
   | 391 | E | E | E | 367 | yes |
   | 400 | E | E | E | 376 | yes |
   | 421 | E | E | E | 397 | yes |
   | 424 | E | E | E | 400 | yes |

   [PASS] all eight anchors read as the expected amino acid by all three routes: every anchor agrees

6. What the cut breaks: disulfide bonds

   A disulfide bond is a covalent link between the sulfur atoms of two
   cysteine residues, written here as the SG atom. It is the strongest
   thing holding a domain's shape together after the backbone, and
   EGFR's extracellular region is held together by a lot of them.

   A cut that leaves one partner inside the fragment and the other
   outside it leaves a cysteine with nothing to bond to. That can make
   the fragment fold differently from the same residues in the intact
   protein, which is the specific failure the unverified entry in the
   decisions log is about. This section lists them. It does not fix
   them: a reported risk is not a repaired one.

   Measured: 24 disulfide bonds in the receptor chain, SG to SG within 2.5 angstroms.

   Checked against the SSBOND records in the deposited file, which
   are the depositors' own statement of which cysteines are bonded:
   [PASS] measured bonds match the file's own SSBOND records: 24 bonds, identical to the file's records

   SEVERED BY THE CUT — one partner inside the fragment, one
   outside it:

   | inside | partner outside | SG-SG distance | partner's domain |
   |---|---|---|---|
   | C470 | C499 | 2.04 A | downstream of domain III (481-645) |

   1 cysteine in the fragment loses the partner it has in the intact protein.
   Each is a place where the fragment could behave differently from
   the intact protein. Whether it actually does is not measured
   here and would be answered by predicting the fragment alone and
   comparing it against the same residues in 6ARU, which is item 11
   in the decisions log's next actions.

   For context: 3 disulfide bonds lie wholly inside the fragment and are unaffected, and 20 lie wholly outside it
   and are removed along with both partners.

7. What the cut breaks: anchors near a cut edge

   The cut creates new chain ends. In the intact protein those
   residues are held by neighbours the fragment no longer contains, so
   they are the part of the model least likely to be right, and an
   anchor beside one of them is the least trustworthy part of the
   result.

   Measured as the closest approach between any non-hydrogen atom of
   the anchor and any non-hydrogen atom of a cut-edge residue. A
   cut-edge residue is one whose sequence neighbour is present in the
   parent structure and absent from the fragment. Gaps the deposited
   structure already had are not cut edges: the cut did not make them.

   Cut edges: [310, 480]

   | anchor | nearest cut edge | distance | within 8 A |
   |---|---|---|---|
   | E344 | 310 | 21.6 A | no |
   | H358 | 310 | 20.9 A | no |
   | D368 | 480 | 16.4 A | no |
   | H370 | 480 | 18.5 A | no |
   | E391 | 480 | 23.9 A | no |
   | E400 | 480 | 12.8 A | no |
   | E421 | 480 | 12.3 A | no |
   | E424 | 480 | 6.6 A | YES |

   1 anchor sits within 8 angstroms of a cut edge: E424
   Treat any design that relies on these as the least supported.
   This does not say they are wrong. It says that if the fragment
   does not hold its shape, these are where it will show first.

8. Writing the fragment

   [PASS] fragment written with the expected residue count: 171 residues in data/structures/6aru_domain3.pdb, expected 171
   The fragment's own numbering range: 286-456, chain A.
   That is 6ARU's numbering, preserved. Add 24 to get ours.

   Anchors read back out of the written file:
     file residue 320 reads E, expected E (our 344) — ok
     file residue 334 reads H, expected H (our 358) — ok
     file residue 344 reads D, expected D (our 368) — ok
     file residue 346 reads H, expected H (our 370) — ok
     file residue 367 reads E, expected E (our 391) — ok
     file residue 376 reads E, expected E (our 400) — ok
     file residue 397 reads E, expected E (our 421) — ok
     file residue 400 reads E, expected E (our 424) — ok
   [PASS] every anchor is present in the written file as the right residue: all eight

   Will the design pipeline accept this file?

   Three things its structure loader does, read out of its own source
   rather than discovered on a paid GPU: it refuses a file containing
   insertion codes, it refuses residues outside the standard twenty,
   and it drops residues missing backbone atoms without stopping.

   [PASS] no insertion codes: none
   [PASS] every residue is one of the standard twenty: 171 standard residues and nothing else
   [PASS] every residue has a complete backbone, so none is dropped: all present

9. The hotspot specification handed to the design run

   Hotspots are the positions on the target the design run is told to
   aim the binder at. They are given in the numbering of the structure
   file supplied alongside them, which here is the trimmed file, which
   carries 6ARU's numbers. Both numberings are written out
   so a reader can check one against the other.

   | role | our numbering | trimmed-file numbering | binder residue |
   |---|---|---|---|
   | E344 — acidic, negative at both pH values | 344 | A320 | histidine |
   | H358 — histidine, neutral at pH 7.4 and positive at pH 6.5 | 358 | A334 | aspartic or glutamic acid |
   | D368 — acidic, negative at both pH values | 368 | A344 | histidine |
   | H370 — histidine, neutral at pH 7.4 and positive at pH 6.5 | 370 | A346 | aspartic or glutamic acid |
   | E391 — acidic, negative at both pH values | 391 | A367 | histidine |
   | E400 — acidic, negative at both pH values | 400 | A376 | histidine |
   | E421 — acidic, negative at both pH values | 421 | A397 | histidine |
   | E424 — acidic, negative at both pH values | 424 | A400 | histidine |

   As one string: A320,A334,A344,A346,A367,A376,A397,A400

   The same eight in our numbering, which is NOT what goes in the
   config file: 344,358,368,370,391,400,421,424

   Wrote design/configs/egfr-domain3-h370.json
   Wrote design/configs/egfr-domain3-h370-notag.json (adds the His tag as an off-target)
   Wrote design/configs/egfr-domain3-h370.md

   WHAT IS CONFIRMED ABOUT THAT CONFIG AND WHAT IS NOT.

   Confirmed. The key names and the hotspot syntax were read out of
   BindCraft2's own repository on 1 October 2026, from its shipped
   example campaign file, its target presets, and its settings
   reference: a campaign file is JSON; targets[] takes name,
   target_path, chains and hotspots; hotspots is one comma-separated
   string in the supplied structure's own numbering, with an optional
   chain-letter prefix. The residue numbers and the chain letter in it
   are computed and checked above by three independent routes.

   Not confirmed. The binder length range and the two design counts are
   defaults rather than decisions, because the molecule category is
   still open and no per-design runtime has been measured. BindCraft2
   also registers a filter metric called Target_Crop_Length, which
   suggests it may crop the target itself; whether that interacts with
   this trim, or renumbers the target in the output, has not been
   checked. It matters because step 10 translates the output's target
   numbering back into ours, and a renumbering would break that
   silently. Both are written down in the companion note.

   The config file itself carries no commentary. BindCraft2 checks a
   campaign file's keys against a fixed catalogue of setting names, so
   a key added for a reader's benefit could make the run refuse the
   file. Everything explanatory is in the markdown file beside it.

   Wrote data/derived/09-hotspot-numbering.csv
   Wrote data/derived/09-severed-disulfides.csv

10. What this changes about the design

   The design run now has a 171-residue target rather
   than a 621-residue one, which
   is what makes a few hundred attempts affordable instead of a few.
   The ratio of run times is not estimated here; it is item 9 in the
   decisions log's next actions, to be measured on the first run by
   timing one trajectory on each.

   1 severed disulfide bond is the price. That is the reason the trimming
   decision stays unverified: nothing here measures whether the
   fragment still holds the shape its anchors sit on.
   1 anchor sits near a cut edge, so any design leaning on
   that position should be ranked below one that does not.

========================================================================
RESULT: PASSED.
Fragment: 171 residues, our 310-480, file numbers 286-456. 1 disulfide bond(s) severed, 1 anchor(s) near a cut edge.
Trimmed target at data/structures/6aru_domain3.pdb; hotspots at design/configs/egfr-domain3-h370.json.
========================================================================
```
