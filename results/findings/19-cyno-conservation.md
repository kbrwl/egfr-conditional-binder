# Cynomolgus monkey EGFR: are the anchors and the sequons conserved?

Computed output of `analysis/19_cyno_conservation.py`. Do not hand-edit.

**This is not a design input and it changes no decision.** The competition page names
human and mouse only, and the decision of 2 October 2026 is to design for human and
mouse only. A pre-launch message had said "mouse and cyno", so whether the face we
aim at survives in cynomolgus monkey (*Macaca fascicularis*, the crab-eating macaque)
was left open rather than answered. This step exists to turn that open question into
a recorded sentence, so that the methods write-up can either make the claim or name
what breaks. No anchor is added, dropped, demoted or promoted because of it.

The answer: 8 of the 8 anchors of the current H370 face hold the same
amino acid in cynomolgus monkey as in human, so none of them loses the charge the
design rule pairs against. The attachment points for sugar chains that sit nearest
those anchors — N352, N361, N413 — are all present in cyno too, including N361, which
is present in human and absent in mouse. The offset between the cyno entry's own
numbering and ours is zero, measured at all 171 positions of domain III rather
than assumed.

EGFR is the epidermal growth factor receptor, the protein we design against. UniProt
is the public archive of protein sequences, and every residue number here is a
position in its full human record P00533.

The cynomolgus monkey sequence is cached at `data/sequences/egfr-macaque-uniprot.fasta`,
with the accession, the address and the date fetched recorded alongside it in
`egfr-macaque-uniprot.source.txt`. Section 1 below says whether that file is tracked
by git, checked by asking git rather than by assuming.

```
========================================================================
CYNOMOLGUS MONKEY EGFR — are the anchors and the sequons conserved?
========================================================================

   THIS IS NOT A DESIGN INPUT AND IT CHANGES NO DECISION. The competition
   page names human and mouse only, and the decision of 2 October 2026 is
   to design for human and mouse only. A pre-launch message had said
   'mouse and cyno', so whether the face we aim at survives in cynomolgus
   monkey was left open rather than answered. This step exists to turn
   that open question into a recorded sentence. Nothing downstream reads
   its output, and no anchor is added, dropped, demoted or promoted
   because of it.

1. The cynomolgus monkey sequence, and where it came from

   UniProt holds no reviewed entry for cynomolgus monkey EGFR: an entry is
   'reviewed' when a curator has read it, and the cyno EGFR entries are
   filled in automatically from a gene prediction instead. The entry used
   here is the better annotated of the two the archive holds: the longer
   of them, the second version of its sequence rather than the first, and
   the one in UniProt's reference proteome for the species. That it is
   unreviewed is the main limit on everything below, which is why a
   reviewed sequence from a closely related monkey is read alongside it in
   section 7.

   cache: data/sequences/egfr-macaque-uniprot.fasta already present (2,152 bytes), not re-fetched

   | field | value |
   |---|---|
   | accession used | A0A2K5WKD8 (unreviewed) |
   | length | 704 aa |
   | address | https://rest.uniprot.org/uniprotkb/A0A2K5WKD8.fasta |
   | date fetched (UTC) | 2026-10-03 (an earlier run) |
   | cached at | data/sequences/egfr-macaque-uniprot.fasta |
   | tracked by git, when this ran | yes, it is tracked by git and travels with the repository |
   | SHA-256 of the cached file | 7d19ee71ad0890ad1681108f808af182add8538948f08c160c609c2f48a3edb6 |

   The header, exactly as UniProt serves it (shown on its own line because
   it contains the character the table above uses as a separator):
     tr|A0A2K5WKD8|A0A2K5WKD8_MACFA Epidermal growth factor receptor OS=Macaca fascicularis OX=9541 PE=4 SV=2

   [PASS] cached sequence file unchanged since it was fetched: fingerprint matches the one recorded at fetch time

2. Lining the two sequences up, and measuring the offset

   Method, the same one step 01 uses for mouse: end-to-end pairwise
   alignment, BLOSUM62 scoring, gap open -11, gap extend -1. Comparing
   position 1 to position 1 straight down two sequences does not work,
   because one extra residue early on shifts everything after it.

   The offset is measured and not assumed. Nothing in a sequence file
   records which numbering convention it follows, and a wrong offset
   shifts every result while raising no error.

   | species | entry | length | positions paired | offset over domain III | share of domain III at that offset |
   |---|---|---|---|---|---|
   | cynomolgus monkey | A0A2K5WKD8 | 704 aa | 704 | ours = theirs +0 | 100.0% of 171 |
   | rhesus macaque | P55245 | 1210 aa | 1210 | ours = theirs +0 | 100.0% of 171 |
   | mouse | Q01279 | 1210 aa | 1208 | ours = theirs +0 | 100.0% of 171 |

   Stated plainly, since this is the number that shifts every later
   result if it is wrong: across domain III our residue number is
   exactly the cyno entry's own residue number, with nothing to add or subtract,
   at every one of the 171 positions (100.0%). Across the whole entry that
   same single offset covers 89.1% of its 704 residues, and the rest is the
   tail discussed next.

   The cyno entry reads straight across, one position after another at
   that single offset, from our position 1 to our position 627. Past
   there the two sequences stop corresponding: this entry is a gene
   prediction, and its last stretch reads the gene differently from the
   human record. Letters the alignment places beyond that point are not
   evidence about the species, and are marked as such below.
   Our extracellular region ends at 645, so its last 18 positions fall in that
   stretch. All 8 anchors sit well inside the readable part, the furthest at
   position 424. The sequon table in section 6 covers the whole
   extracellular region and marks any row that falls beyond the readable
   part, so a reader can see which ones to discount.

   Columns where human has a gap, meaning the cyno entry has extra
   residues there: 0. Columns where the cyno entry has a gap: 506.
   The human record is 506 residues longer than this entry, and that accounts for the gaps exactly.
   The length it is short by is at the far end of the protein, past the
   extracellular region, so it is not our concern.

3. Identity by region

   | Region | Positions | Identical | Identity | Not covered by the entry |
   |---|---|---|---|---|
   | Extracellular region | 25-645 | 600/608 | 98.7% | 13 |
   | Readable part of it | 25-627 | 597/603 | 99.0% | 0 |
   | Domain III | 310-480 | 170/171 | 99.4% | 0 |
   | Candidate epitope | 415-466 | 52/52 | 100.0% | 0 |

   [PASS] domain III identity through the measured map is above 80%: 99.4%. Below this the two sequences are not being compared at matching positions, whatever the alignment says

4. Differences in domain III (310-480), human UniProt numbering
   Notation: S348T means human has S at that position where the other
   species has T.

   human vs cyno:  S348T
   human vs mouse: A313P  S315Y  M318V  V323I  E330D  S348T  N361Y  S364A  R377K  H383R  Q390R  D393E  E412D  R414W  S442G  K467R  (16 of them)

   1 difference between human and cyno across the 171 positions of
   domain III, against 16 between human and mouse. Position 348 differs in both species.

   Inside the fallback epitope 415-466: no differences
   (Mouse has one there, S442G, which is why step 10 rejects any candidate
   touching position 442.)

5. The eight anchors of the current H370 face

   The anchor set is read back from step 08's committed output rather than
   written out here, so there is one copy of it in the project. Six are
   acidic — aspartate or glutamate, negatively charged at both pH values —
   and the binder gets a histidine opposite each. Two are the target
   histidines, H358 and H370, and the binder gets an acidic residue
   opposite those.

   current H370 face

   | anchor | role | cyno | rhesus | mouse | same as human? | what it does to the charge behaviour |
   |---|---|---|---|---|---|---|
   | E344 | acidic | E | E | E | yes | identical |
   | H358 | target histidine | H | H | H | yes | identical |
   | D368 | acidic | D | D | D | yes | identical |
   | H370 | target histidine | H | H | H | yes | identical |
   | E391 | acidic | E | E | E | yes | identical |
   | E400 | acidic | E | E | E | yes | identical |
   | E421 | acidic | E | E | E | yes | identical |
   | E424 | acidic | E | E | E | yes | identical |

   old 415-466 epitope (the recorded fallback)

   | anchor | role | cyno | rhesus | mouse | same as human? | what it does to the charge behaviour |
   |---|---|---|---|---|---|---|
   | D416 | acidic | D | D | D | yes | identical |
   | H418 | target histidine | H | H | H | yes | identical |
   | E421 | acidic | E | E | E | yes | identical |
   | E424 | acidic | E | E | E | yes | identical |
   | H433 | target histidine | H | H | H | yes | identical |
   | E455 | acidic | E | E | E | yes | identical |
   | D458 | acidic | D | D | D | yes | identical |
   | D460 | acidic | D | D | D | yes | identical |

   8 of the 8 current anchors are the same amino acid in cyno as in human.

6. The sequons step 14 found, checked in cyno

   A sequon is the pattern N-X-S/T that a sugar chain attaches to: an
   asparagine, then any amino acid except proline, then a serine or a
   threonine. The rule is step 14's and is called here rather than
   restated, so the two steps cannot disagree about what counts as one.
   A sugar chain is large and mobile, so an attachment point near an
   anchor can cover it some of the time; one present in a species and
   absent in another covers one target and not the other.

   In the status column, 'both' means the pattern is there in human and
   in cyno, 'human only' that cyno has lost it, and 'cyno only' that cyno
   has one where human has none.

   | position | human | cyno, aligned | mouse, aligned (step 14) | cyno status | nearest to a current anchor? |
   |---|---|---|---|---|---|
   | N128 | NKT | NKT | NRT | both | no |
   | N175 | NMS | NMS | NMS | both | no |
   | N196 | NGS | NGS | NGS | both | no |
   | N352 | NAT | NAT | NAT | both | yes |
   | N361 | NCT | NCT | YCT | both | yes |
   | N413 | NRT | NRT | NWT | both | yes |
   | N444 | NIT | NIT | NIT | both | no |
   | N528 | NVS | NVS | NVS | both | no |
   | N568 | NIT | NIT | NIT | both | no |
   | N603 | NNT | NNT | NNT | both | no |
   | N623 | NCT | NCT | NCT | both | no |
   | N639 | NGP | NGS | - | cyno only (past the readable part of the entry) | no |

   The ones that matter here are the attachment points step 14 measured as
   nearest to a current anchor: N352, N361, N413.

     N352: human NAT, cyno NAT, mouse NAT — intact in cyno, intact in mouse
     N361: human NCT, cyno NCT, mouse YCT — intact in cyno, lost in mouse
     N413: human NRT, cyno NRT, mouse NWT — intact in cyno, intact in mouse

   N361 is the one worth a sentence of its own. It is a sequon in human
   and not in mouse, where the asparagine is replaced by a tyrosine, and
   that asymmetry is why step 10 demotes a design leaning on E400, E421
   or E424. In cyno the three residues are NCT, so the sequon is intact:
   cyno matches human here and mouse is the odd one out.

7. Control: the reviewed rhesus macaque entry

   The cyno entry is unreviewed, so a result drawn from it alone rests on
   one automatic gene prediction. Rhesus macaque EGFR is a reviewed entry
   of the full length, and the two monkeys are closely related, so reading
   it alongside says whether anything above is an artefact of that
   prediction. This is a check on the source and not a second answer: the
   question asked is about cyno.

   P55245 (sp|P55245|EGFR_MACMU Epidermal growth factor receptor): 1210 aa, offset ours = theirs +0 across 100.0% of domain III.
   Domain III differences from human: S348T

   [PASS] the two macaque entries agree at every current anchor: they do
   [PASS] the two macaque entries agree about every anchor-adjacent sequon: they do

8. Controls against steps 01 and 14

   This script aligns sequences that step 01 already aligned and finds
   sequons that step 14 already found. Where it recomputes something, it
   is compared against what the earlier step committed, and the run stops
   if the two disagree. Steps 04 and 06 once answered the same question
   two different ways for months before anyone noticed.

   The rule the sequons come from is step 14's own, loaded as a module
   rather than copied, so there is one implementation of it.

   [PASS] cross-check human vs mouse domain III differences, against step 01: 16 residues, identical to 01-domain3-differences.csv
   [PASS] cross-check human sequons in the extracellular region, against step 14: 11 residues, identical to 14-sequons.csv
   [PASS] the mouse substitutions themselves match step 01, not only the positions: all 16 agree

9. The charge rule, tested every time

   Does a substitution keep the charge behaviour the design rule needs?
   Tested on constructed inputs:

     [PASS] E -> E: identical    (no substitution at all)
     [PASS] H -> H: identical    (a target histidine that stays a histidine)
     [PASS] E -> D: charge kept  (acidic stays acidic: glutamate to aspartate)
     [PASS] D -> E: charge kept  (acidic stays acidic: aspartate to glutamate)
     [PASS] E -> Q: charge lost  (glutamate to glutamine: same shape, no charge)
     [PASS] D -> N: charge lost  (aspartate to asparagine: same shape, no charge)
     [PASS] H -> Y: switch lost  (a target histidine replaced, as mouse does at 361)
     [PASS] H -> D: switch lost  (acidic in place of a target histidine is still no switch)

========================================================================
10. What this changes about the design
========================================================================

   Nothing. That is the honest answer and it was the expected one. Cyno is
   not an objective: the three objectives are pH selectivity, mouse
   cross-reactivity and human binding, and the competition page names two
   species. No anchor moves, no ranking term changes, no campaign is
   resized.

   What it records: all 8 anchors of the current H370 face hold the same
   amino acid in cyno as in human, so a binder that works by pairing
   charges against them has the same charges to pair against in cyno.
   The one difference anywhere in domain III, S348T, is not an anchor
   and is the same position that differs in mouse. The attachment
   points nearest the anchors are intact in cyno: N352, N361, N413. That
   includes N361, the sequon human has and mouse does not, which cyno
   has. So the asymmetry behind step 10's demotion of E400, E421 and
   E424 is between human and mouse specifically, rather than between
   human and other species generally. The demotion stands unchanged,
   because mouse cross-reactivity is the objective and cyno is not.
   The methods write-up may say all of this, with the qualifications
   below.

   The qualifications, which belong with the claim wherever it is made:

     - The cyno entry is unreviewed. It is an automatic gene prediction
       that no curator has read. The reviewed rhesus macaque entry agrees
       with it everywhere this step looked, which is reassurance and not
       proof.
     - This is sequence analysis. It says which amino acid sits at each
       position and nothing about whether the face is shaped the same way
       or reachable, which for human is what steps 02 to 06 settled and
       for cyno has not been done at all. No cyno structure was examined.
     - Identical anchors do not make a binder cross-reactive. The
       surrounding surface a binder also touches was not compared position
       by position outside domain III.
     - Nothing was measured about binding. No cyno affinity is predicted,
       claimed or implied.

Wrote data/derived/19-cyno-anchor-conservation.csv, 19-cyno-sequons.csv and
19-cyno-domain3-differences.csv

========================================================================
RESULT: PASSED. 8 of 8 current anchors hold the same amino acid in cyno; 3 of 3 anchor-adjacent
sequons intact in cyno (N352, N361, N413). Changes no decision.
========================================================================
```
