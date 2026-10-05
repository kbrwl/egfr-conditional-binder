# The charge-pair filter

Computed output of `analysis/10_charge_pair_filter.py`. Do not hand-edit.

Takes the candidate complexes a binder-design run produces, works out
which target residue each binder contact position faces, and scores each
candidate by how many of those pairs switch on as the surroundings turn
acidic. The design pipeline has no pH term in its objective, so this is
the step that makes the submission pH-conditional rather than generic.

EGFR is the epidermal growth factor receptor, the protein being designed
against. UniProt is the public sequence archive whose numbering this
project uses. mmCIF is the structure file format the design pipeline
writes.

**Ranked by correct pairs first.** The pairs are what produce the pH
switch, and the design pipeline has no interest in pH. It ranks its own
output by `i_pDAE`, a measure of how confident it is in the interface;
here that value is carried through and used only as the last tie-break,
among candidates whose pair terms are all equal, with the more confident
interface first. It never overrides a pair term.

```
========================================================================
CHARGE-PAIR FILTER
========================================================================

Scores design candidates by how many of their contacts form a charge
pair that switches on as the surroundings turn acidic, and rejects the
ones that break a hard rule. Residue numbers are positions in the
human record P00533 in UniProt, the public sequence archive, unless a
line says they are the structure file's own.

Ranked by correct pairs first, because the pairs are what produce the pH
switch. Two demotions never discard a candidate but rank it below an
equivalent one that does not need the same anchor: leaning on E424, near
the trimmed target's cut edge, and leaning on E400, E421 or E424, the
three anchors nearest the human-only N361 sugar-chain attachment point.
The pipeline's own i_pDAE (higher is better) is carried through and
breaks ties among candidates whose pair terms are all equal. Neither
demotion nor i_pDAE overrides a pair term and neither discards a
candidate.

1. The test cases, built by hand and run every time

   No candidate exists yet. These are synthetic complexes written in
   the same mmCIF format the design pipeline emits, each built to
   exercise one branch of the rule. The target residues are real EGFR
   residues at real positions, numbered the way the trimmed target is
   numbered, so the translation back into our numbering is exercised
   too rather than bypassed.

   four-correct-pairs
     The design target: four pairs that all switch the right way. Two acidic target residues faced by binder histidines, and two target histidines faced by binder acidic residues.
     contacts: E344 faced by binder H, D368 faced by binder H, H370 faced by binder E, H358 faced by binder D
     [PASS] correct_pairs: 4
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 0
     [PASS] unresolved_pairs: 0
     [PASS] verdict: meets the pair target

   binder-his-faces-listed-target-his
     A binder histidine facing H370, one of the two target histidines in our anchor set. Both turn positive at pH 6.5 and push apart, which can cancel a correct pair elsewhere on the same face. Rejected rather than scored, even though two correct pairs are present.
     contacts: E344 faced by binder H, D368 faced by binder H, H370 faced by binder H
     [PASS] correct_pairs: 2
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 1
     [PASS] unresolved_pairs: 0
     [PASS] verdict: rejected
     reasons given: binder histidine faces target histidine H370

   binder-his-faces-unlisted-target-his
     The same fault against H418, a target histidine that is not in our anchor set at all. The rejection has to apply to any target histidine, not only the ones we listed, because the physical problem is identical.
     contacts: E344 faced by binder H, H418 faced by binder H
     [PASS] correct_pairs: 1
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 1
     [PASS] unresolved_pairs: 0
     [PASS] verdict: rejected
     reasons given: binder histidine faces target histidine H418

   contacts-442-among-four
     Four correct pairs, and one contact at position 442, the single human/mouse difference inside the original epitope. Until 5 October 2026 that was a rejection; it is now a demotion, so the candidate is kept and meets the pair target and ranks below the four-correct-pairs case that reaches the same count without touching 442. The pair set is identical to that case, so this isolates the 442 demotion from the E424 and N361 ones. Named to sort alphabetically before four-correct-pairs, so the ranking assertion is a real test of the demotion rather than a pass the name-based final tie-break would have given anyway.
     contacts: E344 faced by binder H, D368 faced by binder H, H370 faced by binder E, H358 faced by binder D, S442 faced by binder A
     [PASS] anchor_group_supported_pairs: 4
     [PASS] correct_pairs: 4
     [PASS] edge_reliant_pairs: 0
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 0
     [PASS] n361_reliant_pairs: 0
     [PASS] species_difference_contacts: 1
     [PASS] supported_pairs: 4
     [PASS] unresolved_pairs: 0
     [PASS] verdict: meets the pair target

   one-correct-pair-bare-interface
     One correct pair and nothing else: a design that touches the target at a single position. It breaks no rule, so it is kept and reported. It exists to be ranked against one-correct-pair-real-interface below, which has the same pair count and a real interface around it. Without the interface-size term the two sort only by name, and this one would win it -- which is how a near-non-binder reached the earlier submission file.
     contacts: D368 faced by binder H
     [PASS] correct_pairs: 1
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 0
     [PASS] total_pairs: 1
     [PASS] unresolved_pairs: 0
     [PASS] verdict: below the pair target

   one-correct-pair-real-interface
     The same single correct pair, with four further contacts around it that carry no charge pair. The pair count is identical, so every term above the interface-size one ties, and this design should rank above the bare one because it actually engages the face. Named to sort alphabetically after the bare case, so the assertion fails if the term is removed.
     contacts: D368 faced by binder H, I365 faced by binder A, L369 faced by binder L, I371 faced by binder V, L372 faced by binder A
     [PASS] correct_pairs: 1
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 0
     [PASS] total_pairs: 5
     [PASS] unresolved_pairs: 0
     [PASS] verdict: below the pair target

   no-correct-pairs
     A well-formed interface with no charge pair anywhere in it. Not rejected, because nothing forbidden happens: it is kept, reported, and ranked last. A candidate is discarded for breaking a rule, never for being weak.
     contacts: I365 faced by binder A, L369 faced by binder L, I371 faced by binder V, L372 faced by binder A
     [PASS] correct_pairs: 0
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 0
     [PASS] unresolved_pairs: 0
     [PASS] verdict: below the pair target

   unresolved-side-chain
     Three pairs that would all be correct, but one binder histidine has no side chain in the structure, so there is no evidence its charge reaches anything. It is counted separately and not as correct, which drops the candidate from three correct pairs to two and below the target. If an unresolved side chain were quietly counted, this case would pass.
     contacts: E344 faced by binder H, E391 faced by binder H, D368 faced by binder H (side chain unresolved)
     [PASS] correct_pairs: 2
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 0
     [PASS] unresolved_pairs: 1
     [PASS] verdict: below the pair target

   e424-pair-among-four
     Four correct pairs, one of them on E424, which sits 6.6 angstroms from the cut in the trimmed target. Still meets the pair target, because the pair is correct, but only three of the four count as supported, so it ranks below the four-correct-pairs case that reaches the same count without E424. E424 is also one of the three anchors nearest N361, so it is demoted on two independent grounds at once, not double-counted within either field.
     contacts: E344 faced by binder H, D368 faced by binder H, H370 faced by binder E, E424 faced by binder H
     [PASS] anchor_group_supported_pairs: 3
     [PASS] correct_pairs: 4
     [PASS] edge_reliant_pairs: 1
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 0
     [PASS] n361_reliant_pairs: 1
     [PASS] supported_pairs: 3
     [PASS] unresolved_pairs: 0
     [PASS] verdict: meets the pair target

   e421-pair-among-four
     Four correct pairs, one of them on E421, one of the three anchors nearest N361, a sugar-chain attachment point that exists in human and not in mouse. Still meets the pair target, because the pair is correct, but only three of the four count toward the anchor-group-supported total, so it ranks below the four-correct-pairs case that reaches the same count using only anchors nearest the shared N352 sequon. E421 is not near the trimmed target's cut edge, so this isolates the new demotion from the E424 edge-reliance one above. Named to sort alphabetically before four-correct-pairs, so the ranking assertion below is a real test of the demotion term rather than a pass that would happen anyway from the name-based final tie-break.
     contacts: E344 faced by binder H, D368 faced by binder H, H370 faced by binder E, E421 faced by binder H
     [PASS] anchor_group_supported_pairs: 3
     [PASS] correct_pairs: 4
     [PASS] edge_reliant_pairs: 0
     [PASS] forbidden_contacts: 0
     [PASS] his_his_pairs: 0
     [PASS] n361_reliant_pairs: 1
     [PASS] supported_pairs: 4
     [PASS] unresolved_pairs: 0
     [PASS] verdict: meets the pair target

   renumbered-target
     Four correct pairs, but the target in the returned file is numbered one place off from the structure the run was handed, which is what a design run that cropped and renumbered its target would return. Every pair would still look right, read against the wrong residues, so the candidate is refused rather than scored.
     contacts: E344 faced by binder H, D368 faced by binder H, H370 faced by binder E, H358 faced by binder D
     [PASS] correct_pairs: 0
     [PASS] numbering_status: does not match the input
     [PASS] verdict: not scored
     reasons given: target numbering cannot be reconciled with the input: target numbering does not match the input: no returned chain resembles the input target by number or by sequence

   Ranking: four correct pairs without E424 against four with it
     [PASS] order: four-correct-pairs then e424-pair-among-four

   Ranking: four correct pairs on the N352 side against four with one on the N361 side
     [PASS] order: four-correct-pairs then e421-pair-among-four

   Ranking: four correct pairs clear of 442 against four that touch it
     [PASS] order: four-correct-pairs then contacts-442-among-four

   Ranking: one correct pair with a real interface against one pair alone
     [PASS] order: one-correct-pair-real-interface then one-correct-pair-bare-interface

   Ranking: four correct pairs touching 442 against one clean pair
     [PASS] order: contacts-442-among-four then one-correct-pair-real-interface

   Ranking: three candidates level on every pair term, ordered by i_pDAE
     [PASS] order: tie-z-confident then tie-a-doubtful then tie-m-no-metric

   Ranking: more pairs but a worse i_pDAE against fewer pairs and a better one
     [PASS] order: pairs-more-doubtful then pairs-fewer-confident

   11 cases, all passed.

2. Reading a campaign folder, end to end

   The classification tests above say nothing about whether this
   script can read what the design pipeline actually emits, which is
   the part most likely to be wrong. So a synthetic campaign folder is
   built in the layout the pipeline documents -- accepted complexes and
   binder-only files together in 3_Ranked/, alongside a !_Ranked.csv
   with its real column names -- and read back through the same code
   path a real folder would take.

   [PASS] the two accepted complexes are found and the binder-only files are skipped: 2
   [PASS] contact pairs are read from them: True
   [PASS] the metrics table is joined to every structure: 2
   [PASS] the pipeline's own i_pDAE is carried through: True
   [PASS] our contact set contains the pipeline's: True
   [PASS] the candidate with a binder histidine facing H370 is rejected: rejected
   [PASS] the candidate with four correct pairs meets the target: meets the pair target
   [PASS] ranking puts the four-pair candidate first: True

   Which chain is which was decided by measurement inside that run,
   not by the letter. The pipeline puts the target on chain A and the
   binder on B, the reverse of BindCraft version 1, so a parser that
   trusted the letter would read the wrong molecule and report a full
   set of plausible nonsense.

3. Real candidates

   results/candidates/egfr-r1-short-r1-short/2_Refolded
   10 candidate complex(es), from every mmCIF below the folder (3_Ranked/ not found)

   Metrics table: !_Refolded.csv, 10 rows, 50 columns.
   Every column is carried through to our output unchanged. None
   of them is used to rank anything here.
   Of those, the ones describing interface quality, confidence or
   charge: i_pDAE, i_pTM, pLDDT, pTM, i_pAE, Unbound_Binder_pLDDT, Target_pLDDT, Interface_Residues_detarget, Interface_Residues, Hotspot_Contact_Fraction, Off_Epitope_Contact_Fraction, Interface_BuriedArea, Interface_BuriedArea_Fraction, Surface_Hydrophobicity, Interface_Hydrophobicity, SS_pLDDT, Binder_pI, Binder_Net_Charge

   Which chain is the target and which the binder, by measurement:

   | chain | residues | translate | read as human EGFR | share |
   |---|---|---|---|---|
   | A | 171 | 171 | 171 | 1.00 |
   | B | 48 | 45 | 1 | 0.02 |

   Target is chain A; binder is B.
   Taken from the measurement above rather than from the letter.
   The pipeline puts the target on A and the binder on B, which is
   the reverse of BindCraft version 1, so a parser that assumed a
   letter would read the wrong molecule and report a full set of
   plausible nonsense.

   Target numbering, first candidate against the structure the run
   was handed (6aru_domain3.pdb):
   target numbering identical: all 171 residues present at the input's own numbers, same amino acids.
   Every candidate is checked the same way, and one that cannot be
   reconciled is refused rather than scored.

   Scored 10 candidates, 105 contact pairs.

   Cross-check against the pipeline's own interface list:
     9 of 10 candidates have our 4.5 A
     contact set containing the pipeline's 4.0 A one,
     which is what a looser cutoff should give. A residue the
     pipeline reports that we do not would be a real disagreement
     and stops the run.

3b. The candidates, ranked by correct pairs

   Ties are broken first by how many correct pairs rest on E424, near the
   trimmed target's cut edge, then by how many rest on E400, E421 or E424,
   the three anchors nearest the human-only N361 sugar-chain attachment
   point (docs/explainers/08-sugar-chains-near-the-anchors.md) -- both
   demotions, never exclusions. Remaining ties are broken by how many of
   those pairs also have their charged groups within reach of each other,
   then by how many distinct target positions are paired, then by the
   pipeline's own i_pDAE with the more confident interface first.

   | rank | design | correct | N361-side | of those, in reach | unresolved | neutral | verdict |
   |---|---|---|---|---|---|---|---|
   | 1 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate7_EGFR_domain3 | 2 | 0 | 2 | 0 | 29 | below the pair target |
   | 2 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate1_EGFR_domain3 | 2 | 0 | 1 | 0 | 37 | below the pair target |
   | 3 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate8_EGFR_domain3 | 1 | 0 | 1 | 0 | 18 | below the pair target |
   | 4 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate2_EGFR_domain3 | 1 | 0 | 1 | 0 | 2 | below the pair target |
   | 5 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate5_EGFR_domain3 | 0 | 0 | 0 | 0 | 3 | below the pair target |
   | 6 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate6_EGFR_domain3 | 0 | 0 | 0 | 0 | 3 | below the pair target |
   | 7 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate9_EGFR_domain3 | 0 | 0 | 0 | 0 | 2 | rejected (binder histidine faces target histidine H433) |
   | 8 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate10_EGFR_domain3 | 0 | 0 | 0 | 0 | 1 | rejected (binder histidine faces target histidine H433) |
   | 9 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate3_EGFR_domain3 | 0 | 0 | 0 | 0 | 1 | below the pair target |
   | 10 | egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_candidate4_EGFR_domain3 | 0 | 0 | 0 | 0 | 1 | below the pair target |

   0 candidates reach 3 correct pairs, of which 0 reach 4.
   2 rejected for breaking a hard rule.
   0 not scored, because the target numbering could not be reconciled with the input.
   8 below the pair target, kept and ranked last rather
   than discarded.

   What this changes about the design. The shortlist is the candidates
   that were not rejected, in this order. Choosing among them is a
   judgement to make by hand using the pipeline's own numbers carried
   through alongside. The aim is the widest gap between pH 6.5 and pH 7.4
   with affinity at pH 6.5 as high as the switch allows.

4. Files written

   data/derived/10-candidate-pairs.csv      105 contact pairs
   data/derived/10-candidate-summary.csv    10 candidates
   results/candidates/shortlist.csv         8 not rejected

========================================================================
RESULT: PASSED. 10 candidates: 0 meet the pair target, 2 rejected, 0 not scored (numbering), 8 below the target and kept.
========================================================================
```
