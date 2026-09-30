# Decisions log

Read this first in any new chat. Update it when something moves between sections.

Last updated: 30 September 2026 (numbering convention, structure-check results)

---

## The task

Anthropic x Adaptyv Protein Design Competition, Challenge 1 (of 5).
Official page: proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
Deadline **4 October 2026, 23:59 AoE**.

Four further challenges follow weekly until **1 November 2026**. Everything built
here — the analysis scripts, the numbering guard, the findings format — should be
reusable against a new target by swapping the sequences and the epitope, not
rewritten. Treat single-use code as a defect.

Other official links, so they never have to be re-fetched:
novelty rules adaptyvbio.com/blog/novelty · Slack join.slack.com/t/proteinbase
(invite `zt-3evw8fs9z-tU9ItWVvw4ySctUuPvIhLQ`) · models
github.com/anthropics/uplifting-biomolecular-modeling. Full copy in
`rules-reference.md`.

Design binders that:

1. bind the extracellular region of human EGFR (residues 25–645, 621 aa)
2. also bind mouse EGFR
3. bind human EGFR at pH 6.5 with **no detectable binding** at pH 7.4

Ranking priority, official: pH selectivity first, mouse cross-reactivity second,
affinity third. Organisers state explicitly that a weak but clearly pH-sensitive
binder may outrank a high-affinity binder that is not pH-sensitive.

**Submission format.** CSV ranked best-first, columns `name`, `sequence`,
`molecule_class` (one of protein / nanobody / scfv / fab_kappa / fab_lambda).
Fabs submitted as `{VH}:{VL}`. Length 10–250 amino acids.
Track 1: 20–40 designs. Tracks 2 and 3: up to 20.

**Molecule categories** (winners announced across these, so pick deliberately):
microbinder <40 aa, minibinder 40–100 aa, large binder >100 aa, nanobody,
antibody.

**Hard rules.** De novo *and* zero-shot — may not start from an existing binder
and modify it. Must have adequate sequence and structural diversity from known
proteins. Embedded instructions or prompt injection in a submission is grounds
for disqualification.

**Selection.** ~1500 designs screened per challenge, split 50% Track 1,
25% Track 2, 25% Track 3. Track 1 takes each participant's own top 20.
Tracks 2 and 3 pool all designs and a Claude workflow selects ~375 per track
based on predicted design quality, design novelty **and method novelty**,
reading the documentation participants submit alongside their sequences.

**Resources.** Claude and Modal credits go to selected Track 1 and 2 teams only;
Track 3 is self-supported. Freely available to everyone regardless of track:
inference-optimised folding and design models at
github.com/anthropics/uplifting-biomolecular-modeling with a simplified
binder-design prompt in the accompanying technical report, and BindCraft2,
released September 2026 for both academic and industry use.

Deliverable is sequences plus written methods. No lab work by us.

---

## Settled

**Numbering convention: full UniProt P00533, always.**
Every residue number in this project is a position in the full human UniProt
record (1210 aa), where residues 1–24 are the signal peptide and the mature
extracellular region runs 25–645. The official challenge constructs are the
mature region, so their position 1 is our position 25:

    UniProt position = challenge-construct position + 24

Nothing in a FASTA file or a PDB file records which convention it uses, so a
mix-up shifts every number by 24 with no error raised — the results just become
wrong. Structure files are worse: PDB entries frequently number by the mature
protein. **Any structure-derived numbering must be offset-checked empirically
before use, never assumed.** Guarded by `analysis/00_numbering_check.py`, which
asserts the convention and exits non-zero if it breaks. Verified 30 September
2026: human construct is character-for-character UniProt 25–645, residue 415 is
T by three independent routes, all eight anchors read correctly in both systems.

One caveat the check surfaced: the mouse construct is 623 aa, two longer than
mouse UniProt 25–645, because of a 2-residue insertion near human-equivalent
position 638. The +24 offset is valid only upstream of that. Domain III sits far
upstream, so it is safe here — but not for mouse positions past ~638.

**Target region: domain III of EGFR, residues ~310–480.**
This is where cetuximab and panitumumab bind, and what the organisers recommend.

**Candidate epitope: residues 415–466.**
52 consecutive positions, one conservative difference (S442G). Six acidic anchor
residues for histidine pairing. Details in `alignment-findings.md`.
Status: leading candidate, **conditional on the structure check passing**.

**pH mechanism: charge pairing, decided per position.**
Histidine is the only amino acid that switches charge between pH 7.4 and 6.5
(neutral above, positive below). An interface is a patterned surface and each
pair behaves independently, so the rule is positional:

- opposite an acidic target residue (D/E) → put **histidine** on the binder
- opposite a target histidine (H418, H433) → put **D or E** on the binder

Both switch on as pH drops. Do **not** place histidine opposite H418 or H433 —
that pair switches off in the tumour and can cancel a correct pair elsewhere.

Two known limits: the switch is partial rather than binary at pH 6.5, so use
three or four pairs, not one; and a histidine's flipping point shifts with its
neighbours, which is why no tool predicts pH selectivity reliably.

**Track: Track 3 (open track).** Not a choice — applications for Tracks 1 and 2
closed 24 September 2026 and we found the competition on 30 September.
Track 3 is self-supported: our own tools, our own compute, no Anthropic credits
or Claude/Modal subscriptions. Up to 20 designs.

**Testing is not guaranteed.** Track 1 teams (~30, selected) hold a reserved
allocation of ~15 designs per challenge. Tracks 2 and 3 have **no guaranteed
allocation**. Roughly 1500 designs are screened per challenge, split 50% Track 1
/ 25% Track 2 / 25% Track 3, so Track 3 competes for ~375 pooled slots. Our
designs may never be made. Plan accordingly: the write-up has standalone value
even if nothing is screened.

**No cash prize.** Non-cash prizes, details unannounced.

**UNVERIFIED — whether winners are judged within or across tracks.** Screening
quotas are per-track (50/25/25). The winner categories described by the
organisers (highest affinity, most cross-reactive, most pH-sensitive) make no
mention of tracks. Best reading: Track 3's disadvantage is getting screened, not
the judging afterwards, meaning a tested design competes with professional lab
submissions. Not stated anywhere. **Ask in the Proteinbase Slack channel.**

**Strategy: compete where no established method exists.**
Affinity optimisation is a solved-ish problem that well-resourced labs will win.
pH selectivity has no validated predictive tool — nobody can compute it, everyone
is reasoning from first principles. Cross-species is a sequence-analysis problem
with a definite answer. Concentrate effort on **pH selectivity and mouse
cross-reactivity**, which are also the two the organisers rank highest. The
organisers confirm this reading explicitly.

(This sentence previously said "objectives 2 and 3". It was written when the
source numbered them affinity first, and it meant mouse and pH. Numbered the way
this repo numbers them — pH, mouse, affinity — it read as mouse and affinity,
which is the opposite of the intent. Objectives are named rather than numbered
throughout for that reason.)

**Target-side histidines are the likely method differentiator.** Off-the-shelf
pipelines add histidines to the binder and score the result. Using the target's
own H418 and H433 as anchors for acidic residues on our side requires reasoning
about the target's protonation, which no standard pipeline does. This is worth
pressing on and worth writing up in detail.

**Strategy: the methods write-up is part of the submission, not paperwork.**
In Track 3 we hold no reserved testing slot, so the submission competes in a
pool. Selection runs all submissions through a Claude workflow weighing
predicted design quality, design novelty and method novelty, reading the
documentation submitted alongside the sequences. Organisers state selection
will not rely on a single in silico metric. A clear written argument for the
epitope choice and the positional charge-pairing rule is therefore the main
available edge, not paperwork. Budget real time for it.

**Sequence provenance verified, both species.** The human sequence (UniProt
P00533, residues 25-645, 621 aa) is character-for-character identical to the
official challenge page. The official mouse construct is 623 aa, two longer than
UniProt Q01279 residues 25-645; the two extra residues sit at human-equivalent
position 638, far downstream of domain III. Aligned against the official
constructs, domain III still shows the same 16 differences, 415-466 still shows
only S442G, and all eight anchors (D416, H418, E421, E424, H433, E455, D458,
D460) are identical in both. Checked 30 September 2026.

**Affinity target is deliberately marginal, not maximal.** The pH-selectivity
requirement demands
*no detectable binding* at pH 7.4 — a threshold, not a ratio. A design binding
at 10 nM / 200 nM (20-fold selective) fails, because 200 nM is plainly
detectable. A design binding at 2 uM / undetectable passes. Organisers confirm a
weak but clearly pH-sensitive binder may outrank a strong non-selective one.
Standard design pipelines maximise affinity by default, which produces exactly
the failing case. **Do not accept the pipeline's built-in objective.** Target a
baseline weak enough that the pH 7.4 state falls under the assay detection floor,
then build the largest switch possible on top of it.

The organisers do not say which assay is used, so **we do not know where the
detection floor actually sits.** The threshold we are designing to is a number we
have not been given. That is why we aim for a wide margin rather than a computed
one, and why the write-up should say the margin is a judgement.

**THE STRUCTURE CHECK PASSED. The 415–466 epitope is real.** Computed
30 September 2026 by `analysis/02`–`06`, all output in `results/findings/`.

*Numbering.* Both 6ARU and 1NQL number by the mature protein: UniProt = PDB + 24,
derived empirically and consistent across 100% of aligned residues in both. This
was the trap the brief warned about, and it was real. Chain A is the receptor in
6ARU, as the competition page said — confirmed by measurement (99.7% identity)
rather than trusted. Chains B and C are the cetuximab light and heavy chains per
the file's own annotation. The whole epitope is resolved in both structures, with
no missing loops, so every position has an answer.

*Exposure (step 03, on the receptor alone with the Fab removed).* **Seven of
eight anchors are exposed or partially exposed.** H418 is buried in 6ARU
(RSA 0.032) — but see the Open section, because 1NQL disagrees. D458 is
classification-sensitive: it lands either side of the boundary depending on which
published normalisation table is used, so its label is not robust. E421 and D458
sit within 0.05 of a cutoff and should be read as near-boundary, not as their
label.

*Geometry (step 05).* **Four anchors cluster within 25 Å** (max span 22.6 Å),
against a minimum of three. Distances measured between side-chain functional
groups — carboxylate carbons and the imidazole ring centroid — not backbones, and
no anchor needed a CB fallback. Four distinct 4-anchor clusters qualify, so the
design is not forced to one contact set. The tightest viable triad is
H433/D458/D460 at 15.6 Å, which keeps a microbinder in play; below ~18 Å only two
anchors fit, which is too few for a stacked switch.

*Cetuximab's real footprint (step 04).* 24 EGFR residues contact the Fab. It
touches 10 of our 52 residues (19.2%) and exactly **one** of our eight anchors,
H433. So we sit adjacent to and partly overlapping the druggable surface without
reproducing it — the position we wanted for the novelty requirement.

*Conformation (step 06).* The epitope is **accessible in both published
conformations**, and slightly *more* accessible in the closed one (mean RSA 0.173
in 6ARU, 0.189 in 1NQL). The comparison is not vacuous: with domain III
superposed (RMSD 1.08 Å, same fold), domains I–II sit 23.4 Å apart between the
two structures, so these genuinely are different global arrangements. EGF contacts
none of our block.

*The second half of the epitope rests against domain IV (step 06, corrected
1 October 2026).* 19 residues in 6ARU and 18 in 1NQL, spanning roughly 446–466,
are within 4.5 Å of residues outside domain III, and the partners are in the
481–524 range, which is domain IV. The two structures give nearly the same list at
nearly the same distances, so this is a standing feature of how the protein folds
rather than something the closed shape introduces — the original concern about the
tether is still answered. Two anchors, **D458 and D460**, are in that group.

This does not change their accessibility numbers: step 03 measured those on the
whole receptor chain with domain IV already present, so its effect is included.
What it adds is that those two anchors sit in a groove between two domains instead
of on an open face, which is a harder shape for a designed binder to fit and makes
those contacts more sensitive to any shift in how the two domains sit together.
Weighed in the cluster comparison under Open.

An earlier version of this entry claimed only one residue (461) was contacted from
outside domain III. That came from a defect in step 06, which looked residue
numbers up in the wrong direction and reported positions 48 away from the real
ones. Corrected and recorded under Resolved.

**New design constraint: avoid contacting position 442.** S442G is the single
human/mouse difference inside our block, and step 04 shows it is in cetuximab's
contact set — a real, used surface. A binder touching 442 risks species-specific
behaviour at the one position where the species differ, which would undermine
mouse cross-reactivity. Computed, not assumed.

**New risk: the epitope carries a sugar chain.** N444 is an N-glycosylation site
*inside* 415–466 — step 03 measured 1.44 Å from N444 to a NAG, which is a
covalent bond, not proximity. All accessibility numbers above are for the bare
protein, so for anchors near N444 they are **upper bounds**. D416 (11.3 Å from
N444) is most affected; E421 (18.2 Å) and H433 (23.0 Å) are within reach of an
extended chain. This cannot be settled from a crystal structure, which resolves
only the innermost ordered sugars.

---

## Open

**Is H418 usable? This is now the highest-value open question.**
The two structures disagree. H418 reads as buried in 6ARU (RSA 0.032) but
partially exposed in 1NQL (0.200), with domain III in the same fold in both
(RMSD 1.08 Å). Because 6ARU is an antibody complex, the burial there may be an
artefact of cetuximab holding that side chain rather than an intrinsic property.
Step 03 excluded H418 on the 6ARU reading alone, so **that exclusion is uncertain,
not settled.**

It matters more than any other open item, because H418 is a *target histidine* and
therefore carries the method-novelty claim. Computed 30 September 2026: if H418 is
usable, the largest cluster becomes **five** anchors — D416, H418, E421, E424,
E455 at 24.0 Å span — which beats every 4-anchor option, contains a target
histidine, and is **entirely outside cetuximab's footprint**. H418–H433 is 25.2 Å,
just over the cutoff, so the two histidines cannot both be reached; it is one or
the other.

Resolve by checking additional EGFR structures and by side-chain rotamer analysis,
not by picking the more convenient of the two readings.

**Which anchor cluster to design against.** Four qualify at size 4; the choice is
a real trade-off rather than an optimisation, and it should be made deliberately
and recorded:

| cluster | span | target His? | cetuximab overlap | within 25 Å of N444 | in the domain IV groove |
|---|---|---|---|---|---|
| D416, E421, E424, E455 | 22.6 Å | none | none | D416, E421 | **none** |
| E424, E455, D458, D460 | 22.8 Å | none | none | none | D458, D460 |
| H433, E455, D458, D460 | 23.6 Å | **H433** | H433 | H433 | D458, D460 |
| D416, E424, E455, D460 | 24.6 Å | none | none | D416 | D460 |
| *D416, H418, E421, E424, E455* | *24.0 Å* | *H418* | *none* | *D416, E421* | *none* |

"Within 25 Å of N444" is distance to the point where the sugar chain attaches, not
to the sugars visible in the structure; step 03 found no anchor within 5 Å of a
visible sugar atom. The two measure different things and both are in the findings.

**Recommendation, revised 1 October 2026.** There is no longer a clean winner among
the four-anchor options. `E424, E455, D458, D460` was previously called the
cleanest, on the grounds that it had neither cetuximab overlap nor a sugar-chain
flag. The domain IV finding above removes that standing: both D458 and D460 sit in
the interdomain groove, so it trades one difficulty for another. Every four-anchor
option now carries at least one of three drawbacks — a sugar chain that might
reach, an anchor in a groove, or no target histidine.

The five-anchor option in the last row avoids the groove, carries a target
histidine, and stays clear of cetuximab, which makes it preferable to any of the
four-anchor options on three of the five columns. It depends entirely on whether
H418 is usable. **Resolving H418 is therefore the most valuable thing
outstanding**, more so than before the domain IV finding. Do not pick a four-anchor
cluster until it is settled. If H418 proves unusable, `D416, E421, E424, E455` is
the best remaining option, on the reasoning that a sugar chain that may sometimes
reach is a softer problem than an anchor sitting in a groove — that comparison is a
judgement, not a computed result.

Clusters 2 and 3 differ by a single residue, E424 against H433, around a shared
core of E455/D458/D460.

**Whether a target histidine is worth its cost.** Pairing an acidic binder residue
against the target's own histidine is the method differentiator, but every cluster
containing one carries a penalty: H433 overlaps cetuximab's footprint (inviting
the "this is cetuximab's epitope" objection), and H418's availability is unproven.
The counter-argument for H433 — that cetuximab contacts it with no pH dependence,
so sharing a residue is not sharing a mechanism — is an argument to make in the
write-up, **not a computed result**, and must not be presented as one.

**Which molecule category to target.** Microbinder <40 aa, minibinder 40–100 aa,
large binder >100 aa, nanobody or antibody. Winners are announced per category, so
a less crowded category may be worth choosing deliberately. Not yet chosen, but
the geometry now constrains it: a 4-anchor cluster spans 22.6–24.6 Å, which points
to a **minibinder (40–100 aa)** as the default. A microbinder (<40 aa) is viable
only against the tightest triad (H433/D458/D460 or E455/D458/D460, both 15.6 Å)
and would carry three pairs rather than four — thinner margin on a switch that is
already partial.

**Glycan shadowing at N444.** A sugar chain is attached inside the epitope.
Accessibility numbers for nearby anchors are upper bounds. Not resolvable from a
crystal structure; usable as a tie-breaker between otherwise equivalent clusters.
Prefer clusters away from N444 where the choice is otherwise even.

**What proportion of receptor is open vs closed in the assay buffer.** Step 06
established the epitope is accessible in *both* published conformations, which
removes the specific fear that domain II covers our face. It does **not** give the
equilibrium in the assay, which is the number that would actually matter, and
which two crystal snapshots cannot supply. Reduced from a potential blocker to an
unquantified source of variance.

**Do pH selectivity and mouse cross-reactivity compete for the same surface?** Both
mouse cross-reactivity constrain the same interface residues. Partly answered: all
eight anchors are species-identical, so the pairing positions themselves are
conserved, and the only conflict found is position 442 — the single species
difference in the block, which is also a cetuximab contact. Constraint recorded in
Settled: do not contact 442.

**Compute environment.** Nothing stood up yet. Colab / Modal / RunPod all viable.
This is the main schedule risk — dependency and CUDA problems eat days.

**Design pipeline choice.** BindCraft is the accessible one-command option.
RFdiffusion + ProteinMPNN gives more control at higher setup cost. Not chosen.

---

## Unverified — do not build on these

- **Domain III boundaries.** 310–480 is the working definition used throughout.
  Confirm against the official annotation. Note that step 06 avoided depending on
  this by measuring occlusion as "contacted by residues outside 310–480" rather
  than by naming which domain does the contacting, so the structure-check results
  do not rest on these boundaries being exactly right.

---

## Resolved — moved out of Unverified

**A defect in step 06, found because two scripts disagreed (1 October 2026).**
Steps 04 and 06 both worked out which epitope residues the cetuximab Fab touches,
from the same file with the same 4.5 Å cutoff, and produced different lists. Step
04 was right; step 06 was wrong.

Cause: step 06 held the numbering as two plain dictionaries, one each way round,
and its contact section indexed the UniProt-to-structure dictionary with a
structure residue number. Because the two ranges overlap, that lookup succeeded
and returned a position 48 away from the correct one, raising no error. Verified by
recomputing both directions and checking against step 04: with the correct
direction the two agree exactly.

What it affected, and what it did not:

- **Wrong:** step 06's contact lists. It reported the Fab touching D416 at 2.74 Å;
  the residue at that distance is S464. It reported one residue contacted from
  outside domain III; the real number is 19 in 6ARU and 18 in 1NQL.
- **Unaffected:** every accessibility number, the domain III superposition (1.08 Å),
  the global conformation comparison (23.4 Å), the H418 ambiguity, and all of step
  04 and step 05. Those parts used the numbering the right way round.
- **D416 is not touched by cetuximab.** Step 04's conclusion stands, and the
  cluster comparison that rests on it is unchanged.

Fixed by moving the calculation into `analysis/egfr_common.py`, which both scripts
now call, and by replacing the two dictionaries with a small class whose methods
are named `uniprot_of` and `pdb_of`, so passing the wrong kind of number returns
nothing rather than a plausible wrong answer. Step 06 now also compares its result
against step 04's table and stops the run if they differ; that check was tested by
corrupting the table on purpose and confirming it fails.

**Cetuximab contact residues — RESOLVED 30 September 2026, and the prior claim is
mostly DISPROVED.** Computed in `analysis/04_cetuximab_contacts.py`: any EGFR
heavy atom within 4.5 Å of any Fab heavy atom in 6ARU, using a spatial index.

The speculation was that Q390R, E412D, R414W and K467R explain cetuximab's failure
on mouse EGFR. Computed result: **1 of the 4 confirmed.**

- **K467R — confirmed in contact** (3.39 Å from Fab chain C).
- **Q390R, E412D, R414W — NOT in contact. Disproved.** None of the three is within
  4.5 Å of the Fab. Recorded as disproved rather than dropped.
- **Two the speculation missed:** R377K (3.61 Å) and S442G (3.26 Å) *are* in the
  contact set.

So the computed explanation is **R377K, S442G and K467R** — three of the sixteen
domain III species differences lie in cetuximab's footprint. This replaces the
remembered version.

What this establishes and what it does not: a differing residue inside the
footprint is consistent with causing the species failure and is far better
evidence than an assumption, but it is not proof. Whether a given substitution
abolishes binding depends on how much that contact contributes energetically,
which a distance calculation does not measure.

**6ARU is the structure of record** — the EGFR extracellular region bound to a
cetuximab Fab mutant, and the entry the competition page references. Preferred
over 1YY9 (still cited in the glossary and in `alignment-findings.md` as the older
pointer) because 6ARU contains the whole extracellular region, which is what the
assay uses, rather than domain III alone.

**1NQL structure choice — VERIFIED 30 September 2026.** Proposed as a tethered
EGFR structure and flagged unverified. Checked against the RCSB entry record: it
is "Structure of the extracellular domain of human epidermal growth factor (EGF)
receptor in an inactive (low pH) complex with EGF", X-ray at 2.8 Å, containing an
EGFR extracellular-region entity (612 residues observed, 99.8% identity to human
EGFR) plus a 53-residue EGF entity. Confirmed suitable. Two circumstances worth
carrying: it was solved at **low pH**, and with **EGF bound**. Both are
crystallisation conditions rather than statements about our assay, and step 06
separates ligand occlusion from conformational occlusion for exactly that reason.

**Is 415–466 exposed in the assayed conformation? — PARTIALLY RESOLVED.** The
epitope is accessible in both published conformations tested, and the specific
fear (domain II covering our face of domain III) is not borne out. What remains
unknown is the open/closed *proportion* in the assay buffer; that has moved to
Open as a variance source rather than a blocker.

---

## Ruled out

**H418 as a confidently usable anchor, on the 6ARU reading alone.** Ruled out by
step 03 (RSA 0.032, buried), then *un*-ruled-out by step 06, which found it
partially exposed in 1NQL. Currently ambiguous, not ruled out. Recorded here
because the ruling was made and then withdrawn, and that sequence should be
visible rather than tidied away.

**Q390R, E412D and R414W as the explanation for cetuximab's species failure.**
Disproved by computation — none is in the contact set. See Resolved above.

**The fallback epitopes 394–411 and 331–347 — not needed, not evaluated.** These
were the contingency if 415–466 failed the structure check. It passed, so they were
never characterised. They remain available and untested if the H418 question or the
cluster choice later makes 415–466 look worse than it does now. `analysis/03`–`05`
take the epitope range as constants at the top of each file and would retarget by
changing those.

---

## Next actions, in order

Steps 1–4 of the previous list are **done** — see Settled, "The structure check
passed". Reproducible via `analysis/00`–`06`.

1. **Resolve the H418 question.** Highest value per unit effort of anything
   outstanding: it decides between a 4-anchor cluster and a 5-anchor cluster that
   carries the method-novelty claim with no cetuximab overlap. Check further EGFR
   structures and examine side-chain rotamers. Do not settle it by preferring the
   more convenient reading.
2. **Choose the anchor cluster and record the reasoning.** Use the trade-off table
   in Open. This is a judgement call between novelty and cleanliness, and the
   write-up is stronger for showing it was made deliberately.
3. **Choose the molecule category.** Geometry points to a minibinder (40–100 aa);
   a microbinder is possible only on the 15.6 Å triad.
4. Stand up a GPU environment and get a binder-design pipeline running end to end
   on any target, however bad the output, before attempting the real one. **This
   remains the main schedule risk** and should start in parallel with 1–3, not
   after them.
5. Generate several hundred candidates against the confirmed patch.
6. **Charge-pair scan.** For each candidate, look at which target residue each
   binder contact position faces, then apply the positional rule: histidine
   opposite acidic targets, acidic opposite the target histidine. Keep versions
   where at least three correct pairs form and the fold still holds. Explicitly
   **reject** any candidate with a histidine facing H418 or H433. Additionally
   reject any candidate contacting position 442 (see Settled).
7. **Hold affinity down, deliberately.** Re-read the marginal-not-maximal rule in
   `CLAUDE.md` before filtering candidates, because the pipeline's default
   objective will fight it.
8. Novelty check, rank, submit up to 20 with written reasoning.
