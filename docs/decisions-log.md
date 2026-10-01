# Decisions log

Read this first in any new chat. Update it when something moves between sections.

Last updated: 2 October 2026 (trim boundary measured, E424 ranking term, target
numbering check, smoke-run instrumentation)

---

## The task

Anthropic x Adaptyv Protein Design Competition, Challenge 1 (of 5).
Official page: proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
Deadline 4 October 2026, 23:59 AoE (Anywhere on Earth: the deadline has not passed
until it has passed in every time zone, which is UTC-12).

Four further challenges follow weekly until 1 November 2026. Everything built
here — the analysis scripts, the numbering guard, the findings format — should be
reusable against a new target by swapping the sequences and the epitope, not
rewritten. Treat single-use code as a defect.

Other official links, so they never have to be re-fetched:
novelty rules adaptyvbio.com/blog/novelty · Slack join.slack.com/t/proteinbase
(invite `zt-3evw8fs9z-tU9ItWVvw4ySctUuPvIhLQ`) · models
github.com/anthropics/uplifting-biomolecular-modeling. Full copy in
`rules-reference.md`.

Design binders that:

1. bind the extracellular region of human EGFR (residues 25–645, 621 aa, where aa
   stands for amino acids)
2. also bind mouse EGFR
3. bind human EGFR at pH 6.5 with no detectable binding at pH 7.4

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
based on predicted design quality, design novelty and method novelty,
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
UniProt is the public archive of protein sequences; P00533 is its record for human
EGFR. Every residue number in this project is a position in that full record, which
is 1210 amino acids long, where residues 1–24 are the signal peptide and the mature
extracellular region runs 25–645. The official challenge constructs are the
mature region, so their position 1 is our position 25:

    UniProt position = challenge-construct position + 24

Nothing in a FASTA sequence file or a PDB structure file records which convention it
uses — PDB being the Protein Data Bank, the public archive of measured
three-dimensional structures — so a
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
residues for histidine pairing. Details in `alignment-findings.md` in this
directory; the computed equivalent is `results/findings/01-alignment.md`.
Status: leading candidate, conditional on the structure check passing.

**pH mechanism: charge pairing, decided per position.**
Histidine is the only amino acid that switches charge between pH 7.4 and 6.5
(neutral above, positive below). An interface is a patterned surface and each
pair behaves independently, so the rule is positional:

- opposite an acidic target residue (D/E) → put histidine on the binder
- opposite a target histidine (H418, H433) → put D or E on the binder

Both switch on as pH drops. Do not place histidine opposite H418 or H433 —
that pair switches off in the tumour and can cancel a correct pair elsewhere.

Two known limits: the switch is partial rather than binary at pH 6.5, so use
three or four pairs rather than one; and a histidine's flipping point shifts with its
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
submissions. Not stated anywhere. Ask in the Proteinbase Slack channel.

**What we are actually doing, honestly stated.** Replaces an earlier section that
claimed we were competing where no established method exists. That claim was
withdrawn on 1 October 2026 after a search of the existing literature; see Ruled out.

1. An epitope chosen by computed human/mouse sequence identity, with the constraint
   against contacting position 442 derived from cetuximab's measured contact set
   rather than assumed.
2. A target histidine chosen using the published experimental mapping of EGFR's
   pH-dependence determinants (Liu et al. 2022, Molecular Therapy — Oncolytics
   27:256–269), which identified H370 and H433 by histidine-to-alanine mutagenesis.
   Our target histidine therefore rests on a published measurement rather than on
   reasoning about which histidine looks promising.
3. A positional charge-pairing rule with an explicit rejection criterion — never a
   histidine on the binder facing a histidine on the target — which has published
   experimental support: Liu et al. tried that arrangement and it produced no effect,
   while the acidic pairing worked.
4. An affinity target deliberately held low, because the objective is a detection
   threshold rather than a ratio.
5. A reproducible record in which every number is regenerated by a script and every
   disproved assumption is recorded rather than dropped. Three of this project's own
   claims were disproved by its own computation.

Histidine engineering for pH-dependent binding is established practice (Sarkar et al.
2002; histidine scanning, Schröter et al. 2015) and is cited as background. Pairing an
acidic binder residue against a target histidine is a published strategy (Wei & Sulea,
mAbs 2024; Strauch et al., PNAS 2014; applied on CTLA-4 and VISTA). We claim neither
as ours.

**On the rules.** Using Liu et al.'s mapping of which EGFR histidines carry pH
dependence is using published information about the target. Starting from their
antibody and modifying it would be starting from an existing binder, which the rules
forbid and which we are not doing. Our designs are generated de novo against the
target surface. We say this explicitly rather than leave it to be inferred.

**One sentence of context, if wanted.** pH selectivity is described by Anthropic, a
co-organiser, as among five problems "at the frontier of today's protein design
capabilities, including challenges such as species cross-reactivity, pH-sensitivity,
and peptide-MHC specificity" (anthropic.com/research/claude-uplifts-biomolecular-modeling,
checked 1 October 2026). That is accurate and citable. It does not stretch further:
no organiser page we could reach compares pH-sensitive design with high-affinity
design.

**On tooling, as a verified fact about a program rather than a claim about people.**
BindCraft2 reports binder charge at a hard-coded pH 7.4 as a suppressed readout
(`bindcraft/filters.py`, `REPORTED_PH = 7.4`, metrics `Binder_pI` and
`Binder_Net_Charge`), with no pH term in its objective and no pH setting in its
235-setting catalogue. AlphaFold2 and AlphaFold-Multimer take sequences and database
paths only; RFdiffusion conditions on backbone, contigs, hotspot residues and a closed
set of radius-of-gyration and contact potentials; ProteinMPNN has no pH, pKa or
protonation argument. This carries no novelty argument on its own and is recorded
because it is checkable.

**BindCraft2's input and output schema, read from its own source (1 October
2026).** Read from `PacesaLab/BindCraft2` on GitHub, main branch — the repository
exists and is public, created 13 September 2026, latest release v1.0.3 of 26
September. Recorded here so it never has to be re-fetched, and because
`analysis/10_charge_pair_filter.py` was written against it before any candidate
existed.

*Input.* A campaign is one JSON file, passed as `bindcraft design <file>`. A
`targets[]` entry takes only `name`, `target_path`, `chains`, `hotspots`,
`coldspots`, `weight` and `objective`. `hotspots` is **one comma-separated
string**, not a list, e.g. `"A54,A56,B12-16"`; the chain-letter prefix is
optional and unprefixed numbers address the first selected chain; ranges use a
hyphen. The documentation says to use the residue numbers from your own
structure, so hotspots are in the supplied file's numbering. `target_path` is
resolved relative to the campaign file's own directory. `binder_lengths` takes
`[low, high]` as an inclusive range. Our campaign file is
`design/configs/egfr-domain3-h370.json`, generated by `analysis/09`, with
`design/configs/egfr-domain3-h370.md` beside it holding the provenance — the
JSON is kept free of commentary because BindCraft2 validates campaign keys
against a fixed catalogue and an invented key could make it refuse the file.

*Its structure loader refuses or drops things.* Insertion codes are rejected
outright; only the standard twenty amino acids are accepted; residues missing
backbone atoms are dropped without stopping the run. `analysis/09` checks the
trimmed file against all three and they pass.

*Output.* `<project_folder>/3_Ranked/!_Ranked.csv`, ranked best-first by
`i_pDAE`, with the accepted complexes at `3_Ranked/<design>_seq<n>[_<target>].cif`
and the binder alone at `<design>_seq<n>_monomer.cif`. **The complexes are mmCIF,
not PDB.** There is no `final_design_stats.csv`; that was BindCraft v1.

**Target residue numbering is preserved from the input file**, each receptor
chain keeping the numbering its own structure gave it. That is what makes
`analysis/10` able to translate BindCraft2's output back into our numbering using
the offset `analysis/02` already measured, and it is why the trimmed target keeps
6ARU's numbers rather than counting from 1.

**The target is chain A and the binder is chain B, which is the reverse of
BindCraft v1.** A parser assuming the binder is chain A would read the wrong
molecule and report a full set of plausible nonsense. `analysis/10` identifies
the chains by measuring which one reads as human EGFR rather than by the letter.

**The one thing the output does not contain, which is the thing we need.**
BindCraft2 writes `Interface_Binder_Residues` and `Interface_Target_Residues` as
two separate lists at a 4.0 Å cutoff. It does not record which binder residue
faces which target residue, and the package contains no contact-pair or
contact-map structure anywhere. So the pairing has to be recomputed from the
`.cif`, which is what `analysis/10` does, at our own 4.5 Å cutoff via
`egfr_common.contact_pairs`. Because 4.5 is the looser cutoff, our contact set
should *contain* BindCraft2's rather than equal it, and step 10 checks that
relationship and stops if a residue BindCraft2 reports is missing from ours.

*Metrics.* 59 registered filter metrics plus 20 `Interface_<X>_Count`
composition counts. `Binder_pI` and `Binder_Net_Charge` are confirmed present,
computed at `REPORTED_PH = 7.4` hard-coded in `bindcraft/filters.py`, with no pH
term in the objective — the claim recorded above is verified. The v1 metrics
`InterfaceHbonds`, `ShapeComplementarity` and `Binder_Energy_Score` do **not**
exist in v2: they were PyRosetta-derived and v2 has no PyRosetta. Do not carry
those names over.

*Still unknown and worth checking before the first real run.* BindCraft2
registers a filter metric called `Target_Crop_Length`, which suggests it may crop
the target itself. Whether that interacts with our trim, or renumbers the target
in the output, has not been checked. It matters because step 10's translation
back into our numbering would break silently if the output were renumbered.

*Resources, from its installation documentation.* No single minimum video memory
is stated. What is stated is a per-worker budget of `2.0 × (3.4 GB + 38 kB × N²)`
for a padded complex of N residues, with 4 GB of the card left as headroom, which
gives 7.1 GB at 64 residues and 11.8 at 256. The implied floor of roughly 11–12 GB
for one worker is our arithmetic, not their statement. About 20 GB of disk. The
only published runtime is a forty-trajectory campaign on a GH200 taking 3620
seconds at one worker and 2021 at seven — on hardware we are not using, and per
trajectory rather than per accepted design, so it does not give us a cost.
BindCraft2 is also under a source-available licence rather than MIT, free for
internal research use and restricted from being offered as a hosted service.

(This sentence previously said "objectives 2 and 3". It was written when the
source numbered them affinity first, and it meant mouse and pH. Numbered the way
this repo numbers them — pH, mouse, affinity — it read as mouse and affinity,
which is the opposite of the intent. Objectives are named rather than numbered
throughout for that reason.)

**The novelty claim is withdrawn.** Recorded under Ruled out, with the three reasons
and the evidence under Resolved. What replaces it is in "How to describe the method
honestly" below.

**How to describe the method honestly (1 October 2026).** With the novelty claim
withdrawn, what is left is still worth submitting, and is defensible without any
claim to be first:

- a target patch chosen by computed cross-species conservation, 52 positions with a
  single difference, so mouse cross-reactivity follows from where we aim rather than
  from later repair
- a positional pairing rule with an explicit rejection criterion, which now has
  published experimental support: Liu et al. 2022 found that a histidine on the
  binder facing EGFR's H433 did nothing while an acidic residue there worked
- a deliberately low affinity target, because the requirement is a threshold and not
  a ratio
- a structure check that reports what it cannot settle: H418 ambiguous between two
  structures, the sugar chain at N444, two anchors in a groove against domain IV, and
  no knowledge of the open-to-closed balance in the assay

What is not yet published as far as we could find: a **de novo designed miniprotein**
that is pH-conditional against EGFR. The pH-dependent EGFR molecules in the
literature are antibodies. De novo pH-sensitive binder design exists on other targets
(Ahn et al. 2025, bioRxiv doi:10.1101/2025.09.29.678932, using RFdiffusion with
histidine-biased ProteinMPNN, switching in the opposite direction and not on EGFR).
So any claim should be about **modality and combination**, stated as "we found no
published example of X" rather than "X has never been done", and it should cite Liu
et al. 2022 prominently rather than leaving a reader to find it.

**Strategy: treat the methods write-up as part of the submission.**
In Track 3 we hold no reserved testing slot, so the submission competes in a
pool. Selection runs all submissions through a Claude workflow weighing
predicted design quality, design novelty and method novelty, reading the
documentation submitted alongside the sequences. Organisers state selection
will not rely on a single in silico metric. A clear written argument for the
epitope choice and the positional charge-pairing rule is therefore the main
available edge. Budget real time for it.

**Sequence provenance verified, both species.** The human sequence (UniProt
P00533, residues 25-645, 621 aa) is character-for-character identical to the
official challenge page. The official mouse construct is 623 aa, two longer than
UniProt Q01279 residues 25-645; the two extra residues sit at human-equivalent
position 638, far downstream of domain III. Aligned against the official
constructs, domain III still shows the same 16 differences, 415-466 still shows
only S442G, and all eight anchors (D416, H418, E421, E424, H433, E455, D458,
D460) are identical in both. Checked 30 September 2026.

**Aim for weak binding, deliberately.** The pH-selectivity requirement demands no
detectable binding at pH 7.4, which is a threshold rather than a ratio: the pH 7.4
state has to fall below what the instrument can see at all. A design binding at
10 nM / 200 nM is twentyfold selective and fails, because 200 nM is easily detected.
A design binding at 2 µM with nothing measurable passes. The organisers confirm that
a weak but clearly pH-sensitive binder may outrank a strong non-selective one.
Standard design pipelines maximise affinity by default, which produces exactly
the failing case. **Do not accept the pipeline's built-in objective.** Target a
baseline weak enough that the pH 7.4 state falls under the assay detection floor,
then build the largest switch possible on top of it.

The organisers do not say which assay is used, so **we do not know where the
detection floor actually sits.** The threshold we are designing to is a number we
have not been given. That is why we aim for a wide margin rather than a computed
one, and why the write-up should say the margin is a judgement.

**The structure check passed; the 415–466 epitope is a real surface.** Computed
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
(relative solvent accessibility, the share of a residue's surface that water can
reach, 0.032) — but see the Open section, because 1NQL disagrees. D458 is
classification-sensitive: it lands either side of the boundary depending on which
published normalisation table is used, so its label is not robust. E421 and D458
sit within 0.05 of a cutoff, so read them as near the boundary rather than as the
label they were given.

*Geometry (step 05).* **Four anchors cluster within 25 Å** (max span 22.6 Å),
against a minimum of three. Distances measured between side-chain functional
groups — carboxylate carbons and the imidazole ring centroid — not backbones, and
no anchor needed a CB fallback. Four distinct 4-anchor clusters qualify, so the
design is not forced to one contact set. The tightest viable triad is
H433/D458/D460 at 15.6 Å, which keeps a microbinder in play; below ~18 Å only two
anchors fit, which is too few for a stacked switch.

*Cetuximab's real footprint (step 04).* 24 EGFR residues contact the Fab. It
touches 10 of our 52 residues (19.2%) and exactly one of our eight anchors,
H433. So we sit adjacent to and partly overlapping the druggable surface without
reproducing it. That is the position we were aiming for: near enough to block the
same functional site, while most of our contact positions are ones cetuximab does
not use.

*Conformation (step 06).* The epitope is accessible in both published
conformations, and slightly more accessible in the closed one: mean relative solvent
accessibility across the block is 0.173 in 6ARU and 0.189 in 1NQL. The comparison is
meaningful because the two structures really are in different shapes. Laying domain
III of one onto domain III of the other matches them to within about one angstrom on
average, so the domain itself has the same fold, and once matched that way the
domain I–II region sits 23.4 Å away from its counterpart, which is a large
rearrangement. EGF contacts none of our block.

*The second half of the epitope rests against domain IV (step 06, corrected
1 October 2026).* 19 residues in 6ARU and 18 in 1NQL, spanning roughly 446–466,
are within 4.5 Å of residues outside domain III, and the partners are in the
481–524 range, which is domain IV. The two structures give nearly the same list at
nearly the same distances, so this is a standing feature of how the protein folds
rather than something the closed shape introduces — the original concern about the
tether is still answered. Two anchors, D458 and D460, are in that group.

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
mouse cross-reactivity. This came out of the contact calculation rather than being
assumed.

**New risk: the epitope carries a sugar chain.** N444 is an N-glycosylation site
*inside* 415–466 — step 03 measured 1.44 Å from N444 to a NAG, which is a
covalent bond rather than two things happening to sit close. All accessibility
numbers above are for the bare
protein, so for anchors near N444 they are upper bounds. D416 (11.3 Å from
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
usable, the largest cluster becomes five anchors — D416, H418, E421, E424,
E455 at 24.0 Å span — which beats every 4-anchor option, contains a target
histidine, and is entirely outside cetuximab's footprint. H418–H433 is 25.2 Å,
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
to a minibinder (40–100 aa) as the default. A microbinder (<40 aa) is viable
only against the tightest triad (H433/D458/D460 or E455/D458/D460, both 15.6 Å)
and would carry three pairs rather than four — thinner margin on a switch that is
already partial.

**Glycan shadowing at N444.** A sugar chain is attached inside the epitope.
Accessibility numbers for nearby anchors are upper bounds. Not resolvable from a
crystal structure; usable as a tie-breaker between otherwise equivalent clusters.
Prefer clusters away from N444 where the choice is otherwise even.

**What proportion of receptor is open vs closed in the assay buffer.** Step 06
established the epitope is accessible in *both* published conformations, which
removes the specific fear that domain II covers our face. It does not give the
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

- **Trimming the target to domain III before running the design pipeline.** The
  full extracellular region is about 620 residues, which makes each attempt slow
  enough to exhaust the free compute budget in a handful of tries; domain III is
  roughly 170, which brings a run into the fast case. The failure mode is real and
  specific: cut too tight and the fragment will not hold its shape in the structure
  predictor, so the binder would be designed against a surface that does not exist
  on the intact protein. Keeping the whole domain rather than carving out the eight
  anchors is the hedge against that, but nothing has tested where the limit sits.
  How much can be left out is a judgement until something measures it. Resolve by
  predicting the trimmed fragment alone and comparing it against the same residues
  in the intact structure before trusting any design made against it.

  **The cut is now made, and what it breaks is measured (1 October 2026).**
  `analysis/09_trim_target.py` writes the fragment and
  `results/findings/09-trimmed-target.md` reports it. The judgement above is
  unchanged and still unverified; what has moved is that its cost is now written
  down instead of guessed at.

  The fragment is 171 residues, our 310–480, every one of them resolved in 6ARU
  with no gap, against 621 for the whole extracellular region. All eight anchors
  of the H370 patch fall inside the boundary.

  *One disulfide bond is severed:* C470 inside the fragment is bonded to C499
  outside it, at 2.04 Å. So one cysteine in the fragment loses the partner it has
  in the intact protein. The other 23 bonds in the receptor chain either sit
  wholly inside the fragment (3) or are removed with both partners (20). The 24
  measured bonds are identical to the SSBOND records the depositors wrote into
  the file, so the measurement is checked rather than merely computed.

  *One anchor sits near a cut edge:* E424 is 6.6 Å from the C-terminal cut at
  480. The other seven are 12.3 Å or further from either end. Designs leaning on
  E424 are therefore the least supported, which is a ranking consideration rather
  than a reason to exclude it.

  Neither finding settles the entry. A fragment can fail to hold its shape for
  reasons a severed-bond count and an edge distance do not capture. It is still
  the first smoke run against the trimmed target that moves this.

  **The boundary was measured on 2 October 2026 and the measurement does not
  separate the candidates.** `analysis/11_trim_boundary.py` predicted three
  fragments with ESMFold and compared each against the same residues in 6ARU:
  310–480 as cut, 310–480 with C470 mutated to serine, and 310–499. Results are in
  `results/findings/11-trim-boundary.md`. Over the residues all three share, the
  fit is 5.00, 4.97 and 4.97 Å. Across the eight anchors fitted on their own it is
  0.64, 0.64 and 0.65 Å. Neither differs by more than the predictor can be
  trusted to resolve. What differs is the cut edge: residues 471–480 sit 4.35 Å
  from 6ARU when they are the fragment's end and 2.95 Å when they are 19 residues
  inside it, the predictor's own confidence there rises from 54 to 78, and the
  predicted C470–C499 bond forms at 2.22 Å. Against that, E424 itself came out
  slightly further from 6ARU in 310–499 (2.71 Å) than in 310–480 (2.36 Å), the
  wrong direction for the residue the extension was meant to help, though by an
  amount inside the noise. Mutating C470 changed nothing measurable.

  *Decision: the boundary stays 310–480.* The extension helps residues that only
  one anchor is near, costs about 11% more target, and shows no benefit on the
  anchor face. The deliberate C470S sequence difference buys nothing and is not
  adopted. The risk E424 carries is handled where it bites, in the ranking in
  `analysis/10`, which places a design that needs E424 below an equivalent one
  that does not.

  *What this does not establish.* The overall fit of about 5 Å is the same for all
  three, so the boundary does not cause it, but this measurement cannot say whether
  it reflects the predictor or the fragment. A fit on eight points is also
  flattering by construction. A real answer needs the smoke run's predictor on our
  hardware. The whole question stays Unverified, with weaker grounds for worry than
  before and none for confidence.

  *If the boundary is ever moved:* `analysis/09` takes its cut from
  `egfr_common.D3_START` and `D3_END`, which are also the definition of domain III
  that steps 01 and 06 use through local copies. Moving them would change step 09's
  output and, once the local copies were synchronised, step 01's difference counts
  and step 06's occlusion numbers. A separate constant for the trim boundary would
  avoid that, and should be introduced rather than moving `D3_END`.

- **Domain III boundaries.** 310–480 is the working definition used throughout.
  Confirm against the official annotation. Note that step 06 avoided depending on
  this by measuring occlusion as "contacted by residues outside 310–480" rather
  than by naming which domain does the contacting, so the structure-check results
  do not rest on these boundaries being exactly right.

---

## Resolved — moved out of Unverified

**The method is not novel. Both halves are published, and the closest paper is on
EGFR and on H433 (1 October 2026).** Searched before claiming novelty, on the
reasoning that a novelty claim in front of readers who know the field is worse than
no claim. Three questions were asked; all three came back against us.

*Histidine engineering for pH-dependent binding is established, roughly twenty
years old.* Cite as background rather than presenting as ours:

- Sarkar CA et al. (2002), "Rational cytokine design ... pH-activated 'histidine
  switching'", Nat Biotechnol 20:908–913, doi:10.1038/nbt725 — origin of the term.
- Igawa T et al. (2010), "Antibody recycling by engineered pH-dependent antigen
  binding improves the duration of antigen neutralization", Nat Biotechnol
  28:1203–1207, doi:10.1038/nbt.1691 — the recycling-antibody work.
- Schröter C et al. (2015), "A generic approach to engineer antibody pH-switches
  using combinatorial histidine scanning libraries and yeast display", mAbs
  7:138–151, doi:10.4161/19420862.2014.985993 — the histidine-scanning method.
- Sulea T et al. (2020), "Structure-based engineering of pH-dependent antibody
  binding for selective targeting of solid-tumor microenvironment", mAbs 12:1682866,
  doi:10.1080/19420862.2019.1682866 — computational histidine placement for tumour
  pH, the closest methodological precedent.
- Review for the histidine rationale: Smith FD et al. (2023), Antibodies 12:55,
  doi:10.3390/antib12030055.

*Pairing an acidic binder residue against a histidine already on the target is also
published, and is stated as a known strategy in a review:* Wei W & Sulea T (2024),
mAbs 16:2404064, doi:10.1080/19420862.2024.2404064. Executed deliberately with
structural validation in at least: Lee PS et al. (2022), "Improved therapeutic index
of an acidic pH-selective antibody", mAbs 14:2024642,
doi:10.1080/19420862.2021.2024642 (engineered Asp/Glu against the single histidine of
CTLA-4); the VISTA programme, Johnston RJ et al. (2019), Nature 574:565–570,
doi:10.1038/s41586-019-1674-5 and Thisted T et al. (2024), Nat Commun 15,
doi:10.1038/s41467-024-47256-x; and, as a design principle for de novo interfaces,
Strauch E-M, Fleishman SJ & Baker D (2014), PNAS 111:675–680,
doi:10.1073/pnas.1313605111.

*The closest paper is on our target and our histidine.* **Liu X et al. (2022), "A
cross-reactive pH-dependent EGFR antibody with improved tumor selectivity and
penetration obtained by structure-guided engineering", Molecular Therapy –
Oncolytics 27:256–269, doi:10.1016/j.omto.2022.11.001, PMC9703009.** They mapped
EGFR's own H370 and H433 in domain III as the determinants of pH-dependent binding by
histidine-to-alanine mutagenesis, then deliberately mutated their antibody's LCDR1
Tyr32 to Glu or Asp to form a new acidic-to-histidine pair with H433. Same target,
same pH pair of 6.5 against 7.4, same mechanism, and human/mouse cross-reactive as
well. Patent family WO2024109709A1 (Huahui Health; Sui, Liu, Tian) covers it.

**One finding from that paper helps us rather than hurting us.** Their Tyr32His
variant — a histidine on the binder facing a histidine on the target — changed
neither pH-dependency nor affinity, while Tyr32Glu and Tyr32Asp improved
pH-dependency substantially. That is published experimental support for the rejection
criterion this project derived from first principles: do not put a histidine on the
binder opposite H418 or H433. The rule can now cite evidence instead of an argument.

**Every histidine in domain III, enumerated (1 October 2026).** Prompted by Liu et al.
naming H370. Computed from our own sequences, so the table is ours rather than taken
from the paper:

| Histidine | Conserved in mouse? | Inside 415–466? | Status |
|---|---|---|---|
| H358 | yes | no | previously missed; never considered |
| H370 | yes | no | previously missed; experimentally implicated by Liu et al. |
| H383 | **no** — H383R | no | **unusable**, see below |
| H418 | yes | yes | in our anchor set; exposure ambiguous |
| H433 | yes | yes | in our anchor set; the one Liu et al. used |

**H383 is ruled out and the reason is recorded so nobody rediscovers it.** Human has
histidine at 383, mouse has arginine. Arginine carries a positive charge at both pH
values, so it cannot switch, and the position is not even the same residue between the
species. A pair built against H383 would fail mouse cross-reactivity and would not
switch in mouse at all. It is one of the sixteen domain III differences.

H358 and H370 were missed because both fall outside 415–466, which was chosen as the
longest well-conserved run rather than by looking for histidines. Whether either
deserves an epitope of its own is being evaluated under the time box in Commitments.

**The antibody in 6ARU is a modified cetuximab, and it does not change our
conclusions (1 October 2026).** 6ARU's title calls it a cetuximab Fab mutant, and
step 04 measured the antibody footprint from that file. Computed in
`analysis/07_fab_mutant_check.py` by comparison against 1YY9, the reference
cetuximab structure (Li S. et al., 2005, Cancer Cell 7:301–311).

Five differences, four of which replace a serine or asparagine with aspartic acid:
light chain S52D and S56D, heavy chain S28D, N31D and R216K. **One is in the
interface: heavy chain N31D, 3.93 Å from H433.** The purpose is unrecorded — the
structure is unpublished (Christie M., Christ D., citation "To Be Published"), the
entry lists no substitutions, and the file's SEQADV records cover only the receptor
chain.

Because a modified residue touches one of our anchors, the footprint was recomputed
on 1YY9. It gives **exactly the same ten residues inside 415–466**, H433 included,
at 3.36 Å against 3.47 Å in the mutant. So step 04's contact set describes cetuximab
and not only this variant, and three results that rest on it stand unchanged: the
disproof of the earlier species-failure claim, the rule not to contact position 442,
and the overlap figures in the cluster comparison.

Worth carrying forward for the novelty question: heavy chain position 31 in 6ARU is
an aspartic acid positioned 3.93 Å from H433, which is structurally an acidic binder
residue paired against a target histidine — the arrangement this project treats as
its distinctive contribution. Whether it was placed there for pH-dependent binding
is unknown. It is a precedent for the arrangement either way.

Also recorded from 6ARU's SEQADV records, which we had not read before: the receptor
chain itself differs from UniProt P00533 at two positions (540 asparagine to lysine,
634 glutamate to arginine) and carries a six-histidine purification tag at the
C-terminus. All three lie outside domain III and outside our epitope.

**The H418 cluster numbers are now computed (1 October 2026).** The 24.0 Å
five-anchor span and the 25.2 Å H418–H433 distance previously existed only in this
log, from a calculation run by hand. Step 05 now runs the clustering under both
assumptions — H418 excluded, following step 03, and H418 included, following step
06 — and writes both to `data/derived/`. Section 7 of
`results/findings/05-anchor-geometry.md` carries the conditional set, labelled as
conditional. The numbers are unchanged; they are now reproducible.

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
on mouse EGFR. Computed result: 1 of the 4 confirmed.

- K467R — confirmed in contact (3.39 Å from Fab chain C).
- Q390R, E412D and R414W — not in contact, and therefore disproved. None of the three is within
  4.5 Å of the Fab. Recorded as disproved rather than dropped.
- Two the speculation missed: R377K (3.61 Å) and S442G (3.26 Å) are in the
  contact set.

So the computed explanation is R377K, S442G and K467R — three of the sixteen
domain III species differences lie in cetuximab's footprint. This replaces the
remembered version.

What this establishes and what it does not: a differing residue inside the
footprint is consistent with causing the species failure and is far better
evidence than an assumption, but it is not proof. Whether a given substitution
abolishes binding depends on how much that contact contributes energetically,
which a distance calculation does not measure.

**6ARU is the structure of record** — the EGFR extracellular region bound to a
cetuximab Fab **mutant**, and the entry the competition page references. The mutant
qualification matters and is resolved above: five residues differ from the reference
cetuximab structure, one of them in the interface, and the footprint was rechecked
against unmodified cetuximab as a result. Preferred
over 1YY9 (still cited in the glossary and in `alignment-findings.md` as the older
pointer) because 6ARU contains the whole extracellular region, which is what the
assay uses, rather than domain III alone.

**1NQL structure choice — VERIFIED 30 September 2026.** Proposed as a tethered
EGFR structure and flagged unverified. Checked against the RCSB entry record: it
is "Structure of the extracellular domain of human epidermal growth factor (EGF)
receptor in an inactive (low pH) complex with EGF", X-ray at 2.8 Å, containing an
EGFR extracellular-region entity (612 residues observed, 99.8% identity to human
EGFR) plus a 53-residue EGF entity. Confirmed suitable. Two circumstances worth
carrying: it was solved at low pH, and with EGF bound. Both are
crystallisation conditions rather than statements about our assay, and step 06
separates ligand occlusion from conformational occlusion for exactly that reason.

**Is 415–466 exposed in the assayed conformation? — PARTIALLY RESOLVED.** The
epitope is accessible in both published conformations tested, and the specific
fear (domain II covering our face of domain III) is not borne out. What remains
unknown is the open/closed *proportion* in the assay buffer; that has moved to
Open as a variance source rather than a blocker.

---

## Ruled out

**Strategy: "compete where no established method exists" — withdrawn 1 October
2026.** The section read:

> Affinity optimisation is a solved-ish problem that well-resourced labs will win.
> pH selectivity has no validated predictive tool — nobody can compute it, everyone
> is reasoning from first principles. Cross-species is a sequence-analysis problem
> with a definite answer.
>
> Target-side histidines are the likely method differentiator. Off-the-shelf
> pipelines add histidines to the binder and score the result. Using the target's own
> H418 and H433 as anchors for acidic residues on our side requires reasoning about
> the target's protonation, which no standard pipeline does.

Withdrawn for three separate reasons, each established by checking rather than
argument:

- The premise is false for this challenge. Liu et al. 2022 built a pH-dependent
  anti-EGFR antibody using EGFR's own H433, at pH 6.5 against 7.4, human/mouse
  cross-reactive, with a patent family. An established method exists and has been
  applied to this target.
- Pairing an acidic binder residue against a target histidine is published as a
  general strategy, not only on EGFR: stated in a 2024 review and executed with
  structures on CTLA-4 and VISTA.
- The tooling sentence was refutable by one grep. BindCraft2 does contain
  pH-dependent code. The accurate version is recorded in Settled and carries no
  novelty argument.

The claim is withdrawn rather than narrowed. It is not relocated to the modality: no
argument is made that a miniprotein rather than an antibody constitutes novelty, on
the grounds that such a carve-out reads badly to anyone who knows the field, and the
people judging this do. A withdrawn strategy is a result and is recorded like any
other.

**H418 as a confidently usable anchor, on the 6ARU reading alone.** Ruled out by
step 03 (relative solvent accessibility 0.032, buried), then reinstated by step 06,
which found it
partially exposed in 1NQL. It is currently ambiguous rather than ruled out.
Recorded here because the ruling was made and then withdrawn, and that sequence is
worth keeping visible.

**Q390R, E412D and R414W as the explanation for cetuximab's species failure.**
Disproved by computation — none is in the contact set. See Resolved above.

**The fallback epitopes 394–411 and 331–347 — never evaluated.** These
were the contingency if 415–466 failed the structure check. It passed, so they were
never characterised. They remain available and untested if the H418 question or the
cluster choice later makes 415–466 look worse than it does now. `analysis/03`–`05`
take the epitope range as constants at the top of each file and would retarget by
changing those.

---

## The H370 epitope evaluation — result and recommendation

Computed 1 October 2026 by `analysis/08_h370_epitope.py`, inside the time box.
Defined from the structure rather than from a sequence window, because residues
adjacent in sequence can point in opposite directions and a window is therefore not
a patch: take H370's imidazole ring centre, collect every residue whose side-chain
tip lies within 25 Å, then keep those that are identical in human and mouse,
reachable by water, and able to carry a charge pair.

**Recommendation: switch the primary epitope to the H370 neighbourhood.** It is
better than 415–466 on every measure we have, and the two reasons that matter most
are not close.

| | 415–466 | H370 neighbourhood |
|---|---|---|
| human/mouse differences in the core run | 1 (S442G) | **0** (365–376, `ISGDLHILPVAF`) |
| candidate anchors, conserved and reachable | 7 | **16** |
| largest cluster on one face | 4 within 22.6 Å, 39° spread | **8 within 24.3 Å, 58° spread** |
| target histidines in that cluster | 1 (H433) | **2 (H358, H370)** |
| cluster anchors cetuximab touches | 1 of 4 (H433) | **0 of 8** |
| anchors in the domain IV groove | 2 (D458, D460) | **0** |
| sugar chain inside the region | yes, N444 | none found within the cluster |
| taken by published work | yes — Liu et al. used H433 | **no — H370 unused** |

The eight-anchor cluster is E344, H358, D368, H370, E391, E400, E421, E424 — six
acidic and two histidine, which is enough for the three or four stacked pairs the
switch needs with redundancy to spare, rather than the bare minimum.

*Answers to the three questions asked.*

(a) **Cetuximab overlap: none in the chosen cluster.** H370 itself is not touched.
Cetuximab's contacts near the conserved run are at 373, 374 and 377, so the tail of
the run is inside the drug's footprint while H370 and D368 are not. Of the 16
candidate anchors, the only one cetuximab touches is H433, which is not in the
recommended cluster.

(b) **Species conservation: the run 365–376 has no human/mouse differences at all.**
It ends where it does because both flanking residues differ: S364A before it and
R377K after. Across the wider spatial neighbourhood, three nearby acidic or histidine
residues are **not** conserved and are excluded for that reason: D393E, H383R and
E412D.

(c) **Clustering: yes, with room to spare.** Eight anchors within 24.3 Å.

*A methodological point that was checked rather than assumed.* The clustering test
used through step 05 asks whether anchors fit inside a ball about 25 Å across. A
binder presents something closer to a flat face, so anchors on opposite sides of that
ball pass the distance test and are still unreachable together. With four anchors the
gap between the test and the thing being tested is small; with eight it is not. So
step 08 adds an angular check — the directions the anchors point away from the
protein's centre — and applies it to **both** epitopes, so the new candidate is not
held to a stricter standard than the incumbent. All the 415–466 clusters pass at
24–44°, and the H370 cluster passes at 58°. The 90° limit is a working convention
chosen here and is cruder than docking a real backbone; it removes clearly
unreachable sets rather than guaranteeing the ones that remain.

*What this does not establish.* The same limits as step 05. Cluster geometry is
necessary and not sufficient: it says nothing about whether a foldable binder can
present the right partners in the right orientations, and nothing about whether the
switch clears the assay's detection floor. H358 is a second target histidine and has
no experimental support of its own — only H370 and H433 were mapped by Liu et al.

*The honest caveat on the main argument.* Preferring H370 because H433 is "taken"
is a judgement about how the work will be read, not a computed result. The computed
part is that H370's neighbourhood has more conserved reachable anchors, more target
histidines, no antibody overlap and no groove or sugar problems.

---

## Pipeline status — the critical path

**State as of 1 October 2026: prepared, blocked on one command only the owner can
run.** No candidate sequence exists yet. Under the hard rule below, this is what
everything else yields to.

*Chosen route: BindCraft2 on Modal.* `design/modal/bindcraft2_smoke.py` builds the
image, checks the GPU, and runs BindCraft2 against its own shipped PD-L1 example —
deliberately not our epitope, so our configuration cannot be the cause of a first
failure.

*What was verified about BindCraft2* (its own installation documentation and
`containers/README.md`, read 1 October 2026): Linux, Python 3.12 or later, a GPU
required with no processor-only install, about 20 GB of disk for setup. The installer
detects the GPU and picks accelerator wheels. It ships Docker, Apptainer and Enroot
recipes but **no prebuilt image in a registry**, so the image has to be built. It
uses JAX with CUDA 13 wheels. **PyRosetta is not mentioned anywhere**, unlike
BindCraft v1, which removes a licensing obstacle. Commands are `bindcraft design`,
`fetch-weights`, `rank`, `filter`, `score`. **No minimum VRAM is stated anywhere in
the documentation** — UNVERIFIED, and the card choice below is a judgement rather
than a requirement.

*Why Modal rather than RunPod or Colab* (modal.com/pricing, read 1 October 2026):
per-second billing with idle time not charged, which suits an exercise that is mostly
failed attempts; **$30 per month of free compute on the Starter plan**; and it takes a
container image definition directly, with image builds billed as processor rather than
GPU time. Colab was rejected because sessions time out and disconnect and the GPU type
is not guaranteed. RunPod is a reasonable second choice but bills a running pod by the
hour, so debugging idle time costs money.

*Costs, at prices read 1 October 2026.* Nvidia L4 $0.80/hour, A10 $1.10/hour, A100
40 GB $2.10/hour, A100 80 GB $2.50/hour. The free $30 is roughly 37 hours on an L4 or
14 on an A100 40 GB. Getting to a first candidate should therefore cost nothing. The
script defaults to A100 40 GB because the folding step is the memory-hungry part and
an out-of-memory failure costs more time than the price difference costs money; drop
to L4 once something works.

*Generating a few hundred designs afterwards is the part that may cost money.*
BindCraft2's documentation does not state a runtime per design, so this cannot be
estimated honestly — UNVERIFIED. Measure it on the first successful run and decide
then.

*The Anthropic repository is not a shortcut.* `uplifting-biomolecular-modeling`
provides optimisation kits that wrap upstream tools including BindCraft on
ColabDesign, BoltzGen and others, targeting an H100 80 GB primarily. It accelerates
existing tools rather than offering an easier route, so it is not the fast path to a
first sequence.

*Blocked on:* `./.venv/bin/modal setup`, which opens a browser and can only be done
by the account holder. Nothing can be spent before that.

*What is ready while it stays blocked (1 October 2026).* Both steps that need no
GPU are written and passing: `analysis/09_trim_target.py` produces the trimmed
target and the campaign file BindCraft2 receives, and
`analysis/10_charge_pair_filter.py` scores whatever it returns. So the blocked
step is now the only blocked step, and the first successful run has somewhere to
send its output rather than needing a filter written afterwards.

---

## Commitments made in advance

These are recorded now so they are commitments rather than judgements made under
deadline pressure later. A later session should treat them as settled.

**HARD RULE — the pipeline outranks everything else.** If by end of day 2 October
2026 the design pipeline has not produced a single candidate sequence against any
target at all, stop all epitope work, abandon the H370 evaluation, and spend
everything remaining on generating sequences.

The reasoning: it is 2 October, the deadline is 4 October 23:59 Anywhere on Earth,
and there is no GPU environment and no candidate sequence. Epitope refinement has
absorbed the effort because it is tractable and produces clean results, while the
pipeline work is the part that eats days on dependency and CUDA problems and has not
started. A perfect epitope with zero sequences is not a submission.

**TIME BOX — the H370 evaluation ends on 2 October 2026.** Done by end of that day or
abandoned, and we proceed with 415–466. Recorded so a later session does not reopen
it.

**DECISION — do not read other entrants' work.** Specifically, do not fetch, read or
search the public repository `sepas1609/EGFR-pH-Conditional-Binder-Design`, which a
literature search surfaced and which describes de novo pH-conditional cross-species
EGFR miniprotein binders for this competition. It has not been fetched or inspected.

The reasoning, recorded so this is a decision rather than an oversight: if we read
another entrant's approach and our designs then resemble it, we cannot honestly claim
independent derivation, and the de novo and zero-shot rules make that awkward to
explain. The only thing reading it would protect is a modality-novelty claim we are
dropping regardless. Nothing to gain, and a real contamination risk.

---

## Next actions, in order

Steps 1–4 of the previous list are done — see Settled, "The structure check
passed". Reproducible via `analysis/00`–`06`.

1. **Resolve the H418 question — now a fallback question, not the critical path.**
   It used to be ranked highest because it was said to carry the method-novelty
   claim. That claim is withdrawn, and 415–466 is no longer the primary epitope, so
   the item is rewritten rather than retired: it still decides something, but less.

   What it decides now: whether the fallback epitope offers a four-anchor cluster
   with no target histidine (D416, E421, E424, E455) or a five-anchor one with a
   target histidine and no cetuximab overlap (D416, H418, E421, E424, E455, span
   24.0 Å). That matters only if the H370 patch fails later.

   One thing worth knowing and untested: the same buried reading in 6ARU that
   excluded H418 from 415–466 also excluded it from the H370 candidate set at step
   08's third filter. If H418 is usable it becomes a seventeenth candidate there,
   and whether it would join or extend the eight-anchor cluster has not been
   computed. That is the only route by which this question touches the primary
   epitope.

   Resolve by checking further EGFR structures and examining side-chain rotamers,
   not by preferring whichever reading is more convenient. Ranked below the
   pipeline either way, under the hard rule above.
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
6. **Charge-pair scan — WRITTEN, 1 October 2026, before any candidate exists.**
   `analysis/10_charge_pair_filter.py`. Deliberately written before generation
   rather than after: a filter written afterwards turns out to need a
   measurement the generation run did not save, and there is no time to
   regenerate. Writing the consumer first is what let the schema above be
   checked in advance.

   It applies the positional rule, rejects any candidate with a binder histidine
   facing **any** target histidine rather than only H418 and H433, rejects any
   candidate contacting 442, and ranks by correct-pair count. Six synthetic test
   cases and one synthetic campaign folder are built and run every time, and
   `--self-test` removes each branch of the rule in turn to confirm a test case
   fails without it. All pass.

   What it still needs: real output. The verdicts are untested against anything
   BindCraft2 actually emitted, only against complexes we built to look like it.
   Run it against the PD-L1 smoke run with `--parse-only` as soon as that exists:
   the verdicts will be meaningless because PD-L1 is not our target, but it
   proves the script can read what the pipeline really writes.
7. **Hold affinity down, deliberately.** Re-read the marginal-not-maximal rule in
   `CLAUDE.md` before filtering candidates, because the pipeline's default
   objective will fight it. Step 10 is built to this: it ranks by charge pairs
   only, never discards a candidate for binding weakly, and carries BindCraft2's
   own `i_pDAE` ordering through untouched so it can be used to pick *downward*
   among candidates that already switch.
8. Novelty check, rank, submit up to 20 with written reasoning.
9. **Have the smoke run record wall-clock time per trajectory and peak GPU
   memory**, and write both into a findings file. Explainer 04 originally quoted a
   250-residue trajectory at five minutes, a 900-residue one at two to three hours,
   and a 32 GB minimum video memory. All three came from BindCraft v1's
   documentation and were removed, because BindCraft2's documentation states
   neither a runtime per design nor a minimum video memory. They are the numbers
   the trimming decision and the card
   choice actually rest on, so they come back as measurements from our own run or
   not at all. Time one trajectory on the untrimmed target and one on domain III,
   so the ratio is measured rather than reasoned.

   *Instrumented, 2 October 2026, not yet run.* `design/modal/bindcraft2_smoke.py`
   now records the card type, peak card memory, wall-clock time for the run and per
   trajectory, and the total residue count of the returned complex, into
   `smoke-report.json` in the results volume. The residue count is what lets a
   PD-L1 figure be carried to EGFR. Peak memory is an upper bound, since JAX can
   reserve most of the card at start-up. Per-trajectory time depends on finding a
   table of trajectories in the output; if none is found the report says so and
   prints every table's row count instead of guessing. The same run now answers
   whether the target comes back at the input's numbering, by comparing the target
   chain of every returned complex against the input structure
   (`analysis/target_numbering.py`), and `analysis/10` refuses to score any
   candidate that fails that comparison. Blocked, like everything else, on
   `modal setup`.
10. **Rewrite `docs/explainers/04-design-pipeline.md` after the first real pipeline
   run**, with what actually happened in place of what was planned. Its figures are
   external, read from the web on 1 October 2026, and are a placeholder for
   measurements we have not taken. Any figure that still matters afterwards moves
   into a findings file produced by a script in `analysis/`, so it is reproducible
   like every other number here.
11. **Test how much of domain III the trim can leave out**, which is the unverified
    decision recorded below. It controls both run time and whether the designs mean
    anything, and nothing has measured it. The trim itself is now written
    (`analysis/09_trim_target.py`) and what it breaks is measured — one severed
    disulfide bond, one anchor near a cut edge — but whether the fragment holds
    its shape is still untested and is what this item is about.

    *Partly done, 2 October 2026.* `analysis/11_trim_boundary.py` predicted three
    candidate fragments with ESMFold and compared them against 6ARU. It did not
    separate them; see the Unverified entry for what it did and did not show. What
    remains is the same comparison with the predictor the design run uses, which
    the smoke run on our own hardware can provide.
