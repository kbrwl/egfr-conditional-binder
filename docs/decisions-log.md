# Decisions log

Read this first in any new chat. Update it when something moves between sections.

Last updated: 3 October 2026. Deadline extended to 6 October (Settled); round one
launched as two length-band campaigns (Pipeline status). The Proteinbase Slack was read and its
answers recorded (assay, construct, selection, submission); the marginal-affinity rule
was withdrawn; Modal was linked and its first check found to have tested nothing;
sugar-chain sequons, the tethered cut and the His-tag screen were written; and the
i_pDAE direction was corrected. Then the first campaign against our own target ran:
it starts and the His-tag detargeting works, but every attempt is stopped at
BindCraft2's target-confidence gate by our own fragment, and whether that is the
fragment or the predictor is being tested now (Pipeline status, Next actions item 2).
No candidate sequence exists yet.

---

## The task

Anthropic x Adaptyv Protein Design Competition, Challenge 1 (of 5).
Official page: proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
Deadline **6 October 2026, 23:59 AoE** (Anywhere on Earth: the deadline has not passed
until it has passed in every time zone, which is UTC-12). Settled 3 October 2026 from the
organisers' official Slack announcement; see Settled.

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

**Deadline: 6 October 2026, 23:59 AoE (3 October 2026).** The organisers announced in the
Proteinbase Slack that challenge 1 is extended to 6 October because some teams received
Modal credits late; challenge 2 starts on time, so the two overlap by about two days. This
replaces the 4 October date and the inferred 5 October one. The competition page was to be
updated to match; re-read it, and if it disagrees, ask. See `docs/rules-reference.md` and
`docs/competition-qa-log.md` section 13a.

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
4. A ranking led by charge-pair count, because the pairs produce the switch, with
   BindCraft2's own confidence ordering used only as a tie-break among equal counts.
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

*Measured on real output, 2 October 2026.* Until then this rested on a reading of the
source. `analysis/target_numbering.py`, run inside the smoke run's container on its two
accepted complexes, reports for both: target numbering identical, all 115 residues present
at the input's own numbers with the same amino acids. BindCraft2 did not crop or renumber
the PD-L1 target. That is one structured target of 115 residues with default settings. It
has not been shown for our 171-residue fragment, which is why `analysis/10` still refuses
any candidate whose numbering does not reconcile.

**The target is chain A and the binder is chain B, which is the reverse of
BindCraft v1.** A parser assuming the binder is chain A would read the wrong
molecule and report a full set of plausible nonsense. `analysis/10` identifies
the chains by measuring which one reads as human EGFR rather than by the letter.
*Measured on the same two complexes:* chain A holds 115 residues, the target, and chain B
holds 144 and 165, the binder.

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

**The H370 anchor face is not occluded by the rest of the receptor — settled
2 October 2026.** `analysis/12_fragment_context_check.py`. Step 09 gives the design
run a 171-residue fragment, so the design run cannot see whether the other 450
residues of the extracellular region fold across the face we are aiming at. If they
did, a design would score well and be unable to reach its target, and nothing in the
run's own numbers would say so.

Step 06 asked this of the old 415–466 epitope and found a real answer: D458 and D460
sit in a groove against domain IV. **That measurement never covered the current
anchors.** Step 08 moved the face to the H370 cluster, and only E421 and E424 fall in
the range step 06 measured, so six of the eight positions the design run is aimed at
had never been checked.

Measured now, on the receptor chain alone, in both conformations:

| | 6ARU (extended) | 1NQL (tethered) |
|---|---|---|
| anchors touched from outside domain III, at 4.5 Å | 0 of 8 | 0 of 8 |
| anchors with any outside neighbour within 8 Å | 2 — E421, E424 | 2 — E421, E424 |
| anchors on open face | 6 | 6 |

The six clean anchors — E344, H358, D368, H370, E391, E400 — have 13.0 Å or more of
clear space to the nearest residue outside domain III, and nothing at all within
8 Å. The method reproduces step 06 exactly when run on the old anchors as a control
(D458 ← Q486 at 2.98 Å, D460 ← K487 at 2.77 Å), and that control is cross-checked
against step 06's committed table, so the run stops if the two ever disagree.

Two measurements were kept apart deliberately, because they answer different
questions: contact at 4.5 Å asks whether a residue is covered over, and the count of
neighbours within 8 and 12 Å asks whether a body the size of a binder can get there.
A residue can pass the first and fail the second by sitting in a cleft.

What this changes: the occlusion risk that cost the old epitope two of its anchors
does not apply to this face, in either published conformation. Truncating to domain
III does not hide the anchors. **What it does not settle:** whether the cut fragment
folds the way the intact protein does — a different question, still open, and only
the design run's own predictor can close it. Step 11 could not resolve it with
ESMFold.

Also worth carrying: E421 and E424 are the only two anchors with any neighbour
nearby, and step 09 separately found E424 sits 6.6 Å from the cut at 480. E424 is
therefore the weakest of the eight on two independent grounds, which supports the
rule step 10 already applies of ranking designs that lean on it below equivalent
designs that do not.

Walked through in plain language in `docs/explainers/06-checking-the-trimmed-target.md`,
which also covers the per-candidate screen in `analysis/13_full_receptor_clash.py`.

**What the organisers have said about the assay, the construct and the submission —
recorded 2 October 2026.** Read from the Proteinbase Slack by a person and installed
here as reported speech; nothing in this entry was computed by us. The attributions,
dates and wording are in `docs/competition-qa-log.md` and the facts are in
`docs/rules-reference.md`; where either disagrees with the competition page, the page
wins. In brief:

- *Assay.* Surface plasmon resonance, with **our designs immobilised and the target
  flowed over them**. Top concentration 1000 nM, reportable range roughly 0.1 nM to
  10 µM. Binding is declared when a K_D can be fitted or the association signal clears
  300% over the negative control. A large K_D shift also qualifies, and no binding at
  pH 7.4 with high affinity at pH 6.5 ranks higher.
- *Buffer.* Matched ionic strength of about 170 mM at both pH values, so any difference
  between them comes from the pH.
- *The target.* The tethered conformation, expressed in HEK293 cells, glycosylated, with
  a C-terminal His tag on both species.
- *The designs.* Built with a C-terminal tail of linker, GFP11, linker and twin-Strep tag,
  immobilised at the C-terminus. Linear chains only; cyclic peptides are not supported.
- *Novelty.* The gate is level 3 of 4, scored automatically on upload, checked by MMseqs2
  against SwissProt, the PDB, the USPTO and EBI patent databases, THPdb, PLAbDab and
  Proteinbase.
- *Selection.* Mainly method novelty, design diversity and a couple of computational
  confidence metrics, the last being the smaller component. Explained in
  `docs/explainers/07-the-assay-and-what-it-changes.md`.
- *Submission.* Up to 20 designs, each chain 10 to 250 residues, **one submission every
  24 hours**. Track 3 needs no confirmation email.

**Allocation follows from that.** Design diversity is a selection axis in its own right,
so the twenty slots should spread across anchor subsets and size categories and not take
the top twenty of one ranked list. And a safe early batch can be uploaded and replaced
later, so there is no reason to hold everything for a single upload.

**The tethered cut: 6ARU stays, with one anchor unjudged — measured 2 October 2026.**
The organisers say the screen uses the tethered form and our fragment is cut from 6ARU,
the extended form. `analysis/16_tethered_fragment.py` cuts the same residues 310–480
from 1NQL, the tethered structure, and compares them (`results/findings/16-tethered-
fragment.md`). The rule for deciding was fixed before the measurement: the cuts are the
same where it matters if no anchor's charged tip moves by 2.0 Å, a third of the 6 Å
reach step 10 allows between two charged groups. That tolerance is a judgement.

The whole fragment fits at 1.08 Å over 171 residues, matching step 06. The 79 residues
within 8 Å of an anchor tip fit at 0.69 Å. The charged tips of seven anchors move 0.3 to
1.2 Å, the largest being E391 at 1.19 Å. So 6ARU stays, for that measured reason.

**E344 is a gap and not a pass.** Its charged tip is not resolved in 1NQL, so the rule
cannot be applied to it. On CB, the nearest atom both structures share, it moves
2.14 Å, and its backbone moves 1.79 Å, about three times any other anchor's. The 1NQL
cut is kept on the shelf at `data/structures/1nql_domain3.pdb` with the same numbering.
H418 reads 0.200 in 1NQL against 0.032 in 6ARU, and under the tethered form the 1NQL
reading is the relevant one; that belongs to the 415–466 fallback and changes nothing
for the H370 face. Explainer: none, because the comparison did not change the cut.

**Sugar chains: all eight anchors are within reach of an attachment point, and N361 is
human-only — measured 2 October 2026.** `analysis/14_glycan_sequons.py`
(`results/findings/14-glycan-sequons.md`, explained in
`docs/explainers/08-sugar-chains-near-the-anchors.md`). The human extracellular region
has 11 sequons, N-X-S/T. Measured from each anchor's charged tip to each sequon's
attachment nitrogen, in both structures:

- E344, H358, D368, H370 and E391 are 10.4 to 14.3 Å from N352 in the tethered
  structure, and a sugar is seen bonded to N352 in both structures.
- E400, E421 and E424 are 16.2 to 16.9 Å from N361, also bonded in both.
- No anchor is beyond 25 Å of a sequon. D368 is 14.3 Å in 1NQL and 15.3 Å in 6ARU, either
  side of the 15 Å band line, so its band is borderline.
- **N361 is a sequon in human and not in mouse.** Mouse has a tyrosine there. N361Y was
  already one of step 01's 16 domain III differences and its effect on sugar chains was
  not drawn until now. The human target carries a chain 16–17 Å from E400, E421 and E424
  that the mouse target lacks, which bears on the mouse cross-reactivity ratio.

The bands are step 05's and rest on a reach of 20–30 Å from memory (Unverified).
Distance to an attachment point is a necessary condition for a chain covering a site and
not evidence of it. The controls reproduce step 03's N444 distance (1.44 Å) and step
05's seven N444 distances to 0.1 Å. **Nothing was dropped at this step:** the work order
expected a few flagged anchors and all eight are flagged, including both target
histidines, so demoting by band would reorder the whole campaign and favour the three
anchors near the human-only chain — that trap is real and ranking by band was rejected
for it.

**Decided 2 October 2026: demote by which sequon is nearest, not by band — the demerit
attaches to N361, not to N352.** The two sequons differ in kind rather than only in
distance. N352 sits in both species, so a chain there costs absolute affinity in human
and mouse alike and leaves the mouse-to-human K_D ratio untouched, which is what mouse
cross-reactivity is scored on; N352 proximity is ignored for ranking. N361 exists in
human and not in mouse, so shielding there moves that ratio directly, which is the one
place a cheap ranking term buys something. The five anchors nearest N352 — E344, H358,
D368, H370, E391 — carry no demerit; the three nearest N361 — E400, E421, E424 — do,
demoted and never excluded, in the same shape as the existing E424 edge-reliance rule.
Implemented in `analysis/10_charge_pair_filter.py` as `N361_SEQUON_RELIANT_POSITIONS`,
with its own test case and ranking assertion. Both target histidines, H358 and H370,
sit on the N352 side, so the half of the pairing rule needing an acidic binder residue
is unaffected, and E424 already carried the edge-reliance demotion on the unrelated
ground of its distance from the trimmed target's cut edge, so it is now demoted on two
independent grounds rather than double-counted within one. This is a tie-break built on
a plausible asymmetry, not a measured effect: nothing measures whether a chain at 16 to
17 Å actually reaches E400, E421 or E424, and the 20–30 Å reach the bands themselves
rest on is still Unverified. Explained in
`docs/explainers/08-sugar-chains-near-the-anchors.md`, section 6.

**The binder's C-terminus is kept out of the interface by a BindCraft2 setting.**
`termini_accessible` points both chain ends away from the target (a loss term and a final
floor of 0.0 on the terminus direction cosine, where +1 is away), so the C-terminal tail
and the chip are not in the way. It is set in both campaign files written by
`analysis/09_trim_target.py`. The work order's fallback of a filter in step 10 is not
needed. It is geometry only. Linear chains only is recorded so no later challenge
rediscovers it; nothing in the current plan conflicts.

**The His tag is avoided inside the design loop and audited afterwards.** Both targets
carry a C-terminal His tag the organisers expect to leave on, and a binder that grips it
looks pH-selective and binds anything with a His tag; our pairing rule makes this project
more exposed than most. `analysis/15_histag_counterscreen.py` was written before any
candidate exists. The predictor is BindCraft2's own AlphaFold2 used through its
detargeting objective, the tag given as the sequence GSHHHHHHGS at weight −0.5.
`bindcraft score` was rejected because it reads existing structures and cannot predict a
complex. "Accepts the tag" is `i_pTM_detarget` of at least 0.4, **BindCraft2's own default
and not a threshold calibrated for this tag.** A blank is never a pass, and a pass is weak
evidence. `design/configs/egfr-domain3-h370-notag.json` is the campaign with the tag as an
off-target; it has **not been run**. Explained in
`docs/explainers/09-the-his-tag-screen.md`.

**The i_pDAE direction was wrong in our own documentation.** BindCraft2 defines it as a
distance-masked interface TM confidence between 0 and 1, higher is better; its source keeps
a list of lower-is-better metrics and i_pDAE is not in it. The glossary said lower was
better, which is true of `i_pAE`, and the tie-break added to step 10 earlier that day sorted
the wrong way round until this was caught and corrected the same afternoon, with a test
that fails on the old direction.

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

**What to do about sugar chains near the H370 anchors — decided, see Settled.** All eight
are within reach of an attachment point. Decided 2 October 2026 to demote by which sequon is
nearest (N361 side demoted, N352 side not) rather than by distance band, implemented in
`analysis/10`. Disclosing the risk in the write-up happens regardless, since the ranking term
is a design preference and the underlying uncertainty — whether a chain actually reaches these
anchors — is not resolved by it. Coldspots and the human/mouse joint-design mode (next item)
were considered and not taken for this purpose: coldspots discourage contact with the flagged
residues themselves without modelling a chain's sweep beyond them, and the joint-design mode
needs a mouse domain III structure we do not hold.

For the old 415–466 epitope the equivalent entry was N444, 11.4 Å from D416 and bonded to a
sugar; it stays as a tie-breaker between otherwise equal clusters if that fallback is ever
used.

**Design against human and mouse together?** BindCraft2 has a multi-target mode for one
binder that engages two orthologs, with hotspots given separately in each structure's own
numbering; its own example is one binder for human and mouse PD-L1. It is aimed at the mouse
cross-reactivity objective, which the organisers define as a mouse-to-human K_D ratio of
about 1, and it could bear on the N361 difference. It needs a structure of mouse domain III,
which we do not hold. All eight anchors are identical in the two species, so the hotspots
carry over. Whether the gain is worth the setup is untested.

**Is E344 usable as a tethered-form anchor?** Its charged tip is unresolved in 1NQL and its
backbone and CB move 1.79 and 2.14 Å between the two forms (Settled). Candidates leaning on
it rest on a position that differs between forms. Could be handled like E424, in the ranking.

**Can the burial of a charge pair be measured well enough to rank on?** Explainer 07 argues
that at roughly 170 mM ionic strength a charge pair in the closed-off core of an interface
contributes much more than the same pair on the exposed rim, because dissolved ions screen
exposed charges. That is a design preference and not a measurement. *Measurable:* yes, with
tools already used here. Per-residue accessibility in the complex and the fraction of area
buried on complex formation both come from Shrake–Rupley, as in steps 03 and 06. Probed on
the three EGFR complexes on disk, the only oppositely charged cross-chain pairs within 4.0 Å
were five (three distinct), with accessibility in the complex from 0.00 to 0.30 and the
fraction buried from 0.30 to 1.00. *Testable:* no. Five pairs are too few to separate buried
from rim, nothing ties a pair's burial to the strength of a pH switch, and step 10's
constructed complexes hold isolated residue pairs with no environment around them. **No
ranking term was built.** Revisit when real complexes from the design run exist; until then
it would be a term that cannot be tested.

**The exact mouse construct.** Residue range and vendor, asked twice in Slack on 1 October
and not answered. Our check was against the sequence on the competition page.

**Is cynomolgus monkey cross-reactivity in scope? — decided 2 October 2026: design for
human and mouse only.** Amir's pre-launch message on 28 September said "mouse and cyno". The
competition page, read on 2 October, names human and mouse only, and so does
`docs/rules-reference.md`. The page is authoritative, and it names two species, so no GPU time
or campaign slot is spent on cyno.

The cheap check this entry deferred has now been done, 3 October 2026:
`analysis/19_cyno_conservation.py`, `results/findings/19-cyno-conservation.md`. All eight
H370-cluster anchors hold the same amino acid in *Macaca fascicularis* as in human, and the
attachment points for sugar chains nearest them — N352, N361 and N413 — are present in cyno
too, N361 included, which is the one human has and mouse lacks. The offset between the cyno
entry's own numbering and ours is zero, measured across all 171 positions of domain III. The
only difference anywhere in domain III is S348T, which is not an anchor and is the same
position that differs in mouse. **This changed no decision and was not expected to.** Cyno is
not an objective, no GPU time or campaign slot went to it, and the demotion of E400, E421 and
E424 over N361 stands, because that asymmetry is between human and mouse and mouse
cross-reactivity is the objective. What it buys is a sentence the methods write-up may use,
with the limits the finding records: the cyno UniProt entry is unreviewed, this is sequence
analysis with no cyno structure examined, and nothing about binding was measured.

**What does novelty level 3 mean?** Simon Dürr said level 3 clears the gate (29 September and
1 October) and on 1 October also described the requirement as the sequence being under 30%
similar to anything existing. The two may not be the same threshold. The definition is at
`adaptyvbio.com/blog/novelty` and has not been read. Novelty is scored automatically on
upload, so it can be tested against the real checker before the deadline.

**Are extra metric columns allowed in the submission CSV?** Asked on 1 October, unanswered.

**Do pH selectivity and mouse cross-reactivity compete for the same surface?** Both
mouse cross-reactivity constrain the same interface residues. Partly answered: all
eight anchors are species-identical, so the pairing positions themselves are
conserved, and the only conflict found is position 442 — the single species
difference in the block, which is also a cetuximab contact. Constraint recorded in
Settled: do not contact 442.

**Compute environment.** Modal, linked on 2 October. The GPU check passes on an NVIDIA
A100-SXM4-40GB, driver 580.95.05, with JAX 0.11.2 on the GPU backend (Pipeline status).

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

- **How far a sugar chain reaches.** Step 05's bands of 15 and 25 Å rest on a reach of 20–30 Å
  from a complex glycan, from memory. They are cautious on purpose. The ordering of the anchors
  by distance is more reliable than the labels.

- **Detargeting against a sequence target — the first half is answered, 2 October 2026.**
  That BindCraft2 accepts the His tag as an off-target alongside a structured target is no
  longer unverified: its reference says `targets[].target_path` takes "PDB, mmCIF or FASTA",
  and the validation run loaded `HisTag (detarget) | 10 residues` beside the structured
  target and designed against both. The tag was avoided, at `i_pTM.HisTag = 0.05` on the one
  trajectory that ran. **What is still unverified is the cost and the rate:** what
  detargeting does to run time and to the share of attempts accepted cannot be read from a
  single trajectory, and no trajectory has yet been accepted at all. The weight of −0.5 and
  the ceiling of 0.4 remain starting values and nothing has tuned them.

- **Domain III boundaries.** 310–480 is the working definition used throughout.
  Confirm against the official annotation. Note that step 06 avoided depending on
  this by measuring occlusion as "contacted by residues outside 310–480" rather
  than by naming which domain does the contacting, so the structure-check results
  do not rest on these boundaries being exactly right.

---

## Resolved — moved out of Unverified

**Which conformation is tested, which assay, and what counts as no detectable binding —
answered by the organisers, 29–30 September (recorded 2 October 2026).** These were gaps in
the official information and are now stated: the tethered conformation, surface plasmon
resonance, and the 300% association rule (Settled; `docs/competition-qa-log.md`). The Open
item on the balance of conformations in the assay is closed by this. What remains of it is
that our fragment is cut from the extended structure, measured under Settled.

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

**Rule: "hold affinity down, deliberately" — withdrawn 2 October 2026.** The rule,
carried in `CLAUDE.md` and built into `analysis/10`, was that affinity should be
deliberately marginal. The reasoning was that the pH-selectivity requirement asks for
no detectable binding at pH 7.4, so a strong binder risks being strong at both pH
values while a weak one that clearly switches passes. Step 10 therefore ranked by
charge pairs alone and carried BindCraft2's `i_pDAE` ordering through untouched so
it could be used to pick downward.

Withdrawn because the organisers answered the question directly. Asked on 30
September whether no binding at pH 7.4 is a hard requirement or whether a large K_D
shift also qualifies, Tudor answered that a large shift also qualifies and that
designs with no binding at pH 7.4 together with high affinity at pH 6.5 rank higher
(`docs/competition-qa-log.md`). The quantity being rewarded is the gap between the
two conditions, with the pH 6.5 end as high as it can be.

The instrument adds a second reason. The target is flowed at a top concentration of
1000 nM and K_D is reportable over roughly 0.1 nM to 10 µM. Occupancy at that
concentration is [analyte] / ([analyte] + K_D) by the definition of K_D, computed in
`analysis/17_occupancy_table.py` (`results/findings/17-occupancy.md`): roughly 91% at
a K_D of 100 nM, 50% at 1 µM, 25% at 3 µM and 9% at 10 µM. So a deliberately weak
design can read as no detectable binding at pH 6.5 as well. That fails human binding
instead of demonstrating pH selectivity.

The aim was right and the route was wrong. Weakness was treated as the way to
guarantee silence at pH 7.4, when what is rewarded is the gap.

**Replacement:** maximise the difference between pH 6.5 and pH 7.4, with affinity at
pH 6.5 as high as the switch allows (`CLAUDE.md`). Charge-pair count still leads the
ranking and `i_pDAE` is now an upward tie-break. **Unchanged:** the positional
pairing rule, rejection of a binder histidine facing any target histidine, rejection
of contact with position 442, the `MIN_CORRECT_PAIRS` floor of 3, `PAIR_TARGET` of 4,
and the E424 edge-reliance demotion. Nothing the organisers said touches them.

Explained in `docs/explainers/07-the-assay-and-what-it-changes.md`.

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

**State as of 2 October 2026, 14:45 IST: Modal linked, the GPU check passes, the PD-L1 smoke
run is running.** No candidate sequence exists yet. Under the hard rule below, this is what
everything else yields to.

*What happened on 2 October.* Modal was linked by the owner (workspace `kunaal11791`, $30
credit). **The first `check` exited successfully having tested nothing.** The file
`design/modal/bindcraft2_smoke.py` defined the image and never passed it to the app, so both
functions ran in Modal's default container, where `/opt/BindCraft2` does not exist. It also
finished in about a minute, too fast for a 20 GB install, which is what pointed at it. Fixed:
the app takes the image, `_run` no longer activates BindCraft2's environment for commands
outside it, and `check` now raises when the GPU, JAX or the `bindcraft` command is missing.
The rerun shows an NVIDIA A100-SXM4-40GB, driver 580.95.05, JAX 0.11.2 on `CudaDevice(id=0)`
with default backend `gpu`, and `bindcraft design --help` working. The image build,
including BindCraft2's `install.sh`, completed.

The PD-L1 smoke run (`design/modal/bindcraft2_smoke.py::smoke`) was launched detached at
14:17 IST with a four-hour ceiling, at most about $8.40 of the $30. Its console output only
appears when the command finishes, because `_run` captures it, so progress is not visible.
**Its end-of-run report is not in this log yet.** Read by hand from the live container at
about 14:46 IST, 29 minutes in, and to be superseded by the script's own
`smoke-report.json` when the run ends:

- *Works end to end.* Trajectories, refolded complexes and ranked designs all exist. Of 7
  trajectories finished, 2 were accepted. The ranked table has 43 columns and every column
  step 10's constructed table assumed is present.
- *Two workers share the one card.* Peak card memory sampled at 1 Hz was 17,950 MiB of
  40,960 MiB, so the memory is an upper bound for two workers and may include memory JAX
  reserved. An A100 40 GB is not shown to be needed; a smaller card is not shown to be
  enough, and no speed ratio between cards was measured.
- *Time per trajectory, design phase only,* for a 115-residue target: 185 s at binder length
  87, 245 s at 67, 310 s at 94, 516 s at 165 and 574 s at 178. The two trajectories that
  included compilation (62 and 144 residues, 207 s and 634 s) are left out because
  compilation inflates them. Validation and refolding are not in these figures.
- *Throughput observed,* including everything: 7 trajectories in about 27 minutes, roughly
  15 per hour. Two of 7 accepted is 29%, from a sample too small to carry to EGFR.
- *Complex size.* 259 and 280 residues (115 target plus 144 and 165 binder). The EGFR
  fragment plus a 60–100 residue binder is 231–271, so these timings are the right order for
  the fragment. The full 621-residue receptor would give 681–721, about two and a half times
  as large, and its slowdown was not measured.
- *i_pDAE is ranked higher-first*: rank 1 has 0.22 and rank 2 has 0.19 in `!_Ranked.csv`,
  matching the source and the correction under Settled.

**What that means for the campaign size.** The campaign file asks for up to 2,000
trajectories and 200 accepted designs. At about 15 trajectories an hour that is more than
130 hours of card time, and the $30 buys about 14 hours of an A100 40 GB at $2.10 an hour,
roughly 200 trajectories and, at the PD-L1 acceptance rate, around 60 accepted designs. The
EGFR rate is unknown and may well be lower. The counts in `analysis/09` were never derived
from a measured run time, as its note says; they should be reset to the budget before the
campaign starts (Next actions, item 4).

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

*No longer blocked.* Modal setup was the one step only the account holder could do, and it
is done. Two read-only helpers, `inspect_cli` and `shell`, run on the processor with no GPU
and were used to read BindCraft2's own documentation and source before configuring a paid
run from them.

*What is ready while it stays blocked (1 October 2026).* Both steps that need no
GPU are written and passing: `analysis/09_trim_target.py` produces the trimmed
target and the campaign file BindCraft2 receives, and
`analysis/10_charge_pair_filter.py` scores whatever it returns. So the blocked
step is now the only blocked step, and the first successful run has somewhere to
send its output rather than needing a filter written afterwards.

*The PD-L1 smoke run was found still live, checked 2 October 2026, after this
section was last written.* `./.venv/bin/modal app list` shows every other app from
that day's work as `stopped`, except `ap-zIp8z3RQdOrgZU0CNvEXIW` (`bindcraft2-smoke`),
state `ephemeral (detached)`, 1 task, created 14:16 IST, no stop time recorded. This
is the same run already described above, launched detached at 14:17 IST with a
four-hour ceiling (`timeout=60*60*4` on the `smoke` function); only a partial check
at 14:46 IST (2 of 7 trajectories accepted, 29 minutes in) was ever captured, and
nobody went back for its final `smoke-report.json`. Whether it is still progressing
or hung, and its actual run time and cost, is being checked now via its logs.
**No EGFR campaign is running or has ever been launched; this is the unrelated
PD-L1 example only.**

*It finished on its own at 15:49:34 IST, before the stop command landed.* Full
report, superseding the 14:46 partial check: card NVIDIA A100-SXM4-40GB, peak
memory 18,218 MiB (confirms the L4 switch above — 24 GB is comfortably enough),
whole run 5540 s (~92 minutes) for 30 trajectories attempted, 184.7 s/trajectory,
81 redesign candidates, **10 of 30 trajectories accepted (33%)**, 2 terminated at
the anneal phase. Complex sizes ranged 127-434 residues. Cost at the A100 rate:
roughly $3.23, well inside the $8.40 ceiling the run was given.

**A bug found in the smoke script's own final verdict, not in the pipeline.**
Its printed conclusion said the target "did NOT come back at the input's
numbers" and warned not to start an EGFR campaign — wrong. Of 95 `.cif` files
the numbering check inspected, 91 (every real designed complex, across
trajectories, refolding and the final ranking) read "identical: all 115
residues present at the input's own numbers". The other 4 — `scaffolds/ARP.cif`,
`Fab.cif`, `VHH.cif`, `scFv.cif` — are BindCraft2's own internal template
structures for other binder modalities (an alternative repeat protein, an
antibody Fab, a nanobody, a single-chain antibody fragment), not anything this
run designed. `smoke()`'s output-copy step globs every `*.cif` under the whole
BindCraft2 install directory with no further filter, so it swept these up along
with the real output, and the numbering check dutifully reported that a
nanobody scaffold does not resemble PD-L1 — true, and meaningless. The real
finding stands uncontaminated: **the target does come back at the input's own
numbers**, for every actual designed complex. Fix before writing the real
EGFR campaign runner: restrict the file copy, or at least the numbering check's
input list, to the run's own output directories (`1_Trajectories/`,
`2_Refolded/`, `3_Ranked/`) rather than globbing the whole install tree.

**The first EGFR campaign runner exists, and the first run against our own target
started up — 2 October 2026, about 17:15 IST.** `design/modal/egfr_campaign.py`, a
new file rather than an edit to the smoke script, so the arrangement known to work
stays available to go back to. The two restrictions above are both implemented, in
two separate places, because they fail differently. Explained in
`docs/explainers/10-running-the-campaign.md`.

*Next actions item 2 is answered for the notag file.* `egfr-domain3-h370-notag.json`
starts against the real fragment. BindCraft2 loaded both targets —
`EGFR_domain3 (target) | 171 residues` with our eight hotspots, and
`HisTag (detarget) | 10 residues` — resolved `termini_accessible` as a design
property, and drew a binder length of 79 from the 30–100 band. Cost $0.20 for one
trajectory, 890 s wall-clock including model compilation. **The plain file has not
been run**, because the notag one started and is the one to prefer.

*Three things checked from BindCraft2's own documentation before spending, which
move Unverified items.* `targets[].target_path` takes "PDB, mmCIF or FASTA", so a
sequence off-target is a supported shape; `crop_fasta_sequence: false` means use the
whole sequence; and `termini_accessible` is valid in a campaign file even though it
is absent from the 134 names `--list-settings` prints, because that list covers only
what `--set` accepts. The last one is recorded because it looks like a defect and is
not.

*The detargeting works.* `i_pTM.HisTag = 0.05` — the binder did not engage the tag.
One trajectory, so this is an existence result rather than a rate.

*Measured on the L4, which had never been measured.* Peak card memory 8,915 MiB of
23,034 MiB, less than half the card, so the earlier arithmetic suggesting roughly
12.4 GB a worker was pessimistic and two workers should fit. The run used one worker
only because BindCraft2 holds fan-out to the trajectory budget. 890 s for one
trajectory with compilation included, against 184.7 s a trajectory on the A100 smoke
run with two workers, so **no clean speed ratio between the cards exists yet**: the
two figures differ in card, worker count, compilation and target size at once.

**The campaign is blocked at BindCraft2's target-confidence gate, and the cause is
the target rather than the binder.** The trajectory was stopped at the screen stage,
the first and cheapest, on `pLDDT.EGFR_domain3 = 0.35` against the filter
`Target_pLDDT >= 0.6`. That is the predictor's confidence in our own 171-residue
fragment inside the complex. BindCraft2 said so itself: "1 trajectories ran and none
were accepted, so the settings rather than the budget are what to change."

*The first reading of this was wrong, and the correction is the useful part.* It was
argued that the gate would stop essentially every attempt, on the grounds that the
same figure on the PD-L1 smoke run was 0.93 with a standard deviation of 0.01 across
ten accepted designs, making it nearly a property of the target. **That inference
does not hold.** The 0.01 was measured across accepted designs only — a set already
filtered by the very gate in question — and was then applied to an unfiltered
screen-stage sample. Survivor variance is not sample variance, and this repository's
own rule that a comparison must be honest about its method is exactly what the
argument broke.

Disproved by measurement the same evening: the relaxed-floor diagnostic ran the same
nominal trajectory against the same fragment and reported `pLDDT.EGFR_domain3 = 0.72`
at screen, where the validation run had reported 0.35 and been rejected. So target
confidence varies substantially between runs of the same trajectory, and the 0.6 gate
rejects **some unknown fraction** of attempts rather than all of them. A pilot against
the committed fragment would have produced candidates, more slowly and more wastefully
than planned, rather than none.

What stands: the pilot was still right not to launch, because its size rested on rates
nothing had measured. What falls: the reason given for it. Measuring the distribution
of target confidence is now what the relaxed-floor diagnostic is for, and it is the
number the boundary decision needs.

*Two explanations were offered and one of them rested on reading the wrong metric.*
**One:** the fragment does not hold its shape, which is the standing Unverified
entry — the cut severs the C470–C499 disulfide. **Two:** the predictor is the
problem. BindCraft2 runs AlphaFold with no evolutionary information (`extra_msa` is
zeroed in `bindcraft/af2.py`), which is far harder for a large beta-solenoid than for
PD-L1's compact and much-studied fold, and ESMFold rates this same fragment 80.6
over the whole fragment and 82.9 across the anchors
(`results/findings/11-trim-boundary.md`), above the usual reliable-backbone
convention of about 70. That much stands, and the disagreement between the two
predictors about this fragment is not resolved.

**What does not stand is the argument from the filter's default name.** It was
claimed that the floor comes from `DEFAULT_DISORDERED_TARGET_PLDDT = 0.6` and is
therefore meant for targets supplied as a bare sequence rather than as a structure
whose coordinates are held. That constant governs `Target_pLDDT`, which is a
different metric from the one the rejections name, so the argument was about
something that was not happening.

**The two metrics, read from `bindcraft/filters.py` and `settings/core/reference.json`
on 2 October 2026, because conflating them produced a wrong diagnosis once already:**

| metric | gated by | thresholds |
|---|---|---|
| `pLDDT.<target>` | `min_plddt_<stage>` (`filters.py:826`) | screen 0.6, refine 0.6, mutate 0.6, anneal 0.65, harden 0.65, final 0.7 |
| `Target_pLDDT` | `min_target_plddt_final` | final stage only; default 0.6 |
| `i_pTM.<target>` | `min_iptm_<stage>` | screen and refine ungated, anneal 0.5, harden 0.5, mutate 0.5, final 0.7 |

Every rejection seen so far names `pLDDT.<target>` or `i_pTM.<target>`, so the
per-stage settings are what bind and `Target_pLDDT` has never been reached.
**The interface floors at the middle stages are 0.5 and not 0.7;** 0.7 applies only at
the final stage. An earlier note in this log read the anneal rejections against 0.7,
which overstated how far short those designs fell.

*Two cheap tests were launched to separate them, about $2 together, authorised by
the owner in chat on 2 October.* Test A runs the same campaign against a fragment
extended to 310–499, giving the severed disulfide its partner back, with **no filter
changed**: if it clears 0.6 on its own, the fragment was the problem and no
rejection rule needs touching. Test B keeps the committed 310–480 fragment and
lowers only `min_target_plddt_final` to 0.25, below the 0.35 observed, to see
whether designs then clear the remaining gates and to get a first rate.

**Test B was mis-specified and is useful anyway.** `min_target_plddt_final` governs
`Target_pLDDT` at the final stage only, so it relaxed a gate no trajectory reaches and
left `min_plddt_screen = 0.6` — the gate that actually rejected the validation
trajectory — untouched. The override was verified to have reached the filter list
(`Target_pLDDT >= 0.25` was printed) which made it look effective; what was not
checked was whether that filter was the one firing. The lesson is narrow and worth
keeping: the rejection message names the metric, and that name is what a diagnostic
should be built against, not a setting whose name resembles it.

Because the changed setting binds only at a stage nothing reached, Test B is in
substance **ten trajectories on the committed fragment under default settings**, which
is the small pilot that was wanted. It was left running for that reason rather than
restarted.

**RESULT, both diagnostics complete, 2 October 2026, 21:07 IST: eleven trajectories
against our own target, zero accepted designs, no candidate sequence.** Test A ran
one trajectory on the extended fragment for $0.48; Test B ran its full ten on the
committed fragment for $2.33, finishing on its trajectory budget rather than on its
spend ceiling. Total spent to date about $6.29 of the $30, leaving about $23.71.

*Where the eleven stopped, and on which gate.* Six at screen on `pLDDT` against its
0.6 floor, three at harden and two at anneal on `i_pTM` against its 0.5 floor.
Nothing reached the mutate or final stages, so **no structure was ever written** and
`3_Ranked/` is empty in both campaigns. Binder lengths were drawn across the whole
30–100 band — 31, 41, 41, 47, 60, 66, 68, 79, 88, 94 — so the failure is not confined
to one size.

*The bottleneck is interface confidence, and the margin is the useful number.*
`i_pTM` readings across every stage of every trajectory ran from 0.20 to 0.62. The
middle-stage floor is 0.5 and the final gate is 0.7, so **even the best trajectory
seen never reached what the final gate requires**, and four of the eleven died on the
middle floor. Complex confidence `pLDDT` was the lesser problem: it ran 0.38 to 0.86
against a 0.6 floor and most trajectories cleared it.

*The His-tag detargeting works, and this is now a rate rather than an anecdote.*
Thirty-two `i_pTM.HisTag` readings across the eleven trajectories ran from 0.01 to
0.08, every one far below the 0.4 ceiling `egfr_common.DETARGET_IPTM_CEILING` sets.
The binder consistently ignores the tag. That half of the design rule is working as
intended, on the only evidence we have.

*The card is settled and oversized rather than marginal.* Peak memory 9,037 MiB of
the L4's 23,034, so 39% of the card for one worker. Two workers would fit, which is
the available throughput lever and was never exercised because BindCraft2 holds
fan-out to the trajectory budget.

*Cost per trajectory, measured: $0.231.* Mean design time 1,040 s over ten
trajectories, taken from the trajectory table's own `Timing` column rather than
wall-clock divided by count, because wall-clock is distorted by worker count and by
the compilation in trajectory 1. The spread matters for sizing: an attempt that dies
at screen averages 487 s and one that reaches harden averages 1,792 s, so a campaign
that fails early is cheaper per attempt than one that fails late. An earlier estimate
in this log of about $72 for 150 trajectories came from Test A's single deep run and
was the worst case; the measured mean gives about **$37 for 150 trajectories**, and
the remaining budget buys about **103**.

**Both rates the pilot existed to measure come back undefined, and that is the
finding.** Candidates per card-hour is zero over eleven trajectories and about 3.6
card-hours — not a small rate but an undefined one. The charge-pair survival rate
cannot be computed at all, because `analysis/10` reads structures and none exist.
Zero of eleven is consistent by the rule of three with a true acceptance rate
anywhere from 0 to about 27%, so eleven attempts cannot pin it down; what it does
give is moderate evidence against PD-L1's 33%, which would produce eleven
consecutive failures about 1% of the time. **A main run sized from these rates
cannot be proposed, because the arithmetic has a zero in it.** What the numbers
support instead is that the obstacle is interface confidence against this epitope,
and that is where a decision is needed rather than more trajectories at the same
settings.

*First trajectory from each, which is one sample each and not a result.* The
extended fragment reached `pLDDT.EGFR_domain3` 0.82 at screen and 0.87 at refine,
clearing the untouched 0.6 gate; the committed fragment under the relaxed gate
reached 0.72 and 0.80. Both passed their screen and refine stages. Consistent with
the severed C470–C499 disulfide mattering, and far too thin to conclude from, given
the run-to-run spread recorded above — the committed fragment alone has now produced
0.35 and 0.72 at the same stage. Interface confidence moved in opposite directions,
rising 0.28 to 0.60 on the extended fragment and falling 0.36 to 0.20 on the
committed one, against a filter needing 0.70. Worth watching rather than believing. The
diagnostic fragment and config are in `data/structures/6aru_domain3_ext499.pdb` and
`design/configs/diagnostic/`, kept apart from the generated campaign files;
`design/configs/diagnostic/README.md` records what each answers. Over the residues
the two fragments share the atom records are byte-identical, so they differ only in
where they stop. **Neither test changes a committed decision.** The boundary is
still 310–480 and the filter floor is still BindCraft2's default; moving either is
the owner's call and would belong in `analysis/09` behind its own constant.

**Every attempt is now recorded, not only the ones that survive.**
`analysis/18_campaign_inventory.py` writes one row per attempt with the outcome at
both gates and the reason, because `analysis/10` reads `3_Ranked/` alone and so
never saw anything BindCraft2 discarded — which is most of what a campaign is paid
for. It imports step 10's rules rather than restating them, the way step 15 already
does. Outputs are `data/derived/18-trajectory-ledger.csv`,
`data/derived/18-file-inventory.csv` (every retained file with a SHA-256
fingerprint, so the record survives in tracked files even where the bulk structures
are not committed) and `results/findings/18-campaign-inventory.md`. Run on the
validation output it recorded the single attempt as rejected at screen, with its
100-round optimisation trace and its final state. It also discovered BindCraft2's
real schema, which had been guessed at: the decision lives in a `terminated` column
that is blank for a surviving attempt and otherwise names the stage that stopped it,
and each attempt's losses table is joined by the folder it sits in rather than by a
column.

**Never quote the 26-trajectory total as one rate — it pools configurations that
were testing different things.** Counted from each campaign's own trajectory
table, 3 October 2026:

| configuration | trajectories | reached the final stage |
|---|---|---|
| eight hotspots (validate, diagnose, floors, ext499 probe) | 20 | 0 |
| four hotspots, the tight4 set | 6 | 2 |

Pooled as "26 attempts, 0 accepted" the project looks stuck. Split, the narrowed
configuration is the only one that has ever reached the last gate, and it did so
on BindCraft2's **default** floors, not the lowered ones. Quote rates per
configuration with their denominator from here on.

*What the split does not prove.* These campaigns differ in more than hotspot
count — the floors and, for one of them, the fragment also moved — so 2 of 6
against 0 of 20 is not a controlled comparison. It is the reason the current
round holds everything except length fixed.

**Novelty, first check on real sequences — 3 October 2026. No significant
similarity to anything in Swiss-Prot.** MMseqs2 18-8cc5c, the ten 48-residue
sequences from the short band against Swiss-Prot release of the same date
(575,748 sequences), sensitivity 7.5.

The raw percent identities look alarming and are not: best hits run 31-42%
identity over the whole query. **Every one of them is statistical noise.** The
best e-value across 2,239 hits is 0.30 and none reaches 0.001. An e-value of
0.30 means a hit that good is expected by chance in a database this size; one of
7.4e+03 means thousands are. A 48-residue query finds 35-40% identity over
40-odd residues by chance routinely, and that is not homology. The organisers'
own wording is that a design is de novo if it hits nothing "with any homology",
which is a significance question rather than a raw-identity one.

*Limits, and they are large.* Swiss-Prot only: the organisers also search the
Protein Data Bank, two patent databases, THPdb, PLAbDab and Proteinbase, none of
which was searched here. Nothing structural was run, and per the novelty blog
the structural half is what decides between level 3 and level 4. And all ten
sequences are ProteinMPNN redesigns of a single backbone, so this is one
backbone checked, not ten.

*What it changes:* nothing yet, but it removes the specific worry that de novo
output would fail the gate for an unanticipated reason. Query and hits kept at
`data/derived/19b-novelty-query.fasta` and `19b-novelty-swissprot-hits.m8`.
**Owed: this was run by hand and should become a script before it is quoted
anywhere final.**

**Round one, the three numbers it existed to measure — 3 October 2026, nine
recorded trajectories per band.**

| | short, 30-60 | long, 60-100 |
|---|---|---|
| trajectories recorded | 9 | 9 |
| reached the end of the stage pipeline | **1** | **0** |
| candidates written | **10** | 0 |
| accepted by BindCraft2 | 0 | 0 |
| card time | $2.51 | $2.77 |
| cost per trajectory | **$0.279** | **$0.308** |
| cost per candidate | **$0.25** | undefined |

*Cost per accepted design is undefined in both bands, again, because the
numerator is zero.* Cost per candidate is defined for the short band and is the
figure that sizes round two, because a candidate is what our own screens can
read and BindCraft2's acceptance is not our criterion.

*The band difference is in completion, not in price.* Per trajectory the short
band is about 10% cheaper, which is less than the length ratio would suggest and
not on its own a reason to prefer it. What separates them is that the short band
carried one trajectory through to candidates and the long band carried none in
nine. On n=9 each that is suggestive and not settled: one completed trajectory
against zero is a difference of one event.

*Charge-pair survival, the longest-standing unmeasured quantity, now has a
denominator:* 0 of 10 candidates reach three correct pairs. Best is two, and
both designs carrying two are rejected for contacting 442. Six are submittable
at one pair or none. All ten come from one backbone, so this is one backbone's
answer rather than ten.

**Round one launched, 3 October 2026.** Two campaigns of ten trajectories each,
`batch_r1_short` (30-60 aa) and `batch_r1_long` (60-100 aa), on the tight4 anchors
and the committed 310-480 fragment, with the lowered floors and
`save_design_sequences` on. Ceiling $3.50 of card time each.

*Why that setting was added.* The tight4 run reached the final stage twice and left
no recoverable sequence: every trajectory folder holds a `_losses.csv` of scalar
metrics and nothing else, on disk and on the volume alike, because
`save_design_sequences` defaults to false. Checked before spending anything.

*What round one measures:* acceptance per band, charge-pair survival with a real
denominator, and cost per accepted design per band. Round two is sized from those.

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

*Status, 2 October 2026, about 14:50 IST.* The condition is met: the PD-L1 smoke run has
produced accepted designs, two by 14:46, against a target that is not ours. That is the
pipeline producing candidate sequences, which is what the rule asks for. The run has not
finished. The CPU-only work done since is independent of it and costs nothing. Meeting the
rule does not license pointing a campaign at EGFR before item 2 of Next actions is done.

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

Reordered on 2 October 2026 after the Slack intake. The smoke run and the CPU-only work are
independent and only one of them costs money.

1. **PD-L1 smoke run — running.** Launched 14:17 IST, detached. When it finishes, report the
   card, the wall-clock time and time per trajectory, the peak memory, the size of the
   returned complex, and the target numbering comparison. The residue count is what lets a
   PD-L1 figure be carried to EGFR. Peak memory is an upper bound, since JAX can reserve most
   of the card at start-up. The same run answers whether the target comes back at the input's
   numbering, by comparing the target chain of every returned complex against the input
   (`analysis/target_numbering.py`); `analysis/10` refuses to score any candidate that
   fails that. Then run `analysis/10 --parse-only` against its output: the verdicts mean
   nothing because PD-L1 is not our target, but it proves the script can read what the
   pipeline really writes. Check the direction `i_pDAE` is ranked in `!_Ranked.csv` against
   the correction under Settled.
2. **Validate both campaign files start up before a real run — DONE for the notag file,
   2 October 2026.** It starts, both targets load, `termini_accessible` resolves and the tag
   is accepted as a sequence off-target and avoided. $0.20. The plain file was not run
   because the notag one started and is the one to prefer; it stays the fallback.
   **A new blocker took its place:** the campaign is stopped at BindCraft2's
   `Target_pLDDT >= 0.6` gate by our own fragment scoring 0.35 (Pipeline status). Two cheap
   tests are in flight to decide whether that is the fragment or the predictor. **Item 4
   does not start until that is settled**, because at the measured behaviour a full campaign
   would accept nothing.
3. **Owner decisions that shape the campaign file:** how to handle sugar chains near the
   anchors, and whether to design against human and mouse together (both under Open); which
   anchor cluster; and the molecule category. Geometry points to a minibinder of 40–100
   residues; a microbinder is possible only on the 15.6 Å triad. Show that the cluster
   choice was made deliberately, because the write-up is stronger for it.
4. **Generate candidates** against the confirmed patch, from `egfr-domain3-h370-notag.json`
   if item 2 shows it starts and from the plain file if not. Designs generated without the
   tag cannot be screened for it afterwards. **Resize the campaign to the budget first:** the
   file asks for up to 2,000 trajectories and 200 accepted designs, which at the measured
   rate is far beyond $30 (Pipeline status). Choose the counts from the card time you are
   willing to spend, and set `max_trajectories` so the run stops by itself.
5. **Test whether the fragment holds its shape,** which is what the cut risks and nothing has
   measured. Predict the 171-residue fragment with the design run's own predictor and
   compare it with 6ARU across the anchors. Step 11 tried with ESMFold and got about 5 Å for
   all three boundaries, the predictor's noise floor. The ratio of run times for 171
   against 621 residues is also unmeasured: time one trajectory on each.
6. **Run the three screens on real output:** `analysis/10_charge_pair_filter.py` (charge
   pairs, then `i_pDAE` as a tie-break with the higher reading first),
   `analysis/13_full_receptor_clash.py` (would each binder fit in the intact receptor) and
   `analysis/15_histag_counterscreen.py` (`--candidates`). All three have only ever seen
   constructed inputs.
7. **Novelty, spread, submit.** Check novelty with MMseqs2 against the databases the
   organisers listed, and test against the real checker, which scores automatically on
   upload. Read `adaptyvbio.com/blog/novelty` first. Spread the twenty slots across anchor
   subsets and size categories. Upload an early safe batch and replace it, one submission per
   24 hours. Deadline is now 6 October 23:59 AoE. Write the reasoning, because selection reads
   the documentation.
8. **Rewrite `docs/explainers/04-design-pipeline.md` after the first real pipeline run,** with
   what happened in place of what was planned. Its figures are external, read on 1 October,
   and any that still matter move into a findings file produced by a script in `analysis/`.
9. **Resolve the H418 question — a fallback question, not the critical path.** It decides
   whether the 415–466 fallback offers a four-anchor cluster with no target histidine or a
   five-anchor one with a target histidine and no cetuximab overlap (D416, H418, E421, E424,
   E455, span 24.0 Å). Under the tethered form 1NQL's reading of 0.200 is the relevant one
   and 6ARU's 0.032 may be an artefact of the antibody. Resolve by checking further
   structures and side-chain rotamers. Do not reopen the epitope choice.
10. **Burial of charge pairs:** revisit when real complexes exist (Open). No ranking term
    until it can be tested.
