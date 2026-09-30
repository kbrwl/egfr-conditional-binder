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
`competition-brief.md`.

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
with a definite answer. Concentrate effort on objectives 2 and 3, which are also
the higher-weighted ones. The organisers confirm this reading explicitly.

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

**Affinity target is deliberately marginal, not maximal.** Objective 3 requires
*no detectable binding* at pH 7.4 — a threshold, not a ratio. A design binding
at 10 nM / 200 nM (20-fold selective) fails, because 200 nM is plainly
detectable. A design binding at 2 uM / undetectable passes. Organisers confirm a
weak but clearly pH-sensitive binder may outrank a strong non-selective one.
Standard design pipelines maximise affinity by default, which produces exactly
the failing case. **Do not accept the pipeline's built-in objective.** Target a
baseline weak enough that the pH 7.4 state falls under the assay detection floor,
then build the largest switch possible on top of it.

---

## Open

**Which of the eight anchors are actually on the surface.**
Acidic anchors D416, E421, E424, E455, D458, D460, plus target histidines H418
and H433 (which anchor an acidic residue on our side). Some are certainly
buried. Needs solvent accessibility calculation on the 3D structure.

**Whether the usable anchors form a spatial cluster.**
They must sit within roughly 25 Å for one binder to reach them. If they are
scattered across faces, fall back to run 394–411 or 331–347.

**Cetuximab's true contact set.** See "Unverified" below.

**Which molecule category to target.** Microbinder <40 aa, minibinder 40–100 aa,
large binder >100 aa, nanobody or antibody. Winners are announced per category,
so a less crowded category may be worth choosing deliberately. Not yet chosen.

**Is 415-466 accessible in the conformation being assayed?** EGFR's extracellular
region switches between a tethered (closed) and an extended (open) form, and
domain III is partly occluded by domain II in the tethered state. The assay uses
the full extracellular region, so the dominant conformation in that buffer
decides what a binder can physically reach. If our block is buried in the
tethered form, a correct design measures as nothing. Unknown. Check against 6ARU
and the tethered structures before committing to the epitope.

**Do objectives 2 and 3 compete for the same surface?** Both pH selectivity and
mouse cross-reactivity constrain the same interface residues. Not yet checked
whether the charge-pairing positions we want are all inside the conserved block.

**Compute environment.** Nothing stood up yet. Colab / Modal / RunPod all viable.
This is the main schedule risk — dependency and CUDA problems eat days.

**Design pipeline choice.** BindCraft is the accessible one-command option.
RFdiffusion + ProteinMPNN gives more control at higher setup cost. Not chosen.

---

## Unverified — do not build on these

- **Cetuximab contact residues.** The claim that Q390R, E412D, R414W and K467R
  explain cetuximab's failure on mouse EGFR was stated from memory. It is
  plausible and consistent with the alignment, but it has not been computed.
  Resolve by measuring contacts in structure 6ARU.
  **6ARU is the structure of record for this project** — the EGFR extracellular
  region bound to a cetuximab Fab mutant, and the entry the competition page
  itself references. Preferred over 1YY9 (cited in the glossary and in
  `alignment-findings.md`) because 6ARU contains the whole extracellular region,
  which is what the assay uses, rather than domain III alone.
- **Domain III boundaries.** 310–480 is the working definition used throughout.
  Confirm against the official annotation.

---

## Ruled out

Nothing yet.

---

## Next actions, in order

1. Download structure 6ARU — EGFR extracellular region with a cetuximab Fab
   mutant, the reference the competition itself uses. Note chain A is the receptor.
2. Compute solvent accessibility for every residue in 415–466. Keep the exposed ones.
3. Compute cetuximab's real contact set at a 4.5 Å cutoff.
4. Check the surviving anchors cluster within ~25 Å.
5. Stand up a GPU environment and get a binder-design pipeline running end to end
   on any target, however bad the output, before attempting the real one.
6. Generate several hundred candidates against the confirmed patch.
7. Charge-pair scan. For each candidate, look at which target residue each
   binder contact position faces, then apply the positional rule: histidine
   opposite acidic targets, acidic opposite H418/H433. Keep versions where at
   least three correct pairs form and the fold still holds. Explicitly reject
   any candidate with a histidine facing H418 or H433.
8. Novelty check, rank, submit top 20 with written reasoning.

Steps 1–4 are cheap and can run anywhere. Step 5 is the schedule risk and should
start in parallel, not after.
