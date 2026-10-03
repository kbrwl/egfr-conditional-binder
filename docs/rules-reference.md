# Rules reference

The competition's own facts: dates, rules, formats, limits, tracks and links.
Recorded so they do not have to be looked up again.

This file holds only what the organisers state. Our readings, inferences,
positions and strategy are in `decisions-log.md`. If a line here would be
embarrassing were it wrong, it is an inference and belongs in that file instead.

Captured 30 September 2026. Anthropic x Adaptyv Protein Design Competition,
Challenge 1 of 5.

Extended 2 October 2026 with what the organisers said in the Proteinbase Slack,
which the competition page does not carry. Those additions are marked *(Slack)*,
and the person, date and wording behind each are in `competition-qa-log.md`.
Where a Slack statement and the competition page disagree, the page wins.

---

## Links

| What | Where |
|---|---|
| Challenge 1 (EGFR) page | https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr |
| Competition main page | https://proteinbase.com/competitions/anthropic-adaptyv-2026 |
| Novelty rules (Adaptyv blog) | https://www.adaptyvbio.com/blog/novelty |
| Anthropic folding and design models | https://github.com/anthropics/uplifting-biomolecular-modeling |
| Proteinbase Slack | https://join.slack.com/t/proteinbase/shared_invite/zt-3evw8fs9z-tU9ItWVvw4ySctUuPvIhLQ |

---

## Dates

Challenge 1 deadline: **6 October 2026, 23:59 Anywhere on Earth** *(Slack, official)*.

"Anywhere on Earth" is a defined deadline convention: the deadline has not passed
until it has passed in every time zone, which is UTC−12, so about midday UTC on
7 October.

**Settled 3 October 2026.** The original date was 4 October. An earlier Slack message
(Tudor, 30 September) was read as a one-day extension to 5 October. The official
announcement in the Proteinbase Slack supersedes that: challenge 1 is extended to
**6 October 23:59 AoE** because some teams received Modal credits late. The competition
page was still to be updated to match when this was recorded; if the page and this file
disagree, re-read the page. Challenge 2 starts on time, so the two overlap by about two days.

Four further challenges follow weekly until 1 November 2026.

---

## The target and the objectives

Design protein binders against the extracellular region of human EGFR (epidermal
growth factor receptor), UniProt P00533 residues 25–645, 621 amino acids. The
organisers recommend targeting domain III, where the approved antibodies cetuximab
and panitumumab bind.

Three objectives, in the organisers' stated order of importance:

1. pH selectivity: bind human EGFR at pH 6.5 with no detectable binding at pH 7.4.
2. Mouse cross-reactivity: the same sequence must also bind mouse EGFR.
3. Affinity against human EGFR.

The organisers state that a weak but clearly pH-sensitive binder may outrank a
high-affinity binder that is not pH-sensitive.

---

## The assay and the construct *(Slack)*

None of this is on the competition page. All of it was stated by the organisers in
Slack between 28 September and 2 October 2026; see `competition-qa-log.md` for who
said what and when.

**Measurement**

- The main assay is surface plasmon resonance (SPR). Bio-layer interferometry (BLI)
  is used only to cross-check hits.
- **The designs are immobilised and the target is flowed over them.** The target is
  captured through a C-terminal twin-Strep tag, described as the default and
  changeable.
- Top analyte concentration is 1000 nM. The reportable affinity range is about
  0.1 nM to 10 µM.
- A design counts as binding at a given pH if its trace can be fitted to give a K_D,
  or, where no fit is possible, if the association signal rises more than 300%
  above the negative control. Failing both is "no detectable binding".
- Designs that bind at both pH values with a large K_D shift still qualify. Designs
  with no binding at pH 7.4 and high affinity at pH 6.5 rank higher.
- Mouse cross-reactivity means a mouse-to-human K_D ratio of about 1.

**Buffer**

- At pH 7.4: 10 mM HEPES, 150 mM NaCl, 0.2% Tween-20, 3 mM EDTA.
- At pH 6.5: MES replaces HEPES.
- Total ionic strength is matched between the two conditions with NaCl, at roughly
  170 mM.

**The target**

- The **tethered** (closed, inactive) conformation.
- Expressed in HEK293 cells and **glycosylated**.
- Both the human and the mouse protein carry a C-terminal His tag. The tags are
  unlikely to be cleaved.
- The organisers will weight in silico scoring more heavily where a design might be
  binding a tag, and will run follow-up neutralisation assays.

**The designs**

- Built with a C-terminal tail of linker, GFP11, linker, twin-Strep tag, and
  immobilised at the C-terminus.
- **Cyclic peptides are not supported.** Linear chains only.

**Novelty**

- The gate is level 3 of 4, scored automatically on upload.
- The organisers' check is MMseqs2 against SwissProt, the Protein Data Bank, the
  USPTO patent database, the EBI patent database, THPdb, PLAbDab and Proteinbase.
- One organiser also described the requirement as the sequence being under 30%
  similar to anything existing. Whether this is the same threshold as level 3 is not
  established. **It is not the same**: read against the blog on 3 October 2026,
  under 30% on sequence gives level 4 only when structural similarity is also low.

*The levels, from the blog.* Two measurements combine: sequence similarity by
MMseqs2, and structural similarity by predicting the fold and matching it with
FoldSeek and TM-align. "Moderate" structural similarity means over 70% of the
sequence is covered by domains matching a known structure at TM-score 0.5 or
better; "high" means the same coverage at 0.8 or better.

| level | what it takes |
|---|---|
| 4, de novo | sequence 30% or under **and** structural below moderate |
| 3, partly novel | exactly **one** of: sequence over 30%, or moderate structural |
| 2, familiar | sequence over 70%; or high structural; or sequence over 30% **and** moderate structural |
| 1, essentially known | sequence over 70% **and** at least moderate structural |

*What follows for us.* The gate is level 3, so levels 3 and 4 both pass. A
designed binder under 30% on sequence that still folds into a common shape, such
as a helical bundle, scores level 3 and passes. Failing the gate takes landing on
level 2, which from a de novo design means sequence similarity over 30% together
with moderate structural similarity. So the thing to watch is sequence similarity,
not the fold. This also removes the worry that a conventional-looking fold alone
could disqualify a design.

**Selection**

- Mainly method novelty, design diversity and a couple of in silico confidence
  metrics. In silico scoring is described as the smaller component, and worse
  confidence scores on a rigid structure or an obstructed epitope will not be
  penalised much, if at all.
- The selection prompts will be shared after the competition.

**Submission mechanics**

- Up to 20 designs for Track 3 (40 for Track 1), each chain 10 to 250 amino acids.
- Track 3 needs no confirmation email: designs can be submitted directly.

**CONFLICT, unresolved as of 3 October 2026 — how many times we may submit.**

| source | says |
|---|---|
| the competition page, read 3 October | "One submission per participant per challenge" |
| a participant quoting the submission page in Slack, 1 October | "one submission every 24 hours" |

Our standing rule is that the page beats Slack, and the page is the stricter of
the two. **Plan on having exactly one submission.** The wording may only mean one
CSV file rather than several, and may still permit replacing it — but that reading
is not established, and acting on it risks spending the entire entry on a draft.

This matters because a strategy of "submit something safe early and improve it
later" is only sound under the 24-hour reading. Under the page's reading it would
end the entry. Resolve by asking the organisers, or by checking whether the
submission portal itself offers a replace option, before submitting anything.

---

## Submission format

A CSV file, ranked best-first, with these columns:

| Column | Contents |
|---|---|
| `name` | your identifier for the design |
| `sequence` | the amino acid sequence |
| `molecule_class` | one of `protein`, `nanobody`, `scfv`, `fab_kappa`, `fab_lambda` |

- Fabs are submitted in one field as `{VH}:{VL}` — heavy chain, colon, light chain.
  A Fab is the gripping arm of an antibody; VH and VL are its two chains, heavy and
  light.
- Length limit: 10–250 amino acids.
- Row order is the ranking.
- Designs are submitted alongside written documentation, which is read during
  selection.

---

## Molecule categories

Judged separately from one another:

| Category | Length |
|---|---|
| Microbinder | under 40 amino acids |
| Minibinder | 40–100 |
| Large binder | over 100 |
| Nanobody | — |
| Antibody | — |

---

## Rules

- **De novo.** Designed from scratch. An existing binder may not be taken and
  modified.
- **Zero-shot.** Submitted without having been tested in a lab and refined.
- **Diversity.** Adequate sequence and structural diversity from known proteins is
  required. Novelty is a judged criterion.
- **No embedded instructions.** Embedded instructions or prompt injection in a
  submission is grounds for disqualification.

---

## Tracks

Applications for Tracks 1 and 2 closed 24 September 2026. We are in Track 3, the
open track.

| | Track 1 | Track 2 | Track 3 |
|---|---|---|---|
| Support | Anthropic credits | Anthropic credits | self-supported |
| Designs allowed | 20–40 | up to 20 | up to 20 |
| Testing allocation | about 15 reserved per challenge | none guaranteed | none guaranteed |
| Screening slots | 50% of about 1500 | about 375 pooled | about 375 pooled |

About 1500 designs are screened per challenge, split 50% / 25% / 25%. Track 1 teams
hold a reserved allocation. Tracks 2 and 3 compete for pooled slots.

### Selection

For Tracks 2 and 3, all submissions are run through a Claude workflow that weighs
predicted design quality, design novelty and method novelty, reading the
documentation submitted alongside the sequences. The organisers state that
selection will not rely on a single computational metric.

---

## Prizes

Non-cash. Details unannounced at the time of capture.

The winner categories the organisers describe are highest affinity, most
cross-reactive and most pH-sensitive. They do not say whether these are judged
within tracks or across all of them.

---

## Available to all tracks

- Inference-optimised folding and design models at
  `github.com/anthropics/uplifting-biomolecular-modeling`, with a simplified
  binder-design prompt in the accompanying technical report.
- BindCraft2, released September 2026, for academic and industry use.

Claude and Modal credits go to selected Track 1 and 2 teams only.

---

## Not specified by the organisers

Recorded because these are gaps in the official information, not our opinions about
them. What we decided to do about each is in `decisions-log.md`.

- The exact residue boundaries of domain III.
- The exact mouse EGFR construct: residue range and vendor. Asked twice in Slack on
  1 October and not answered.
- Whether extra metric columns are tolerated in the submission CSV beyond the three
  in the template. Asked on 1 October and not answered.
- **Whether cynomolgus monkey cross-reactivity is in scope for this challenge.** A
  pre-launch message from Amir at Anthropic on 28 September described the challenge
  as wanting binding "to mouse and cyno". The competition page, read on 2 October,
  names human and mouse only. One of the two is out of date.
- Whether the competition page will be updated to carry the Slack statements above.
  Tudor said it would be updated to name the tethered form; on 2 October it was not.

Answered since capture and moved up into the assay section: which assay is used,
what counts as "no detectable binding", and which conformation is tested.
