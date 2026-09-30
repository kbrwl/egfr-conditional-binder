# Rules reference

The competition's own facts: dates, rules, formats, limits, tracks and links.
Recorded so they do not have to be looked up again.

This file holds only what the organisers state. Our readings, inferences,
positions and strategy are in `decisions-log.md`. If a line here would be
embarrassing were it wrong, it is an inference and belongs in that file instead.

Captured 30 September 2026. Anthropic x Adaptyv Protein Design Competition,
Challenge 1 of 5.

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

Challenge 1 deadline: **4 October 2026, 23:59 Anywhere on Earth**.

"Anywhere on Earth" is a defined deadline convention: the deadline has not passed
until it has passed in every time zone, which is UTC−12, so about midday UTC on
5 October.

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

- Which assay is used, and therefore what counts as "no detectable binding".
- Which conformation of the extracellular region predominates in the assay buffer.
- The exact residue boundaries of domain III.
