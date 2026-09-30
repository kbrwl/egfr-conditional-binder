# Competition brief — the official facts

Reference copy of the rules, dates, formats and links, so they never have to be
re-fetched or half-remembered. **Official statements only.** Our reasoning,
decisions and interpretations live in `decisions-log.md`; this file is the
external ground truth those decisions are built on.

Captured 30 September 2026.

---

## Links

| What | Where |
|---|---|
| Challenge 1 (EGFR) page | https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr |
| Competition main page | https://proteinbase.com/competitions/anthropic-adaptyv-2026 |
| Novelty rules (Adaptyv blog) | https://www.adaptyvbio.com/blog/novelty |
| Anthropic folding/design models | https://github.com/anthropics/uplifting-biomolecular-modeling |
| Proteinbase Slack | https://join.slack.com/t/proteinbase/shared_invite/zt-3evw8fs9z-tU9ItWVvw4ySctUuPvIhLQ |

---

## Dates

**Challenge 1 deadline: 4 October 2026, 23:59 Anywhere on Earth.**

"Anywhere on Earth" means the deadline has not passed until it has passed in
every time zone on the planet — effectively UTC−12, so roughly midday UTC on
5 October. Do not rely on that margin.

Four further challenges follow weekly until **1 November 2026**. Consequence for
how we build: the analysis pipeline should retarget by swapping inputs, not by
being rewritten. Five challenges in five weeks rewards reusable tooling and
punishes one-off scripts.

---

## The design problem

Design protein binders against the **extracellular region of human EGFR**
(epidermal growth factor receptor), UniProt P00533 residues 25–645, 621 amino
acids. The organisers recommend targeting **domain III**, where the approved
antibodies cetuximab and panitumumab bind.

Three objectives, **ranked in this order of importance by the organisers**:

1. **pH selectivity.** Bind human EGFR at pH 6.5 with **no detectable binding**
   at pH 7.4. *This is a threshold, not a ratio.*
2. **Mouse cross-reactivity.** The same sequence must also bind mouse EGFR.
3. **Affinity** against human EGFR.

The organisers state explicitly that **a weak but clearly pH-sensitive binder may
outrank a high-affinity binder that is not pH-sensitive.**

That sentence is the single most consequential line in the brief, because it
inverts the default objective of every standard design pipeline. See
`decisions-log.md` → "Affinity target is deliberately marginal, not maximal".

Why pH 6.5 versus 7.4: blood and healthy tissue sit at pH 7.4, while solid
tumours drift down to roughly 6.5. A binder that only grips in the acidic
environment hits the tumour and releases healthy tissue.

---

## Submission format

A CSV, **ranked best-first**, with columns:

| Column | Contents |
|---|---|
| `name` | your identifier for the design |
| `sequence` | the amino acid sequence |
| `molecule_class` | one of `protein`, `nanobody`, `scfv`, `fab_kappa`, `fab_lambda` |

- Fabs are submitted as a single field, `{VH}:{VL}` — heavy chain, colon, light
  chain.
- **Length 10–250 amino acids.**
- Row order is the ranking, so ordering is a judgement we are making on the
  record.

Designs are submitted **alongside written documentation**, which is read during
selection. See "Selection" below.

---

## Molecule categories

Judged separately, so the category is a deliberate choice rather than a
by-product of design:

| Category | Length |
|---|---|
| Microbinder | under 40 aa |
| Minibinder | 40–100 aa |
| Large binder | over 100 aa |
| Nanobody | — |
| Antibody | — |

---

## Hard rules

- **De novo.** Designed from scratch. You may **not** take an existing binder and
  modify it.
- **Zero-shot.** No experimental feedback loop — designs are submitted without
  having been tested and refined.
- **Adequate sequence and structural diversity from known proteins** is required.
  Novelty is a judged criterion, not just a filter. See the Adaptyv novelty blog.
- **Embedded instructions or prompt injection in a submission is grounds for
  disqualification.** Selection is performed partly by a Claude workflow reading
  submitted documentation; attempting to influence it is disqualifying. Our
  documentation states facts and reasoning and never addresses the reader as an
  instruction.

---

## Tracks, and our position in them

We are in **Track 3 (open track)**. Applications for Tracks 1 and 2 closed
24 September 2026, before we found the competition. This was not a choice.

| | Track 1 | Track 2 | Track 3 |
|---|---|---|---|
| Support | Anthropic credits | Anthropic credits | **self-supported** |
| Designs allowed | 20–40 | up to 20 | **up to 20** |
| Testing allocation | ~15 reserved per challenge | none guaranteed | **none guaranteed** |
| Screening slots | 50% of ~1500 | ~375 pooled | **~375 pooled** |

Track 3 means our own tools and our own compute, with no Anthropic credits.

**Testing is not guaranteed.** Roughly 1500 designs are screened per challenge,
split 50% / 25% / 25%. Track 1 teams hold a reserved allocation. Tracks 2 and 3
compete for pooled slots, so our designs may never be physically made.

### Selection

For Tracks 2 and 3, all submissions are run through a **Claude workflow** that
weighs:

- predicted design quality
- **design novelty**
- **method novelty**

reading **the documentation submitted alongside the sequences**. The organisers
state selection will not rely on a single in-silico metric.

**Therefore the write-up is part of the submission, not paperwork.** In a pooled
track with no reserved slot, a clear written argument for the epitope choice and
the design rule is the main available edge.

---

## Prizes

Non-cash. Details unannounced at time of capture.

**UNVERIFIED — whether winners are judged within tracks or across all tracks.**
Screening quotas are explicitly per-track. The winner categories described by the
organisers (highest affinity, most cross-reactive, most pH-sensitive) make no
mention of tracks. Best available reading is that Track 3's disadvantage is
*getting screened*, not the judging afterwards — meaning a design that does get
tested competes directly with professional lab submissions. **This is not stated
anywhere and remains an inference.** Worth asking in the Proteinbase Slack.

---

## Freely available to all tracks

- Inference-optimised folding and design models at
  `github.com/anthropics/uplifting-biomolecular-modeling`, with a simplified
  binder-design prompt in the accompanying technical report.
- BindCraft2, released September 2026, for academic and industry use.

Claude and Modal credits go to selected Track 1 and 2 teams only.

---

## What the organisers do *not* specify

Recorded because these are the gaps we had to make our own decisions about, and
the write-up should be honest that they are decisions rather than requirements:

- Which assay, and therefore **what the detection floor actually is**. Objective 1
  is defined against "detectable", so the threshold we are designing to is a
  number we do not have. We target a wide margin rather than a computed one.
- Which **conformation** of the extracellular region predominates in the assay
  buffer (tethered vs extended), which determines what is physically reachable.
- The exact **domain III boundaries**. We use 310–480 throughout as a working
  definition.
