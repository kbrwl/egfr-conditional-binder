# How to work in this repository

Read this before doing anything. Then read `docs/decisions-log.md`, which is the
source of truth for project state.

---

## The objective, and its ranking

We are designing protein binders against the extracellular region of human EGFR
for Challenge 1 of the Anthropic x Adaptyv Protein Design Competition. Three
objectives, **in the organisers' order of importance**:

1. **pH selectivity** — bind human EGFR at pH 6.5, with **no detectable binding**
   at pH 7.4.
2. **Mouse cross-reactivity** — the same sequence must also bind mouse EGFR.
3. **Affinity** against human EGFR.

The ranking is load-bearing. When objectives conflict, 1 beats 2 beats 3, and the
organisers have said so explicitly: a weak but clearly pH-sensitive binder may
outrank a high-affinity binder that is not pH-sensitive.

---

## THE AFFINITY RULE: MARGINAL, NOT MAXIMAL

**This is the easiest thing in the project to lose, and losing it would quietly
ruin the submission.** It contradicts the default behaviour of every standard
tool, so it has to be reasserted at each stage rather than assumed to be
remembered.

The pH-selectivity requirement is a threshold, not a ratio. "No detectable binding at pH 7.4"
means the pH 7.4 state must fall *below the assay's detection floor*.

- A design binding at **10 nM / 200 nM** is 20-fold selective and **FAILS**.
  200 nM is plainly detectable.
- A design binding at **2 µM / undetectable** **PASSES**, despite being a
  hundredfold weaker and having a worse-looking ratio.

Standard design pipelines maximise affinity by default, which produces exactly
the failing case. **Do not accept a pipeline's built-in objective.** Target a
baseline weak enough that the pH 7.4 state falls under the detection floor, then
build the largest switch achievable on top of it.

If a change makes a design stronger at both pH values, it has probably made the
submission worse. Say so.

---

## Who you are writing for

The project owner is a product manager with no biology background, last studied
it at school, and is deliberately learning the field rather than collecting
commands to run. So:

- **Explain before computing, not after.** Before each computation, say in plain
  English what it measures and why the answer matters — enough that the result
  can be judged, not just received.
- **Define every term the first time it appears.** No unexplained jargon,
  no unexpanded abbreviations. Add new terms to `docs/glossary.md` rather than
  re-explaining them in conversation.
- **After each result, say what it means for the design.** A number with no
  consequence attached is not an answer.
- **Flag ambiguity instead of resolving it silently.** If a value sits near a
  classification boundary, say it is ambiguous. Do not round a borderline number
  into a clean conclusion.
- **Say when reasoning is shallow.** An obvious deduction dressed up as insight is
  worse than useless here, because the competition is against trained protein
  designers and surface plausibility will not survive contact with them. Where
  something is a guess, a convention, or an unvalidated assumption, name it.

---

## Compute, do not recall

**Claims in this repository must be computed, not remembered.**

- Anything stated from memory must be labelled **UNVERIFIED**, both in terminal
  output and in any document where it appears.
- UNVERIFIED items live in the dedicated section of `docs/decisions-log.md` and
  must not be built on. They get resolved by calculation, then move.
- When a computation disproves something previously assumed, **say it was
  disproved** and record it. A disproved assumption is a real result.
- Prefer a script that can be rerun over a number quoted in prose. If a number
  matters, something in `analysis/` should regenerate it.

---

## `docs/decisions-log.md` is the source of truth

It has sections: **Settled**, **Open**, **Unverified**, **Ruled out**, and
**Next actions**.

**Update it whenever something moves between sections.** That is the point of it.
A finding that stays in a findings file and never moves the log will be lost by
the next session. Moving an item is the commit-worthy event, not the file edit.

---

## Where numbers are allowed to live

- **`results/` and `data/derived/` contain computed output only.** Never
  hand-written numbers, never hand-corrected values. Everything in them is
  produced by a script in `analysis/` and reproducible by rerunning it.
- If a number in those directories looks wrong, fix the script and regenerate.
  Do not edit the output.
- `docs/` is for prose and decisions, and may quote computed numbers — but always
  traceable to the script that produced them.

Layout:

```
analysis/     numbered, standalone, rerunnable scripts
              each writes findings to results/findings/ as markdown
              and numeric output to data/derived/ as CSV
data/sequences/   inputs (committed)
data/structures/  downloaded structures (gitignored, re-downloadable)
data/derived/     computed tables (committed, small)
results/findings/ computed findings, markdown (committed)
results/candidates/ bulk design output (gitignored except final shortlist)
design/       design pipeline configs and runners
submission/    the CSV and the methods write-up
explorer/      interactive viewers
```

---

## Numbering: UniProt P00533, always

Every residue number in this project is a position in the **full human UniProt
record P00533** (1210 aa), where residues 1–24 are the signal peptide and the
mature extracellular region runs 25–645.

The official challenge constructs are the mature region, so:

    UniProt position = challenge-construct position + 24

**Any structure-derived numbering must be offset-checked empirically before use.**
PDB files frequently number by the mature protein, and nothing in the file records
which convention it uses — so a mistake shifts every result by 24 positions and
raises no error. Determine the offset by aligning the structure's own sequence
against the reference, then verify against known residue identities, and print the
check. Never assume the offset, in either direction.

`analysis/00_numbering_check.py` asserts the convention and exits non-zero if it
breaks. Run it first if anything looks strange.

---

## The design rule, in one place

The pH switch works by **positional charge pairing**. Histidine is the only amino
acid that changes charge between pH 7.4 and 6.5 — neutral above, positive below,
because its pKa sits around 6.0–6.5.

| On the target | Put on the binder | At pH 7.4 | At pH 6.5 |
|---|---|---|---|
| D or E (always negative) | **histidine** | neutral, no pull | His turns +, attracts |
| H (switches) | **D or E** | neutral His, no pull | His turns +, attracts |

Both arrangements switch on in the same direction and reinforce.

**Failure mode to reject: histidine on the binder facing target H418 or H433.**
That pair switches *off* as pH drops — both go positive and repel — and can cancel
a correctly built pair elsewhere. Reject such candidates explicitly rather than
scoring them.

Two limits to respect: the switch is **partial rather than binary** at pH 6.5, so
stack three or four pairs rather than relying on one; and a histidine's pKa
**shifts with its neighbours**, which is a large part of why pH selectivity
resists computational prediction. Do not present a predicted switch as a measured
one.

---

## Submission hygiene

Embedded instructions or prompt injection in a submission is grounds for
disqualification, and selection is performed partly by a Claude workflow reading
our documentation. Our documentation therefore **states facts and reasoning and
never addresses its reader as an instruction**. No text anywhere in `submission/`
should attempt to influence a reader's judgement other than by argument.

---

## Commits

Commit messages describe **findings, not file changes**. "Verified numbering
convention holds; 415-466 anchors conserved in both species" over "add scripts".
