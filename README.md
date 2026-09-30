# A pH-conditional EGFR binder, designed by positional charge pairing

Submission work for **Challenge 1 of the Anthropic x Adaptyv Protein Design
Competition** (Track 3, open track). Deadline 4 October 2026.

This repository contains the reasoning, the computations and the code behind a
set of designed proteins intended to grip a cancer target **only in the acidic
environment of a tumour**, and to let go in healthy tissue.

It is written to be readable by someone who knows neither protein design nor the
author. Every term is defined in [docs/glossary.md](docs/glossary.md), and every
number in `results/` and `data/derived/` was computed by a script in
[analysis/](analysis/) rather than quoted from memory.

---

## The problem, in plain English

**EGFR** (epidermal growth factor receptor) is a protein that sits threaded
through the outer membrane of a cell, with one part outside and one part inside.
It works like a doorbell: when a signalling molecule lands on the outside part,
the inside part switches on and tells the cell to divide. In many cancers this
doorbell is stuck down, so the divide signal never stops. Blocking it slows the
tumour, and two antibody drugs that do exactly that — cetuximab and panitumumab —
are already approved and in clinical use.

They have a problem. EGFR is not only on tumours. It is also on healthy skin and
gut cells, and a drug that cannot tell the difference attacks those too.
Cetuximab causes a severe acne-like rash, sometimes bad enough that patients stop
treatment.

**So: can a binder be built that only works inside a tumour?**

There is a physical difference to exploit. Healthy tissue and blood sit at
**pH 7.4**. Solid tumours drift down to around **pH 6.5**, because they burn
sugar inefficiently, dump acid, and drain poorly. That gap of 0.9 pH units is the
only handle available, and it is the whole basis of this design.

The competition asks for a binder that satisfies three objectives, **in this
order of importance**:

1. **pH selectivity** — binds human EGFR at pH 6.5, with **no detectable binding**
   at pH 7.4.
2. **Mouse cross-reactivity** — the same sequence must also bind mouse EGFR
   (so the design can be tested in mice before humans).
3. **Affinity** — how tightly it grips human EGFR.

---

## The counter-intuitive part: stronger is worse

Objective 1 is a **threshold, not a ratio**. "No detectable binding at pH 7.4"
means the pH 7.4 state has to fall below what the measuring instrument can see at
all. It is not a claim about the gap between two numbers.

The consequence inverts normal practice:

| Design | Ratio | Verdict |
|---|---|---|
| binds at 10 nM (pH 6.5) / 200 nM (pH 7.4) | 20-fold selective | **FAILS** — 200 nM is plainly detectable |
| binds at 2 µM (pH 6.5) / nothing measurable | worse on paper | **PASSES** |

The second design is a hundred times weaker and wins, because the requirement is
about crossing a floor rather than opening a gap.

Every standard protein design pipeline maximises binding strength by default,
which produces precisely the failing case. **The affinity target here is
deliberately marginal, not maximal** — weak enough that pH 7.4 disappears under
the detection floor, with the largest achievable switch built on top. This is
recorded prominently in [CLAUDE.md](CLAUDE.md) because it is the easiest thing in
the project to lose and it would ruin the submission quietly, not loudly.

---

## The approach: positional charge pairing

### Why histidine

Of the twenty amino acids, **histidine is the only one that changes electrical
charge between pH 7.4 and pH 6.5.** It is neutral above roughly pH 6.5 and
positive below it, because its switching point (its **pKa**) happens to sit at
about 6.0–6.5 — right in the window between healthy tissue and tumour. Nothing
else in the toolkit switches in that range. So any pH switch has to be built out
of histidine.

### Why *positional*

A protein interface is a **patterned** surface, not a uniformly charged one. Each
position on the binder faces one specific position on the target, and
electrostatic attraction falls away quickly with distance — a charge twenty
ångströms away barely affects a given pair. So contact pairs behave more or less
independently, and the design rule is per-position rather than global. "Sprinkle
histidines across the interface and score it" is not the same thing and is not
what we do.

The rule has **two** halves:

| On the target | Put on the binder | At pH 7.4 | At pH 6.5 |
|---|---|---|---|
| D or E (acidic, always negative) | **histidine** | neutral, no pull | histidine turns +, attracts |
| **H (histidine — switches)** | **D or E** | neutral histidine, no pull | histidine turns +, attracts |

Both arrangements switch on in the same direction as pH falls, so they reinforce
each other.

**The failure mode we reject explicitly:** a histidine on the binder placed
opposite a histidine on the target. That pair switches *off* as pH drops — both
become positive and repel — and can cancel out a correctly built pair elsewhere.
Candidates with this arrangement are rejected, not scored.

### What we think is novel about this

The second row of that table. Standard pipelines add histidines to the *binder*
and score the result. Reasoning about the **target's own protonation** — using
EGFR's histidines as anchor points for acidic residues on our side — requires
modelling a charge state on the molecule you did not design, which off-the-shelf
tooling does not do.

We are claiming this as a **method** contribution rather than a result, and
stating plainly what it is not: it is a design rule derived from first principles,
not a validated predictor. Two known limits are carried openly throughout:

- **The switch is partial, not binary.** At pH 6.5 a histidine is perhaps a third
  to a half protonated, not fully. Each pair contributes a fraction of a charge
  change, which is why three or four pairs are needed rather than one.
- **A histidine's pKa shifts with its neighbours.** The same amino acid flips at a
  different pH depending on what surrounds it. This is a large part of why pH
  selectivity resists reliable computational prediction — including ours.

---

## Target selection: why this patch of the protein

The target is **domain III** of the EGFR extracellular region (roughly residues
310–480), the region both approved antibodies bind and the one the organisers
recommend.

Within it, the candidate epitope is **residues 415–466**, chosen by computed
human/mouse comparison rather than by inspection
([results/findings/alignment-findings.md](results/findings/alignment-findings.md)):

| Region | Human/mouse identity | Differences |
|---|---|---|
| Whole extracellular region (25–645) | ~87% | — |
| Domain III (310–480) | 90.6% | 16 in 171 positions |
| **Candidate epitope (415–466)** | **98.1%** | **1 in 52 positions** (S442G) |

That one difference is serine to glycine — two of the smallest amino acids, so the
local shape barely changes. Fourteen of domain III's sixteen human/mouse
differences fall *before* position 415, which is why this block is the obvious
place to aim if objective 2 (mouse cross-reactivity) is to be satisfied by
construction rather than by luck.

The block contains **eight** anchor points for the pairing rule, all identical in
both species:

- **acidic, pair with binder histidine:** D416, E421, E424, E455, D458, D460
- **target histidines, pair with binder D or E:** H418, H433

Also present and conserved, but useless for the switch because they do not change
charge between the two pH values: R427, R429, K431, R451, K454.

### Why this was not yet enough

Everything above is **sequence** analysis — comparing letters in a row. Letters
tell you what a residue *is*, not which direction it *points*. Domain III folds
into a solenoid, a spiral-staircase shape, in which residues that sit next to each
other in the sequence can point in opposite directions. Some of the eight anchors
are certainly buried inside the protein core and therefore unusable, and a binder
is a single small object that can only touch one patch — its contact points have
to sit within roughly 25 ångströms of one another.

So the epitope was explicitly held as **conditional on a structure check** before
any design work: are the anchors actually on the surface, and do they cluster?

That check is [analysis/02](analysis/) through [analysis/06](analysis/), and its
outcome is recorded in the current-state section below.

---

## Current state

See [docs/decisions-log.md](docs/decisions-log.md) for the authoritative and
current version, split into Settled / Open / Unverified / Ruled out. In summary:

### The structure check passed. The epitope is real.

Computed 30 September 2026 from PDB entries 6ARU and 1NQL.

| Question | Answer |
|---|---|
| Does the PDB numbering match ours? | **No — off by 24.** Both structures number by the mature protein. Derived empirically, 100% consistent. |
| Are the anchors on the surface? | **7 of 8** exposed or partially exposed. H418 buried in 6ARU. |
| Do they cluster within one binder's reach? | **Yes — 4 anchors within 22.6 Å**, against a minimum of 3. |
| What does cetuximab actually touch? | 24 residues; **10 of our 52**, and only **1 of our 8 anchors** (H433). |
| Is the epitope covered when EGFR closes? | **No.** Slightly *more* accessible in the closed structure. |

Four things worth drawing out, because they changed the plan rather than
confirming it:

1. **The numbering trap was real.** 6ARU numbers by the mature protein, so
   UniProt = PDB + 24. Had that been assumed rather than measured, every result
   here would have been computed on the wrong residues and nothing would have
   errored.

2. **A prior assumption was disproved.** The working explanation for why
   cetuximab fails on mouse EGFR named four residue differences. Computation put
   only one of them (K467R) in the contact set, and found two the assumption had
   missed (R377K, S442G). Recorded as disproved rather than quietly corrected.

3. **The method-novelty claim narrowed.** It depends on pairing acidic residues
   against the target's *own* histidines. Of the two available, H418 is buried in
   6ARU — so the claim rests on H433, which is also the single anchor cetuximab
   touches. The mechanistically distinctive choice is therefore the one most
   exposed to a "that's cetuximab's epitope" objection. That tension is real and
   is documented rather than smoothed over.

4. **Two new constraints emerged that no sequence analysis could have found.**
   Position 442 is both the only human/mouse difference in the block and a
   cetuximab contact, so the design must avoid it. And N444 — inside the block —
   carries a covalently attached sugar chain, which means accessibility numbers
   for nearby anchors are upper bounds rather than measurements.

### Still open

The highest-value question is **whether H418 is usable**: 6ARU says buried, 1NQL
says partially exposed, with domain III in the same fold in both (RMSD 1.08 Å).
It matters because if H418 is usable the best cluster grows to five anchors,
carries a target histidine, and sits entirely outside cetuximab's footprint. Two
experimental structures disagree, so this is recorded as ambiguous rather than
resolved in whichever direction would be more convenient.

Also open: which of the four qualifying anchor clusters to design against (a
genuine trade-off, tabulated in the decisions log); which molecule size category
to enter; the design pipeline and compute environment; and the proportion of open
to closed receptor in the assay buffer, which two crystal structures cannot
supply.

---

## Reproducing this

No GPU is needed for anything in `analysis/`. The structure files download from
the public RCSB archive.

```bash
git clone <this repo>
cd egfr-conditional-binder

python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt

# Run in order. 00 and 01 are regression tests: if either fails,
# the environment is wrong and nothing downstream should be trusted.
./.venv/bin/python analysis/00_numbering_check.py
./.venv/bin/python analysis/01_alignment.py
./.venv/bin/python analysis/02_structure_prep.py
./.venv/bin/python analysis/03_solvent_accessibility.py
./.venv/bin/python analysis/04_cetuximab_contacts.py
./.venv/bin/python analysis/05_anchor_geometry.py
./.venv/bin/python analysis/06_tethered_occlusion.py
```

Each script is standalone and rerunnable, writes a human-readable findings file
to `results/findings/` and a machine-readable table to `data/derived/`, and prints
its checks rather than asserting silently.

| Script | What it answers |
|---|---|
| `00_numbering_check.py` | Does the residue-numbering convention hold? (guards every later number) |
| `01_alignment.py` | Where do human and mouse EGFR differ? (regression test) |
| `02_structure_prep.py` | What is the structure's numbering offset, and which chain is the receptor? |
| `03_solvent_accessibility.py` | Which anchors are on the surface rather than buried? |
| `04_cetuximab_contacts.py` | What does the existing drug actually touch? |
| `05_anchor_geometry.py` | Do the surviving anchors cluster tightly enough for one binder? |
| `06_tethered_occlusion.py` | Is the epitope reachable in the closed conformation too? |

---

## Repository layout

```
docs/
  decisions-log.md      source of truth: settled / open / unverified / ruled out
  competition-brief.md  official rules, dates, formats, links
  glossary.md           every term used, defined in plain English
analysis/               numbered standalone scripts (see table above)
data/
  sequences/            UniProt records and official challenge constructs
  structures/           downloaded structures (gitignored, re-downloadable)
  derived/              computed tables, CSV (committed)
results/
  findings/             computed findings, markdown (committed)
  candidates/           bulk design output (gitignored except final shortlist)
design/                 design pipeline configs and runners
submission/             the submission CSV and the methods write-up
explorer/               interactive alignment viewer
```

---

## On honesty in this repository

The competition is judged partly on method novelty, by a process that reads the
submitted documentation. That makes it tempting to overstate. Two deliberate
choices against that:

1. **Anything asserted from memory is labelled UNVERIFIED** and kept in a separate
   section of the decisions log, where it cannot be built on until a computation
   resolves it. When a computation *disproves* a prior assumption, the
   disproof is recorded as a result rather than quietly dropped.
2. **Predicted is not measured.** No design here has been tested in a laboratory —
   the competition requires that (it is "zero-shot"). Where a claim rests on a
   convention, an unvalidated assumption, or a cutoff someone chose rather than
   measured, the text says so.

The author is a product manager learning this field, not a trained protein
designer, and the repository is written to make that verifiable rather than
disguised.
