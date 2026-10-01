# Explainer 04 — the design pipeline, and the machinery it needs

What actually produces binder sequences, why it needs a rented graphics card, what
Modal is, and where the H370 patch from explainer 03 gets consumed.

Written for a reader with no biology or machine-learning-infrastructure background.

---

## How this document differs from the other explainers

Explainers 01 to 03 describe completed work and quote numbers from
`results/findings/`. This one describes work that has **not been done yet**: the
design pipeline is the project's open blocker as of 1 October 2026.

It exists because the next session needs to understand what it is setting up before
setting it up, and because the external figures below were checked on 1 October 2026
and will otherwise be recalled rather than looked up.

**Every number in this document is external and unverified against our own runs.**
Sources are named inline with the date they were read. Treat them as planning
estimates. When the pipeline has actually run, rewrite this explainer with what
happened, and move the figures that matter into a findings file produced by a script.

---

## The gap this fills

The structure work answered where to aim: eight residues on one face of EGFR. It
produced no binder.

A **binder**, in this project, is a small protein of our own design — roughly 70 to
100 amino acids — that sticks to a chosen patch of a larger protein. It cannot be
written by hand. With 20 amino acids to choose from at each of about 80 positions,
the number of possible sequences exceeds the number of atoms in the observable
universe, and nearly all of them fold into something useless. The practical method
is to let machine-learning models propose shapes and sequences, and then check their
own proposals.

---

## What the pipeline does, in three jobs

### 1. Invent a backbone

The **backbone** is the protein's skeleton: the chain of atoms that defines its
shape, before deciding which side chains hang off it. The job is to invent a shape
whose surface is complementary to the target patch, the way a key's shape is
complementary to a lock.

The widely used tool for this step alone is **RFdiffusion**.

### 2. Choose the sequence

Given a shape, which 80 amino acids will fold into it? This is **inverse folding**:
the normal direction is sequence to structure, and this runs the other way, structure
to sequence. The standard tool is **ProteinMPNN**.

### 3. Check it

Feed the proposed sequence into a structure predictor, **AlphaFold2**, and ask two
questions. Does this sequence fold into the shape we asked for? Does it sit on the
target where we wanted?

If the prediction disagrees with the design, the design is discarded. This is a
**self-consistency check**, which is the model marking its own homework. It is weak
evidence about whether the binder works in a tube, and it is good at catching
obvious nonsense, which is what it is for.

### BindCraft2 bundles all three

**BindCraft2** runs those three steps in one automated loop, repeatedly. That is why
this project chose it over wiring three tools together by hand with three days
available. The version matters: statements about BindCraft **v1** do not all carry
over, and this document means v2 wherever it says BindCraft2.

One attempt through the loop is a **trajectory**. Most trajectories produce nothing:
they get rejected at one of the checks. In a published sample run against the test
target PD-L1, a trajectory that finished in under five minutes still ended with no
accepted designs (BindCraft v1 run log published by the US National Institutes of
Health high-performance computing group, read 1 October 2026). Many trajectories are
run to get a few survivors.

---

## Why this needs a rented graphics card

AlphaFold2 runs inside that loop, many times per trajectory. It needs a graphics
processing unit, a **GPU** — not for drawing anything, but because the chips built to
render game frames happen to be good at the matrix arithmetic that neural networks
are made of.

BindCraft2 requires an NVIDIA card with CUDA support (CUDA is NVIDIA's programming
interface for using the card for general computation). How much memory that card
needs is not stated anywhere: the check recorded in `docs/decisions-log.md` found no
minimum video memory in BindCraft2's documentation, and the entry is marked
unverified. The card choice is therefore a judgement rather than a requirement.

`design/modal/bindcraft2_smoke.py` defaults to an A100 with 40 GB on the reasoning
that the folding step is the memory-hungry part, and that a run failing for want of
memory costs more time than a larger card costs money. It drops to a cheaper card
once something is known to work. A laptop has nothing comparable to any of these, so
the card is rented by the hour either way.

---

## What Modal is

**Modal** is serverless GPU rental. "Serverless" means no machine is rented that sits
idle costing money. You write a Python file describing what should run and on what
hardware, run `modal run yourfile.py`, and Modal starts a container somewhere with
that GPU attached, runs the job, returns the results, and shuts it down. Billing is
per second of actual execution.

The file is short, and roughly says: build an image with these Python packages,
attach a GPU of a named type, run this function, write output to this storage
volume. `design/modal/bindcraft2_smoke.py` is that file for this project.

Why Modal over the alternatives:

- **Google Colab** gives free GPUs, and it disconnects, hands out whatever card is
  spare, and makes a reproducible run difficult. This project's central rule is that
  everything regenerates by rerunning a script, and Colab works against that.
- **A plain cloud virtual machine** on AWS or Google Cloud means installing drivers,
  CUDA and Python environments by hand, and paying while it idles, including the
  hours you forget to shut it down.
- **Modal** costs more per GPU-hour than the cheapest providers and removes the
  setup, which is the thing that consumes a three-day budget.

---

## The setup

1. `pip install modal`, then `modal setup`. This opens a browser, you log in, and it
   writes an authentication token to the machine. That is the browser step: Modal
   needs to know which account to bill.
2. **Attach a payment card, then immediately set a spending limit.** Published rates
   as of 1 October 2026: about $5 per month in credits with no card attached, $30 per
   month with one. Setting the limit at $30 makes overspending impossible rather than
   unlikely. Do this before the first run.
3. Build the image from BindCraft2's own installer, once, and let it cache.
   `design/modal/bindcraft2_smoke.py` in this repository already does this: it
   clones BindCraft2, runs its `install.sh` at image build time so the work is not
   repeated on every run, and fetches model weights into a separate volume so a
   rebuild does not re-download them. BindCraft2 ships container recipes but no
   prebuilt image in a registry, so the image has to be built. Note that
   BindCraft **v1** depends on PyRosetta, which is slow to install and carries a
   licence; BindCraft2 does not mention it anywhere, so that obstacle is gone.
4. **Smoke run on the shipped example first.** BindCraft2 includes a PD-L1 example
   target. Run it unmodified and get one design out. A **smoke run** is a test that
   proves the chain works end to end — authentication, image build, GPU allocation,
   file output — before any of our own choices enter it. Pointing it at EGFR first
   means a failure cannot be attributed: epitope, settings or CUDA initialisation
   would all look the same.

---

## The decision that controls cost and time

Run time depends on the size of the whole complex: target residues plus binder
residues. There is no published figure to plan against. BindCraft2's documentation
does not state a runtime per design, which `docs/decisions-log.md` records as
unverified, so what follows is reasoning about direction and scale rather than
arithmetic.

The full EGFR extracellular region is about 620 residues. With an 80-residue binder
that is roughly 700. Domain III is 171 residues, so trimming to it and adding the
same binder gives about 250 — under four tenths of the size. Because the cost is
driven by the size of the whole complex rather than by the binder alone, that is a
large reduction in run time.

**How large is unknown.** It is the difference between a budget that covers this
exercise and one that does not, and no number here establishes which. The smoke run
is what produces it: time one trajectory on the untrimmed target and one on the
trimmed one, and the ratio stops being a guess.

**So trim the target.** Give BindCraft2 only the part of EGFR around our patch.
Domain III contains or neighbours all eight anchors, which makes it the natural
unit to keep.

The trim is a real decision with a real failure mode. Cut too tight and the fragment
will not hold its shape in the model, and the design will have been made against a
surface that does not exist on the intact protein. Keep the whole domain rather than
carving out the eight anchors.

---

## Where the epitope work is consumed

BindCraft2 takes **hotspots**: the residues on the target it is told to aim at. That
input is the output of `analysis/08`, written out:

> **E344, H358, D368, H370, E391, E400, E421, E424**

Everything in explainer 03 exists to produce that one line of configuration. The
funnel from 154 residues to 16 to 8, the face test, the conservation filter — all of
it was deciding what goes in that field. Hotspot choice is also where binder
campaigns commonly fail, because a badly chosen set gets every trajectory rejected
without saying why.

---

## What the pipeline will not do

BindCraft2 does not design for pH. The check recorded in `docs/decisions-log.md`
found that it reports binder charge at a hard-coded pH 7.4 as a suppressed readout —
`bindcraft/filters.py`, `REPORTED_PH = 7.4`, the metrics `Binder_pI` and
`Binder_Net_Charge` — with no pH term in its objective and no pH setting in its
235-setting catalogue. It does contain pH-dependent code, including a side-chain
pKa table, so the accurate statement is that pH is a reported readout rather than
something the tool optimises against. An earlier version of this project's own
notes said the tooling had no pH handling at all; that was refuted by one search of
the source and is recorded under Ruled out.

It will produce a strong, well-folded, completely pH-blind binder. Under this
project's objective ordering, that is a failing submission. So the pipeline is the
first half of the work and the second half is ours:

- generate several hundred candidates against the eight anchors
- for each candidate, look at which target residue every binder contact position
  faces, and apply the positional pairing rule: histidine on the binder opposite the
  acidic anchors, aspartic or glutamic acid opposite H358 and H370
- reject any candidate with a histidine on the binder facing a histidine on the
  target
- then deliberately weaken the survivors, because the requirement is no *detectable*
  binding at pH 7.4 and BindCraft2's built-in objective pushes the opposite way

That filtering and rescoring step is where the submission is decided. It needs
candidates to operate on, which is why the pipeline is the blocker rather than a
convenience.

---

## Planning figures

All external, all read 1 October 2026, none verified against our own runs.

| Figure | Value | Source |
|---|---|---|
| Modal free credits, card attached | $30 per month | Modal pricing pages, as reported by several pricing trackers |
| Modal free credits, no card | about $5 per month | same |
| Concurrent GPU jobs on the free tier | 10 | same |
| Cost per accepted design | about $2.90 | Adaptyv's own published tool comparison, averaged over 7 targets |

**What this table got wrong, and in which direction.** It originally carried three
more rows: a minimum of 32 GB of video memory with 48 GB recommended, a 250-residue
trajectory at about five minutes on an H100, and a 900-residue trajectory at two to
three hours. All three were removed on 1 October 2026 after checking them against
`docs/decisions-log.md`. Each was quoted from BindCraft **v1**'s documentation while
this project uses **BindCraft2**, whose documentation states neither a minimum video
memory nor a runtime per design — both already recorded there as unverified. The
error ran one way: it made the hardware requirement look like a fixed threshold and
the trimming decision look like arithmetic, when the first is a judgement and the
second is a direction of unknown magnitude. The rows left above survived the check.
The removed figures come back only as measurements from our own smoke run.

Two consequences worth stating. Twenty designs sits comfortably inside $30 with a
trimmed target, and the untrimmed case is unknown rather than merely worse. And
BindCraft2 cannot split a single job across several cards, though several jobs can
run at once writing into the same output folder, which is how wall-clock time is
compressed on the free tier.

---

## What this does not settle

- None of the figures above has been reproduced here. The first real run replaces
  them.
- The self-consistency checks inside BindCraft2 say whether a design is internally
  coherent. They say nothing about whether it binds in the assay.
- Nothing in this pipeline addresses the pH requirement, which is the objective the
  competition ranks first.
- How much of the trimmed domain III is needed for the fragment to hold its shape is
  a judgement, not a measurement, until something tests it.

---

## Terms introduced here

For `docs/glossary.md`: AlphaFold2, backbone, BindCraft2, binder, CUDA, GPU, hotspot,
inverse folding, Modal, ProteinMPNN, RFdiffusion, self-consistency check, serverless,
smoke run, trajectory.
