# Explainer 10 — running our own design campaign on a rented graphics card

How the design run that produces our candidate sequences is set up, what it costs,
what it measures, and what is kept from it. The runner is
`design/modal/egfr_campaign.py`. It replaces nothing: the earlier smoke-test runner
`design/modal/bindcraft2_smoke.py` stays where it is.

Written for a reader with no biology or software background. It stands on its own:
every term is explained where it first appears.

---

## How this document differs from the others

**This covers a runner that has produced no candidate yet.** It was written while the
first validation run was still going. Every number below is either a measurement from
a previous run against a different protein, or a price, or arithmetic — none of it is
a result from our own campaign. The document is rewritten when there is real output,
and the sections that are placeholders say so where they appear.

It is written now for the reason the explainer rules give: a tool is understood before
its results are believed. When the candidate sequences arrive they will come with
numbers attached, and those numbers are only worth anything if the thing that produced
them is understood first.

---

## 1. What this is, in one paragraph

We are designing a small protein, called a **binder**, that sticks to a chosen patch
of **EGFR** (epidermal growth factor receptor), a protein on the surface of human
cells. The designing is done by a program called **BindCraft2**, which needs a
**GPU** (graphics processing unit — a specialised processor originally built for
rendering graphics, now the standard hardware for this kind of modelling). We do not
own one, so we rent one by the second from a company called **Modal**. Until now the
only thing we had run on that rented hardware was BindCraft2's own worked example,
against an unrelated protein. `egfr_campaign.py` is the first code that points the
whole machine at *our* target with *our* settings.

---

## 2. Why it is a separate file from the smoke test

On 2 October a run called the **smoke test** was carried out. The name is from
electronics: switch on the new circuit and see whether smoke comes out. Its purpose
was to prove the pipeline works end to end — the software installs, the graphics card
is visible, the model weights download, a design run completes, and a sequence comes
out the other side. It deliberately used BindCraft2's own shipped example target,
**PD-L1** (programmed death-ligand 1, an unrelated protein), so that if anything
broke, our configuration could not be the cause.

It worked. So the question became how to run the same pipeline against our own target
without disturbing the arrangement that is known to work.

The alternative was to edit the smoke-test file so it could do both. That was rejected
for a specific reason: the smoke test is the thing we go back to when something new
fails, to find out whether the problem is ours or the environment's. A file that has
been edited to do two jobs is no longer a reliable answer to that question. So
`egfr_campaign.py` is a new file, and the parts of the smoke test that are already
known to work on real output — the memory sampling, the trajectory counting, the
numbering check — were carried over unchanged rather than rewritten.

---

## 3. The two runs, and why the order is not optional

### Step one: validate

One attempt. The point is to find out whether our campaign file starts up at all
against the real target fragment. This is listed in `docs/decisions-log.md` under
Next actions as item 2, and the entry says plainly that it has never been done and is
a prerequisite before spending real money.

There are two specific things that have never been exercised:

- **`termini_accessible`**, a setting that tells BindCraft2 to point both ends of the
  binder away from the target. We need it because our designs are glued to the
  measuring instrument by their tail end, so that tail must not be part of the
  sticking surface.
- **The His-tag off-target.** A **His tag** is a short run of histidine residues added
  to a protein so it can be caught and purified. The target used in the competition's
  measurement carries one and the organisers expect to leave it on. Our whole design
  approach builds pockets that grip histidines, so there is a real risk of designing
  something that grips the tag rather than the protein. The fix is to hand BindCraft2
  the tag as a thing to *avoid* while it designs. Whether BindCraft2 accepts a plain
  sequence in that role, as opposed to a full three-dimensional structure, had never
  been tested.

A validation run costs cents. Discovering either problem partway through a run costing
several dollars is the thing it exists to prevent.

### Step two: the pilot

150 to 200 attempts. This is the run that produces our first real candidate sequences
and, more importantly, measures the two rates that tell us how big the main run should
be. Explained in section 5.

---

## 4. What a "trajectory" is, and why most of them are thrown away

A **trajectory** is one attempt at designing a binder. BindCraft2 starts from nothing,
grows a protein shape against the target patch, works out a sequence of amino acids
that would fold into that shape, then predicts what that sequence would actually do
and checks the prediction against its own quality thresholds.

Most attempts fail those thresholds and are discarded. On the PD-L1 smoke run, 10 of
30 attempts were accepted — about a third. That is the normal shape of this work: you
pay for every attempt and keep a minority.

The word **accepted** below always means "cleared BindCraft2's own quality
thresholds". It does not mean the design is any good for *our* purpose, because
BindCraft2 has no interest in pH. It reports the binder's electrical charge at a
fixed pH 7.4 as a readout it then ignores, and there is no pH term anywhere in what
it is trying to optimise. Selecting for the pH switch is entirely our own job,
afterwards, and that is what `analysis/10_charge_pair_filter.py` does.

---

## 4a. The gates an attempt has to pass, and the two that stop ours

An attempt is not judged once at the end. It passes through named stages — screen,
refine, anneal, harden, mutate, then a final check — and at most stages it is
measured against a minimum and dropped if it falls short. Knowing which minimum
stopped an attempt is the whole of knowing why it failed, so these are worth
setting out. Read from BindCraft2's `bindcraft/filters.py` and
`settings/core/reference.json` on 2 October 2026.

Three numbers do nearly all the work:

| what it measures, in plain terms | its name | the minimum, by stage |
|---|---|---|
| how confident the model is in the whole predicted complex | `pLDDT.<target>` | screen 0.6, refine 0.6, mutate 0.6, anneal 0.65, harden 0.65, final 0.7 |
| how confident it is that the binder and target really sit together as drawn | `i_pTM.<target>` | screen and refine have no minimum, anneal 0.5, harden 0.5, mutate 0.5, final 0.7 |
| how confident it is in the target's own shape | `Target_pLDDT` | final stage only, default 0.6 |

**pLDDT** is the predicted local distance difference test, the model's confidence
in where it has put each part, from 0 to 1. **i_pTM** is the interface predicted
TM-score, its confidence specifically about the join between the two molecules
rather than about either piece alone. Both are the model's opinion of itself, so a
high score is not a measurement of binding — a prediction can be confident and
wrong.

Two things in that table caused a wrong diagnosis here and are worth pointing at.

**The first two rows are different metrics with similar names.** `pLDDT.<target>` is
checked at every stage; `Target_pLDDT` is checked only at the very end. A setting was
changed on the assumption that it governed the gate doing the rejecting, and it did
not — it governed the final-stage one, which no attempt had reached. The lesson is
that the rejection message names the metric that stopped the attempt, and that name
is what to act on, rather than a setting whose name resembles it.

**The middle-stage interface minimum is 0.5, not 0.7.** The 0.7 applies at the final
check. An attempt reading 0.56 at anneal has passed, not nearly failed, and reading
the middle stages against 0.7 makes designs look further off than they are.

## 5. The two numbers the pilot exists to measure

The campaign file currently asks for up to 2,000 attempts and 200 accepted designs.
Those numbers were never derived from a measured run time — the notes beside the file
(`design/configs/egfr-domain3-h370.md`) say so. They are placeholders. The pilot
replaces them with something derived from measurement.

**Candidates per card-hour.** How many accepted designs we get per hour of rented
graphics card. Everything about cost follows from this. We cannot borrow the figure
from the smoke run, for two reasons that both push it in unknown directions: that run
used a different and more expensive card (an A100 40GB, where we now use an L4), and
it was working on a 115-residue target where ours is 171 residues. Run time grows with
the size of the thing being modelled.

**Charge-pair survival rate.** Of the designs BindCraft2 accepts, the share that pass
*our* test. Our test counts **charge pairs**: places where a residue on the binder sits
across from a residue on the target such that the pair attracts at pH 6.5 and does not
at pH 7.4. We require at least three. This rate has never been measured on real output
because no real output has ever existed — it has only ever been run on complexes built
by hand to test the code.

Multiply the two together and you get designs-that-pass-our-filter per hour, which is
the number that converts a budget into a campaign size. That conversion is the whole
point of the pilot.

---

## 6. Getting our own files onto a rented computer

This is mundane and it is also the thing most likely to go wrong silently, so it is
worth a section.

When you rent a machine from Modal, your code runs inside a **container** — a sealed,
prepackaged computer environment, built from a recipe called an **image**. Our image
contains Ubuntu Linux, the graphics-card libraries, and BindCraft2 installed from
source. Building it takes a while and about 20 GB, so it is built once and cached.

The smoke test needed nothing of ours inside that container, because it ran
BindCraft2's own example. Our campaign needs four files that are not in the image:

| file | what it is |
|---|---|
| `data/structures/6aru_domain3.pdb` | the trimmed target: 171 residues cut out of a measured EGFR structure |
| `data/sequences/his-tag-offtarget.fasta` | the His tag, written as a plain sequence |
| `design/configs/egfr-domain3-h370.json` | the campaign settings |
| `design/configs/egfr-domain3-h370-notag.json` | the same, plus the tag as a thing to avoid |

The mechanism for putting them there is `add_local_file`. **This was checked against
the installed software rather than recalled,** which matters: an older way of doing
the same thing, `copy_local_file`, has been removed from the version we have (Modal
1.6.0). Recalling the old name would have produced an error at the worst moment. This
is the repository's "compute rather than recall" rule applied to a software interface
instead of a measurement.

There is a choice in how the files are added. The default attaches them when a
container starts; the alternative bakes them permanently into the image. The default
was chosen, because baking them in means that editing a settings file forces a rebuild
of the 20 GB image behind it.

### The part that is easy to get wrong

**The files have to sit in the same arrangement inside the container as they do in our
project folder.** This is not tidiness. Our campaign file refers to the target as
`../../data/structures/6aru_domain3.pdb` — "go up two folders, then down into data".
BindCraft2 works that out relative to wherever the campaign file itself is sitting
(confirmed by reading its own reference documentation). So if the four files were
dumped together into one folder, that instruction would point at nothing.

The alternative was to rewrite the paths into absolute ones inside the container. That
was rejected because the campaign file is generated by `analysis/09_trim_target.py` and
documented as not to be hand-edited, and because it would mean the file BindCraft2
reads is not the file in our repository. The self-test asserts both that the mirrored
arrangement resolves correctly and that flattening it breaks — so the thing we rely on
has been seen to fail when removed.

---

## 7. What was checked before spending anything

Four questions were answered by reading BindCraft2's own documentation inside the
container, using a helper that runs without a graphics card and therefore costs almost
nothing. Each would otherwise have been discovered by a failed run that we paid for.

1. **Can a plain sequence be used as a thing to avoid?** Yes. BindCraft2's reference
   says the target path takes "PDB, mmCIF or FASTA" — the first two are structure file
   formats, the third is a plain sequence. So handing it the His tag as a sequence is a
   supported arrangement. `docs/decisions-log.md` lists this under Unverified; this
   narrows the question from "is this even allowed" to "what does it do to the run",
   which only running it answers.
2. **Will it chop up our ten-residue tag?** No. The setting `crop_fasta_sequence` takes
   a window, and `false` means use the whole sequence. Our file sets `false`.
3. **Is `termini_accessible` a real setting?** This one was briefly alarming. It does
   not appear in the list of 134 settings the program prints, and the decisions log
   warns that an invented key could make BindCraft2 refuse the campaign file outright.
   It is valid. The printed list covers only settings that can be overridden from the
   command line, which is a narrower set than the keys a campaign file may contain, and
   BindCraft2's own README documents `"termini_accessible": true` as the campaign-file
   spelling. Recorded here because the next person to check will have the same scare.
4. **Can the number of attempts be capped from the command line?** Yes,
   `--set max_trajectories=N`. This matters because it lets the pilot be bounded
   without editing the committed campaign file, whose own counts describe the eventual
   main run and are still unconfirmed.

---

## 8. How the money is bounded

The budget is a $30 monthly allowance from Modal, of which roughly $3.23 went on the
PD-L1 smoke run. The card we use, an L4, costs $0.80 an hour, so a dollar buys
4,500 seconds.

Two independent limits, because each has a failure the other does not cover.

**The attempt cap** stops the run once it has spent the number of attempts we asked
for. This is the limit that should fire, and when it does the run ends by itself.

**A time limit on the container** is the backstop for the case where an attempt hangs,
or the rate turns out far worse than measured. It is set from a spend ceiling stated in
dollars rather than in seconds, so the limit is expressed in the same unit as the
budget.

The time limit kills the container, and a killed container loses everything on its
local disk — which for a multi-hour run would mean losing output we had already paid
for. So a background task copies the run's output folder to permanent storage every
five minutes. A run cut short still yields both rates, because both are rates and a
rate can be computed from however many attempts finished.

### Being able to see what is happening

The smoke test held all its output until the run finished, which meant that while it
was going, a healthy run and a hung one looked identical. The decisions log records
the consequence: that run was left unattended and nobody went back for its final
report. For a one-attempt validation this does not matter. For a pilot lasting hours
it is the difference between noticing a stall in ten minutes and noticing it in six,
and the gap is paid for in rented card time. The runner now prints each line as it
arrives. Four self-test cases cover it, including a check that a failed command still
reports its failure — without which a dead run would read as a finished one.

---

## 9. A bug from the smoke run that this runner does not repeat

Worth describing because it shows how a check can produce a confident wrong answer.

BindCraft2 preserves the target's residue numbering in its output. We depend on that:
`analysis/10` translates the numbers in a returned structure back into our own
numbering, and if the numbering had silently changed, it would produce confident
conclusions about the wrong residues. So the smoke run compared the target in every
returned structure against the file that went in.

The copying step gathered every structure file under the whole BindCraft2 installation
folder. That swept up four template structures BindCraft2 ships for other kinds of
binder — an antibody fragment, a nanobody, and two others — which have nothing to do
with the run. The check dutifully reported that a nanobody template does not resemble
PD-L1, which is true and meaningless, and the run's printed conclusion said the target
had **not** come back correctly and warned against starting an EGFR campaign.

That conclusion was wrong. All 91 genuinely designed structures read as identical.

The new runner guards against it in two separate places: it copies only the run's own
output folder, and the numbering check reads only the three subfolders BindCraft2
writes its results into. Two guards rather than one, because they fail differently — a
wrong copy loses files we paid for, and a wrong check list produces a confident false
verdict. The self-test includes a folder containing both real output and those template
files, and asserts the templates do not reach the check.

---

## 10. What is kept

Everything the run produces is kept, and that is deliberate rather than incidental.

The runner copies the entire project folder, which includes every attempt, not only the
accepted ones. Failed attempts are evidence: they are what the acceptance rate is
computed from, and a record that only contained successes would make the pipeline look
better than it is and would make the rate uncomputable.

A distinction worth being clear about, because the two are easy to confuse:

- **Kept** means the files exist, on permanent storage attached to the rented machine
  and copied back to this one. Everything is kept.
- **Committed to version control** means the files are tracked in the project's
  history. Bulk design output is deliberately not, because it is large; `.gitignore`
  admits only the final shortlist. That is a decision about repository size and says
  nothing about whether the data is retained.

### The gap that was open here, and what closes it

`analysis/10` reads only the folder of designs BindCraft2 accepted. So the attempts
BindCraft2 discarded were retained as files and then never tested, and there was no
per-attempt record of what was tried or why it did not survive. For a run where most
attempts are discarded, that is most of what we pay for going unrecorded.

`analysis/18_campaign_inventory.py` closes it. It writes one row per attempt, covering
every attempt, with the outcome at each of two gates and the reason:

- **Gate one, BindCraft2's own thresholds.** Did the attempt fold into a confident
  interface? Passing says nothing about pH, because BindCraft2 has no pH term in what
  it optimises.
- **Gate two, our charge-pair floor.** At least three correct pairs, no histidine on
  the binder facing a histidine on the target, no contact at position 442.

It adds no rule of its own. The charge-pair rules are imported from `analysis/10`
rather than restated, by the same mechanism `analysis/15` already uses, so the two
cannot drift apart — this repository has already been bitten once by the same quantity
being computed twice and disagreeing.

Two things it produces are worth knowing about:

- `data/derived/18-file-inventory-<campaign>.csv` lists every retained file with its size and a
  **SHA-256 fingerprint** — a short code derived from the file's contents, such that
  two files with the same code hold the same bytes. This is what lets a structure be
  shown later to be the one this inventory described rather than assumed to be.
  Because `data/derived/` is tracked in version control, the *record* survives even
  where the large structure files themselves are not committed.
- A coverage table, deliberately including the count of designs BindCraft2 accepted
  that carry no verdict from us. That combination should never occur, and reporting
  it as zero every time is the only way a non-zero would ever be noticed.

**One thing this makes newly answerable.** If attempts BindCraft2 discards turn out to
pass our charge-pair gate, then its thresholds are throwing away the thing we are
paying for, and the main run should be configured differently. Nothing could see that
before, because nothing looked at the discarded attempts at all.

---

## 11. What this does not establish

- **Nothing about whether the designs are any good.** The pilot measures rates and
  produces first candidates. Whether those candidates switch with pH is not settled by
  any of this, and cannot be settled by a prediction.
- **Nothing about the trimmed fragment holding its shape.** We hand the design run a
  171-residue piece cut out of a much larger protein. Whether that piece keeps the
  shape it has in the intact protein is still unresolved
  (`docs/decisions-log.md`, Unverified). A binder designed against a shape that does
  not exist would score well and fail in the laboratory.
- **The validation run's timing is an upper bound, not an estimate.** The first attempt
  in any run carries the one-off cost of preparing the models for the specific card,
  which the smoke run saw directly: its two attempts that included that preparation
  took 207 and 634 seconds against 185 to 574 for the rest.
- **Whether the cheaper card is fast enough has not been measured.** The switch from the
  A100 to the L4 was made on memory alone: the smoke run peaked at 18,218 MiB, which
  fits inside the L4's 24,564 MiB. No speed comparison between the two has ever been
  made. There is also an open question about how many parallel workers fit: BindCraft2's
  own stated memory formula works out at roughly 12.4 GB per worker for a complex of our
  size, which on a 24 GB card with headroom allows one worker where the A100 ran two.
  That is arithmetic from their documentation rather than a measurement, and if it holds
  it roughly halves throughput on top of any difference in clock speed. The pilot is
  what settles this.

---

## 12. Where to look

| | |
|---|---|
| the runner | `design/modal/egfr_campaign.py` |
| its self-test | `./.venv/bin/python design/modal/egfr_campaign.py` — no account, no card, no network |
| the earlier smoke test | `design/modal/bindcraft2_smoke.py` |
| the campaign settings, and the reasoning behind them | `design/configs/egfr-domain3-h370.md` |
| the charge-pair test applied afterwards | `analysis/10_charge_pair_filter.py`, explained in `docs/explainers/05-input-and-filter.md` |
| why the His tag is avoided | `docs/explainers/09-the-his-tag-screen.md` |
| what the trimmed target is and what cutting it cost | `docs/explainers/06-checking-the-trimmed-target.md` |
| where the project stands | `docs/decisions-log.md` |
