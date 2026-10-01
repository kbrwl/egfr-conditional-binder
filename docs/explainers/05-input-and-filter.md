# Explainer 05 — preparing the target, and filtering what comes back

What `analysis/09_trim_target.py` and `analysis/10_charge_pair_filter.py` do, why
each one is harder than it sounds, and the things that would have failed silently.

Written for a reader with no biology or software background.

---

## How this document differs from the others

Explainers 01 to 03 describe completed work with results. Explainer 04 describes
machinery not yet set up. This one sits between them: **both scripts are finished and
tested, and neither has processed a real result yet.**

It is written now rather than after the run for one reason. Script 10 is the part of
this submission that is not simply "we ran BindCraft" — it is the method we are
claiming. If it is only understood after it produces a shortlist, the shortlist gets
approved on trust. A tool has to be understood before its output is believed, not
after.

Rewrite this explainer once both scripts have run on real candidates, with what
actually came out.

---

## Where these two sit

The design tool, BindCraft2, needs a file describing what to aim at, and returns a
pile of candidate binders. Script 09 makes the first. Script 10 processes the second.

```
09  →  BindCraft2 (on a rented GPU)  →  10  →  shortlist
```

Neither script needs a graphics card. Both were written and tested before the design
tool was set up, deliberately: if the filter is written after the candidates exist,
it reliably turns out to need something the generation run did not save, and
everything has to be regenerated.

---

# Script 09 — making the input

## What it makes

BindCraft2's run time scales with the size of the whole complex, meaning the target
plus the binder being designed. The full EGFR extracellular region is about 620
residues. Script 09 cuts out domain III, the 171-residue region from 310 to 480 that
contains the eight anchors, and writes that out as the target file.

Three things make the cut less simple than taking a slice.

## Complication 1 — cutting renumbers everything

In the full protein, histidine 370 is the 370th residue in the chain. Cut away the
first 309 and that same histidine becomes residue 61 of the new file.

BindCraft2 reads the new file. So if it is handed "370" as a target, it will aim at
whatever happens to sit at position 370 of the fragment, which is a different residue
entirely, 309 places further along. The design would be built against the wrong
surface and nothing would announce it.

The script therefore writes the target list in both numbering systems and checks the
mapping **by residue identity rather than by arithmetic**. It confirms that the
position it believes is H370 really does come back as a histidine, that 344 comes
back as glutamic acid, and so on for all eight. Arithmetic can be off by one and
still look correct; a histidine that turns out to be a leucine cannot.

This is not a hypothetical caution. The same class of error produced the
disagreement between steps 04 and 06 earlier in this project.

## Complication 2 — proteins are stapled, and cutting breaks staples

A **disulfide bond** is a chemical staple between two cysteine residues. The two can
be far apart in the sequence and touching in the folded shape, and the staple is part
of what holds the fold together.

Cutting the protein can leave one end of a staple inside the fragment and the other
outside, leaving a cysteine with nothing to bond to.

That happens once here. The chain has 24 such staples. Three sit wholly inside the
fragment and are unaffected. Twenty are removed with both ends, which is harmless.
One is severed: **C470 is inside the fragment and its partner C499 is outside**,
measured at 2.04 Å apart in the intact structure.

The 24 measured bonds match the records the structure's depositors wrote into the
file, so this is checked against independent evidence rather than resting on the
script's own distance threshold.

**Why this matters less than it sounds.** The trimmed fragment is only ever an input
to a modelling program. It is never manufactured, and the real tests are run against
real EGFR. The only question is whether the fragment holds its shape well enough in
the model to design against. An unpaired cysteine makes that slightly less certain,
and the open question is how much.

## Complication 3 — the cut edge is the least trustworthy part

A fragment's ends are where a structure prediction is least reliable, because in the
real protein the chain continues and in the fragment it stops.

Seven of the eight anchors sit 12.3 Å or further from either end. One does not:
**E424 is 6.6 Å from the cut at 480.** Designs that lean heavily on E424 are leaning
on the part of the model least likely to be right, which is why that affects how
candidates are ranked rather than appearing only as a note.

## The open decision

Extending the boundary past 499 would make the severed staple whole and give E424
margin, at a cost of about 20 residues and the risk that the start of the next domain
does not behave well on its own. The way to settle it is to predict each candidate
fragment on its own and compare it against the same residues in the intact structure.
Recorded under Unverified in `docs/decisions-log.md`.

---

# Script 10 — the filter

This is the one worth understanding properly.

## The problem it solves

BindCraft2 will return a few hundred candidates. They will be well-folded, they will
bind, and not one of them will care about pH, because BindCraft2 has no pH term in
what it optimises. It reports binder charge at a fixed pH 7.4 as a readout and
optimises against none of it.

Under this project's objective ordering — switching behaviour first, binding strength
last — a strong pH-blind binder is a failing submission. Script 10 is what converts
the pile into a pH-conditional shortlist. It is the part of the method that is ours.

## The rule it applies

For each candidate, the script reads the three-dimensional structure, works out which
residue of our binder sits opposite which residue of EGFR, and classifies each
facing pair.

| On EGFR | On our binder | Verdict |
|---|---|---|
| aspartic or glutamic acid | histidine | good pair — switches on as acidity rises |
| histidine | aspartic or glutamic acid | good pair — switches on as acidity rises |
| histidine | histidine | **candidate rejected** |
| anything else | anything | neutral; counted, not scored |

Both good arrangements work the same way. At pH 7.4 the histidine is neutral and
nothing happens; at pH 6.5 it becomes positive and is pulled toward the negative
partner. Two different-looking arrangements, one mechanism.

The rejection matters as much as the acceptance. Two histidines facing each other
both become positive at pH 6.5 and push apart, and that repulsion can cancel a good
pair elsewhere on the same face. The rejection applies to **any** histidine on EGFR
the binder faces, not only H358 and H370 — a histidine we never listed causes the
same physical problem.

A candidate is also rejected if it touches position 442, a standing rule from earlier
in the project.

The design target is three or four good pairs. One is not enough, because the
histidine switch is partial rather than all-or-nothing at pH 6.5.

## The discovery that made this harder

BindCraft2 does report the interface. It reports it as two separate lists: the binder
residues involved, and the target residues involved, at a 4.0 Å threshold.

That is **a guest list without a seating plan**. It tells you who came to dinner. It
does not tell you who sat opposite whom. Every rule in the table above needs the
seating plan.

A search of the whole package confirmed that no setting produces it — there is no
verbosity flag and no option that writes contact pairs. So script 10 reconstructs the
seating plan itself, from the three-dimensional coordinates of each returned
structure.

This has one consequence that governs the entire generation run: **the structure
files are not optional output, they are the only usable output.** The summary tables
are the guest list. A generation run that returns tables and no structures has to be
done again from scratch.

## Why it does not rank by binding strength

The competition requires no *detectable* binding at pH 7.4. That is a threshold, not
a ratio. A weak binder that clearly switches passes; a strong binder that switches by
the same proportion can fail, because its pH 7.4 signal is still above the
instrument's floor.

So script 10 ranks by the number of good pairs, carries BindCraft2's strength
measurements through without ranking on them, and leaves them available for choosing
*downward*. This is deliberately contrary to what the design tool optimises for, and
it is the easiest place in this project to ruin the submission by doing the normal
thing.

---

## Three things that would have failed silently

Worth recording, because none would have produced an error message.

**The structures were not being collected.** The script that runs BindCraft2 on the
rented machine copied results back by file type, and the list did not include the
format BindCraft2 writes complexes in. It would have returned every summary table and
not one structure — the guest list, and nothing else. Found and fixed before any
money was spent. The copy-back also flattened the output into a single folder, which
would have destroyed the directory structure script 10 reads and allowed files with
the same name to overwrite each other.

**The chains are the other way round.** BindCraft2 puts the target on chain A and the
binder on chain B, which is the reverse of version 1. Script 10 therefore identifies
the chains by checking which one reads as human EGFR, never by the letter.

**The target may be cropped again.** BindCraft2 has an internal measurement called
`Target_Crop_Length`, which suggests it may trim the target itself. If it does that
and renumbers while doing so, script 10's translation back to our numbering breaks
without complaint, producing confident verdicts about the wrong residues. Unresolved.
The run against BindCraft2's own example target is the cheapest way to find out, and
script 10 refuses to score any candidate whose numbering it cannot reconcile with the
input rather than assuming nothing moved.

---

## How both scripts are tested

Each one breaks its own checks on purpose, because a check that has never fired is
not known to work.

Script 09 shifts its numbering by one, confirms all eight anchors are caught, and
confirms it stops before writing anything.

Script 10 disables each rule in turn and confirms the test case covering it then
fails. One result is worth repeating: switching off the rule about unresolved side
chains flipped a candidate from "below target" to "meets target". That is precisely
the silent inflation the rule exists to prevent, and it only became visible because
the rule was deliberately broken.

Script 10 also builds a fake campaign folder in BindCraft2's documented layout, with
real column names and decoy files included, and reads it end to end on every run.

---

## What this does not settle

- Neither script has processed real output. Everything above is behaviour against
  constructed inputs.
- Whether BindCraft2 crops and renumbers the target is unresolved.
- Whether the boundary should extend past 499 is unresolved.
- A good pair counted here is a geometric arrangement. Whether it produces the
  intended switch in a tube is what the experiment decides.
- Script 10 scores what BindCraft2 returns. It cannot recover a design that was never
  generated, which is why the target list script 09 produces matters more than
  anything in this document.

---

## Terms introduced here

For `docs/glossary.md`: chain, contact pair, disulfide bond, interface, mmCIF,
residue numbering, side chain.
