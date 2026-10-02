# Explainer 05 — preparing what the design tool is given, and filtering what it returns

What a design program called BindCraft2 is, where it runs, what it needs, what it
cannot do for us, and the two scripts (`analysis/09_trim_target.py` and
`analysis/10_charge_pair_filter.py`) that cover those two gaps. A third script,
`analysis/11_trim_boundary.py`, measured one decision inside the first.

Written for a reader with no biology or software background. It stands on its own:
every term is explained where it first appears.

---

## How this document differs from the others

Explainers 01 to 03 describe completed work with results. Explainer 04 describes
machinery not yet set up. This one covers **finished scripts that have not yet
processed a real result**. Scripts 09 and 10 are written and tested against
constructed inputs, and no design run has produced anything for them to work on.
Script 11 has results, from a weaker predictor than the one the design run uses, and
they are reported below with that limit attached.

It is written now because script 10 is the part of this submission that is our own
method, and a method has to be understood before its output is believed. When the
scripts have run on real output, this explainer is rewritten with what came out.

---

## 1. What we are trying to build

**EGFR** (epidermal growth factor receptor) is a protein that sits on the surface of
cells, and it is the target of this competition challenge. We are designing a **binder**:
a small protein of our own, roughly 70 to 100 amino acids long, that sticks to a
chosen patch of EGFR. Amino acids are the building blocks of every protein, and a
protein is a chain of them folded into a particular three-dimensional shape.

The competition asks for a binder with a particular behaviour. It has to stick to
EGFR at **pH 6.5**, the mildly acidic condition found in the surroundings of a
tumour, and show no detectable binding at **pH 7.4**, the condition in healthy
tissue. pH is the scale of acidity: 7 is neutral, and lower numbers are more acidic.
Two further requirements are that the same sequence also binds mouse EGFR, and that
it binds human EGFR with some strength. The competition ranks them in that order, so
the pH behaviour matters most.

**How a binder can be made to care about pH.** Of the twenty amino acids, only one,
**histidine**, changes its electrical charge across that pH range. At pH 7.4 a
histidine is electrically neutral. At pH 6.5 it picks up a positive charge. Opposite
charges attract. So if a binder has a histidine sitting directly across from a
negatively charged amino acid on EGFR, the two ignore each other at pH 7.4 and
attract each other at pH 6.5. That is a switch.

Two arrangements work, and both switch on in the same direction as the pH falls:

- a histidine on the binder facing a negatively charged amino acid on EGFR (aspartic
  acid or glutamic acid, written D and E)
- a D or E on the binder facing a histidine that EGFR already has

Each such pair is a weak effect, and the switch is partial at pH 6.5, so a design
needs three or four of them stacked.

**The anchors.** Earlier work (explainers 01 to 03) chose where on EGFR to build
these pairs. It found one face of the protein carrying eight suitable residues, which
this project calls **anchors**: E344, H358, D368, H370, E391, E400, E421 and E424.
The letter is the amino acid and the number is its position along the chain. Six are
D or E and two are histidines. The binder is to be aimed at those eight.

One earlier rule has been withdrawn and is not applied below. This project first
aimed for deliberately weak binding; the organisers have since described how
binding is scored, and the design target is now the largest gap between pH 6.5 and
pH 7.4 with affinity at pH 6.5 as high as the switch allows. The reasoning is in
`docs/explainers/07-the-assay-and-what-it-changes.md`.

---

## 2. What builds the binder: BindCraft2

**BindCraft2** is a free program for designing binders. It is the program this
project uses to produce candidate binder sequences. (A version of the same tool,
BindCraft version 1, exists and behaves differently in places. This document means
version 2 throughout.)

**What goes in.** Two things, supplied in a settings file called a **campaign
file**:

1. A **structure file** of the target. A structure file is a plain text list of the
   x, y and z coordinates of every atom in a protein, as measured in a laboratory.
   The one used here comes from 6ARU, an entry in the Protein Data Bank, which is the
   public archive of measured protein structures.
2. A list of **hotspots**: the positions on the target the binder should be aimed at.
   For us these are the eight anchors.

**What it does.** It runs a loop, over and over. Each pass through the loop is a
**trajectory**. A trajectory has three steps:

1. Invent a **backbone**, the skeleton of the protein, which fixes its overall shape.
   The aim is a shape whose surface fits the target patch, as a key fits a lock.
2. Choose the amino acid sequence that would fold into that shape. This runs the
   usual direction backwards, from structure to sequence, and is called inverse
   folding.
3. Check the result. A structure predictor, AlphaFold2, is asked whether that
   sequence really folds into the invented shape and sits on the target where
   intended. A design that fails the check is discarded. The check is the program
   marking its own work, so it is good at catching nonsense and weak evidence that a
   binder works in a test tube.

Most trajectories end in rejection, so BindCraft2 runs many to get a few survivors.

**What comes back.** A folder of results. The survivors are called **candidates**.
Each candidate is saved as a **complex**: a structure file holding the target and the
designed binder together in the position the program predicts they would sit. The
folder also holds a table of numbers describing each candidate, including the
program's own confidence scores. BindCraft2 saves complexes in a file format called
**mmCIF** (macromolecular crystallographic information file), the newer standard
format for protein structures.

**Why this project uses it.** The three steps above are normally done by three
separate tools wired together by hand. BindCraft2 bundles them into one automated
loop, and with three days available, wiring and debugging three tools was the larger
risk.

---

## 3. Where it runs, and what Modal is

BindCraft2 needs a **GPU** (graphics processing unit). A GPU is a chip built for
drawing graphics, and the arithmetic that draws a game frame is the same kind of
arithmetic that a neural network such as AlphaFold2 is made of. AlphaFold2 runs
inside BindCraft2's loop many times per trajectory. On an ordinary laptop processor
it would be too slow to be useful, and the cards that suit it are data-centre
hardware, so this project does not own one.

**Modal** is a service that rents one by the second. You write a Python file saying
what program to run and what hardware to run it on. Modal starts a machine somewhere
with that GPU attached, runs the job, hands back the results and shuts the machine
down. Nothing is billed while nothing runs. The two are separate things: BindCraft2
is the work, and Modal is where the work happens.

The file that does this is `design/modal/bindcraft2_smoke.py`. As of this writing it
has not been run, because Modal has to be linked to an account through a login in a
web browser, and only the account holder can do that. Until that happens no
candidates exist, which is why both scripts below were tested on constructed inputs.

---

## 4. The two gaps, and why there are two scripts

BindCraft2 does the central job of inventing binders. Two things around it are
missing, and each is a script.

```
target structure  →  [09]  →  campaign file  →  BindCraft2 on Modal  →  candidates  →  [10]  →  shortlist
```

**Gap 1: the input does not exist yet.** BindCraft2 needs a target structure file and
a hotspot list. The structure file we hold covers all of EGFR's outer region, about
620 amino acids, which makes every trajectory slow and expensive. And the hotspot
list has to be written in the numbering of whichever file is supplied, which turns
out to be easy to get wrong. Script 09 produces the right input.

**Gap 2: the output cannot answer our question.** BindCraft2 produces binders that
fold and bind, and ranks them by how confident it is in the interface. It has no
interest in pH. It reports a binder's charge at a fixed pH of 7.4 as a readout and
does not use it when deciding what to keep. Left alone it would produce strong,
well-folded binders that do not care about pH, which fail the competition's first
requirement. Script 10 is the step that selects, from what BindCraft2 returns, the
candidates that carry the charge switch from section 1.

---

# Script 09: making the input

## What it makes

Some terms first, because the rest of this section depends on them.

A **residue** is one amino acid at one position in a protein chain. "Residue 370"
means the 370th amino acid counting from the start. A **chain** is one continuous
string of residues, and a structure file can hold several chains. A **domain** is a
section of a protein that folds into its own compact blob, rather like one bead on a
string. EGFR's outer region contains four, and the anchors sit in the third, called
**domain III**, which this project defines as residues 310 to 480.

BindCraft2's run time grows with the total size of the complex: the target plus the
binder. Script 09 therefore cuts domain III out of the full structure and writes it
as a smaller structure file. The cut-out piece is called the **fragment** from here
on. It contains 171 residues, which is 28% of the 621 in EGFR's outer region, and it
contains all eight anchors. The target file BindCraft2 receives is this fragment, and
script 09 also writes the campaign file listing the hotspots.

Three things make the cut less simple than taking a slice.

## Complication 1: cutting a protein renumbers it

In the full protein, histidine 370 is the 370th residue. The structure file we hold
counts differently from our project's numbering (all residue numbers in this project
count from the start of the full human sequence in UniProt, the public sequence
archive), and the cut fragment keeps the structure file's own numbers. Hotspot
numbers handed to BindCraft2 have to use the numbering of the file BindCraft2
actually reads. If it is handed "370" in our numbering, it will aim at whatever sits
at position 370 in that file, which is a different residue. The design would then be
built against the wrong surface and nothing would announce it.

Script 09 writes the hotspot list in both numbering systems. It also checks the
mapping by **residue identity** rather than by arithmetic: it confirms that the
position it believes is H370 really is a histidine, that the position it believes is
E344 really is glutamic acid, and so on for all eight, by three independent routes.
Arithmetic can be off by one and still look correct. A histidine that turns out to
be a leucine cannot.

The same class of error produced a disagreement between two earlier analysis steps
in this project.

## Complication 2: cutting can break a disulfide bond

A **disulfide bond** is a chemical link between two cysteine residues (cysteine is
one of the twenty amino acids). The two can be far apart along the chain and touching
in the folded shape, and the link is part of what holds the fold together. Think of a
staple through two pages of a book.

When the fragment is cut out, a staple can have one end inside the fragment and the
other end left behind in the part that was discarded. That cysteine then has no
partner.

EGFR's outer region has 24 disulfide bonds. Of those:

- 3 lie wholly inside the fragment and are unaffected
- 20 lie wholly outside it and are removed together with both ends, which is harmless
- 1 is severed: cysteine 470 is inside the fragment and its partner, cysteine 499, is
  outside it. In the intact protein the two sulfur atoms are 2.04 Å apart. (An
  **ångström**, written Å, is a ten-billionth of a metre, and an atom is roughly one to
  two ångströms across.)

The count of 24 matches the records the original researchers wrote into the
structure file, so it rests on their statement as well as on our own distance
measurement.

**Why it matters.** The fragment is only ever an input to a modelling program. It is
never manufactured, and the experiments are run against real EGFR. The question is
whether an unpaired cysteine makes the fragment fold differently in the model from
the same residues in the intact protein. A design made against a surface that exists
only in the model would be useless.

## Complication 3: the cut edges are the least reliable part

Where a fragment is cut, the chain stops. In the real protein the chain continues,
and the residues that were next to the cut are held in place by neighbours the
fragment no longer has. So a predicted structure is least reliable at the cut ends.

Seven of the eight anchors are 12.3 Å or more from the nearest cut end. One is not:
**E424 is 6.6 Å from the cut at residue 480.** A design that depends heavily on E424
is depending on the part of the model least likely to be right. Script 10 acts on
this when it ranks candidates (see below).

## The decision: where the cut should end

Extending the cut from 480 to 499 would put both cysteines of the severed bond inside
the fragment and move E424 further from the end. It would cost about 20 more
residues of target and brings in the start of the next domain, which might not hold
its shape on its own.

Script 11 measured this instead of arguing it. It predicted the three-dimensional
structure of three versions of the fragment from their sequences alone: 310 to 480 as
cut; the same with cysteine 470 changed to serine, which removes the unpaired
cysteine at the cost of one deliberate difference from the real sequence; and 310 to
499. It then compared each prediction against the same residues in the measured
structure 6ARU. Results are in `results/findings/11-trim-boundary.md`.

The predictor was **ESMFold**, which works from a single sequence and is less
accurate than the predictor the design run uses. It was chosen because it needs no
GPU, and GPU access is still blocked. Every absolute number below carries that
limit. What can be read is the difference between the three versions, since all
three went through the same instrument.

The comparison is expressed as **RMSD** (root-mean-square deviation): the average
distance between matching atoms in two structures after one has been moved on top of
the other as well as it can be. A small number means the same shape.

**The three versions cannot be separated on anything the design depends on.** Over
the residues all three contain, they sit 5.00, 4.97 and 4.97 Å from the measured
structure. Across the eight anchors, lined up on those eight alone, they sit 0.64,
0.64 and 0.65 Å. Those differences are far smaller than this predictor can resolve,
and a fit on only eight points is flattering by construction, so 0.64 Å should not be
read as the anchors being modelled to better than an ångström.

The one place the versions differ is the cut edge. The last ten residues (471 to 480)
sit 4.35 Å from the measured structure when they are the end of the fragment, and
2.95 Å when they sit 19 residues inside it. The predictor's own confidence score
there, called **pLDDT**, rises from 54 to 78 on a 0 to 100 scale. The predicted
470–499 bond forms at 2.22 Å. Against that, E424 itself came out slightly further from
the measured structure in the longer version (2.71 Å against 2.36 Å), which is the
wrong direction for the residue the extension was meant to help, though the difference
is small enough to be noise. Changing cysteine 470 to serine made no measurable
difference.

**The cut stays at 480.** The extension improves residues that only E424 is near,
costs about 11% more target, and shows no benefit at the anchors. The serine change
adds a deliberate difference from the real sequence for no measurable gain, so it is
not adopted. E424's risk is handled in the ranking instead.

What this leaves unsettled: the overall fit of about 5 Å is the same for all three
versions, so the boundary does not cause it, but this measurement cannot say whether
it comes from the predictor or from the fragment. Only the design run's own predictor
can. The question of whether the fragment keeps its shape stays open.

---

# Script 10: the filter

This is the part worth understanding properly.

## What it looks at, and what it decides

For each candidate that BindCraft2 returns, script 10 reads the complex file and
works out which residue of the binder sits directly opposite which residue of EGFR.
Residues that face each other across the join are in contact, and the join itself is
called the **interface**. A single binder residue facing a single target residue is a
**contact pair**.

It then judges each pair against the rule from section 1:

| On EGFR | On our binder | Verdict |
|---|---|---|
| aspartic or glutamic acid | histidine | good pair: switches on as acidity rises |
| histidine | aspartic or glutamic acid | good pair: switches on as acidity rises |
| histidine | histidine | **candidate rejected** |
| anything else | anything | neutral: counted, not scored |

The histidine-against-histidine rejection matters as much as the good pairs. Two
histidines facing each other both turn positive at pH 6.5 and push apart, which can
cancel a good pair elsewhere on the same face. The rule applies to **any** histidine
on EGFR the binder faces, including ones not among the eight anchors, because the
physical problem is identical.

A candidate is also rejected if it contacts residue 442. That is the only position
where human and mouse EGFR differ inside the original candidate patch (residues 415
to 466), and it is also one of the residues the antibody cetuximab touches, so
contacting it risks behaving differently in the two species.

The decision for each candidate is one of: rejected, meets the target (three or more
good pairs; four is the design aim), or below the target. A candidate below the
target is kept and ranked last. It is dropped only for breaking a rule, never for
being weak.

**Ranking.** Candidates are ordered by how many of their good pairs do not rest on
E424, then by the total number of good pairs, then by how many have their charged
groups close enough to reach each other. A design that gets to the same number of
pairs without E424 therefore ranks above one that needs it. A pair on E424 still
counts as good, because it is physically a good pair. It only costs a design its
place against an equivalent design that does not lean on the least reliable part of
the model.

## Why BindCraft2 cannot do this itself

BindCraft2 does report the interface, as two separate lists: the binder residues
involved, and the target residues involved, using a threshold of 4.0 Å. It does not
record which binder residue faces which target residue.

That is a guest list without a seating plan. It says who came to dinner and says
nothing about who sat opposite whom, and every row of the table above needs the
seating plan. The source code of BindCraft2 was read to confirm that no setting
produces it. So script 10 works out the pairing itself, from the three-dimensional
coordinates in each returned complex file.

The consequence governs the whole generation run: **the complex files are the only
usable output.** The summary tables are the guest list. A run that returned tables
and no complex files would have to be done again from scratch.

## Why pair count leads the ranking

Script 10 ranks by good pairs first, because the charge pairs are what produce the
pH switch. BindCraft2's own confidence numbers are carried through and used only as
a tie-break among candidates with equal pair counts, preferring the more confident
interface. They do not override the pair count, since BindCraft2 has no interest in
pH and a high-confidence interface with no charge pairs has no switch.

## What it checks before scoring anything

The target in a returned complex has to be numbered the way the file BindCraft2 was
given was numbered. Otherwise script 10's translation back into our numbering is
reading the wrong residues. So for every candidate it compares the target against the
input file, residue by residue, and sorts the result into four outcomes: identical;
cropped with the numbers kept (accepted, since every position still means what it
did); renumbered; or not matching. The last two are refused, not scored, and the
reason is printed.

---

## Three things that would have failed silently

None of these would have produced an error message.

**The complex files were not being collected.** The script that runs BindCraft2 on
Modal copied results back by file type, and the list did not include mmCIF, the
format BindCraft2 saves complexes in. It would have returned every summary table and
not one complex file: the guest list and nothing else. The copy also flattened the
output into a single folder, which would have destroyed the folder layout script 10
reads and let files with the same name overwrite each other. Both were found and fixed
before any money was spent.

**The chains are the other way round.** A chain is one string of residues, and a
complex file holds at least two: the target and the binder. BindCraft2 puts the target
on chain A and the binder on chain B. Version 1 did the opposite. A script that
trusted the letter would read the wrong molecule and report plausible nonsense.
Script 10 identifies the target by checking which chain reads as human EGFR.

**The target might be cropped and renumbered.** BindCraft2 has an internal measurement
called `Target_Crop_Length`, which suggests it might trim the target. If it did that
and renumbered the residues, script 10 would translate them into our numbering
wrongly and report confident verdicts about the wrong residues.

What is known: the source code defines that measurement as a count of residues that
are not padding, and padding is a device for making batches the same length, so the
reading is that nothing is cropped. The documentation and output code also say the
input's residue numbers are kept. Both are readings of source code and neither has
been checked against what a run returns. What is built: the check described just
above, which refuses to score a candidate whose numbering does not reconcile, and
which was tested by deliberately renumbering a target and by switching the check off
to confirm that test then fails. The first run on Modal, against BindCraft2's own
example target, is set up to make the same comparison on real output. The question is
open until that run happens.

---

## How both scripts are tested

Each script breaks its own checks on purpose, because a check that has never been
seen to fail is not known to work.

Script 09 shifts its numbering by one, confirms that all eight anchors are caught,
and confirms that it stops before writing anything. Script 11 does the same with its
mapping from predicted residues to ours.

Script 10 builds constructed complex files, each designed to exercise one branch of
the rule, and runs them every time. There are eight, including four good pairs, a
histidine facing a histidine, a contact at 442, a design resting on E424, and a target
numbered one place off. It then switches each rule off in turn and confirms that the
matching test fails. One result is worth recording. Switching off the rule about
side chains that are missing from a structure turned a candidate from "below target"
into "meets target". That is the silent inflation the rule exists to prevent, and it
was visible only because the rule was deliberately broken.

Script 10 also builds a fake campaign folder in BindCraft2's documented layout, with
its real column names and with the binder-only files that sit among the real
candidates, and reads it end to end on every run.

---

## What this does not settle

- Neither script 09 nor script 10 has processed real output. Everything above is
  behaviour against constructed inputs.
- Whether BindCraft2 crops and renumbers the target is open. The source says it does
  not, and script 10 refuses a candidate if it does.
- Whether the 171-residue fragment holds the shape the anchors sit on is open. Script
  11 did not separate the three boundaries it tried, and it used a weaker predictor
  than the design run.
- A good pair, as counted here, is a geometric arrangement in a predicted structure.
  Whether it produces the intended switch in a test tube is what the experiment
  decides.
- Script 10 selects from what BindCraft2 returns. It cannot recover a design that was
  never generated, which is why the target and hotspot list from script 09 matter more
  than anything else in this document.

---

## Terms introduced here

For `docs/glossary.md`: anchor, campaign, candidate, chain, complex, contact pair,
disulfide bond, ESMFold, fragment, interface, mmCIF, pLDDT, residue, residue numbering,
RMSD, side chain, structure file.
