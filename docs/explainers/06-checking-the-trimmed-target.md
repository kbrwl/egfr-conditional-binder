# Explainer 06 — checking the cut-down target against the whole receptor

Whether cutting one piece out of EGFR and designing against that piece alone can
hide the surface we are aiming at, and the two scripts that answer it:
`analysis/12_fragment_context_check.py` and `analysis/13_full_receptor_clash.py`.

Written for a reader with no biology or software background. It stands on its own:
every term is explained where it first appears.

---

## How this document differs from the others

The two scripts here are at different stages, and the difference matters when
reading their numbers.

Script 12 has real results. It measures published, laboratory-determined structures
of EGFR that we already hold on disk, so its numbers are measurements of the real
molecule and are reported as such.

Script 13 has no real results. It screens the output of the design run, and no
design run has produced anything yet, so every number it has ever printed comes from
inputs built on purpose to test it. This explainer says so wherever its numbers
appear, and is rewritten when there is real output.

Numbers below are quoted from `results/findings/12-fragment-context.md` and
`results/findings/13-full-receptor-clash.md`. Where this document and those files
disagree, those files are correct.

---

## 1. The result

The face of EGFR we are aiming at is not covered up by the rest of EGFR. Six of our
eight target positions have nothing from the rest of the protein within 8 ångströms
of them, and 13 ångströms or more of clear space to the nearest piece. None of the
eight is in contact with the rest of the protein at all. The same holds in both of
the two shapes EGFR is known to adopt.

An **ångström**, written Å, is a ten-billionth of a metre. An atom is roughly one to
two ångströms across, and the small protein we are designing is about 20 to 30
ångströms across, so 13 ångströms of clear space is room for part of it to sit.

What this changes about the design: the main risk created by cutting a piece out of
EGFR turns out not to apply to the face we chose. Designing against the cut-out
piece is now supported by a measurement rather than by an assumption, and the
campaign can proceed on that basis. A separate risk from the cut, whether the piece
holds its shape once separated, is untouched by this and stays open. Section 7 says
what would settle it.

A second thing came out of this. Of the eight positions, E421 and E424 are the only
two with any neighbouring piece of protein close by. Earlier work found that E424
also sits 6.6 Å from one of the cut edges, which is the least reliable part of a
predicted structure. E424 is therefore the weakest of the eight positions on two
separate counts, which supports a rule `analysis/10_charge_pair_filter.py` already
applies: a design that leans on E424 ranks below an equivalent design that does not.

---

## 2. What we are aiming at, and the risk in how we aim

**EGFR** (epidermal growth factor receptor) is a protein on the surface of cells and
the target of this competition challenge. A **protein** is a chain of building
blocks called **amino acids**, folded into one particular three-dimensional shape. A
**residue** is one amino acid at one position in the chain, so "residue 370" means
the 370th one counting from the start. Every residue number in this project counts
from the start of the full human sequence in **UniProt**, the public archive of
protein sequences.

We are designing a **binder**: a small protein of our own, roughly 70 to 100
residues long, that sticks to a chosen patch of EGFR. The competition requires it to
stick at pH 6.5, the mildly acidic condition around a tumour, and to show no
detectable sticking at pH 7.4, the condition in healthy tissue. pH is the scale of
acidity, where 7 is neutral and lower numbers are more acidic.

The switch is built from pairs of oppositely charged residues. Of the twenty amino
acids, only **histidine** changes its electrical charge across that pH range: it is
electrically neutral at pH 7.4 and picks up a positive charge at pH 6.5. Placing a
histidine on our binder directly across from a permanently negative residue on EGFR
gives two parts that ignore each other at pH 7.4 and attract at pH 6.5.

Earlier work (explainer 03) chose eight positions on EGFR to build these pairs
against. This project calls them **anchors**, and they are E344, H358, D368, H370,
E391, E400, E421 and E424. The letter is the amino acid and the number is the
position. They sit on one face of the protein, within reach of a single small
binder.

**Where the risk comes in.** EGFR's outer region is 621 residues long. The eight
anchors all sit inside one section of it, residues 310 to 480, which folds into its
own compact blob. A section like that is called a **domain**, and this one is EGFR's
third, so it is called **domain III**.

`analysis/09_trim_target.py` cuts domain III out and hands that 171-residue piece to
the design program as the target, in place of the full 621 residues. The piece is
called the **fragment**. The reason is cost: the design program's run time grows with
the size of the target, and a smaller target buys more attempts for the same money,
which explainer 05 covers.

The consequence is that the design program never sees the other 450 residues. That
creates a specific hazard. The other 450 residues might fold across the very face we
are aiming at. If they do, a binder could be designed against a surface that is
covered over in the real molecule. One part of a molecule covering another part so
that nothing can reach it is called **occlusion**.

Occlusion is dangerous because it is silent. Every score the design program reports
would look correct, because the program cannot see the thing doing the covering. The
binder would then measure as nothing in the laboratory, and no number anywhere would
distinguish that from a binder that was simply designed badly.

---

## 3. Why this had to be measured again

An earlier step, `analysis/06_tethered_occlusion.py`, asked this question and found a
real answer. Of the eight anchors it examined, D458 and D460 sit in a groove between
two domains rather than on an open face, pressed against the domain that follows
ours.

Those are not our current anchors. Step 06 ran when the project was aiming at a
different patch of EGFR, residues 415 to 466. Explainer 03 describes the measurement
that moved the target to the face around H370, and only two of the current eight
anchors, E421 and E424, fall inside the range step 06 examined. The other six —
E344, H358, D368, H370, E391 and E400 — had never been checked for occlusion.
Script 12 checks them.

---

## 4. Script 12: is the face reachable?

### What it measures, and why two measurements

The script makes two different measurements for each anchor, because one of them
alone gives the wrong answer in a case that matters.

The first asks whether anything from outside domain III comes within 4.5 ångströms
of the anchor. Two residues that close are touching, and the measurement answers
"is this residue covered over". It uses the same shared function step 06 used, so
the two steps cannot mean different things by the word.

The second counts how many residues from outside domain III fall within 8 ångströms
and within 12 ångströms. This answers a different question: "can something the size
of a binder get here". A binder is a body roughly 20 to 30 ångströms across, so an
anchor can have nothing touching it and still sit at the bottom of a narrow cleft
that nothing that size can reach into. The first measurement would call such a site
open. This project calls the second measurement **clearance**, which is our own term
rather than standard vocabulary.

Neither 8 nor 12 ångströms is a boundary in the molecule. Nothing changes as an atom
crosses either line, so the counts are reported and read as a description of how
crowded a position is, rather than used to pass or fail an anchor.

The alternative to both of these was to name the domains doing the covering and
measure against those. That was rejected because it would make the answer depend on
exactly where each domain is taken to begin and end, and those boundaries are
something this project is unsure of. Asking what lies "outside domain III" needs only
one boundary, the one `analysis/09` already uses for the cut.

### Which structures, and why two of them

A **structure file** is a list of the measured positions of every atom in a protein,
determined in a laboratory and published. The script uses two:

- **6ARU**, in which EGFR is in its extended shape, with its domains splayed out.
  This is the file `analysis/09` cuts the fragment from.
- **1NQL**, in which EGFR is in its tethered shape, folded back on itself.

EGFR's outer region switches between those two shapes. The competition tests the
whole outer region in solution, and we do not know which shape predominates there,
so an anchor face that is open in one shape and covered in the other would be a
risk. Both are measured for that reason. In each file the measurement uses the EGFR
chain alone, leaving out whatever else was bound in that crystal: 6ARU has an
antibody fragment attached and 1NQL has a signalling molecule attached, and
including them would mostly measure the difference between those two attachments,
which is a circumstance of growing each crystal rather than a description of our
test.

### The result

For each anchor: what touches it at 4.5 Å, how far away the nearest piece of protein
from outside domain III is, and how many such residues fall within 8 and 12 Å.

In 6ARU, the extended shape:

| anchor | touched at 4.5 Å | nearest piece from outside domain III | within 8 Å | within 12 Å |
|---|---|---|---|---|
| E344 | nothing | 24.1 Å | 0 | 0 |
| H358 | nothing | 23.8 Å | 0 | 0 |
| D368 | nothing | 13.0 Å | 0 | 0 |
| H370 | nothing | 17.0 Å | 0 | 0 |
| E391 | nothing | 17.4 Å | 0 | 0 |
| E400 | nothing | 13.8 Å | 0 | 0 |
| E421 | nothing | 4.6 Å | 2 | 9 |
| E424 | nothing | 7.7 Å | 1 | 7 |

In 1NQL, the tethered shape, the same eight: nearest distances 23.6, 19.1, 14.8,
17.4, 19.3, 14.0, 7.1 and 8.5 Å in the same order, with zero residues within 8 Å
for the first six, three for E421 and none for E424.

The row that decides the question is not any single one. What decides it is that the
"touched at 4.5 Å" column reads "nothing" for all eight anchors in both shapes, and
that the "within 8 Å" column reads zero for six of them. Those six are on an open
face. E421 and E424 are clear of contact and have a neighbouring piece close enough
to be worth recording.

### How the measurement is checked

A measurement that has never been seen to fail is not known to work, so the script
is checked two ways.

First, it runs the old anchors from step 06 through the same code as a control. If it
cannot reproduce step 06's finding, its answer about the new anchors is not to be
believed either. It reproduces it: D458 is touched by Q486 at 2.98 Å and D460 by
K487 at 2.77 Å, the same partners at the same distances step 06 reported. That
result is then compared automatically against the table step 06 committed, and the
run stops with an error if the two ever differ. This follows the rule in `CLAUDE.md`
that came out of steps 04 and 06 once disagreeing about the same quantity.

Second, the script can break its own numbering on purpose. Run as
`python analysis/12_fragment_context_check.py --break-numbering 1`, it shifts every
residue number by one place and then checks whether each position holds the amino
acid our numbering says it should. All eight checks fail, reporting among other
things that the position believed to be histidine 370 holds a leucine, and the script
exits with an error. Checking by amino acid identity rather than by arithmetic is
deliberate: arithmetic can be off by one and still look correct, and a histidine that
turns out to be a leucine cannot.

### What it changes about the design

The occlusion risk that cost the old 415–466 patch two of its anchors does not apply
to the face around H370, in either published shape. Cutting domain III out does not
hide the anchors, so `analysis/09`'s decision to hand the design run a 171-residue
fragment is supported for this target. The entry for it has moved into the Settled
section of `docs/decisions-log.md`.

Two things it does not license. It says the anchors are reachable, and says nothing
about whether the separated fragment folds the way the intact protein does, which is
section 7. And it is a measurement of the target, which leaves a gap that section 5
fills.

---

## 5. Script 13: would each individual design fit?

### Why the answer in section 4 is not enough

Section 4 establishes a property of the target: the face is reachable. A binder
could still be designed that approaches that reachable face at an angle, or is large
enough, that part of it ends up where another domain of EGFR already is. The design
program would score it exactly as well as a binder that lies neatly on the face,
because the design program cannot see the other domains.

So section 4 licenses the campaign and `analysis/13_full_receptor_clash.py` screens
the candidates it produces.

### How it places a binder back into the whole receptor

The design program returns each surviving candidate as a **complex**: a structure
file holding the target and the designed binder together, in the positions the
program predicts. The target in that file is the 171-residue fragment.

To ask whether the binder would fit in the intact protein, the script lays the
candidate's fragment on top of the same residues of the full 6ARU structure, then
moves the binder by exactly the same amount. Laying one structure on another by
rotating and sliding it until it matches as closely as possible is called
**superposition**, and one that may only rotate and slide, never bend, is a
**rigid-body superposition**. The result places the designed binder in the
coordinate frame of the whole receptor, where its distance to every other domain can
be measured.

It then counts how many binder atoms come close to the 438 residues the design
program never saw. Contacts with domain III itself are counted separately and are
not a fault: those are the intended join between binder and target, and counting
them serves only as a check that the binder is where it is supposed to be.

Two distances separate the outcomes:

- Below **2.5 ångströms**, two atoms would have to occupy the same space. A carbon
  atom is about 1.7 ångströms in radius and atoms in ordinary contact sit 3.5 to 4.5
  ångströms apart, so anything below 2.5 is an overlap of the atoms themselves. Two
  atoms being in the same place is a different kind of problem from two being close:
  it cannot happen. This is called a **steric clash**.
- Between 2.5 and 4.5 ångströms, the two are touching. A real receptor can shift a
  little to accommodate that, so the script calls this outcome marginal and reports
  it rather than treating it as a failure.

A candidate with no atom within 4.5 ångströms of anything outside domain III is
called clear.

### Why a rigid measurement, and what it cannot see

This asks where things are, which needs no structure prediction and no **GPU**
(graphics processing unit, the specialised hardware a structure predictor needs). So
it runs on every candidate in seconds, instead of competing for the rented machine
the design run needs.

The alternative was to predict each candidate's structure against the whole 621-residue
outer region. That would answer a richer question, including whether EGFR would move
to make room for a binder that overlaps it slightly, and it costs considerably more
per candidate. How much more is UNVERIFIED: a structure predictor's cost is
understood to grow with roughly the square of the chain length, which would put it
near thirteen times the work, but that figure is from memory and is not measured
here. It is item 9 in the decisions log's next actions, to be measured on the first
run. Either way, the richer check suits a final shortlist and this one suits
screening every candidate, so both have a place.

What the rigid measurement cannot see is exactly that movement. This is why a
marginal verdict is reported as marginal and carried forward with its numbers,
rather than being counted as a failure.

### Why it flags rather than discards

A clashing candidate stays on the list, labelled with the overlap measured. The
decision to discard a design belongs in `analysis/10_charge_pair_filter.py`'s
ranking, so that every decision about what to submit is made in one place and can be
reviewed in one place.

The script does refuse to score a candidate in two situations, and refusing differs
from discarding: a refusal means no verdict was reached, and the reason is printed.

- The superposition is poor, above 2.5 ångströms. Domain III is one rigid blob, so a
  good match should come in well under that. A poor match means the movement applied
  to the binder is unreliable, so every distance computed from it would be
  unfounded, and a confident clash figure derived from a bad placement is worse than
  no figure.
- No chain in the file reads as human EGFR, so which molecule is the target cannot
  be established. The design run's output really does contain binder-only files
  alongside the complexes, and this is how they are handled.

Which chain is the target is established by measurement rather than by the chain's
letter. The design program puts the target on chain A and the binder on chain B, and
version 1 of the same program did the opposite, so a script that trusted the letter
would read the wrong molecule and report plausible nonsense. Each chain is scored by
how many of its residues read as the amino acid human EGFR has at that position.

### What it has been run on

No design run has produced anything, so the script builds its own candidates and
checks all four on every run. Each uses the real trimmed fragment as its target, so
the superposition being exercised is the real one, with a compact cluster of 27
alanine residues standing in for a binder. The stand-in carries no designed sequence
because this step measures where a binder's volume is and nothing else about it. The
files are written in the same format the design run writes, so they go through the
same parser real output will.

| constructed candidate | what it is | required verdict | result |
|---|---|---|---|
| binder on open face | placed out from the anchor face, where a real binder would sit | clear | clear, superposition 0.00 Å |
| binder inside domain IV | placed in the space the rest of the receptor occupies | clashing | clashing, 55 atom overlaps |
| target of the wrong shape | a target whose shape is not 6ARU's | refused | refused, superposition 5.98 Å over 171 atoms |
| binder-only file | no target chain present | refused | refused, no chain reads as EGFR |

The script can also switch each of its rules off, and the matching case must then
fail. Switching off the overlap count turns the domain IV case from clashing into
marginal; switching off the superposition limit turns the wrong-shape case from
refused into clear. Both report a failure and exit with an error, so both rules are
known to be doing something.

One thing is worth recording about how the third case was built. It first moved the
whole target 69 ångströms away and expected a refusal. It came back clear, with a
superposition of 0.00 ångströms, which is the correct answer: a superposition exists
in order to undo a rigid movement, so a target moved bodily across the box lands
back in place exactly. Only a change of shape can make a superposition fail. The
case now distorts the target's shape instead, and the reasoning is written into the
script and the glossary so the test is not rebuilt the same way later.

### What it changes about the design

It adds a screen that turns the remaining truncation risk from an assumption into a
measurement made on each candidate. A clear verdict means a design is not ruled out
by where its atoms sit, and says nothing about whether it binds or whether it
switches with pH, which `analysis/10` and the laboratory decide.

Until the design run produces output, this is a prepared check rather than a result.

---

## 6. One change to shared code

Steps 12 and 13 need quantities steps 09 and 10 already compute: the list of eight
anchors, how to read a candidate's structure file, and how to tell the target chain
from the binder chain. Those three now live in `analysis/egfr_common.py` and all four
steps call them.

`CLAUDE.md` requires this. Steps 04 and 06 once computed the same quantity twice and
reported different answers, and a second copy of a calculation is a second chance to
get a different answer from the same input. Steps 09 and 10 keep their own function
names, which now pass the work through to the shared versions, and both still pass
their own tests.

---

## 7. What this does not settle

- **Whether the 171-residue fragment holds the shape its anchors sit on.** This is
  the remaining risk from cutting domain III out, and nothing here touches it.
  Script 12 measures the intact protein, so it cannot speak to what the separated
  piece does. `analysis/11_trim_boundary.py` tried and could not resolve it: all
  three boundary options it tested sat about 5 ångströms from the measured
  structure, which is the noise floor of the predictor it used, so the comparison
  separates nothing. Only the design run's own structure predictor can settle it,
  and the first run on rented hardware can do so by predicting the fragment and
  comparing it against 6ARU across the eight anchors.
- **Script 13 has processed no real output.** Its four results are from constructed
  inputs.
- **Accessibility of the current anchors in the tethered shape has not been
  computed.** Section 4 measures geometry in both shapes. The separate measurement of
  how much of each residue's surface water can reach was computed for the extended
  shape only, by `analysis/03_solvent_accessibility.py`, on the whole receptor chain,
  so it already includes any covering by the rest of the protein in that shape. The
  function to reuse for the tethered shape is `receptor_only_sasa` in
  `analysis/06_tethered_occlusion.py`. Given 14 ångströms or more of clear space
  around six of the anchors in the tethered shape, this is unlikely to change the
  conclusion, which makes it a gap to close rather than a doubt to act on.
- **Both structures are crystals of a protein that moves.** 6ARU and 1NQL are two
  shapes EGFR is known to take, and the balance between them in the competition's
  test remains unknown. The conclusion in section 1 holds for both shapes we can
  measure, which is as far as two measured structures reach.
- **A clear verdict from script 13 is a statement about space.** Whether a design
  binds, and whether it switches with pH, is what `analysis/10` and the experiment
  decide.

---

## Terms introduced here

For `docs/glossary.md`, all four added in the same commit as this document:
clearance, occlusion, rigid-body superposition, steric clash.

Terms used here and defined in earlier explainers: ångström, anchor, amino acid,
binder, complex, domain, EGFR, fragment, GPU, histidine, protein, residue, residue
numbering, structure file, superposition, UniProt.
