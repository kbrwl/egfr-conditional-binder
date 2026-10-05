# The assay we are actually designing for, and what it changes

**This explainer covers work with no results of its own yet.** Every number in it
is external, stated by the organisers in the Proteinbase Slack between 28
September and 2 October 2026, and recorded with attribution in
`docs/competition-qa-log.md`. Nothing here was computed by this project. Three
checks that these facts demand are written down at the end and have not been run.
This file is rewritten when they have been.

The point of reading it: until now this project designed against the competition
brief, which says what the objectives are. It said almost nothing about how the
objectives would be measured. The answers to that turn out to move four of our
own decisions, one of them in the opposite direction from where we had it.

**Four words carried over from earlier explainers,** repeated here so this file
can be read on its own. A **residue** is one amino acid at one position in a
protein chain, so "residue 370" means the 370th one counting from the start. An
**epitope** is the patch on the target that a binder grips. An **anchor** is this
project's word for one of the eight charged residues in that patch we are aiming
at: E344, H358, D368, H370, E391, E400, E421 and E424. And an **angstrom**,
written Å, is a ten-billionth of a metre, which is roughly the width of one atom;
two residues touching sit about 4 Å apart.

**The three objectives** are referred to below by name rather than by number,
which `CLAUDE.md` requires. They are **human binding**, binding the extracellular
region of human EGFR; **mouse cross-reactivity**, the same sequence also binding
mouse EGFR; and **pH selectivity**, binding at pH 6.5 with no detectable binding
at pH 7.4. The repository numbers them two different ways — the decisions log
lists them as a design task and `CLAUDE.md` lists them in the organisers' order of
importance, which puts pH selectivity first — so a number alone is ambiguous here.

---

## 1. The experiment, in order

Our designed binder is made as a protein with a tail added to its C-terminal end.
The **C-terminus** is one of the two ends of a protein chain; the other is the
N-terminus. The tail is a linker, then **GFP11**, then a **twin-Strep tag**, each
separated by a short flexible linker.

That tail is how the binder gets stuck to a sensor chip. So our binder is the
part that sits still.

EGFR is then washed over the chip in solution. It is produced in **HEK293 cells**,
a human cell line, which means it carries human-like sugar chains. It is in its
**tethered** shape, the closed inactive one. It carries a C-terminal **His tag**
and is caught on the surface, where needed, by a C-terminal twin-Strep tag.

The instrument is **surface plasmon resonance**, abbreviated SPR. It shines light
at the underside of a metal-coated chip and watches how the reflection changes as
mass piles up on the other side. More protein stuck to the surface, bigger
change. It needs no dye and no radioactive label, and it watches in real time, so
it sees both the sticking-on and the falling-off.

This is run twice: once in a buffer at pH 7.4 and once at pH 6.5. Everything else
is held the same, including the salt.

---

## 2. What the instrument can and cannot see

The highest concentration of EGFR they flow is **1000 nM**, and the affinity
range they can report runs from about **0.1 nM to 10 µM**.

**Affinity** is reported as a **KD**, the dissociation constant. It is the
concentration at which half the binder on the chip is occupied. Smaller number,
tighter grip. 1 nM is a strong therapeutic-grade binder. 10 µM is two molecules
that brush past each other.

A design counts as **binding** at a pH if its trace can be fitted to give a KD,
or, failing that, if the rising part of the trace clears 300% of the **negative
control**. A negative control is the same measurement run with something known
not to bind, so it shows what the instrument reads when nothing is happening;
everything else is judged against it. If neither test holds, the organisers write
down **no detectable binding**.

The practical floor matters more than the ceiling here, and it is arithmetic. The
fraction of binder on the chip that is occupied at any moment is the analyte
concentration divided by itself plus the K_D. At the top concentration of 1000 nM
that gives roughly 91% occupancy for a K_D of 100 nM, 50% for 1 µM, 25% for 3 µM
and 9% for 10 µM. So a design in the high nanomolar range produces a strong,
easily fitted signal, and one at the far end of the reportable range produces a
small one. Where exactly the instrument stops returning a fittable curve depends
on the size of the design and the density on the chip as well as the occupancy,
and the organisers have not said. What the arithmetic does establish is the
direction: deliberate weakening walks a design toward the end of the range where
the signal is smallest, and far enough along it the design reads as no detectable
binding at pH 6.5 as well as at pH 7.4, which fails human binding.

The table is computed by `analysis/17_occupancy_table.py` rather than typed in by
hand, since the same figures also appear in `CLAUDE.md` and `docs/decisions-log.md`
and three hand-typed copies of one calculation is exactly the risk `CLAUDE.md` warns
against. See `results/findings/17-occupancy.md`.

---

## 3. The decision this reverses

`CLAUDE.md` carries a rule this project set early, that affinity should be
deliberately marginal rather than maximal. The reasoning was that pH selectivity
asks for no detectable binding at pH 7.4, so a strong binder risks being strong
at both pHs, while a weak one that clearly switches passes. `analysis/10` is
built to that rule: it ranks by charge pairs alone, never discards a candidate
for binding weakly, and carries BindCraft2's own confidence ordering through
untouched so it can be used to pick *downward*.

Asked on 30 September whether no binding at 7.4 is a hard requirement or whether
a large KD shift also counts, the competition lead at Adaptyv answered that a large
KD shift also qualifies, and that designs showing no binding at pH 7.4 **together
with high affinity at 6.5** rank higher.

So the quantity being rewarded is the size of the gap between the two pH
conditions, with the pH 6.5 end as high as it can be. Our rule aimed at the right
objective by the wrong route: we treated weakness as the way to guarantee silence
at pH 7.4, when what the organisers reward is the gap itself. A deliberately
weakened design gives up ranking at the top and risks falling off the bottom of
the instrument's range at the same time.

**What changes.** The rule in `CLAUDE.md` is withdrawn and replaced: aim for the
largest achievable difference between pH 6.5 and pH 7.4, with affinity at 6.5 as
high as the switch allows. `analysis/10`'s downward use of the confidence
ordering is removed. Ranking stays led by charge-pair count, because the charge
pairs are what produce the switch, with the confidence ordering used upward as
the tie-break rather than downward.

**What does not change.** The pairing rule itself, the rejection of a binder
histidine facing a target histidine, and the exclusion of position 442 all stand.
Nothing in the organisers' answers touches them.

---

## 4. Salt, and where the charge pairs should sit

The buffer is 10 mM HEPES, 150 mM NaCl, 0.2% Tween-20 and 3 mM EDTA at pH 7.4,
with the HEPES swapped for MES at pH 6.5. Total **ionic strength** is matched
between the two conditions at roughly 170 mM using extra sodium chloride.

Ionic strength is a measure of how many charged particles are dissolved in the
solution. Every charged group on a protein attracts a loose cloud of oppositely
charged ions from the solution, and that cloud cancels out much of the group's
pull on anything further away. The effect is called **screening**.

The matching is good news: it means any difference we see between the two pH
conditions comes from the pH rather than from the salt, which is exactly the
comparison our mechanism is built on.

The level is the consideration. At 150 mM, which is roughly the salt
concentration of blood, screening is substantial.

A **salt bridge** is what this project's charge pairs are: one positively charged
group and one negatively charged group close enough to hold each other, usually
within about 4 Å. A salt bridge formed between two charges sitting out in the
open, both still fully surrounded by water and by dissolved ions, contributes
modestly to the grip. The same pair contributes far more when it sits in the
closed-off core of an interface where water has been squeezed out, because there
is no cloud of ions in there to cancel it. As a rough picture: two magnets held
apart in a bucket of iron filings pull on each other far less than the same two
magnets in clean air.

This is a design preference rather than a measurement, and it is stated here as
such: among candidates with equal charge-pair counts, prefer the ones whose pairs
sit in the buried part of the interface rather than on the solvent-exposed rim.
Whether `analysis/10` can measure that distinction is an open question, listed at
the end.

---

## 5. The His tag, which sits directly on our mechanism

Both the human and the mouse target carry a C-terminal **His tag**, a short run
of histidine residues, usually six, added so the protein can be caught and
purified on a metal column. The competition lead said on 30 September that cleaving the tags is
awkward and the screen will probably use the His-tagged protein as it is.

an Anthropic organiser of Anthropic put the problem plainly: a binder that targets
the tag would look pH-selective and would in fact bind anything carrying a His
tag.

This project is more exposed to that than most entries, for a structural reason.
Our positional rule puts acidic residues on the binder facing the target's own
histidines, H358 and H370.

The word doing the work here is **protonated**, and it is the hinge of the whole
mechanism, so it is worth stating plainly. A proton is a hydrogen atom stripped of
its electron, and it carries one positive charge. Protonated means a group has
picked one up. Histidine is the only one of the twenty amino acids that switches
between picking one up and letting it go right around the pH range we care about:
at pH 7.4 most copies of a histidine are uncharged, and at pH 6.5 a useful
fraction have gained a proton and become positively charged. That gained charge is
what our acidic residues are there to grab, and it is why the binder sticks at the
lower pH and lets go at the higher one.

The problem is that an acidic pocket shaped to hold a protonated histidine will
hold the histidines of a tag just as readily. The tag is also easier for the
pocket to reach: it is a floppy unstructured tail free to thread in, where H358
and H370 are held in place by the fold around them.

The organisers' countermeasure is to lean harder on the in silico scoring in
cases that look like tag binding, since the predicted structure has no tag in it,
and to run follow-up neutralisation assays.

**What this adds.** A counter-screen, written before any candidate exists, in the
same spirit as `analysis/10`. Model each surviving candidate against a free
six-histidine peptide and reject or demote any candidate whose acidic pocket
accepts it. Two things come out of that. A design that passes is less likely to
be a tag binder dressed up as a switch, and the methods write-up can say we
checked, which matters when the selection step is reading our documentation.

---

## 6. The tethered form, and what our structure of record is

Asked which of EGFR's two shapes to design against, the competition lead answered on 29
September: the tethered form, the closed inactive one. He restated it on 30
September, saying the wet-lab screen is against glycosylated, tethered EGFR, and
that this is why domain III was suggested.

This project's structure of record is **6ARU**, which is EGFR in the extended
open form, held in a complex with the antibody cetuximab. The fragment
`analysis/09` cuts and hands to the design run comes from 6ARU. So does the
surface-exposure reading in `analysis/03`, the cetuximab contact set in
`analysis/04`, and the comparison in `analysis/11`.

Two measurements already in hand keep the damage small. `analysis/12` measured
the H370 anchor face in both 6ARU and **1NQL**, which is the tethered form, and
found none of the eight anchors contacted from outside domain III in either. And
domain III itself has nearly the same fold in both, recorded in the decisions log
at 1.08 Å **RMSD**.

RMSD, root-mean-square deviation, is the average distance between two structures'
matching atoms once one has been laid on top of the other as well as it will go.
Small number means the same shape. 1.08 Å is close to the width of a single atom,
so domain III is doing nearly the same thing in the open and closed forms even
though the rest of the receptor is not. The epitope call survives.

What does not automatically survive is the fragment we hand the design run. If
the protein in the assay is tethered, the model the designs are shaped against
should be the tethered one, and the cheapest way to settle it is to cut the same
310–480 fragment out of 1NQL and compare the two cut fragments directly.

**One side effect worth recording.** H418 reads as buried in 6ARU, with a
**relative solvent accessibility** of 0.032, and partially exposed in 1NQL, at
0.200.

Relative solvent accessibility, abbreviated RSA, is how much of a residue's
surface is reachable by water, given as a fraction of how much would be reachable
if that residue were floating free on its own. It runs from 0 to 1. Around 0.03
means almost entirely covered over by the rest of the protein, so a binder could
not touch it. Around 0.20 means a fifth of it is exposed, which is enough to
reach. The gap between the two readings is therefore the difference between an
unusable anchor and a usable one.
The decisions log already flags that the exclusion rested on the 6ARU reading
alone and that 6ARU being an antibody complex may be the cause. If the tethered
form is what gets tested, the 1NQL reading is the relevant one. This does not
affect the H370 epitope we are designing against; it strengthens the 415–466
fallback if we ever need it.

---

## 7. Glycans are real, and nothing has checked the current anchors

The target is expressed in HEK293 cells and is glycosylated. This is confirmed,
not assumed.

A **glycan** is a branched chain of sugars attached to the protein after it is
built. **N-linked** glycans attach to an asparagine that sits inside a
three-residue pattern called a **sequon**, written N-X-S/T: asparagine, then
almost any residue, then serine or threonine. A sequon can therefore be found by
reading the sequence alone, with no structure needed.

Worked through, using the single-letter codes for amino acids: a stretch reading
N-A-T is a sequon, because it is asparagine, then anything, then threonine, so a
sugar chain is expected on that asparagine. A stretch reading N-A-V is not, because
valine is neither serine nor threonine. N-P-T is not either, because proline in
the middle position blocks the attachment even though the pattern otherwise fits.

Glycans are large, they wave about, and a crystal structure shows only the first
one or two sugars that happen to hold still. So the real territory a glycan
covers is bigger than any structure shows, and a binding site next to one may be
unreachable in practice while looking open in the model.

`analysis/03` measured this for the old 415–466 epitope and found N444 sitting
inside it, 1.44 Å from a visible sugar. At that distance the two are chemically
bonded, since a bond between two atoms runs about 1.4 to 2 Å and anything merely
sitting nearby would be 3 Å or more away. **That measurement was never repeated
for the H370 anchors.** Step 08 moved the face, and nothing since has asked which
sequons lie near the current eight positions.

Of the checks this explainer asks for, this is the one to run first: it needs no
graphics card and no design output, and it could remove an anchor we are about to
spend the whole campaign budget aiming at.

---

## 8. How the designs get selected, and what that means for the write-up

The competition lead, 2 October: selection is mainly method novelty, design diversity, and a
couple of in silico confidence metrics. He added that this is deliberately less
biased toward any single metric than previous competitions, which selected on
interface confidence scores alone, and that the selection prompts will be
published afterwards. On 30 September he had already said that in silico scoring
is only a smaller component, and that weak confidence scores on a rigid structure
or an obstructed epitope would not be penalised much, if at all.

For Track 3, where up to 20 designs go into a pool and a Claude workflow reads
the accompanying documentation, three things follow.

**The write-up carries most of the weight.** This is already the plan, and it is
now the plan for a stated reason rather than an assumed one.

**Twenty near-identical designs waste the allocation.** Diversity is named as a
selection axis in its own right. A submission should spread across anchor subsets
and across size categories rather than filling all twenty slots with the top of
one ranked list.

**Our confidence scores matter less than we assumed.** The project has been
treating interface confidence as something to be defended. It is one input among
several and the organisers have said so.

One more piece of mechanics changes how the last two days should be spent:
submissions are allowed once every 24 hours, with up to 20 designs each time, and
novelty is scored automatically on upload. An early safe batch can be replaced by
a better one. There is no reason to hold everything for a single upload.

---

## 9. The novelty gate

Level 3 of 4 clears the submission gate, stated twice by the novelty lead at Adaptyv.
On 1 October he also described the requirement as the sequence being under 30%
similar to anything existing. Those two statements may not describe the same
threshold and the blog at `adaptyvbio.com/blog/novelty` should be read before
either is relied on. Recorded as unresolved in `docs/competition-qa-log.md`.

The check itself is MMseqs2, a fast sequence-similarity search, run against
SwissProt, the Protein Data Bank, the USPTO and EBI patent databases, THPdb,
PLAbDab and Proteinbase. Our own novelty step should mirror that list rather than
invent one.

---

## 10. What this does not settle

- Nothing here is measured by us. All of it is reported speech from the
  organisers, and the competition page overrides it wherever they disagree.
- Whether any sequon sits near the eight H370 anchors. Unmeasured.
- Whether the 310–480 fragment cut from the tethered structure differs usefully
  from the one cut from 6ARU. Unmeasured.
- Whether an acidic pocket built for H358 and H370 also accepts a free
  six-histidine peptide. Unmeasured, and the thing most likely to invalidate a
  design that otherwise looks good.
- Whether `analysis/10` can tell a buried charge pair from a solvent-exposed one
  well enough to rank on it. Open question, not yet a plan.
- The exact mouse construct used in the screen, which was asked twice in Slack
  and never answered.
- Whether cynomolgus monkey cross-reactivity is in scope for this challenge. The
  pre-launch announcement said mouse and cyno; the recorded objectives say human
  and mouse.

---

## Terms introduced here

For `docs/glossary.md`: SPR, BLI, analyte, ligand (in the assay sense), KD,
association, immobilisation, His tag, twin-Strep tag, GFP11, linker, HEK293,
glycan, N-linked glycosylation, sequon, ionic strength, screening, MES, HEPES,
tethered form, extended form, MMseqs2, novelty level, neutralisation assay,
negative control, protonation, salt bridge, RMSD, relative solvent accessibility,
angstrom.
