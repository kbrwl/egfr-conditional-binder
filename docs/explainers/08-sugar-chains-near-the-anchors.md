# Explainer 08 — sugar chains near the anchors

Whether the target in the assay carries sugar chains close to the eight positions we
are aiming at, and what that does to the design. The script is
`analysis/14_glycan_sequons.py`.

Written for a reader with no biology or software background. It stands on its own:
every term is explained where it first appears.

---

## 1. The result

All eight target positions are within reach of a place where a sugar chain attaches
to EGFR, and five of them are close to one such place.

The attachment point is an asparagine at position 352, written N352. In the tethered
structure, which the organisers say is the form tested, five of the eight positions
lie between 10.4 and 14.3 ångströms from it: E344 at 12.0, H358 at 10.4, D368 at
14.3, H370 at 10.6 and E391 at 11.3. A sugar is seen bonded to N352 in both of the
structures we hold. The other three positions, E400, E421 and E424, lie 16.2 to
16.9 ångströms from a second attachment point, N361.

An **ångström**, written Å, is a ten-billionth of a metre. An atom is roughly one to
two ångströms across, and the small protein we are designing is about 20 to 30
ångströms across, so 10 to 14 Å is within reach of a sugar chain that sticks out
from the surface.

Two further facts matter more than the distances.

**N361 exists in human EGFR and not in mouse.** Mouse has a tyrosine at that position
(Y) where human has the asparagine, which removes the attachment point. This
difference was already on our list: step 01 recorded N361Y among the 16 differences
between human and mouse in domain III. It was recorded as a change of one amino acid
and its effect on sugar chains was never drawn. The human target carries a chain at
N361 that the mouse target lacks, 16 to 17 Å from E400, E421 and E424. The
competition scores mouse cross-reactivity as a ratio of the mouse to the human
binding strength that should be about 1, so a difference between the two targets
near our anchors is a risk to that ratio.

**The two histidine anchors are among the closest.** H358 and H370 are the target's
own histidines, the positions our pairing rule is built around, and they are the
two nearest anchors to N352 at 10.4 and 10.6 Å.

What this changes about the design: it adds a risk we had not measured, and it does
not remove any anchor. A sugar chain is large and mobile, and what the measurement
shows is that a chain *could* cover part of each position some of the time, which is
different from blocking it. Section 6 sets out what remains to be decided and by
whom.

---

## 2. Why this had to be measured

**EGFR** (epidermal growth factor receptor) is the protein on the surface of cells
that this competition targets. A **protein** is a chain of building blocks called
**amino acids**, and a **residue** is one amino acid at one position in the chain,
so "residue 352" means the 352nd counting from the start. Every residue number in
this project counts from the start of the full human sequence in **UniProt**, the
public archive of protein sequences. The eight **anchors** are E344, H358, D368,
H370, E391, E400, E421 and E424: the positions on EGFR our binder is designed to pair
against. The letter is the amino acid and the number is the position.

On 30 September the organisers confirmed that the target in the assay is made in
**HEK293 cells**, a human cell line, and is **glycosylated**
(`docs/competition-qa-log.md`). A **glycan** is a branched tree of sugars attached to
a protein after it is built. Cells make their proteins carrying glycans, and a
protein made in a human cell line carries human-like ones.

Glycans matter to a binder because they are large and they move. A crystal
structure, which is a list of the measured positions of every atom in a protein,
shows only the first sugar or two of a chain, because only those hold still. The real
chain is longer than the structure shows, and it sweeps through the space around its
attachment point. A binder aimed beside one can find the site covered in the real
molecule while the model shows it open, and the design program cannot see the
difference, since the fragment we hand it contains no sugars.

An earlier step, `analysis/03_solvent_accessibility.py`, checked this for the
original target patch, residues 415 to 466. It found N444 inside that patch, bonded to
a sugar 1.44 Å away. Step 08 then moved the target to the face around H370, and the
sugar check was never repeated. Six of the eight current anchors had never been
measured against an attachment point.

---

## 3. What was measured

### Finding the attachment points by reading the sequence

N-linked glycans attach to the asparagine in a three-residue pattern called a
**sequon**, written N-X-S/T. It is an asparagine (N), then any amino acid except
proline (X), then a serine (S) or threonine (T). The pattern can be found by reading
the sequence alone, with no structure needed.

Worked through: N-A-T is a sequon, because the third residue is a threonine. N-A-V is
not, because valine is neither serine nor threonine. N-P-T is not either, because
proline in the middle position blocks the attachment even though the rest of the
pattern fits. Two can overlap: in N-N-S-T both the first and the second asparagine
are followed by a pattern that fits.

The human extracellular region, residues 25 to 645, contains 11 sequons: N128, N175,
N196, N352, N361, N413, N444, N528, N568, N603 and N623. A sequon is a necessary
condition for attachment and does not show that a chain is attached. A sugar seen
bonded to the asparagine in a structure is evidence that it is. The pattern N-X-C
also supports attachment, at a much lower rate, and is not counted, which makes the
count of 11 a lower bound.

### Measuring from the attachment nitrogen, and not from the visible sugar

Step 03 measured each position's distance to the sugar atoms that are actually
present in the structure file, with a 5 Å cutoff, and found none within reach. That
measurement answers whether an anchor touches a sugar we have coordinates for. It
cannot answer whether the full chain could reach it, because the structure shows only
the innermost sugars of a chain that continues.

This step measures to **ND2**, the nitrogen at the end of the asparagine's side
chain, which is the atom the chain is bonded to. The distance does not depend on how
much of the chain happens to be visible. The two measurements answer different
questions and do not contradict each other.

Each anchor is measured from the charged tip of its side chain, because that is the
atom that forms a charge pair, using the same rule as `analysis/05_anchor_geometry.py`.
The rule now lives in `analysis/egfr_common.py` so the two steps cannot measure from
different points.

### Two structures

The distances were measured in both structures we hold: **6ARU**, EGFR in its
extended shape with its domains splayed out, which is our structure of record, and
**1NQL**, EGFR in its tethered shape, folded back on itself. The organisers say the
screen uses the tethered shape, and a sequon belonging to another domain could sit
closer to an anchor in one shape than the other.

### The bands

The distances are read against two bands that step 05 set for the same purpose:

- under 15 Å: likely shadowed at least some of the time
- under 25 Å: within reach of an extended chain
- beyond that: probably clear

These bands rest on a reach of 20 to 30 Å for a complex glycan, which is from
memory and UNVERIFIED here. They are cautious on purpose, and "within reach" means
covered some of the time and not blocked. Nothing in the molecule changes at 15 or 25
Å, so a value close to either is the same thing as one just across it.

---

## 4. The result in full

Distance from each anchor to the nearest attachment nitrogen, and which one:

| anchor | 6ARU, extended | 1NQL, tethered | band, using the closer |
|---|---|---|---|
| E344 | N352, 12.4 Å | N352, 12.0 Å | likely shadowed |
| H358 | N352, 10.8 Å | N352, 10.4 Å | likely shadowed |
| D368 | N352, 15.3 Å | N352, 14.3 Å | likely shadowed |
| H370 | N352, 11.1 Å | N352, 10.6 Å | likely shadowed |
| E391 | N352, 12.0 Å | N352, 11.3 Å | likely shadowed |
| E400 | N361, 17.2 Å | N361, 16.5 Å | within reach |
| E421 | N413, 16.3 Å | N361, 16.2 Å | within reach |
| E424 | N361, 17.6 Å | N361, 16.9 Å | within reach |

No anchor is probably clear. The row worth reading with care is D368, which sits at
14.3 Å in the tethered structure and 15.3 Å in the extended one, either side of the
15 Å line. Its band depends on which structure is used and it should be read as
borderline.

The tethered structure gives the shorter distance for every anchor. The differences
between the two structures are small, at most 1.1 Å (D368, 15.3 against 14.3), and are
within what structures at 3.2 Å and 2.8 Å resolution can resolve.

### Is a sugar really there?

A sugar within 2.0 Å of the attachment nitrogen is a bond, since a bond between a
nitrogen and a carbon runs about 1.45 Å.

| sequon | 6ARU | 1NQL |
|---|---|---|
| N352 | bonded, 1.44 Å | bonded, 1.43 Å |
| N361 | bonded, 1.44 Å | bonded, 1.36 Å |
| N413 | bonded, 1.44 Å | not bonded, 4.97 Å |
| N444 | bonded, 1.44 Å | bonded, 1.45 Å |
| N528, N568, N603 | not bonded | bonded |

The two sequons nearest our anchors, N352 and N361, carry a visible sugar in both
structures. N413 carries one in 6ARU and not in 1NQL. That is a statement about which
sugars were ordered enough to be resolved in each crystal, and not about whether a
chain is attached, which is the same reason a missing sugar is weak evidence of an
empty site.

### Human against mouse

Of the 11 sequons, 10 are present in both species after the two sequences are lined
up. The exception is N361, which reads NCT in human and YCT in mouse. The sequences
were aligned before comparing, because a single inserted residue shifts every
position after it; there is no gap anywhere near position 361, so the alignment there
is unambiguous.

---

## 5. How the measurement was checked

A measurement that has never been seen to fail is not known to work, so the script is
checked in four ways.

**It reproduces earlier steps.** The old 415 to 466 anchors go through the same code.
N444 to its nearest sugar atom comes out at 1.44 Å, the figure step 03 committed. The
seven old-anchor distances to the N444 side-chain carbon agree with the table step 05
committed to within 0.1 Å. The run stops with an error if either differs.

**The sequon rule is tested on constructed sequences,** including N-A-V, N-P-T and
overlapping patterns, and each part of the rule is switched off in turn on purpose.
In each case the matching test fails. The test that checks the human and mouse
comparison is broken the same way.

**The numbering is confirmed by amino acid and not by arithmetic.** Every asparagine
the script finds, and every anchor, is checked to be that amino acid in each
structure. With the numbering shifted by one on purpose, the checks fail: the position
believed to be H370 reads as a leucine.

**Step 05's outputs did not change** when its functional-point rule and the two bands
moved into the shared module, confirmed by comparing every file it writes.

---

## 6. What this changes, and what is left to decide

It changes what we know and does not change the design on its own. Nothing was
dropped and `analysis/10_charge_pair_filter.py` was not changed.

The work order that asked for this check expected a few anchors to be flagged and
asked for step 10's ranking to demote them, in the way it already demotes a design
that leans on E424 near the cut edge. All eight are flagged, including both
histidines. A demotion by band would reorder every candidate, and it would favour
pairs on E400, E421 and E424, which sit near the human-only chain at N361. That
swaps one risk for another. It is a judgement for the project owner and not a
default to set quietly.

Options that exist, with what each costs:

1. **Disclose and proceed.** Record the risk in the write-up. The organisers
   recommended domain III knowing the target is glycosylated, which is some reason to
   think the region is usable. It costs nothing and removes nothing.
2. **Demote by band in step 10.** Prefers E400, E421 and E424 over the five anchors
   near N352, at the price of leaning on the three anchors near the human-only chain.
3. **Mark the footprint with coldspots.** BindCraft2 accepts coldspots, residues to
   keep clear, and a participant in the Proteinbase Slack suggested treating glycan
   sites that way. A coldspot discourages contact with those residues themselves, and
   a chain's sweep reaches beyond them, so this reduces direct contact and does not
   model the chain.
4. **Design against human and mouse together.** BindCraft2 has a multi-target mode for
   one binder that engages two orthologs, with hotspots given separately in each
   structure's numbering. It is aimed at the mouse objective and could also bear on N361.
   It needs a structure of mouse domain III, which we do not hold.

---

## 7. What this does not settle

- **Whether a chain covers any anchor.** Distance to an attachment point is a
  necessary condition for covering and not evidence of it. A mobile chain covers a
  site some of the time. Only an experiment settles this.
- **The 15 and 25 Å bands.** They rest on a reach figure from memory. The ranking
  among anchors, 10.4 Å against 14.3 Å, is more reliable than the labels.
- **Which sequons are occupied in the assay.** Expression in HEK293 cells makes
  occupancy likely at the sequons that are occupied in the crystals, and it is not
  measured for this construct.
- **The exact mouse construct.** The organisers were asked twice for the residue range
  and vendor and did not answer. The N361Y finding is from the reference mouse
  sequence on the competition page.
- **Sequons that are not N-X-S/T.** The rarer N-X-C pattern is not counted.
- **Anything about the shape the glycan takes.** No sugar chain was modelled.

---

## Terms introduced here

For `docs/glossary.md`: glycan, N-glycosylation, sequon and HEK293 are already
there from earlier work or from explainer 07. New in this explainer, and defined in
the same commit: ND2, attachment point, coldspot.

Terms used here and defined in earlier explainers: ångström, amino acid, anchor,
EGFR, protein, residue, structure file, tethered and extended conformations,
UniProt.
