# Explainer 01 — the structure check

What `analysis/00` through `analysis/08` did, why, and what came back.

A later document, `02-what-was-already-published.md`, covers the literature search that followed and
withdrew this project's novelty claim. Read this one for the measurements and that
one for what they turned out to be worth.
Written for a reader with no biology background who intends to understand the
work rather than trust it.

Covers work done 30 September and 1 October 2026.
Numbers here are quoted from `results/findings/`; if the two ever disagree, the
findings files are correct and this one needs updating.

---

## How to read this

Each section names the files involved, explains the question in plain language,
gives the technical terms for what was done, reports what came back, and says
what it changed.

The technical terms are given deliberately. Knowing that a measurement is called
"relative solvent accessibility" is how you read other people's work later, ask a
useful question, and recognise the term when a design tool prints it. Full
definitions live in `docs/glossary.md`.

---

## The question all of this was answering

We picked residues 415 to 466 of human EGFR as the patch to aim at. That choice
came from comparing the human and mouse sequences letter by letter: those 52
positions are identical in both species except one.

Comparing letters tells you what each residue is. It cannot tell you which
direction it points. A protein chain folds up, and a residue can end up on the
outside where a binder can reach it, or buried in the middle where nothing can.

Domain III, the part of EGFR we are aiming at, folds into a **solenoid** — a
repeating spiral, like a staircase. In that kind of fold, two residues next to
each other in the sequence often point in opposite directions.

So the whole structure check exists to answer: of the eight anchor residues we
chose from sequence, which ones are actually reachable, and can one binder touch
enough of them at once?

**Anchor** is our word for a residue on EGFR that the pH switch would pair
against. There are two kinds:

- six acidic residues (D416, E421, E424, E455, D458, D460), which carry a
  negative charge and would be paired with a histidine on our binder
- two histidines already on the target (H418, H433), which would be paired with
  an acidic residue on our binder

---

## Step 00 — the counting problem

**Files**
`analysis/00_numbering_check.py` → `results/findings/00-numbering-check.md`,
`data/derived/numbering-check.csv`

**The problem.** Two different counting systems exist for the same protein.

The public sequence archive, **UniProt**, counts from the very start of the chain
as the cell first builds it. But the first 24 amino acids are a **signal
peptide** — a temporary address label telling the cell where to send the protein,
which gets cut off and discarded. The mature protein that sits on the cell
surface begins at UniProt position 25.

The competition publishes the mature region. So their position 1 is our position
25:

    UniProt position = challenge-construct position + 24

Nothing in a sequence file records which system it uses. Read one as the other
and every residue number in the project shifts by 24, with no error raised.

**What the script does.** Checks lengths, confirms the competition's human
sequence is exactly UniProt residues 25–645, confirms position 415 is a threonine
by three independent routes, and confirms all eight anchors read correctly in
both systems and in the mouse sequence. Exits with an error if any check fails.

**Result.** All checks pass. One caveat recorded: the mouse construct is 623
amino acids, two longer, because mouse carries a two-residue insertion near
position 638. The +24 offset is only valid before that point. Our epitope sits
well before it.

---

## Step 01 — where human and mouse differ

**Files**
`analysis/01_alignment.py` → `results/findings/01-alignment.md`,
`data/derived/01-domain3-differences.csv`,
`data/derived/01-identical-runs.csv`

The hand-written interpretation of these numbers is `docs/alignment-findings.md`.
It was written before the scripts existed and lives in `docs/` rather than
`results/findings/`, because that directory holds computed output only.

**What it does.** Lines the two sequences up position by position and lists the
differences. Lining two sequences up is called an **alignment**, and it has to be
done properly rather than by simple side-by-side comparison, because one sequence
may have extra residues that shift everything after them.

The method is **global pairwise alignment** using **BLOSUM62**, a standard
scoring table that knows some swaps are chemically mild (serine for glycine, both
tiny) and others are drastic.

**Result.** Domain III has 16 human/mouse differences. Inside 415–466 there is
exactly one: **S442G**, human serine against mouse glycine.

This step is now also a regression test. We already know the answer, so if it
ever disagrees, the environment is broken and nothing downstream should be
trusted.

---

## Step 02 — reading the 3D structure

**Files**
`analysis/02_structure_prep.py` → `results/findings/02-structure-prep.md`,
`data/derived/02-chain-inventory.csv`, `data/derived/02-numbering-offset.csv`,
`data/structures/6aru_receptor_only.pdb`

**What a structure file is.** The **Protein Data Bank (PDB)** is a public archive
of experimentally measured protein structures. A structure file is a plain text
table listing the x, y and z coordinates of every atom. Our structure of record
is entry **6ARU**: the EGFR extracellular region with a cetuximab **Fab mutant**
stuck to it. A Fab is the gripping end of an antibody, cut away from the rest. The
antibody in this file is a modified cetuximab rather than the approved drug — its
own title says "Fab mutant" — and step 07 works out which residues differ and
whether any of them touches EGFR.

A structure file has **chains** — separate molecules in the same file, each given
a letter. Here, one chain is EGFR and two are the antibody's halves.

**The numbering trap.** Structure files carry their own residue numbering, and
for secreted proteins they usually count from the mature protein. So their
residue 392 is UniProt's 416. Nothing in the file says which convention it uses.

The script does not assume. It pulls the amino acid sequence out of the
structure, aligns it against the reference sequence, and derives the shift from
the data. Then it checks the result against residues whose identity we already
know.

**Result.** The offset is +24, consistent across 100% of 609 matched residues.
All eight anchors read as expected: D at 416, H at 418, and so on. Chain A is
EGFR at 99.7% identity to the human sequence, confirmed by measurement rather
than taken from the competition page. Chains B and C are the antibody's light and
heavy chains, per the file's own `COMPND` annotation.

The script also writes out a **receptor-only** copy with the antibody removed,
which step 03 needs.

---

## Step 03 — which anchors are on the surface

**Files**
`analysis/03_solvent_accessibility.py` →
`results/findings/03-solvent-accessibility.md`,
`data/derived/03-epitope-rsa.csv`, `data/derived/03-anchor-verdicts.csv`

**What is measured.** Roll a ball the size of a water molecule over the outside
of the protein and record, for each residue, how much of its surface that ball
can reach. The result is an area in square ångströms, called **solvent-accessible
surface area (SASA)**. The specific published algorithm used is
**Shrake–Rupley**; the implementation is `Bio.PDB.SASA.ShrakeRupley`.

An **ångström** is a ten-billionth of a metre. For a sense of scale measured from
our own structure: the distance between two bonded carbon atoms in 6ARU averages
1.53 ångströms, across 559 such pairs. That is a bond length — the gap between two
atoms' centres — rather than the width of an atom, which is roughly twice as
much.

**Why raw area is not enough.** Amino acids differ hugely in size. 50 square
ångströms of exposed tryptophan is mostly buried; 50 of exposed glycine is wide
open. So the raw number is divided by the largest area that residue type could
possibly expose, giving **relative solvent accessibility (RSA)**, a fraction from
roughly 0 to 1. RSA is the number to read.

The reference maxima come from **Tien et al. 2013**, a published table of those
per-residue maxima. Several such tables exist and disagree by up to 20%, so which
one was used has to be stated for a number to be reproducible.

**Why the receptor alone.** Cetuximab sits directly on the surface we care about.
Measuring on the complex would report our epitope as buried when it is only
covered by an antibody that will not be anywhere near the assay.

**Result.**

| Anchor | RSA | Reading |
|---|---|---|
| H433 | 0.648 | two-thirds in open air |
| E424 | 0.429 | clearly reachable |
| E455 | 0.348 | clearly reachable |
| D458 | 0.248 | on the boundary |
| E421 | 0.210 | a slice available |
| D460 | 0.186 | a slice available |
| D416 | 0.117 | about a ninth showing |
| H418 | 0.032 | effectively inside the protein |

Seven of eight survive. H418 is lost — and H418 was one of our two target
histidines.

**Two honest qualifications the script records.** The cutoffs (0.25 exposed, 0.05
buried) are conventions the field settled on for convenience; nothing physical
happens at 0.25, so residues within 0.05 of a cutoff are flagged as borderline.
And running the numbers through both published reference tables shows D458
changes category depending on which one you pick, so its status is ambiguous
rather than settled.

**One risk found that was not in the plan.** EGFR is a **glycoprotein**: sugar
chains are attached at specific points, a process called **N-glycosylation**.
N444 is one of those points, and it sits inside our block. The measured 1.44
ångströms from N444 to the sugar is a chemical bond rather than two things
happening to sit close. Since accessibility was computed on the bare protein,
numbers for anchors near N444 are upper bounds. A crystal structure only resolves
the innermost, most ordered sugars; real chains reach further and move.

---

## Step 04 — what the existing drug touches

**Files**
`analysis/04_cetuximab_contacts.py` →
`results/findings/04-cetuximab-contacts.md`,
`data/derived/04-cetuximab-contacts.csv`, `data/derived/04-epitope-overlap.csv`

**The definition.** There is no natural line where one protein stops touching
another, so a convention is needed. The one used: an EGFR residue is in contact
if any of its **heavy atoms** sits within 4.5 ångströms of any heavy atom of
either antibody chain. "Heavy atom" means any atom except hydrogen — hydrogens
are too light to appear in most crystal structures and are simply absent from the
file.

Distances are found with `Bio.PDB.NeighborSearch`, which builds a spatial index.
The obvious alternative is comparing every atom against every other atom, which
for 5,000 atoms is 25 million comparisons. Same answer, far slower.

**What it settled.** A previous working assumption held that four human/mouse
differences — Q390R, E412D, R414W, K467R — explain why cetuximab fails on mouse
EGFR. That assumption was stated from memory and never computed.

**One of the four was right.** K467R is genuinely in contact at 3.39 ångströms.
Q390R, E412D and R414W are not within 4.5 ångströms of the antibody at all.
Recorded as disproved rather than quietly dropped.

The computation found two the assumption missed: **R377K** and **S442G**.

**S442G is the useful one.** It is the single position inside our 52-residue
block where the species differ, and it turns out to sit on a surface a real
therapeutic antibody uses. That produces a hard design rule: **do not let our
binder contact position 442.** A binder that touches it risks behaving
differently in the two species at the one position where they differ, which would
undermine the mouse cross-reactivity objective.

**Overlap with our block.** Cetuximab touches 10 of our 52 residues, and exactly
one of our eight anchors: H433. So we sit adjacent to the drug's surface and
partly share it, without reproducing it.

---

## Step 05 — can one binder reach them all

**Files**
`analysis/05_anchor_geometry.py` → `results/findings/05-anchor-geometry.md`,
`data/derived/05-anchor-distance-matrix.csv`,
`data/derived/05-anchor-clusters.csv`,
`explorer/anchor-viewer.html`, `explorer/anchor-view.pml`

**The question.** A binder is one object presenting one face. Everything it grips
has to fit on that face. So: do the surviving anchors sit close enough together for
one face to reach them all?

The threshold used is 25 ångströms, fixed before the measurement as a working
estimate of how far one small binder can reach. It is a convention chosen in
advance so it could not be tuned to whatever answer came out, not a measured
property of binders. A cluster that only just fits should be read as borderline.

**How distance is measured.** Between the parts of each side chain that actually
form the charge pair — the carboxylate carbon for aspartic acid and glutamic
acid, the centre of the imidazole ring for histidine — rather than between
backbone atoms. The **imidazole ring** is the five-membered ring that makes
histidine's charge switch work.

**Result.** Four anchors fit within 25 ångströms, spanning 22.6. The epitope
survives.

Four different groups of four qualify, and they are not equivalent:

| Group | Span | Target histidine? | Overlaps cetuximab? | Near the N444 sugar? | In the domain IV groove? |
|---|---|---|---|---|---|
| D416, E421, E424, E455 | 22.6 Å | none | none | D416, E421 | none |
| E424, E455, D458, D460 | 22.8 Å | none | none | none | D458, D460 |
| H433, E455, D458, D460 | 23.6 Å | H433 | H433 | H433 | D458, D460 |
| D416, E424, E455, D460 | 24.6 Å | none | none | D416 | D460 |
| *D416, H418, E421, E424, E455* | *24.0 Å* | *H418* | *none* | *D416, E421* | *none* |

The last row is only available if H418 turns out to be usable. See step 06.

**The bind.** With H418 excluded, H433 is the only remaining target histidine, so
it is the only place left to use the half of the pairing rule that puts an acidic
residue on the binder. It is also the single anchor cetuximab touches. So that
choice is the one most open to a "that is just cetuximab's epitope" objection.

There is a counter-argument: cetuximab grips H433 with no pH dependence at all,
so sharing one residue with it is not sharing a mechanism. That is an argument to
make in the write-up. Nothing computes it, and it must not be presented as a
result.

A separate point, found later and more important: **this arrangement is published,
and on H433 specifically.** Liu et al. 2022 mapped EGFR's own H370 and H433 as the
residues responsible for pH-dependent antibody binding, then deliberately put an
acidic residue against H433. See "What the prior-art check found" below.

**What this does not establish.** Geometric reachability is necessary and not
sufficient. It does not show that a foldable binder exists which presents the
right partner residues in the right orientations, nor that the resulting switch
is large enough to clear the assay's detection floor at pH 7.4.

---

## Step 06 — EGFR has two shapes

**Files**
`analysis/06_tethered_occlusion.py` →
`results/findings/06-tethered-occlusion.md`,
`data/derived/06-conformation-comparison.csv`,
`data/derived/06-intra-chain-occlusion.csv`

**The worry.** The outer part of EGFR is not a fixed shape. It folds two ways:
**extended** (open, domains splayed out) and **tethered** (closed, folded back on
itself and self-inhibited, with domain II pressed against domain III). These two
shapes are called **conformations**.

The competition assays the whole outer region floating in buffer. If our patch is
covered in the closed shape, a perfectly good binder measures as nothing, and we
would have no way to tell that apart from a bad design.

**The difficulty.** You cannot watch a protein move. You get frozen snapshots.
So the test compares one open structure (6ARU) against one closed structure
(1NQL), measuring how exposed our 52 residues are in each.

That comparison only means something if three things hold.

**One — strip the partners.** 6ARU has an antibody attached; 1NQL has EGF
attached. Measuring with those in place would mostly measure the difference
between an antibody and a growth factor. Both are measured on the bare receptor
chain.

**Two — is domain III itself folded the same way in both?** If the domain had
rearranged internally, an exposure difference could be the domain reshaping
rather than its surroundings changing.

The test is **superposition**: lay one copy onto the other, rotating and sliding
until it fits as closely as possible, then measure the leftover gap. It uses one
atom per residue, the **alpha carbon (CA)** — the central carbon every amino acid
has, which traces the backbone.

The leftover gap is reported as **RMSD (root-mean-square deviation)**: take the
gap between each matched pair, square them, average, take the square root. 171
alpha carbons matched, RMSD **1.08 ångströms**. For scale, measured from 6ARU: a
glycine spans about 3.0 ångströms at its widest, an alanine about 3.4, and
consecutive alpha carbons sit about 3.8 apart. So the leftover gap is roughly a
third of the width of a small amino acid — comfortably less than one residue, which
is what "same fold" means here.

**Three — are the two structures actually in different overall shapes?** If both
happened to be open, the comparison would say nothing about the closed form and a
reassuring answer would be empty. With domain III pinned in place, how far does
the rest of the molecule sit from its counterpart?

| Region | Gap after superposing domain III |
|---|---|
| Domains I–II (25–309) | 23.44 Å |
| Domain III (310–480) | 1.08 Å |
| Domain IV (481–620) | 2.12 Å |

The front half swings about 23 ångströms between the two. Genuinely different
arrangements, so the comparison is real.

### Result 1 — the original worry is answered

Average RSA across our 52 residues is 0.173 open and 0.189 closed. Slightly more
exposed when closed. All eight anchors stay reachable. Domain II folding across
our face does not happen.

### Result 2 — domain IV, which nobody was looking for

The step also asked a separate question: which residues in 415–466 have another
part of *the same receptor chain*, from outside domain III, within 4.5 ångströms?
That is "covered by another bit of itself", phrased so the answer does not depend
on where domain boundaries are drawn.

Answer: 19 residues in the open structure, 18 in the closed one, spanning roughly
446 to 466. The partners doing the covering sit at positions 481–524, which is
**domain IV**, the domain immediately after ours.

Nearly the same list in both structures at nearly the same distances, so this is
a standing feature of how the protein folds rather than something closing
introduces. The tethering question stays answered.

What changes is the shape of our surface. **D458 and D460 sit in a groove between
domain III and domain IV rather than on an open face.** Their exposure numbers do
not change — step 03 measured on the whole chain with domain IV already present,
so its effect was always included. What is new is that a binder reaching those
two has to fit into a cleft, which is harder to design than a flat patch, and it
makes those contacts sensitive to any shift in how the two domains sit against
each other.

### Result 3 — H418 comes back

Step 03 called H418 buried at 0.032, reading only the open structure. Step 06
finds it at 0.200 in the closed structure, partially exposed. Domain III has the
same fold in both, so the difference is not the domain rearranging.

The likely explanation: 6ARU has cetuximab clamped onto it, and the antibody may
be physically pinning that side chain down. The burial could be an artefact of
the antibody rather than a property of the receptor on its own.

**H418 is ambiguous rather than dead.** This matters more than anything else
outstanding, because H418 is a target histidine, and the five-anchor group it
enables — D416, H418, E421, E424, E455 at 24.0 ångströms — avoids the groove,
avoids cetuximab's footprint, and carries a target histidine. No four-anchor
option manages all three. H418 and H433 are 25.2 ångströms apart, just over the
reach cutoff, so it is one histidine or the other and never both.

Those two figures are computed in section 7 of `results/findings/05-anchor-geometry.md`,
which runs the clustering twice: once with H418 excluded, as step 03 has it, and
once with H418 included, as step 06 suggests it might be. Both sets are written to
`data/derived/05-anchor-clusters.csv` and
`data/derived/05-anchor-clusters-with-h418.csv`. The second set is conditional on
the unresolved question and labelled that way.

D458 also disagrees between the structures: 0.248 open against 0.446 closed.

### Limits

- Two snapshots cannot give the proportion of open to closed in the assay buffer,
  which is the number that would actually matter.
- 1NQL was crystallised at low pH and with EGF bound. Both are conditions for
  getting a crystal to form rather than descriptions of the assay, though the low
  pH is worth noticing given pH is the variable the whole design turns on.
- The two structures differ in construct, resolution and crystallisation
  conditions as well as conformation, so some of the measured difference comes
  from those.
- Exposure was computed on bare protein and ignores the sugar chain at N444.

---

## Step 07 — the antibody in our structure is not quite cetuximab

**Files**
`analysis/07_fab_mutant_check.py` → `results/findings/07-fab-mutant-check.md`,
`data/derived/07-fab-differences.csv`

**The problem.** 6ARU's own title calls the antibody a cetuximab Fab **mutant**.
Step 04 measured "cetuximab's" contacts using that file, and three results rest on
those contacts: the disproof of the earlier species-failure claim, the rule not to
contact position 442, and the overlap figures in step 05's cluster table. If a
modified residue sits in the interface, those describe a modified antibody rather
than the drug.

**Why it had to be computed.** There was nothing to look it up in. The file's
`SEQADV` records — the place a depositor lists differences from a reference
sequence — exist only for the receptor chain, and show two conflicts at UniProt 540
and 634 plus a six-histidine purification tag, all outside domain III. There are
none for either antibody chain. The Protein Data Bank entry names the mutant
without listing substitutions. And the primary citation is "To Be Published"
(Christie M., Christ D., deposited 2017, released 2018), so there is no paper.

So the mutations were found by comparison against **1YY9**, the reference structure
of cetuximab (Li S. et al., 2005, *Cancer Cell* 7:301–311).

**Result: five differences, one of them in the interface.**

| Chain | Position | 1YY9 | 6ARU | In the interface? |
|---|---|---|---|---|
| Light | 52 | S | D | no |
| Light | 56 | S | D | no |
| Heavy | 28 | S | D | no |
| Heavy | 31 | N | D | **yes — 3.93 Å from H433** |
| Heavy | 216 | R | K | no |

Four of the five replace a serine or asparagine with **aspartic acid**, adding
negative charge across the gripping region. The purpose is unrecorded, because the
structure is unpublished.

**The follow-up that settled it.** Since heavy chain 31 touches H433, and H433 is
one of our anchors, the question became whether unmodified cetuximab touches the
same place. Recomputing the footprint on 1YY9 gives **exactly the same ten residues
inside 415–466**, H433 included, at 3.36 ångströms against 3.47 in the mutant.

So step 04's footprint describes cetuximab and not merely this variant. The 442
rule stands, and so do the overlap figures.

**One thing worth noticing.** Heavy chain position 31 in 6ARU is an aspartic acid
sitting 3.93 ångströms from H433 — an acidic residue on the binder positioned
against a histidine on the target. Whether it was put there for pH-dependent binding
is unknown, since the structure is unpublished. Either way it is a precedent for the
arrangement, and it turned out to be the first of several: see below.

---

## The defect found on 1 October, and what it changed

**Files**
`analysis/egfr_common.py` (new), plus corrections in `06_tethered_occlusion.py` and
a refactor of `04_cetuximab_contacts.py` onto the shared helper, so both callers now
run the same contact calculation

Steps 04 and 06 both computed which residues the antibody touches, from the same
file with the same 4.5 ångström cutoff, and produced different answers.

**The cause.** Step 06 held the numbering as two plain dictionaries, one pointing
each way, and its contact section looked a structure-file number up in the
dictionary that expects UniProt numbers. The two number ranges overlap, so the
lookup succeeded and returned a position 48 away from the correct one — the
24-residue offset applied in the wrong direction, doubling. No error was raised.

**Which was right.** Step 04. Recomputing both directions reproduces step 04's
list exactly. The residue 2.74 ångströms from the antibody is S464, not D416.
D416 is not touched by cetuximab, and step 05's cluster table is numerically
unchanged.

**What the fix exposed.** With the numbering corrected, step 06's count of
epitope residues covered from outside domain III went from 1 to 19 — which is
the domain IV finding described above. That result only appeared because the bug
was fixed.

**The prevention, now a standing rule in `CLAUDE.md`.** Anything computed twice
lives in one function in `analysis/egfr_common.py` and both callers use it. The
two dictionaries are replaced by a `Numbering` class whose methods are named
`uniprot_of` and `pdb_of`, so handing it the wrong kind of number returns nothing
instead of a plausible wrong answer. Step 06 now compares its result against step
04's committed table and stops the run on disagreement. That guard was tested by
corrupting the input on purpose and confirming it fails.

---

## What the prior-art check found

Before claiming the method was new, we searched for whether it already existed. It
does. This is the most consequential thing in this document.

**Putting histidines on a binder to make binding pH-dependent is about twenty years
old.** It has a name in the literature — histidine switching, after Sarkar et al.
2002 — a standard method for finding the right positions (histidine scanning,
Schröter et al. 2015), clinical-stage molecules built on it, and review articles
covering it. It is background to cite, not an idea to claim.

**Putting an acidic residue on the binder against a histidine already on the target
is also published**, and is described as a known strategy in a 2024 review (Wei &
Sulea, *mAbs* 16:2404064). It has been done deliberately, with crystal structures, on
CTLA-4 (Lee et al. 2022) and on VISTA (Johnston et al. 2019, taken to a clinical
antibody by Thisted et al. 2024). Designing a de novo interface around a target
histidine's two charge states was stated as a principle by the Baker lab in 2014
(Strauch et al., *PNAS* 111:675–680).

**And the closest paper is on our target, our histidine and our pH pair.** Liu et al.
2022 (*Molecular Therapy – Oncolytics* 27:256–269) built a pH-dependent anti-EGFR
antibody. They found EGFR's own H370 and H433 were responsible for the pH-dependence
by mutating each histidine to alanine, then deliberately changed their antibody's
Tyr32 to glutamate or aspartate to pair with H433. Binding at pH 6.5 against 7.4, and
human/mouse cross-reactive as well. There is a patent family over it.

So the novelty claim is withdrawn. It was a memory-based belief, exactly the kind
this project's rules say must be checked, and checking it was worth more than
keeping it.

**One result in that paper helps us.** They also tried Tyr32**histidine** — a
histidine on the binder facing the histidine on the target. It changed nothing, while
the acidic versions worked well. That is published experimental evidence for the
rejection criterion this project worked out from first principles: never put a
histidine opposite H418 or H433. The rule now rests on a measurement rather than an
argument.

**And it pointed at two anchors we had missed.** Because that paper named H370, every
histidine in domain III was listed from our own sequences: H358, H370, H383, H418,
H433. H358 and H370 are identical in human and mouse and were never considered,
because both fall outside 415–466. H383 differs between the species, so it is no use
for cross-reactivity. Whether H358 or H370 deserves an epitope of its own is open,
and H370 has experimental evidence behind it.

**What is left that is defensible.** No published example was found of a **de novo
designed miniprotein** that is pH-conditional against EGFR — the EGFR molecules in
the literature are antibodies, and de novo pH-sensitive design has been published on
other targets (Ahn et al. 2025). So the honest claim is about modality and
combination, phrased as "we found no published example" rather than "this has never
been done", and citing Liu et al. 2022 openly rather than leaving a reader to find
it.

---

## Where this leaves the project

**Established.** Seven of eight anchors are reachable. Four groups of four sit
close enough for one binder. A hard rule not to contact position 442. Cetuximab's
real footprint is known, and checked against the unmodified antibody.

The tethering risk is reduced without being eliminated: the epitope is accessible in
both published shapes, and the specific fear about domain II folding across it is
not what the numbers show, but two crystal snapshots cannot give the balance of open
to closed in the assay buffer, which is the number that would actually matter.

**Withdrawn.** The claim that the method was novel, entirely rather than narrowed.
See `02-what-was-already-published.md`, which covers the search, what it found and what it gave back.

**Superseded.** 415-466 is no longer the primary epitope. The prior-art search named
EGFR's H370 as the other experimentally implicated histidine, and evaluating it
produced a better epitope on every measure we have: eight conserved reachable anchors
on one face instead of four, two target histidines instead of one, no overlap with
cetuximab's footprint, and none of the domain IV groove or sugar-chain complications.
Computed in `analysis/08_h370_epitope.py`; comparison table in
`docs/decisions-log.md`. 415-466 stays fully characterised and is the fallback.

**The open question that matters most.** Whether H418 is usable. It decides
between a four-anchor group that carries at least one drawback whichever you
pick, and a five-anchor group that carries none of them. It is settled by
checking further EGFR structures and by examining **rotamers** — the different
orientations a side chain can adopt around its rotatable bonds — rather than by
preferring whichever reading is more convenient.

Full current state is in `docs/decisions-log.md`.

---

## Terms introduced here

All defined in `docs/glossary.md`: alignment, alpha carbon, ångström, anchor,
BLOSUM62, chain, conformation, epitope, extended, Fab, glycoprotein, heavy atom,
imidazole ring, N-glycosylation, PDB, receptor-only, relative solvent
accessibility, RMSD, rotamer, signal peptide, solenoid, solvent-accessible
surface area, Shrake–Rupley, superposition, tethered, Tien et al. 2013, UniProt.

Checked against the glossary on 1 October 2026: all present.
