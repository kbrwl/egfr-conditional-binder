# Glossary

Every term defined so far in this project, in plain English. Add to it rather
than re-explaining in chat.

---

## Protein basics

**Amino acid** — the building block of a protein. Twenty different kinds exist.
Each has an agreed single-letter code, often not the first letter of its name:
R is arginine, K is lysine, W is tryptophan.

**Protein** — a chain of amino acids that folds into one specific, reproducible
3D shape. The sequence determines the shape; the shape determines what it does.
A protein sequence written out is just a string of letters, e.g. `MKTAYIAKQRQ...`

**Residue** — one amino acid at one position in a chain. "Residue 390" means the
390th link counting from the start. Used interchangeably with "position".

**Domain** — a chunk of a protein that folds into its own self-contained blob.
One protein chain can contain several, strung together like beads.

**Sequence identity** — the percentage of positions where two proteins have the
same amino acid, after lining them up. 87% identity means they agree at 87
positions out of every 100.

**MSA (multiple sequence alignment)** — the same idea as an alignment but with many
sequences at once, usually the same protein from many different species. Useful
because a position that has stayed the same across millions of years of evolution is
usually one the protein cannot afford to change, which tells you it matters. Several
structure-prediction tools take one as input. We have not needed one yet; our
comparison is between exactly two sequences, human and mouse.

**Alignment** — lining up two sequences position by position so you can compare
them, inserting gaps where one has extra residues. The comparison is only
meaningful after alignment; you cannot just compare position 1 to position 1.

**BLOSUM62** — the standard scoring table used during alignment. It knows that
some substitutions are chemically mild (serine to glycine) and others are
drastic (glycine to tryptophan), and scores accordingly.

**FASTA** — the standard plain-text file format for sequences. A header line
starting with `>`, then the letters.

**Signal peptide** — a short leader sequence at the very start of a protein whose
only job is to direct the cell where to send it. It is cut off and discarded once
delivery is done, so it is absent from the finished protein. For human EGFR it is
residues 1–24. This is the entire reason two numbering systems exist for the same
molecule, and the reason every number in this project has to say which one it
means.

**Mature protein** — the protein after the signal peptide has been cut off: what
actually exists on the cell. Counting from the mature protein's first residue
gives numbers 24 lower than counting from the UniProt record. Databases and
structure files disagree about which they use, and neither records the choice,
so the offset has to be measured rather than assumed. See
`analysis/00_numbering_check.py`.

---

## Binding

**Binding** — two protein surfaces fitting against each other and sticking,
held by many individually weak forces: matching bumps and hollows, opposite
charges attracting, water-avoiding patches clumping together, hydrogen bonds.
Thirty weak contacts at once makes a grip that lasts.

**Binder** — any protein designed to stick to a chosen target. What we are
designing.

**Affinity / K_D (dissociation constant)** — how tightly a binder grips,
measured as a concentration, where a lower number means a tighter grip.
Micromolar (µM) is
a weak grip. Nanomolar (nM) is drug-grade, a thousand times tighter.
Picomolar (pM) is a vise.

**Epitope** — the specific patch on the target that a binder touches. Choosing
the epitope is choosing where to aim.

**De novo** — designed from scratch rather than copied or tweaked from an
existing molecule. The competition requires it.

**Zero-shot** — submitted without ever having been tested in a lab and refined.
A separate requirement from de novo, and easy to conflate with it. De novo is
about where the design came from: invented rather than copied. Zero-shot is about
*what happened to it afterwards* (no experimental feedback loop). A design could
be de novo but not zero-shot — invented from scratch, then improved over three
rounds of lab measurement. The competition forbids that: one shot, predicted
only. It is also why the write-up carries so much weight, since nobody has
measurements to show.

**Molecule size categories** — the competition judges these separately, so the
category is a deliberate choice rather than a by-product. Microbinder under 40
amino acids; minibinder 40–100; large binder over 100; plus
nanobody and antibody as their own categories. Overall length limit
10–250.

**Detection floor** — the weakest binding an assay can still see. Anything
weaker is reported as nothing at all. This is what makes the pH-selectivity requirement a
**threshold rather than a ratio**: "no detectable binding at pH 7.4" is a claim
about landing underneath the floor, not about the gap between two numbers. A
binder at 10 nM / 200 nM has a 20-fold ratio and still fails, because 200 nM is
comfortably visible. A binder at 2 µM / nothing has a worse ratio on paper and
passes. Aiming for stronger binding therefore works against us here, which is
the opposite of what every standard design pipeline does by default.

**Antibody** — a Y-shaped protein the immune system makes. The two tips of the Y
are the gripping surfaces; the body can generate billions of variants of just
those tips.

**VH and VL** — the two chains that make up the gripping end of an antibody: VH is
the heavy chain's variable part, VL the light chain's. "Variable" because these are
the regions the immune system varies to produce different grips, while the rest
stays the same. The competition asks for a Fab to be submitted as the two written
together with a colon between them, `{VH}:{VL}`.

**scFv / Fab** — two standard ways of packaging just the gripping part of an
antibody into a smaller molecule. Both are allowed submission formats.

**Ångström (Å)** — a unit of distance used for atoms. One ten-billionth of a
metre. Two atoms within about 4.5 Å are considered in contact. A binder's
contact points need to sit within roughly 25 Å of each other.

---

## Our specific target

**EGFR (epidermal growth factor receptor)** — a protein threaded through the cell
membrane, with part outside the cell and part inside. Works like a doorbell: when
a signalling protein called EGF lands on the outside part, two EGFRs pair up, the
inside part switches on, and the cell is told to divide. In many cancers this is
stuck on, so the divide signal never stops. Blocking it slows the tumour.

**Extracellular region** — the part of EGFR outside the cell, roughly residues
25 to 645. Folds into four sub-blobs labelled domain I, II, III and IV.

**Tethered and extended conformations** — EGFR's extracellular region does not
hold one fixed shape. It folds shut on itself (**tethered**, also called closed
or autoinhibited) or opens out (**extended**). In the tethered form, domain II
folds back and lies against domain III, physically covering part of it. This
matters directly: the competition assays the *whole* extracellular region, so
whichever form dominates in the assay buffer decides what a binder can actually
reach. An epitope that is exposed in an open structure but covered in the closed
one produces a design that is correct on paper and measures as nothing in the
tube. Checked in `analysis/06_tethered_occlusion.py`.

**Domain III** — roughly residues 310–480. One of the two domains that grip EGF,
and the patch both approved antibody drugs target. Folds into a solenoid, a
spiral-staircase shape, which means residues adjacent in the sequence can point
in completely opposite directions.

**Cetuximab (Erbitux)** — an approved antibody for colorectal and head-and-neck
cancer. Clamps onto domain III, blocking the EGF landing site. Does not bind
mouse EGFR, because the residues its tips touch differ between species.

**Panitumumab (Vectibix)** — a second approved antibody grabbing an overlapping
patch on domain III.

---

## The pH mechanism

**pH** — the acidity scale, 0 to 14. Seven is neutral. Each whole unit is a
tenfold change in acidity. Blood sits at 7.4. Tumours drift down to about 6.5
because they burn sugar inefficiently and dump acid, with poor drainage.

**Conditional binding** — binding that depends on conditions. Here: grip at pH
6.5, let go at pH 7.4. The point is sparing healthy tissue. EGFR also sits on
healthy skin and gut cells, and cetuximab hitting those causes a severe
acne-like rash, bad enough that patients sometimes stop treatment.

**Histidine (H)** — the one amino acid whose side chain changes charge in the
window between healthy tissue and tumour: neutral above about pH 6.5, positive
below it. No other standard amino acid switches in that range. This makes it the
only practical building block for a pH switch.

**Acidic residues — aspartic acid (D) and glutamic acid (E)** — carry a negative
charge at both pH values we care about. Place a histidine across from one of
these and nothing happens at pH 7.4 (neutral histidine, no attraction), but at
pH 6.5 the histidine turns positive and snaps onto it. That is the switch.

**Basic residues — arginine (R) and lysine (K)** — carry a positive charge at
both pH values. Useful for plain affinity, useless for the switch, since
nothing about them changes between 7.4 and 6.5.

**Contact pair** — one position on the binder facing one position on the target.
An interface is a patterned surface rather than a uniformly charged one, and
electrostatic attraction falls off fast with distance, so each pair behaves
more or less independently. A positive charge twenty ångströms away barely
affects a given pair. This is why the design rule is positional rather than
global: histidine opposite acidic target residues, acidic residues opposite
target histidines. Both switch on together as pH drops.

**pKa** — the pH at which an amino acid is half charged and half neutral;
its flipping point. Histidine's is around 6.0–6.5, which is the whole reason
it works here. Two consequences worth remembering. First, at pH 6.5 a histidine
is only partly positive, not fully, so the switch is gradual — stack three or
four pairs rather than relying on one. Second, pKa shifts with surroundings:
a histidine next to a negative residue holds its charge more easily and flips at
a higher pH, one next to a positive residue flips lower. The same amino acid
behaves differently depending on its neighbours, which is a large part of why
pH selectivity resists reliable computational prediction.

---

## Structure and tools

**PDB (Protein Data Bank)** — free public archive of experimentally measured 3D
protein structures. A structure file is a plain text table of x, y, z
coordinates for every atom.

**1YY9** — the archive entry for cetuximab bound to EGFR domain III.

**6ARU** — the archive entry used as this project's structure of record: the EGFR
extracellular region bound to a cetuximab Fab mutant. Preferred over 1YY9 because
it contains the whole extracellular region rather than domain III alone, which
lets us also ask whether neighbouring domains cover our epitope. It is also the
entry the competition page itself references.

**CA and CB (alpha carbon, beta carbon)** — two specific atoms named by their
position in a residue. Every amino acid has a CA: it is the atom on the protein's
backbone, the continuous chain running through the whole molecule, so a list of CA
positions is a compact way of describing a protein's overall shape. CB is the first
atom of the side chain, the part that differs between amino acid types and sticks out
from the backbone. Where a side chain is too mobile to locate in a structure, CB is
sometimes used as a stand-in for it, which is approximate because the part that
actually forms a charge pair can be several angstroms further out.

**Heavy atom** — any atom in a structure except hydrogen. Hydrogen is too light
to register in most crystal structures, so it is usually simply absent from the
file. "Heavy-atom distance" therefore means "distance between the atoms we
actually have coordinates for", and is the standard basis for deciding whether
two residues touch.

**Solvent accessibility** — how much of a residue is exposed to surrounding
water. Calculated by rolling a water-sized ball over the protein surface and
measuring how much of each residue it can touch. High means a binder can reach
it; near zero means buried and useless to us.

**SASA (solvent-accessible surface area)** — the measured version of the above,
reported as an area in square ångströms.

**Shrake–Rupley** — the specific published algorithm that performs that
calculation, by scattering test points over a sphere around each atom and
counting how many are not blocked by neighbouring atoms. The implementation used
here is `Bio.PDB.SASA.ShrakeRupley`.

**RSA (relative solvent accessibility)** — raw SASA divided by the largest area
that amino acid type could possibly expose if it were fully unobstructed. The
result runs from about 0 to about 1. This normalisation is necessary because raw
area is not comparable between residue types: tryptophan is a far bigger amino
acid than glycine, so 50 Å² of exposed tryptophan is mostly buried while 50 Å² of
exposed glycine is wide open. Read the relative figure; the raw area on its own
is misleading.

**Tien et al. 2013 maxima** — the reference table of those theoretical maximum
areas, one per amino acid type, from Tien, Meyer, Sydykova, Spielman and Wilke,
"Maximum allowed solvent accessibilites of residues in proteins", PLOS ONE 8(11):
e80635. Several such tables exist and they disagree slightly, so which one was
used has to be stated for a number to be reproducible.

**Exposed / partial / buried** — the conventional RSA bands: exposed at 0.25 and
above, buried at 0.05 and below, partially exposed in between. Worth being clear
that these cutoffs are conventions agreed for convenience rather than physical
constants — nothing changes
in the molecule at 0.25. A residue at 0.24 and one at 0.26 are essentially the
same thing, so values near a boundary should be read as ambiguous rather than
rounded into a verdict.

**Superposition** — rotating and sliding one structure until it sits on top of
another as closely as possible, so the two can be compared. Necessary because
coordinates are arbitrary: the same protein solved twice can be described by
completely different numbers in the file.

**RMSD (root-mean-square deviation)** — the average distance left between
corresponding atoms after superposition, in ångströms. It is the standard measure
of how well two structures match. Low means they agree; a couple of ångströms
across a domain is a good match.

**AlphaFold** — software that predicts a protein's 3D structure from its
sequence alone. Outputs coordinates plus a confidence score.

**ipTM (interface predicted TM-score)** — a confidence score that structure
prediction software reports for a complex of two proteins, running from 0 to 1, where
higher is better. It answers "how sure is the model that these two pieces sit
together the way it has drawn them", as opposed to how sure it is about each piece on
its own. Commonly used as the first filter on designed binders, on the reasoning that
if the software cannot confidently place the binder on the target, the design is
unlikely to work. It is a statement about the model's confidence and not a
measurement of binding.

**pAE interaction (predicted aligned error at the interface)** — a confidence
number from structure prediction. Lower is better: it means the model is
confident the two pieces sit where you claim they do. Used to filter
candidate binders.

**BindCraft** — an open pipeline that takes a target structure plus a list of
target residues and generates candidate binder sequences. Runs on a rented GPU.

**RFdiffusion / ProteinMPNN** — the other standard design route. More control,
more setup difficulty. RFdiffusion proposes 3D shapes; ProteinMPNN works out
which amino acid sequence would fold into a given shape.
