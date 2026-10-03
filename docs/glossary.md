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

**Backbone** — the protein's skeleton: the chain of atoms running through every
residue in order, which defines the shape, before deciding which side chains hang
off it. Design tools invent a backbone first and choose the sequence afterwards,
because the shape is what has to be complementary to the target.

**Residue** — one amino acid at one position in a chain. "Residue 390" means the
390th link counting from the start. Used interchangeably with "position".

**Chain** — one continuous string of amino acids. A structure file can hold several,
labelled A, B, C and so on, because a crystal often contains the target and whatever
was bound to it. 6ARU holds EGFR in one chain and the two halves of an antibody arm
in others, which is why `analysis/02` writes a receptor-only file before measuring
anything: a surface measured with the antibody still attached reads as buried where
the antibody is sitting.

**Side chain** — the part of a residue that hangs off the backbone and makes one
amino acid different from another. The backbone is identical in all twenty; the side
chain is what carries the charge, the water-avoiding bulk or the hydrogen-bonding
group, so it is the part that does the chemistry. It matters for the charge-pair
rule that a side chain can be absent from a structure file while the backbone is
present: the atoms were not resolved, meaning the measurement was not sharp enough
to place them. A histidine with no side chain in the file is a histidine with no
evidence its charge reaches anything, which is why
`analysis/10_charge_pair_filter.py` counts those separately rather than as correct
pairs.

**Fragment** — a piece cut out of a larger protein structure and used on its own.
Here it means the 171 residues of domain III (our 310 to 480) that
`analysis/09_trim_target.py` cuts out of EGFR's outer region and hands to the design
tool, because the tool's run time grows with the size of what it is given. Cutting
creates two problems that the full protein does not have: residue numbers change, and
a disulfide bond can lose one of its two ends.

**Structure file** — a plain text list of the x, y and z coordinates of every atom in
a protein, as measured in a laboratory or predicted by a program. The file records
what each atom is and where it sits, and nothing in it says which numbering
convention its residue numbers follow. Common formats are PDB and mmCIF.

**Complex** — two or more molecules held together, saved as one structure file. Here
it means a designed binder together with the target in the position the design tool
predicts they would sit. Each returned candidate is a complex file, and
`analysis/10_charge_pair_filter.py` works out the pairing between binder and target
residues from its coordinates.

**Domain** — a chunk of a protein that folds into its own self-contained blob.
One protein chain can contain several, strung together like beads.

**Residue numbering** — which number each position in a protein is called. There is
no single answer, which is the problem: a sequence archive numbers from the first
amino acid the gene encodes, while a structure file carries whatever numbering its
depositors chose, often counting from the start of the mature protein after the
leader sequence is cut off, and nothing inside the file records which convention it
follows. This project's numbering is positions in the full human UniProt record
P00533; both structures used here are offset by 24 from it. Getting an offset wrong
raises no error and returns a plausible residue in the wrong place, so the offset is
measured by aligning sequences and then checked by asking what amino acid each
position turns out to be. See `analysis/00_numbering_check.py` and the "Residue
numbering" section of `CLAUDE.md`.

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

**UniProt** — the main public archive of protein sequences. Each protein has an
accession, a short identifier: human EGFR is P00533, mouse EGFR is Q01279. Every
residue number in this project is a position in a UniProt record, counted from the
very first residue of the chain as the cell builds it. Saying "UniProt numbering"
is how you specify which of the two possible counting systems you mean.

**Solenoid** — a fold shape: the chain wraps round and round in a repeating spiral,
like a staircase or a coiled spring, rather than packing into a compact ball.
Domain III of EGFR is one. The consequence that matters for design is that going
one step along the sequence moves you one step around the spiral, so residues that
are neighbours in the sequence can point in opposite directions, and a stretch that
looks like a single patch on paper can be spread across several faces in reality.

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

**Interface** — the pair of surfaces where two proteins touch, counted as the
residues of each that come within a chosen distance of the other. The distance is a
convention rather than a physical boundary, so an interface is only defined once the
cutoff is stated: this project uses 4.5 angstroms between non-hydrogen atoms, the
choice steps 04, 06 and 10 all share, and BindCraft2 reports its own at 4.0. The
looser cutoff gives the larger set, so our interface should contain BindCraft2's
rather than match it, and `analysis/10_charge_pair_filter.py` checks that
relationship rather than expecting the two to agree.

**Affinity / K_D (dissociation constant)** — how tightly a binder grips,
measured as a concentration, where a lower number means a tighter grip.
Micromolar (µM) is
a weak grip. Nanomolar (nM) is drug-grade, a thousand times tighter.
Picomolar (pM) is a vise. The competition's instrument reports K_D over roughly
0.1 nM to 10 µM, flowing the target at a top concentration of 1000 nM, so a grip
near the 10 µM end, where only about 9% of the chip's sites are occupied at the top
concentration, may return no number at all rather than a weak one.

**Anchor** — our own term for a residue on the target that the pH switch is built
against. Not standard vocabulary; if you use it with someone outside the project
you would need to explain it. There are two kinds here: an acidic residue on EGFR,
which gets a histidine opposite it on the binder, and a histidine already on EGFR,
which gets an acidic residue opposite it. Eight were identified from sequence
(D416, E421, E424, E455, D458, D460, H418, H433) before any check of whether they
are reachable.

**Candidate anchor** — an anchor that has passed the three filters that decide
whether it could be used at all: it is a residue type that can carry a charge pair
(acidic, or a histidine), it is identical in human and mouse, and water can reach it.
"Candidate" is the reminder that passing those filters is not the same as being
usable, because the set still has to be reachable by one binder face. The H370
measurement produced 16 candidate anchors and 8 that survive the face test.

**Hotspot** — the residues on the target that a design tool is told to aim at,
given to it as a list. It is the field where this project's epitope work ends up: the
eight anchors from `analysis/08` are what goes in it. The tool's own word, so it is
worth knowing; it means the same thing as the chosen anchors here.

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
weaker is reported as nothing at all. For this competition the instrument flows the
target at a top concentration of 1000 nM and reports K_D over roughly 0.1 nM to
10 µM. "No detectable binding" is the organisers' verdict when a trace can neither
be fitted to give a K_D nor clear 300% above the negative control. The floor cuts
both ways: a design with nothing at pH 7.4 and a strong signal at pH 6.5 ranks
highest, and a design made too weak can fall under the floor at pH 6.5 as well. At a
K_D of 10 µM only about 9% of the chip's sites are occupied at 1000 nM, from
1000 / (1000 + 10000).

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

**Glycoprotein** — a protein with sugar chains attached to it. Most proteins on the
outside of a cell are glycoproteins, EGFR included. The sugars are large, flexible
and branched, and they matter here for a practical reason: a structure file usually
shows only the first one or two sugars of a chain that continues well past them, so
a surface can look more open in the file than it is in reality.

**N-glycosylation** — the specific way those sugar chains are attached: through the
nitrogen atom of an asparagine side chain, which is why it is "N". It happens only
at asparagines in particular sequence contexts, so the attachment points are fixed
and identifiable. That context is the **sequon**, N-X-S/T, which has its own entry.
N444, inside the old 415-466 epitope, is one of them. The 1.44 ångström distance
measured from N444 to the neighbouring sugar is a chemical bond, which is how we
know the chain is attached there rather than passing nearby. The target in the
assay is confirmed glycosylated, expressed in HEK293 cells, and a structure shows
only the first ordered sugar or two, so the real reach of a chain is longer than any
distance measured from a structure.

**Tethered and extended conformations** — EGFR's extracellular region does not
hold one fixed shape. It folds shut on itself (**tethered**, also called closed
or autoinhibited) or opens out (**extended**). In the tethered form, domain II
folds back and lies against domain III, physically covering part of it. This
matters directly: the competition assays the *whole* extracellular region, so
whichever form dominates in the assay buffer decides what a binder can actually
reach. An epitope that is exposed in an open structure but covered in the closed
one produces a design that is correct on paper and measures as nothing in the
tube. Checked in `analysis/06_tethered_occlusion.py`.

The organisers have since answered which form the screen uses: the **tethered**
one, stated twice and recorded in `docs/competition-qa-log.md`. Our structure of
record, 6ARU, is the extended form. `analysis/12_fragment_context_check.py`
measured the eight H370 anchors in both forms and found none of them contacted from
outside domain III in either, and domain III has nearly the same fold in both at
1.08 ångström RMSD, so the epitope choice survives the answer. Whether the fragment
handed to the design run should be cut from the tethered structure instead is open.

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

**Histidine switching** — the established technique of putting histidine residues
into a binding surface so that binding depends on pH. Named in Sarkar et al. 2002.
Roughly twenty years old and the background our own work sits on; we cite it rather
than claim it.

**Histidine scanning** — the standard laboratory method for finding where to put
those histidines: make many variants of a binder, each with a histidine in a
different position, and measure which ones actually switch with pH. Schröter et al.
2015 is the canonical method paper. It is a search rather than a prediction, which is
worth knowing because it is what people do when the prediction is unreliable — and
pH behaviour is unreliable to predict.

**Alanine substitution**, also called **alanine scanning** when it is done across
many positions in turn — replacing a residue with alanine, one of the smallest and
chemically dullest amino acids, to find out whether the original residue mattered. If
the effect you were studying disappears, that residue was involved. This is how Liu
et al. 2022 established that EGFR's H370 and H433 are the histidines responsible for
pH-dependent antibody binding: they removed each histidine in turn and saw which
removal removed the pH effect. It is a direct measurement of involvement, which is
why it is worth more than inferring importance from a structure.

**Recycling antibody** — an antibody engineered to release its target inside the
acidic compartments of a cell, so the target is destroyed while the antibody survives
to be used again. A well-known application of pH-dependent binding running in the
opposite direction to ours: release when acidic, rather than grip when acidic.
Igawa et al. 2010.

**Imidazole ring** — the five-membered ring at the end of a histidine side chain,
containing two nitrogen atoms. It is the part that picks up or loses a positive
charge as pH changes, so it is the actual working component of the pH switch. The
charge sits across the whole ring rather than on one atom, which is why distances
involving a histidine are measured from the ring's centre point rather than from a
single chosen atom.

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

**mmCIF** — the newer structure file format, which has replaced the PDB format
as the archive's own. The same information in a different syntax, without the
PDB format's fixed column widths and so without its limits on how many atoms and
chains a file can hold. It matters here because BindCraft2 writes its designed
complexes as mmCIF rather than PDB, so anything reading its output has to parse
both.

**Disulfide bond** — a covalent link between the sulfur atoms of two cysteine
residues, written as a bond between their SG atoms and measuring about 2.05
angstroms. After the backbone itself it is the strongest thing holding a
protein's shape together, and it can join residues that are far apart in the
sequence. It is why cutting a domain out of a protein is risky: a cysteine whose
partner is outside the cut has nothing to bond to, and an unpaired cysteine can
make the fragment behave differently from the same residues in the intact
protein. `analysis/09_trim_target.py` reports every bond its cut severs.

**SSBOND record** — the line in a structure file where the depositors wrote down
which cysteines they say are bonded. Useful as independent evidence against a
distance calculation: if measuring the sulfur-to-sulfur distances gives a
different set of bonds from the ones the file claims, one of the two is wrong and
worth knowing about.

**Salt bridge** — two oppositely charged groups holding each other, usually taken
to need about 4 angstroms between the charged groups themselves. The thing a
charge pair is trying to form at pH 6.5. Note the difference from a contact: two
residues can be counted as touching on a heavy-atom cutoff while their charged
groups point away from each other and never reach, which is why
`analysis/10_charge_pair_filter.py` reports the charge-group distance separately
rather than assuming a contact is a salt bridge.

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

**Receptor-only structure** — our term for a copy of a structure file with
everything except the target protein deleted: no antibody, no growth factor, no
water, no sugars. Made because accessibility has to be measured on the target as it
would be in the assay. Measuring on the file as deposited would report the patch
under the bound antibody as covered, when it is only covered by something that will
not be present. Written by `analysis/02_structure_prep.py`.

**Rotamer** — one of the specific orientations a side chain can adopt by rotating
about its own single bonds. The backbone holds still; the side chain can swing to
point in noticeably different directions. Only some orientations are comfortable,
so a side chain tends to sit in one of a handful of preferred rotamers. This matters
for H418: a side chain caught in one rotamer can look buried while the same residue
in another rotamer is exposed, which is one explanation for two structures
disagreeing about how accessible it is.

**Superposition** — rotating and sliding one structure until it sits on top of
another as closely as possible, so the two can be compared. Necessary because
coordinates are arbitrary: the same protein solved twice can be described by
completely different numbers in the file.

**Rigid-body superposition** — a superposition that may only rotate and slide a
structure, never bend it. `analysis/13_full_receptor_clash.py` uses one to put a
designed binder back into the full receptor: it lays the candidate's target on the
same residues of the intact structure and moves the binder by the same amount. Note
what this means for testing it — a superposition exists precisely to undo a rigid
movement, so a target shifted bodily across the box lands back in place with an
RMSD of zero, correctly. Only a change of *shape* can make a superposition fail.

**Steric clash** — two atoms closer together than their sizes allow, meaning they
would have to occupy the same space. A carbon atom is about 1.7 ångströms in radius
and heavy atoms in ordinary contact sit 3.5 to 4.5 ångströms apart, so
`analysis/13_full_receptor_clash.py` counts anything below 2.5 ångströms as an
overlap rather than a close approach. Unlike a contact, a clash is not a matter of
degree: two sets of atoms cannot both be there.

**Clearance** — our own term, not standard vocabulary, for how much empty space
sits above a residue. Distinct from accessibility, and the distinction matters: a
residue can have no neighbours within contact distance and still sit at the bottom
of a cleft nothing the size of a protein can reach into. Accessibility asks whether
water can touch a residue; clearance asks whether a binder can. `analysis/12`
measures it as the number of residues from outside domain III falling within 8 and
12 ångströms, alongside the nearest-neighbour distance.

**Occlusion** — one part of a molecule covering another part up, so that a binder
cannot reach it. The risk it creates is silent: a design aimed at an occluded
surface scores perfectly well in a model that cannot see the covering, and then
measures as nothing in the assay, with no number to distinguish that from a badly
designed binder. `analysis/06` found two anchors of the old 415–466 epitope
occluded by domain IV; `analysis/12` found none of the current H370 anchors
occluded, in either measured conformation.

**ESMFold** — a protein structure predictor that works from a single sequence,
with no multiple sequence alignment, meaning no column of related proteins from
other species to read conservation out of. That makes it far cheaper to run than
AlphaFold2 and less accurate in absolute terms. `analysis/11_trim_boundary.py`
uses it through a public web interface because it needs no graphics card or account
and the project's graphics-card access was still blocked. Its absolute error is
unknown here, so the script reads differences between candidates predicted the same
way and not any single number.

**pLDDT (predicted local distance difference test)** — a predictor's own
confidence in each residue's position, from 0 to 100. Above about 70 is usually
read as a reliable backbone and below about 50 as essentially no information; those
are the field's conventions, not something measured here. It is the model's opinion
of itself, so a prediction can be confident and wrong, and it is reported beside the
comparison against the measured structure and never in place of it.

**RMSD (root-mean-square deviation)** — the average distance left between
corresponding atoms after superposition, in ångströms. It is the standard measure
of how well two structures match. Low means they agree; a couple of ångströms
across a domain is a good match.

**Span** — of a set of anchors, the distance between its two furthest-apart members,
in ångströms. It answers whether the whole set fits inside the reach of one small
binder. Our working limit is 25 ångströms.

**Angular spread** — of a set of anchors, the widest angle between any two of them,
measured as the directions they point away from the centre of the protein. Anchors
on one face of a protein point roughly the same way, so their spread is small;
anchors wrapped around opposite sides point apart, so it is large. It is the
measurement the face test is built on.

**Face test** — our own term, not standard vocabulary, for checking that a set of
anchors sits on one face of the protein rather than wrapped around it. It exists
because span on its own is not enough: two anchors on opposite sides of a 25-ångström
ball are within that distance of each other and still impossible for one binder to
touch at the same time. A set whose angular spread is under 90 degrees is treated as
reachable by one face. The 90 degrees is a working convention chosen in this project
rather than a measured property of real binders, and it is cruder than fitting an
actual protein backbone against the surface, so it should be read as a filter that
removes the clearly impossible rather than a guarantee about what it lets through.

**AlphaFold** — software that predicts a protein's 3D structure from its
sequence alone. Outputs coordinates plus a confidence score. **AlphaFold2** is the
version the design pipeline uses, run inside the loop to check each proposed binder
rather than to study a known protein.

**Inverse folding** — working out which amino acid sequence would fold into a shape
you already have. The usual direction is sequence to structure; this runs the other
way, structure to sequence. It is the second step of binder design, once a backbone
has been invented. ProteinMPNN is the standard tool for it.

**Self-consistency check** — feeding a designed sequence into a structure predictor
and asking whether it folds into the shape it was designed as, and lands where it
was meant to. The model is marking its own homework, so it is weak evidence about
whether the binder works in a tube. It is good at catching obvious nonsense, which
is what it is for.

**Trajectory** — one attempt through a design loop: invent a backbone, choose a
sequence, check it, keep or discard. Most trajectories end with nothing, so many are
run to get a few survivors. Run time per trajectory scales with the size of the
whole complex, target plus binder, which is why the target gets trimmed before a run.

**GPU (graphics processing unit)** — the chip built to render game frames, used here
because the same matrix arithmetic that draws graphics is what neural networks are
made of. Structure prediction needs one; a laptop processor is not a substitute. The
cards involved are data-centre hardware, so they are rented by the hour.

**CUDA** — NVIDIA's programming interface for using one of its graphics cards for
general computation rather than graphics. Software that needs a GPU usually means it
needs an NVIDIA card with CUDA, which is why the card choice is not free.

**Serverless** — renting compute that exists only while a job is running. No machine
sits idle costing money between runs, and billing is per second of execution. The
trade is a higher price per hour against paying for nothing in between.

**Modal** — the serverless GPU provider this project uses. A Python file describes
what to run and on what hardware; `modal run` starts a container with that GPU
attached, runs the job, returns the output and shuts down. Chosen over a plain cloud
virtual machine because it removes driver and environment setup, which is the part
that consumes a short time budget.

**Smoke run** — a first test on a known-good example, to prove the whole chain works
end to end before any of your own choices enter it. Authentication, image build, GPU
allocation and file output are all exercised at once. Running your own target first
means a failure cannot be attributed, because a bad epitope, a wrong setting and a
broken install all look the same.

**Design stage** — one named phase of a single design attempt. BindCraft2 runs an
attempt through screen, refine, anneal, harden and mutate, then a final check, and
most stages measure the attempt against a minimum and drop it if it falls short. The
stage an attempt was dropped at, and the metric named alongside it, together are the
reason it failed, which is why `analysis/18_campaign_inventory.py` records both for
every attempt.

**`Target_pLDDT`, and how it differs from `pLDDT.<target>`** — two confidence
readings with similar names that are checked at different times, and conflating them
produced a wrong diagnosis in this project on 2 October 2026.
`pLDDT.<target>` is the model's confidence in the whole predicted complex and is
checked at **every** stage, against 0.6 at screen rising to 0.7 at the end.
`Target_pLDDT` is its confidence in the target's own shape and is checked at the
**final stage only**, against a default of 0.6. A rejection message names whichever
metric stopped the attempt, and that name is what to act on; a setting whose name
merely resembles it may govern a gate the attempt never reached. See also
**design stage**.

**Container** — a sealed, prepackaged computing environment that carries its own
operating system files, libraries and installed software, so a job runs the same way
wherever it is started. The rented machine runs our job inside one. Nothing of ours
is inside it unless we put it there, which is why the campaign's own files have to be
attached explicitly.

**Container image** — the recipe a container is built from, and the stored result of
building it. Ours installs Linux, the graphics-card libraries and BindCraft2 from
source, which takes about 20 GB, so it is built once and reused. Files can either be
baked into the image, which makes editing one force a rebuild, or attached when a
container starts, which does not. This project attaches them.

**Card-hour** — one hour of one rented graphics card. This is the unit the price is
quoted in (an L4 is $0.80 a card-hour) and the unit a campaign's size is worked out
in. Two cards running for one hour and one card running for two hours are both two
card-hours and cost the same, so it measures spend rather than elapsed time.

**Candidates per card-hour** — how many designs clear BindCraft2's own quality
thresholds per hour of rented card. The figure that turns a budget into a campaign
size. It cannot be carried across from one target to another, because run time grows
with the size of the complex being modelled, nor across cards, because their speeds
differ.

**Charge-pair survival rate** — of the designs BindCraft2 accepts, the share that
also reach this project's floor of three correct charge pairs
(`analysis/10_charge_pair_filter.py`). It exists as a separate number because
BindCraft2's thresholds and ours measure different things: BindCraft2 has no pH term
anywhere in what it optimises, so the pH switch is selected for entirely afterwards.
A low rate means we are paying for designs that are confident interfaces and not pH
switches.

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
target residues and generates candidate binder sequences. Runs on a rented GPU. It
bundles the three design steps — invent a backbone, choose a sequence, check it —
into one automated loop, which is why it was chosen over wiring the tools together
by hand. **BindCraft2** is the version this project actually uses, and the
difference matters: v1 depends on PyRosetta and v2 does not mention it anywhere,
which removes a licensing obstacle. See the decisions log for what was verified
about it.

**Campaign** — BindCraft2's word for one design run against one target,
configured by a single JSON file and writing one output folder. Our campaign file
is `design/configs/egfr-domain3-h370.json`, generated by
`analysis/09_trim_target.py`.

**Detargeting, and off-target** — BindCraft2's way of steering a binder away from
something. A target given a negative weight is an **off-target**: the binder is
pushed away from it during design, and a design is rejected if its predicted
interface confidence with the off-target reaches a ceiling, 0.4 by default. The
output tables record `i_pTM_detarget`, which is the same interface confidence
(`ipTM`, above, written `i_pTM` in BindCraft2's tables) measured against the
off-target, with lower the better here because the aim is avoidance. A prediction
that finds no interface with an off-target is weak evidence that none exists.
`analysis/15_histag_counterscreen.py` uses it to avoid the His tag on the target.

**i_pDAE** — the confidence score BindCraft2 ranks its own output by. Its output
documentation describes it as a distance-masked interface TM confidence between 0
and 1, computed over contacts within 8 ångströms, and **higher is better**; its
source confirms the direction, since its own list of lower-is-better metrics
(`i_pAE`, `i_pTM_detarget` and others) does not include it. An earlier version of
this entry said lower was better and expanded the name as an error measure, which
the documentation does not support; the full form of the name is not given there.
It has no equivalent in BindCraft v1. It is worth naming because it is the ordering
this project demotes below the charge-pair count: it measures how sure the model is
about the interface, which is a proxy for binding strength.
`analysis/10_charge_pair_filter.py` reranks by charge-pair count first and uses
i_pDAE only as a tie-break among candidates with equal pair counts, preferring the
higher reading, which means the more confident interface.

**Interface_Target_Residues / Interface_Binder_Residues** — the two columns
BindCraft2 writes listing which residues lie in the interface, each a list of
residue identities and numbers, measured at a 4.0 angstrom cutoff. They are two
separate lists rather than a pairing: BindCraft2 does not record which binder
residue faces which target residue, which is exactly what the charge-pair rule
needs, so step 10 recomputes the pairing from the structure file.

**RFdiffusion / ProteinMPNN** — the other standard design route. More control,
more setup difficulty. RFdiffusion proposes 3D shapes; ProteinMPNN works out
which amino acid sequence would fold into a given shape.

---

## The assay, the construct and the submission

These are the conditions the designs are actually measured in. All of them come
from the organisers rather than from anything this project computed; the
attributions and dates are in `docs/competition-qa-log.md` and the consequences
are worked through in `docs/explainers/07-the-assay-and-what-it-changes.md`.

**SPR (surface plasmon resonance)** — the instrument the competition measures
binding with. One partner is fixed to a metal-coated sensor chip, the other is
washed over it in solution, and the instrument shines light at the underside of
the chip and watches the reflection change as mass builds up on the other side.
It needs no dye and no radioactive label and it watches in real time, so it sees
both the sticking-on and the falling-off.

**BLI (bio-layer interferometry)** — a second instrument measuring the same thing
by a different optical trick. The organisers use it only to cross-check hits that
SPR has already found.

**Ligand and analyte** — in an assay of this kind, the **ligand** is the partner
fixed to the chip and the **analyte** is the partner flowed over it in solution.
Worth stating because the arrangement here runs the opposite way to the obvious
guess: our designed binder is the ligand, stuck to the chip, and EGFR is the
analyte flowed over it. The word ligand is used in a second, unrelated sense
elsewhere in biology, for a small molecule a protein binds; the assay sense is the
one meant here.

**Immobilisation** — fixing a protein to the sensor surface so it cannot move. Our
designs are immobilised at their C-terminal end, the end of the protein chain that
carries the added tail. That end is therefore not free to take part in binding,
which is a design constraint rather than a detail.

**Association** — the rising part of an SPR trace, recorded while the analyte is
flowing on and accumulating. The falling part, after the flow switches to plain
buffer, is dissociation.

**Negative control** — the same measurement run with something known not to bind,
so it records what the instrument reads when nothing is happening. Every other
trace is judged against it. The organisers count a design as binding if its
association signal rises more than 300% above the negative control, in cases where
the trace cannot be fitted to give a K_D.

**No detectable binding** — the organisers' verdict when a trace can neither be
fitted to give a K_D nor clear 300% above the negative control. It is a statement
about what the instrument could see, so it depends on the instrument's range: top
analyte concentration 1000 nM and a reportable window of roughly 0.1 nM to 10 µM.
A binder near the weak end of that window may read as no detectable binding at
pH 6.5 as well as at pH 7.4, which fails human binding rather than demonstrating pH
selectivity. At a K_D of 10 µM, the top of the reportable window, only about 9% of
the chip's sites are occupied at the 1000 nM top concentration. See also
**detection floor**.

**His tag** — a short run of histidine residues, usually six, added to the end of a
protein so it can be caught and purified on a metal column. Both the human and the
mouse target carry one at the C-terminus and the organisers expect to screen
without removing them. This sits directly on our mechanism: our design rule builds
acidic pockets to grip the target's own histidines, and such a pocket will grip the
histidines of a tag as readily, while the tag is a floppy exposed tail and the real
target histidines are held in a fold.

**Twin-Strep tag** — a pair of short peptide tags that stick tightly to a
streptavidin surface, used to catch a protein on a chip. It is how the target is
captured and also the last element of the tail added to our designs.

**GFP11** — the eleventh strand of green fluorescent protein, about sixteen
residues. Alone it does nothing; supplied with the other ten strands it completes
them and the pair glows, which is how the amount of protein on the chip gets
measured. It is a beta strand, a flat extended piece of chain that pairs up
edge-to-edge with other beta strands, so a design whose own edge is an exposed beta
strand has a fused strand dangling beside it.

**Linker** — a short, deliberately floppy run of residues joining two parts of a
construct so that neither constrains the other's shape. The designs carry a
linker–GFP11–linker–twin-Strep tail on the C-terminus.

**HEK293** — a human embryonic kidney cell line, widely used to produce proteins.
It matters here because a protein made in human cells carries human-like sugar
chains, so the target in the assay is glycosylated rather than bare.

**Glycan** — a branched tree of sugars attached to a protein after the chain is
built. Glycans are large, they wave about, and a crystal structure resolves only
the first one or two sugars that happen to hold still, so the area a glycan really
covers is larger than any structure shows. A binding site beside one can look open
in a model and be unreachable in practice.

**Attachment point (of a sugar chain)** — the asparagine a glycan is bonded to,
found at a sequon. Distances to the attachment point do not depend on how much of
the chain a structure happens to show, which is why `analysis/14` measures to it and
not to the visible sugars. A sugar seen within about 2 ångströms of its nitrogen is a
bond, since a nitrogen-to-carbon bond runs about 1.45 ångströms.

**ND2** — the nitrogen atom at the end of an asparagine's side chain, which a sugar
chain is bonded to. It is the atom `analysis/14` measures each anchor's distance to.

**Coldspot** — a residue on the target that a design tool is told to keep clear of,
the opposite of a hotspot. BindCraft2 accepts them. A coldspot discourages contact
with that residue itself and does not model what lies beyond it, such as a mobile
sugar chain that sweeps past its own attachment point.

**Sequon** — the three-residue pattern that tells you where an N-linked glycan will
be attached, written N-X-S/T: asparagine, then almost any residue, then serine or
threonine. Proline in the middle position blocks attachment despite fitting the
pattern. Because it is a pattern in the sequence, a sequon can be found by reading
the sequence alone with no structure. So N-A-T is a sequon and N-A-V is not,
valine being neither serine nor threonine.

**Ionic strength** — a measure of how many charged particles are dissolved in a
solution. The assay buffer sits at roughly 170 mM, matched between the two pH
conditions, which is near the value in blood.

**Screening (electrostatic)** — the cancelling-out of a charged group's pull by the
cloud of oppositely charged ions that collects around it in solution. The higher
the ionic strength, the more screening. The consequence for this project: a charge
pair sitting out in the open, surrounded by water and ions, contributes much less
grip than the same pair buried in the closed-off core of an interface where the
water has been squeezed out. The word is also used loosely to mean testing a
library of candidates, which is a different thing.

**HEPES and MES** — two buffers, substances that hold a solution at a steady pH.
HEPES works around pH 7.4 and MES around pH 6.5, so the assay uses HEPES for the
neutral condition and swaps in MES for the acidic one. Everything else in the
buffer is held the same and the ionic strength is matched, so a difference between
the two conditions comes from the pH rather than from the salt.

**Neutralisation assay** — a follow-up experiment asking whether a binder actually
blocks what the target does, rather than merely sticking to it. The organisers will
run these partly to catch designs that bind an added tag rather than the target.

**MMseqs2** — a fast tool for searching a protein sequence against large databases
of known sequences. The organisers run it to decide whether a design is new,
against SwissProt, the Protein Data Bank, the USPTO and EBI patent databases,
THPdb, PLAbDab and Proteinbase.

**Novelty level** — the organisers' score for how unlike anything existing a
submitted sequence is, on a scale of 4. Level 3 clears the submission gate, and the
score is computed automatically on upload, so a design can be tested against the
real checker before the deadline. One organiser also described the requirement as
being under 30% similar to anything existing; whether those two statements describe
the same threshold is recorded as unresolved in `docs/competition-qa-log.md`.

**Protonation** — gaining a proton, which is a hydrogen atom stripped of its
electron and carries one positive charge. A group that has gained one is
**protonated**. Histidine is the only one of the twenty amino acids that switches
between protonated and not across the pH range this competition measures: mostly
uncharged at pH 7.4, and a useful fraction positively charged at pH 6.5. That
gained charge is what the acidic residues on our binder are placed to grab, and it
is the whole of the switch. See also **pKa** and **histidine switching**.
