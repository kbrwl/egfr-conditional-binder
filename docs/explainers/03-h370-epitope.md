# Explainer 03 — the epitope we are actually aiming at

How the target patch moved from residues 415–466 to a patch built around H370,
what was measured to decide that, and which residues the binder now has to grip.

Written for a reader with no biology background.
Covers `analysis/08_h370_epitope.py`, run 1 October 2026.
Numbers are quoted from `results/findings/08-h370-epitope.md` and
`data/derived/08-h370-neighbourhood.csv`. Where those disagree with this
document, they are right and this one needs updating.

---

## First, the thing that confuses people about this change

H370 is **not** the epitope. An **epitope** is the patch of the target's surface
that a binder actually touches — a group of residues, not one residue. A single
residue is far too small for a protein to grip; it would be like gripping a wall
by one brick.

H370 is the **centre point** we measured outwards from. The epitope is the set of
residues that came back from that measurement: **eight of them**, sitting together
on one face of the protein, with H370 among them.

Those eight are:

> **E344, H358, D368, H370, E391, E400, E421, E424**

Six acidic residues and two histidines. That is the patch. The rest of this
document is how we got to it and why we believe it is better than the old one.

---

## Why H370 at all

The pH switch this whole project is built on needs a **histidine**. Histidine is
the only one of the twenty amino acids whose electrical charge changes between
pH 7.4 (healthy tissue, blood) and pH 6.5 (inside a tumour): neutral above,
positive below. Everything that makes our binder switch on in a tumour and stay
silent elsewhere hangs off that single chemical fact.

Human EGFR's outer region contains many histidines. Guessing which one produces a
real, measurable pH effect would be a guess. We do not have to guess, because
somebody measured it.

Liu et al. 2022 (*Molecular Therapy — Oncolytics* 27:256–269) took EGFR and changed
each histidine, one at a time, into alanine — a small, chemically inert residue.
That technique is called **alanine substitution**, or alanine scanning: you break
one thing at a time and see what stops working. Two histidines turned out to carry
EGFR's pH-dependent antibody binding:

- **H433** — which sits inside our old epitope, and which that paper then built its
  own molecule against
- **H370** — the other one, outside our old epitope, and, as far as we could find,
  not used by any published binder

Both are experimentally supported. Only one is taken. That is the reason H370 was
worth a measurement.

---

## How the candidate set was built, and why not from the sequence

The obvious move would be to pick a window of the sequence around position 370 and
rerun the earlier steps on it. That is what we did for 415–466, and it is the wrong
tool here.

A protein chain folds. Two residues sitting next to each other in the written
sequence can end up pointing in opposite directions once it folds, and a residue
forty positions away in the sequence can end up right beside H370 in space. The
part of EGFR we are aiming at, **domain III**, folds into a **solenoid** — a
repeating spiral, like a spiral staircase — which makes this especially severe.

So a sequence window is a list of letters. An epitope is a patch of surface. They
are different objects.

The candidate set was therefore built from the three-dimensional structure instead:

1. Take the centre of H370's **imidazole ring** — the five-membered ring on
   histidine's side chain where its charge actually sits.
2. Collect every residue whose side-chain tip lies within **25 ångströms** of that
   point. An ångström is a ten-billionth of a metre; 25 of them is roughly how far
   one small designed protein can reach with a single face. That number was fixed
   in advance, before measuring, so it could not be tuned to produce a nice answer.
3. Filter what comes back.

We do also report the sequence run, because it is real and useful: **365–376**,
twelve residues, `ISGDLHILPVAF`. It contains both H370 and D368. See the species
section below for why it matters.

---

## The funnel: 154 → 16 → 8

**154 residues** fall within 25 ångströms of H370. That is just everything nearby,
most of it useless to us. Three filters cut it down.

### Filter 1 — can this residue carry a charge pair?

Our switch works by pairing charges across the interface, so a residue is only
useful if it carries a charge, or can be made to. Three residue types qualify:

- **D** (aspartic acid) and **E** (glutamic acid) — acidic, negatively charged at
  both pH values. Against these we put a histidine on *our* binder: neutral at
  pH 7.4 so nothing happens, positive at pH 6.5 so it is pulled in.
- **H** (histidine) — the switch itself. Against these we put an acidic residue
  on our binder: again, nothing at pH 7.4, attraction at pH 6.5.

Everything else is scenery.

### Filter 2 — is it the same in mouse?

The competition requires the same sequence to bind **both** human and mouse EGFR.
The cheapest way to get that is to aim at a patch where the two species are already
identical, so cross-reactivity falls out of where we aim rather than needing repair
later.

Three nearby candidates were thrown out here, and it is worth naming them so the
exclusion is visible rather than silent:

| Position | Human | Mouse | Why it is out |
|---|---|---|---|
| 393 | D | E | different residue in mouse |
| 383 | H | R | different in mouse, and R cannot switch at all |
| 412 | E | D | different residue in mouse |

### Filter 3 — can water reach it?

A residue buried inside the protein cannot be touched by a binder, whatever its
chemistry. The measurement is **solvent-accessible surface area (SASA)**: roll a
water-sized ball over the structure and record how much of each residue it can
touch. Divide that by the largest value that residue type could ever have and you
get **relative solvent accessibility (RSA)**, a number between 0 and 1 that lets
different residue types be compared. Our working cutoffs: RSA ≥ 0.25 is *exposed*,
≤ 0.05 is *buried*, in between is *partial*.

This filter removed **H418** (RSA 0.032, buried in structure 6ARU), which is the
same verdict step 03 reached and the same ambiguity step 06 raised. H418 is
conserved in mouse and would otherwise have qualified.

### What survived: 16 candidate anchors

**Anchor** is this project's word for a residue on EGFR that the pH switch pairs
against. Sixteen survive all three filters:

| Anchor | Distance from H370 | RSA | Exposure | Cetuximab touches it? |
|---|---|---|---|---|
| H370 | 0.0 Å | 0.134 | partial | — |
| D368 | 5.5 Å | 0.122 | partial | — |
| H433 | 9.7 Å | 0.648 | exposed | **yes** |
| D347 | 12.4 Å | 0.706 | exposed | — |
| E344 | 14.2 Å | 0.480 | exposed | — |
| D379 | 14.4 Å | 0.183 | partial | — |
| H358 | 16.6 Å | 0.500 | exposed | — |
| D460 | 16.9 Å | 0.186 | partial | — |
| D458 | 16.9 Å | 0.248 | partial | — |
| E400 | 17.3 Å | 0.147 | partial | — |
| E391 | 18.9 Å | 0.283 | exposed | — |
| E455 | 21.2 Å | 0.348 | exposed | — |
| D388 | 21.2 Å | 0.426 | exposed | — |
| D416 | 21.5 Å | 0.117 | partial | — |
| E421 | 21.9 Å | 0.210 | partial | — |
| E424 | 22.0 Å | 0.429 | exposed | — |

Sixteen candidates, three target histidines among them (H358, H370, H433), and
exactly one residue that the existing drug cetuximab touches (H433).

---

## Being honest about "7 versus 16"

The decisions log compares 7 candidate anchors at 415–466 against 16 here, and that
comparison needs a caveat the raw numbers do not carry.

The two sets were built by different methods. The 7 came from a 52-residue sequence
window; the 16 came from a 25-ångström sphere in space. The sphere is simply a
bigger and better-shaped region, and it **contains seven of the old anchors** —
D416, E421, E424, H433, E455, D458, D460 all reappear in the list above.

So "7 to 16" is not nine brand-new discoveries. It is a better-defined search area
that keeps the old candidates and finds nine more.

The number that actually improved, and the one the design depends on, is the next
one.

---

## The eight that can be reached at once

Having candidates is not enough. A binder is one object with one face, and that
face has to press against all of its anchors at the same time. So the real question
is: how many of the sixteen can a single face cover?

### The distance test, and why it is not sufficient on its own

The test used through step 05 asks whether a set of anchors fits inside a ball 25
ångströms across. That test has a flaw that only becomes serious as sets get
larger: a binder presents something closer to a flat face than a ball, so two
anchors on *opposite sides* of that ball pass the distance test while being
completely unreachable together.

With four anchors the gap between the test and the thing it is standing in for is
small. With eight it is not. So a second test was added.

### The face test

For each anchor, take the direction it points away from the protein's centre.
Anchors on one face all point roughly the same way; anchors wrapped around the
protein point apart. Measure the widest angle between any two directions in the
set: the **angular spread**. A set under 90° of spread is treated as reachable by
one face.

The 90° limit is a working convention chosen here. It is not a measured property
of real binders, and it is cruder than docking an actual protein backbone against
the surface. Read it as a filter that removes the clearly impossible, not as a
guarantee about what it lets through.

### Result

The largest set that passes both tests:

> **E344, H358, D368, H370, E391, E400, E421, E424**
> span 24.3 Å · angular spread 58° (widest pair: H358 and E424)
> six acidic, two histidine · includes H370

Eight anchors, against a minimum of three. Three or four pairs are what the switch
needs; eight gives redundancy rather than the bare minimum, which matters because
the switch is partial rather than all-or-nothing at pH 6.5 and because some pairs
will be lost when a real foldable backbone is fitted to the surface.

For comparison, the same face test run on the old epitope's clusters, judged on
identical terms:

| Cluster | Span | Spread |
|---|---|---|
| D416, E421, E424, E455 | 22.6 Å | 39° |
| E424, E455, D458, D460 | 22.8 Å | 24° |
| H433, E455, D458, D460 | 23.6 Å | 30° |
| D416, E424, E455, D460 | 24.6 Å | 38° |
| D416, H418, E421, E424, E455 *(only if H418 is usable)* | 24.0 Å | 44° |

Every one of those passes the face test. None of them exceeds five anchors.

---

## Where the two epitopes differ, measure by measure

| | 415–466 | H370 patch |
|---|---|---|
| human/mouse differences in the core run | 1 (S442G) | **0** (365–376) |
| candidate anchors, conserved and reachable | 7 | 16 |
| largest cluster on one face | **4** | **8** |
| target histidines in that cluster | 1 (H433) | **2** (H358, H370) |
| cluster anchors cetuximab touches | 1 of 4 | **0 of 8** |
| anchors in the domain IV groove | 2 (D458, D460) | **0** |
| sugar chain inside the region | yes, at N444 | none within the cluster |
| already used by published work | yes, H433 | no |

The 415–466 column describes the cluster H433, E455, D458, D460 (span 23.6 Å),
which is the one of that epitope's four 4-anchor clusters that contains a target
histidine. Step 05 found that the clusters containing a target histidine are
exactly the ones overlapping cetuximab's footprint, because H433 is the single
anchor cetuximab touches. The other three 4-anchor clusters avoid cetuximab
entirely and contain no histidine, so they cannot drive a pH switch at all. The
comparison uses the histidine-containing cluster because that is the only kind
of cluster this project can build on.

The two rows that decide it are the cluster size and the histidine count. The rest
are the complications that disappear as a side effect.

### Species conservation

The run 365–376 has **no** human/mouse differences at all. It ends exactly where it
does because both residues immediately outside it differ: position 364 is serine in
human and alanine in mouse, position 377 is arginine in human and lysine in mouse.
The run is bounded by the differences, which is a good sign that it is a genuinely
conserved stretch rather than an arbitrary slice.

### Cetuximab overlap

Cetuximab is the approved anti-EGFR antibody whose measured footprint we computed
in step 04. Overlapping it invites the obvious objection that we have redrawn an
existing drug's epitope.

None of the eight chosen anchors is touched by cetuximab. H370 itself is not
touched. Cetuximab's nearest contacts are at positions 373, 374 and 377, which
means the *tail* of the conserved run sits inside the drug's footprint while H370
and D368 sit outside it. The only candidate anchor cetuximab touches is H433, and
H433 is not in the chosen cluster.

### The two complications that went away

Two anchors in the old epitope, D458 and D460, sit in a groove where domain III
meets domain IV, which makes them awkward to reach from outside. And N444, inside
the old block, carries an **N-glycosylation** site — a flexible sugar chain
attached to the protein, which is largely invisible in crystal structures and
blocks more surface than the structure suggests. Neither problem exists in the new
cluster.

---

## What this does not settle

Worth stating plainly, because geometry is persuasive and shallow.

- **Reachability is necessary, not sufficient.** Eight anchors sitting on one face
  says a binder *could* touch them. It says nothing about whether a protein exists
  that folds stably and presents the right partner residues at the right angles.
  That is what the design pipeline has to answer.
- **It says nothing about the detection floor.** The hardest requirement is no
  *detectable* binding at pH 7.4, and we do not know what instrument the organisers
  use or where its floor sits.
- **The 90° convention is ours.** It is cruder than real docking.
- **Exposure comes from one structure.** RSA values here are measured in 6ARU. H418
  was buried in 6ARU and partially exposed in 1NQL, which is a reminder that these
  numbers are one snapshot rather than a property of the protein.
- **415–466 is not deleted.** It stays fully characterised and is the fallback if
  the new patch fails later. Everything in explainer 01 about it still holds.

---

## Terms introduced here

For `docs/glossary.md`: alanine substitution, angular spread, candidate anchor,
domain III, face test, imidazole ring, N-glycosylation, relative solvent
accessibility, solenoid, solvent-accessible surface area.
