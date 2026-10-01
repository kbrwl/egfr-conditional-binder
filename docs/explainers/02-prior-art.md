# Explainer 02 — the prior-art check, and what it cost us

What we searched for, what we found, what we withdrew, and why the withdrawal is
the most useful thing in this repository.

Written for a reader with no biology background. Covers work done 1 October 2026.
Numbers and citations here are traceable to `docs/decisions-log.md` under Resolved
and Ruled out.

---

## The short version

This project believed it had found a method nobody else was using. It had not. The
idea is about twenty years old, it is a named technique with its own literature, and
the single closest paper did it **on our target, using the same histidine, at the same
two pH values, with the same cross-species requirement**, in 2022, and patented it.

We found that out by looking, before claiming anything. That is the point of this
document.

---

## Why look at all

The competition is judged partly on method novelty, by people who work in protein
design. A novelty claim that a reader can refute from memory is worse than making no
claim: it costs the credibility of everything else in the submission, including the
parts that are sound.

The project's own working rule says that anything asserted from memory is labelled
UNVERIFIED and must not be built on until a computation or a check resolves it. The
novelty claim was exactly that — a belief, never checked. So it got checked.

---

## What we were claiming

Two things, one of which we knew was ordinary and one of which we thought was ours.

**The ordinary half.** Put **histidine** residues on the designed binder. Histidine
is the only one of the twenty amino acids whose electrical charge changes between
pH 7.4 (blood, healthy tissue) and pH 6.5 (inside a tumour): neutral above, positive
below. So a contact that depends on a histidine being positive is a contact that
switches on in a tumour and off in healthy tissue. We never thought this was new.

**The half we thought was ours.** Put an **acidic** residue — aspartate or glutamate,
which carry a negative charge at both pH values — on the binder, positioned opposite a
histidine that is **already present on the target protein**. The switch works in the
same direction: at pH 7.4 the target's histidine is neutral and there is no
attraction; at pH 6.5 it becomes positive and is pulled towards our negative residue.

The appeal was that this reasons about the **target's** chemistry rather than the
binder's. We believed standard tools do not do that.

---

## What the search found

Three questions, asked separately. All three came back against us.

### Is histidine engineering for pH-dependent binding established? Yes, thoroughly.

It has a name in the literature — **histidine switching**, after Sarkar et al. 2002 —
and a standard laboratory method for finding the right positions, called **histidine
scanning**: make many variants with histidines in different places and measure which
ones switch (Schröter et al. 2015). There is computational work on placing them for
tumour pH specifically (Sulea et al. 2020), clinical-stage molecules built on the
principle, and review articles covering the field.

There is also a well-known application in a different direction: **recycling
antibodies**, engineered to let go of their target inside the cell's acidic
compartments so the antibody survives to work again (Igawa et al. 2010).

This is background to cite. It was never ours.

### Is the narrower idea published? Yes, repeatedly.

Pairing an acidic binder residue against a histidine already on the target is
described as a **known strategy** in a 2024 review (Wei & Sulea, *mAbs* 16:2404064).
It has been done deliberately, with crystal structures to confirm the arrangement, on
at least two other targets:

- **CTLA-4** (Lee et al. 2022): engineered aspartate and glutamate residues on the
  antibody, placed against the single histidine of the target's outer region.
- **VISTA** (Johnston et al. 2019, taken to a clinical antibody by Thisted et al.
  2024): a target whose outer surface carries many histidines, with antibodies
  engineered to bind only when those are protonated.

And designing a from-scratch interface around a target histidine's two charge states
was stated as a design principle by the Baker laboratory in 2014 (Strauch, Fleishman
& Baker, *PNAS* 111:675–680), though in the opposite direction — their binder
released at acidic pH rather than gripping.

### Is there published pH-selective work on EGFR? Yes, and this is the one that hurts.

**Liu X et al. (2022), "A cross-reactive pH-dependent EGFR antibody with improved
tumor selectivity and penetration obtained by structure-guided engineering",
*Molecular Therapy – Oncolytics* 27:256–269.**

They did the following:

1. Found which of EGFR's own histidines actually carry pH-dependent binding, by
   changing each one to **alanine** — a small, inert amino acid — and seeing which
   change removed the effect. Two came out: **H370 and H433**.
2. Deliberately changed a residue on their antibody (Tyr32) to **glutamate or
   aspartate**, to form a new acidic-to-histidine pair with **H433**.
3. Measured binding at **pH 6.5 against pH 7.4**.
4. Confirmed the antibody binds **both human and mouse** EGFR.
5. Patented it (family WO2024109709A1).

That is this project's strategy — same target, same histidine, same pH pair, same
mechanism, same cross-species objective — published four years ago.

---

## What we withdrew

The novelty claim, entirely. Not narrowed, not relocated.

In particular we do **not** argue that designing a small protein rather than an
antibody makes the method novel. That carve-out would be transparent to anyone who
knows the field, and the people judging this do know it.

The old strategy section said we were "competing where no established method exists".
It is kept in the decisions log under **Ruled out**, quoted in full, with the reasons
it was withdrawn. A withdrawn strategy is a result and is recorded like any other
result, rather than quietly deleted.

Two other sentences went at the same time, for a related reason. "Affinity
optimisation is a solved-ish problem that well-resourced labs will win" and "nobody
can compute it" were both claims about **what other people do and can do**. We have no
way to know either. A claim about a *tool* is checkable from its documentation; a claim
about *people* is not. Those were deleted rather than re-sourced.

---

## What the search gave back

Three things, and they are worth more than the claim we lost.

### 1. Our rejection criterion now rests on a measurement

This project worked out from first principles that you must **never** put a histidine
on the binder facing a histidine on the target: at pH 6.5 both become positive and
repel, so that pair switches *off* as pH drops and can cancel a correctly built pair
elsewhere.

Liu et al. tested exactly that arrangement. They made the Tyr32→**histidine** variant
alongside the acidic ones. It changed neither the pH-dependency nor the affinity,
while the acidic versions improved pH-dependency substantially.

So a rule we had derived by argument is now a rule supported by a published
experiment. That is a strict improvement.

### 2. It found two anchors we had missed, and one of them is better than what we had

Because that paper named H370, we listed every histidine in domain III from our own
sequences rather than taking it on trust:

| Histidine | Same in mouse? | In our old epitope? | Status |
|---|---|---|---|
| H358 | yes | no | missed; never considered |
| H370 | yes | no | missed; **experimentally implicated** |
| H383 | **no** — mouse has arginine | no | **unusable** |
| H418 | yes | yes | exposure ambiguous between two structures |
| H433 | yes | yes | the one Liu et al. used |

H383 is worth recording as a dead end so nobody rediscovers it: mouse has **arginine**
there, which is positively charged at both pH values and so cannot switch at all. A
pair built against it would fail cross-reactivity outright.

H358 and H370 were missed because both sit outside 415–466, and that block was chosen
as the longest well-conserved stretch rather than by looking for histidines. That was
a real blind spot in the method, not bad luck.

Evaluating H370 then produced a **better epitope than the one we had**: more conserved
reachable anchors, two target histidines instead of one, no overlap with cetuximab's
footprint, and none of the complications that affect 415–466. Details in
`docs/decisions-log.md`.

### 3. We are using a measurement where others will be guessing

EGFR's outer region contains many histidines. Which of them actually produce pH
dependence is not something you can tell by looking at a structure — it was determined
experimentally, one alanine substitution at a time. That result is published, and we
are using it.

That is not novelty. It is reading the literature. But it is a real advantage over
choosing a histidine because it looked exposed.

---

## What is left that is honest

No novelty claim. What the submission argues instead:

1. An epitope chosen by **computed** human/mouse conservation, with the rule against
   contacting position 442 derived from cetuximab's **measured** contact set rather
   than assumed.
2. A target histidine chosen from the **published experimental mapping** of which EGFR
   histidines carry pH dependence.
3. A positional pairing rule with an explicit rejection criterion that has
   **published experimental support**.
4. An affinity target deliberately held **low**, because the requirement is a
   detection threshold rather than a ratio.
5. A record in which every number is regenerated by a script, and in which every
   disproved assumption is kept. **Three of this project's own claims were disproved
   by its own computation**: the cetuximab species-failure explanation, the "only one
   residue is occluded" finding, and this novelty claim.

### On the rules

Using Liu et al.'s mapping of which EGFR histidines matter is using published
information about **the target**. Starting from their **antibody** and modifying it
would be starting from an existing binder, which the competition forbids and which we
are not doing. Our designs are generated from scratch against the target surface. We
state that explicitly rather than leave it to be inferred.

### One deliberate omission

A literature search also surfaced a public repository by another entrant pursuing a
similar approach in this competition. **We did not open it.** Reading a competitor's
approach would mean that any resemblance between our designs and theirs could not
honestly be called independent, and the de novo and zero-shot rules make that awkward
to explain. The only thing reading it would have protected is a modality-novelty claim
we are dropping anyway. The decision and its reasoning are recorded in the decisions
log so it is a choice rather than an oversight.

---

## Why this document is in the submission

A reader who knows this field will think of Liu et al. 2022 within a minute of reading
our epitope choice. There are two ways that can go: they find it themselves and
conclude we did not look, or they find it cited by us, with the part that supports our
design rule drawn out and the part that defeats our novelty claim stated plainly.

The second is better, and it is also just true.

---

## Terms introduced here

In `docs/glossary.md`: alanine substitution, histidine scanning, histidine switching,
recycling antibody, prior art. Plus, from earlier: anchor, epitope, Fab, histidine,
pH, UniProt.
