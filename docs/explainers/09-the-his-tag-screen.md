# Explainer 09 — the His-tag screen, and keeping the binder's tail out of the way

Why a binder could look like a pH switch while really gripping a purification tag,
what we found BindCraft2 offers against that, and what the two campaign files now
do. The script is `analysis/15_histag_counterscreen.py`. The campaign files are
written by `analysis/09_trim_target.py`.

Written for a reader with no biology or software background. It stands on its own:
every term is explained where it first appears.

---

## How this document differs from the others

**This covers a finished script that has processed no real candidate.** Nothing has
been designed yet, so every result below comes from inputs built to test the script.
The configuration that avoids the tag has also never been run on real hardware. The
document is rewritten when there is real output.

It is written now because the screen exists to be understood before its output is
believed, and because the methods write-up will want to say the check was made.

---

## 1. The result

A binder can look pH-selective while binding something unrelated to our target, and
this project's design rule makes that more likely than usual. A screen for it now
exists, and the design run can be told to avoid the problem as it designs, not only
check for it afterwards.

Three things came out of reading BindCraft2's documentation and source on 2 October.

1. **BindCraft2 can steer a binder away from a chosen sequence while it designs.** A
   target given a negative weight is avoided. The His tag can be given as a short
   sequence, which BindCraft2 co-folds with the binder as a floppy disordered region.
2. **BindCraft2 has a setting that keeps the binder's tail out of the binding
   surface,** which the assay construct needs. It is switched on in both campaign
   files.
3. **BindCraft2 can design one binder against several targets at once,** with an
   example of one binder that engages both the human and the mouse form of a protein.
   That is relevant to the mouse objective and is listed as an option, with its cost,
   at the end.

What this changes about the design: the campaign should be run from the file that
includes the tag as an off-target, if that file starts. A design generated without it
cannot be screened for tag binding afterwards without a prediction run that does not
yet exist.

---

## 2. The problem

**EGFR** (epidermal growth factor receptor) is the protein on the surface of cells
that this competition targets. We are designing a **binder**, a small protein that
sticks to a chosen patch of it. The competition requires the binder to stick at
pH 6.5, the mildly acidic condition around a tumour, and not at pH 7.4, the condition
in healthy tissue.

The switch rests on **histidine**, the one amino acid whose electrical charge changes
across that range: neutral at pH 7.4, and a useful fraction positively charged at
pH 6.5. Our design places acidic residues on the binder, which carry a negative
charge, directly across from the target's own histidines, H358 and H370. At pH 6.5
the histidine gains a charge and the two attract; at pH 7.4 they do not.

Both the human and the mouse target in the assay carry a **His tag** at their
C-terminal end. A His tag is a short run of histidine residues, usually six, added so
the protein can be caught on a metal column during purification. The organisers said
on 30 September that removing the tags is awkward and the screen will probably use the
protein as it is (`docs/competition-qa-log.md`).

A pocket on our binder shaped to grip a protonated histidine will grip the histidines
of a tag as readily. A binder that does so would also change its grip with pH, since
the tag's histidines gain a charge in the same way. It would look like exactly the
switch we want, and it would stick to anything carrying a His tag. Amir at Anthropic
said such a binder should not be selected, and Tudor at Adaptyv said the organisers
will give more weight to the computational scoring where a design might be hitting a
tag.

The project is more exposed than most to this, for a structural reason. The target's
own histidines are held in place by the fold of the protein, while the tag is a free
tail that can thread into a pocket, so a pocket built for H358 and H370 may fit the
tag better than it fits the real thing.

---

## 3. How the screen was designed, and what was rejected

The work order asked for the predictor and the meaning of "accepts the tag" to be
decided before the script was written. They were decided from BindCraft2's own
documentation and source.

### The predictor

**BindCraft2's own AlphaFold2, used through its detargeting objective.** AlphaFold2 is
the structure predictor that validates every design inside BindCraft2's loop. A
**detargeted** target is one the binder is pushed away from. Its output tables then
record the binder's predicted interface confidence with that target, as
`i_pTM_detarget`, together with `i_pAE_detarget` and `Interface_Residues_detarget`.
`i_pTM` is the predictor's confidence in the interface of a predicted complex, from 0
to 1.

Three alternatives were rejected.

- **`bindcraft score`.** This was the obvious candidate, since BindCraft2 has a
  command with that name. Its own help says it scores one existing design structure
  on checks taken off coordinates, and that the prediction readings are reported from
  the design's own record and are not recomputed. It cannot predict a binder against a
  peptide.
- **Placing a peptide by geometry.** A six-residue floppy tail has no fixed pose to
  place.
- **A separate prediction run after the campaign.** It would need a prediction runner
  that does not exist yet. The in-loop version adds no run and steers designs away
  from the tag, where the separate run would only find them afterwards.

### "Accepts it"

A candidate accepts the tag if `i_pTM_detarget` is at least 0.4.

That figure is BindCraft2's own default ceiling for rejecting an off-target
(`max_detarget_iptm`). **It is the tool's default and not a threshold calibrated for
this tag.** Calibrating one needs binders known to grip a His tag and binders known
not to, and neither exists. The threshold rests on the tool's own choice and on
nothing measured here, and the findings file says so.

The screen asks whether the binder is predicted to bind the tag at all, and does not
ask whether it does so through its acidic pocket. A binder that grips the tag by any
residues confounds the assay in the same way.

### A blank is not a pass

A candidate with no off-target reading is reported "not recorded", which means no
evidence either way. BindCraft2's documentation says not to treat a blank as zero.
This is the situation for every design generated from a campaign file that does not
include the tag.

### A pass is weak evidence

A prediction that finds no interface does not show there is none. AlphaFold2 can
miss a real one, especially for a floppy tail. The screen can show that a risk is
present and cannot show that it is absent, and the methods write-up should say so.

---

## 4. What the campaign files now do

`analysis/09_trim_target.py` writes two campaign files from one set of values, so
they cannot disagree about which positions are aimed at.

- `design/configs/egfr-domain3-h370.json` is the plain campaign: the human fragment,
  the eight hotspots, and the termini setting below.
- `design/configs/egfr-domain3-h370-notag.json` is the same campaign with the His tag
  added as an off-target.

The off-target is the sequence GSHHHHHHGS at weight −0.5. The six histidines are the
tag. The serine and glycine on each side are there because BindCraft2 samples a window
of at least 10 residues from a sequence target by default, and six is shorter than
that, and because a tag sits at the end of a chain with something before it. The
weight is the value in BindCraft2's own example of a cross-reactive campaign, where a
magnitude of 0.5 means the push away is half the strength of the pull towards a
binding target. It is a starting value and has not been tuned.

### The termini setting

Our designs are held on the sensor chip by a tail added to the binder's C-terminal
end: a linker, then GFP11, then a linker, then a twin-Strep tag, with the binder
attached to the chip through that end. So the C-terminal end of the binder must not
be part of the binding surface, or the tail and the chip would be in the way.

BindCraft2's property `termini_accessible` points both chain ends of the binder away
from the target. It adds a term to the design objective steering the ends away, and a
final filter requiring the direction of each end, measured as a cosine where +1 is
away from the target and −1 is toward it, to be at least 0.0. The 0.0 floor, a right
angle, is the tool's default. The work order's fallback was a filter in step 10, which
is not needed. It is geometry only and says nothing about whether the tail expresses
or folds.

---

## 5. How the screen is tested

No candidate exists, so the screen is tested on constructed inputs, run every time.

**The rule, on constructed readings:** 0.12, 0.39, 0.40, 0.71, an empty cell, a
missing column, and a value that is not a number. The first two must come out clear,
0.40 and 0.71 must come out as accepting the tag, and the last three must come out as
not recorded. The reading of exactly 0.40 is included on purpose, since the threshold
is "at least".

**A constructed campaign folder,** in BindCraft2's documented layout with its real
column names, is read through the same reader step 10 uses. This exercises the path
real output will take.

**The campaign files themselves** are checked: the tag is present as a sequence
target with a negative weight, its sequence contains six consecutive histidines, the
sequence file exists, the termini setting is on in both files, and the two files give
the EGFR target and hotspots identically.

**Each rule is switched off in turn.** With the threshold removed, the readings of 0.40
and 0.71 come out clear and two tests fail. With a blank treated as zero, the three
blank-like cases come out clear and three tests fail. In both cases the tests are
known to be doing work.

---

## 6. What this does not settle

- **The detargeting variant has not been run.** Whether BindCraft2 accepts a short
  sequence as an off-target alongside a structured target, and what it does to the run
  time and the hit rate, are untested. A start-up failure would cost little, and the
  plain campaign file is available.
- **The weight of −0.5 and the ceiling of 0.4 are starting values.** Neither is
  calibrated for this tag.
- **No candidate has been screened.** The evidence is the rule's tests and a check of
  the files.
- **A tag binder might not look like a tag binder to AlphaFold2.** See section 3.
- **Mouse cross-reactivity is not addressed by any of this.** BindCraft2's
  multi-target mode, one binder against human and mouse together with hotspots given in
  each structure's own numbering, is the nearest tool for it. It needs a structure of
  mouse domain III, which we do not hold, and the human and mouse forms differ at N361
  in a way that affects the sugar chains (explainer 08). Using it is a decision for the
  project owner and is recorded as an open question in the decisions log.

---

## Terms introduced here

For `docs/glossary.md`: His tag, detargeting, off-target, i_pTM. His tag is already
defined from explainer 07. Detargeting, off-target and i_pTM are added in the same
commit as this document.

Terms used here and defined in earlier explainers: amino acid, binder, EGFR,
histidine, protein, residue, tethered and extended conformations.
