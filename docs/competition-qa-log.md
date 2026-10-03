# Organiser answers from the Proteinbase Slack

**What this file is.** A standing record of what the competition organisers have
actually said, in their own Slack workspace, in answer to participants' questions.
It exists because the competition page and the FAQ do not carry most of it, and
because four more challenges follow this one, so the same questions will recur.

**What this file is not.** It is not a findings file. Nothing here was computed by
this project. Every line is external information, read on the date given, from the
person named. Where it disagrees with the official competition page, the page wins
and this file gets corrected. Where it bears on a decision, the decision lives in
`docs/decisions-log.md` and the reasoning lives in an explainer, not here.

Read on 2 October 2026 from the `proteinbase.slack.com` workspace, channels
`#anthropic_adaptyv_competition`, `#general`, `#design-methods` and `#feedback`.

Names and roles, so attributions below are readable:

- **Tudor-Stefan Cotet** (`tudor@adaptyvbio.com`) — Adaptyv, running the competition.
  Source of nearly all the assay detail.
- **Simon Dürr** (`simond@adaptyvbio.com`) — Adaptyv, owns the novelty scoring.
- **Amir Shanehsazzadeh** (`ashanehsazzadeh@anthropic.com`) — Anthropic.
- **Bailey Bova** (`bailey@anthropic.com`) — Anthropic, community and access.

---

## 1. How binding will actually be measured

**Method.** Surface plasmon resonance, abbreviated SPR, as the main assay.
Bio-layer interferometry, abbreviated BLI, only to cross-check hits.
(Tudor, 30 September.)

Both are instruments that measure two proteins sticking together without
attaching a dye or a radioactive label to either one. One partner is fixed to a
sensor surface, the other is washed over it in solution, and the instrument
watches the mass accumulating on the surface in real time.

**Which partner is fixed.** **The designs are immobilised and the target is
flowed over them.** The target is caught on the surface through a C-terminal
twin-Strep tag, described as the default and changeable. (Tudor, 30 September.)

This is the reverse of the arrangement most people assume. The protein in
solution, called the **analyte**, is EGFR. The protein stuck to the chip, called
the **ligand**, is our binder.

**Concentrations.** Top analyte concentration 1000 nM, with a usable affinity
window of roughly 0.1 nM to 10 µM. (Tudor, 30 September.)

Affinity is reported as **KD**, the dissociation constant, in units of
concentration. Lower means tighter. A KD of 1 nM is a strong binder; 10 µM is
barely binding at all. The window above is the range in which the instrument can
return a number.

## 2. What counts as binding, and as not binding

A design counts as **binding** at a given pH if either of these holds
(Tudor, 30 September):

1. the sensor trace can be fitted to give a KD, or
2. where no fit is possible, the association signal rises more than 300% above
   the negative control.

If both fail, the organisers declare **no detectable binding**.

"Association" is the rising part of the trace, while the target is flowing on.
So the second test is a crude height check for designs whose curves are too ugly
to fit.

## 3. How the pH objective is ranked

Asked directly whether no detectable binding at 7.4 is a hard requirement, or
whether a large KD shift with binding at both pHs also qualifies, Tudor answered
on 30 September:

> Designs that show a large KD shift will also qualify, but designs that show no
> binding at pH 7.4 and a high affinity at 6.5 will be ranked higher.

Two separate things follow. Binding at both pHs is not disqualifying if the
shift is large. And the top of the ranking wants **high affinity at pH 6.5**
together with nothing at 7.4.

Measurement pH values are 6.5 and 7.4. Tudor noted on 28 September that this is
deliberately harder than measuring at 6.0 or 5.0 would be, and pointed at
PubMed 36458200 as the source of the conditions.

## 4. What counts as mouse cross-reactive

A KD ratio of mouse over human of **around 1**. (Tudor, 30 September.) He
pointed at the same paper, Mol Ther Oncolytics
`S2372-7705(22)00136-X`, as the comparison.

This is stricter than "binds both". A design that binds human at 50 nM and mouse
at 2 µM binds both and is not cross-reactive by this definition.

## 5. Buffer and ionic strength

Base buffer, used for the pH 7.4 condition (Tudor, 29 September):

| component | concentration | what it is |
|---|---|---|
| HEPES | 10 mM | the thing holding pH steady at 7.4 |
| NaCl | 150 mM | ordinary salt, sets the ionic strength |
| Tween-20 | 0.2% | detergent, stops proteins sticking to plastic and to each other |
| EDTA | 3 mM | mops up stray metal ions |

For the acidic condition HEPES is replaced with **MES** (Tudor, 30 September).
MES is a buffer whose useful range covers pH 6.5, where HEPES has none.

**Total ionic strength is matched between the two conditions with NaCl, at
approximately 170 mM.** (Tudor, 29 September.)

That matching is deliberate and it matters to us. Ionic strength is how many
charged particles are floating in the solution. A high one surrounds every
charged group on a protein with a cloud of counter-ions that cancels much of its
pull, which is called **screening**. Because the two conditions are matched, any
difference we see between pH 6.5 and 7.4 is from the pH and not from the salt.
Because the level is physiological rather than low, a charge pair sitting out in
the open on the surface contributes less than the same pair would in a
low-salt buffer.

Tudor added that these are the organisers' standard platform conditions,
described in `nature.com/articles/s41587-026-03187-0`.

## 6. The target as it will actually exist

| property | value | source |
|---|---|---|
| conformation | **tethered** (the closed, inactive form) | Tudor, 29 September |
| expression system | HEK293 cells | Tudor, 30 September |
| glycosylation | **glycosylated** | Tudor, 30 September |
| tags, both species | C-terminal His tag | Tudor, 30 September |
| capture tag | C-terminal twin-Strep | Tudor, 30 September |

On 30 September Tudor restated the whole of it in one line: the wet-lab screen is
against the glycosylated, tethered EGFR, and that is why domain III was suggested
as the place to aim.

**On the tethered form.** EGFR's outer region has two well-known shapes. In the
**tethered** or closed form, domain II folds back and clips onto domain IV, and
the receptor is inactive. In the **extended** or open form it unfolds and the
receptor can pair up with another copy and switch the cell on. Different surfaces
are exposed in each. Asked which one to target, Tudor said the tethered form, and
said the competition page would be updated to say so.

**On glycosylation.** A **glycan** is a tree of sugars attached to the protein
after it is made. **N-linked** glycans attach to an asparagine that sits in a
three-residue pattern called a **sequon**, written N-X-S/T, meaning asparagine,
then almost any residue, then serine or threonine. Expression in HEK293, a human
cell line, means the sugars are there and are human-like. They are large, they
move about, and a crystal structure shows only the first ordered sugar or two, so
their real reach is larger than any structure suggests.

## 7. The His tag, and why it is a trap for a pH-sensitive binder

This came up as a joke and turned serious. Jonathan Ouyang asked on 30 September
what stops someone simply designing a binder against the His tag, and noted that
histidine is exactly the residue that responds to pH.

Tudor's answers, same day:

- the tags are unlikely to be cleaved, and the screening will probably use
  His-tagged EGFR as-is
- follow-up neutralisation assays will be run, which is often easier than
  cleaving tags
- in silico scoring will be weighted more heavily in cases where a design might
  be hitting a tag, because the structure prediction has no tag in it

Amir Shanehsazzadeh added:

> The tag or really any component of the system that is there for screening
> should definitely not be targeted. Good call out given it's Histidine. Any such
> binder might look pH-selective but it would bind to anything with a His tag.

A **His tag** is a short run of histidine residues, usually six, added to a
protein so it can be purified by sticking it to a metal column. The problem for
this project specifically: our pH mechanism puts acidic residues on the binder
opposite the target's own histidines. A pocket shaped to grip a histidine will
grip the histidines in a tag just as readily, and the tag is floppy and exposed
where the real target histidines are held in a fold.

Chengbo Chen raised the same shape of problem for glycans, and Jonathan Ouyang
suggested treating glycan sites as coldspots, that is, positions a design is told
to avoid.

## 8. What our binder gets attached to

Asked what tags and linkers go on the designs, Tudor answered on 30 September:

> We use C-terminal tags by default and add the following to the constructs:
> linker - GFP11 - linker - TwinStrep tag. Immobilization is on the C-terminus.

Also: **cyclic peptides are not supported in this competition, linear only.**

**GFP11** is the eleventh strand of green fluorescent protein, about sixteen
residues. On its own it does nothing; when the other ten strands are supplied
separately it completes them and the pair glows, which is how the amount of
protein on the surface gets measured. It is a beta strand, a flat extended piece
of chain that likes to pair up edge-to-edge with other beta strands.

Two consequences for a design. The C-terminal end of our binder is the end that
gets tied down, so it should not be part of the binding surface. And a binder
whose own edge is an exposed beta strand has a fused beta strand dangling beside
it.

## 9. Novelty: the gate and how it is checked

**The gate is level 3 of 4.** Simon Dürr said so twice, on 29 September
("Level 3 will pass the submission gate") and on 1 October, where he described it
as the sequence needing to be under 30% similar to anything existing. The scale
is defined at `adaptyvbio.com/blog/novelty`.

*The two statements are not obviously the same thing and this should be checked
against the blog before relying on either. Recorded here as stated, unresolved.*

**How they check it** (Tudor, 1 October): in previous competitions they ran
MMseqs2, a fast sequence-similarity search tool, against SwissProt, the Protein
Data Bank, the USPTO patent database, the EBI patent database, THPdb, PLAbDab and
Proteinbase itself, and called a design de novo if it hit nothing with any
homology in any of them.

Novelty is scored **automatically on upload**, so a design can be tested against
the real checker before the deadline rather than guessed at.

Simon noted that the pipeline currently gives antibodies and nanobodies special
treatment that other rigid scaffolds such as affibodies and DARPins do not get,
and that they will try to fix this for the next challenge.

## 10. How designs get selected, which is the whole game for Track 3

Asked whether a strategy relying on pH-dependent conformational change would be
penalised by the selection modelling, Tudor answered on 2 October:

> Selection will be mainly based on method novelty, design diversity, and a
> couple of in silico/confidence metrics. [...] Compared to previous competitions
> (e.g. where we only selected by ipSAE or ipAE), this selection strategy will
> likely be less biased against more creative design methods. We'll share the
> selection prompts after the competition.

And on 30 September, to Yannick Kiefl:

> In silico scoring will only be a smaller component in the overall selection
> criteria. [...] worse confidence scores on the rigid structure or on obstructed
> epitopes won't be penalized much (if at all).

Taken together: the written method and the spread of the designs carry more
weight than the confidence numbers. A set of twenty near-identical binders spends
twenty slots on one idea.

## 11. Submission mechanics

From the submission page, quoted in `#general` on 1 October by a participant:

- up to 20 proteins per submission
- each chain between 10 and 250 amino acids
- **one submission every 24 hours**
- novelty score 3 of 4 or higher, checked automatically after upload

The 24-hour cadence means an early safe batch can be replaced by a better one
later, rather than the whole entry resting on one upload.

**Track 3 needs no confirmation email.** Tudor, 30 September: Track 3
pre-registration existed only to anticipate how much compute the organisers would
need for selection, and designs can be submitted directly.

**Deadline.** Amir confirmed on 29 September that it was Sunday 4 October,
anywhere on Earth. Tudor's 30 September message announced an extension for the Modal
credit delays, read at the time as one day (5 October). **Superseded: the official
announcement gives 6 October 2026, 23:59 anywhere on Earth** (section 13a). The second
challenge still starts on time, so challenges 1 and 2 overlap by about two days.

## 12. Public data the organisers have released

Mentioned in `#general`, all public and none of it another current entrant's
working material:

- the EGFR round-1 competition results collection on Proteinbase
- a negative data bundle, which is the set of designs that did not work — rarer
  and often more useful than the successes
- the Nipah binder competition results, 171 additional proteins in the last batch

Anthony Gitter reported on 18 April that his test of the *unselected* EGFR
round-1 designs came back negative across the board except for the EGF positive
control.

In `#design-methods`, Nick Boyd published a write-up on why his team won the
de novo category of the Nipah competition, concluding that it came down mainly to
which epitope they targeted: `blog.escalante.bio`.

## 13. Practical notes

**Model access.** Amir announced on 28 September that Claude Sonnet 5.5 performs
roughly on par with Opus 5.5 at de novo binder design, and carries lighter bio
classifiers, making it more usable for this work. Several participants in
`#feedback` report Opus 5.5 refusing ordinary protein-design prompts.

**Accelerated models.** Amir, 29 September: Anthropic has published accelerated
versions of popular open-source design and folding models, averaging 4–5x
speed-ups. `github.com/anthropics/uplifting-biomolecular-modeling`.

## 13a. Deadline extended to 6 October (recorded 3 October 2026)

Official announcement in the Proteinbase Slack, relayed to this project on 3 October:
challenge 1 is extended to **6 October 2026, 23:59 Anywhere on Earth**, because some
teams received their Modal credits late. The competition page is to be updated to match.
Challenge 2 starts on time, so the challenges overlap by about two days. This supersedes
the earlier reading of Tudor's 30 September message as a one-day extension to 5 October.
Message author and exact wording were not captured at the time of recording.

## 14. Asked and still unanswered as of 2 October

- The exact mouse EGFR construct, residue range and vendor. Olga Lavinda asked
  twice, on 1 October, and got no reply. Our own construct check on 30 September
  was against the sequence on the competition page, so this is a confirmation we
  do not have.
- Whether extra metric columns are tolerated in the submission CSV beyond the
  three in the template. Asked by Jose Farias on 1 October, unanswered.
- **Whether cynomolgus monkey cross-reactivity is in scope for challenge 1.**
  Amir's pre-launch message on 28 September said the challenge wanted
  cross-reactive binding "to mouse and cyno". The official objectives as recorded
  in `docs/rules-reference.md` name human and mouse only. One of the two is out
  of date and it should be checked on the competition page.
