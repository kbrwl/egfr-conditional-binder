# How to work in this repository

Read this first, then `docs/decisions-log.md`, which records where the project
currently stands.

Abbreviations used here: EGFR is the epidermal growth factor receptor, the protein
we are designing against. UniProt is the public archive of protein sequences. PDB is
the Protein Data Bank, the public archive of measured three-dimensional protein
structures.

---

## What we are trying to do

Design protein binders against the part of human EGFR that sits outside the cell,
for Challenge 1 of the Anthropic x Adaptyv Protein Design Competition. There are
three objectives, listed here in the order of importance the organisers give them:

1. pH selectivity: bind human EGFR at pH 6.5, with no detectable binding at pH 7.4.
2. Mouse cross-reactivity: the same sequence must also bind mouse EGFR.
3. Affinity against human EGFR.

The order matters when the objectives pull against each other, and the organisers
have said so directly: a weak but clearly pH-sensitive binder may outrank a
high-affinity binder that is not pH-sensitive.

Refer to the objectives by name rather than by number. The source we originally
worked from numbered them the other way round, and a sentence in the decisions log
that said "objectives 2 and 3" ended up meaning the opposite of what was intended
once the numbering changed under it.

---

## The affinity rule: aim for weak binding

This is the easiest thing in the project to lose track of, and losing it would spoil
the submission without anything looking wrong. It contradicts what every standard
design tool does by default, so it needs restating at each stage.

The pH-selectivity requirement is a threshold rather than a ratio. "No detectable
binding at pH 7.4" means the pH 7.4 state has to fall below what the measuring
instrument can see at all. It is not a statement about the gap between two numbers.

- A binder at 10 nM at pH 6.5 and 200 nM at pH 7.4 is twentyfold selective and
  fails, because 200 nM is easily detected.
- A binder at 2 µM at pH 6.5 with nothing measurable at pH 7.4 passes, even though
  it is a hundred times weaker and its ratio looks worse.

Standard design pipelines maximise binding strength by default, which produces
exactly the failing case. Do not accept a pipeline's built-in objective. Aim for a
baseline weak enough that the pH 7.4 state disappears under the detection floor,
then build the largest switch achievable on top of that.

If a change makes a design bind more strongly at both pH values, it has probably
made the submission worse. Say so when that happens.

We also do not know where the detection floor sits, because the organisers have not
said which assay they use. The margin we aim for is a judgement rather than a
calculation, and should be described that way.

---

## Who you are writing for

The project owner is a product manager with no biology background, who last studied
it at school, and who is reading every line of this repository in order to
understand it rather than take it on trust.

- Explain what a computation measures and why the answer matters before running it,
  in enough detail that the result can be judged rather than just received.
- Expand every abbreviation the first time it appears in each file. Once per file,
  not once per project, because files get read on their own. Add new terms to
  `docs/glossary.md` rather than explaining them again in conversation.
- When you use a term of art, give the plain meaning in the same sentence. "An
  average deviation of 1.08 angstroms" tells this reader nothing on its own; that
  the two structures sit on top of each other to within about one angstrom, which
  for a protein means the same shape, tells them something.
- Keep the rigour. This is not a request for simplification or approximation. Keep
  every number, caveat and qualification, and say the same thing in plainer words.
- Where a technical decision was made, say what the alternative was and why it was
  rejected.
- A number with no consequence attached is not a finished answer. After each result,
  say what it changes about the design.
- Flag ambiguity rather than resolving it quietly. If a value sits near a
  classification boundary, or two sources disagree, say so instead of rounding it
  into a clean conclusion.
- Say when reasoning is thin. Where something is a guess, a convention or an
  unvalidated assumption, name it. The competition is judged against trained protein
  designers, and an obvious deduction presented as insight is worse than no comment.

### How to write

Plain sentences, ordinary words, the way a knowledgeable colleague explains
something across a desk.

Avoid the contrast construction "X, not Y" unless it prevents a specific expensive
mistake. It is doing real work in "a threshold rather than a ratio", above, because
that distinction stops a particular error. Elsewhere it is only emphasis, and the
repository had accumulated a great deal of it.

Avoid dramatisation. Say what happened and what follows from it. Reserve bold and
capitals for warnings a reader genuinely must not miss; if most paragraphs have
emphasis, none of it registers. Lead with the result rather than building up to it.

---

## Compute rather than recall

Claims in this repository have to be computed.

- Anything stated from memory must be labelled UNVERIFIED, in terminal output and in
  any document where it appears.
- UNVERIFIED items live in their own section of `docs/decisions-log.md` and must not
  be built on. They are resolved by calculation, then moved.
- When a computation disproves something previously assumed, record that it was
  disproved. A disproved assumption is a real result and should not be quietly
  dropped.
- Prefer a script that can be rerun to a number written into prose. If a number
  matters, something in `analysis/` should regenerate it.

This applies to generated output too. The interactive viewer was once written
against a library version number that did not exist, because it was recalled rather
than checked, and no assertion was watching that file.

---

## Anything computed twice must be computed once

Steps 04 and 06 both worked out which residues an antibody touches, from the same
file with the same settings, and reported different answers for months' worth of
conclusions. One of them was wrong.

The rule that came out of it:

- If two scripts need the same quantity, it lives in one function in
  `analysis/egfr_common.py` and both call it.
- If that is impractical, the second script must compare its result against the
  first script's committed output and stop the run when they differ. Use
  `egfr_common.cross_check_residue_set`.
- A check that has never been seen to fail is not known to work. Break the input on
  purpose once and confirm it fails.

Related: pass residue numbers around as `egfr_common.Numbering`, whose methods are
named `uniprot_of` and `pdb_of`, rather than as bare dictionaries. Handing that class
the wrong kind of number gets nothing back. Handing a bare dictionary the wrong kind
of number returned a plausible answer 48 positions away, which is how the
disagreement above happened.

---

## `docs/decisions-log.md` records the current state

Its sections are Settled, Open, Unverified, Resolved, Ruled out and Next actions.

Update it whenever something moves between sections. That is what it is for. A
finding that stays in a findings file and never reaches the log will be lost by the
next session, so moving the item is the part that matters.

---

## Where numbers are allowed to live

`results/` and `data/derived/` hold computed output only. No hand-written numbers,
no hand-corrected values. Everything in them is produced by a script in `analysis/`
and comes back by rerunning it. If a number in those directories looks wrong, fix
the script and regenerate; do not edit the output. The findings files say "do not
hand-edit" and that applies to you as much as anyone.

`docs/` holds prose and decisions, and may quote computed numbers, as long as they
are traceable to the script that produced them.

```
analysis/            numbered standalone scripts, plus egfr_common.py
                     each writes prose findings to results/findings/
                     and tables to data/derived/ as CSV
data/sequences/      inputs (tracked)
data/structures/     downloaded structures (not tracked, re-downloadable)
data/derived/        computed tables (tracked, small)
results/findings/    computed findings in markdown (tracked)
results/candidates/  bulk design output (not tracked, except the final shortlist)
design/              design pipeline configs and runners
submission/          the submission CSV and the methods write-up
explorer/            interactive viewers
```

---

## Residue numbering

Every residue number in this project is a position in the full human UniProt record
P00533, which is 1210 amino acids long. Residues 1 to 24 are the signal peptide, a
short leader that is cut off when the protein is made, so the mature protein outside
the cell runs from 25 to 645.

The official competition constructs are that mature region, so:

    UniProt position = challenge-construct position + 24

Structure files are the greater risk. They carry their own numbering and nothing in
the file records which convention it follows; files of secreted proteins often count
from the mature protein. Both structures used here are offset by 24. Measure the
offset by aligning the structure's own sequence against the reference, check the
result against residues whose identity you already know, and print the check. Never
assume it, in either direction, because getting it wrong shifts every result and
raises no error.

`analysis/00_numbering_check.py` asserts the convention and exits with an error if it
breaks. Run it first if anything looks strange.

---

## The design rule

The pH switch works by pairing charges position by position. Histidine is the only
amino acid that changes charge between pH 7.4 and pH 6.5: neutral above, positive
below, because the pH at which it is half charged sits around 6.0 to 6.5.

| On the target | Put on the binder | At pH 7.4 | At pH 6.5 |
|---|---|---|---|
| D or E, acidic and always negative | histidine | neutral, no attraction | histidine turns positive and is attracted |
| H, which switches | D or E | neutral histidine, no attraction | histidine turns positive and is attracted |

Both arrangements switch on in the same direction as pH falls, so they reinforce
each other.

Reject any candidate with a histidine on the binder facing target H418 or H433. That
pair switches off as pH drops, because both become positive and repel, and it can
cancel out a correctly built pair elsewhere. Reject such candidates rather than
scoring them.

Also reject any candidate that contacts position 442. It is the only difference
between human and mouse inside our epitope, and it is in cetuximab's contact set, so
touching it risks species-specific behaviour at the one position where the species
differ.

Two limits to respect. The switch is partial rather than complete at pH 6.5, so
stack three or four pairs instead of relying on one. And the pH at which a histidine
flips shifts depending on its neighbours, which is much of why pH selectivity
resists reliable prediction. Do not present a predicted switch as a measured one.

---

## Submission hygiene

Embedded instructions or prompt injection in a submission is grounds for
disqualification, and selection is performed partly by a Claude workflow that reads
our documentation. Our documentation therefore states facts and reasoning and does
not address its reader as an instruction. Nothing in `submission/` should try to
influence a reader's judgement by any means other than argument.

---

## Commits

Commit messages describe findings. "Verified numbering convention holds; 415-466
anchors conserved in both species" rather than "add scripts".
