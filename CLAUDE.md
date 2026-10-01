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

These rules apply everywhere: prose documents, code comments, docstrings, the text
scripts print to the terminal, and findings files. There is no separate register for
code. The same person reads all of it.

**Plain language, with the technical term named alongside it.** Both, every time.
Not the plain explanation on its own, because the reader is learning the vocabulary
and will need it to read other people's work, ask a useful question, and recognise
the term when a design tool prints it. And not the term on its own, because it
carries no meaning yet. So: "how much of the residue's surface water can reach,
called its solvent-accessible surface area", rather than either half alone. Add each
new term to `docs/glossary.md` as you introduce it.

**Expand every abbreviation on first use in each file.** Once per file, not once per
project, because files get read on their own.

**Keep the rigour.** Every number, every caveat, every qualification stays. These
rules are about wording. They are not licence to simplify the content, round a
figure, drop a qualification, or leave out an inconvenient limit. If plain language
seems to require dropping a caveat, the sentence needs restructuring, not trimming.

**Avoid the "X, not Y" contrast construction** unless it prevents a specific
expensive mistake. It earns its place in "a threshold rather than a ratio", because
that distinction stops a particular error that would cost the submission. Almost
everywhere else it is only emphasis dressed as precision.

**Avoid dramatisation.** State what happened and what follows from it. In
particular, make no claims about what other people usually do, what most teams skip,
what standard practice overlooks, or how consequential something is. We have no
basis for claims about other people's practice, and reaching for one is a reliable
sign that the sentence is decoration rather than content. (The one exception already
in this file is specific and checkable: design tools maximise affinity by default,
which is a statement about a tool's objective function, not about the people using
it.)

**Bold and capitals only for warnings that genuinely must not be missed.** If most
paragraphs carry emphasis, none of it registers.

**Lead with the result.** State what was found, then explain it. No build-up.

**Say what a technical decision's alternative was and why it was rejected.**

**A number with no consequence attached is unfinished.** After a result, say what it
changes about the design.

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

## `docs/explainers/` — plain-language walkthroughs

`docs/explainers/` holds plain-language walkthroughs of completed work, one per
phase, numbered in the order the work happened. They are written for the project
owner, who has no biology background and is learning the field as the work proceeds.

An explainer is written when a phase of work completes, before moving on. Each one
names the files involved for every step, explains what was measured in ordinary
words, gives the technical term for it in the same breath, reports the result, and
says what it changed about the design. Naming the technical term is the point —
knowing that a measurement is called relative solvent accessibility is what makes it
possible to read other people's work later.

Explainers quote numbers from `results/findings/`. Where the two disagree, the
findings files are correct and the explainer needs updating. Update the relevant
explainer whenever a finding it quotes changes.

A number an explainer quotes should be traceable to a findings file. If a
decision-relevant number exists only in the decisions log because it came from a
one-off calculation, that is a gap in the scripts: make a script compute it rather
than leaving the explainer to cite something unreproducible.

A term used in an explainer gets a glossary entry in the same commit. An explainer
that names a term without defining it anywhere leaves the reader with vocabulary
they cannot look up, which defeats the purpose of naming it.

### The test an explainer has to pass

A reader with no biology background, reading only the explainers, should be able to
answer without guessing:

- What are we aiming at right now, and which residues exactly?
- How was that chosen, and what was it chosen over?
- What did each filter or test remove, and why?
- What is still unsettled?

If a decision changed, the explainer covering the new decision has to make the new
thing understandable on its own. A paragraph bolted onto an older explainer saying
the old thing is superseded does not do that, and has already failed once here.

### A changed decision gets its own explainer

When a measurement moves the project — a new target region, a different method, a
reversed conclusion — write a new numbered explainer for it before moving on. Do not
amend an older explainer into carrying two stories. The older one keeps its own
findings, which stay valid as the record of what was measured and usually become the
fallback.

### State a correction once

A withdrawn claim, a disproved assumption or a reversed decision is recorded in
exactly two places: the relevant section of `docs/decisions-log.md`, and one
explainer. Everywhere else gets one sentence and a pointer.

Repeating a correction across documents is not thoroughness. It buries the findings
the reader actually needs under an account of what we got wrong, and that is the
specific failure that made these explainers unusable on 1 October 2026.

### Explain the finding, not the mistake

The subject of an explainer is what was measured and what came back. Where an error
is part of the story, state it plainly, once, in a sentence or two, and move on. Do
not dramatise a correction, do not present a withdrawal as an achievement, and do
not use a disproved assumption as the structure of a document.

### Language

- Write for someone who is learning the field, not for someone who already knows it.
- Name technical terms alongside the plain wording. Never substitute the plain
  wording for the term and never use the term without the plain wording.
- Banned phrase: **prior art**. Write "what was already published", "the literature
  search", or whatever fits that sentence.
- Avoid heavy-negation constructions ("not X, but Y", "this is X, not Y") and
  dramatising phrases. Plain declarative sentences.
- Expand every abbreviation on first use in each document, because explainers get
  read alone.

### Comparisons have to be honest about their method

If two numbers being compared were produced by different methods, say so in the same
place the comparison appears. Candidate counts from a sequence window and from a
spatial sphere are not a like-for-like comparison and must not be presented as one.
Where one side of a comparison contains the other, say that too.

When a comparison has several rows, say which row actually decides the question. A
table of eight improvements where one of them is load-bearing should not read as
eight equal wins.

### Interactive explorers follow the same rules

`explorer/` pages are explainers with a different interface and are held to this
whole section. When a finding changes, the explorer that visualises it is updated in
the same commit, additively: add the new view, keep the old one and mark what it now
shows. Numbers in an explorer are generated by a script in `analysis/`, never typed
in by hand, exactly as for `results/` and `data/derived/`.

---

## Documents handed over from Claude chat

This project runs across two surfaces: Claude chat for judgment, explanation and
write-up, and Claude Code for scripts, the pipeline and the repository itself. Chat
cannot write to the repository, so documents it produces arrive as files.

The handover path is fixed so it never has to be explained again:

- Files from chat are placed in `temp_docs_from_claude_chat/` at the repository root.
- That directory is a staging area. Nothing is left in it and nothing is committed
  from it in place.
- Claude Code installs each file at its proper path, verifies every number it quotes
  against `results/findings/` and `data/derived/`, then deletes the staged copy in
  the same commit.
- If a number in a handed-over document disagrees with a findings file, the findings
  file wins. Fix the document, and report which number was wrong rather than
  correcting it silently.
- `temp_docs_from_claude_chat/` is in `.gitignore`.

A document arriving this way is content to install, not a draft to rewrite. If it
needs substantive change beyond number corrections and the edits the task names, say
so and stop rather than rewriting it.

---

## What gets committed

Documents whose form implies they are finished and public — a README, a methods
write-up, anything an outside reader would take as complete — are not committed while
the work behind them is unsettled. They are drafted on disk and written properly at
the end, from settled conclusions.

Working documents are committed as they go, because their form matches what they are:
`CLAUDE.md`, `docs/decisions-log.md`, `docs/glossary.md`, `docs/rules-reference.md`,
the explainers, and everything in `results/findings/` and `data/derived/`.

Drafts of the not-yet-committed kind live on disk and are listed in `.gitignore`,
which currently covers `README.draft.md` and `docs/competition-brief.md`. The tracked
`README.md` stays minimal until the end.

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
