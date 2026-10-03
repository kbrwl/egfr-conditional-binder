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

1. **pH selectivity**: bind human EGFR at pH 6.5, with no detectable binding at
   pH 7.4.
2. **Mouse cross-reactivity**: the same sequence must also bind mouse EGFR.
3. **Human binding**: bind human EGFR at all, with affinity as the measure.

The order matters when the objectives pull against each other, and the organisers
have said so directly: a weak but clearly pH-sensitive binder may outrank a
high-affinity binder that is not pH-sensitive.

**Refer to the objectives by name rather than by number.** The source we originally
worked from numbered them the other way round, and a sentence in the decisions log
that said "objectives 2 and 3" ended up meaning the opposite of what was intended
once the numbering changed under it. The decisions log and this file also order them
differently — the log lists them as a design task, this file in the organisers' order
of importance — so a number alone is ambiguous even now.

**The three names are "pH selectivity", "mouse cross-reactivity" and "human
binding".** Use those exact words everywhere: in prose, in code comments, in
docstrings, in terminal output and in findings files. "Affinity" names the
measurement, not the objective, so "fails human binding" rather than "fails
affinity". `docs/rules-reference.md` quotes the organisers' own list, which still
says "Affinity against human EGFR" there because that is a record of their wording,
not ours.

---

## Read this first

You have read this file. **You have not read `docs/`.** Nothing in `docs/`,
`results/findings/` or `data/derived/` is in your context unless you open it.

**Always, at the start of every session:** open `docs/decisions-log.md`. It records
where the project stands, in six sections — Settled, Open, Unverified, Resolved,
Ruled out and Next actions — and it is authoritative for all of them. If a decision
is not written there, it has not been made. Do not infer one.

**Then, before the specific kinds of work below, open the file named:**

| If you are about to… | Open first |
|---|---|
| Say anything about what the competition requires — the assay, the deadline, tags, buffers, novelty, submission format, how designs are selected | `docs/rules-reference.md`, then `docs/competition-qa-log.md` |
| Set, change or defend an affinity target, or decide how candidates are ranked | "The affinity rule" below, then `docs/explainers/07-the-assay-and-what-it-changes.md` |
| Write or change any script in `analysis/` | `analysis/egfr_common.py`, then "Anything computed twice must be computed once" below |
| Quote any number in prose, a comment or terminal output | the file in `results/findings/` that produced it |
| Write or change an explainer | "`docs/explainers/`" below, then `docs/glossary.md` |
| Touch the trim, the fragment handed to the design run, or the campaign configuration | `results/findings/09-trimmed-target.md` and `11-trim-boundary.md`, then `docs/explainers/06-checking-the-trimmed-target.md` |
| Touch the charge-pair filter, the rejection rules or the shortlist ordering | "The design rule" below, then the header comments in `analysis/10_charge_pair_filter.py` and `docs/explainers/05-input-and-filter.md` |
| Change the anchor set, or weigh a risk attached to one anchor | "The design rule" below, then `results/findings/14-glycan-sequons.md` and `docs/explainers/08-sugar-chains-near-the-anchors.md` |
| Change how the His tag is avoided, in the design loop or in the audit | the tag constants in `analysis/egfr_common.py`, then `analysis/15_histag_counterscreen.py` |
| Spend anything on Modal, or start any campaign against EGFR | Next actions in `docs/decisions-log.md`, then "What this does not settle" in `docs/explainers/07-the-assay-and-what-it-changes.md` |
| Install a file from `temp_docs_from_claude_chat/` | "Documents handed over from Claude chat" below |
| State what the organisers have said about anything | `docs/competition-qa-log.md`. Never from memory, never from this file |
| Write or change anything in `submission/` | "Submission hygiene" below |
| Wonder whether something has already been settled | `docs/decisions-log.md`, then ask Kunal |

**These are not suggestions.** Working from this file alone will produce output that
contradicts decisions you cannot see, and has already done so once. If a task
touches a row above and you have not opened that file in this session, open it
before writing anything.

---

## Where things live

| File or directory | Contains | Status |
|---|---|---|
| `docs/decisions-log.md` | Where the project stands. Settled, Open, Unverified, Resolved, Ruled out, Next actions | **Authoritative** |
| `docs/rules-reference.md` | Official competition facts, from the competition page and FAQ | **Binding** |
| `docs/competition-qa-log.md` | What the organisers have said in the Proteinbase Slack, with speaker and date. Grows across all five challenges | **Binding**, superseded by the competition page where they disagree |
| `docs/glossary.md` | Every term any explainer names. One entry per term, added in the same commit as the explainer that introduces it | |
| `docs/explainers/NN-*.md` | Plain-language walkthroughs, one per phase, numbered in the order the work happened | |
| `results/findings/NN-*.md` | Computed output of the matching `analysis/` script | **Never hand-edited** |
| `data/derived/*.csv` | Computed tables | **Never hand-edited** |
| `analysis/egfr_common.py` | Anything two scripts both need: the numbering class, the anchor set, the contact calculation | |
| `analysis/NN_*.py` | The analysis steps, numbered in the order they run | |
| `design/configs/` | Campaign files for the design run, generated by `analysis/09` | Generated |
| `design/modal/` | The Modal runner and its instrumentation | |
| `submission/` | The submission CSV and the methods write-up | |
| `explorer/` | Interactive explainers. Held to the explainer rules below | |
| `temp_docs_from_claude_chat/` | Staging only. Nothing stays here | Gitignored |

Where two of these disagree: a findings file beats a document, the competition page
beats `docs/competition-qa-log.md`, and `docs/decisions-log.md` beats everything
else about what we have decided.

---

## These come back to Kunal, in chat

Do not decide them in the repository. Flag and stop.

- Changing the target epitope or the anchor set
- Changing the affinity direction, the pairing rule, or any rejection criterion
- Starting a real EGFR campaign on Modal, as opposed to a smoke run
- Anything that contradicts an entry under Commitments in `docs/decisions-log.md`
- What goes into the submission, how many designs, and in what order
- Posting anything in the Proteinbase Slack. **The workspace is read-only for us**
- Reading another entrant's work. There is a standing decision not to, recorded
  under Commitments, and it is not yours to revisit

---

## Blocked — do not invent an answer

| Question | Why it is blocked | What to do |
|---|---|---|
| The exact mouse EGFR construct used in the screen — residue range, vendor, catalogue number | Asked twice in the Proteinbase Slack on 1 October 2026 and never answered | Use the sequence from the competition page, and say in any document that the screening construct is unconfirmed |
| Whether cynomolgus monkey cross-reactivity is in scope | The pre-launch announcement said mouse and cyno; the recorded objectives say human and mouse | Check the competition page. Do not assume either way. Decided 2 October 2026 to design for human and mouse only; see Settled in `docs/decisions-log.md` |
| What novelty level 3 requires exactly | Adaptyv stated it twice in terms that may not match — "level 3 of 4" and "under 30% similar to anything existing" | Read `adaptyvbio.com/blog/novelty` before relying on either |

Record anything newly blocked in the Open or Unverified section of
`docs/decisions-log.md`, which this table summarises rather than replaces.

---

## The affinity rule: maximise the pH gap

Design for the largest achievable difference between binding at pH 6.5 and binding
at pH 7.4, with affinity at pH 6.5 as high as the switch allows. This replaces an
earlier rule of holding affinity down deliberately; the withdrawal is recorded in
`docs/decisions-log.md` under Ruled out and explained in
`docs/explainers/07-the-assay-and-what-it-changes.md`.

What the organisers said, on 30 September (`docs/competition-qa-log.md`): designs
with a large K_D shift qualify even if they bind at both pH values, and designs with
no binding at pH 7.4 together with high affinity at pH 6.5 rank higher. The
instrument flows the target at a top concentration of 1000 nM and reports K_D over
roughly 0.1 nM to 10 µM.

What follows for a design:

- No binding at pH 7.4 is the best outcome and not the only qualifying one.
- A binder made deliberately weak gives up the top of the ranking. The fraction
  of binder occupied at the assay's top concentration is [analyte] divided by
  itself plus K_D — by definition, since K_D is the concentration at which half
  the sites are occupied. At 1000 nM that gives roughly 91% occupancy for a K_D
  of 100 nM, falling to 9% at 10 µM (`analysis/17_occupancy_table.py`,
  `results/findings/17-occupancy.md`), so a design weakened toward the far end
  of the reportable range can read as no detectable binding at pH 6.5 as well,
  which fails human binding.
- Charge-pair count leads the ranking, because the charge pairs are what produce the
  switch. BindCraft2's own confidence ordering is a tie-break among candidates with
  equal pair counts, used to prefer the higher-confidence interface.
- A change that raises binding at pH 7.4 as much as at pH 6.5 has not widened the
  gap. Say so when that happens.

We do not know how the organisers score a "large" K_D shift, as a ratio or as an
absolute difference, or where they draw the line. The design target above is a
judgement and should be described as one.

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
expensive mistake. It earns its place when it separates two things that get
confused at a cost, such as a binder that is weak and a binder that is selective:
the first can read as no detectable binding at either pH. Almost everywhere else it
is only emphasis dressed as precision.

**Avoid dramatisation.** State what happened and what follows from it. In
particular, make no claims about what other people usually do, what most teams skip,
what standard practice overlooks, or how consequential something is. We have no
basis for claims about other people's practice, and reaching for one is a reliable
sign that the sentence is decoration rather than content. (The one exception already
in this file is specific and checkable: BindCraft2 ranks its own output by an
interface-confidence score, which is a statement about a tool's objective function,
not about the people using it.)

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

A finished artefact gets an explainer even when the phase around it is still open,
if the reader will otherwise have to trust its output without understanding it. A
tool is understood before its results are believed, not after. Such an explainer
states in its own opening that it covers work with no results yet, and is rewritten
when there are.

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
`docs/competition-qa-log.md`, the explainers, and everything in `results/findings/`
and `data/derived/`.

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

**The current anchor set is the eight residues of the H370 patch:** E344, H358,
D368, H370, E391, E400, E421 and E424. Six are acidic and two — H358 and H370 — are
the target histidines the second row of the table applies to. These replaced the
415–466 set at step 08; that older set, whose target histidines were H418 and H433,
is the recorded fallback and is not what anything is currently aimed at.

**Reject any candidate with a histidine on the binder facing any histidine on the
target.** That pair switches off as pH drops, because both become positive and
repel, and it can cancel out a correctly built pair elsewhere. The rule applies to
every target histidine the binder faces, including ones outside the anchor set,
because the physical problem is identical. `analysis/10_charge_pair_filter.py`
implements it that way and has a test case against H418 to prove it. Reject such
candidates rather than scoring them.

**Also reject any candidate that contacts position 442.** S442G is the single
human/mouse difference inside the original 415–466 epitope and it sits in
cetuximab's measured contact set, so touching it risks species-specific behaviour at
the one position where the species differ. The rule stands whichever epitope is
current.

**The target carries a C-terminal histidine tag in the assay, and it will probably
not be cleaved.** An acidic pocket built to hold a target histidine will hold the
histidines of a tag just as readily, and a design that binds the tag looks
pH-selective while binding anything His-tagged. This is handled twice over: the tag
is given to BindCraft2 as an off-target inside the design loop
(`HIS_TAG_SEQUENCE`, `HIS_TAG_WEIGHT` and `DETARGET_IPTM_CEILING` in
`analysis/egfr_common.py`), and surviving candidates are audited afterwards by
`analysis/15_histag_counterscreen.py`. Both the weight and the ceiling are starting
values rather than tuned ones. Background: explainer 07.

**The binder's own C-terminus is kept out of the interface** by BindCraft2's
`termini_accessible` setting, because that end carries the tail the design is
immobilised by. It is set in both campaign files written by `analysis/09`.

**Every anchor sits within reach of a sugar-chain attachment point, and the demotion
this gets in the ranking depends on which one.** `analysis/14_glycan_sequons.py`
found two attachment points: N352, 10.4–14.3 Å from the five anchors E344, H358,
D368, H370 and E391, present in both species; and N361, 16.2–16.9 Å from E400, E421
and E424, a sequon in human and not in mouse (mouse has tyrosine there). Decided
2 October 2026: a chain at N352 costs both species alike and is not demoted, because
it leaves the mouse-to-human K_D ratio untouched; a chain at N361 is demoted,
because whatever shielding it causes happens on the human target only, which moves
that ratio directly. `analysis/10_charge_pair_filter.py` demotes — never excludes —
a design leaning on E400, E421 or E424, the same shape as the pre-existing
edge-reliance demotion of E424 near the trimmed target's cut edge, and the two stack
for E424 itself. This is a tie-break on a plausible asymmetry, not a measured
effect: nothing measures whether a chain at this distance actually reaches these
three anchors. Explainer: `docs/explainers/08-sugar-chains-near-the-anchors.md`.

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
