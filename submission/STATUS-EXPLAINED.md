# Where this project actually is — 3 October 2026

Written because the work has sprawled across several sessions and the
running account of it has been unreadable. No jargon goes unexplained here.
Everything below is either a measured number from this repository or is
labelled as a judgement.

---

## Read this bit first: you have ten sequences

You said there is not a single binder sequence in front of you. There are ten.
They were produced today and nobody put them where you could see them. Here
they are. Each is a chain of 48 amino acids, written in the standard one-letter
code where each letter is one amino acid.

| # | sequence (48 amino acids) | charge pairs | fits? | verdict |
|---|---|---|---|---|
| 1 | `DPEEEERFWKWIHEQLREAPELPPHERASFASRALTEFQDPANRSRFL` | 2 | yes | **REJECTED** (rule broken) |
| 7 | `DPAEEEEFWKYVHEKLREAPPLPPHERASLASRLLTEFQQPEKRSEFL` | 2 | yes | **REJECTED** (rule broken) |
| 2 | `SPEEEEKFWKYVHEKLRKMPELPPHERAEFASRLLTEFQKPENRSKFL` | 1 | yes | usable |
| 8 | `SPEEEEKFWKYIHEKLREAPELPPHERAEFASRLLTEFQKPENRSKFL` | 1 | yes | usable |
| 3 | `SPEEEERFWKYIHEQLRKAPPLPPHERAEFASRLLTEFQDPANRSRFL` | 0 | yes | usable |
| 4 | `SPEEEEQFWKYVHERLREMPPLPPHERAEFASRLLTEFQQPEKRSEFL` | 0 | yes | usable |
| 5 | `DPEEEEQFWKYIHEQLRAMPPLPPHERASFASRKLTEFQNPENRSRFL` | 0 | no | **REJECTED** (collides) |
| 6 | `DPKEEEEFWKYVHEKLRKMPPLPPHERASLASRLLTEFQDPKNRSKFL` | 0 | yes | usable |
| 9 | `SPEEEEKFWKWIHEKLRAMPPLPPHERAEFASRALTEFQDPAKRSEFL` | 0 | yes | **REJECTED** (rule broken) |
| 10 | `SPEEEERFWKYVHELLREAPELPPHERAELASRLLTEFQDPANRSRFL` | 0 | yes | **REJECTED** (rule broken) |

**Five of those ten are usable right now.** They could be submitted today.
Four break a charge rule we set ourselves and one collides with part of EGFR
that was not in the picture when it was designed. All explained below.

They are also saved as a plain file at `submission/round1-sequences.fasta`.

---

## 1. What we are trying to make, in plain words

**The target.** EGFR is a protein that sits on the surface of human cells, poking
out of them. It is involved in cell growth, and in many cancers there is too much
of it. It is one of the best studied drug targets there is.

**The thing we are designing.** A "binder" is a small protein that sticks to a
chosen patch on a bigger protein, the way a key fits one lock. We are designing a
small one: around 50 amino acids, where an amino acid is one link in a protein
chain.

**The twist, and the whole point of the competition entry.** Our binder is
supposed to stick to EGFR *only when the surroundings are slightly acidic*, and let
go when they are not.

Why that is worth doing: the fluid inside a solid tumour is more acidic than the
fluid in healthy tissue. Tumour roughly pH 6.5, healthy tissue roughly pH 7.4,
where pH is the standard acidity scale and lower means more acidic. Healthy cells
carry EGFR too. So a binder that cannot tell the two apart hits healthy tissue as
well. One that only switches on in acid concentrates itself where the acidity is.

**How we make it switch.** Of the twenty amino acids, exactly one changes its
electrical charge between those two acidities: histidine, written H. It is neutral
at pH 7.4 and becomes positive at pH 6.5. Every other amino acid is the same at
both.

So the design works like this. Opposite charges attract. We find spots on EGFR
that carry a permanent negative charge, and we place a histidine on our binder
facing each one. At pH 7.4 our histidine is neutral, nothing much happens, and the
binder drifts off. At pH 6.5 it turns positive, is pulled towards the negative spot
on EGFR, and the binder sticks. Stack several of these pairs and the effect adds up.

There is a second version of the same trick. Where EGFR itself has a histidine, we
put a negative amino acid on our binder facing it. Same result, same direction: in
acid, EGFR's histidine turns positive and is pulled to our negative one.

Those are the **"charge pairs"** referred to everywhere in this project. The count
of them is the single most important quality number we have, because they are what
produces the switch.

**One rule that follows from this, and it matters later.** Never put a histidine on
our binder facing a histidine on EGFR. In acid they both turn positive, and two
positives push each other apart. That one pair actively works against every correct
pair elsewhere. Any design that does it gets thrown out.

---

## 2. Which patch of EGFR we aim at, and why that was a real decision

EGFR's outside portion is big, roughly 620 amino acids. We cannot aim at all of
it; a 50-amino-acid binder covers a small patch. Choosing the patch was most of the
early work.

**We aim at a patch built around EGFR's own histidine number 370**, written H370.
Numbers like 370 are positions counted along EGFR's chain from one end.

Why that patch:

- A published paper (Liu and colleagues, 2022) tested which of EGFR's histidines
  are responsible for acid-dependent binding, by removing them one at a time and
  seeing what broke. Two mattered: H370 and H433. Building on a position already
  shown by experiment to drive acid-dependent behaviour beats guessing.
- Of those two, H433 sits underneath where cetuximab binds. Cetuximab is an
  approved cancer antibody. Designing into the footprint of an existing drug is a
  worse position to be in than designing somewhere unoccupied. H370 is not touched
  by it.
- The surrounding spots had to be *conserved*, meaning identical in mouse as well
  as human, because the competition requires the same binder to work on mouse EGFR
  too. They are.

**The eight anchor positions** chosen around H370: E344, H358, D368, H370, E391,
E400, E421, E424. E and D are the permanently negative amino acids; H is histidine.
So six of them get a histidine on our binder, and two of them get a negative.

**Then we cut it down to four.** That is the single most consequential change in the
whole project and it came from a measurement, not an opinion. Explained in section 4.

---

## 3. The machine that makes the designs, and the two different scoreboards

**BindCraft2** is the open-source tool that actually invents the binders. We rent a
graphics card by the hour from a company called Modal and run it there, because
this kind of work needs that hardware.

Roughly, for each attempt it: invents a shape, predicts how that shape and EGFR
would sit together, scores how convincing the result looks, nudges it, repeats.
One such attempt is called a **trajectory**. Each trajectory costs real money,
around 25 to 30 cents.

A trajectory runs through five named stages in order — screen, refine, anneal,
harden, mutate. **At the end of each stage there is a minimum score. Miss it and
the trajectory is killed on the spot and the money is spent for nothing.** Most
trajectories die partway. That is normal, not a malfunction.

**Here is the thing that has confused this project for two days, and it is worth
getting straight.** There are two completely separate scoreboards:

**Scoreboard one — BindCraft2's own.** Its main number is called i_pTM, which is
the tool's confidence that the two proteins really would sit together the way it
drew them. It runs 0 to 1. BindCraft2's own default is that a design must reach
0.7 to be "accepted". *This number has nothing to do with acidity.* BindCraft2 has
no concept of pH anywhere in it.

**Scoreboard two — ours.** We count charge pairs, as described in section 1, and we
apply our own rejection rules. This is the scoreboard that measures the thing we
actually care about.

**"0 accepted designs" has always meant scoreboard one.** It never meant "no
sequences exist". Those are different statements and conflating them is why this
has felt like total failure. BindCraft2 declining to stamp a design as high
confidence is not the same as the design being useless to us — especially since
the competition organisers have said in writing that their in-silico scoring is
"only a smaller component" of how entries are judged, and that weak confidence
scores on an awkward target will not be penalised much.

---

## 4. Everything that has been run, what it cost, and what it taught

In order. Every number here is read from the saved output of the run itself.

### The first five campaigns — 26 attempts, before today

| campaign | attempts | what it was testing | result |
|---|---|---|---|
| validate | 1 | does the setup start at all | it starts |
| probe (ext499) | 1 | would a longer piece of EGFR behave better | no improvement |
| diagnose | 10 | eight anchors, one filter loosened | nothing reached the end |
| floors | 8 | eight anchors, confidence floors lowered | nothing reached the end |
| tight4 | 6 | **four anchors instead of eight** | **two reached the last stage** |

Cost of those: roughly $8.

**The finding that mattered, and it was being hidden by how we reported it.** Added
together, those five campaigns read "26 attempts, 0 accepted", which sounds like
nothing works. Split apart, they say something quite different:

- **Eight anchors: 0 out of 20 attempts reached the final stage.**
- **Four anchors: 2 out of 6 did.**

The likely reason: eight anchor points are spread over about 24 angstroms, which
is a wider stretch than a 50-amino-acid binder can cover well. Asking for all
eight produces a binder that reaches for everything and grips nothing. Four
neighbouring ones are an achievable ask.

Pooling those campaigns into one number made the project look stuck when the
narrowed version was working. That is now reported split, permanently.

### A mistake that got caught, and is worth you knowing about

Four conclusions drawn earlier from one or two attempts turned out to be wrong when
more attempts were run. The clearest: it was argued that the target-confidence
filter would block essentially every attempt, based on how little that number
varied across ten *accepted* designs from an unrelated test. But a set of accepted
designs has already been filtered by the very thing being measured, so it cannot
tell you the spread of the unfiltered population. The same piece of EGFR later
scored 0.72 where it had scored 0.35.

The rule that came out of it, and it governs everything since: **on this target,
fewer than about ten attempts settles nothing.** The spread is too wide.

### Today's run, which is where your ten sequences came from

Two campaigns side by side, ten attempts each. Everything held identical between
them except **the length of the binder**, so that any difference could only be
caused by length. Both used the four-anchor patch that had worked.

| | short: 30-60 amino acids | long: 60-100 amino acids |
|---|---|---|
| attempts | 10 | 10 |
| ran all the way through | **1** | **0** |
| sequences produced | **10** | **0** |
| passed BindCraft2's own bar | 0 | 0 |
| cost | $2.89 | $3.34 |
| peak memory used on the card | 8,894 MB | 17,227 MB |

**The short band works and the long band does not.** Ten against ten is a real
result, not a fluke. The long band never produced a single structure at any point.
Its confidence scores while working ran 0.21 to 0.25; the short band's ran 0.73 to
0.83. A 60-to-100 amino acid binder is not finding this patch.

**A bonus finding nobody was looking for.** The short band uses half the memory.
The rented card has 23,034 MB. Two short jobs fit on one card at once; two long
ones do not. So the short band is not just better, it can be run at double the rate
for the same rental.

### One fix that was essential and nearly missed

The earlier tight4 run reached the last stage twice and **saved nothing**. The
setting that writes sequences to disk was switched off by default. Two successful
attempts, and the output was thrown away. That setting was turned on for today's
run, which is the only reason you have sequences at all.

### Four faults in our own scoring code, found today

The scoring scripts had never once seen real design output. The moment they did,
four separate faults appeared. Each one would have stopped us at the deadline:

1. Our structure reader crashed on every file BindCraft2 produced, because
   BindCraft2 leaves out a column the reader demanded. Nothing could be scored.
2. A column-name collision crashed the results writer.
3. Files describing the binder against the purification tag were being scored as
   though they were EGFR, producing ten meaningless rows.
4. **The worst one.** The tag safety check reported "not recorded" for all ten
   designs while the reading was sitting right there in the table under a
   different name. A safety check that fails by reporting *nothing* is the most
   dangerous kind, because ten blanks read like ten passes.

Two further faults in the record-keeping script were fixed after that: it was
counting a completed attempt as an accepted design, and each campaign's records
were overwriting the previous campaign's.

---

## 5. How good are the ten sequences, honestly

They went through all three of our checks. This is the first time those checks
have ever run on real output.

**Check one — charge pairs.** Our target was three or more correct pairs.

- **Nought of ten reach three.** The best two have two pairs each.
- **Those best two are both disqualified anyway**, because they touch EGFR
  position 442. We banned that position early: it is the one spot in the region
  where human and mouse EGFR genuinely differ, so touching it risks a binder that
  behaves differently in the two species, and mouse cross-reactivity is a scored
  requirement.
- Two more are disqualified for the histidine-facing-histidine error from
  section 1, against EGFR's H433.
- That leaves six by this check alone. One of those six then fails the
  physical-fit check below, so **five are usable overall**.

**Check two — does it physically fit.** We designed against a cut-out piece of
EGFR. This check puts the design back against the whole receptor to see whether it
would collide with parts that were not in the picture. **Nine of ten fit. One
collides.**

**Check three — does it stick to the purification tag instead.** The EGFR used in
the competition's lab test carries a tag made of histidines. Our binder is
*designed* to have pockets that grab histidines. So there is a real danger of
building something that grabs any tagged protein and looks acid-sensitive while
being useless. **All ten pass comfortably**: scores 0.08 to 0.11 against a limit
of 0.4. This is the strongest result we have and it holds up on real output.

**Novelty.** The competition requires designs not to resemble existing proteins.
Searching all ten against Swiss-Prot, the curated database of known proteins, found
**no meaningful resemblance to anything**. The raw similarity percentages look
highish, 31 to 42%, but every single match is statistical noise — short random
sequences coincidentally resemble things in a database of 575,748 entries all the
time. No match came anywhere near significance.

### So what is the real verdict

**The good:** sequences exist, they are the right length, they are novel, they
avoid the tag trap, they almost all physically fit, and the hard safety rules are
visibly catching real violations rather than sitting unused.

**The bad, stated plainly:** the charge-pair counts are poor. We wanted three and
got a maximum of two, and the two that managed two are disqualified on other
grounds. The pH switch is the entire idea of this entry, and on this evidence it is
weakly built.

**The important caveat in both directions:** *all ten sequences came from a single
successful attempt.* They are ten variations on one backbone, not ten independent
designs. So "0 out of 10 reached three pairs" is really one backbone's answer
repeated ten times. It is much weaker evidence than it sounds — in both directions.

---

## 6. Money

| | |
|---|---|
| free monthly allowance from Modal | $30 |
| spent before today | about $8 |
| spent today on the two campaigns | **$6.23** |
| left | **about $13 to $15** |

The allowance does not refill until 1 November. There are no competition credits
coming: those went to Track 1 and Track 2 entrants and we are Track 3.

At the short band's measured rate of about $0.28 an attempt, $13 buys roughly **45
more attempts** — and closer to 90 if two are run per card, which the memory
measurement says is possible for the short band.

Today's spend came in under its $7 ceiling, and the ceiling was enforced by a
hard timeout rather than by anyone watching.

---

## 7. The deadline

**6 October, 23:59 "anywhere on earth".** Extended from 4 October because some
teams got their Modal credits late. Anywhere-on-earth means it has not passed until
it has passed in every timezone, so in practice around midday UK time on 7 October.

**That is three days away, not hours.** Submission rules: up to 20 designs per
submission, one submission per 24 hours, so about three more windows. An early
submission can be replaced by a better one later.

**Nothing has been submitted yet.** That is the single most important open item,
and it is waiting on you, not on more work.

---

## 8. What I think you should decide

Three things. My recommendation on each, with the reasoning, so you can disagree
with the reasoning rather than just the conclusion.

### A. Submit the five usable sequences now, today?

**I would say yes, with the five that pass everything.** The reasoning: it costs nothing, a submission can be replaced
in 24 hours, and it flushes out any problem with the submission format while there
is still time to fix it. Going into the last day having never once tested the
upload is an avoidable risk.

The argument against, which is real: all five are variations on one backbone, so
it spends a submission window on a single idea, and diversity is one of the three
things entries are judged on.

**My actual suggestion: submit them today as a safety net, and plan to replace
them.** A held submission that is never bettered scores zero.

### B. How to spend the remaining money

**Short band only.** The long band is finished as an idea — nought out of ten. I
would spend the next round on getting *different backbones*, because the single
biggest weakness in what we have is that everything descends from one.

### C. The one thing I would change in the design itself

Our best designs reach two charge pairs, not three. We have never tried telling
BindCraft2 to concentrate harder on the four anchor points. There is a setting for
that. It was ruled out earlier for a specific technical reason involving the tag
being supplied as a bare sequence, and that reason should be re-checked now rather
than taken on trust, because two pairs versus three is exactly the gap it might close.

---

## 9. Straight answers to what you actually asked

**"There is not a single binder sequence in front of me."**
There are ten, at the top of this file and in `submission/round1-sequences.fasta`.
They were produced today at about 16:00. They existed for several hours before
anyone showed them to you, and that is a reporting failure rather than a work one.

**"So much time and money and nothing to show."**
What the money bought, stated fairly:
- The eight-anchor patch does not work and the four-anchor one does. That took
  26 attempts to establish and it is the finding the rest rests on.
- Long binders do not work on this patch and short ones do. Twenty attempts.
- Ten real sequences that pass three of our four quality checks.
- A scoring pipeline that now actually runs. It had six faults, all invisible
  until real output existed, several of which would have been fatal on the last day.
- The acid-switch idea is reachable but weakly: two pairs achieved against a
  target of three.

What it did not buy: a design we can be confident switches strongly with acidity.
That is the honest gap.

**"You keep asking my input and I have no idea what is going on."**
Fair. The three decisions in section 8 are the only ones that need you, and
section 8 gives a recommendation for each. If you would rather not decide: say
"submit the five and carry on", and that is enough to proceed.

---

## 10. Where things live, if you want to look

| what | where |
|---|---|
| **the ten sequences** | `submission/round1-sequences.fasta` |
| their full scores | `data/derived/10-candidate-summary.csv` |
| physical-fit check | `data/derived/13-candidate-clashes.csv` |
| tag safety check | `data/derived/15-tag-screen.csv` |
| every attempt ever made, per campaign | `data/derived/18-trajectory-ledger-*.csv` |
| the full decision history | `docs/decisions-log.md` (long) |
| methods write-up, draft | `submission/methods.draft.md` |
| competition facts | `docs/rules-reference.md` |

---

## 11. One thing I want to flag about my own reporting

Several times today I described progress using the phrase "0 accepted designs"
without immediately saying that it refers to BindCraft2's internal confidence bar
and not to whether sequences exist. Those are different things, and running them
together is a large part of why this has felt like hours of work producing
nothing. The sequences were real and scoreable for some time before that was
stated plainly. Worth recording so it does not repeat.

