#!/usr/bin/env python3
"""
01_alignment.py — compare human and mouse EGFR. This is a regression test.

WHAT THIS MEASURES AND WHY IT MATTERS
-------------------------------------
One of the competition's objectives, mouse cross-reactivity, is that a single
binder sequence must bind both human and mouse EGFR. The cheapest way to get
there is to aim at a patch of the target where the two species already have the
same residues (a residue is one building block of a protein), rather than
designing for human first and fixing up the mouse case afterwards. If the binder
never touches a position where the two species differ, it works on both by
construction instead of by luck.

Finding such a patch means knowing exactly where the two proteins differ, and
that needs an alignment: the two sequences laid side by side position by
position, with gaps inserted where one of them has extra residues. Comparing
position 1 to position 1 straight down the sequences does not work, because a
single insertion early on shifts everything after it and makes every later
comparison wrong.

Method: global pairwise alignment, meaning the two sequences are aligned end to
end rather than only at their best-matching fragments. Scoring uses BLOSUM62, a
table of how chemically similar each pair of residues is, with a penalty of -11
for opening a gap and -1 for extending one. BLOSUM62 records that some swaps are
chemically mild (serine to glycine) and others drastic, so the alignment prefers
biologically sensible arrangements. The gap penalties are the standard defaults:
opening a gap is expensive and extending one already open is cheap, which
reflects that a real insertion is usually one event several residues long.

This script is a regression test. The expected answer is already known from
earlier work, so the script asserts it:
  - exactly 16 differences in domain III (residues 310-480)
  - exactly one difference inside 415-466, namely S442G: human has serine at
    position 442 where mouse has glycine
If either assertion fails, the environment or the input data is wrong, and
nothing downstream should be trusted until that is fixed.

Outputs:
  results/findings/01-alignment.md
  data/derived/01-domain3-differences.csv
  data/derived/01-identical-runs.csv

Run standalone:  python analysis/01_alignment.py
Exit code 0 means the expectations above were met.
"""

import sys
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "data" / "sequences" / "egfr-uniprot-full.fasta"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"

# All ranges are positions in the full human record P00533 in UniProt, the public
# protein sequence archive.
ECD_START, ECD_END = 25, 645      # the mature part of EGFR that sits outside the cell
D3_START, D3_END = 310, 480       # domain III, the working definition we use
EPI_START, EPI_END = 415, 466     # candidate epitope: the patch we aim the binder at

MIN_RUN = 10                       # report identical runs of at least this length

# The answers this script asserts, carried over from earlier work.
EXPECTED_D3_DIFFS = [
    "A313P", "S315Y", "M318V", "V323I", "E330D", "S348T", "N361Y", "S364A",
    "R377K", "H383R", "Q390R", "D393E", "E412D", "R414W", "S442G", "K467R",
]
EXPECTED_EPITOPE_DIFFS = ["S442G"]


def read_fasta(path):
    records, header, chunks = {}, None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks)
    return records


def find(records, token):
    for header, seq in records.items():
        if token in header:
            return seq
    raise KeyError(f"no record matching {token!r}")


def main():
    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    failures = []

    def check(label, ok, detail):
        emit(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")
        if not ok:
            failures.append(label)

    records = read_fasta(FASTA)
    human = find(records, "EGFR_HUMAN")
    mouse = find(records, "EGFR_MOUSE")

    emit("=" * 72)
    emit("HUMAN vs MOUSE EGFR ALIGNMENT  (regression test)")
    emit("=" * 72)
    emit()
    emit(f"human P00533: {len(human)} aa (amino acids, the building blocks)")
    emit(f"mouse Q01279: {len(mouse)} aa")
    emit("method: end-to-end pairwise alignment, BLOSUM62 table, gap open -11, gap extend -1")
    emit("numbering: positions in the full human UniProt record throughout")
    emit()

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignment = aligner.align(human, mouse)[0]

    # Walk the alignment column by column, numbering everything by its human
    # position. A column where the human side is a gap has no human position to
    # assign, so it is left out of the numbering and counted only in the tally of
    # mouse insertions reported below.
    rows = []            # (human_pos, human_aa, mouse_aa, identical)
    human_pos = 0
    human_gaps_in_mouse = 0
    for h_char, m_char in zip(alignment[0], alignment[1]):
        if h_char != "-":
            human_pos += 1
            rows.append((human_pos, h_char, m_char, h_char == m_char))
        else:
            human_gaps_in_mouse += 1

    def identity(lo, hi):
        sel = [r for r in rows if lo <= r[0] <= hi]
        n_same = sum(1 for r in sel if r[3])
        return n_same, len(sel), 100.0 * n_same / len(sel)

    emit("1. Identity by region")
    emit()
    emit("   | Region | Positions | Identical | Identity |")
    emit("   |---|---|---|---|")
    for label, lo, hi in [
        ("Extracellular region", ECD_START, ECD_END),
        ("Domain III", D3_START, D3_END),
        ("Candidate epitope", EPI_START, EPI_END),
    ]:
        same, total, pct = identity(lo, hi)
        emit(f"   | {label} | {lo}-{hi} | {same}/{total} | {pct:.1f}% |")
    emit()
    emit(f"   Columns where human has a gap (mouse insertions): {human_gaps_in_mouse}")
    emit()

    # Domain III differences, in human UniProt numbering.
    d3_diffs = [(r[0], r[1], r[2]) for r in rows
                if D3_START <= r[0] <= D3_END and not r[3]]
    d3_labels = [f"{h}{p}{m}" for p, h, m in d3_diffs]

    emit("2. Domain III differences, human UniProt numbering")
    emit("   Notation: Q390R means human has Q at position 390 where mouse has R.")
    emit()
    for i in range(0, len(d3_labels), 8):
        emit("   " + "  ".join(d3_labels[i:i + 8]))
    emit()
    before = sum(1 for p, _, _ in d3_diffs if p < EPI_START)
    emit(f"   {len(d3_labels)} differences total; {before} fall before position "
         f"{EPI_START}, {len(d3_labels) - before} at or after.")
    emit()

    epi_diffs = [f"{h}{p}{m}" for p, h, m in d3_diffs
                 if EPI_START <= p <= EPI_END]
    emit(f"3. Differences inside the candidate epitope {EPI_START}-{EPI_END}")
    emit(f"   {epi_diffs if epi_diffs else 'none'}")
    emit()

    # Stretches where the two species match residue for residue. These are found
    # across the whole extracellular region first and filtered to domain III only
    # afterwards, so a stretch that crosses the domain III boundary is not cut short.
    runs = []
    start = None
    for r in rows:
        if not (ECD_START <= r[0] <= ECD_END):
            continue
        if r[3]:
            if start is None:
                start = r[0]
        else:
            if start is not None and r[0] - start >= MIN_RUN:
                runs.append((start, r[0] - 1))
            start = None
    if start is not None and ECD_END - start + 1 >= MIN_RUN:
        runs.append((start, ECD_END))

    d3_runs = [(a, b) for a, b in runs if a <= D3_END and b >= D3_START]
    d3_runs.sort(key=lambda ab: ab[1] - ab[0], reverse=True)

    emit(f"4. Identical runs of {MIN_RUN}+ consecutive residues overlapping domain III")
    emit()
    emit("   | Positions | Length | Sequence |")
    emit("   |---|---|---|")
    for a, b in d3_runs:
        seq = human[a - 1:b]
        emit(f"   | {a}-{b} | {b - a + 1} | {seq} |")
    emit()

    emit("5. Regression assertions")
    check("domain III difference count", len(d3_labels) == 16,
          f"{len(d3_labels)} (expected 16)")
    check("domain III difference identities", d3_labels == EXPECTED_D3_DIFFS,
          "match expected list" if d3_labels == EXPECTED_D3_DIFFS
          else f"does not match the expected list: got {d3_labels}")
    check(f"epitope {EPI_START}-{EPI_END} differences",
          epi_diffs == EXPECTED_EPITOPE_DIFFS,
          f"{epi_diffs} (expected {EXPECTED_EPITOPE_DIFFS})")
    emit()

    emit("6. What this means for the design")
    emit()
    emit(f"   The {EPI_START}-{EPI_END} block is 52 positions with a single")
    emit("   difference, S442G: human serine, mouse glycine. Both are among the")
    emit("   smallest amino acids, so the local shape barely changes, which is about")
    emit("   as close to species-identical as a real surface patch gets.")
    emit()
    emit("   Fourteen of the sixteen domain III differences fall before 415, which")
    emit("   is the argument for aiming here rather than elsewhere in domain III: a")
    emit("   binder confined to this block satisfies mouse cross-reactivity by")
    emit("   construction.")
    emit()
    emit("   What this result does not settle: it is sequence analysis, so it tells")
    emit("   us which residue sits at each position and nothing about which")
    emit("   direction that residue points. Domain III folds into a solenoid, a")
    emit("   spiral-staircase shape in which residues next to each other in the")
    emit("   sequence can point opposite ways. Whether these anchors are actually")
    emit("   reachable is decided by steps 02-06.")
    emit()

    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "01-domain3-differences.csv").open("w") as fh:
        fh.write("uniprot_pos,challenge_pos,human_aa,mouse_aa,label,inside_epitope\n")
        for p, h, m in d3_diffs:
            fh.write(f"{p},{p - 24},{h},{m},{h}{p}{m},"
                     f"{EPI_START <= p <= EPI_END}\n")
    with (DERIVED / "01-identical-runs.csv").open("w") as fh:
        fh.write("start_uniprot,end_uniprot,length,sequence\n")
        for a, b in d3_runs:
            fh.write(f"{a},{b},{b - a + 1},{human[a - 1:b]}\n")
    emit("Wrote data/derived/01-domain3-differences.csv")
    emit("Wrote data/derived/01-identical-runs.csv")

    emit()
    emit("=" * 72)
    if failures:
        emit("RESULT: REGRESSION TEST FAILED -- STOP. Do not run 02 onward.")
        for f in failures:
            emit(f"  - {f}")
    else:
        emit("RESULT: REGRESSION TEST PASSED. Environment reproduces prior work.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "01-alignment.md").write_text(
        "# Human vs mouse EGFR alignment (computed)\n\n"
        "Computed output of `analysis/01_alignment.py`. Do not hand-edit.\n\n"
        "This is a regression test: the expected answer was known before the script\n"
        "was written, and the script checks itself against it. If it ever disagrees,\n"
        "the environment or the input data has changed and nothing computed after it\n"
        "should be trusted until that is sorted out. See `alignment-findings.md` in\n"
        "this directory for what the numbers mean for the design.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
