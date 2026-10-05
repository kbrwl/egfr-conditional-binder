"""20_submission_check.py -- refuse to let a malformed submission be uploaded.

Runs against `submission/challenge1-egfr-submission.csv` and must pass before any
upload. It is a separate program from `analysis/21_build_submission.py`, which
writes that file, because a builder that validates its own output tests only that
it is self-consistent. This one re-derives what the file should contain from the
screens' own tables and compares.

Abbreviations, expanded because this file gets read on its own:
  EGFR  epidermal growth factor receptor, the protein we design against.
  CSV   comma-separated values, the plain-text table format the organisers want.
  aa    amino acids, the unit a protein's length is counted in.

WHAT IT CHECKS, AND WHY EACH ONE IS HERE
-----------------------------------------
The first four come straight from the organisers' stated format
(`docs/rules-reference.md`). A file that fails any of them is rejected by the
platform, or worse, accepted and read as something we did not mean:

  1. exactly three columns, named and ordered `name`, `sequence`,
     `molecule_class`
  2. `molecule_class` drawn from single_chain / nanobody / scfv / fab_kappa /
     fab_lambda. **This set was wrong until 6 October 2026**: it read `protein`
     as the first value, taken from the execution plan rather than from the
     platform, and the platform rejected the upload. The check therefore passed
     a file the platform refused, because it encoded the same unverified
     assumption the file did. A check built from an assumption validates the
     assumption, not the file
  3. every sequence 10 to 250 aa, written only in the twenty standard amino
     acids. A lowercase letter, a gap character left in from an alignment, or
     an X standing for "unknown" all read as a sequence we did not design
  4. no duplicate sequences and no duplicate names. A duplicate sequence spends
     two of twenty slots on one design

The fifth is ours rather than theirs:

  5. the rows are in the order step 10 ranks them, with the designs step 13
     found clashing removed. This is the check that the uploaded file is the
     screens' output rather than something edited afterwards. It re-reads
     `results/candidates/shortlist.csv` and
     `data/derived/13-candidate-clashes.csv` and rebuilds the expected sequence
     order independently of the builder

The sixth reports rather than passes or fails:

  6. the pairwise identity matrix within the submission. Design diversity is one
     of the axes selection is scored on, so how alike our own designs are is a
     number we should know and state rather than discover afterwards. There is
     no threshold here because the organisers publish none; the figure goes in
     the write-up as measured

WHAT THE IDENTITY FIGURE MEANS, AND WHAT IT DOES NOT
-----------------------------------------------------
This is identity within our own submission, between our designs. It is not the
novelty check: that compares each design against public databases and is scored
by the platform on upload. A set of designs can be entirely novel against the
world and nearly identical to each other, which is exactly the risk when every
sequence comes from one backbone, and that is the case this number is here to
make visible.

Sequences of equal length are compared position by position. Unequal lengths are
aligned first with Needleman-Wunsch, a standard global alignment, because round
two draws binder lengths from a band and the comparison must not silently become
a comparison of prefixes.
"""
import argparse
import csv
import itertools
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUBMISSION = ROOT / "submission" / "challenge1-egfr-submission.csv"
SHORTLIST = ROOT / "results" / "candidates" / "shortlist.csv"
CLASHES = ROOT / "data" / "derived" / "13-candidate-clashes.csv"

COLUMNS = ["name", "sequence", "molecule_class"]
# Read from the platform's own upload form on 6 October 2026, after it rejected
# `protein`. Not taken from any document of ours.
MOLECULE_CLASSES = {"single_chain", "nanobody", "scfv", "fab_kappa",
                    "fab_lambda"}
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LENGTH, MAX_LENGTH = 10, 250

# Needleman-Wunsch scores. Identity is counted over aligned non-gap columns, so
# these only decide where the gaps fall, not what is then reported.
MATCH, MISMATCH, GAP = 1, -1, -2


def align(a, b):
    """Global alignment of two sequences. Returns the aligned pair."""
    rows, cols = len(a) + 1, len(b) + 1
    score = [[0] * cols for _ in range(rows)]
    for i in range(1, rows):
        score[i][0] = score[i - 1][0] + GAP
    for j in range(1, cols):
        score[0][j] = score[0][j - 1] + GAP
    for i in range(1, rows):
        for j in range(1, cols):
            score[i][j] = max(
                score[i - 1][j - 1] + (MATCH if a[i - 1] == b[j - 1] else MISMATCH),
                score[i - 1][j] + GAP,
                score[i][j - 1] + GAP)
    out_a, out_b, i, j = [], [], len(a), len(b)
    while i or j:
        if i and j and score[i][j] == score[i - 1][j - 1] + (
                MATCH if a[i - 1] == b[j - 1] else MISMATCH):
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i -= 1; j -= 1
        elif i and score[i][j] == score[i - 1][j] + GAP:
            out_a.append(a[i - 1]); out_b.append("-"); i -= 1
        else:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1
    return "".join(reversed(out_a)), "".join(reversed(out_b))


def identity(a, b):
    """Percent identity over aligned columns where neither side is a gap."""
    if len(a) != len(b):
        a, b = align(a, b)
    columns = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    if not columns:
        return 0.0
    return 100.0 * sum(1 for x, y in columns if x == y) / len(columns)


def expected_sequences():
    """Rebuild the order the submission should be in, from the screens' tables."""
    if not SHORTLIST.exists() or not CLASHES.exists():
        return None
    clashing = {row["design"] for row in csv.DictReader(CLASHES.open())
                if row["verdict"] == "clashing"}
    return [row["Binder_Sequence"].strip()
            for row in csv.DictReader(SHORTLIST.open())
            if row["design"] not in clashing]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--submission", type=Path, default=SUBMISSION)
    args = parser.parse_args()

    failures, notes = [], []
    if not args.submission.exists():
        sys.exit(f"missing {args.submission}. "
                 "Run analysis/21_build_submission.py first.")

    raw = args.submission.read_text()
    if "\r" in raw:
        failures.append("the file contains carriage returns; it should be "
                        "newline-terminated only")

    rows = list(csv.reader(args.submission.open(newline="")))
    if not rows:
        sys.exit("the submission file is empty")
    header, body = rows[0], [r for r in rows[1:] if r]

    print("=" * 72)
    # relative_to raises for any path outside the repository, which made
    # --submission unusable for a file in /tmp and turned a deliberate
    # break-test into a crash that looked like a rejection.
    try:
        shown = args.submission.relative_to(ROOT)
    except ValueError:
        shown = args.submission
    print(f"SUBMISSION CHECK  {shown}")
    print("=" * 72)
    print()
    print(f"1. Columns. Expect exactly {COLUMNS}.")
    if header != COLUMNS:
        failures.append(f"header is {header}, expected {COLUMNS}")
        print(f"   [FAIL] {header}")
    else:
        print(f"   [PASS] {header}")
    print()

    print(f"2. Rows. {len(body)} designs, each with three fields.")
    for index, row in enumerate(body, start=1):
        if len(row) != 3:
            failures.append(f"row {index} has {len(row)} fields, expected 3")
    print(f"   [{'PASS' if all(len(r) == 3 for r in body) else 'FAIL'}] "
          f"every row has three fields")
    print()

    names = [r[0] for r in body if len(r) == 3]
    sequences = [r[1].strip() for r in body if len(r) == 3]
    classes = [r[2] for r in body if len(r) == 3]

    print("3. Molecule class. One of "
          + " / ".join(sorted(MOLECULE_CLASSES)) + ".")
    bad = {c for c in classes if c not in MOLECULE_CLASSES}
    if bad:
        failures.append(f"molecule_class values not permitted: {sorted(bad)}")
    print(f"   [{'FAIL' if bad else 'PASS'}] "
          + ", ".join(f"{c}×{classes.count(c)}" for c in sorted(set(classes))))
    print()

    print(f"4. Sequences. {MIN_LENGTH}-{MAX_LENGTH} aa, twenty standard amino "
          "acids only.")
    for name, seq in zip(names, sequences):
        if not MIN_LENGTH <= len(seq) <= MAX_LENGTH:
            failures.append(f"{name} is {len(seq)} aa, outside "
                            f"{MIN_LENGTH}-{MAX_LENGTH}")
        odd = sorted(set(seq) - STANDARD_AA)
        if odd:
            failures.append(f"{name} contains non-standard characters: {odd}")
    lengths = sorted(set(len(s) for s in sequences))
    print(f"   [{'PASS' if not failures else 'see below'}] lengths "
          + (f"all {lengths[0]} aa" if len(lengths) == 1
             else f"{min(lengths)}-{max(lengths)} aa"))
    print()

    print("5. No duplicates.")
    dup_names = {n for n in names if names.count(n) > 1}
    dup_seqs = {s for s in sequences if sequences.count(s) > 1}
    if dup_names:
        failures.append(f"duplicate names: {sorted(dup_names)}")
    if dup_seqs:
        failures.append(f"{len(dup_seqs)} duplicate sequence(s)")
    print(f"   [{'FAIL' if dup_names or dup_seqs else 'PASS'}] "
          f"{len(set(names))} distinct names, {len(set(sequences))} distinct "
          f"sequences, from {len(names)} rows")
    print()

    print("6. Order. Rows must follow step 10's ranking with step 13's "
          "clashing designs removed.")
    expected = expected_sequences()
    if expected is None:
        notes.append("the screens' tables were not found, so the order could "
                     "not be checked independently")
        print("   [SKIP] shortlist.csv or 13-candidate-clashes.csv missing")
    else:
        want = expected[:len(sequences)]
        if sequences != want:
            failures.append("the submission order does not match step 10's "
                            "ranking after step 13's removals")
            for index, (got, wanted) in enumerate(zip(sequences, want), start=1):
                if got != wanted:
                    print(f"   [FAIL] row {index} is not what the ranking gives")
                    break
        else:
            print(f"   [PASS] all {len(sequences)} rows in ranked order")
    print()

    print("7. Pairwise identity within the submission.")
    print("   Reported, not graded. This is how alike our own designs are, "
          "which\n   bears on the design-diversity axis selection is scored on. "
          "It is not\n   the novelty check, which compares against public "
          "databases.")
    print()
    if len(sequences) > 1:
        print("        " + "".join(f"{i + 1:>6}" for i in range(len(sequences))))
        for i, a in enumerate(sequences):
            cells = "".join(
                "     -" if i == j else f"{identity(a, b):>6.0f}"
                for j, b in enumerate(sequences))
            print(f"   {i + 1:>3}  {cells}")
        pairs = [identity(a, b) for a, b in itertools.combinations(sequences, 2)]
        print()
        print(f"   minimum {min(pairs):.0f}%   maximum {max(pairs):.0f}%   "
              f"mean {sum(pairs) / len(pairs):.0f}%   over "
              f"{len(pairs)} pairs")
        if min(pairs) > 70:
            print()
            print("   Every pair is above 70% identical to every other. These "
                  "are variants\n   of one backbone, not independent designs, "
                  "and the write-up should say\n   so rather than let the "
                  "count of seven imply seven tries.")
    else:
        print("   only one sequence; nothing to compare")
    print()

    print("=" * 72)
    if failures:
        print(f"SUBMISSION CHECK FAILED: {len(failures)} problem(s).")
        for failure in failures:
            print(f"  - {failure}")
        print("Do not upload this file.")
    else:
        print(f"SUBMISSION CHECK PASSED. {len(sequences)} designs, "
              "format and order both correct.")
    for note in notes:
        print(f"  note: {note}")
    print("=" * 72)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
