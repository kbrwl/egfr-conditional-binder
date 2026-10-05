"""21_build_submission.py -- write the submission CSV from the screens' own output.

The submission file is generated here and never typed by hand. That is the same
rule `results/` and `data/derived/` live under, and it applies with more force to
the one file that actually gets uploaded: a sequence transcribed by hand into a
submission is a sequence nobody screened.

Abbreviations, expanded because this file gets read on its own:
  EGFR  epidermal growth factor receptor, the protein we design against.
  CSV   comma-separated values, the plain-text table format the organisers want.
  aa    amino acids, the unit a protein's length is counted in.

WHAT GOES IN, AND IN WHAT ORDER
-------------------------------
Two inputs, both computed:

  results/candidates/shortlist.csv        every candidate step 10 did not reject,
                                          already in step 10's ranking order
  data/derived/13-candidate-clashes.csv   whether each candidate would fit on the
                                          intact receptor rather than the trimmed
                                          fragment it was designed against

Step 10 rejects on the hard rules -- a binder histidine facing a target histidine
is the only one still standing -- and ranks what survives. Step 13 is a separate
question: a design can pair correctly against the 171-residue fragment and still
be unable to occupy that space on the whole 621-residue receptor, because the
fragment does not carry the rest of the protein around it. A design that cannot
physically be there cannot bind, so a hard overlap is a removal and not a
demotion. It is the one exclusion applied here that step 10 does not already make.

The order is step 10's and is not recomputed. This script selects and renames; it
does not rank. If the order looks wrong, the ranking rule is in
`analysis/10_charge_pair_filter.py` and that is where it changes.

NAMES
-----
`egfr-ph-h370-NN`, numbered from the top of the ranking. The name carries the
target face rather than the campaign or the trajectory, because the organisers
read the name and our internal trajectory identifiers mean nothing to them. The
mapping back to the design that produced each row is written beside the CSV in
`submission/submission-provenance.csv`, so no row in the uploaded file is
untraceable.

MOLECULE CLASS
--------------
`protein` for all current rows. The organisers' permitted values are protein,
nanobody, scfv, fab_kappa and fab_lambda, and ours are single-chain designed
miniproteins of 48 aa, which is the minibinder category (40-100 aa) and the
`protein` class. This is asserted rather than assumed: any row outside 40-100 aa
stops the run so the class is reconsidered deliberately.

WHAT THIS DOES NOT DO
---------------------
It does not validate the file it writes. `analysis/20_submission_check.py` does
that, as a separate program, because a builder that passes its own output proves
nothing. Run 20 after this and before any upload.
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHORTLIST = ROOT / "results" / "candidates" / "shortlist.csv"
CLASHES = ROOT / "data" / "derived" / "13-candidate-clashes.csv"
SUBMISSION = ROOT / "submission" / "challenge1-egfr-submission.csv"
PROVENANCE = ROOT / "submission" / "submission-provenance.csv"

# The organisers take up to 20 designs in Track 3.
MAX_DESIGNS = 20

# The minibinder band, 40-100 aa, which is what the `protein` class covers here.
MINIBINDER_RANGE = (40, 100)


def clashing_designs():
    """Designs with at least one hard overlap against the intact receptor.

    Only the complexes are read. Step 13 also reports the binder-alone
    `_monomer` files and marks them "not scored", because no chain in those
    files reads as human EGFR so there is nothing to place the binder against.
    Treating a "not scored" monomer as a clash would remove every candidate.
    """
    if not CLASHES.exists():
        sys.exit(f"missing {CLASHES}. Run analysis/13_full_receptor_clash.py first.")
    out = set()
    for row in csv.DictReader(CLASHES.open()):
        if row["verdict"] == "clashing":
            out.add(row["design"])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-designs", type=int, default=MAX_DESIGNS)
    args = parser.parse_args()

    if not SHORTLIST.exists():
        sys.exit(f"missing {SHORTLIST}. Run analysis/10_charge_pair_filter.py first.")

    shortlist = list(csv.DictReader(SHORTLIST.open()))
    clashing = clashing_designs()

    print(f"shortlist           {len(shortlist)} candidates step 10 did not reject")
    print(f"clashing on intact  {len(clashing)} removed here and nowhere else")

    chosen, dropped = [], []
    for row in shortlist:
        if row["design"] in clashing:
            dropped.append(row)
            continue
        chosen.append(row)
        if len(chosen) >= args.max_designs:
            break

    for row in dropped:
        print(f"  removed {row['design']}: hard overlap against the intact receptor")

    if not chosen:
        sys.exit("nothing survived both screens; no submission written")

    SUBMISSION.parent.mkdir(parents=True, exist_ok=True)
    # Newlines are written as a bare \n. Python's csv module defaults to the
    # \r\n that RFC 4180 specifies, which is correct for the format and wrong
    # for this use: the upload is read by a web form, every other text file in
    # this repository is \n, and a stray carriage return inside a sequence field
    # is the kind of thing that survives a visual check and fails a parser.
    with SUBMISSION.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "sequence", "molecule_class"])
        for index, row in enumerate(chosen, start=1):
            sequence = row["Binder_Sequence"].strip()
            length = len(sequence)
            low, high = MINIBINDER_RANGE
            if not low <= length <= high:
                sys.exit(f"{row['design']} is {length} aa, outside the "
                         f"{low}-{high} aa minibinder band this script assigns "
                         f"molecule_class 'protein' for. Choose the class "
                         f"deliberately rather than letting this default stand.")
            writer.writerow([f"egfr-ph-h370-{index:02d}", sequence, "protein"])

    with PROVENANCE.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "design", "correct_pairs", "total_contact_pairs",
                         "species_difference_contact", "i_pDAE", "length"])
        for index, row in enumerate(chosen, start=1):
            contacted = row["target_positions_contacted"].split()
            writer.writerow([f"egfr-ph-h370-{index:02d}", row["design"],
                             row["correct_pairs"], row["total_pairs"],
                             "yes" if "442" in contacted else "no",
                             row["i_pDAE"], row["length"]])

    print(f"\nwrote {SUBMISSION.relative_to(ROOT)}  {len(chosen)} designs")
    print(f"wrote {PROVENANCE.relative_to(ROOT)}  the mapping back to each design")
    print("\nNow run analysis/20_submission_check.py before uploading anything.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
