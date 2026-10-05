"""25_submission_note.py -- the note that goes up with the submission.

A short record of what is being uploaded, when, what it replaces, and the
caveats a reader of the submission should have in front of them. Generated from
the submission files themselves rather than written by hand, so it cannot
describe a set that is not the set being uploaded.

Abbreviations, expanded because this file gets read on its own:
  EGFR  epidermal growth factor receptor, the protein we design against.
  CSV   comma-separated values, the plain-text table the organisers want.
  aa    amino acids, the unit a protein's length is counted in.

The organisers allow one submission every 24 hours, with the most recent
designated unless another is chosen, so this note records which upload it
replaces. That cadence was confirmed by an organiser on 4 October 2026 and
replaced a standing assumption that we had exactly one attempt.
"""
import argparse
import csv
import datetime as dt
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUBMISSION = ROOT / "submission" / "challenge1-egfr-submission.csv"
PROVENANCE = ROOT / "submission" / "submission-provenance.csv"
NOVELTY = ROOT / "data" / "derived" / "24-structural-novelty.csv"
NOTE = ROOT / "submission" / "SUBMISSION-NOTE.md"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--replaces", default="none — this is the first upload")
    ap.add_argument("--billed-usd", type=float, required=True,
                    help="total billed compute for the whole project, from "
                         "the provider's own billing report")
    ap.add_argument("--date", default=dt.date.today().isoformat())
    args = ap.parse_args()

    for p in (SUBMISSION, PROVENANCE):
        if not p.exists():
            sys.exit(f"missing {p}")

    rows = list(csv.DictReader(SUBMISSION.open()))
    prov = {r["name"]: r for r in csv.DictReader(PROVENANCE.open())}
    nov = {}
    if NOVELTY.exists():
        nov = {r["design"]: r for r in csv.DictReader(NOVELTY.open())}

    backbones = {r.get("design", "").rsplit("_candidate", 1)[0]
                 for r in prov.values()}
    lengths = sorted({len(r["sequence"]) for r in rows})
    pairs = [int(prov[r["name"]]["correct_pairs"]) for r in rows if r["name"] in prov]

    out = []
    out.append("# Submission note — Challenge 1, EGFR\n")
    out.append(f"**Date:** {args.date}  ")
    out.append(f"**Replaces:** {args.replaces}  ")
    out.append(f"**Designs:** {len(rows)}  ")
    out.append(f"**Distinct backbones:** {len(backbones)}  ")
    out.append(f"**Lengths:** "
               + (f"all {lengths[0]} aa" if len(lengths) == 1
                  else f"{min(lengths)}–{max(lengths)} aa") + "  ")
    out.append(f"**Total billed compute for the whole project:** "
               f"USD {args.billed_usd:.2f}\n")

    out.append("## Files uploaded\n")
    out.append("| file | what it is |")
    out.append("|---|---|")
    out.append("| `challenge1-egfr-submission.csv` | the ranked designs, "
               "best first |")
    out.append("| `methods.md` | the written methods, with Appendix A |\n")

    out.append("## The designs\n")
    out.append("| rank | name | correct pairs | contacts (4.5 Å) | novelty | "
               "flags |")
    out.append("|---|---|---|---|---|---|")
    for i, r in enumerate(rows, start=1):
        p = prov.get(r["name"], {})
        n = nov.get(r["name"], {})
        flags = []
        if p.get("species_difference_contact") == "yes":
            flags.append("contacts 442")
        out.append(f"| {i} | {r['name']} | {p.get('correct_pairs','?')} | "
                   f"{p.get('total_contact_pairs','?')} | "
                   f"level {n.get('novelty_level','not scored')} | "
                   f"{', '.join(flags) or '—'} |")
    out.append("")

    out.append("## What a reader should know before weighing this\n")
    out.append(f"- **No design reaches the three correct charge pairs our own "
               f"rule asks for.** The best reach {max(pairs) if pairs else 0}. "
               "The rule was set as a floor for a switch that is partial rather "
               "than complete at pH 6.5, and this submission does not meet it.")
    out.append(f"- **All {len(rows)} descend from {len(backbones)} backbone"
               f"{'s' if len(backbones) != 1 else ''}.** Pairwise sequence "
               "identity within the set is reported by "
               "`analysis/20_submission_check.py`. A set this alike is "
               "variants of one idea, not independent attempts.")
    out.append("- **The pH switch is predicted, never measured.** Every pair "
               "in Appendix A is an inference from which residues face each "
               "other in a predicted structure.")
    out.append("- **Novelty is our own reading of a published scale.** The "
               "platform scores it automatically on upload and that score is "
               "the one that counts.")
    out.append("")
    out.append("## Before uploading\n")
    out.append("Run `analysis/20_submission_check.py`. It must pass.")
    out.append("")

    NOTE.write_text("\n".join(out) + "\n")
    print(f"wrote {NOTE.relative_to(ROOT)}")
    print(f"  {len(rows)} designs, {len(backbones)} backbone(s), "
          f"billed USD {args.billed_usd:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
