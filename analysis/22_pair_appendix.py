"""22_pair_appendix.py -- write Appendix A of the methods document.

The methods document's central claim is that for any submitted design we can say
which binder residue faces which target residue and what each does at the two pH
values the assay uses. Asserting that is cheap. This writes the table that shows
it, one block per submitted design, generated from the screen's own output so it
cannot drift from what was actually measured.

Abbreviations, expanded because this file gets read on its own:
  EGFR  epidermal growth factor receptor, the protein we design against.
  aa    amino acids, the unit a protein's length is counted in.
  His   histidine, the one amino acid whose charge changes between pH 7.4 and
        pH 6.5: neutral above, positive below, because the pH at which half its
        copies carry the extra charge sits around 6.0 to 6.5.

WHAT EACH ROW SAYS, AND WHY THE PH COLUMNS ARE WRITTEN THE WAY THEY ARE
------------------------------------------------------------------------
A correct pair is one of two arrangements, and they behave differently enough
that the table states each separately rather than collapsing them into one
phrase:

  target acidic (D or E), binder His    at pH 7.4 the binder's His is neutral
                                        and there is no attraction; at pH 6.5 it
                                        turns positive and is drawn to the
                                        target's permanent negative charge
  target His, binder acidic (D or E)    at pH 7.4 the target's His is neutral
                                        and there is no attraction; at pH 6.5 it
                                        turns positive and is drawn to the
                                        binder's permanent negative charge

Both switch on as pH falls, which is why they reinforce rather than cancel. The
rejected arrangement -- His on the binder facing His on the target -- switches
off as pH falls, because both turn positive and repel.

Correct pairs are listed in full, one row each. Neutral contacts are summarised
as a count rather than listed, because a design makes tens of them and listing
them would bury the pairs that do the work. The count is still given, because it
is the quantity the ranking's third term uses and a reader should be able to see
how much interface sits around the pairs.

ON THE WORD "CONTACTS"
----------------------
The per-design summary line reports contacts at the 4.5 angstrom heavy-atom
cutoff: every residue-residue contact the design makes with the target face,
charged or not. An earlier draft of the methods document called these "charged
contacts", which was wrong -- for the leading design, 2 of its 31 contacts are
charge pairs and 29 are neutral. The label used here is the corrected one.

WHAT THIS TABLE IS NOT
----------------------
Every distance here is measured on a predicted structure, and every pH statement
is an inference from which amino acids face each other, not a measurement of
binding at either pH. A pair listed as correct is a pair built the right way
round; whether the switch it implies survives in a real assay is exactly what the
competition is testing and is not settled by anything in this file.
"""
import argparse
import collections
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAIRS = ROOT / "data" / "derived" / "10-candidate-pairs.csv"
PROVENANCE = ROOT / "submission" / "submission-provenance.csv"
APPENDIX = ROOT / "submission" / "appendix-a-pair-tables.md"

THREE = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly",
         "H": "His", "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn",
         "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val",
         "W": "Trp", "Y": "Tyr"}

# The two correct arrangements, and what each does at the two pH values.
BEHAVIOUR = {
    "target_acidic": ("binder His turns positive, drawn to the target's "
                      "negative charge",
                      "binder His neutral, no attraction"),
    "target_his": ("target His turns positive, drawn to the binder's "
                   "negative charge",
                   "target His neutral, no attraction"),
}


def arrangement(target_aa, binder_aa):
    if target_aa in ("D", "E") and binder_aa == "H":
        return "target_acidic"
    if target_aa == "H" and binder_aa in ("D", "E"):
        return "target_his"
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=APPENDIX)
    args = parser.parse_args()

    for path in (PAIRS, PROVENANCE):
        if not path.exists():
            sys.exit(f"missing {path}")

    by_design = collections.defaultdict(list)
    for row in csv.DictReader(PAIRS.open()):
        by_design[row["design"]].append(row)

    submitted = list(csv.DictReader(PROVENANCE.open()))

    out = []
    out.append("## Appendix A — the pair table, per submitted design\n")
    out.append(
        "For every submitted design, which binder residue faces which target\n"
        "residue, and what each pair does at the two pH values the assay uses.\n"
        "Generated by `analysis/22_pair_appendix.py` from\n"
        "`data/derived/10-candidate-pairs.csv`, so it cannot drift from what was\n"
        "measured.\n")
    out.append(
        "Histidine is the only amino acid whose charge changes across this\n"
        "range: neutral at pH 7.4, positive at pH 6.5. Aspartate and glutamate\n"
        "carry a negative charge at both. A pair is correct when it is built so\n"
        "that attraction appears as the pH falls, which happens in two ways —\n"
        "a binder histidine facing a target acidic residue, or a binder acidic\n"
        "residue facing a target histidine. Both are listed in full below.\n")
    out.append(
        "Contacts are counted at a 4.5 Å heavy-atom cutoff and include every\n"
        "residue-residue contact with the target face, charged or not. Neutral\n"
        "contacts are summarised as a count rather than listed one by one.\n")
    out.append(
        "**Every distance here is measured on a predicted structure, and every\n"
        "pH statement is an inference from which residues face each other rather\n"
        "than a measurement of binding.**\n")

    for entry in submitted:
        rows = by_design.get(entry["design"], [])
        correct = [r for r in rows
                   if r["classification"] == "correct"
                   and r["charge_atoms_resolved"] == "True"]
        neutral = [r for r in rows if r["classification"] == "neutral"]
        out.append(f"\n### {entry['name']}\n")
        if correct:
            out.append("| binder residue | target residue | verdict | "
                       "distance (Å) | behaviour at pH 6.5 | at pH 7.4 |")
            out.append("|---|---|---|---|---|---|")
            for r in sorted(correct, key=lambda r: int(r["target_pos"])):
                kind = arrangement(r["target_aa"], r["binder_aa"])
                low, high = BEHAVIOUR[kind]
                dist = r["charge_group_dist"] or r["min_dist"]
                out.append(
                    f"| {r['binder_aa']}{r['binder_pos']} "
                    f"({THREE.get(r['binder_aa'], r['binder_aa'])}) "
                    f"| {r['target_aa']}{r['target_pos']} "
                    f"({THREE.get(r['target_aa'], r['target_aa'])}) "
                    f"| correct pair | {float(dist):.1f} | {low} | {high} |")
        else:
            out.append("No correct charge pair. This design is submitted on the "
                       "interface it makes rather than on the switch, and is "
                       "ranked accordingly.")
        flags = []
        if entry["species_difference_contact"] == "yes":
            flags.append("contacts 442")
        summary = (f"\ncorrect pairs: {len(correct)} · "
                   f"contacts at 4.5 Å: {entry['total_contact_pairs']} "
                   f"(of which {len(neutral)} neutral) · "
                   f"length: {entry['length']} aa")
        if flags:
            summary += " · flags: " + ", ".join(flags)
        out.append(summary)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out) + "\n")
    total = sum(1 for e in submitted)
    print(f"wrote {args.out.relative_to(ROOT)} for {total} submitted designs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
