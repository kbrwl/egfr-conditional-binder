#!/usr/bin/env python3
"""
00_numbering_check.py — guard against the 24-residue offset trap.

WHAT THIS CHECKS AND WHY IT MATTERS
-----------------------------------
This project refers to residues by their position in the full UniProt record
P00533 (human EGFR), which is 1210 amino acids long. In that record, residues
1-24 are the signal peptide -- a short leader sequence that is cut off and
thrown away when the protein is manufactured, so it is not present in the
mature protein that actually sits on the cell surface.

The official competition constructs are the MATURE extracellular region:
UniProt residues 25-645. Their position 1 is our position 25.

  UniProt position = challenge-construct position + 24

Nothing in a FASTA file records which convention it uses. If we ever read a
challenge-construct position as a UniProt position, every residue number in the
project shifts by 24 and every result is silently wrong -- no error, no crash,
just wrong. This script asserts the convention holds, so that failure is loud.

It verifies:
  1. Both challenge constructs have the expected length (human 621, mouse 623).
  2. The human challenge construct is an exact substring of UniProt 25-645.
  3. UniProt position 415 is threonine (T), reached three independent ways.
  4. All eight pH anchors read as expected in BOTH numbering systems.

Run standalone:  python analysis/00_numbering_check.py
Exit code 0 = convention holds. Non-zero = STOP, do not run anything else.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UNIPROT_FASTA = ROOT / "data" / "sequences" / "egfr-uniprot-full.fasta"
CHALLENGE_FASTA = ROOT / "data" / "sequences" / "egfr-challenge-constructs.fasta"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"

# The offset between the two numbering systems. Signal peptide is residues 1-24,
# so the mature protein's residue 1 is UniProt residue 25.
SIGNAL_PEPTIDE_LEN = 24

# Mature extracellular region, in UniProt numbering.
ECD_START, ECD_END = 25, 645

# The eight pH-switch anchors, in UniProt numbering, with expected identity.
# Acidic (pair with a histidine on the binder): D416 E421 E424 E455 D458 D460
# Target histidines (pair with D or E on the binder): H418 H433
ANCHORS = {
    416: "D", 418: "H", 421: "E", 424: "E",
    433: "H", 455: "E", 458: "D", 460: "D",
}


def read_fasta(path):
    """Return {header_without_'>': sequence} preserving file order."""
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
            return header, seq
    raise KeyError(f"no record whose header contains {token!r} in {list(records)}")


def main():
    failures = []
    lines = []

    def emit(text=""):
        print(text)
        lines.append(text)

    def check(label, ok, detail):
        mark = "PASS" if ok else "FAIL"
        emit(f"  [{mark}] {label}: {detail}")
        if not ok:
            failures.append(f"{label}: {detail}")
        return ok

    uniprot = read_fasta(UNIPROT_FASTA)
    challenge = read_fasta(CHALLENGE_FASTA)

    _, human_full = find(uniprot, "EGFR_HUMAN")
    _, mouse_full = find(uniprot, "EGFR_MOUSE")
    _, human_chal = find(challenge, "EGFR_HUMAN_challenge")
    _, mouse_chal = find(challenge, "EGFR_MOUSE_challenge")

    emit("=" * 72)
    emit("NUMBERING CONVENTION CHECK")
    emit("=" * 72)
    emit()
    emit("Convention: all residue numbers in this project are positions in the")
    emit("FULL UniProt record. UniProt position = challenge position + 24.")
    emit()

    emit("1. Record lengths")
    check("human UniProt P00533 full length", len(human_full) == 1210,
          f"{len(human_full)} aa (expected 1210)")
    check("mouse UniProt Q01279 full length", len(mouse_full) > 0,
          f"{len(mouse_full)} aa (no fixed expectation asserted)")
    check("human challenge construct length", len(human_chal) == 621,
          f"{len(human_chal)} aa (expected 621)")
    check("mouse challenge construct length", len(mouse_chal) == 623,
          f"{len(mouse_chal)} aa (expected 623)")
    emit()

    # Slice the UniProt record to the mature ECD. Python is 0-indexed and slice
    # ends are exclusive, so UniProt residues 25..645 are [24:645].
    human_ecd = human_full[ECD_START - 1:ECD_END]
    emit("2. Does the official human construct equal UniProt 25-645?")
    check("slice length", len(human_ecd) == 621, f"{len(human_ecd)} aa")
    identical = human_ecd == human_chal
    check("character-for-character identity", identical,
          "identical" if identical else "DIFFERS -- provenance broken")
    if not identical:
        diffs = [(i + ECD_START, a, b)
                 for i, (a, b) in enumerate(zip(human_ecd, human_chal)) if a != b]
        for pos, a, b in diffs[:10]:
            emit(f"        UniProt {pos}: full={a} challenge={b}")
    emit()

    emit("3. Is UniProt residue 415 a threonine, by three independent routes?")
    via_full = human_full[415 - 1]
    via_ecd = human_ecd[415 - ECD_START]
    via_chal = human_chal[415 - SIGNAL_PEPTIDE_LEN - 1]
    check("full UniProt record, index 414", via_full == "T", f"{via_full}")
    check("UniProt 25-645 slice", via_ecd == "T", f"{via_ecd}")
    check("challenge construct, index 390", via_chal == "T", f"{via_chal}")
    emit()

    emit("   CAVEAT: the mouse challenge construct is 623 aa, two longer than")
    emit("   mouse UniProt 25-645 (621 aa). Mouse carries a 2-residue insertion")
    emit("   near human-equivalent position 638. The +24 offset is therefore only")
    emit("   valid UPSTREAM of that insertion. Domain III (310-480) and our")
    emit("   epitope (415-466) sit well upstream, so the offset is safe here --")
    emit("   but do not reuse it for mouse positions past ~638 without rechecking.")
    emit()

    emit("4. The eight pH anchors, read in both numbering systems")
    emit("   (identity must also match between human and mouse -- these anchors")
    emit("    are the ones we claim are cross-species conserved)")
    rows = []
    for pos in sorted(ANCHORS):
        want = ANCHORS[pos]
        got_full = human_full[pos - 1]
        got_chal = human_chal[pos - SIGNAL_PEPTIDE_LEN - 1]
        got_mouse = mouse_chal[pos - SIGNAL_PEPTIDE_LEN - 1]
        ok = got_full == want and got_chal == want
        conserved = got_mouse == want
        kind = "acidic (binder gets HIS)" if want in "DE" else "target HIS (binder gets D/E)"
        check(f"anchor {want}{pos}",
              ok and conserved,
              f"UniProt={got_full} challenge={got_chal} mouse={got_mouse} | {kind}")
        rows.append((pos, pos - SIGNAL_PEPTIDE_LEN, want, got_full, got_chal,
                     got_mouse, conserved, kind))
    emit()

    DERIVED.mkdir(parents=True, exist_ok=True)
    out_csv = DERIVED / "numbering-check.csv"
    with out_csv.open("w") as fh:
        fh.write("uniprot_pos,challenge_pos,expected_aa,human_uniprot_aa,"
                 "human_challenge_aa,mouse_challenge_aa,conserved,role\n")
        for r in rows:
            fh.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]},{r[6]},\"{r[7]}\"\n")
    emit(f"Wrote {out_csv.relative_to(ROOT)}")

    emit()
    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED -- STOP.")
        for f in failures:
            emit(f"  - {f}")
    else:
        emit("RESULT: ALL CHECKS PASSED.")
        emit("The numbering convention holds. UniProt pos = challenge pos + 24.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "00-numbering-check.md").write_text(
        "# Numbering convention check\n\n"
        "Computed output of `analysis/00_numbering_check.py`. Do not hand-edit.\n\n"
        "All residue numbers in this project are positions in the full UniProt\n"
        "record (human P00533). The official challenge constructs are the mature\n"
        "extracellular region, UniProt 25-645, so:\n\n"
        "    UniProt position = challenge-construct position + 24\n\n"
        "```\n" + "\n".join(lines) + "\n```\n"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
