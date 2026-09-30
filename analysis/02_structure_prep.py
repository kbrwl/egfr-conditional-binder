#!/usr/bin/env python3
"""
02_structure_prep.py — download 6ARU, resolve its numbering, split out the receptor.

WHAT THIS DOES AND WHY IT MATTERS
---------------------------------
Steps 00 and 01 compared letters in a row. That tells us what each residue is,
but not which direction it points, and domain III folds into a solenoid -- a
spiral-staircase shape -- in which residues that sit next to each other in the
sequence can face opposite ways. To find out which of our eight anchors a binder
could actually reach, we need measured 3D coordinates.

6ARU is an entry in the Protein Data Bank (PDB), the public archive of measured
3D protein structures. It holds the EGFR extracellular region bound to a mutant
cetuximab Fab -- a Fab is the gripping arm of an antibody, cut free of the rest.
We use 6ARU rather than 1YY9 because it contains the whole extracellular region
rather than domain III alone, which also lets step 06 ask whether neighbouring
domains cover our epitope. It is the entry the competition page references too.

THE NUMBERING OFFSET
--------------------
A PDB file carries its own residue numbers, and there is no rule about which
convention they follow. Structures of secreted proteins are very often numbered
by the mature protein -- residue 1 = residue 25 in UniProt, the public protein
sequence archive -- because the mature protein is what was actually crystallised.
Nothing in the file records the choice.

Assume wrongly and every number in steps 03-06 shifts by 24 positions, with no
error raised anywhere. We would compute the solvent-accessible surface area (how
much of a residue's surface water can reach) of the wrong residues, confidently,
and design against a surface that does not exist.

So this script does not assume. It extracts the receptor chain's observed
sequence from the coordinates, aligns that against the human UniProt sequence,
and works the offset out position by position from that alignment. It then
verifies the result by checking that our eight anchors really are the amino acids
(aa) we expect. If any disagrees, the script stops.

It also identifies which chain is which, from the file's own annotation rather
than from assumption, and reports which residue ranges are actually resolved.
Crystal structures routinely have missing loops -- regions too mobile to see --
and a missing residue means no coordinates and therefore no answer for that
position.

Outputs:
  data/structures/6aru.pdb                  (downloaded, gitignored)
  data/structures/6aru_receptor_only.pdb    (Fab removed, for step 03)
  results/findings/02-structure-prep.md
  data/derived/02-chain-inventory.csv
  data/derived/02-numbering-offset.csv

Run standalone:  python analysis/02_structure_prep.py
Exit code 0 = offset resolved and anchors verified.
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

import requests
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "data" / "sequences" / "egfr-uniprot-full.fasta"
STRUCT_DIR = ROOT / "data" / "structures"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"

PDB_ID = "6ARU"
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
CIF_URL = f"https://files.rcsb.org/download/{PDB_ID.lower()}.cif"

EPI_START, EPI_END = 415, 466
ANCHORS = {416: "D", 418: "H", 421: "E", 424: "E",
           433: "H", 455: "E", 458: "D", 460: "D"}

THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}


def read_human():
    seq, keep = [], False
    for line in FASTA.read_text().splitlines():
        if line.startswith(">"):
            keep = "EGFR_HUMAN" in line
        elif keep:
            seq.append(line.strip())
    return "".join(seq)


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return f"already present ({dest.stat().st_size:,} bytes)"
    dest.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return f"downloaded ({len(resp.content):,} bytes)"


def parse_compnd(pdb_path):
    """Read the file's own COMPND annotation: molecule name -> chain IDs.

    COMPND is the header record where the depositors wrote down what each chain is,
    so taking the mapping from here means we never have to guess which chain is which.
    """
    text, current, mapping = pdb_path.read_text().splitlines(), None, {}
    raw = []
    for line in text:
        if line.startswith("COMPND"):
            raw.append(line[10:].strip())
        elif line.startswith("ATOM"):
            break
    joined = " ".join(raw)
    # COMPND is a set of "TOKEN: value;" pairs, grouped per molecule.
    mol = None
    for field in joined.split(";"):
        field = field.strip()
        if field.upper().startswith("MOLECULE:"):
            mol = field.split(":", 1)[1].strip()
        elif field.upper().startswith("CHAIN:") and mol:
            chains = [c.strip() for c in field.split(":", 1)[1].split(",")]
            for c in chains:
                mapping[c] = mol
    return mapping


class ChainSelect(Select):
    def __init__(self, keep):
        self.keep = set(keep)

    def accept_chain(self, chain):
        return chain.id in self.keep

    def accept_residue(self, residue):
        # Keep only standard amino acids: drop water, ions and glycans (attached
        # sugar chains), so the file handed to the solvent-accessible surface area
        # (SASA) step in 03 is bare protein surface.
        return is_aa(residue, standard=True)


def main():
    out, failures = [], []

    def emit(text=""):
        print(text)
        out.append(text)

    def check(label, ok, detail):
        emit(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")
        if not ok:
            failures.append(f"{label}: {detail}")

    emit("=" * 72)
    emit(f"STRUCTURE PREPARATION — PDB {PDB_ID}")
    emit("=" * 72)
    emit()

    pdb_path = STRUCT_DIR / f"{PDB_ID.lower()}.pdb"
    emit("0. Acquire structure")
    emit(f"   {PDB_URL}")
    emit(f"   {download(PDB_URL, pdb_path)} -> data/structures/{pdb_path.name}")
    try:
        cif_path = STRUCT_DIR / f"{PDB_ID.lower()}.cif"
        emit(f"   {download(CIF_URL, cif_path)} -> data/structures/{cif_path.name}"
             "  (mmCIF backup)")
    except Exception as exc:                                   # noqa: BLE001
        emit(f"   mmCIF backup unavailable ({exc}) — not required")
    emit()

    # Report what the file says about itself.
    header_lines = [ln for ln in pdb_path.read_text().splitlines()[:40]
                    if ln.startswith(("HEADER", "TITLE", "EXPDTA", "REMARK   2"))]
    emit("1. What the file says about itself")
    for ln in header_lines[:8]:
        emit(f"   {ln.rstrip()}")
    emit()

    compnd = parse_compnd(pdb_path)
    emit("   Chain annotation as recorded in the file's own COMPND records,")
    emit("   which is how the depositors described each chain:")
    for c in sorted(compnd):
        emit(f"     chain {c}: {compnd[c]}")
    emit()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, str(pdb_path))
    model = structure[0]
    if len(structure) > 1:
        emit(f"   Note: {len(structure)} models present; using the first.")
        emit()

    human = read_human()

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    # Inventory every chain: observed residues, sequence, and identity to human EGFR.
    emit("2. Chain inventory, measured from the coordinates")
    emit()
    emit("   'Observed' means residues with actual coordinates. Residues listed in")
    emit("   the file's sequence but too mobile to locate are absent here.")
    emit()
    chain_data = {}
    for chain in model:
        residues = [r for r in chain if is_aa(r, standard=True)]
        if not residues:
            continue
        seq = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in residues)
        numbers = [r.id[1] for r in residues]
        # Local alignment to score how EGFR-like this chain is.
        local = Align.PairwiseAligner()
        local.mode = "local"
        local.open_gap_score = -11
        local.extend_gap_score = -1
        local.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aln = local.align(human, seq)[0]
        matches = sum(1 for a, b in zip(aln[0], aln[1]) if a == b and a != "-")
        pct = 100.0 * matches / len(seq)
        chain_data[chain.id] = dict(residues=residues, seq=seq, numbers=numbers,
                                    egfr_identity=pct)
        emit(f"   chain {chain.id}: {len(residues):4d} observed residues, "
             f"numbered {min(numbers)}-{max(numbers)}, "
             f"{pct:5.1f}% identity to human EGFR")
        emit(f"              annotation: {compnd.get(chain.id, '(none in file)')}")
    emit()

    # The receptor is whichever chain matches human EGFR, taken from the identity
    # scores measured above rather than from the chain letter.
    ranked = sorted(chain_data.items(), key=lambda kv: -kv[1]["egfr_identity"])
    receptor_id = ranked[0][0]
    receptor_pct = ranked[0][1]["egfr_identity"]
    fab_ids = [cid for cid, d in chain_data.items() if cid != receptor_id]

    emit("3. Which chain is the receptor?")
    emit()
    emit("   Answered by measuring each chain's identity to the human EGFR sequence.")
    check("receptor identified by EGFR identity",
          receptor_pct > 90,
          f"chain {receptor_id} at {receptor_pct:.1f}% — "
          f"next best {ranked[1][0]} at {ranked[1][1]['egfr_identity']:.1f}%")
    check("competition page claim 'chain A is the receptor'",
          receptor_id == "A",
          f"computed receptor is chain {receptor_id} — "
          f"{'claim confirmed' if receptor_id == 'A' else 'CLAIM CONTRADICTED'}")
    emit()
    emit(f"   Fab chains (everything not the receptor): {', '.join(sorted(fab_ids))}")
    for cid in sorted(fab_ids):
        emit(f"     chain {cid}: {compnd.get(cid, '(none in file)')} "
             f"— {len(chain_data[cid]['residues'])} residues")
    emit()
    emit("   Which of the two Fab chains is the heavy one and which the light comes")
    emit("   from the file's COMPND annotation above. Steps 03-04 do not depend on")
    emit("   that: contacts are computed against 'either Fab chain', and both chains")
    emit("   are removed together for the receptor-only file.")
    emit()

    # ---- The offset, measured from the alignment rather than assumed ----
    emit("4. The numbering offset, measured rather than assumed")
    emit()
    emit("   Method: align the receptor chain's observed sequence against the full")
    emit("   human UniProt sequence, then for every observed residue compute")
    emit("      offset = UniProt position - PDB residue number")
    emit("   A single dominant value means the file follows one consistent")
    emit("   convention. We measure it because nothing in the file records which")
    emit("   convention that is, and guessing wrong shifts every later result.")
    emit()

    rec = chain_data[receptor_id]
    aln = aligner.align(human, rec["seq"])[0]

    # Walk the alignment, tracking index into human and into the observed chain.
    hi, si = 0, 0
    pdb_to_uniprot = {}
    offsets = Counter()
    for a, b in zip(aln[0], aln[1]):
        if a != "-":
            hi += 1
        if b != "-":
            si += 1
        if a != "-" and b != "-":
            pdb_num = rec["numbers"][si - 1]
            pdb_to_uniprot[pdb_num] = hi
            offsets[hi - pdb_num] += 1

    total = sum(offsets.values())
    emit("   | offset | aligned residues | share |")
    emit("   |---|---|---|")
    for off, n in offsets.most_common(5):
        emit(f"   | {off:+d} | {n} | {100.0 * n / total:.1f}% |")
    emit()

    dominant_offset, dominant_n = offsets.most_common(1)[0]
    share = 100.0 * dominant_n / total
    check("offset is consistent across the chain", share > 95,
          f"{share:.1f}% of aligned residues share offset {dominant_offset:+d}")

    if dominant_offset == 0:
        interpretation = ("PDB numbering already matches UniProt numbering. "
                          "No conversion needed.")
    elif dominant_offset == 24:
        interpretation = ("PDB numbers by the mature protein, so UniProt = PDB + 24. "
                          "This is the convention the brief warned about.")
    else:
        interpretation = (f"PDB uses an unexpected convention. "
                          f"UniProt = PDB + {dominant_offset}.")
    emit(f"   Interpretation: {interpretation}")
    emit()

    uniprot_to_pdb = {u: p for p, u in pdb_to_uniprot.items()}

    # ---- Verify against the eight anchors ----
    emit("5. Verification: are the eight anchors the residues we expect?")
    emit()
    emit("   This is the check that catches an offset error. If the offset were")
    emit("   wrong by 24, these positions would read as the wrong amino acids, and")
    emit("   the script stops rather than letting steps 03-06 run on them.")
    emit()
    anchor_rows = []
    for pos in sorted(ANCHORS):
        want = ANCHORS[pos]
        if pos not in uniprot_to_pdb:
            check(f"anchor {want}{pos}", False,
                  "not resolved in the structure — no coordinates for it")
            anchor_rows.append((pos, want, "", "", False, "unresolved"))
            continue
        pdb_num = uniprot_to_pdb[pos]
        res = next(r for r in rec["residues"] if r.id[1] == pdb_num)
        got = THREE_TO_ONE.get(res.get_resname(), "X")
        ok = got == want
        check(f"anchor {want}{pos}", ok,
              f"PDB {receptor_id}/{res.get_resname()}{pdb_num} reads {got}"
              f"{'' if ok else f' — expected {want}'}")
        anchor_rows.append((pos, want, pdb_num, got, ok, "resolved"))
    emit()

    # ---- Resolved ranges and gaps ----
    emit("6. Which residues are actually resolved?")
    emit()
    observed_uniprot = sorted(pdb_to_uniprot.values())
    segments = []
    start = prev = observed_uniprot[0]
    for u in observed_uniprot[1:]:
        if u != prev + 1:
            segments.append((start, prev))
            start = u
        prev = u
    segments.append((start, prev))
    emit(f"   Receptor chain {receptor_id}, in UniProt numbering:")
    emit(f"   spans {observed_uniprot[0]}-{observed_uniprot[-1]}, "
         f"{len(observed_uniprot)} residues observed in "
         f"{len(segments)} continuous segment(s)")
    for a, b in segments:
        emit(f"     {a}-{b}  ({b - a + 1} residues)")
    if len(segments) > 1:
        emit("   Unresolved gaps (mobile loops, no coordinates):")
        for (a1, b1), (a2, _) in zip(segments, segments[1:]):
            emit(f"     {b1 + 1}-{a2 - 1}  ({a2 - b1 - 1} residues missing)")
    emit()

    epi_missing = [p for p in range(EPI_START, EPI_END + 1)
                   if p not in uniprot_to_pdb]
    emit(f"   Inside our epitope {EPI_START}-{EPI_END}:")
    check("epitope fully resolved", not epi_missing,
          "all 52 residues have coordinates" if not epi_missing
          else f"{len(epi_missing)} missing: {epi_missing}")
    if epi_missing:
        emit("   A missing residue has no coordinates, so steps 03 and 05 have no")
        emit("   answer for that position. A missing anchor cannot be assessed, and")
        emit("   must not be counted as usable on the strength of the sequence alone.")
    emit()

    # ---- Write receptor-only structure ----
    emit("7. Writing receptor-only structure for step 03")
    emit()
    emit("   Step 03 measures solvent accessibility on the receptor by itself, so")
    emit("   the Fab chains come out here. The alternative, measuring on the whole")
    emit("   complex, was rejected because cetuximab sits on the very surface we")
    emit("   care about: it would report our epitope as buried when it is only")
    emit("   covered by an antibody that will not be present in our assay, and we")
    emit("   would discard a workable epitope for that reason.")
    emit()
    io = PDBIO()
    io.set_structure(structure)
    receptor_only = STRUCT_DIR / f"{PDB_ID.lower()}_receptor_only.pdb"
    io.save(str(receptor_only), select=ChainSelect([receptor_id]))
    kept = sum(1 for r in parser.get_structure("chk", str(receptor_only))[0][receptor_id]
               if is_aa(r, standard=True))
    check("receptor-only file written", kept == len(rec["residues"]),
          f"{kept} residues written to data/structures/{receptor_only.name}")
    emit()

    # ---- CSVs ----
    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "02-chain-inventory.csv").open("w") as fh:
        fh.write("chain,role,file_annotation,observed_residues,"
                 "pdb_num_min,pdb_num_max,pct_identity_to_human_egfr\n")
        for cid, d in sorted(chain_data.items()):
            role = "receptor" if cid == receptor_id else "fab"
            ann = compnd.get(cid, "").replace(",", ";")
            fh.write(f"{cid},{role},\"{ann}\",{len(d['residues'])},"
                     f"{min(d['numbers'])},{max(d['numbers'])},"
                     f"{d['egfr_identity']:.2f}\n")
    with (DERIVED / "02-numbering-offset.csv").open("w") as fh:
        fh.write("uniprot_pos,pdb_resnum,offset\n")
        for pdb_num, uni in sorted(pdb_to_uniprot.items(), key=lambda kv: kv[1]):
            fh.write(f"{uni},{pdb_num},{uni - pdb_num}\n")
    with (DERIVED / "02-anchor-verification.csv").open("w") as fh:
        fh.write("uniprot_pos,expected_aa,pdb_resnum,observed_aa,match,status\n")
        for r in anchor_rows:
            fh.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]}\n")
    emit("   Wrote data/derived/02-chain-inventory.csv")
    emit("   Wrote data/derived/02-numbering-offset.csv")
    emit("   Wrote data/derived/02-anchor-verification.csv")

    emit()
    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED — STOP, do not run step 03.")
        for f in failures:
            emit(f"  - {f}")
    else:
        emit("RESULT: PASSED.")
        emit(f"Receptor is chain {receptor_id}. "
             f"UniProt position = PDB residue number {dominant_offset:+d}.")
        emit("All eight anchors read as the expected amino acids; steps 03-06 can proceed.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "02-structure-prep.md").write_text(
        f"# Structure preparation: PDB {PDB_ID}\n\n"
        "Computed output of `analysis/02_structure_prep.py`. Do not hand-edit.\n\n"
        "Works out, by measurement rather than assumption, how the residue numbers\n"
        "in this structure file relate to ours. PDB is the Protein Data Bank, the\n"
        "public archive of measured three-dimensional structures; UniProt is the\n"
        "public sequence archive our numbering follows. Also identifies which chain\n"
        "is the receptor and which are the antibody, reports which residues the\n"
        "structure actually resolves, and writes out a receptor-only file for the\n"
        "accessibility step.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
