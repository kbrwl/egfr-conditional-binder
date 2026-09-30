#!/usr/bin/env python3
"""
04_cetuximab_contacts.py — what does cetuximab actually touch?

WHAT THIS MEASURES AND WHY IT MATTERS
-------------------------------------
Cetuximab is an approved antibody drug that binds EGFR domain III. It is known
not to bind mouse EGFR. A previous working assumption in this project was that
four specific human/mouse differences -- Q390R, E412D, R414W and K467R -- sit in
cetuximab's grip and explain that failure. That assumption was stated from memory
and never computed. It is flagged UNVERIFIED in the decisions log. This script
resolves it by measurement.

A "contact" needs a definition, because there is no natural boundary where one
protein stops touching another. The convention used here: an EGFR residue is in
contact if ANY of its heavy atoms sits within 4.5 angstroms of ANY heavy atom of
either Fab chain. "Heavy atom" means any atom except hydrogen -- hydrogens are
too light to appear in most crystal structures and are simply absent from this
file, so heavy-atom distances are what we actually have. 4.5 A is the standard
cutoff for "these two residues are touching".

Distance queries use Bio.PDB.NeighborSearch, which builds a spatial index, rather
than a brute-force double loop over every atom pair. With ~5000 protein atoms a
double loop is ~25 million comparisons; the index makes it near-instant. Same
answer, different cost.

THREE QUESTIONS THIS ANSWERS
----------------------------
(a) Do any of the 16 human/mouse domain III differences fall in the contact set,
    and specifically Q390R, E412D, R414W, K467R? If yes, that is a computed
    explanation for cetuximab's species failure. If no, the prior speculation is
    DISPROVED and must be recorded as such.

(b) How much does cetuximab's footprint overlap our 415-466 block? We want to be
    NEARBY -- close enough that our binder blocks the same functional site -- but
    NOT on identical residues, because the competition rules forbid starting from
    an existing binder, and landing on exactly cetuximab's footprint invites the
    novelty objection even for a genuinely de novo design.

(c) Which residues in 415-466 does cetuximab NOT touch? Those are where our
    design can differentiate itself.

Outputs:
  results/findings/04-cetuximab-contacts.md
  data/derived/04-cetuximab-contacts.csv
  data/derived/04-epitope-overlap.csv

Run standalone:  python analysis/04_cetuximab_contacts.py
"""

import sys
from pathlib import Path

from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

ROOT = Path(__file__).resolve().parents[1]
COMPLEX_PDB = ROOT / "data" / "structures" / "6aru.pdb"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"
OFFSET_CSV = DERIVED / "02-numbering-offset.csv"
DIFFS_CSV = DERIVED / "01-domain3-differences.csv"

CUTOFF = 4.5          # angstroms, heavy atom to heavy atom
RECEPTOR_CHAIN = "A"  # established by measurement in step 02
EPI_START, EPI_END = 415, 466
ANCHORS = {416: "D", 418: "H", 421: "E", 424: "E",
           433: "H", 455: "E", 458: "D", 460: "D"}

# The four differences the earlier speculation named.
SPECULATED = ["Q390R", "E412D", "R414W", "K467R"]

THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}


def load_offset():
    if not OFFSET_CSV.exists():
        raise SystemExit("Run analysis/02_structure_prep.py first.")
    pdb_to_uni = {}
    for line in OFFSET_CSV.read_text().splitlines()[1:]:
        uni, pdb, _ = line.split(",")
        pdb_to_uni[int(pdb)] = int(uni)
    return pdb_to_uni


def load_differences():
    if not DIFFS_CSV.exists():
        raise SystemExit("Run analysis/01_alignment.py first.")
    diffs = []
    for line in DIFFS_CSV.read_text().splitlines()[1:]:
        pos, _chal, h, m, label, _inside = line.split(",")
        diffs.append((int(pos), h, m, label))
    return diffs


def main():
    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    pdb_to_uni = load_offset()
    differences = load_differences()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(COMPLEX_PDB))
    model = structure[0]

    receptor = model[RECEPTOR_CHAIN]
    fab_chains = [c for c in model
                  if c.id != RECEPTOR_CHAIN
                  and any(is_aa(r, standard=True) for r in c)]

    emit("=" * 72)
    emit("CETUXIMAB CONTACT SET — resolving an UNVERIFIED claim")
    emit("=" * 72)
    emit()
    emit(f"Structure: 6ARU, receptor chain {RECEPTOR_CHAIN}, "
         f"Fab chains {', '.join(c.id for c in fab_chains)}")
    emit(f"Definition: an EGFR residue is in contact if any heavy atom is within")
    emit(f"            {CUTOFF} A of any heavy atom of either Fab chain.")
    emit("Method: Bio.PDB.NeighborSearch spatial index (not a double loop).")
    emit("Numbering: UniProt = PDB + 24, established empirically by step 02.")
    emit()

    # Build the search index over Fab heavy atoms.
    fab_atoms = [a for c in fab_chains for r in c if is_aa(r, standard=True)
                 for a in r if a.element != "H"]
    emit(f"Fab heavy atoms indexed: {len(fab_atoms)}")
    ns = NeighborSearch(fab_atoms)

    # For each receptor residue, find the closest Fab atom within the cutoff.
    contacts = {}
    for res in receptor:
        if not is_aa(res, standard=True):
            continue
        uni = pdb_to_uni.get(res.id[1])
        if uni is None:
            continue
        best = None
        for atom in res:
            if atom.element == "H":
                continue
            for near in ns.search(atom.coord, CUTOFF):
                d = atom - near
                if best is None or d < best[0]:
                    best = (d, near.get_parent().get_parent().id,
                            atom.get_id(), near.get_id())
        if best is not None:
            contacts[uni] = dict(
                aa=THREE_TO_ONE.get(res.get_resname(), "X"),
                pdb_num=res.id[1], min_dist=best[0], fab_chain=best[1],
                egfr_atom=best[2], fab_atom=best[3])

    emit(f"EGFR residues in contact with the Fab: {len(contacts)}")
    emit()

    # ---- Full contact list ----
    emit("1. Full contact list")
    emit()
    emit("   | UniProt | aa | PDB# | min dist A | Fab chain | EGFR atom | Fab atom |")
    emit("   |---|---|---|---|---|---|---|")
    for uni in sorted(contacts):
        c = contacts[uni]
        emit(f"   | {uni} | {c['aa']} | {c['pdb_num']} | {c['min_dist']:.2f} | "
             f"{c['fab_chain']} | {c['egfr_atom']} | {c['fab_atom']} |")
    emit()
    lo, hi = min(contacts), max(contacts)
    emit(f"   Footprint spans UniProt {lo}-{hi}.")
    by_chain = {}
    for c in contacts.values():
        by_chain[c["fab_chain"]] = by_chain.get(c["fab_chain"], 0) + 1
    emit("   Residues by nearest Fab chain: "
         + ", ".join(f"chain {k}: {v}" for k, v in sorted(by_chain.items())))
    emit()

    # ---- (a) species differences in the contact set ----
    emit("2. QUESTION (a): do human/mouse differences fall in the contact set?")
    emit()
    contact_set = set(contacts)
    in_contact, not_in_contact = [], []
    for pos, h, m, label in differences:
        (in_contact if pos in contact_set else not_in_contact).append(
            (pos, h, m, label))

    emit(f"   Of the 16 domain III differences, {len(in_contact)} are in contact:")
    if in_contact:
        for pos, h, m, label in in_contact:
            c = contacts[pos]
            emit(f"     {label}  —  {c['min_dist']:.2f} A from Fab chain "
                 f"{c['fab_chain']}")
    else:
        emit("     none")
    emit()
    emit(f"   The other {len(not_in_contact)} are not: "
         + ", ".join(lbl for _, _, _, lbl in not_in_contact))
    emit()

    emit("   Testing the four specifically speculated differences:")
    verdicts = {}
    for label in SPECULATED:
        pos = int(label[1:-1])
        hit = pos in contact_set
        verdicts[label] = hit
        if hit:
            c = contacts[pos]
            emit(f"     {label}: IN CONTACT — {c['min_dist']:.2f} A from chain "
                 f"{c['fab_chain']}")
        else:
            emit(f"     {label}: NOT in contact"
                 + (f" (nearest contact residue is "
                    f"{min(contact_set, key=lambda p: abs(p - pos))})"
                    if contact_set else ""))
    emit()

    n_confirmed = sum(verdicts.values())
    emit("   VERDICT ON THE PRIOR SPECULATION:")
    if n_confirmed == len(SPECULATED):
        emit("   CONFIRMED. All four speculated differences are in cetuximab's")
        emit("   contact set. This is now a computed explanation for why")
        emit("   cetuximab fails on mouse EGFR, and replaces the speculation.")
    elif n_confirmed == 0:
        emit("   DISPROVED. None of the four speculated differences is in")
        emit("   cetuximab's contact set. The memory-based explanation was wrong")
        emit("   and must be recorded as disproved, not quietly dropped.")
    else:
        emit(f"   PARTIALLY CONFIRMED — {n_confirmed} of {len(SPECULATED)} are in")
        emit("   contact. The speculation was directionally right but not")
        emit("   accurate as stated. It must be restated to the computed set")
        emit("   rather than kept as originally written.")
        emit()
        emit("   Confirmed in contact: "
             + ", ".join(k for k, v in verdicts.items() if v))
        emit("   Not in contact: "
             + ", ".join(k for k, v in verdicts.items() if not v))
    emit()
    if in_contact:
        emit("   Note what this does and does not establish. Showing that a")
        emit("   differing residue is inside the footprint is consistent with it")
        emit("   causing the species failure, and is far better evidence than an")
        emit("   assumption. It is not proof: whether a given substitution")
        emit("   actually abolishes binding depends on how much that contact")
        emit("   contributes, which this calculation does not measure.")
    emit()

    # ---- (b) overlap with our block ----
    emit(f"3. QUESTION (b): overlap between cetuximab's footprint and {EPI_START}-{EPI_END}")
    emit()
    epi_contacts = sorted(p for p in contact_set if EPI_START <= p <= EPI_END)
    epi_size = EPI_END - EPI_START + 1
    emit(f"   Cetuximab touches {len(epi_contacts)} of the {epi_size} residues in")
    emit(f"   our block ({100.0 * len(epi_contacts) / epi_size:.1f}%).")
    emit(f"   Of cetuximab's {len(contact_set)} total contacts, "
         f"{len(epi_contacts)} "
         f"({100.0 * len(epi_contacts) / len(contact_set):.1f}%) are in our block.")
    emit()
    if epi_contacts:
        emit("   Which ones:")
        for p in epi_contacts:
            c = contacts[p]
            mark = "  <-- ANCHOR" if p in ANCHORS else ""
            emit(f"     {c['aa']}{p}: {c['min_dist']:.2f} A from chain "
                 f"{c['fab_chain']}{mark}")
    else:
        emit("   None.")
    emit()

    anchors_in_contact = [p for p in sorted(ANCHORS) if p in contact_set]
    emit(f"   Anchors inside cetuximab's footprint: "
         f"{len(anchors_in_contact)} of 8"
         + (f" — {', '.join(f'{ANCHORS[p]}{p}' for p in anchors_in_contact)}"
            if anchors_in_contact else ""))
    emit()
    emit("   Reading this for the design:")
    frac = len(epi_contacts) / epi_size
    if frac == 0:
        emit("   Zero overlap. Our block is a genuinely distinct surface. Good for")
        emit("   the novelty requirement — but raises the question of whether we")
        emit("   are still blocking the functionally important site at all.")
    elif frac < 0.35:
        emit("   Partial overlap. This is the position we wanted: adjacent to and")
        emit("   partly sharing the druggable surface, so a binder here plausibly")
        emit("   blocks the same function, while most of our contact positions are")
        emit("   cetuximab-independent. Novelty is defensible on residue identity.")
    else:
        emit("   Heavy overlap. Our block largely reproduces cetuximab's footprint,")
        emit("   which invites the novelty objection even for a de novo design.")
        emit("   Consider shifting the design's contact centre toward the")
        emit("   non-overlapping residues identified below.")
    emit()

    # ---- (c) what cetuximab does NOT touch ----
    emit(f"4. QUESTION (c): residues in {EPI_START}-{EPI_END} cetuximab does NOT touch")
    emit()
    untouched = [p for p in range(EPI_START, EPI_END + 1) if p not in contact_set]
    emit(f"   {len(untouched)} of {epi_size} residues are cetuximab-free.")
    emit("   These are where our design can differentiate.")
    emit()
    for i in range(0, len(untouched), 13):
        emit("   " + " ".join(str(p) for p in untouched[i:i + 13]))
    emit()
    free_anchors = [p for p in sorted(ANCHORS) if p in untouched]
    emit(f"   Anchors cetuximab does NOT touch: {len(free_anchors)} of 8"
         + (f" — {', '.join(f'{ANCHORS[p]}{p}' for p in free_anchors)}"
            if free_anchors else ""))
    emit()
    emit("   This is the most directly useful output of this step: an anchor that")
    emit("   is both reachable (step 03) and outside cetuximab's footprint gives")
    emit("   pH-switch capability on a surface no approved drug occupies.")
    emit()

    # ---- CSVs ----
    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "04-cetuximab-contacts.csv").open("w") as fh:
        fh.write("uniprot_pos,pdb_resnum,aa,min_dist_a,fab_chain,egfr_atom,"
                 "fab_atom,in_epitope,is_anchor,is_species_difference\n")
        diff_pos = {p for p, _, _, _ in differences}
        for uni in sorted(contacts):
            c = contacts[uni]
            fh.write(f"{uni},{c['pdb_num']},{c['aa']},{c['min_dist']:.3f},"
                     f"{c['fab_chain']},{c['egfr_atom']},{c['fab_atom']},"
                     f"{EPI_START <= uni <= EPI_END},{uni in ANCHORS},"
                     f"{uni in diff_pos}\n")
    with (DERIVED / "04-epitope-overlap.csv").open("w") as fh:
        fh.write("uniprot_pos,is_anchor,cetuximab_contact,min_dist_a\n")
        for p in range(EPI_START, EPI_END + 1):
            d = f"{contacts[p]['min_dist']:.3f}" if p in contacts else ""
            fh.write(f"{p},{p in ANCHORS},{p in contact_set},{d}\n")
    emit("Wrote data/derived/04-cetuximab-contacts.csv")
    emit("Wrote data/derived/04-epitope-overlap.csv")

    emit()
    emit("=" * 72)
    emit(f"RESULT: {len(contact_set)} EGFR residues contact the Fab. "
         f"{len(in_contact)} of 16 species differences are among them.")
    emit(f"Prior speculation: {n_confirmed} of 4 confirmed.")
    emit(f"Overlap with our block: {len(epi_contacts)}/{epi_size} residues.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "04-cetuximab-contacts.md").write_text(
        "# Cetuximab contact set (computed)\n\n"
        "Computed output of `analysis/04_cetuximab_contacts.py`. Do not hand-edit.\n\n"
        "Resolves the UNVERIFIED claim that Q390R, E412D, R414W and K467R explain\n"
        "cetuximab's failure on mouse EGFR. Contact definition: any EGFR heavy atom\n"
        f"within {CUTOFF} A of any heavy atom of either Fab chain, in PDB 6ARU.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
