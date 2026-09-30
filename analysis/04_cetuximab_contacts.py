#!/usr/bin/env python3
"""
04_cetuximab_contacts.py — what does cetuximab actually touch?

WHAT THIS MEASURES AND WHY IT MATTERS
-------------------------------------
Cetuximab is an approved antibody drug that binds EGFR domain III. It is known
not to bind mouse EGFR. An earlier working assumption in this project was that
four specific human/mouse differences -- Q390R, E412D, R414W and K467R -- sit in
cetuximab's grip and explain that failure. That assumption came from memory and
was never measured, so it is flagged UNVERIFIED in the decisions log. This script
settles it by measurement.

The measurement is made on structure 6ARU from the Protein Data Bank (PDB), the
public archive of measured 3D protein structures. Residue numbers here are
positions in the human record in UniProt, the public protein sequence archive,
which is the numbering used everywhere in this project.

A "contact" has to be defined, because there is no natural boundary where one
protein stops touching another. The convention used here: an EGFR residue is in
contact if any of its heavy atoms sits within 4.5 angstroms of any heavy atom of
either Fab chain. A Fab is the gripping arm of an antibody, separated from the
rest of the antibody; the 6ARU file holds two Fab chains alongside the EGFR
chain. "Heavy atom" means any atom except hydrogen -- hydrogens are too light to
appear in most crystal structures and are simply absent from this file, so
heavy-atom distances are what we actually have to work with. 4.5 A is the
standard cutoff for "these two residues are touching".

Distance queries use Bio.PDB.NeighborSearch, which builds a spatial index, rather
than a brute-force double loop over every atom pair. With ~5000 protein atoms a
double loop is ~25 million comparisons; the index makes it near-instant. Both
routes give the same answer, so this is purely a question of cost.

THREE QUESTIONS THIS ANSWERS
----------------------------
(a) Do any of the 16 human/mouse domain III differences fall in the contact set,
    and in particular Q390R, E412D, R414W, K467R? If yes, that is a measured
    explanation for cetuximab's species failure. If no, the earlier speculation
    is disproved and gets recorded in the decisions log as disproved.

(b) How much does cetuximab's footprint overlap our 415-466 block? We want to be
    nearby -- close enough that our binder blocks the same functional site -- but
    off the identical residues, because the competition rules forbid starting
    from an existing binder, and landing on exactly cetuximab's footprint invites
    the novelty objection even for a genuinely de novo design.

(c) Which residues in 415-466 does cetuximab leave untouched? Those are where our
    design can differentiate itself.

Outputs:
  results/findings/04-cetuximab-contacts.md
  data/derived/04-cetuximab-contacts.csv
  data/derived/04-epitope-overlap.csv

Run standalone:  python analysis/04_cetuximab_contacts.py
"""

import sys
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

import egfr_common as common

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

# The four differences the earlier, memory-based assumption named.
SPECULATED = ["Q390R", "E412D", "R414W", "K467R"]

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

    numbering = common.load_numbering()
    differences = load_differences()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(COMPLEX_PDB))
    model = structure[0]

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
    emit(f"            {CUTOFF} A of any heavy atom of either Fab chain. A Fab is")
    emit("            the gripping arm of an antibody, on its own; a heavy atom is")
    emit("            any atom except hydrogen, and hydrogens are absent from this")
    emit("            structure file.")
    emit("Method: Bio.PDB.NeighborSearch builds a spatial index, so we avoid")
    emit("        comparing all ~25 million atom pairs one by one. Same answer.")
    emit("Numbering: UniProt position = PDB position + 24, measured by step 02.")
    emit()

    # The calculation lives in egfr_common.contacts_to_partner so that step 06,
    # which asks the same question, gets the same answer. Step 06 used to keep its
    # own copy of this code and reported a different contact set.
    contacts = common.contacts_to_partner(
        model, RECEPTOR_CHAIN, [c.id for c in fab_chains], numbering,
        cutoff=CUTOFF)

    emit(f"EGFR residues in contact with the Fab: {len(contacts)}")
    emit()

    # ---- Full contact list ----
    emit("1. Full contact list")
    emit()
    emit("   One row per EGFR residue cetuximab touches. 'aa' is the amino acid in")
    emit("   one-letter code, and 'min dist' is the closest heavy-atom approach")
    emit("   between that residue and the Fab.")
    emit()
    emit("   | UniProt | aa | PDB# | min dist A | Fab chain | EGFR atom | Fab atom |")
    emit("   |---|---|---|---|---|---|---|")
    for uni in sorted(contacts):
        c = contacts[uni]
        emit(f"   | {uni} | {c['aa']} | {c['pdb_resnum']} | {c['min_dist']:.2f} | "
             f"{c['partner_chain']} | {c['receptor_atom']} | {c['partner_atom']} |")
    emit()
    lo, hi = min(contacts), max(contacts)
    emit(f"   Footprint spans UniProt {lo}-{hi}.")
    by_chain = {}
    for c in contacts.values():
        by_chain[c["partner_chain"]] = by_chain.get(c["partner_chain"], 0) + 1
    emit("   Residues by nearest Fab chain: "
         + ", ".join(f"chain {k}: {v}" for k, v in sorted(by_chain.items())))
    emit()

    # ---- (a) species differences in the contact set ----
    emit("2. Question (a): do human/mouse differences fall in the contact set?")
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
                 f"{c['partner_chain']}")
    else:
        emit("     none")
    emit()
    emit(f"   The other {len(not_in_contact)} are not: "
         + ", ".join(lbl for _, _, _, lbl in not_in_contact))
    emit()

    emit("   Testing the four differences the earlier assumption named:")
    verdicts = {}
    for label in SPECULATED:
        pos = int(label[1:-1])
        hit = pos in contact_set
        verdicts[label] = hit
        if hit:
            c = contacts[pos]
            emit(f"     {label}: in contact — {c['min_dist']:.2f} A from chain "
                 f"{c['partner_chain']}")
        else:
            emit(f"     {label}: not in contact"
                 + (f" (nearest contact residue is "
                    f"{min(contact_set, key=lambda p: abs(p - pos))})"
                    if contact_set else ""))
    emit()

    n_confirmed = sum(verdicts.values())
    emit("   Verdict on the earlier assumption:")
    if n_confirmed == len(SPECULATED):
        emit("   Confirmed. All four of the named differences are in cetuximab's")
        emit("   contact set. That gives us a measured explanation for why")
        emit("   cetuximab fails on mouse EGFR, and it replaces the assumption.")
    elif n_confirmed == 0:
        emit("   Disproved. None of the four named differences is in cetuximab's")
        emit("   contact set. The explanation we were carrying from memory was")
        emit("   wrong, and the decisions log records it as disproved.")
    else:
        emit(f"   Partly confirmed: {n_confirmed} of {len(SPECULATED)} are in")
        emit("   contact. The assumption pointed in the right direction but was")
        emit("   inaccurate as written, so the decisions log now carries the")
        emit("   measured set in its place.")
        emit()
        emit("   Confirmed in contact: "
             + ", ".join(k for k, v in verdicts.items() if v))
        emit("   Not in contact: "
             + ", ".join(k for k, v in verdicts.items() if not v))
    emit()
    if in_contact:
        emit("   What this establishes, and what it does not. A differing residue")
        emit("   sitting inside the footprint is consistent with it causing the")
        emit("   species failure, and is far better evidence than an assumption.")
        emit("   It is not proof: whether a given substitution actually abolishes")
        emit("   binding depends on how much that particular contact contributes")
        emit("   to the grip, which this calculation does not measure. So we can")
        emit("   cite it as the likely explanation, but not as a settled one.")
    emit()

    # ---- (b) overlap with our block ----
    emit(f"3. Question (b): overlap between cetuximab's footprint and {EPI_START}-{EPI_END}")
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
                 f"{c['partner_chain']}{mark}")
    else:
        emit("   None.")
    emit()

    anchors_in_contact = [p for p in sorted(ANCHORS) if p in contact_set]
    emit(f"   Anchors inside cetuximab's footprint: "
         f"{len(anchors_in_contact)} of 8"
         + (f" — {', '.join(f'{ANCHORS[p]}{p}' for p in anchors_in_contact)}"
            if anchors_in_contact else ""))
    emit()
    emit("   What this means for the design:")
    frac = len(epi_contacts) / epi_size
    if frac == 0:
        emit("   Zero overlap. Our block is a genuinely separate surface, which")
        emit("   satisfies the novelty requirement easily — but it leaves open")
        emit("   whether we are still blocking the functionally important site.")
    elif frac < 0.35:
        emit("   Partial overlap. Our block sits alongside the druggable surface and")
        emit("   shares part of it, so a binder here plausibly blocks the same")
        emit("   function, while most of our contact positions are ones cetuximab")
        emit("   never touches. That makes the novelty claim defensible on residue")
        emit("   identity alone, and we can keep the block as it stands.")
    else:
        emit("   Heavy overlap. Our block largely reproduces cetuximab's footprint,")
        emit("   which invites the novelty objection even for a de novo design. The")
        emit("   response is to shift the design's contact centre toward the")
        emit("   non-overlapping residues listed below.")
    emit()

    # ---- (c) what cetuximab leaves untouched ----
    emit(f"4. Question (c): residues in {EPI_START}-{EPI_END} cetuximab leaves untouched")
    emit()
    untouched = [p for p in range(EPI_START, EPI_END + 1) if p not in contact_set]
    emit(f"   {len(untouched)} of {epi_size} residues are cetuximab-free.")
    emit("   These are the positions our design can build on to look unlike")
    emit("   cetuximab while still sitting on the same surface.")
    emit()
    for i in range(0, len(untouched), 13):
        emit("   " + " ".join(str(p) for p in untouched[i:i + 13]))
    emit()
    free_anchors = [p for p in sorted(ANCHORS) if p in untouched]
    emit(f"   Anchors cetuximab leaves untouched: {len(free_anchors)} of 8"
         + (f" — {', '.join(f'{ANCHORS[p]}{p}' for p in free_anchors)}"
            if free_anchors else ""))
    emit()
    emit("   These free anchors are the most directly useful output of this step. An")
    emit("   anchor that is both reachable (step 03) and outside cetuximab's")
    emit("   footprint lets us build the pH switch on a surface no approved drug")
    emit("   occupies, so those anchors get first preference in the design.")
    emit()

    # ---- CSVs ----
    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "04-cetuximab-contacts.csv").open("w") as fh:
        fh.write("uniprot_pos,pdb_resnum,aa,min_dist_a,fab_chain,egfr_atom,"
                 "fab_atom,in_epitope,is_anchor,is_species_difference\n")
        diff_pos = {p for p, _, _, _ in differences}
        for uni in sorted(contacts):
            c = contacts[uni]
            fh.write(f"{uni},{c['pdb_resnum']},{c['aa']},{c['min_dist']:.3f},"
                     f"{c['partner_chain']},{c['receptor_atom']},{c['partner_atom']},"
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
        "Settles a claim that had been carried as UNVERIFIED, meaning stated from\n"
        "memory rather than computed: that the differences Q390R, E412D, R414W and\n"
        "K467R explain why cetuximab does not bind mouse EGFR.\n\n"
        "A residue counts as in contact when any of its heavy atoms, meaning any\n"
        "atom except hydrogen, comes within\n"
        f"{CUTOFF} angstroms of any heavy atom of either chain of the antibody\n"
        "fragment, in structure 6ARU. Hydrogen is left out because it is too light\n"
        "to appear in most measured structures.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
