#!/usr/bin/env python3
"""
07_fab_mutant_check.py — which residues of the antibody in 6ARU are mutated, and
does any of them touch EGFR?

WHY THIS MATTERS
----------------
Structure 6ARU is titled "Structure of Cetuximab Fab mutant in complex with EGFR
extracellular domain". Fab means the gripping arm of an antibody, separated from
the rest of it. The antibody in that file is therefore a modified version of
cetuximab, not the approved drug as it is sold.

Step 04 measured cetuximab's contacts using this file, and three results elsewhere
in the project rest on those contacts:

  - the disproof of an earlier claim about why cetuximab does not bind mouse EGFR
  - the rule not to let our binder contact position 442
  - the overlap figures in the candidate-cluster comparison in step 05

If a mutated residue sits in the interface, the footprint we measured is the
footprint of the modified antibody rather than of cetuximab, and those three
results need qualifying. If no mutated residue is in the interface, the footprint
stands as a description of cetuximab itself.

WHY THIS HAS TO BE COMPUTED
---------------------------
There is no document to look the answer up in:

  - the file's own SEQADV records, which are where a depositor lists differences
    from a reference sequence, cover only the receptor chain. There are none for
    either antibody chain.
  - the entry record at the Protein Data Bank (PDB, the public archive of measured
    three-dimensional structures) names the entry a mutant but lists no
    substitutions.
  - the primary citation is "To Be Published" (Christie M., Christ D., deposited
    2017-08-23, released 2018-08-29), so there is no paper to read.

So the mutations are found by comparing this antibody against a reference
structure of unmodified cetuximab: PDB entry 1YY9, from Li S. et al. (2005),
"Structural basis for inhibition of the epidermal growth factor receptor by
cetuximab", Cancer Cell 7:301-311, doi:10.1016/j.ccr.2005.03.003. That entry is
titled as the Fab fragment of cetuximab/Erbitux/IMC-C225 with no mutant
qualifier.

ONE LIMIT ON THIS METHOD, WHICH THE OUTPUT REPEATS
--------------------------------------------------
1YY9's entry record does not state in so many words that its Fab is unmodified. It
is being treated as the reference because it is the structure that accompanied the
description of how cetuximab works. Any difference this script finds is therefore
a difference between two deposited structures. It is good evidence of a deliberate
mutation where the difference is a single clear substitution, and weaker evidence
where the two files simply disagree. The interface test is what decides whether a
difference matters for our conclusions, and that test does not depend on which of
the two files is the "correct" cetuximab.

Outputs:
  results/findings/07-fab-mutant-check.md
  data/derived/07-fab-differences.csv

Run standalone:  python analysis/07_fab_mutant_check.py
"""

import sys
from pathlib import Path

import requests
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

import egfr_common as common

ROOT = Path(__file__).resolve().parents[1]
STRUCT = ROOT / "data" / "structures"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"

MUTANT_ID = "6ARU"     # cetuximab Fab mutant with EGFR, our structure of record
REFERENCE_ID = "1YY9"  # cetuximab Fab with EGFR, Li et al. 2005

CONTACT_CUTOFF = common.CONTACT_CUTOFF   # 4.5 angstroms, as used in step 04
EPI_START, EPI_END = common.EPI_START, common.EPI_END


def download(pdb_id):
    dest = STRUCT / f"{pdb_id.lower()}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest, "already present"
    dest.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(
        f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=120)
    response.raise_for_status()
    dest.write_bytes(response.content)
    return dest, f"downloaded ({len(response.content):,} bytes)"


def chain_sequence(chain):
    residues = common.protein_residues(chain)
    seq = "".join(common.THREE_TO_ONE.get(r.get_resname(), "X")
                  for r in residues)
    return residues, seq


def identity_to(seq_a, seq_b):
    """Percentage of matching positions when two sequences are lined up.

    Uses local alignment, which finds the best-matching stretch rather than
    forcing the full length of both to correspond. That is the right choice here
    because the two files resolve different amounts of each chain.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignment = aligner.align(seq_a, seq_b)[0]
    matches = sum(1 for x, y in zip(alignment[0], alignment[1])
                  if x == y and x != "-")
    return 100.0 * matches / max(min(len(seq_a), len(seq_b)), 1)


def compare(residues_m, seq_m, residues_r, seq_r):
    """List positions where the two chains differ.

    Returns (differences, aligned_count). Each difference is
    (mutant_resnum, mutant_aa, reference_resnum, reference_aa).
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignment = aligner.align(seq_m, seq_r)[0]

    i_m = i_r = aligned = 0
    differences = []
    for char_m, char_r in zip(alignment[0], alignment[1]):
        if char_m != "-":
            i_m += 1
        if char_r != "-":
            i_r += 1
        if char_m != "-" and char_r != "-":
            aligned += 1
            if char_m != char_r:
                differences.append((residues_m[i_m - 1].id[1], char_m,
                                    residues_r[i_r - 1].id[1], char_r))
    return differences, aligned


def main():
    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    emit("=" * 72)
    emit(f"IS THE ANTIBODY IN {MUTANT_ID} MUTATED AT THE INTERFACE?")
    emit("=" * 72)
    emit()
    emit(f"Compares the two antibody chains of {MUTANT_ID}, which is titled a")
    emit(f"cetuximab Fab mutant, against {REFERENCE_ID}, the reference structure of")
    emit("cetuximab (Li et al. 2005, Cancer Cell 7:301-311). Fab means the gripping")
    emit("arm of an antibody, separated from the rest of it.")
    emit()
    emit("The mutations are not recorded anywhere we could read them: the file's own")
    emit("SEQADV records, which list differences from a reference sequence, exist")
    emit("only for the receptor chain; the Protein Data Bank entry names the mutant")
    emit("but lists no substitutions; and the citation is 'To Be Published'.")
    emit()

    parser = PDBParser(QUIET=True)
    structures, paths = {}, {}
    for pdb_id in (MUTANT_ID, REFERENCE_ID):
        path, status = download(pdb_id)
        paths[pdb_id] = path
        structures[pdb_id] = parser.get_structure(pdb_id, str(path))
        emit(f"   {pdb_id}: {status}")
    emit()

    human = common.human_sequence()

    # Sort each file's chains into receptor and antibody by measured identity to
    # human EGFR, rather than by trusting chain letters.
    inventory = {}
    emit("1. Which chain is which, decided by measured identity to human EGFR")
    emit()
    for pdb_id in (MUTANT_ID, REFERENCE_ID):
        model = structures[pdb_id][0]
        rows = []
        for chain in model:
            residues, seq = chain_sequence(chain)
            if len(residues) < 20:
                continue
            pct = identity_to(human, seq)
            rows.append((chain.id, residues, seq, pct))
        rows.sort(key=lambda r: -r[3])
        receptor = rows[0]
        antibody = rows[1:]
        inventory[pdb_id] = dict(receptor=receptor, antibody=antibody)
        emit(f"   {pdb_id}:")
        emit(f"     chain {receptor[0]}: {len(receptor[1]):3d} residues, "
             f"{receptor[3]:5.1f}% identity to human EGFR  -> receptor")
        for cid, residues, seq, pct in antibody:
            emit(f"     chain {cid}: {len(residues):3d} residues, "
                 f"{pct:5.1f}% identity to human EGFR  -> antibody")
    emit()

    # Pair each antibody chain in the mutant with its counterpart in the
    # reference, by which pairing scores highest. Chain letters are not assumed to
    # mean the same thing in both files.
    emit("2. Pairing the antibody chains between the two files")
    emit()
    emit("   Matched by sequence similarity rather than by chain letter, because")
    emit("   the same letter need not mean the same molecule in two files.")
    emit()
    pairings = []
    used = set()
    for cid_m, res_m, seq_m, _ in inventory[MUTANT_ID]["antibody"]:
        best = None
        for cid_r, res_r, seq_r, _ in inventory[REFERENCE_ID]["antibody"]:
            if cid_r in used:
                continue
            score = identity_to(seq_m, seq_r)
            if best is None or score > best[0]:
                best = (score, cid_r, res_r, seq_r)
        if best is None:
            continue
        used.add(best[1])
        pairings.append((cid_m, res_m, seq_m, best[1], best[2], best[3], best[0]))
        emit(f"   {MUTANT_ID} chain {cid_m} <-> {REFERENCE_ID} chain {best[1]}: "
             f"{best[0]:.1f}% identical")
    emit()

    # Build the set of antibody residues that touch the receptor in 6ARU, so a
    # difference can be tested for whether it is at the interface. This is the
    # same 4.5 angstrom heavy-atom rule step 04 uses; "heavy atom" means any atom
    # except hydrogen, which is absent from these files.
    mutant_model = structures[MUTANT_ID][0]
    receptor_cid = inventory[MUTANT_ID]["receptor"][0]
    receptor_atoms = common.heavy_atoms(
        common.protein_residues(mutant_model[receptor_cid]))
    search = NeighborSearch(receptor_atoms)
    numbering = common.load_numbering()

    interface = {}
    for cid_m, res_m, _, _, _, _, _ in pairings:
        for residue in res_m:
            closest = None
            for atom in residue:
                if atom.element == "H":
                    continue
                for near in search.search(atom.coord, CONTACT_CUTOFF):
                    distance = atom - near
                    partner = numbering.uniprot_of(near.get_parent().id[1])
                    if closest is None or distance < closest[0]:
                        closest = (distance, partner)
            if closest is not None:
                interface[(cid_m, residue.id[1])] = closest

    emit("3. Antibody residues that touch the receptor in " + MUTANT_ID)
    emit()
    emit(f"   Any antibody heavy atom within {CONTACT_CUTOFF} angstroms of any")
    emit("   receptor heavy atom, the same rule step 04 uses.")
    emit(f"   Interface residues found: {len(interface)}")
    emit()

    # Differences, and whether each is at the interface.
    emit("4. Differences between the two antibodies")
    emit()
    all_rows, at_interface = [], []
    for cid_m, res_m, seq_m, cid_r, res_r, seq_r, _ in pairings:
        differences, aligned = compare(res_m, seq_m, res_r, seq_r)
        emit(f"   chain {cid_m} against {REFERENCE_ID} chain {cid_r}: "
             f"{len(differences)} difference(s) over {aligned} aligned positions")
        if not differences:
            emit("     none")
        for resnum_m, aa_m, resnum_r, aa_r in differences:
            key = (cid_m, resnum_m)
            contact = interface.get(key)
            if contact:
                mark = (f"  AT INTERFACE: {contact[0]:.2f} A from EGFR residue "
                        f"{contact[1]}")
                at_interface.append((cid_m, resnum_m, aa_m, aa_r, contact))
            else:
                mark = ""
            emit(f"     {cid_m} {resnum_m}: {aa_r} in {REFERENCE_ID} -> "
                 f"{aa_m} in {MUTANT_ID}{mark}")
            all_rows.append((cid_m, resnum_m, aa_m, cid_r, resnum_r, aa_r,
                             bool(contact),
                             f"{contact[0]:.3f}" if contact else "",
                             contact[1] if contact else ""))
        emit()

    # Verdict.
    emit("5. Verdict")
    emit()
    if not all_rows:
        emit("   The two antibodies have no sequence differences over the residues")
        emit("   both files resolve. Either the mutation lies in a part that one of")
        emit("   the files does not resolve, or it is not a substitution.")
    elif not at_interface:
        emit(f"   {len(all_rows)} difference(s) found, and none of them is within")
        emit(f"   {CONTACT_CUTOFF} angstroms of the receptor.")
        emit()
        emit("   What this means for the project: step 04's contact set is a")
        emit("   description of cetuximab's own footprint, because every residue")
        emit("   doing the touching is the same in the modified antibody as in the")
        emit("   reference structure. The disproof of the earlier species-failure")
        emit("   claim, the rule not to contact position 442, and the overlap")
        emit("   figures in step 05's cluster comparison all stand as they are.")
    else:
        emit(f"   {len(at_interface)} of {len(all_rows)} difference(s) sit within")
        emit(f"   {CONTACT_CUTOFF} angstroms of the receptor:")
        for cid, resnum, aa_m, aa_r, contact in at_interface:
            emit(f"     chain {cid} {resnum}: {aa_r} -> {aa_m}, "
                 f"{contact[0]:.2f} A from EGFR residue {contact[1]}")
        emit()
        emit("   What this means for the project: step 04's contact set describes")
        emit("   the modified antibody rather than cetuximab itself. Three results")
        emit("   that rest on it need qualifying in the decisions log: the disproof")
        emit("   of the earlier species-failure claim, the rule not to contact")
        emit("   position 442, and the overlap figures in step 05's cluster")
        emit("   comparison. Check in particular whether the EGFR residues listed")
        emit("   above fall inside 415-466.")
        touching_epitope = [c[4][1] for c in at_interface
                            if c[4][1] and EPI_START <= c[4][1] <= EPI_END]
        if touching_epitope:
            emit()
            emit(f"   Mutated antibody residues contact these residues inside our")
            emit(f"   epitope: {sorted(set(touching_epitope))}")
        else:
            emit()
            emit("   None of the mutated positions contacts a residue inside")
            emit(f"   {EPI_START}-{EPI_END}, so the overlap figures for our own block")
            emit("   are unaffected even though the footprint as a whole is not")
            emit("   purely cetuximab's.")
    emit()

    # If a mutated residue touches our epitope, the follow-up question is whether
    # the unmodified antibody touches the same place. That decides whether the
    # overlap figures in step 05's cluster comparison describe cetuximab or only
    # this variant.
    emit("5b. Does the unmodified antibody touch the same residues?")
    emit()
    emit(f"   Recomputes the footprint on {REFERENCE_ID}, which has no engineered")
    emit("   substitution at the interface, and compares the two inside our")
    emit(f"   epitope {EPI_START}-{EPI_END}. Both use the same rule: any receptor")
    emit(f"   heavy atom within {CONTACT_CUTOFF} angstroms of any antibody heavy")
    emit("   atom.")
    emit()
    ref_model = structures[REFERENCE_ID][0]
    ref_receptor_cid = inventory[REFERENCE_ID]["receptor"][0]
    ref_receptor_residues = inventory[REFERENCE_ID]["receptor"][1]
    ref_numbering, ref_pct = common.numbering_from_alignment(
        ref_receptor_residues, human)
    ref_offset, ref_share, _ = ref_numbering.dominant_offset()
    emit(f"   {REFERENCE_ID} receptor chain {ref_receptor_cid}: numbering offset")
    emit(f"   UniProt = PDB {ref_offset:+d}, covering "
         f"{100.0 * ref_share:.1f}% of matched residues, derived by alignment.")
    emit()

    ref_antibody_cids = [cid for cid, _, _, _ in inventory[REFERENCE_ID]["antibody"]]
    ref_contacts = common.contacts_to_partner(
        ref_model, ref_receptor_cid, ref_antibody_cids, ref_numbering,
        cutoff=CONTACT_CUTOFF, restrict_to=(EPI_START, EPI_END))
    mutant_contacts = common.contacts_to_partner(
        mutant_model, receptor_cid,
        [cid for cid, _, _, _, _, _, _ in pairings], numbering,
        cutoff=CONTACT_CUTOFF, restrict_to=(EPI_START, EPI_END))

    in_both = sorted(set(ref_contacts) & set(mutant_contacts))
    only_mutant = sorted(set(mutant_contacts) - set(ref_contacts))
    only_ref = sorted(set(ref_contacts) - set(mutant_contacts))

    emit(f"   Inside {EPI_START}-{EPI_END}:")
    emit(f"     touched in both: {in_both}")
    emit(f"     touched only in {MUTANT_ID} (the mutant): {only_mutant}")
    emit(f"     touched only in {REFERENCE_ID} (the reference): {only_ref}")
    emit()

    anchors_ref = [p for p in sorted(common.ANCHORS) if p in ref_contacts]
    anchors_mut = [p for p in sorted(common.ANCHORS) if p in mutant_contacts]
    emit(f"   Anchors touched in {MUTANT_ID}: "
         f"{[f'{common.ANCHORS[p]}{p}' for p in anchors_mut] or 'none'}")
    emit(f"   Anchors touched in {REFERENCE_ID}: "
         f"{[f'{common.ANCHORS[p]}{p}' for p in anchors_ref] or 'none'}")
    emit()
    for pos in sorted(set(anchors_mut) | set(anchors_ref)):
        label = f"{common.ANCHORS[pos]}{pos}"
        d_m = (f"{mutant_contacts[pos]['min_dist']:.2f} A"
               if pos in mutant_contacts else "not touched")
        d_r = (f"{ref_contacts[pos]['min_dist']:.2f} A"
               if pos in ref_contacts else "not touched")
        emit(f"     {label}: {MUTANT_ID} {d_m}; {REFERENCE_ID} {d_r}")
    emit()
    if set(anchors_mut) != set(anchors_ref):
        emit("   The two structures disagree about which of our anchors the")
        emit("   antibody touches. The figure to use for the cluster comparison is")
        emit(f"   the one from {REFERENCE_ID}, because that describes cetuximab")
        emit("   rather than an engineered variant. Record both.")
    else:
        emit("   Both structures touch the same anchors, so the overlap figures in")
        emit("   step 05's cluster comparison hold for cetuximab and not only for")
        emit("   this variant.")
    emit()

    emit("6. Limits")
    emit()
    emit("   The two structures were solved separately, at different resolutions")
    emit("   and in different crystals, so a residue right at the cutoff can fall")
    emit("   on either side for reasons unconnected to the mutation. Treat a")
    emit("   difference of a few tenths of an angstrom as noise and a residue")
    emit("   appearing in one list but not the other as worth checking rather than")
    emit("   settled.")
    emit()
    emit(f"   {REFERENCE_ID}'s entry record does not state that its Fab is")
    emit("   unmodified. It is used as the reference because it is the structure")
    emit("   published with the account of how cetuximab works. A difference found")
    emit("   here is a difference between two deposited structures, which is good")
    emit("   evidence of a deliberate change where it is a single clean")
    emit("   substitution and weaker where the files merely disagree.")
    emit()
    emit("   The comparison covers only residues both files resolve. A mutation in")
    emit("   a stretch missing from either structure would not appear here.")
    emit()
    emit("   The interface test uses the positions in this one crystal structure.")
    emit("   It answers whether a mutated residue touches the receptor in this")
    emit("   structure, which is the question that matters for step 04, and does")
    emit("   not measure how much any contact contributes to binding.")
    emit()

    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "07-fab-differences.csv").open("w") as fh:
        fh.write("mutant_chain,mutant_resnum,mutant_aa,reference_chain,"
                 "reference_resnum,reference_aa,at_interface,distance_a,"
                 "contacts_egfr_uniprot_pos\n")
        for row in all_rows:
            fh.write(",".join(str(x) for x in row) + "\n")
    emit("Wrote data/derived/07-fab-differences.csv")

    emit()
    emit("=" * 72)
    emit(f"RESULT: {len(all_rows)} antibody difference(s), "
         f"{len(at_interface)} at the interface.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "07-fab-mutant-check.md").write_text(
        f"# Is the antibody in {MUTANT_ID} mutated at the interface?\n\n"
        f"Computed output of `analysis/07_fab_mutant_check.py`. Do not hand-edit.\n\n"
        f"{MUTANT_ID} is titled a cetuximab Fab mutant, and step 04 measured the\n"
        "antibody footprint using it. This checks whether any modified residue is\n"
        "in the interface, which decides whether that footprint describes cetuximab\n"
        "or only this modified version of it.\n\n"
        f"Reference for comparison: {REFERENCE_ID}, Li S. et al. (2005), \"Structural\n"
        "basis for inhibition of the epidermal growth factor receptor by\n"
        "cetuximab\", Cancer Cell 7:301-311, doi:10.1016/j.ccr.2005.03.003.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
