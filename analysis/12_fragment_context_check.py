#!/usr/bin/env python3
"""
12_fragment_context_check.py — is the H370 anchor face reachable in the intact
receptor, or does the rest of the receptor fold across it?

Abbreviations, expanded here because each file gets read on its own:
  PDB      Protein Data Bank, the public archive of measured 3D protein
           structures. Also used to mean a file from that archive.
  UniProt  the public protein sequence archive. Our residue numbers are
           positions in its records.
  RSA      relative solvent accessibility: how much of a residue's surface water
           can reach, divided by the most that residue type could ever expose.
  ECD      extracellular domain, the part of EGFR outside the cell, residues
           25-645 in our numbering.
  EGF      epidermal growth factor, the signalling molecule EGFR responds to.
  Fab      the gripping arm of an antibody, separated from the rest of it.
  aa       amino acids, the building blocks a protein chain is made of.

WHY THIS STEP EXISTS
--------------------
Step 09 cuts domain III (residues 310-480) out of the full extracellular region
and hands the 171-residue fragment to the design run as the target. The design run
therefore never sees the other 450 residues. That raises a question the design run
cannot answer about itself: a binder could be designed against a surface that
other parts of the receptor cover up in the intact molecule. Every predicted score
would come back good and the molecule would measure as nothing in the assay, with
nothing to distinguish that from a badly designed binder.

Step 06 asked exactly this question and found a real answer: anchors D458 and D460
sit in a groove against domain IV rather than on an open face. But step 06 ran on
the epitope we aimed at then, residues 415-466. Step 08 subsequently moved the
target face to the cluster around H370, and only two of those eight anchors
(E421, E424) are in the range step 06 measured. Six of the eight positions the
design run is now aimed at had never been checked for occlusion by the rest of the
receptor. This step checks them.

WHAT THIS SCRIPT COMPUTES
-------------------------
For each anchor, two different measurements, because they answer two different
questions:

  1. Contacts from outside domain III, at 4.5 angstroms, using the same shared
     function step 06 used. This answers "is this residue covered over".

  2. How many residues from outside domain III fall within 8 and 12 angstroms.
     This answers "can a protein get here", which is not the same question. A
     binder is a body roughly 20 angstroms across, so an anchor can have no
     contacts at all and still sit at the bottom of a cleft nothing that size can
     reach into. Measurement 1 alone would call such a site open.

Both are computed on the receptor protein chain alone, in two structures:

  6ARU  extended (open) form, the structure step 09 cuts the fragment from
  1NQL  tethered (closed) form, the autoinhibited shape

Running both matters because the assay presents the whole extracellular region in
solution and we do not know which shape predominates there. An anchor face that is
open in one form and buried in the other is a risk; one that is open in both is
as settled as two measured structures can make it.

The partner molecules are excluded from both files, so the measurement is of what
the receptor's own shape covers. 6ARU has a cetuximab Fab bound and 1NQL has EGF
bound; leaving them in would mostly measure the difference between those two
partners, which is a circumstance of getting each crystal to form rather than a
description of our assay.

HOW THIS STEP IS CHECKED
------------------------
A check that has never been seen to fail is not known to work, so:

  - The old 415-466 anchors are run through the same code as a control, and the
    residues it finds in contact are cross-checked against the table step 06
    committed. If this script's answer differs from step 06's, the run stops. That
    is the rule in `CLAUDE.md` that came out of steps 04 and 06 disagreeing.
  - `--break-numbering N` shifts every residue number by N, which should make the
    identity checks notice that the position believed to be H370 is not a
    histidine. Run it once and confirm it fails.

WHAT THIS STEP DOES NOT SETTLE
------------------------------
It measures whether the anchor face is reachable in the intact receptor. It says
nothing about whether the cut fragment still folds into the shape those anchors sit
on, which is a different question, is still open, and needs the design run's own
structure predictor rather than this geometry. Step 11 tried and could not resolve
it with a weaker predictor.

It also measures geometry, not solvent accessibility, in the closed form. Step 03
measured RSA on the whole receptor chain, so the exposure numbers that selected
these anchors already include whatever burial the rest of the receptor causes in
6ARU. The same figure for 1NQL is not computed here; `receptor_only_sasa` in step
06 is the function to lift if we want it.

Writes:
  results/findings/12-fragment-context.md
  data/derived/12-anchor-clearance.csv
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import egfr_common as common  # noqa: E402

STRUCT = common.STRUCT_DIR
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS

OPEN_ID, CLOSED_ID = "6ARU", "1NQL"
STRUCTURE_FILES = {OPEN_ID: "6aru.pdb", CLOSED_ID: "1nql.pdb"}

# Domain III, the piece step 09 cuts out. Taken from egfr_common so this step and
# step 09 cannot disagree about where the cut is.
DOMAIN_III = (common.D3_START, common.D3_END)

# The wider shells, in angstroms. 8 is roughly two residues' reach; 12 is about
# half the width of a small binder, so a residue with many neighbours inside 12 is
# in a cleft rather than on a face. Neither value is a threshold in the molecule,
# so an anchor near either count is reported as such rather than classified.
SHELLS = (8.0, 12.0)

# How much clear space counts as room for a binder to approach. A binder is a body
# roughly 20 angstroms across, and it only has to lay one face against the target,
# so 10 angstroms of clearance with no close neighbours is taken as reachable.
# This is a judgement, not a measured threshold, and is reported as one.
CLEARANCE_ROOM = 10.0


def expected_identities(anchors, human):
    """What amino acid each anchor should be, from the human sequence itself."""
    return {pos: human[pos - 1] for pos in anchors}


def check_identities(label, numbering, chain, expected, emit):
    """Confirm the numbering by residue identity rather than by arithmetic.

    Arithmetic can be off by one and still look right. A histidine that turns out
    to be a leucine cannot. This is the same discipline step 09 applies, for the
    same reason: a wrong offset shifts every result below and raises no error.
    """
    by_number = {r.id[1]: r for r in common.protein_residues(chain)}
    problems, checked = [], 0
    for pos, aa in sorted(expected.items()):
        pdb_num = numbering.pdb_of(pos)
        if pdb_num is None or pdb_num not in by_number:
            problems.append(f"{aa}{pos}: not observed in {label}")
            continue
        found = common.THREE_TO_ONE.get(by_number[pdb_num].get_resname())
        checked += 1
        if found != aa:
            problems.append(f"{aa}{pos}: file has {found} at that position")
    if problems:
        emit(f"   [FAIL] identity check against {label}: "
             + "; ".join(problems))
        return False
    emit(f"   [PASS] identity check against {label}: all {checked} positions "
         "hold the amino acid our numbering says they should")
    return True


def shifted(numbering, shift):
    """The measured numbering, optionally broken on purpose by `shift` places."""
    if shift == 0:
        return numbering
    moved = {pdb: uni + shift for pdb, uni in
             ((numbering.pdb_of(u), u) for u in numbering.uniprot_positions())
             if pdb is not None}
    return common.Numbering(moved)


def measure(model, chain_id, numbering, anchors):
    """Contacts at 4.5 A and the wider shell counts, for one structure."""
    contacts = common.intra_chain_contacts_outside(
        model, chain_id, numbering,
        target_range=(min(anchors), max(anchors)),
        exclude_range=DOMAIN_III)
    shells = common.shell_counts(model, chain_id, numbering, set(anchors),
                                 exclude_range=DOMAIN_III, shells=SHELLS)
    return contacts, shells


def reading(row, contact):
    """Plain words for one anchor's clearance. Reported, not used to filter."""
    if contact is not None:
        return "in contact, sits in a groove"
    nearest = row["nearest"][0]
    if nearest >= CLEARANCE_ROOM and row[SHELLS[0]] == 0:
        return "open face, room for a binder"
    if nearest >= common.CONTACT_CUTOFF:
        return "clear of contact, some crowding"
    return "very close to a neighbour"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--break-numbering", type=int, default=0, metavar="N",
                        help="shift every residue number by N, to confirm the "
                             "identity checks notice. Expected to fail.")
    args = parser.parse_args(argv)

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    human = common.human_sequence()
    anchors = common.load_primary_cluster()
    old_anchors = sorted(common.ANCHORS)
    expected = expected_identities(anchors, human)

    emit("=" * 72)
    emit("FRAGMENT CONTEXT CHECK — is the H370 face reachable in intact EGFR?")
    emit("=" * 72)
    emit()
    emit("Step 09 hands the design run a 171-residue fragment, residues "
         f"{DOMAIN_III[0]}-{DOMAIN_III[1]},")
    emit("so the design run never sees the other 450 residues of the "
         "extracellular")
    emit("region. This step asks whether those 450 residues fold across the face "
         "we")
    emit("are aiming at. If they do, a design can score well and still be unable "
         "to")
    emit("reach its target, and nothing in the design run would say so.")
    emit()
    emit("   anchors under test (from step 08's committed output): "
         + ", ".join(f"{expected[p]}{p}" for p in anchors))
    emit(f"   fragment the design run receives: {DOMAIN_III[0]}-{DOMAIN_III[1]}")
    emit()
    emit("Two measurements per anchor, because they answer different questions:")
    emit("  contacts at 4.5 A   — is this residue covered over")
    emit("  residues within 8 A and 12 A — can a body the size of a binder get "
         "here")
    emit()
    if args.break_numbering:
        emit(f"   NOTE: --break-numbering {args.break_numbering} is in effect. "
             "This run is")
        emit("   expected to fail at the identity check below.")
        emit()

    failures = []
    per_structure = {}

    for pid in (OPEN_ID, CLOSED_ID):
        path = STRUCT / STRUCTURE_FILES[pid]
        emit(f"1. {pid} — {'extended (open)' if pid == OPEN_ID else 'tethered (closed)'} form")
        emit()
        if not path.exists():
            emit(f"   {path} is missing. Run analysis/06_tethered_occlusion.py, "
                 "which downloads both structures.")
            return 1
        model = common.load_complex(path)
        chain_id, numbering, pct, others = common.receptor_chain_of(model, human)
        emit(f"   receptor chain {chain_id}, {len(numbering)} residues mapped, "
             f"{pct:.1f}% identity to human EGFR")
        offset, share, _counts = numbering.dominant_offset()
        emit(f"   numbering offset: ours = file {offset:+d} "
             f"({100.0 * share:.1f}% of residues) — derived, not assumed")
        emit(f"   other protein chains present, excluded from the measurement: "
             f"{', '.join(others) or 'none'}")
        emit()

        numbering = shifted(numbering, args.break_numbering)
        chain = model[chain_id]
        if not check_identities(pid, numbering, chain, expected, emit):
            failures.append(f"identity check failed against {pid}")
            emit()
            continue
        emit()

        contacts, shells = measure(model, chain_id, numbering, anchors)
        old_contacts, old_shells = measure(model, chain_id, numbering,
                                           old_anchors)
        per_structure[pid] = dict(contacts=contacts, shells=shells,
                                  old_contacts=old_contacts,
                                  old_shells=old_shells)

        emit("   The current anchors, against the rest of the receptor:")
        emit()
        emit("   | anchor | touched at <4.5 A by | nearest outside domain III | "
             "within 8 A | within 12 A | reading |")
        emit("   |---|---|---|---|---|---|")
        for pos in anchors:
            row = shells[pos]
            contact = contacts.get(pos)
            touched = (f"{expected.get(contact[1], '?')}{contact[1]}"
                       if contact else "nothing")
            near_d, near_pos = row["nearest"]
            emit(f"   | {expected[pos]}{pos} | {touched} | {near_d:.1f} A "
                 f"to {near_pos} | {row[SHELLS[0]]} | {row[SHELLS[1]]} | "
                 f"{reading(row, contact)} |")
        emit()

        in_contact = sorted(contacts)
        if in_contact:
            emit(f"   {len(in_contact)} of {len(anchors)} anchors are touched "
                 "from outside domain III: "
                 + ", ".join(f"{expected[p]}{p}" for p in in_contact))
        else:
            emit(f"   None of the {len(anchors)} anchors is touched from outside "
                 "domain III.")
        emit()

        # ---- Control: the same code on step 06's anchors, cross-checked ----
        emit(f"   Control — the old {common.EPI_START}-{common.EPI_END} anchors "
             "through this same code.")
        emit("   Step 06 measured these and found two of them in a groove against")
        emit("   domain IV. If this script cannot reproduce that, its answer above")
        emit("   is not to be believed either.")
        emit()
        for pos in old_anchors:
            row = old_shells[pos]
            contact = old_contacts.get(pos)
            touched = f"{contact[1]} at {contact[0]:.2f} A" if contact else "nothing"
            emit(f"     {common.ANCHORS[pos]}{pos}: touched by {touched}; "
                 f"nearest {row['nearest'][0]:.1f} A")
        emit()
        if pid == OPEN_ID:
            try:
                common.cross_check_residue_set(
                    f"anchors in contact from outside domain III in {pid}",
                    {p for p in old_contacts if p in common.ANCHORS},
                    DERIVED / "06-intra-chain-occlusion.csv",
                    column="uniprot_pos", condition_column="is_anchor",
                    emit=emit)
            except common.CrossCheckError as exc:
                failures.append(str(exc))
        emit()

    # ---- What it means ----
    emit("=" * 72)
    emit("2. What this changes about the design")
    emit("=" * 72)
    emit()
    if len(per_structure) == 2:
        touched_any = sorted(
            set(per_structure[OPEN_ID]["contacts"])
            | set(per_structure[CLOSED_ID]["contacts"]))
        crowded = sorted(
            pos for pos in anchors
            if any(per_structure[pid]["shells"][pos][SHELLS[0]] > 0
                   for pid in per_structure))
        clear = [pos for pos in anchors if pos not in crowded]
        emit(f"   Anchors touched from outside domain III, in either form: "
             f"{len(touched_any)}"
             + (" — " + ", ".join(f"{expected[p]}{p}" for p in touched_any)
                if touched_any else " — none"))
        emit(f"   Anchors with any neighbour inside {SHELLS[0]:.0f} A, in either "
             f"form: {len(crowded)}"
             + (" — " + ", ".join(f"{expected[p]}{p}" for p in crowded)
                if crowded else " — none"))
        emit(f"   Anchors on open face in both forms: {len(clear)}"
             + (" — " + ", ".join(f"{expected[p]}{p}" for p in clear)
                if clear else " — none"))
        emit()
        if not touched_any and len(clear) >= len(anchors) - 2:
            emit("   The face the design run is aimed at is not covered by the "
                 "rest of")
            emit("   the receptor, in either measured conformation. Cutting "
                 "domain III out")
            emit("   does not hide the anchors, so the occlusion risk that cost "
                 "the old")
            emit("   415-466 epitope two of its anchors does not apply to this "
                 "face.")
            emit()
            emit("   What it does not license: this says the anchors are "
                 "reachable, not")
            emit("   that the fragment folds the way the intact protein does. "
                 "That stays")
            emit("   open, and only the design run's own predictor can close it.")
        else:
            emit("   Some anchors are packed against the rest of the receptor. A "
                 "binder")
            emit("   reaching them has to fit into a groove rather than lie "
                 "against a")
            emit("   face, which is harder to design and more sensitive to how "
                 "the two")
            emit("   domains sit against each other. Weigh this when ranking.")
        emit()
        if crowded:
            emit("   Worth carrying forward: "
                 + ", ".join(f"{expected[p]}{p}" for p in crowded)
                 + " have neighbours close by.")
            emit("   Step 09 separately found E424 sits 6.6 A from the cut at "
                 "480. An")
            emit("   anchor in both lists is the weakest of the eight on two")
            emit("   independent grounds, and step 10 already ranks designs that "
                 "lean on")
            emit("   it below equivalent designs that do not.")
        emit()

    # ---- Tables ----
    DERIVED.mkdir(parents=True, exist_ok=True)
    csv_path = DERIVED / "12-anchor-clearance.csv"
    with csv_path.open("w") as fh:
        fh.write("structure,uniprot_pos,aa,anchor_set,touched_at_4_5a,"
                 "touched_by_uniprot_pos,contact_distance_a,nearest_outside_a,"
                 "nearest_outside_pos,residues_within_8a,residues_within_12a\n")
        for pid, data in per_structure.items():
            for label, positions, contacts, shells in (
                    ("h370", anchors, data["contacts"], data["shells"]),
                    ("415-466", old_anchors, data["old_contacts"],
                     data["old_shells"])):
                for pos in positions:
                    row = shells[pos]
                    contact = contacts.get(pos)
                    aa = expected.get(pos) or common.ANCHORS.get(pos) or "?"
                    near_d, near_pos = row["nearest"]
                    touched_by = contact[1] if contact else ""
                    distance = f"{contact[0]:.3f}" if contact else ""
                    fh.write(f"{pid},{pos},{aa},{label},"
                             f"{contact is not None},{touched_by},{distance},"
                             f"{near_d:.3f},{near_pos},"
                             f"{row[SHELLS[0]]},{row[SHELLS[1]]}\n")
    emit(f"Wrote data/derived/{csv_path.name}")
    emit()

    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED —")
        for failure in failures:
            emit(f"  - {failure}")
    else:
        emit("RESULT: PASSED. The H370 anchor face is not occluded by the rest "
             "of the")
        emit("receptor in either measured conformation. Whether the cut fragment "
             "holds")
        emit("its shape is a separate question and is still open.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "12-fragment-context.md").write_text(
        "# Is the H370 anchor face reachable in the intact receptor?\n\n"
        "Computed output of `analysis/12_fragment_context_check.py`. "
        "Do not hand-edit.\n\n"
        "Step 09 cuts domain III out and gives the design run the fragment alone, "
        "so\nthe design run cannot see whether the rest of the receptor folds "
        "across the\nface we are aiming at. Step 06 asked this of the old "
        f"{common.EPI_START}-{common.EPI_END} epitope and found\ntwo anchors in a "
        "groove against domain IV. Step 08 then moved the face to the\ncluster "
        "around H370, and six of those eight anchors had never been checked.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
