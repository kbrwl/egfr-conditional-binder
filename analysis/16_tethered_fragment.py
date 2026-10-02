#!/usr/bin/env python3
"""
16_tethered_fragment.py — does the fragment cut from the tethered structure differ
from the one cut from the extended structure, where the anchors are?

Abbreviations, expanded here because each file gets read on its own:
  EGFR     epidermal growth factor receptor, the protein we design against.
  UniProt  the public protein sequence archive. Every residue number in this
           project is a position in its human record P00533.
  PDB      the Protein Data Bank, the public archive of measured 3D structures.
           Also used to mean a file from that archive.
  CA       the alpha carbon, one atom present in every residue, used here as a
           stand-in for the residue's position along the backbone.
  RMSD     root-mean-square deviation: the average distance left between matched
           atoms once two structures have been laid on top of each other. A small
           number means the same shape in the same place.
  EGF      epidermal growth factor, the signalling molecule EGFR responds to.
  Fab      the gripping arm of an antibody, separated from the rest of it.

WHY THIS STEP EXISTS
--------------------
The organisers said the screen uses the tethered form of EGFR, the closed inactive
one (docs/competition-qa-log.md). Our structure of record, 6ARU, is the extended
form held by a cetuximab Fab, and the 171-residue fragment step 09 hands to the
design run is cut from it. If the target in the assay is tethered, the model the
designs are shaped against should be the tethered one, unless the two shapes are
the same where the anchors are.

Two things already in hand keep this from being a crisis. Step 12 measured the
H370 face in both structures and found no anchor contacted from outside domain III
in either. And step 06 found domain III itself superposes between the two at 1.08
angstroms RMSD over all 171 residues. What those do not say is whether the
anchors' charged groups sit in the same places, which is what the pairing rule is
built on. This step measures that directly, and decides on the measurement and not
on the argument.

WHAT THIS SCRIPT COMPUTES
-------------------------
The same residues, 310 to 480 in our numbering, are taken from 6ARU (extended) and
from 1NQL (tethered), and 1NQL is laid on 6ARU by fitting all the shared CA atoms.
Then, in that common frame:

  - the whole fragment: RMSD, and which residues deviate most;
  - the anchor face: every residue with an atom within 8 A of any anchor's charged
    tip, so that a comparison of eight points is not the only evidence;
  - each anchor: how far its CA moved, and how far its charged tip moved. The tip
    is the atom the pairing rule uses, so it is the quantity that decides.

THE RULE FOR DECIDING, FIXED BEFORE THE MEASUREMENT
---------------------------------------------------
The two cuts are called the same where it matters if no anchor's charged tip moves
by TIP_TOLERANCE or more once the structures are fitted. That is 2 angstroms, a
third of the 6 angstrom reach step 10 allows between two charged groups. It is a
judgement and not a measured threshold, and the moves are printed so that anyone
can apply their own. The two structures are 3.2 A (6ARU) and 2.8 A (1NQL)
resolution, so differences of a few tenths of an angstrom are inside what they can
resolve.

If the cuts are the same, 6ARU stays, for that stated reason. If any anchor tip
moves by the tolerance or more, the tethered cut is the one to design against, and
this script has already written it.

HOW THIS STEP IS CHECKED
------------------------
  - Fitting 1NQL onto 6ARU over the whole of domain III must reproduce the 1.08 A
    that step 06 committed. The run stops if it does not.
  - The deciding rule is tested on constructed displacements, run every time: no
    movement and a 1.9 A movement must give 'same', and a 2.0 A and a 5 A movement
    must give 'differs'. `--break-rule decision` removes the threshold, and the
    matching test must then fail.
  - Every anchor is confirmed to be the amino acid our numbering says in both
    structures, by identity and not by arithmetic. `--break-numbering N` shifts the
    numbering and should fail.

Writes:
  results/findings/16-tethered-fragment.md
  data/derived/16-anchor-displacement.csv
  data/structures/1nql_domain3.pdb      (the tethered cut, kept on the shelf)

Run standalone:  python analysis/16_tethered_fragment.py
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, Select, Superimposer

sys.path.insert(0, str(Path(__file__).resolve().parent))

import egfr_common as common  # noqa: E402

STRUCT = common.STRUCT_DIR
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS

FILES = {"6ARU": "6aru.pdb", "1NQL": "1nql.pdb"}
FRAGMENT = (common.D3_START, common.D3_END)

# A third of the 6 angstrom reach step 10 allows between two charged groups. A
# judgement, not a measured threshold; see the module docstring.
TIP_TOLERANCE = 2.0

# How close an atom must be to an anchor's charged tip for its residue to count as
# part of the anchor face. 8 A is roughly two residues' reach.
FACE_RADIUS = 8.0

# The fewest shared CA atoms for a superposition to mean anything.
MIN_ATOMS = 30

BROKEN_RULE = None
BREAKABLE_RULES = {
    "decision": (
        "deciding rule",
        "an anchor whose charged tip had moved by several angstroms would still "
        "be called the same, and the tethered cut would never be chosen"),
}


def decide(tip_moves):
    """'same' if every anchor tip moved less than the tolerance, else 'differs'."""
    if BROKEN_RULE == "decision":
        return "same"
    return "same" if all(m < TIP_TOLERANCE for m in tip_moves) else "differs"


def run_tests(emit):
    """The deciding rule on constructed displacements."""
    emit("   The deciding rule, on constructed displacements of anchor tips (A):")
    emit()
    cases = [
        ([0.0] * 8, "same", "no movement"),
        ([0.3, 1.9, 0.5, 0.2, 0.1, 0.8, 1.2, 0.4], "same",
         "the largest at 1.9 A, just under the tolerance"),
        ([0.3, 2.0, 0.5, 0.2, 0.1, 0.8, 1.2, 0.4], "differs",
         "the largest at exactly the tolerance"),
        ([0.3, 0.2, 5.0, 0.2, 0.1, 0.8, 1.2, 0.4], "differs",
         "one anchor 5 A out, the rest still"),
    ]
    failures = []
    for moves, want, why in cases:
        got = decide(moves)
        ok = got == want
        emit(f"     [{'PASS' if ok else 'FAIL'}] {why}: {got}")
        if not ok:
            failures.append(f"deciding rule: {why} gave {got}, expected {want}")
    emit()
    return failures


class FragmentSelect(Select):
    """Keep one chain and exactly the residues whose file numbers are listed."""

    def __init__(self, chain_id, keep):
        self.chain_id = chain_id
        self.keep = set(keep)

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[1] in self.keep


def load(pid, human, shift, emit):
    path = STRUCT / FILES[pid]
    if not path.exists():
        raise SystemExit(f"{path} is missing. Run analysis/06_tethered_occlusion.py, "
                         "which downloads both structures.")
    model = common.load_complex(path)
    chain_id, numbering, pct, _others = common.receptor_chain_of(model, human)
    offset, share, _ = numbering.dominant_offset()
    emit(f"   {pid}: receptor chain {chain_id}, {len(numbering)} residues mapped, "
         f"{pct:.1f}% identity; ours = file {offset:+d} ({100 * share:.1f}%), derived")
    if shift:
        numbering = common.Numbering(
            {numbering.pdb_of(u): u + shift for u in numbering.uniprot_positions()
             if numbering.pdb_of(u) is not None})
    by_file = {r.id[1]: r for r in common.protein_residues(model[chain_id])}
    by_ours = {}
    for file_num, res in by_file.items():
        ours = numbering.uniprot_of(file_num)
        if ours is not None:
            by_ours[ours] = res
    return dict(id=pid, model=model, chain=chain_id, numbering=numbering,
                by_file=by_file, by_ours=by_ours)


def shared_point(res_a, res_b, aa):
    """The same atom, or atoms, from both structures, so that what is compared is a
    movement and not the distance between two different atoms.

    Returns (point in a, point in b, which atoms, whether it is a fallback).
    The preferred point is the charged tip step 05 and step 14 use. If the tip is
    not resolved in both structures, the first atom of the side chain (CB) is used
    in both, and failing that the alpha carbon (CA). A fallback is reported as one:
    a movement of CB understates a movement of the tip.

    The first version of this function took each structure's own best atom. Where
    1NQL did not resolve E344's tip, that compared 6ARU's CD against 1NQL's CB and
    reported the 3 A between two different atoms as a movement.
    """
    wanted = common.FUNCTIONAL_ATOM[aa]
    if all(a in res_a for a in wanted) and all(a in res_b for a in wanted):
        if aa == "H":
            return (np.mean([res_a[a].coord for a in wanted], axis=0),
                    np.mean([res_b[a].coord for a in wanted], axis=0),
                    "imidazole centroid", False)
        return res_a[wanted[0]].coord, res_b[wanted[0]].coord, wanted[0], False
    for atom in ("CB", "CA"):
        if atom in res_a and atom in res_b:
            return res_a[atom].coord, res_b[atom].coord, atom, True
    raise ValueError("no atom shared by both structures")


def read_step06_rmsd():
    """The domain III RMSD step 06 committed, from its findings file. Raises if the
    row is not there, so a reworded file cannot make this control silently pass."""
    text = (FINDINGS / "06-tethered-occlusion.md").read_text()
    m = re.search(r"domain III \(310-480\)\s*\|\s*(\d+)\s*\|\s*([0-9.]+) A", text)
    if not m:
        raise common.CrossCheckError(
            "could not find the domain III row in 06-tethered-occlusion.md")
    return int(m.group(1)), float(m.group(2))


def main(argv=None):
    global BROKEN_RULE

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--break-rule", choices=list(BREAKABLE_RULES), default=None,
                        help="switch the deciding rule off, to confirm its test "
                             "then fails. Expected to fail.")
    parser.add_argument("--break-numbering", type=int, default=0, metavar="N",
                        help="shift every residue number by N, to confirm the "
                             "identity checks notice. Expected to fail.")
    args = parser.parse_args(argv)
    BROKEN_RULE = args.break_rule
    write_files = BROKEN_RULE is None and not args.break_numbering

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    human = common.human_sequence()
    anchors = common.load_primary_cluster()
    aa = {p: human[p - 1] for p in anchors}
    failures = []

    emit("=" * 72)
    emit("TETHERED FRAGMENT — is the cut from 1NQL different where the anchors are?")
    emit("=" * 72)
    emit()
    emit("The organisers say the screen uses the tethered form. Our fragment is cut")
    emit("from 6ARU, the extended form. Step 12 found no anchor covered in either,")
    emit("and step 06 found domain III superposes at 1.08 A. Neither says whether the")
    emit("anchors' charged groups sit in the same places, which is what the pairing")
    emit("rule is built on. This measures it, and decides on the measurement.")
    emit()
    if BROKEN_RULE or args.break_numbering:
        emit("   NOTE: a break flag is in effect. This run is expected to fail and")
        emit("   writes no files.")
        emit()

    emit("1. The deciding rule, tested every time")
    emit()
    failures.extend(run_tests(emit))

    emit("2. The two structures")
    emit()
    six = load("6ARU", human, args.break_numbering, emit)
    one = load("1NQL", human, args.break_numbering, emit)
    emit()

    for st in (six, one):
        problems = []
        for pos in anchors:
            res = st["by_ours"].get(pos)
            found = common.THREE_TO_ONE.get(res.get_resname()) if res else None
            if found != aa[pos]:
                problems.append(f"{aa[pos]}{pos} reads as {found}")
        if problems:
            emit(f"   [FAIL] identity check, anchors, {st['id']}: "
                 + "; ".join(problems))
            failures.append(f"identity check failed against {st['id']}")
        else:
            emit(f"   [PASS] identity check, anchors, {st['id']}: all "
                 f"{len(anchors)} hold the amino acid our numbering says")
    emit()
    if failures:
        return finish(out, failures, emit, write_files, None)

    # ---- 3. the fragment, residue by residue ---------------------------------
    wanted = range(FRAGMENT[0], FRAGMENT[1] + 1)
    have = {st["id"]: {p for p in wanted if p in st["by_ours"]
                       and "CA" in st["by_ours"][p]} for st in (six, one)}
    shared = sorted(have["6ARU"] & have["1NQL"])
    emit(f"3. The fragment, residues {FRAGMENT[0]}-{FRAGMENT[1]} "
         f"({FRAGMENT[1] - FRAGMENT[0] + 1} positions)")
    emit()
    for pid in ("6ARU", "1NQL"):
        missing = sorted(set(wanted) - have[pid])
        emit(f"   {pid}: {len(have[pid])} residues present"
             + (f"; missing {common_ranges(missing)}" if missing else ""))
    emit(f"   shared by both and used for the fit: {len(shared)}")
    emit()
    if len(shared) < MIN_ATOMS:
        failures.append(f"only {len(shared)} shared residues; need {MIN_ATOMS}")
        return finish(out, failures, emit, write_files, None)

    sup = Superimposer()
    sup.set_atoms([six["by_ours"][p]["CA"] for p in shared],
                  [one["by_ours"][p]["CA"] for p in shared])
    rot, tran = sup.rotran
    rmsd = float(sup.rms)

    def moved(xyz):
        return np.asarray(xyz, dtype=float) @ rot + tran

    ca_dev = {p: float(np.linalg.norm(moved(one["by_ours"][p]["CA"].coord)
                                      - six["by_ours"][p]["CA"].coord))
              for p in shared}

    emit(f"4. The whole fragment: 1NQL laid on 6ARU over {len(shared)} CA atoms")
    emit()
    emit(f"   RMSD {rmsd:.2f} A")
    committed_n, committed = read_step06_rmsd()
    ok = abs(round(rmsd, 2) - committed) < 0.006 and len(shared) == committed_n
    emit(f"   [{'PASS' if ok else 'FAIL'}] control against step 06: "
         f"{rmsd:.2f} A over {len(shared)} residues here, {committed:.2f} A over "
         f"{committed_n} there")
    if not ok:
        failures.append(f"control: {rmsd:.2f} A over {len(shared)} here, "
                        f"{committed:.2f} A over {committed_n} in step 06")
    worst = sorted(ca_dev.items(), key=lambda kv: -kv[1])[:6]
    emit("   The six residues that deviate most (CA, after the fit): "
         + ", ".join(f"{human[p - 1]}{p} {d:.1f} A" for p, d in worst))
    emit()

    # ---- 5. the anchor face ---------------------------------------------------
    tips6, tips1, labels = {}, {}, {}
    for pos in anchors:
        t6, t1, label, fb = shared_point(six["by_ours"][pos], one["by_ours"][pos],
                                         aa[pos])
        tips6[pos], tips1[pos] = np.asarray(t6, dtype=float), moved(t1)
        labels[pos] = (label, fb)
    fallbacks = [p for p in anchors if labels[p][1]]
    # The face is drawn round the tips as 6ARU has them, not the shared points.
    face_tips = [np.asarray(common.functional_point(six["by_ours"][q], aa[q])[0],
                            dtype=float) for q in anchors]

    face = sorted(p for p in shared if any(
        np.linalg.norm(a.coord - tip) <= FACE_RADIUS
        for tip in face_tips for a in six["by_ours"][p] if a.element != "H"))
    face_rmsd = float(np.sqrt(np.mean([ca_dev[p] ** 2 for p in face])))

    emit("5. The anchor face")
    emit()
    emit(f"   {len(face)} residues have an atom within {FACE_RADIUS:.0f} A of an "
         "anchor's charged tip in 6ARU.")
    emit(f"   Their CA RMSD under the whole-fragment fit: {face_rmsd:.2f} A "
         f"(whole fragment {rmsd:.2f} A)")
    emit()

    # fit on the anchors alone, with the usual caveat
    sup8 = Superimposer()
    sup8.set_atoms([six["by_ours"][p]["CA"] for p in anchors],
                   [one["by_ours"][p]["CA"] for p in anchors])
    emit(f"   Fitting on the {len(anchors)} anchors' CA atoms alone gives "
         f"{float(sup8.rms):.2f} A. A fit on {len(anchors)} points is flattering "
         "by construction,")
    emit("   so this is shown and not relied on; the whole-fragment fit above is "
         "the one the")
    emit("   decision uses.")
    emit()

    emit("6. Each anchor, in the common frame")
    emit()
    emit("   | anchor | CA moved | atom compared | that atom moved | basis |")
    emit("   |---|---|---|---|---|")
    tip_moves, measured = [], []
    for pos in anchors:
        d = float(np.linalg.norm(tips1[pos] - tips6[pos]))
        label, fb = labels[pos]
        if not fb:
            tip_moves.append(d)
            measured.append(pos)
        lacks = [st["id"] for st in (six, one)
                 if not all(a in st["by_ours"][pos]
                            for a in common.FUNCTIONAL_ATOM[aa[pos]])]
        emit(f"   | {aa[pos]}{pos} | {ca_dev[pos]:.2f} A | {label} | {d:.2f} A | "
             + ("charged tip, resolved in both" if not fb else
                f"the tip is not resolved in {' or '.join(lacks)}, so the nearer "
                f"atom {label} is compared in both") + " |")
    emit()
    if fallbacks:
        emit("   Not measurable at the charged tip: "
             + ", ".join(f"{aa[p]}{p}" for p in fallbacks)
             + ". The atom compared is the same in both structures, so the figure is "
             "a real movement")
        emit("   of that atom. The tip sits further out along the side chain and "
             "usually moves at least")
        emit("   as much, though a side chain can rotate so that it moves less, "
             "and this does not")
        emit("   establish which. These anchors are excluded from the deciding "
             "rule, which is about")
        emit("   tips, and listed as a gap.")
        emit()
    emit(f"   Charged tip measured for {len(measured)} of {len(anchors)} anchors. "
         f"Largest move among them: {max(tip_moves):.2f} A "
         f"({aa[measured[int(np.argmax(tip_moves))]]}"
         f"{measured[int(np.argmax(tip_moves))]}); mean "
         f"{float(np.mean(tip_moves)):.2f} A; tolerance {TIP_TOLERANCE:.1f} A")
    emit()

    # ---- 7. decision ---------------------------------------------------------
    verdict = decide(tip_moves)
    emit("=" * 72)
    emit("7. Decision, and what it changes about the design")
    emit("=" * 72)
    emit()
    if verdict == "same":
        emit("   The two cuts are the same where the anchors' charged tips can be")
        emit(f"   compared, by the rule fixed above: no tip moves by "
             f"{TIP_TOLERANCE:.1f} A or more. The fragment stays")
        emit("   cut from 6ARU. The reason is the measurement and not the argument "
             "that domain III")
        emit("   is a rigid blob.")
        if fallbacks:
            gap = ", ".join(f"{aa[p]}{p} (CA moved {ca_dev[p]:.2f} A, "
                            f"{labels[p][0]} moved {float(np.linalg.norm(tips1[p] - tips6[p])):.2f} A)"
                            for p in fallbacks)
            emit()
            emit(f"   The verdict covers {len(measured)} of {len(anchors)} anchors. "
                 f"Not covered: {gap}.")
            emit("   The charged tip of each is unresolved in one structure, so the "
                 "rule could not be")
            emit("   applied to it. This is a stated gap in the verdict and not a "
                 "finding of no movement.")
        emit()
        emit("   The tethered cut is not needed for the design run. It is on the "
             "shelf at")
        emit("   data/structures/1nql_domain3.pdb in case the design run's own "
             "predictor")
        emit("   disagrees.")
    else:
        far = [f"{aa[p]}{p} ({m:.2f} A)" for p, m in zip(measured, tip_moves)
               if m >= TIP_TOLERANCE]
        emit("   The cuts DIFFER where the anchors are: " + ", ".join(far)
             + f" moved by {TIP_TOLERANCE:.1f} A or more.")
        emit("   By the rule fixed above the tethered cut is the one to design "
             "against.")
        emit("   It is written to data/structures/1nql_domain3.pdb and keeps the "
             "file's own")
        emit("   numbering, which is ours minus 24 as for 6ARU, so the hotspot "
             "numbers in")
        emit("   design/configs/egfr-domain3-h370.json carry over. Step 09's "
             "disulfide and")
        emit("   cut-edge checks have to be rerun on it before it replaces the "
             "6ARU cut.")
    emit()
    emit("   H418 is a separate matter and changes nothing for the H370 face: it "
         "reads as")
    emit("   buried in 6ARU (RSA 0.032) and partly exposed in 1NQL (0.200), and 6ARU "
         "being an")
    emit("   antibody complex is the suspected cause. If the tethered form is what "
         "is tested,")
    emit("   the 1NQL reading is the relevant one. That belongs to the 415-466 "
         "fallback.")
    emit()
    emit("   What this does not settle: it compares two crystal structures "
         "of a molecule that")
    emit("   moves, at 3.2 A and 2.8 A resolution. It says nothing about whether "
         "the cut")
    emit("   fragment keeps its shape once separated, which is the open question "
         "from step 11.")
    emit()

    if write_files:
        DERIVED.mkdir(parents=True, exist_ok=True)
        with (DERIVED / "16-anchor-displacement.csv").open("w") as fh:
            fh.write("uniprot_pos,aa,ca_moved_a,compared_atom,atom_moved_a,"
                     "tip_resolved_in_both\n")
            for pos in anchors:
                m = float(np.linalg.norm(tips1[pos] - tips6[pos]))
                fh.write(f"{pos},{aa[pos]},{ca_dev[pos]:.3f},{labels[pos][0]},"
                         f"{m:.3f},{not labels[pos][1]}\n")
        keep = [one["numbering"].pdb_of(p) for p in wanted
                if one["numbering"].pdb_of(p) in one["by_file"]]
        io = PDBIO()
        io.set_structure(one["model"].get_parent())
        io.save(str(STRUCT / "1nql_domain3.pdb"),
                select=FragmentSelect(one["chain"], keep))
        emit("Wrote data/derived/16-anchor-displacement.csv and "
             "data/structures/1nql_domain3.pdb")
        emit()

    summary = (f"verdict {verdict} for {len(measured)} of {len(anchors)} anchors; "
               f"largest charged-tip move {max(tip_moves):.2f} A; fragment RMSD "
               f"{rmsd:.2f} A")
    return finish(out, failures, emit, write_files, summary)


def common_ranges(positions):
    """'310-312, 340' from [310, 311, 312, 340]."""
    if not positions:
        return "none"
    runs, start, prev = [], positions[0], positions[0]
    for p in positions[1:]:
        if p != prev + 1:
            runs.append((start, prev))
            start = p
        prev = p
    runs.append((start, prev))
    return ", ".join(f"{a}" if a == b else f"{a}-{b}" for a, b in runs)


def finish(out, failures, emit, write_files, summary):
    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED —")
        for failure in failures:
            emit(f"  - {failure}")
    else:
        emit(f"RESULT: PASSED. {summary}")
    emit("=" * 72)
    if write_files and not failures:
        FINDINGS.mkdir(parents=True, exist_ok=True)
        (FINDINGS / "16-tethered-fragment.md").write_text(
            "# Is the fragment cut from the tethered structure different?\n\n"
            "Computed output of `analysis/16_tethered_fragment.py`. "
            "Do not hand-edit.\n\n"
            "The organisers say the screen uses the tethered form of EGFR. Our "
            "fragment is\ncut from 6ARU, the extended form. This cuts the same "
            "residues, 310 to 480, from\n1NQL, the tethered form, lays it on the "
            "6ARU cut, and decides on the measured\nmovement of each anchor's "
            "charged tip.\n\n"
            "```\n" + "\n".join(out) + "\n```\n")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
