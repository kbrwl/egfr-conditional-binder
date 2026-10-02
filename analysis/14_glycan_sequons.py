#!/usr/bin/env python3
"""
14_glycan_sequons.py — which anchors sit near a place where a sugar chain attaches?

Abbreviations, expanded here because each file gets read on its own:
  EGFR     epidermal growth factor receptor, the protein we design against.
  UniProt  the public protein sequence archive. Every residue number in this
           project is a position in its human record P00533.
  PDB      the Protein Data Bank, the public archive of measured 3D structures.
           Also used to mean a file from that archive.
  ECD      extracellular domain, the part of EGFR outside the cell, residues
           25-645 in our numbering.
  HEK293   a human cell line. A protein made in it carries human-like sugar
           chains.
  NAG      N-acetylglucosamine, the first sugar of an N-linked sugar chain.
  ND2      the nitrogen atom at the end of an asparagine's side chain, which is
           the atom a sugar chain is bonded to.
  CB       the first carbon of a side chain, one atom in from the tip.
  N, X, S, T, P   single-letter codes for asparagine, any amino acid, serine,
           threonine and proline.

WHY THIS STEP EXISTS
--------------------
The organisers confirmed on 30 September that the target in the assay is made in
HEK293 cells and is glycosylated (docs/competition-qa-log.md). A glycan is a
branched tree of sugars attached to the protein after it is built. Glycans are
large and they move about, so a binder aimed next to one may find the site covered
in the real molecule while the model shows it open.

Step 03 measured this for the old 415-466 epitope and found N444 inside it, bonded
to a sugar, 1.44 angstroms away. Step 08 then moved the target face to the cluster
around H370 and the measurement was never repeated. Six of the eight current
anchors have never been checked against a place where a sugar chain attaches.

WHAT THIS SCRIPT COMPUTES
-------------------------
1. Every N-X-S/T pattern in the human extracellular sequence, found by reading the
   sequence. This pattern is called a sequon: an asparagine (N), then any amino
   acid except proline (X), then a serine or threonine (S/T). N-linked glycans
   attach to the asparagine of a sequon. So N-A-T is a sequon and N-A-V is not,
   because valine is neither serine nor threonine, and N-P-T is not either,
   because proline in the middle position blocks the attachment even though the
   rest of the pattern fits. The pattern N-X-C also supports attachment at a much
   lower rate and is not counted here, which makes the count a lower bound.

2. The same sequons against the mouse sequence. A sequon present in one species
   and absent in the other would put a sugar chain on one target and not the other,
   which breaks mouse cross-reactivity at the one place we cannot afford it. The
   two sequences are aligned first, because a single inserted residue shifts every
   position after it.

3. For each anchor, the distance to the attachment nitrogen (ND2) of every
   sequon, in two measured structures:
     6ARU  EGFR in its extended shape, the structure of record
     1NQL  EGFR in its tethered shape, which the organisers say is the one tested
   The distance is to the attachment nitrogen and not to the sugars that happen to
   be visible. A crystal structure resolves only the first sugar or two, because
   only those hold still, so the real reach of a chain is longer than what the file
   shows. Step 03's 5 angstrom test against visible sugars therefore cannot settle
   the question, and says so.

   Distances use the bands step 05 set, from egfr_common: under 15 A the anchor is
   likely shadowed at least some of the time, and under 25 A it is within reach of
   an extended chain. UNVERIFIED: the 20-30 A reach those bands rest on is from
   memory and is not measured here. They are cautious on purpose, and 'within
   reach' means covered some of the time rather than blocked.

   Whether a sequon is actually occupied is not knowable from the sequence. A
   sequon is a necessary condition for attachment and a sugar seen bonded to it in
   a structure is evidence of occupancy. The absence of a visible sugar is weak
   evidence that nothing is there, for the reason above.

WHAT AN ANCHOR NEAR A SEQUON MEANS FOR THE DESIGN
-------------------------------------------------
Nothing is dropped here. An anchor near a glycan is reported with its distance, and
the ranking in step 10 is the place to demote a design that leans on one, in the
same way step 10 already demotes a design that leans on E424 near the cut edge.
Every decision about what to submit is made in one place.

HOW THIS STEP IS CHECKED
------------------------
  - The old 415-466 anchors are run through the same code as a control. Their
    distances to the N444 CB are compared against the table step 05 committed, and
    the distance from N444 to its sugar against the figure step 03 committed. The
    run stops if either differs.
  - The sequon rule is tested on constructed sequences and each branch of it is
    broken on purpose: `--break-rule proline`, `hydroxyl`, `overlap` and
    `mouse-compare`. The matching test must then fail, and `--self-test` runs all
    four.
  - Every asparagine this script finds is confirmed to be an asparagine in each
    structure, by residue identity and not by arithmetic. `--break-numbering N`
    shifts the numbering and should fail.

Writes:
  results/findings/14-glycan-sequons.md
  data/derived/14-sequons.csv
  data/derived/14-anchor-glycan-distance.csv

Run standalone:  python analysis/14_glycan_sequons.py
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import egfr_common as common  # noqa: E402

STRUCT = common.STRUCT_DIR
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS

STRUCTURE_FILES = {"6ARU": "6aru.pdb", "1NQL": "1nql.pdb"}
SHAPE = {"6ARU": "extended", "1NQL": "tethered"}

ECD = (common.ECD_START, common.ECD_END)

# A sugar bonded to an asparagine sits about 1.45 angstroms from its nitrogen. 2.0
# leaves room for the resolution of the structure without admitting a sugar that
# merely sits nearby, which would be 3 angstroms or more.
BOND_MAX = 2.0

# Set by --break-rule, which disables one branch on purpose so the tests can be
# confirmed to catch it. A check that has never been seen to fail is not known to
# work, and that applies to a test case as much as to an assertion.
BROKEN_RULE = None

# Each branch of the sequon rule, the test that is supposed to catch its removal,
# and what goes wrong if nothing does.
BREAKABLE_RULES = {
    "proline": (
        "proline blocks",
        "an N-P-T pattern would be counted as a place a sugar attaches, although "
        "proline in the middle position prevents it"),
    "hydroxyl": (
        "serine or threonine",
        "an N-A-V pattern would be counted, although valine cannot support "
        "attachment, so every asparagine would look like a sequon"),
    "overlap": (
        "overlapping sequons",
        "two sequons sharing residues, as in N-N-S-T, would be counted as one, "
        "and an attachment point would go unmeasured"),
    "mouse-compare": (
        "mouse comparison",
        "a sequon present in one species and absent in the other would be reported "
        "as shared, hiding exactly the difference that would break cross-reactivity"),
}


def rule_live(name):
    return BROKEN_RULE != name


# ---------------------------------------------------------------------------
# Finding sequons, and comparing them between species
# ---------------------------------------------------------------------------

def find_sequons(sequence, first=1):
    """Every N-X-S/T in `sequence`, as (position of the N, the three residues).

    `first` is the number of the first residue, so a slice of a longer sequence
    keeps its own numbering. Sequons are allowed to overlap, as in N-N-S-T where
    both the first and the second asparagine are followed by a pattern that fits.
    """
    found = []
    skip_until = 0
    for i in range(len(sequence) - 2):
        if sequence[i] != "N":
            continue
        if not rule_live("overlap") and i < skip_until:
            continue
        if sequence[i + 1] == "P" and rule_live("proline"):
            continue
        if sequence[i + 2] not in "ST" and rule_live("hydroxyl"):
            continue
        found.append((first + i, sequence[i:i + 3]))
        skip_until = i + 3
    return found


def column_map(human, mouse):
    """Our position in the human sequence -> the aligned position in the mouse one,
    for columns where both have a residue."""
    alignment = common.make_aligner().align(human, mouse)[0]
    human_i = mouse_i = 0
    mapping = {}
    for h_char, m_char in zip(alignment[0], alignment[1]):
        if h_char != "-":
            human_i += 1
        if m_char != "-":
            mouse_i += 1
        if h_char != "-" and m_char != "-":
            mapping[human_i] = mouse_i
    return mapping


def compare_species(human, mouse, lo, hi):
    """Sequons in the human range lo-hi against the aligned mouse sequence.

    Returns rows of {pos, human, mouse, status} where `pos` is a human position,
    `human` and `mouse` are the three residues at that place in each species, and
    `status` is one of 'both', 'human only' or 'mouse only'. A mouse sequon whose
    asparagine has no aligned human residue is reported 'mouse only' with a
    position of None rather than dropped.
    """
    h2m = column_map(human, mouse)
    m2h = {m: h for h, m in h2m.items()}
    human_sequons = {p: t for p, t in find_sequons(human[lo - 1:hi], first=lo)}

    mapped = [h2m[p] for p in range(lo, hi + 1) if p in h2m]
    mouse_lo, mouse_hi = (min(mapped), max(mapped)) if mapped else (1, len(mouse))
    mouse_sequons = {p: t for p, t in
                     find_sequons(mouse[mouse_lo - 1:mouse_hi], first=mouse_lo)}

    def mouse_triplet(h_pos):
        letters = []
        for k in range(3):
            m = h2m.get(h_pos + k)
            letters.append(mouse[m - 1] if m else "-")
        return "".join(letters)

    def human_triplet(h_pos):
        return human[h_pos - 1:h_pos + 2]

    rows = []
    for pos, tri in sorted(human_sequons.items()):
        m_pos = h2m.get(pos)
        shared = m_pos in mouse_sequons and rule_live("mouse-compare")
        if not rule_live("mouse-compare"):
            shared = True
        rows.append(dict(pos=pos, human=tri, mouse=mouse_triplet(pos),
                         status="both" if shared else "human only"))
    for m_pos, tri in sorted(mouse_sequons.items()):
        h_pos = m2h.get(m_pos)
        if h_pos in human_sequons or not rule_live("mouse-compare"):
            continue
        rows.append(dict(pos=h_pos,
                         human=human_triplet(h_pos) if h_pos else "---",
                         mouse=tri, status="mouse only"))
    return rows


# ---------------------------------------------------------------------------
# Reading a structure
# ---------------------------------------------------------------------------

def sugar_atoms(model):
    """Every atom of every sugar residue in the model, in any chain."""
    return [a for chain in model for r in chain
            if r.get_resname().strip() in common.GLYCAN_NAMES
            for a in r if a.element != "H"]


def read_structure(pid, human, shift, emit):
    """Everything this step needs from one structure file.

    Returns a dict with the numbering, the receptor residues keyed by file number,
    and the sugar atoms in the file.
    """
    path = STRUCT / STRUCTURE_FILES[pid]
    if not path.exists():
        raise SystemExit(f"{path} is missing. Run analysis/06_tethered_occlusion.py, "
                         "which downloads both structures.")
    model = common.load_complex(path)
    chain_id, numbering, pct, others = common.receptor_chain_of(model, human)
    offset, share, _ = numbering.dominant_offset()
    emit(f"   {pid} ({SHAPE[pid]}): receptor chain {chain_id}, "
         f"{len(numbering)} residues mapped, {pct:.1f}% identity to human EGFR")
    emit(f"        numbering offset: ours = file {offset:+d} "
         f"({100.0 * share:.1f}% of residues), derived and not assumed")
    sugars = sugar_atoms(model)
    emit(f"        sugar atoms in the file: {len(sugars)}; other protein chains "
         f"left out of the measurement: {', '.join(others) or 'none'}")
    if shift:
        moved = {numbering.pdb_of(u): u + shift for u in numbering.uniprot_positions()
                 if numbering.pdb_of(u) is not None}
        numbering = common.Numbering(moved)
    by_file = {r.id[1]: r for r in common.protein_residues(model[chain_id])}
    return dict(id=pid, model=model, numbering=numbering, residues=by_file,
                sugars=sugars)


def residue_at(structure, uniprot_pos):
    file_num = structure["numbering"].pdb_of(uniprot_pos)
    return structure["residues"].get(file_num) if file_num is not None else None


def check_identity(structure, expected, label, emit):
    """Every position we believe is a given amino acid reads as that amino acid.

    By residue identity, because arithmetic can be off by one and still look
    right, and a histidine that turns out to be a leucine cannot. Positions the
    structure does not contain are reported as such and are not a failure.
    """
    problems, unresolved, checked = [], [], 0
    for pos, aa in sorted(expected.items()):
        res = residue_at(structure, pos)
        if res is None:
            unresolved.append(pos)
            continue
        checked += 1
        found = common.THREE_TO_ONE.get(res.get_resname())
        if found != aa:
            problems.append(f"{aa}{pos} reads as {found}")
    if problems:
        emit(f"   [FAIL] identity check, {label}, {structure['id']}: "
             + "; ".join(problems))
        return False, unresolved
    emit(f"   [PASS] identity check, {label}, {structure['id']}: all {checked} "
         "positions hold the amino acid our numbering says"
         + (f"; {len(unresolved)} not present in the file" if unresolved else ""))
    return True, unresolved


def attachment_evidence(structure, asn):
    """Is a sugar bonded to this asparagine in the structure?

    Returns (distance from the attachment nitrogen to the nearest sugar atom,
    whether that is close enough to be a bond). (None, False) when the nitrogen
    or any sugar is missing.
    """
    if asn is None or "ND2" not in asn or not structure["sugars"]:
        return None, False
    nitrogen = asn["ND2"].coord
    best = min(float(np.linalg.norm(a.coord - nitrogen)) for a in structure["sugars"])
    return best, best <= BOND_MAX


def residue_to_sugar_min(structure, asn):
    """Closest approach of any atom of the residue to any sugar atom, which is how
    step 03 measured N444 and is the number the control compares."""
    if asn is None or not structure["sugars"]:
        return None
    return min(float(np.linalg.norm(a.coord - s.coord))
               for a in asn if a.element != "H" for s in structure["sugars"])


def band(distance):
    if distance < common.GLYCAN_NEAR:
        return "likely shadowed"
    if distance < common.GLYCAN_PLAUSIBLE:
        return "within reach of an extended chain"
    return "probably clear"


# ---------------------------------------------------------------------------
# The tests, run every time
# ---------------------------------------------------------------------------

def run_tests(emit):
    """Constructed sequences, each exercising one branch, with the answer written
    out rather than computed."""
    failures = []

    emit("   Sequon rule, on constructed sequences:")
    emit()
    cases = [
        ("KNATQ", [(2, "NAT")], "serine or threonine / a plain sequon"),
        ("KNAVQ", [], "serine or threonine: valine is not S or T"),
        ("KNPTQ", [], "proline blocks attachment in the middle position"),
        ("KNASQ", [(2, "NAS")], "serine or threonine / serine accepted"),
        ("KNNSTQ", [(2, "NNS"), (3, "NST")], "overlapping sequons are both counted"),
        ("NAT", [(1, "NAT")], "a sequon at the very start"),
        ("QQNA", [], "an asparagine too near the end to form a sequon"),
    ]
    for sequence, want, why in cases:
        got = find_sequons(sequence)
        ok = got == want
        emit(f"     [{'PASS' if ok else 'FAIL'}] {sequence:7s} -> {got}   ({why})")
        if not ok:
            label = {"KNPTQ": "proline blocks", "KNAVQ": "serine or threonine",
                     "KNNSTQ": "overlapping sequons"}.get(sequence, "sequon rule")
            failures.append(f"{label}: {sequence} gave {got}, expected {want}")
    emit()

    emit("   Mouse comparison, on constructed sequences:")
    emit()
    human = "AAKNATQQNKSAA"
    # position:   1 2 3 4 5 6 7 8 9 10 11 12 13
    #             A A K N A T Q Q N K  S  A  A     sequons: 4 (NAT) and 9 (NKS)
    lost = "AAKNAAQQNKSAA"          # T6 -> A: the sequon at 4 is lost in mouse
    gained_h = "AAKNAAQQNKSAA"      # human without the sequon at 4
    gained_m = "AAKNATQQNKSAA"      # mouse with it
    tests = [
        ("sequon lost in mouse", human, lost,
         {4: "human only", 9: "both"}),
        ("sequon gained in mouse", gained_h, gained_m,
         {4: "mouse only", 9: "both"}),
        ("identical sequences", human, human, {4: "both", 9: "both"}),
    ]
    for name, h, m, want in tests:
        got = {r["pos"]: r["status"] for r in compare_species(h, m, 1, len(h))}
        ok = got == want
        emit(f"     [{'PASS' if ok else 'FAIL'}] {name}: {got}")
        if not ok:
            failures.append(f"mouse comparison: {name} gave {got}, expected {want}")
    emit()
    return failures


# ---------------------------------------------------------------------------
# Controls against the earlier steps
# ---------------------------------------------------------------------------

def read_step03_n444():
    """The distance step 03 committed for N444 to its sugar, read from its
    findings file. Raises if the line is not there, so a reworded findings file
    cannot make this control silently pass."""
    text = (FINDINGS / "03-solvent-accessibility.md").read_text()
    match = re.search(r"N444:\s*([0-9.]+)\s*A from NAG", text)
    if not match:
        raise common.CrossCheckError(
            "could not find the N444 distance in 03-solvent-accessibility.md")
    return float(match.group(1))


def read_step05_table():
    """Step 05's committed table of anchor distances to the N444 CB, as
    {anchor position: distance}."""
    text = (FINDINGS / "05-anchor-geometry.md").read_text()
    if "distance to N444 CB" not in text:
        raise common.CrossCheckError(
            "could not find the N444 table in 05-anchor-geometry.md")
    section = text[text.index("distance to N444 CB"):]
    rows = re.findall(r"\|\s*[DEH](\d+)\s*\|\s*([0-9.]+) A\s*\|", section)
    table = {}
    for pos, dist in rows:
        table.setdefault(int(pos), float(dist))
    if not table:
        raise common.CrossCheckError("the N444 table in step 05 has no rows")
    return table


# ---------------------------------------------------------------------------

def main(argv=None):
    global BROKEN_RULE

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--break-rule", choices=list(BREAKABLE_RULES), default=None,
                        help="switch one branch off, to confirm the matching test "
                             "then fails. Expected to fail.")
    parser.add_argument("--break-numbering", type=int, default=0, metavar="N",
                        help="shift every residue number by N, to confirm the "
                             "identity checks notice. Expected to fail.")
    parser.add_argument("--self-test", action="store_true",
                        help="switch each branch off in turn and confirm a test "
                             "notices.")
    args = parser.parse_args(argv)
    if args.self_test:
        return run_self_test()
    BROKEN_RULE = args.break_rule

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    human = common.human_sequence()
    mouse = common.mouse_sequence()
    current = common.load_primary_cluster()
    old = sorted(common.ANCHORS)
    anchor_aa = {p: human[p - 1] for p in set(current) | set(old)}

    emit("=" * 72)
    emit("GLYCAN SEQUONS — which anchors sit near a place a sugar chain attaches?")
    emit("=" * 72)
    emit()
    emit("The target in the assay is made in human cells and is glycosylated. A")
    emit("sugar chain is large and mobile, so a binder aimed beside one may find the")
    emit("site covered in the real molecule while the model shows it open. Step 03")
    emit("checked this for the old 415-466 epitope and found N444 inside it, bonded")
    emit("to a sugar. It was never repeated for the current face around H370.")
    emit()
    if BROKEN_RULE:
        emit(f"   NOTE: --break-rule {BROKEN_RULE} is in effect: "
             f"{BREAKABLE_RULES[BROKEN_RULE][1]}.")
        emit("   This run is expected to fail.")
        emit()
    if args.break_numbering:
        emit(f"   NOTE: --break-numbering {args.break_numbering} is in effect. "
             "This run is")
        emit("   expected to fail at the identity checks.")
        emit()

    failures = []

    # ---- 1. The tests -------------------------------------------------------
    emit("1. The sequon rule, tested every time")
    emit()
    failures.extend(run_tests(emit))

    # ---- 2. Sequons in the human sequence, and against mouse ---------------
    lo, hi = ECD
    sequons = compare_species(human, mouse, lo, hi)
    human_sequons = [r for r in sequons if r["status"] != "mouse only"]
    emit(f"2. Sequons in the human extracellular sequence, residues {lo}-{hi}")
    emit()
    emit(f"   {len(human_sequons)} found by reading the sequence. Compared against "
         "the mouse sequence")
    emit("   after aligning the two:")
    emit()
    emit("   | position | human | mouse, aligned | status |")
    emit("   |---|---|---|---|")
    for r in sequons:
        pos = r["pos"] if r["pos"] is not None else "no aligned position"
        emit(f"   | N{pos} | {r['human']} | {r['mouse']} | {r['status']} |")
    emit()
    different = [r for r in sequons if r["status"] != "both"]
    if different:
        emit(f"   {len(different)} of {len(sequons)} differ between the species: "
             + ", ".join(f"N{r['pos']} ({r['status']})" for r in different))
        emit("   A sequon present in one species and absent in the other puts a sugar")
        emit("   chain on one target and not the other. Where one of these lies near")
        emit("   an anchor it is a cross-reactivity risk at that anchor, and section 4")
        emit("   says which.")
    else:
        emit(f"   All {len(sequons)} are present in both species.")
    emit()

    # ---- 3. Structures ------------------------------------------------------
    emit("3. The two structures")
    emit()
    structures = {}
    for pid in ("6ARU", "1NQL"):
        structures[pid] = read_structure(pid, human, args.break_numbering, emit)
    emit()

    sequon_positions = [r["pos"] for r in human_sequons]
    asn_expected = {p: "N" for p in sequon_positions}
    anchor_expected = dict(anchor_aa)
    emit("   Identity checks, by amino acid and not by arithmetic:")
    for pid, st in structures.items():
        ok, _ = check_identity(st, asn_expected, "every sequon asparagine", emit)
        if not ok:
            failures.append(f"identity check failed against {pid} for the "
                            "sequon asparagines")
        ok, _ = check_identity(st, anchor_expected, "every anchor", emit)
        if not ok:
            failures.append(f"identity check failed against {pid} for the anchors")
    emit()

    # ---- the sugar evidence per sequon ----
    emit("   Is a sugar seen bonded to each sequon in the structures? A sugar within")
    emit(f"   {BOND_MAX} A of the attachment nitrogen is a bond. Its absence is weak")
    emit("   evidence of nothing there, because a structure shows only the sugars that")
    emit("   hold still.")
    emit()
    emit("   | sequon | 6ARU nitrogen to nearest sugar | 1NQL nitrogen to nearest sugar |")
    emit("   |---|---|---|")
    evidence = {}
    for r in human_sequons:
        cells = []
        for pid, st in structures.items():
            asn = residue_at(st, r["pos"])
            d, bonded = attachment_evidence(st, asn)
            evidence[(pid, r["pos"])] = (d, bonded)
            if asn is None:
                cells.append("not in the file")
            elif d is None:
                cells.append("nitrogen not resolved")
            else:
                cells.append(f"{d:.2f} A" + (" — bonded" if bonded else ""))
        emit(f"   | N{r['pos']} | {cells[0]} | {cells[1]} |")
    emit()

    # ---- 4. Anchor distances ------------------------------------------------
    def distances(structure, anchors):
        """For each anchor, distances to every measurable sequon's nitrogen and CB."""
        result = {}
        for a in anchors:
            res = residue_at(structure, a)
            if res is None:
                result[a] = None
                continue
            point, label, fallback = common.functional_point(res, anchor_aa[a])
            rows = []
            for s in sequon_positions:
                asn = residue_at(structure, s)
                if asn is None or "ND2" not in asn:
                    continue
                d_nd2 = float(np.linalg.norm(np.asarray(point) - asn["ND2"].coord))
                d_cb = (float(np.linalg.norm(np.asarray(point) - asn["CB"].coord))
                        if "CB" in asn else None)
                rows.append((s, d_nd2, d_cb))
            result[a] = dict(point_label=label, fallback=fallback, rows=rows)
        return result

    emit("4. Distance from each anchor to every sequon's attachment nitrogen")
    emit()
    emit("   Measured from the charged tip of each anchor's side chain, the same")
    emit("   point step 05 uses, to the ND2 of every sequon. Bands, from step 05:")
    emit(f"     under {common.GLYCAN_NEAR:.0f} A   likely shadowed at least some of "
         "the time")
    emit(f"     under {common.GLYCAN_PLAUSIBLE:.0f} A   within reach of an extended "
         "chain")
    emit("     beyond      probably clear")
    emit()

    summary_rows = []
    dist_by = {}
    for set_name, anchors in (("current H370 face", current),
                              ("old 415-466 epitope", old)):
        emit(f"   {set_name}")
        emit()
        emit("   | anchor | nearest sequon, 6ARU | nearest sequon, 1NQL | "
             "closest overall | assessment |")
        emit("   |---|---|---|---|---|")
        for a in anchors:
            cells, best = [], None
            for pid, st in structures.items():
                if (pid, set_name) not in dist_by:
                    dist_by[(pid, set_name)] = distances(st, anchors)
                info = dist_by[(pid, set_name)][a]
                if info is None or not info["rows"]:
                    cells.append("not measurable")
                    continue
                s, d, _cb = min(info["rows"], key=lambda t: t[1])
                cells.append(f"N{s} at {d:.1f} A")
                if best is None or d < best[0]:
                    best = (d, s, pid)
                summary_rows.append((set_name, pid, a, s, d, band(d),
                                     sum(1 for t in info["rows"]
                                         if t[1] < common.GLYCAN_PLAUSIBLE)))
            if best:
                emit(f"   | {anchor_aa[a]}{a} | {cells[0]} | {cells[1]} | "
                     f"N{best[1]} at {best[0]:.1f} A ({best[2]}) | {band(best[0])} |")
            else:
                emit(f"   | {anchor_aa[a]}{a} | {cells[0]} | {cells[1]} | — | — |")
        emit()

    # ---- which anchors, in one line each ----
    flagged = {}
    for set_name, pid, a, s, d, b, n_within in summary_rows:
        if set_name != "current H370 face":
            continue
        if d < common.GLYCAN_PLAUSIBLE:
            cur = flagged.get(a)
            if cur is None or d < cur[0]:
                flagged[a] = (d, s, pid, b)
    emit("   The current anchors, one line each:")
    emit()
    for a in current:
        if a in flagged:
            d, s, pid, b = flagged[a]
            species = next((r["status"] for r in sequons if r["pos"] == s), "?")
            bonded = [p for p in structures if evidence.get((p, s), (None, False))[1]]
            emit(f"     {anchor_aa[a]}{a}: {b}. Nearest attachment point N{s}, "
                 f"{d:.1f} A away in {pid}; sequon {species} across species"
                 + (f"; a sugar is seen bonded there in {', '.join(bonded)}"
                    if bonded else "; no sugar seen bonded there"))
        else:
            emit(f"     {anchor_aa[a]}{a}: probably clear. No sequon within "
                 f"{common.GLYCAN_PLAUSIBLE:.0f} A in either structure.")
    emit()

    # ---- 5. Controls --------------------------------------------------------
    emit("5. Controls against steps 03 and 05")
    emit()
    emit("   The old 415-466 anchors go through this same code. If it cannot")
    emit("   reproduce what those steps committed, the answer above is not to be")
    emit("   believed either.")
    emit()
    six = structures["6ARU"]

    n444 = residue_at(six, 444)
    emit(f"   N444 in the human sequence: {human[443:446]} — "
         + ("a sequon" if 444 in sequon_positions else "NOT found as a sequon"))
    if 444 not in sequon_positions:
        failures.append("control: N444 was not found as a sequon in the human "
                        "sequence")
    committed = read_step03_n444()
    mine = residue_to_sugar_min(six, n444)
    if mine is None:
        failures.append("control: N444 or the sugars could not be read in 6ARU")
    else:
        ok = abs(round(mine, 2) - committed) < 0.006
        emit(f"   [{'PASS' if ok else 'FAIL'}] N444 to its nearest sugar atom in "
             f"6ARU: {mine:.2f} A here, {committed:.2f} A in step 03")
        if not ok:
            failures.append(f"control: N444 to sugar is {mine:.2f} A here and "
                            f"{committed:.2f} A in step 03")

    table = read_step05_table()
    info6 = dist_by.get(("6ARU", "old 415-466 epitope"), {})
    mismatches, compared = [], 0
    for a, committed_d in sorted(table.items()):
        info = info6.get(a)
        row = next((t for t in (info or {}).get("rows", []) if t[0] == 444), None)
        if row is None or row[2] is None:
            mismatches.append(f"{anchor_aa.get(a, '?')}{a}: not measurable here")
            continue
        compared += 1
        if abs(round(row[2], 1) - committed_d) > 0.051:
            mismatches.append(f"{anchor_aa[a]}{a}: {row[2]:.1f} A here, "
                              f"{committed_d:.1f} A in step 05")
    if mismatches:
        emit("   [FAIL] anchors to the N444 CB against step 05's table: "
             + "; ".join(mismatches))
        failures.append("control: " + "; ".join(mismatches))
    else:
        emit(f"   [PASS] anchors to the N444 CB against step 05's table: all "
             f"{compared} distances agree to 0.1 A")
    emit()

    # ---- 6. What it changes -------------------------------------------------
    emit("=" * 72)
    emit("6. What this changes about the design")
    emit("=" * 72)
    emit()
    near = [a for a in current if a in flagged and flagged[a][0] < common.GLYCAN_NEAR]
    plausible = [a for a in current if a in flagged
                 and flagged[a][0] >= common.GLYCAN_NEAR]
    clear = [a for a in current if a not in flagged]
    emit(f"   Likely shadowed (under {common.GLYCAN_NEAR:.0f} A): "
         + (", ".join(f"{anchor_aa[a]}{a}" for a in near) or "none"))
    emit(f"   Within reach of an extended chain (under {common.GLYCAN_PLAUSIBLE:.0f}"
         " A): "
         + (", ".join(f"{anchor_aa[a]}{a}" for a in plausible) or "none"))
    emit("   Probably clear: "
         + (", ".join(f"{anchor_aa[a]}{a}" for a in clear) or "none"))
    emit()
    emit("   Nothing is dropped. These are reported so that step 10's ranking can")
    emit("   demote a design leaning on a flagged anchor, the way it already demotes")
    emit("   one leaning on E424 near the cut edge. A crystal structure cannot")
    emit("   settle whether a mobile chain actually covers a site, so 'within reach'")
    emit("   means covered some of the time and not blocked.")
    emit()

    # ---- files --------------------------------------------------------------
    if BROKEN_RULE is None and not args.break_numbering:
        DERIVED.mkdir(parents=True, exist_ok=True)
        with (DERIVED / "14-sequons.csv").open("w") as fh:
            fh.write("uniprot_pos,human_residues,mouse_residues_aligned,"
                     "species_status,6ARU_nitrogen_to_sugar_a,6ARU_bonded,"
                     "1NQL_nitrogen_to_sugar_a,1NQL_bonded\n")
            for r in sequons:
                cells = []
                for pid in ("6ARU", "1NQL"):
                    d, bonded = evidence.get((pid, r["pos"]), (None, False))
                    cells += ["" if d is None else f"{d:.3f}", str(bonded)]
                fh.write(f"{r['pos']},{r['human']},{r['mouse']},{r['status']},"
                         + ",".join(cells) + "\n")
        with (DERIVED / "14-anchor-glycan-distance.csv").open("w") as fh:
            fh.write("anchor_set,structure,uniprot_pos,aa,nearest_sequon,"
                     "nitrogen_distance_a,band,sequons_within_25a\n")
            for set_name, pid, a, s, d, b, n_within in summary_rows:
                fh.write(f"{set_name},{pid},{a},{anchor_aa[a]},{s},{d:.3f},{b},"
                         f"{n_within}\n")
        emit("Wrote data/derived/14-sequons.csv and "
             "data/derived/14-anchor-glycan-distance.csv")
        emit()

    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED —")
        for failure in failures:
            emit(f"  - {failure}")
    else:
        emit(f"RESULT: PASSED. {len(human_sequons)} sequons in the extracellular "
             f"region; {len(near)} anchor(s) likely shadowed, {len(plausible)} "
             f"within reach of an extended chain, {len(clear)} probably clear.")
    emit("=" * 72)

    if BROKEN_RULE is None and not args.break_numbering:
        FINDINGS.mkdir(parents=True, exist_ok=True)
        (FINDINGS / "14-glycan-sequons.md").write_text(
            "# Glycan sequons near the H370 anchors\n\n"
            "Computed output of `analysis/14_glycan_sequons.py`. Do not hand-edit.\n\n"
            "The target in the assay is glycosylated. Step 03 found N444 bonded to a\n"
            "sugar inside the old 415-466 epitope, and the measurement was never\n"
            "repeated after step 08 moved the face to the cluster around H370. This\n"
            "finds every sequon, N-X-S/T, by reading the sequence, compares them with\n"
            "mouse, and measures each anchor's distance to each attachment nitrogen in\n"
            "the extended and the tethered structure.\n\n"
            "```\n" + "\n".join(out) + "\n```\n")
    return 1 if failures else 0


def run_self_test():
    """Remove each branch in turn and confirm a test notices."""
    print("=" * 72)
    print("SELF-TEST: is each branch of the sequon rule covered by a test?")
    print("=" * 72)
    print()
    failures = []
    for rule, (expected_case, consequence) in BREAKABLE_RULES.items():
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--break-rule", rule],
            capture_output=True, text=True,
            cwd=str(Path(__file__).resolve().parent))
        caught = [ln.strip() for ln in proc.stdout.splitlines()
                  if ln.strip().startswith("- ")]
        named = any(expected_case in ln for ln in caught)
        ok = proc.returncode != 0 and named
        print(f"   [{'PASS' if ok else 'FAIL'}] {rule}")
        print(f"     if unnoticed: {consequence}")
        print(f"     exit code {proc.returncode}, {len(caught)} assertion(s) failed")
        for ln in caught[:3]:
            print(f"       {ln}")
        if not ok:
            failures.append(f"{rule}: "
                            + ("removing this branch broke no test"
                               if proc.returncode == 0
                               else f"tests failed but not the one meant to cover "
                                    f"it ({expected_case})"))
        print()
    print("=" * 72)
    if failures:
        print(f"SELF-TEST FAILED: {len(failures)} branch(es) not covered.")
        for failure in failures:
            print(f"  - {failure}")
    else:
        print(f"SELF-TEST PASSED. All {len(BREAKABLE_RULES)} branches are covered by "
              "a test that")
        print("fails when the branch is removed.")
    print("=" * 72)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
