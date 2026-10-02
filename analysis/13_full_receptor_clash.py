#!/usr/bin/env python3
"""
13_full_receptor_clash.py — would each designed binder still fit if the rest of
EGFR were there?

Abbreviations, expanded here because each file gets read on its own:
  PDB      Protein Data Bank, the public archive of measured 3D protein
           structures. Also used to mean a file from that archive.
  mmCIF    macromolecular crystallographic information file, the newer standard
           format for protein structures, and the one the design run writes.
  UniProt  the public protein sequence archive. Our residue numbers are
           positions in its records.
  ECD      extracellular domain, the part of EGFR outside the cell, residues
           25-645 in our numbering.
  CA       the alpha carbon, one atom present in every residue, often used as a
           single stand-in for the whole residue's position.
  RMSD     root-mean-square deviation: the average distance left between matched
           atoms once two structures have been laid on top of each other. A small
           number means the same shape in the same place.
  aa       amino acids, the building blocks a protein chain is made of.

WHY THIS STEP EXISTS
--------------------
The design run is given the 171-residue fragment step 09 cut out, and it never
sees the other 450 residues of the extracellular region. So it cannot tell the
difference between a binder that lies against an open face and a binder that
occupies space another part of the receptor is already in. Both score the same.
The assay uses the whole extracellular region, so the second kind would measure as
nothing, and nothing in the design run's own numbers would say why.

Step 12 answers this for the target: the eight anchors are not covered by the rest
of the receptor in either measured conformation, so the face is reachable. That is
a property of the target and it is settled. This step answers the remaining half,
which is a property of each individual design: a reachable face does not stop one
particular binder from being too large, or approaching at the wrong angle, and
putting part of itself where domain II or domain IV already is.

So step 12 licenses the campaign and this step screens the candidates.

WHAT THIS SCRIPT COMPUTES
-------------------------
For each candidate complex the design run returns:

  1. Identify which chain is the target and which is the binder, by measurement
     rather than by chain letter, using the same shared function step 10 uses. The
     design run puts the target on A and the binder on B, which is the opposite of
     BindCraft version 1, so trusting the letter would read the wrong molecule.

  2. Lay the candidate's target on top of the same residues of the full 6ARU
     structure, and apply that same movement to the binder. This puts the designed
     binder into the coordinate frame of the intact receptor. Done on CA atoms of
     the residues present in both.

  3. Refuse the candidate if that superposition is poor. A high RMSD means the
     movement is not a reliable way to place the binder, so every distance
     computed from it would be unfounded. Reporting a confident verdict from a bad
     superposition is the failure this check exists to prevent.

  4. Count how many binder atoms come within reach of receptor atoms from OUTSIDE
     domain III — the part the design run could not see. Contacts with domain III
     itself are counted too, but only as a sanity check that the binder is where it
     is supposed to be: those are the intended interface and are not a fault.

A candidate is called clashing, marginal or clear. Clashing candidates are not
dropped from the shortlist by this script; they are flagged with the overlap
measured, because a decision to discard a design belongs in one place and that
place is step 10's ranking.

WHY A RIGID SUPERPOSITION IS ENOUGH
-----------------------------------
This asks whether the binder's volume is already occupied, which is a question
about where things are, not about what they would do. It needs no structure
predictor and no GPU (graphics processing unit, the hardware a predictor needs), so
it runs on every candidate in seconds rather than competing for the rented machine.

The alternative was to predict each candidate against the whole extracellular
region. That would answer a richer question, including whether the receptor would
shift to accommodate the binder, and it costs much more per candidate: the target
goes from 171 residues to 621. UNVERIFIED, from memory and not measured here: a
structure predictor's cost grows with roughly the square of the chain length, which
would make that about thirteen times the work per candidate. The ratio is item 9 in
the decisions log's next actions, to be measured on the first run. Either way the
richer check is worth doing for a final shortlist rather than for screening, and
this step is the screen.

What a rigid superposition cannot see: the receptor could move to accommodate a
binder that overlaps it slightly, so a marginal verdict is genuinely marginal
rather than a soft failure. A hard overlap is a different matter, because two sets
of atoms cannot be in the same place.

HOW THIS STEP IS CHECKED
------------------------
No design run has returned anything yet, so every number this script has produced
is from constructed inputs. It builds its own test complexes, in the format the
design run writes, and runs them on every invocation:

  - a binder placed out from the anchor face, which must come back clear
  - a binder placed inside the volume domain IV occupies, which must come back
    clashing
  - a target whose coordinates do not match 6ARU, which must be refused rather
    than scored

`--break-rule` switches off each rule in turn, and the matching test must then
fail. A check that has never been seen to fail is not known to work.

Writes:
  results/findings/13-full-receptor-clash.md
  data/derived/13-candidate-clashes.csv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFIO, Superimposer
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

sys.path.insert(0, str(Path(__file__).resolve().parent))

import egfr_common as common  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
STRUCT = common.STRUCT_DIR
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS
CANDIDATES = ROOT / "results" / "candidates"

FULL_RECEPTOR = STRUCT / "6aru.pdb"
TRIMMED_TARGET = STRUCT / "6aru_domain3.pdb"

DOMAIN_III = (common.D3_START, common.D3_END)

# Two atoms this close cannot both be there: it is an overlap of the atoms
# themselves, not a close approach. Heavy atoms in contact sit around 3.5-4.0 A
# apart, and a carbon atom is about 1.7 A in radius, so below 2.5 A the two are
# interpenetrating.
CLASH_HARD = 2.5

# Touching, which may be tolerable or may mean the binder is straddling a groove.
# Same value step 10 uses for a contact, so the two steps mean the same thing by
# the word.
CLASH_CONTACT = common.CONTACT_CUTOFF

# Above this, the superposition is not a reliable way to place the binder, so no
# distance computed from it can be believed. Domain III is one rigid blob in both
# structures, so a good match should be well under 2 A; 2.5 leaves room for the
# design run having nudged the target without having moved it somewhere else.
MAX_SUPERPOSITION_RMSD = 2.5

# Fewest matched CA atoms for a superposition to mean anything.
MIN_SUPERPOSITION_ATOMS = 30

VERDICT_CLEAR = "clear"
VERDICT_MARGINAL = "marginal"
VERDICT_CLASHING = "clashing"
VERDICT_REFUSED = "not scored"

BREAKABLE_RULES = {
    "hard-clash": "stop counting atom overlaps below "
                  f"{CLASH_HARD} A, so a binder buried in domain IV reads clear",
    "superposition": "accept any superposition however poor, so a target that is "
                     "not 6ARU gets scored anyway",
}
BROKEN_RULE = None


def rule_live(name):
    return BROKEN_RULE != name


# ---------------------------------------------------------------------------
# Placing a candidate into the intact receptor
# ---------------------------------------------------------------------------

def receptor_context():
    """The full receptor chain, split into the part the design run saw and the
    part it did not.

    Returns (model, chain id, numbering, {our pos: residue} inside domain III,
    [residues outside domain III]).
    """
    if not FULL_RECEPTOR.exists():
        raise SystemExit(
            f"{FULL_RECEPTOR} is missing. Run analysis/06_tethered_occlusion.py, "
            "which downloads it.")
    model = common.load_complex(FULL_RECEPTOR)
    human = common.human_sequence()
    chain_id, numbering, pct, others = common.receptor_chain_of(model, human)
    inside, outside = {}, []
    for residue in common.protein_residues(model[chain_id]):
        pos = numbering.uniprot_of(residue.id[1])
        if pos is None:
            continue
        if DOMAIN_III[0] <= pos <= DOMAIN_III[1]:
            inside[pos] = residue
        else:
            outside.append(residue)
    return dict(model=model, chain=chain_id, numbering=numbering, identity=pct,
                partners=others, inside=inside, outside=outside)


def superpose_onto_receptor(target_residues, target_numbering, context):
    """Lay a candidate's target on the same residues of the intact receptor.

    Returns (Superimposer, matched atom count, rmsd) or (None, n, None) when there
    is too little to match on.
    """
    fixed, moving = [], []
    for residue in target_residues:
        pos = target_numbering.uniprot_of(residue.id[1])
        if pos is None or pos not in context["inside"]:
            continue
        reference = context["inside"][pos]
        if "CA" not in residue or "CA" not in reference:
            continue
        fixed.append(reference["CA"])
        moving.append(residue["CA"])
    if len(moving) < MIN_SUPERPOSITION_ATOMS:
        return None, len(moving), None
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    return sup, len(moving), float(sup.rms)


def count_overlaps(binder_coords, residues, hard=CLASH_HARD,
                   contact=CLASH_CONTACT):
    """How many binder atoms come within reach of a set of receptor residues.

    Returns (hard overlaps, contacts, closest distance, the residue that is
    closest). Counted per binder atom rather than per pair, so one atom wedged
    against several residues counts once.
    """
    atoms = [a for r in residues for a in r if a.element != "H"]
    if not atoms or len(binder_coords) == 0:
        return 0, 0, float("inf"), None
    coords = np.array([a.coord for a in atoms])
    hard_n = contact_n = 0
    closest, closest_residue = float("inf"), None
    for point in binder_coords:
        distances = np.linalg.norm(coords - point, axis=1)
        index = int(np.argmin(distances))
        best = float(distances[index])
        if best < closest:
            closest = best
            closest_residue = atoms[index].get_parent()
        if best < contact:
            contact_n += 1
        if best < hard and rule_live("hard-clash"):
            hard_n += 1
    return hard_n, contact_n, closest, closest_residue


def verdict_for(hard, contacts):
    if hard > 0:
        return VERDICT_CLASHING
    if contacts > 0:
        return VERDICT_MARGINAL
    return VERDICT_CLEAR


def analyse_candidate(name, path, context, human, emit):
    """Place one candidate into the intact receptor and measure the overlap."""
    numbering = common.load_numbering()
    model = common.load_complex(path)
    target_chain, binder_chains, _evidence = common.identify_chains(
        model, numbering, human, emit=None)

    blank = dict(design=name, path=str(path), verdict=VERDICT_REFUSED,
                 reason="", matched_atoms=0, rmsd=None, hard=0, contacts=0,
                 closest=None, closest_pos=None, domain3_contacts=0)

    if target_chain is None:
        blank["reason"] = ("no chain in this file reads as human EGFR, so which "
                           "molecule is the target cannot be established")
        return blank
    if not binder_chains:
        blank["reason"] = "no binder chain present alongside the target"
        return blank

    target_residues = common.protein_residues(model[target_chain])
    sup, matched, rmsd = superpose_onto_receptor(target_residues, numbering,
                                                 context)
    blank["matched_atoms"] = matched
    if sup is None:
        blank["reason"] = (f"only {matched} residues could be matched against "
                           f"6ARU, fewer than the {MIN_SUPERPOSITION_ATOMS} "
                           "needed to place the binder")
        return blank
    blank["rmsd"] = rmsd
    if rmsd > MAX_SUPERPOSITION_RMSD and rule_live("superposition"):
        blank["reason"] = (f"the target does not lie on 6ARU: RMSD {rmsd:.2f} A "
                           f"over {matched} CA atoms, above the "
                           f"{MAX_SUPERPOSITION_RMSD} A limit. The binder cannot "
                           "be placed, so no clash figure would mean anything")
        return blank

    binder_atoms = [a for cid in binder_chains
                    for r in common.protein_residues(model[cid])
                    for a in r if a.element != "H"]
    if not binder_atoms:
        blank["reason"] = "the binder chain holds no atoms"
        return blank

    # Move the binder by the same rotation and translation that laid the target
    # on the receptor. Working on a copy of the coordinates leaves the parsed
    # structure untouched, so nothing later in the run sees moved atoms.
    rotation, translation = sup.rotran
    coords = np.array([a.coord for a in binder_atoms]) @ rotation + translation

    hard, contacts, closest, closest_residue = count_overlaps(
        coords, context["outside"])
    _h3, domain3_contacts, _c3, _r3 = count_overlaps(
        coords, list(context["inside"].values()))

    closest_pos = (context["numbering"].uniprot_of(closest_residue.id[1])
                   if closest_residue is not None else None)
    return dict(design=name, path=str(path), verdict=verdict_for(hard, contacts),
                reason="", matched_atoms=matched, rmsd=rmsd, hard=hard,
                contacts=contacts,
                closest=None if closest == float("inf") else closest,
                closest_pos=closest_pos, domain3_contacts=domain3_contacts,
                binder_atoms=len(binder_atoms))


# ---------------------------------------------------------------------------
# Constructed test cases, because no design run has returned anything yet
# ---------------------------------------------------------------------------

BINDER_ATOMS = (("N", "N", (-1.2, 0.0, 0.0)),
                ("CA", "C", (0.0, 0.0, 0.0)),
                ("C", "C", (1.2, 0.0, 0.0)),
                ("O", "O", (1.8, 1.0, 0.0)),
                ("CB", "C", (0.0, 1.5, 0.0)))


def write_test_complex(path, target_residues, binder_centre, name,
                       distort=0.0, include_target=True):
    """Write a target plus a blob of binder atoms, in the format the run writes.

    The target is the real trimmed fragment, so the superposition this exercises is
    the real one. The binder is a compact cluster of alanine residues rather than a
    designed sequence, because this step measures where a binder's volume is and
    nothing else about it.

    `distort` moves every target atom by that many angstroms in a random direction,
    which is how a returned target that is not the structure we supplied would look.
    Moving the target as a whole would not do: a superposition exists precisely to
    undo a rigid movement, so a translated target lands back on the receptor with an
    RMSD of zero, correctly. Only a change of shape can make the superposition fail,
    and only a failed superposition should stop a candidate being scored.

    `include_target=False` writes the binder alone, which is what the binder-only
    files that sit among the run's real output look like.

    Written as mmCIF through Biopython's own writer, the same way step 10 builds
    its test cases, so the tests go through the parser real output will go through
    rather than through a hand-rolled format that only these tests ever exercise.
    """
    structure = Structure("test")
    model = Model(0)
    structure.add(model)
    target_chain, binder_chain = Chain("A"), Chain("B")
    if include_target:
        model.add(target_chain)
    model.add(binder_chain)

    # A fixed seed, so a failing test fails the same way twice and can be looked at.
    rng = np.random.default_rng(0)
    if include_target:
        for residue in target_residues:
            copied = residue.copy()
            for atom in copied:
                coord = np.array(atom.coord, dtype=float)
                if distort:
                    step = rng.normal(size=3)
                    coord = coord + distort * step / np.linalg.norm(step)
                atom.set_coord(coord)
            target_chain.add(copied)

    # A compact globule: alanine residues on a small lattice 3.8 A apart, which is
    # roughly how far apart neighbouring residues sit in a real protein. 27 of them
    # is a body about 11 A across, small for a binder but large enough that an
    # overlap is unambiguous.
    centre = np.array(binder_centre, dtype=float)
    spacing, per_side = 3.8, 3
    resnum, serial = 1, 1
    for i in range(per_side):
        for j in range(per_side):
            for k in range(per_side):
                base = centre + spacing * np.array([i - 1, j - 1, k - 1],
                                                   dtype=float)
                residue = Residue((" ", resnum, " "), "ALA", "")
                for atom_name, element, step in BINDER_ATOMS:
                    residue.add(Atom(atom_name, base + np.array(step),
                                     20.0, 1.0, " ", atom_name, serial,
                                     element=element))
                    serial += 1
                binder_chain.add(residue)
                resnum += 1

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    return path


def test_geometry(context, anchors):
    """Where to put a binder so that it must be clear, and so that it must clash.

    Both are read out of the real structure rather than written in by hand, so the
    test does not quietly stop exercising the thing it was built for if the anchor
    set changes.
    """
    anchor_coords = np.array([context["inside"][p]["CA"].coord
                              for p in anchors if p in context["inside"]])
    face = anchor_coords.mean(axis=0)
    all_coords = np.array([a.coord for r in context["inside"].values()
                           for a in r if a.element != "H"])
    centre = all_coords.mean(axis=0)
    outward = face - centre
    outward /= np.linalg.norm(outward)
    # 14 A out from the anchor face: the near side of a blob this size then sits a
    # few angstroms off the face, which is where a binder goes, and well away from
    # anything outside domain III.
    clear_centre = face + 14.0 * outward
    # Inside the volume the rest of the receptor occupies.
    outside_coords = np.array([a.coord for r in context["outside"]
                               for a in r if a.element != "H"])
    clash_centre = outside_coords.mean(axis=0)
    return clear_centre, clash_centre


def run_tests(context, human, emit):
    """Build the constructed candidates and check each comes back as it must."""
    import tempfile

    anchors = common.load_primary_cluster()
    if not TRIMMED_TARGET.exists():
        emit(f"   {TRIMMED_TARGET} is missing. Run analysis/09_trim_target.py "
             "first.")
        return ["trimmed target missing, so the tests could not be built"]

    fragment = common.protein_residues(
        common.load_complex(TRIMMED_TARGET)[
            next(iter(common.load_complex(TRIMMED_TARGET))).id])
    clear_centre, clash_centre = test_geometry(context, anchors)

    cases = [
        dict(name="binder-on-open-face", centre=clear_centre,
             expect=VERDICT_CLEAR,
             why="a binder out from the anchor face, where a real one would sit"),
        dict(name="binder-inside-domain-iv", centre=clash_centre,
             expect=VERDICT_CLASHING,
             why="a binder in the volume the rest of the receptor occupies"),
        dict(name="target-wrong-shape", centre=clear_centre, distort=6.0,
             expect=VERDICT_REFUSED,
             why="a target whose shape is not 6ARU's, so the binder cannot be "
                 "placed"),
        dict(name="binder-only-file", centre=clear_centre, include_target=False,
             expect=VERDICT_REFUSED,
             why="a binder-only file, which the run's output really does contain"),
    ]

    failures = []
    emit("   | constructed candidate | what it is | expected | got | "
         "overlaps | RMSD |")
    emit("   |---|---|---|---|---|---|")
    with tempfile.TemporaryDirectory() as tmp:
        for case in cases:
            name = case["name"]
            path = write_test_complex(
                Path(tmp) / f"{name}.cif", fragment, case["centre"], name,
                distort=case.get("distort", 0.0),
                include_target=case.get("include_target", True))
            result = analyse_candidate(name, path, context, human, emit)
            got, expected = result["verdict"], case["expect"]
            rmsd = ("—" if result["rmsd"] is None else f"{result['rmsd']:.2f} A")
            emit(f"   | {name} | {case['why']} | {expected} | {got} | "
                 f"{result['hard']} | {rmsd} |")
            if result["reason"]:
                emit(f"   >  {name}: {result['reason']}")
            if got != expected:
                failures.append(f"{name}: expected {expected}, got {got}")
    return failures


# ---------------------------------------------------------------------------
# Real candidates, when there are any
# ---------------------------------------------------------------------------

def find_candidate_structures(folder):
    """Every candidate complex under a design-run output folder.

    Binder-only files sit among the complexes in the run's output, and they have no
    target chain, so they are refused by the chain check rather than filtered by
    filename here.
    """
    if folder is None or not folder.exists():
        return []
    return sorted(p for p in folder.rglob("*")
                  if p.suffix.lower() in (".cif", ".pdb"))


def main(argv=None):
    global BROKEN_RULE

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--candidates", type=Path, default=None, metavar="DIR",
                        help="a design-run output folder to screen. Without it, "
                             "only the constructed tests run.")
    parser.add_argument("--break-rule", choices=list(BREAKABLE_RULES),
                        default=None,
                        help="switch one rule off, to confirm the matching test "
                             "then fails. Expected to fail.")
    args = parser.parse_args(argv)
    BROKEN_RULE = args.break_rule

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    human = common.human_sequence()

    emit("=" * 72)
    emit("FULL-RECEPTOR CLASH CHECK — would each binder fit in the intact "
         "receptor?")
    emit("=" * 72)
    emit()
    emit("The design run sees only the 171-residue fragment, so it cannot tell a "
         "binder")
    emit("that lies on an open face from one that occupies space another part of "
         "the")
    emit("receptor is already in. Both score the same. This step puts each "
         "candidate")
    emit("back into the intact receptor and measures the overlap.")
    emit()
    emit("   Step 12 established that the anchor face itself is reachable. That "
         "is a")
    emit("   property of the target. This is a property of each design: a "
         "reachable")
    emit("   face does not stop one particular binder approaching at an angle "
         "that")
    emit("   puts part of it inside domain II or domain IV.")
    emit()
    if BROKEN_RULE:
        emit(f"   NOTE: --break-rule {BROKEN_RULE} is in effect — "
             f"{BREAKABLE_RULES[BROKEN_RULE]}.")
        emit("   This run is expected to fail.")
        emit()

    context = receptor_context()
    emit(f"1. The intact receptor: {FULL_RECEPTOR.name}")
    emit()
    emit(f"   receptor chain {context['chain']}, "
         f"{context['identity']:.1f}% identity to human EGFR")
    emit(f"   residues the design run saw (domain III, "
         f"{DOMAIN_III[0]}-{DOMAIN_III[1]}): {len(context['inside'])}")
    emit(f"   residues it did not see, which this step measures against: "
         f"{len(context['outside'])}")
    emit(f"   other chains in the file, left out of the measurement: "
         f"{', '.join(context['partners']) or 'none'}")
    emit()
    emit(f"   An overlap below {CLASH_HARD} A is two sets of atoms in the same "
         "place.")
    emit(f"   Between {CLASH_HARD} and {CLASH_CONTACT} A is touching, which the "
         "receptor")
    emit("   may be able to move to accommodate, so it is reported as marginal "
         "rather")
    emit("   than counted as a failure.")
    emit()

    emit("2. Constructed tests, run every time")
    emit()
    emit("   No design run has returned anything yet, so these are the only "
         "candidates")
    emit("   this script has ever scored. Each one exercises one branch.")
    emit()
    failures = run_tests(context, human, emit)
    emit()

    results = []
    folder = args.candidates or (CANDIDATES if CANDIDATES.exists() else None)
    structures = find_candidate_structures(folder)
    emit("3. Real candidates")
    emit()
    if not structures:
        emit("   None found"
             + (f" under {folder}" if folder else "")
             + ". Nothing from the design run exists yet, which is why section 2")
        emit("   is the whole of this script's evidence. Rerun this step with "
             "--candidates")
        emit("   pointing at the run's output folder once there is one.")
    else:
        emit(f"   {len(structures)} file(s) under {folder}")
        emit()
        emit("   | design | verdict | overlaps <2.5 A | contacts | closest | "
             "nearest residue | interface to domain III |")
        emit("   |---|---|---|---|---|---|---|")
        for path in structures:
            result = analyse_candidate(path.stem, path, context, human, emit)
            results.append(result)
            closest = ("—" if result["closest"] is None
                       else f"{result['closest']:.2f} A")
            emit(f"   | {result['design']} | {result['verdict']} | "
                 f"{result['hard']} | {result['contacts']} | {closest} | "
                 f"{result['closest_pos'] or '—'} | "
                 f"{result['domain3_contacts']} |")
            if result["reason"]:
                emit(f"   >  {result['design']}: {result['reason']}")
        emit()
        counts = {}
        for result in results:
            counts[result["verdict"]] = counts.get(result["verdict"], 0) + 1
        for verdict in (VERDICT_CLEAR, VERDICT_MARGINAL, VERDICT_CLASHING,
                        VERDICT_REFUSED):
            emit(f"   {verdict}: {counts.get(verdict, 0)}")
    emit()

    emit("=" * 72)
    emit("4. What this changes about the design")
    emit("=" * 72)
    emit()
    emit("   A clashing candidate is not dropped here. It is flagged with the "
         "overlap")
    emit("   measured, and the decision to discard belongs in step 10's ranking, "
         "so")
    emit("   that every decision about what to submit is made in one place.")
    emit()
    emit("   A clear verdict does not mean the design works. It means the design "
         "is")
    emit("   not ruled out by where its atoms are. Whether it binds, and whether "
         "it")
    emit("   switches with pH, are what step 10 and the experiment decide.")
    emit()
    emit("   This is a rigid measurement: the receptor is treated as unable to "
         "move.")
    emit("   A real receptor can shift to accommodate a small overlap, which is "
         "why a")
    emit("   marginal verdict is reported as marginal rather than as a failure.")
    emit()

    DERIVED.mkdir(parents=True, exist_ok=True)
    csv_path = DERIVED / "13-candidate-clashes.csv"
    with csv_path.open("w") as fh:
        fh.write("design,verdict,hard_overlaps,contacts,closest_a,"
                 "closest_uniprot_pos,domain3_contacts,superposition_rmsd_a,"
                 "matched_ca_atoms,reason\n")
        for result in results:
            rmsd = "" if result["rmsd"] is None else f"{result['rmsd']:.3f}"
            closest = "" if result["closest"] is None else f"{result['closest']:.3f}"
            reason = result["reason"].replace(",", ";")
            fh.write(f"{result['design']},{result['verdict']},{result['hard']},"
                     f"{result['contacts']},{closest},"
                     f"{result['closest_pos'] or ''},"
                     f"{result['domain3_contacts']},{rmsd},"
                     f"{result['matched_atoms']},{reason}\n")
    emit(f"Wrote data/derived/{csv_path.name}"
         + ("" if results else " (header only; no real candidates yet)"))
    emit()

    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED —")
        for failure in failures:
            emit(f"  - {failure}")
    else:
        emit("RESULT: PASSED. Every constructed case came back as it had to: a "
             "binder")
        emit("on the open face reads clear, one inside domain IV reads clashing, "
             "and a")
        emit("target that is not 6ARU is refused rather than scored.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "13-full-receptor-clash.md").write_text(
        "# Would each designed binder still fit in the intact receptor?\n\n"
        "Computed output of `analysis/13_full_receptor_clash.py`. "
        "Do not hand-edit.\n\n"
        "The design run is given the 171-residue fragment and never sees the "
        "other\n450 residues, so it scores a binder that lies on an open face and "
        "one that\noccupies space domain II or domain IV is already in exactly the "
        "same. This\nstep superposes each candidate back into the full 6ARU "
        "structure and measures\nthe overlap.\n\n"
        "Step 12 settled that the anchor face is reachable, which is a property of "
        "the\ntarget. This is a property of each design.\n\n"
        "No design run has returned anything yet, so the evidence below is from "
        "the\nconstructed cases in section 2.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
