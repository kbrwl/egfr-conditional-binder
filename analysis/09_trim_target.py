#!/usr/bin/env python3
"""
09_trim_target.py — cut the receptor down to domain III and write the files the
design run receives.

WHAT THIS IS FOR
----------------
A binder-design run has to fold the whole complex at every step: the target plus
the binder being designed. The cost of folding grows faster than the number of
residues, so the size of the target sets how many attempts we can afford. The
EGFR extracellular region -- the part of the receptor outside the cell, residues
25-645 in our numbering -- is 621 residues. Domain III, the part our epitope sits
on, is 171. Handing the design run the smaller piece is the difference between a
handful of attempts and a few hundred on the same compute budget.

EGFR is the epidermal growth factor receptor, the protein we are designing
against. An epitope is the patch of the target surface a binder is aimed at. An
anchor is a residue inside that patch which can form one half of a pH-sensitive
charge pair.

WHAT THIS COSTS, STATED UP FRONT
--------------------------------
Trimming is a judgement, not a measurement, and it is recorded as unverified in
`docs/decisions-log.md` for that reason. Cutting a domain out of a protein can
leave it unable to hold its shape, and a fragment that does not hold its shape in
the structure predictor presents a surface that does not exist on the intact
protein. Every design made against it would then be aimed at nothing. Nothing in
this script tests that. What it does instead is write down exactly what the cut
breaks, so the risk is visible rather than hidden:

  - every disulfide bond the cut severs. A disulfide bond is a covalent link
    between the sulfur atoms of two cysteine residues, which is the strongest
    thing holding a protein domain's shape together after the backbone itself.
    Cutting a domain out can leave a cysteine whose partner is now outside the
    fragment, and an unpaired cysteine can make the fragment behave differently
    from the intact protein.
  - which anchors sit near a cut edge. The two new chain ends created by the cut
    are the least reliable part of the model, because in the intact protein they
    are held by residues the fragment no longer contains. An anchor next to one
    of them is the least trustworthy part of the result.

This script reports both. It does not try to repair either. A reported risk is
not a repaired one, and the first smoke run against the trimmed target is what
would move the unverified entry, not anything here.

THE NUMBERING, WHICH IS THE THING MOST LIKELY TO GO WRONG
---------------------------------------------------------
Three numbering systems are in play at once, and the design run has to be given
the right one.

  UniProt         every residue number in this project is a position in the
                  human record P00533 in UniProt, the public sequence archive.
                  H370 means UniProt position 370.
  the 6ARU file   the structure file numbers by the mature protein, so its
                  numbers are ours minus 24. H370 is residue 346 in that file.
                  Measured, not assumed, by `analysis/02_structure_prep.py`.
  the trimmed file  this script preserves 6ARU's numbers rather than renumbering
                  from 1, for the reason below.

The design pipeline reads hotspot residues -- the positions on the target it is
told to aim the binder at -- in the numbering of the structure file it is handed.
Give it UniProt numbers and it will aim 24 residues away, at a different surface,
and raise no error. That is the same class of defect that produced the step 04
against step 06 disagreement, which went unnoticed because a wrong lookup
returned a plausible answer rather than failing.

So this script:
  - passes numbers around as `egfr_common.Numbering`, never as bare dictionaries
  - emits the hotspot list in both numberings, side by side
  - checks the mapping by residue identity rather than by arithmetic: position
    344 has to come back as glutamic acid, 358 as histidine, 368 as aspartic
    acid, 370 as histidine, and so on for all eight, read independently from the
    UniProt sequence and from the trimmed structure file
  - can be told to break its own numbering on purpose, so that the check is known
    to work rather than merely present. Run with `--self-test`.

WHY THE TRIMMED FILE KEEPS 6ARU'S NUMBERS RATHER THAN COUNTING FROM 1
---------------------------------------------------------------------
Renumbering the fragment 1 to 171 was considered and rejected. It would create a
fourth numbering system, with nothing in the file recording that it exists, in a
project whose most expensive defect so far was a numbering mix-up. Keeping 6ARU's
numbers means every residue in the trimmed file can be found in the parent file
under the same number, and the one conversion anybody has to do -- add 24 to get
our numbering -- is the conversion `analysis/02` already measured and committed.

Outputs:
  data/structures/6aru_domain3.pdb        the trimmed target, gitignored as a
                                          *.pdb because this script regenerates it
  design/configs/egfr-domain3-h370.json   the hotspot specification, committed
  data/derived/09-hotspot-numbering.csv   the same hotspots in both numberings
  data/derived/09-severed-disulfides.csv  what the cut breaks
  results/findings/09-trimmed-target.md

Reads data/derived/02-numbering-offset.csv and data/derived/08-h370-clusters.csv.
Run analysis/02 and analysis/08 first.

Run standalone:  python analysis/09_trim_target.py
                 python analysis/09_trim_target.py --self-test
Exit code 0 = the fragment was written and every check passed.
"""

import argparse
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch

import egfr_common as common

ROOT = Path(__file__).resolve().parents[1]
RECEPTOR_PDB = ROOT / "data" / "structures" / "6aru_receptor_only.pdb"
PARENT_PDB = ROOT / "data" / "structures" / "6aru.pdb"
TRIMMED_PDB = ROOT / "data" / "structures" / "6aru_domain3.pdb"
CONFIG_DIR = ROOT / "design" / "configs"
CAMPAIGN_NAME = "egfr-domain3-h370"
HOTSPOT_CONFIG = CONFIG_DIR / f"{CAMPAIGN_NAME}.json"
HOTSPOT_NOTE = CONFIG_DIR / f"{CAMPAIGN_NAME}.md"

# Binder length range handed to the design run, in amino acids. The molecule
# category is still open in docs/decisions-log.md; 40-100 is the minibinder band
# and the anchor cluster spans 24.3 angstroms, which is what points there. This
# is a default rather than a decision and the companion note says so.
BINDER_LENGTHS = (60, 100)

# How many accepted designs to stop at, and the cap on attempts. Neither is
# derived from a measured runtime, because no per-design runtime for this
# pipeline has been measured by us or published by it.
FINAL_DESIGNS = 200
MAX_TRAJECTORIES = 2000
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS

STRUCTURE_LABEL = "6ARU"

# The boundary of the cut, in our numbering. Taken from egfr_common rather than
# written again here, because it is also what steps 01 and 06 use and a second
# copy could drift away from them. 310-480 is this project's working definition
# of domain III and is recorded as unverified against the official annotation in
# docs/decisions-log.md; step 06 was written so that it does not depend on these
# numbers being exactly right, and this script does depend on them.
TRIM_START, TRIM_END = common.D3_START, common.D3_END

# The eight anchors of the primary epitope, with the amino acid each one has to
# turn out to be. The positions themselves are not decided here: they are read
# back from step 08's committed cluster file, and this dictionary is the written
# expectation that the read is checked against, so that if step 08's answer ever
# changes this script stops rather than silently trimming for a different set.
#
# Six are acidic -- aspartic acid (D) or glutamic acid (E), which carry a
# negative charge at both pH values -- so the binder gets a histidine opposite
# them. Two are histidines, which are neutral at pH 7.4 and positive at pH 6.5,
# so the binder gets an acidic residue opposite them. Both arrangements switch on
# in the same direction as the surroundings turn acidic.
EXPECTED_ANCHORS = {344: "E", 358: "H", 368: "D", 370: "H",
                    391: "E", 400: "E", 421: "E", 424: "E"}

# An anchor closer than this to one of the two ends the cut creates is reported
# as least trustworthy. 8 angstroms is roughly two residues' reach, and is a
# working convention chosen here rather than a measured property of anything.
CUT_EDGE_WARN = 8.0

# Two cysteine sulfur atoms this close are bonded to each other. A real disulfide
# bond measures about 2.05 angstroms; 2.5 is generous enough to catch a bond in a
# structure solved at modest resolution without catching two sulfurs that merely
# sit near one another, the next-nearest pairs being several angstroms further.
DISULFIDE_MAX = 2.5


class RangeSelect(Select):
    """Keep exactly the residues whose structure-file numbers are listed."""

    def __init__(self, chain_id, keep_resnums):
        self.chain_id = chain_id
        self.keep = set(keep_resnums)

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[1] in self.keep


def load_primary_cluster():
    """The anchor positions, read back from step 08's committed output.

    Step 08 decided which residues cluster on one face around H370 and wrote the
    answer to `data/derived/08-h370-clusters.csv`. That answer is not recomputed
    here and not retyped here: if two scripts need the same quantity, one of them
    computes it and the other reads it, which is the rule in `CLAUDE.md` that
    came out of steps 04 and 06 disagreeing.

    Returns the positions of the largest cluster that contains H370.
    """
    path = DERIVED / "08-h370-clusters.csv"
    if not path.exists():
        raise SystemExit(
            f"{path} is missing. Run analysis/08_h370_epitope.py first; this "
            "step takes the anchor set from its output rather than keeping its "
            "own copy of it.")
    best = None
    for line in path.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        # cluster_size,max_internal_span_a,"344 358 ...",contains_h370
        _before, _, rest = line.partition('"')
        anchors, _, after = rest.partition('"')
        members = [int(tok) for tok in anchors.split()]
        contains_centre = after.strip(", ").strip() == "True"
        if not contains_centre:
            continue
        if best is None or len(members) > len(best):
            best = members
    if best is None:
        raise SystemExit(
            f"{path} contains no cluster that includes H370. Rerun "
            "analysis/08_h370_epitope.py.")
    return sorted(best)


def load_numbering_shifted(shift=0):
    """The numbering analysis/02 measured, optionally broken on purpose.

    With `shift` at zero this is exactly `egfr_common.load_numbering()`. With any
    other value every position is moved by that many places, which is what the
    self-test uses to confirm the identity checks notice. The alternative, a
    test that reasons about what would happen, was rejected: a check that has
    never been seen to fail is not known to work.
    """
    numbering = common.load_numbering()
    if shift == 0:
        return numbering
    path = DERIVED / "02-numbering-offset.csv"
    shifted = {}
    for line in path.read_text().splitlines()[1:]:
        uni, pdb, _offset = line.split(",")
        shifted[int(pdb)] = int(uni) + shift
    return common.Numbering(shifted)


def read_ssbond_records(path):
    """The disulfide bonds the depositors recorded in the parent structure file.

    A PDB file can carry SSBOND records, which are the depositors' own statement
    of which cysteines are bonded to which. Reading them is independent evidence
    against the distances this script measures, so the two can be compared
    instead of trusting either alone.

    Returns {(lower structure number, higher structure number)} for the receptor
    chain, or None if the parent file is not on disk -- it is gitignored and
    re-downloaded by analysis/02, so its absence is normal rather than an error.
    """
    if not path.exists():
        return None
    pairs = set()
    for line in path.read_text().splitlines():
        if not line.startswith("SSBOND"):
            continue
        # Columns are fixed-width: chain and sequence number for each partner.
        try:
            chain_1, num_1 = line[15], int(line[17:21])
            chain_2, num_2 = line[29], int(line[31:35])
        except ValueError:
            continue
        pairs.add((chain_1, chain_2, min(num_1, num_2), max(num_1, num_2)))
    return pairs


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Trim the receptor to domain III for the design run.")
    parser.add_argument(
        "--break-numbering", type=int, default=0, metavar="N",
        help="Shift the measured numbering by N positions on purpose, to "
             "confirm the identity checks catch it. Writes no output files.")
    parser.add_argument(
        "--self-test", action="store_true",
        help="Run this script again with a deliberately broken numbering and "
             "confirm it exits non-zero. A check that has never been seen to "
             "fail is not known to work.")
    args = parser.parse_args(argv)

    if args.self_test:
        return run_self_test()

    out, failures = [], []

    def emit(text=""):
        print(text)
        out.append(text)

    def check(label, ok, detail):
        emit(f"   [{'PASS' if ok else 'FAIL'}] {label}: {detail}")
        if not ok:
            failures.append(f"{label}: {detail}")

    dry_run = args.break_numbering != 0

    emit("=" * 72)
    emit(f"TRIMMING THE TARGET TO DOMAIN III ({TRIM_START}-{TRIM_END})")
    emit("=" * 72)
    emit()
    emit("Residue numbers below are positions in the human record P00533 in")
    emit("UniProt, the public protein sequence archive, unless the line says")
    emit(f"otherwise. The structure is {STRUCTURE_LABEL} from the Protein Data")
    emit("Bank, the public archive of measured three-dimensional structures,")
    emit("with the antibody already removed by analysis/02.")
    emit()
    if dry_run:
        emit("!! RUNNING WITH THE NUMBERING DELIBERATELY SHIFTED BY "
             f"{args.break_numbering:+d}.")
        emit("!! This is the self-test. The identity checks below are supposed")
        emit("!! to fail, and no output file will be written.")
        emit()

    if not RECEPTOR_PDB.exists():
        raise SystemExit(f"{RECEPTOR_PDB} is missing. Run "
                         "analysis/02_structure_prep.py first.")

    numbering = load_numbering_shifted(args.break_numbering)

    structure = PDBParser(QUIET=True).get_structure("receptor",
                                                    str(RECEPTOR_PDB))
    model = structure[0]
    chain = next(iter(model))
    residues = common.protein_residues(chain)
    by_uniprot, uniprot_of_pdb = {}, {}
    for residue in residues:
        pos = numbering.uniprot_of(residue.id[1])
        if pos is not None:
            by_uniprot[pos] = residue
            uniprot_of_pdb[residue.id[1]] = pos

    offset, share, _counts = numbering.dominant_offset()
    emit("1. The numbering this script is working in")
    emit()
    emit("   Measured by analysis/02 and read back from its committed table,")
    emit("   rather than assumed here:")
    emit(f"     our number = {STRUCTURE_LABEL} residue number {offset:+d}, "
         f"holding for {100 * share:.1f}% of the chain")
    emit(f"     chain {chain.id}, {len(residues)} residues with coordinates, "
         f"file numbers {min(r.id[1] for r in residues)}-"
         f"{max(r.id[1] for r in residues)}")
    emit()

    # ---- the anchor set, taken from step 08 rather than decided here ----
    emit("2. The anchor set, read back from step 08")
    emit()
    emit("   Step 08 decided which residues around H370 cluster on one face. It")
    emit("   is not recomputed here and not retyped here. What this script adds")
    emit("   is a written expectation of what that set should be, so that if")
    emit("   step 08's answer changes, this run stops.")
    emit()
    cluster = load_primary_cluster()
    expected_positions = sorted(EXPECTED_ANCHORS)
    check("anchor set matches step 08's largest H370 cluster",
          cluster == expected_positions,
          f"step 08 gives {cluster}"
          + ("" if cluster == expected_positions
             else f", this script expected {expected_positions}"))
    anchors = {pos: EXPECTED_ANCHORS[pos] for pos in expected_positions}
    emit(f"   {len(anchors)} anchors: "
         + ", ".join(f"{aa}{pos}" for pos, aa in sorted(anchors.items())))
    emit()

    # ---- the boundary assertion ----
    emit(f"3. Do all eight anchors fall inside {TRIM_START}-{TRIM_END}?")
    emit()
    emit("   If any does not, the cut would remove a residue the design run has")
    emit("   to aim at. The script stops rather than widening the boundary,")
    emit("   because a boundary that moves to fit the answer is not a boundary.")
    emit()
    outside = [pos for pos in anchors if not TRIM_START <= pos <= TRIM_END]
    for pos in sorted(anchors):
        inside = TRIM_START <= pos <= TRIM_END
        emit(f"     {anchors[pos]}{pos}: "
             f"{'inside' if inside else 'OUTSIDE THE BOUNDARY'}")
    check(f"all anchors inside {TRIM_START}-{TRIM_END}", not outside,
          "all eight within the boundary" if not outside
          else f"{len(outside)} outside: {sorted(outside)} — STOP, do not trim")
    emit()
    if outside:
        emit("   Stopping. Widening the boundary to include these would make the")
        emit("   fragment larger and would also mean the domain III definition")
        emit("   in egfr_common.py is wrong, which is a separate thing to fix.")
        return finish(out, failures, emit, wrote_files=False,
                      write_findings=not dry_run)

    # ---- which residues survive the cut ----
    keep_uniprot = sorted(p for p in by_uniprot
                          if TRIM_START <= p <= TRIM_END)
    keep_pdb = {by_uniprot[p].id[1] for p in keep_uniprot}
    emit("4. What the fragment contains")
    emit()
    segments = contiguous_segments(keep_uniprot)
    emit(f"   {len(keep_uniprot)} residues with coordinates, in "
         + plural(len(segments), "continuous piece", "continuous pieces") + ":")
    for first, last in segments:
        emit(f"     our numbering {first}-{last}  "
             f"({STRUCTURE_LABEL} numbering "
             f"{by_uniprot[first].id[1]}-{by_uniprot[last].id[1]}), "
             f"{last - first + 1} residues")
    requested = TRIM_END - TRIM_START + 1
    missing = [p for p in range(TRIM_START, TRIM_END + 1)
               if p not in by_uniprot]
    emit()
    if missing:
        emit(f"   The boundary asks for {requested} residues and the structure")
        emit(f"   resolves {len(keep_uniprot)} of them. "
             + plural(len(missing), "residue has", "residues have")
             + " no coordinates in this file, so they cannot be written into")
        emit("   the fragment and the design run will see a gap where they are.")
        emit(f"     missing: {compact_ranges(missing)}")
    else:
        emit(f"   The boundary asks for {requested} residues and the structure")
        emit("   resolves every one of them, so the fragment has no gap the")
        emit("   trimming introduced and none it inherited.")
    emit()
    emit("   The fragment is bare protein. The sugar chains and ions in the")
    emit("   deposited file were already dropped by analysis/02, including the")
    emit("   one attached at N444, so the design run sees a surface with no")
    emit("   sugars on it. That is a known simplification rather than a")
    emit("   measurement, and it applies to the untrimmed target equally.")
    emit()
    emit(f"   For scale: the extracellular region is "
         f"{common.ECD_END - common.ECD_START + 1} residues and this fragment "
         f"is {len(keep_uniprot)}, which is "
         f"{100.0 * len(keep_uniprot) / (common.ECD_END - common.ECD_START + 1):.0f}% "
         "of it.")
    emit()

    # ---- the identity check, three independent routes ----
    emit("5. Is the mapping right? Checked by residue identity, not arithmetic")
    emit()
    emit("   An offset error does not raise an error; it returns a plausible")
    emit("   residue 24 positions away. The only way to catch it is to ask what")
    emit("   amino acid each position actually turns out to be, and compare")
    emit("   that against what it has to be. Three independent routes:")
    emit()
    emit("     written    the expectation recorded at the top of this script")
    emit("     sequence   read out of the human UniProt sequence")
    emit(f"     structure  read out of the {STRUCTURE_LABEL} coordinates, via the")
    emit("                offset analysis/02 measured")
    emit()
    human = common.human_sequence()
    emit("   | our pos | written | sequence | structure | file pos | agree |")
    emit("   |---|---|---|---|---|---|")
    hotspot_rows = []
    for pos in sorted(anchors):
        want = anchors[pos]
        from_sequence = human[pos - 1]
        residue = by_uniprot.get(pos)
        if residue is None:
            from_structure, pdb_num = "-", None
        else:
            from_structure = common.THREE_TO_ONE.get(residue.get_resname(), "X")
            pdb_num = residue.id[1]
        agree = (want == from_sequence == from_structure)
        emit(f"   | {pos} | {want} | {from_sequence} | {from_structure} | "
             f"{pdb_num if pdb_num is not None else 'unresolved'} | "
             f"{'yes' if agree else 'NO'} |")
        hotspot_rows.append(dict(uniprot_pos=pos, expected_aa=want,
                                 sequence_aa=from_sequence,
                                 structure_aa=from_structure,
                                 trimmed_resnum=pdb_num, agree=agree))
    emit()
    disagreeing = [r for r in hotspot_rows if not r["agree"]]
    check("all eight anchors read as the expected amino acid by all three routes",
          not disagreeing,
          "every anchor agrees" if not disagreeing
          else f"{len(disagreeing)} disagree: "
               + ", ".join(f"{r['uniprot_pos']} expected {r['expected_aa']}, "
                           f"sequence {r['sequence_aa']}, structure "
                           f"{r['structure_aa']}" for r in disagreeing))
    emit()
    if disagreeing:
        emit("   Stopping before anything is written. A hotspot list emitted on")
        emit("   a broken mapping would aim the design run at the wrong surface")
        emit("   and nothing downstream would notice.")
        return finish(out, failures, emit, wrote_files=False,
                      write_findings=not dry_run)

    if dry_run:
        emit("   The self-test expected these checks to fail and they passed,")
        emit("   which means the checks do not catch a shifted numbering.")
        failures.append("self-test: shifted numbering was not detected")
        return finish(out, failures, emit, wrote_files=False,
                      write_findings=False)

    # ---- what the cut breaks: disulfide bonds ----
    emit("6. What the cut breaks: disulfide bonds")
    emit()
    emit("   A disulfide bond is a covalent link between the sulfur atoms of two")
    emit("   cysteine residues, written here as the SG atom. It is the strongest")
    emit("   thing holding a domain's shape together after the backbone, and")
    emit("   EGFR's extracellular region is held together by a lot of them.")
    emit()
    emit("   A cut that leaves one partner inside the fragment and the other")
    emit("   outside it leaves a cysteine with nothing to bond to. That can make")
    emit("   the fragment fold differently from the same residues in the intact")
    emit("   protein, which is the specific failure the unverified entry in the")
    emit("   decisions log is about. This section lists them. It does not fix")
    emit("   them: a reported risk is not a repaired one.")
    emit()
    bonds, free = find_disulfides(residues, uniprot_of_pdb)
    severed = [b for b in bonds
               if (TRIM_START <= b["a_uniprot"] <= TRIM_END)
               != (TRIM_START <= b["b_uniprot"] <= TRIM_END)]
    internal = [b for b in bonds
                if TRIM_START <= b["a_uniprot"] <= TRIM_END
                and TRIM_START <= b["b_uniprot"] <= TRIM_END]
    emit("   Measured: "
         + plural(len(bonds), "disulfide bond", "disulfide bonds")
         + " in the receptor chain, SG to SG within "
         + f"{DISULFIDE_MAX} angstroms.")
    emit()

    ssbond = read_ssbond_records(PARENT_PDB)
    if ssbond is None:
        emit(f"   The parent file {PARENT_PDB.name} is not on disk, so the")
        emit("   depositors' own SSBOND records could not be read and the")
        emit("   measured bonds stand unchecked. Rerun analysis/02 to download")
        emit("   it and this comparison will appear.")
    else:
        receptor_ssbond = {(lo, hi) for c1, c2, lo, hi in ssbond
                           if c1 == c2 == chain.id}
        measured = {(min(b["a_resnum"], b["b_resnum"]),
                     max(b["a_resnum"], b["b_resnum"])) for b in bonds}
        only_measured = sorted(measured - receptor_ssbond)
        only_recorded = sorted(receptor_ssbond - measured)
        emit("   Checked against the SSBOND records in the deposited file, which")
        emit("   are the depositors' own statement of which cysteines are bonded:")
        check("measured bonds match the file's own SSBOND records",
              not only_measured and not only_recorded,
              f"{len(measured)} bonds, identical to the file's records"
              if not only_measured and not only_recorded
              else f"measured but not recorded {only_measured}; "
                   f"recorded but not measured {only_recorded}")
        if only_recorded:
            emit("   A bond the file records but we did not measure usually means")
            emit("   one of its two cysteines has no coordinates in this file.")
    emit()

    if severed:
        emit("   SEVERED BY THE CUT — one partner inside the fragment, one")
        emit("   outside it:")
        emit()
        emit("   | inside | partner outside | SG-SG distance | partner's domain |")
        emit("   |---|---|---|---|")
        for b in sorted(severed, key=lambda b: b["a_uniprot"]):
            inside_pos, outside_pos = b["a_uniprot"], b["b_uniprot"]
            if not TRIM_START <= inside_pos <= TRIM_END:
                inside_pos, outside_pos = outside_pos, inside_pos
            emit(f"   | C{inside_pos} | C{outside_pos} | "
                 f"{b['distance']:.2f} A | {where(outside_pos)} |")
        emit()
        emit("   " + plural(len(severed), "cysteine in the fragment loses",
                              "cysteines in the fragment lose")
             + " the partner it has in the intact protein.")
        emit("   Each is a place where the fragment could behave differently from")
        emit("   the intact protein. Whether it actually does is not measured")
        emit("   here and would be answered by predicting the fragment alone and")
        emit("   comparing it against the same residues in 6ARU, which is item 11")
        emit("   in the decisions log's next actions.")
    else:
        emit("   None. Every disulfide bond with a partner inside the fragment")
        emit("   has its other partner inside it too, so the cut leaves no")
        emit("   cysteine without the partner it has in the intact protein.")
    emit()
    emit("   For context: "
         + plural(len(internal), "disulfide bond lies", "disulfide bonds lie")
         + " wholly inside the fragment and "
         + ("is" if len(internal) == 1 else "are") + " unaffected, and "
         + f"{len(bonds) - len(internal) - len(severed)} lie wholly outside it")
    emit("   and are removed along with both partners.")
    free_inside = [p for p in free if TRIM_START <= p <= TRIM_END]
    if free_inside:
        emit()
        emit("   Separately, these cysteines inside the fragment are already")
        emit("   unpaired in the intact structure, so the cut is not what makes")
        emit(f"   them so: {sorted(free_inside)}")
    emit()

    # ---- what the cut breaks: anchors near a cut edge ----
    emit("7. What the cut breaks: anchors near a cut edge")
    emit()
    emit("   The cut creates new chain ends. In the intact protein those")
    emit("   residues are held by neighbours the fragment no longer contains, so")
    emit("   they are the part of the model least likely to be right, and an")
    emit("   anchor beside one of them is the least trustworthy part of the")
    emit("   result.")
    emit()
    emit("   Measured as the closest approach between any non-hydrogen atom of")
    emit("   the anchor and any non-hydrogen atom of a cut-edge residue. A")
    emit("   cut-edge residue is one whose sequence neighbour is present in the")
    emit("   parent structure and absent from the fragment. Gaps the deposited")
    emit("   structure already had are not cut edges: the cut did not make them.")
    emit()
    edges = cut_edges(by_uniprot, keep_uniprot)
    emit(f"   Cut edges: {sorted(edges)}")
    emit()
    emit(f"   | anchor | nearest cut edge | distance | within {CUT_EDGE_WARN:.0f} A |")
    emit("   |---|---|---|---|")
    near_edge = []
    edge_rows = {}
    for pos in sorted(anchors):
        nearest, distance = nearest_edge(by_uniprot[pos],
                                        [by_uniprot[e] for e in sorted(edges)],
                                        uniprot_of_pdb)
        flagged = distance <= CUT_EDGE_WARN
        edge_rows[pos] = (nearest, distance, flagged)
        if flagged:
            near_edge.append(pos)
        emit(f"   | {anchors[pos]}{pos} | {nearest} | {distance:.1f} A | "
             f"{'YES' if flagged else 'no'} |")
    emit()
    if near_edge:
        emit("   " + plural(len(near_edge), "anchor sits", "anchors sit")
             + f" within {CUT_EDGE_WARN:.0f} angstroms of a cut edge: "
             + ", ".join(f"{anchors[p]}{p}" for p in near_edge))
        emit("   Treat any design that relies on these as the least supported.")
        emit("   This does not say they are wrong. It says that if the fragment")
        emit("   does not hold its shape, these are where it will show first.")
    else:
        emit(f"   None. Every anchor is further than {CUT_EDGE_WARN:.0f} "
             "angstroms from both ends the cut created, which is the better")
        emit("   outcome: the part of the fragment the design run aims at is the")
        emit("   part furthest from the damage.")
    emit()

    # ---- write the fragment ----
    emit("8. Writing the fragment")
    emit()
    io = PDBIO()
    io.set_structure(structure)
    TRIMMED_PDB.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(TRIMMED_PDB), select=RangeSelect(chain.id, keep_pdb))
    written = PDBParser(QUIET=True).get_structure("trim", str(TRIMMED_PDB))[0]
    written_chain = next(iter(written))
    written_residues = common.protein_residues(written_chain)
    check("fragment written with the expected residue count",
          len(written_residues) == len(keep_uniprot),
          f"{len(written_residues)} residues in "
          f"data/structures/{TRIMMED_PDB.name}, expected {len(keep_uniprot)}")
    file_numbers = [r.id[1] for r in written_residues]
    emit(f"   The fragment's own numbering range: "
         f"{min(file_numbers)}-{max(file_numbers)}, chain {written_chain.id}.")
    emit(f"   That is {STRUCTURE_LABEL}'s numbering, preserved. Add {offset} to "
         "get ours.")
    emit()

    # Re-read the anchors out of the file that was actually written, rather than
    # trusting that the file contains what we selected. This is the check that
    # the hotspot numbers handed to the design run point at real residues in the
    # real file.
    emit("   Anchors read back out of the written file:")
    reread = {r.id[1]: common.THREE_TO_ONE.get(r.get_resname(), "X")
              for r in written_residues}
    mismatched = []
    for pos in sorted(anchors):
        pdb_num = by_uniprot[pos].id[1]
        got = reread.get(pdb_num)
        ok = got == anchors[pos]
        emit(f"     file residue {pdb_num} reads {got or 'ABSENT'}, "
             f"expected {anchors[pos]} (our {pos}) — "
             f"{'ok' if ok else 'MISMATCH'}")
        if not ok:
            mismatched.append(pos)
    check("every anchor is present in the written file as the right residue",
          not mismatched,
          "all eight" if not mismatched else f"wrong or absent: {mismatched}")
    emit()

    # ---- will the design pipeline accept this file at all? ----
    # BindCraft2's structure loader refuses some things outright and silently
    # drops others, both of which are worth knowing before a GPU is paid for.
    # Read from bindcraft/protein.py in PacesaLab/BindCraft2 on 1 October 2026:
    # insertion codes are rejected with an error telling you to renumber; only
    # the twenty standard amino acids are accepted; and residues missing backbone
    # atoms are dropped without stopping the run. A dropped residue would shift
    # nothing, because numbering is preserved, but it would remove a hotspot.
    emit("   Will the design pipeline accept this file?")
    emit()
    emit("   Three things its structure loader does, read out of its own source")
    emit("   rather than discovered on a paid GPU: it refuses a file containing")
    emit("   insertion codes, it refuses residues outside the standard twenty,")
    emit("   and it drops residues missing backbone atoms without stopping.")
    emit()
    insertion_codes = [r.id for r in written_chain if r.id[2] != " "]
    check("no insertion codes", not insertion_codes,
          "none" if not insertion_codes
          else f"{len(insertion_codes)} present: {insertion_codes[:5]} — the "
               "loader will refuse this file")
    nonstandard = [r.get_resname() for r in written_chain
                   if r.id[0] == " " and r.get_resname() not in common.THREE_TO_ONE]
    check("every residue is one of the standard twenty", not nonstandard,
          f"{len(written_residues)} standard residues and nothing else"
          if not nonstandard else f"non-standard present: {sorted(set(nonstandard))}")
    backbone = ("N", "CA", "C")
    incomplete = [r.id[1] for r in written_residues
                  if any(atom not in r for atom in backbone)]
    anchor_incomplete = [pos for pos in anchors
                         if by_uniprot[pos].id[1] in incomplete]
    check("every residue has a complete backbone, so none is dropped",
          not incomplete,
          "all present" if not incomplete
          else f"{len(incomplete)} missing one of N, CA or C: {incomplete[:8]}"
               + (f" — including the anchors {anchor_incomplete}"
                  if anchor_incomplete else ""))
    emit()

    # ---- the hotspot specification ----
    emit("9. The hotspot specification handed to the design run")
    emit()
    emit("   Hotspots are the positions on the target the design run is told to")
    emit("   aim the binder at. They are given in the numbering of the structure")
    emit("   file supplied alongside them, which here is the trimmed file, which")
    emit(f"   carries {STRUCTURE_LABEL}'s numbers. Both numberings are written out")
    emit("   so a reader can check one against the other.")
    emit()
    emit("   | role | our numbering | trimmed-file numbering | binder residue |")
    emit("   |---|---|---|---|")
    for pos in sorted(anchors):
        aa = anchors[pos]
        pdb_num = by_uniprot[pos].id[1]
        binder = ("histidine" if aa in "DE"
                  else "aspartic or glutamic acid")
        role = ("acidic, negative at both pH values" if aa in "DE"
                else "histidine, neutral at pH 7.4 and positive at pH 6.5")
        emit(f"   | {aa}{pos} — {role} | {pos} | "
             f"{written_chain.id}{pdb_num} | {binder} |")
    emit()
    hotspot_string = ",".join(
        f"{written_chain.id}{by_uniprot[pos].id[1]}" for pos in sorted(anchors))
    emit(f"   As one string: {hotspot_string}")
    emit()
    emit("   The same eight in our numbering, which is NOT what goes in the")
    emit("   config file: "
         + ",".join(str(p) for p in sorted(anchors)))
    emit()

    config_text = write_hotspot_config(
        chain_id=written_chain.id, anchors=anchors, by_uniprot=by_uniprot,
        offset=offset, keep_uniprot=keep_uniprot, severed=severed,
        near_edge=near_edge)
    emit(f"   Wrote design/configs/{HOTSPOT_CONFIG.name}")
    emit(f"   Wrote design/configs/{HOTSPOT_NOTE.name}")
    emit()
    emit("   WHAT IS CONFIRMED ABOUT THAT CONFIG AND WHAT IS NOT.")
    emit()
    emit("   Confirmed. The key names and the hotspot syntax were read out of")
    emit("   BindCraft2's own repository on 1 October 2026, from its shipped")
    emit("   example campaign file, its target presets, and its settings")
    emit("   reference: a campaign file is JSON; targets[] takes name,")
    emit("   target_path, chains and hotspots; hotspots is one comma-separated")
    emit("   string in the supplied structure's own numbering, with an optional")
    emit("   chain-letter prefix. The residue numbers and the chain letter in it")
    emit("   are computed and checked above by three independent routes.")
    emit()
    emit("   Not confirmed. The binder length range and the two design counts are")
    emit("   defaults rather than decisions, because the molecule category is")
    emit("   still open and no per-design runtime has been measured. BindCraft2")
    emit("   also registers a filter metric called Target_Crop_Length, which")
    emit("   suggests it may crop the target itself; whether that interacts with")
    emit("   this trim, or renumbers the target in the output, has not been")
    emit("   checked. It matters because step 10 translates the output's target")
    emit("   numbering back into ours, and a renumbering would break that")
    emit("   silently. Both are written down in the companion note.")
    emit()
    emit("   The config file itself carries no commentary. BindCraft2 checks a")
    emit("   campaign file's keys against a fixed catalogue of setting names, so")
    emit("   a key added for a reader's benefit could make the run refuse the")
    emit("   file. Everything explanatory is in the markdown file beside it.")
    emit()
    del config_text

    # ---- CSVs ----
    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "09-hotspot-numbering.csv").open("w") as fh:
        fh.write("uniprot_pos,trimmed_resnum,chain,aa,role,binder_partner,"
                 "nearest_cut_edge_uniprot,distance_to_cut_edge_a,"
                 "within_cut_edge_warning\n")
        for pos in sorted(anchors):
            aa = anchors[pos]
            nearest, distance, flagged = edge_rows[pos]
            role = "acidic" if aa in "DE" else "target histidine"
            partner = "H" if aa in "DE" else "D or E"
            fh.write(f"{pos},{by_uniprot[pos].id[1]},{written_chain.id},{aa},"
                     f"{role},{partner},{nearest},{distance:.3f},{flagged}\n")
    with (DERIVED / "09-severed-disulfides.csv").open("w") as fh:
        fh.write("inside_uniprot,outside_uniprot,inside_resnum,outside_resnum,"
                 "sg_sg_distance_a,outside_region\n")
        for b in sorted(severed, key=lambda b: b["a_uniprot"]):
            inside_pos, outside_pos = b["a_uniprot"], b["b_uniprot"]
            inside_num, outside_num = b["a_resnum"], b["b_resnum"]
            if not TRIM_START <= inside_pos <= TRIM_END:
                inside_pos, outside_pos = outside_pos, inside_pos
                inside_num, outside_num = outside_num, inside_num
            fh.write(f"{inside_pos},{outside_pos},{inside_num},{outside_num},"
                     f"{b['distance']:.3f},{where(outside_pos)}\n")
    emit("   Wrote data/derived/09-hotspot-numbering.csv")
    emit("   Wrote data/derived/09-severed-disulfides.csv")
    emit()

    emit("10. What this changes about the design")
    emit()
    emit(f"   The design run now has a {len(keep_uniprot)}-residue target rather")
    emit(f"   than a {common.ECD_END - common.ECD_START + 1}-residue one, which")
    emit("   is what makes a few hundred attempts affordable instead of a few.")
    emit("   The ratio of run times is not estimated here; it is item 9 in the")
    emit("   decisions log's next actions, to be measured on the first run by")
    emit("   timing one trajectory on each.")
    emit()
    if severed:
        emit("   " + plural(len(severed), "severed disulfide bond is",
                              "severed disulfide bonds are")
             + " the price. That is the reason the trimming")
        emit("   decision stays unverified: nothing here measures whether the")
        emit("   fragment still holds the shape its anchors sit on.")
    if near_edge:
        emit("   " + plural(len(near_edge), "anchor sits", "anchors sit")
             + " near a cut edge, so any design leaning on")
        emit("   that position should be ranked below one that does not.")
    if not severed and not near_edge:
        emit("   No disulfide bond is severed and no anchor sits near a cut edge,")
        emit("   which is the best available outcome from a cut of this size. It")
        emit("   does not make the trimming decision verified: a fragment can")
        emit("   still fail to hold its shape for reasons a bond count and an")
        emit("   edge distance do not capture.")
    emit()

    return finish(out, failures, emit, wrote_files=True,
                  summary=(f"Fragment: {len(keep_uniprot)} residues, "
                           f"our {TRIM_START}-{TRIM_END}, file numbers "
                           f"{min(file_numbers)}-{max(file_numbers)}. "
                           f"{len(severed)} disulfide bond(s) severed, "
                           f"{len(near_edge)} anchor(s) near a cut edge."))


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------

def plural(n, singular, plural_form):
    """'1 disulfide bond' rather than '1 disulfide bond(s)'."""
    return f"{n} {singular if n == 1 else plural_form}"


def contiguous_segments(positions):
    """Turn a sorted list of positions into (first, last) runs."""
    segments, start, prev = [], positions[0], positions[0]
    for pos in positions[1:]:
        if pos != prev + 1:
            segments.append((start, prev))
            start = pos
        prev = pos
    segments.append((start, prev))
    return segments


def compact_ranges(positions):
    """'310-312, 319' rather than a list of every number."""
    return ", ".join(f"{a}-{b}" if a != b else str(a)
                     for a, b in contiguous_segments(sorted(positions)))


def where(pos):
    """Which part of the receptor a position falls in, in words."""
    if pos < common.ECD_START:
        return "signal peptide"
    if pos < TRIM_START:
        return f"domains I-II ({common.ECD_START}-{TRIM_START - 1})"
    if pos <= TRIM_END:
        return f"domain III ({TRIM_START}-{TRIM_END})"
    if pos <= common.ECD_END:
        # Named by position rather than by domain, because this project has not
        # confirmed where domain IV begins and ends and does not need to here.
        return f"downstream of domain III ({TRIM_END + 1}-{common.ECD_END})"
    return "inside the cell or across the membrane"


def find_disulfides(residues, uniprot_of_pdb):
    """Every pair of cysteines bonded to each other, and every unpaired one.

    Measured between SG atoms, the sulfur at the end of a cysteine side chain,
    which is the atom the bond is actually between. Returns (bonds, unpaired),
    where bonds is a list of dicts and unpaired is a list of our residue numbers
    for cysteines with an SG atom and no partner within the cutoff.
    """
    sulfurs = []
    for residue in residues:
        if common.THREE_TO_ONE.get(residue.get_resname()) != "C":
            continue
        if "SG" not in residue:
            continue            # side chain unresolved; nothing to measure
        pos = uniprot_of_pdb.get(residue.id[1])
        if pos is None:
            continue
        sulfurs.append((pos, residue.id[1], residue["SG"]))

    bonds, bonded = [], set()
    for (pos_a, num_a, sg_a), (pos_b, num_b, sg_b) in itertools.combinations(
            sulfurs, 2):
        distance = float(np.linalg.norm(sg_a.coord - sg_b.coord))
        if distance > DISULFIDE_MAX:
            continue
        lo, hi = sorted(((pos_a, num_a), (pos_b, num_b)))
        bonds.append(dict(a_uniprot=lo[0], a_resnum=lo[1],
                          b_uniprot=hi[0], b_resnum=hi[1],
                          distance=distance))
        bonded.update({pos_a, pos_b})
    unpaired = [pos for pos, _num, _sg in sulfurs if pos not in bonded]
    return sorted(bonds, key=lambda b: b["a_uniprot"]), sorted(unpaired)


def cut_edges(by_uniprot, keep_uniprot):
    """Fragment residues the cut created a new chain end at.

    A kept residue is a cut edge when the residue next to it in the sequence has
    coordinates in the parent structure and is not in the fragment. A gap the
    deposited structure already had is not a cut edge, because the cut is not
    what made it.
    """
    kept = set(keep_uniprot)
    edges = set()
    for pos in kept:
        for neighbour in (pos - 1, pos + 1):
            if neighbour in by_uniprot and neighbour not in kept:
                edges.add(pos)
    return edges


def nearest_edge(anchor_residue, edge_residues, uniprot_of_pdb):
    """Closest approach from an anchor to any cut-edge residue.

    Returns (our number of the nearest edge residue, distance in angstroms).
    Measured between non-hydrogen atoms, because hydrogens are absent from these
    files and "heavy atom" therefore means the atoms we have positions for.
    """
    if not edge_residues:
        return None, float("inf")
    edge_atoms = common.heavy_atoms(edge_residues)
    search = NeighborSearch(edge_atoms)
    best = (float("inf"), None)
    for atom in anchor_residue:
        if atom.element == "H":
            continue
        # Widen until something is found, so the answer is a real distance
        # rather than "further than the cutoff".
        radius = CUT_EDGE_WARN
        found = []
        while not found and radius <= 200:
            found = search.search(atom.coord, radius)
            radius *= 2
        for near in found:
            distance = atom - near
            if distance < best[0]:
                best = (distance,
                        uniprot_of_pdb.get(near.get_parent().id[1]))
    return best[1], best[0]


def write_hotspot_config(chain_id, anchors, by_uniprot, offset, keep_uniprot,
                         severed, near_edge):
    """The campaign file BindCraft2 is given, plus a companion note about it.

    Two files rather than one, and the reason matters. BindCraft2 validates the
    keys in a campaign file against a fixed catalogue of setting names. A comment
    field invented for a reader's benefit is not in that catalogue and could make
    the run refuse the file, so the JSON holds nothing but keys that were read
    out of BindCraft2's own shipped example and target-preset files. Everything a
    person needs in order to decide whether to trust it goes in the markdown file
    beside it.

    The alternative, embedding the provenance as underscore-prefixed keys, was
    rejected for that reason: a config file that explains itself and then fails
    to load is worse than one that loads and is explained next door.
    """
    hotspots = ",".join(f"{chain_id}{by_uniprot[pos].id[1]}"
                        for pos in sorted(anchors))
    first = min(by_uniprot[p].id[1] for p in keep_uniprot)
    last = max(by_uniprot[p].id[1] for p in keep_uniprot)

    campaign = {
        "campaign_name": CAMPAIGN_NAME,
        "project_folder": f"results/candidates/{CAMPAIGN_NAME}",
        "modality": "binder",
        "targets": [{
            "name": "EGFR_domain3",
            # Resolved relative to this file's own directory, per BindCraft2's
            # reference documentation. Change it if the campaign is launched
            # from somewhere the repository is laid out differently, which is
            # the case on Modal.
            "target_path": "../../data/structures/6aru_domain3.pdb",
            "chains": chain_id,
            "hotspots": hotspots,
        }],
        "binder_lengths": list(BINDER_LENGTHS),
        "number_of_final_designs": FINAL_DESIGNS,
        "max_trajectories": MAX_TRAJECTORIES,
    }
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    HOTSPOT_CONFIG.write_text(json.dumps(campaign, indent=2) + "\n")

    rows = []
    for pos in sorted(anchors):
        aa = anchors[pos]
        rows.append(
            f"| {aa}{pos} | {pos} | {chain_id}{by_uniprot[pos].id[1]} | "
            f"{'acidic' if aa in 'DE' else 'target histidine'} | "
            f"{'histidine' if aa in 'DE' else 'aspartic or glutamic acid'} | "
            f"{'yes' if pos in near_edge else 'no'} |")

    HOTSPOT_NOTE.write_text(f"""# Notes on `{HOTSPOT_CONFIG.name}`

Generated by `analysis/09_trim_target.py`. Do not hand-edit either file; rerun
that script.

`{HOTSPOT_CONFIG.name}` is the campaign file handed to BindCraft2. It is kept
free of commentary because BindCraft2 checks the keys in a campaign file against
a fixed catalogue of setting names, and a key invented for a reader could make
the run refuse the file. Everything explanatory is here instead.

## The hotspots, in both numberings

Hotspots are the positions on the target BindCraft2 is told to aim the binder at.
They are given in the numbering of the structure file supplied alongside them.
That file is `data/structures/6aru_domain3.pdb`, which preserves {STRUCTURE_LABEL}'s
own numbering, so **our UniProt numbering = the number in the config {offset:+d}**.

Putting UniProt numbers in the config would aim the run {abs(offset)} residues
away, at a different surface, and raise no error. That is the single most likely
way for this file to be wrong.

| anchor | our numbering | in the config | role on the target | what the binder gets | near a cut edge |
|---|---|---|---|---|---|
{chr(10).join(rows)}

Six anchors are acidic -- aspartic acid (D) or glutamic acid (E), negatively
charged at both pH values -- so the binder gets a histidine opposite them. Two are
histidines, neutral at pH 7.4 and positive at pH 6.5, so the binder gets an acidic
residue opposite them. Both arrangements switch on as the surroundings turn acidic.

## The fragment

Our numbering {TRIM_START}-{TRIM_END}, {len(keep_uniprot)} residues with
coordinates, chain {chain_id}, config numbering {first}-{last}.
{plural(len(severed), 'disulfide bond is', 'disulfide bonds are')} severed by the
cut and {plural(len(near_edge), 'anchor sits', 'anchors sit')} within
{CUT_EDGE_WARN:.0f} angstroms of a cut edge. Both are detailed in
`results/findings/09-trimmed-target.md`.

## What is confirmed about this config and what is not

Read out of BindCraft2's own repository (`PacesaLab/BindCraft2`, main branch) on
1 October 2026, from `examples/pdl1_custom_target.json`, `settings/target/*.json`
and `docs/source/reference.md`:

- a campaign file is JSON, passed as `bindcraft design <file>`
- `targets[]` entries take `name`, `target_path`, `chains`, `hotspots`,
  `coldspots`, `weight` and `objective`, and nothing else
- `hotspots` is one comma-separated string, not a list. A chain-letter prefix is
  optional and unprefixed numbers address the first selected chain. The
  documentation says to use the residue numbers from your own structure, which is
  what this file does.
- `binder_lengths` takes `[low, high]` as an inclusive range
- `target_path` is resolved relative to the campaign file's own directory
- `number_of_final_designs` is the count of accepted designs to stop at, and
  `max_trajectories` caps the attempts

Computed and checked by `analysis/09`, by three independent routes: the residue
numbers, the chain letter, and the amino acid each hotspot turns out to be.

**Not confirmed.** `binder_lengths` of {BINDER_LENGTHS[0]}-{BINDER_LENGTHS[1]} is a
default rather than a decision: the molecule category is still open in
`docs/decisions-log.md`, and the minibinder band is 40-100 amino acids. The counts
{FINAL_DESIGNS} and {MAX_TRAJECTORIES} are not derived from a measured runtime,
because no per-design runtime for BindCraft2 has been measured by us or published
by it.

**One thing to check before the first real run.** BindCraft2 registers a filter
metric named `Target_Crop_Length`, which suggests it may crop the target itself.
If it does, cropping could interact with the trim this script performs, and could
renumber the target in the output. Neither has been checked. It matters because
`analysis/10_charge_pair_filter.py` translates the target numbering in
BindCraft2's output back into ours, and a renumbering there would break that
silently.

## What this file deliberately does not do

Nothing here constrains binding strength. The requirement is no *detectable*
binding at pH 7.4, which is a threshold and not a ratio, so a weak binder that
clearly switches beats a strong one. BindCraft2 maximises confidence and interface
quality by default and this file does not stop it. That is handled downstream by
`analysis/10_charge_pair_filter.py`, which ranks by correct charge pairs and
deliberately does not rank by binding strength.
""")
    return campaign


def run_self_test():
    """Break the numbering on purpose and confirm the script fails.

    `CLAUDE.md`: a check that has never been seen to fail is not known to work.
    This runs the script again as a separate process with the numbering shifted
    by one position, so the exit code is genuinely observed rather than reasoned
    about. Every anchor should then read as the wrong amino acid and the run
    should stop before writing anything.
    """
    print("=" * 72)
    print("SELF-TEST: does the identity check catch a shifted numbering?")
    print("=" * 72)
    print()
    print("Running this script again with the measured numbering shifted by one")
    print("position. Every anchor should read as the wrong amino acid, the run")
    print("should stop before writing any file, and the exit code should be")
    print("non-zero. If it is zero, the check is decoration.")
    print()
    before = TRIMMED_PDB.stat().st_mtime if TRIMMED_PDB.exists() else None
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--break-numbering", "1"],
        capture_output=True, text=True, cwd=str(Path(__file__).resolve().parent))
    tail = [ln for ln in proc.stdout.splitlines() if "[FAIL]" in ln
            or "[PASS]" in ln or "RESULT" in ln]
    for line in tail:
        print("  " + line.strip())
    print()
    print(f"Exit code: {proc.returncode}")
    after = TRIMMED_PDB.stat().st_mtime if TRIMMED_PDB.exists() else None
    untouched = before == after
    print(f"Output file left untouched: {'yes' if untouched else 'NO'}")
    print()
    if proc.returncode != 0 and untouched:
        print("SELF-TEST PASSED. The check catches a shifted numbering, stops")
        print("the run, and writes nothing.")
        return 0
    print("SELF-TEST FAILED. A shifted numbering did not stop the run"
          + ("" if untouched else " and output was overwritten") + ".")
    print("Fix the check before trusting any hotspot list this script emits.")
    return 1


def finish(out, failures, emit, wrote_files, summary="", write_findings=True):
    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED — "
             "do not use any output of this run.")
        for failure in failures:
            emit(f"  - {failure}")
    else:
        emit("RESULT: PASSED.")
        if summary:
            emit(summary)
        if wrote_files:
            emit(f"Trimmed target at data/structures/{TRIMMED_PDB.name}; "
                 f"hotspots at design/configs/{HOTSPOT_CONFIG.name}.")
    emit("=" * 72)

    if not write_findings:
        # The self-test run deliberately breaks its own input. Letting it write
        # the findings file would replace a real result with a broken one.
        return 1 if failures else 0

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "09-trimmed-target.md").write_text(
        "# Trimming the target to domain III\n\n"
        "Computed output of `analysis/09_trim_target.py`. Do not hand-edit.\n\n"
        "Cuts the EGFR extracellular region down to domain III and writes the\n"
        "files a binder-design run receives: the trimmed structure, and the\n"
        "hotspot list -- the positions on the target the run is told to aim at --\n"
        "expressed in the trimmed file's own numbering rather than ours.\n\n"
        "EGFR is the epidermal growth factor receptor, the protein being designed\n"
        "against. UniProt is the public sequence archive whose numbering this\n"
        "project uses throughout. The Protein Data Bank, abbreviated PDB, is the\n"
        "public archive of measured three-dimensional structures.\n\n"
        "The cut is made to bring the design run's cost down far enough to afford\n"
        "a few hundred attempts rather than a few. Whether a fragment this size\n"
        "holds the shape its anchors sit on is **not** measured here and is\n"
        "recorded as unverified in `docs/decisions-log.md`. What is measured here\n"
        "is what the cut breaks: which disulfide bonds it severs, and which\n"
        "anchors sit near a cut edge.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
