#!/usr/bin/env python3
"""
10_charge_pair_filter.py — the step that turns generic binders into pH-conditional
ones.

WHAT THIS IS FOR
----------------
The design pipeline has no concept of pH. BindCraft2 reports the binder's charge
at a hard-coded pH 7.4 as a readout it does not act on -- `REPORTED_PH = 7.4` in
its `bindcraft/filters.py`, surfaced as the metrics `Binder_pI` and
`Binder_Net_Charge` -- and has no pH term in its objective. It will return several
hundred well-folded, strongly binding, completely pH-blind candidates. Under this
project's objective ordering, where pH selectivity comes first, that is a failing
submission.

This script is what converts those candidates into pH-conditional ones. It is the
part of the method that is actually ours.

Abbreviations, expanded here because this file gets read on its own:
  EGFR     epidermal growth factor receptor, the protein we design against.
  UniProt  the public protein sequence archive. Every residue number in this
           project is a position in its human record P00533.
  PDB      the Protein Data Bank, the public archive of measured
           three-dimensional structures. Also a file format.
  mmCIF    the newer structure file format, which is what the design pipeline
           writes. Same information, different syntax.
  pH       how acidic something is. Lower means more acidic. Tumour tissue sits
           around 6.5 and blood around 7.4, which is the difference the design
           has to detect.
  aa       amino acids, the building blocks a protein chain is made of.

THE RULE BEING APPLIED
----------------------
Histidine is the only amino acid that changes charge between pH 7.4 and pH 6.5:
neutral above, positive below, because the pH at which half of its copies carry
the extra charge sits around 6.0 to 6.5. So a charge pair that uses a histidine
switches on as the surroundings turn acidic, and one that does not, does not.

  | On the target        | On the binder | Verdict                             |
  |----------------------|---------------|-------------------------------------|
  | D or E, acidic       | histidine     | correct -- switches on as pH falls  |
  | histidine            | D or E        | correct -- switches on as pH falls  |
  | histidine            | histidine     | REJECT THE WHOLE CANDIDATE          |
  | anything else        | anything      | neutral: counted, not scored        |

The rejection applies to **any** target histidine the binder faces, not only the
two in our anchor set. A binder histidine facing a target histidine we never
listed is the same physical problem: both turn positive at pH 6.5 and push apart,
which can cancel a correct pair elsewhere on the same face. That rule has
published experimental support rather than resting on argument -- Liu et al. 2022
tried exactly that arrangement on EGFR's H433 and it changed nothing, while the
acidic pairing improved pH dependence substantially.

A candidate contacting position 442 is also rejected, which is a standing rule in
`docs/decisions-log.md`. S442G is the single human/mouse difference inside the
original epitope and sits in cetuximab's measured contact set, so touching it
risks species-specific behaviour at the one position where the species differ.

WHY PAIR COUNT LEADS, AND CONFIDENCE ONLY BREAKS TIES
-----------------------------------------------------
The charge pairs are what produce the pH switch, so the number of correct pairs
leads the ranking. The pipeline has no interest in pH, which means an interface it
is confident about can carry no switch at all.

The organisers said on 30 September that designs with a large shift in
dissociation constant qualify even when they bind at both pH values, and that
designs with no binding at pH 7.4 and high affinity at pH 6.5 rank higher
(`docs/competition-qa-log.md`). The quantity rewarded is the gap between the two
conditions with the pH 6.5 end as high as the switch allows, so nothing here prefers
weaker binding. An earlier version of this script did; that was withdrawn, see Ruled
out in `docs/decisions-log.md`.

The pipeline's own confidence ordering, `i_pDAE` (higher is better), is carried
through to the output and used as the last term of the ranking. Among candidates
whose pair terms are all equal, the one with the more confident interface comes
first. It never overrides a pair term, it is not folded into a score, and it does
not discard a candidate. The pipeline ranks its own output by `i_pDAE` alone; this
ordering puts the pair terms ahead of it. A candidate with no `i_pDAE` value sorts
after one that has it.

Note what is not claimed. Counting pairs says a switch is built, not that it
works. How far a histidine's flipping point moves depends on its neighbours,
which is much of why pH selectivity resists reliable prediction, and nothing here
measures how large a shift the organisers will count as large, which they have not
said.

TESTING WITH NO CANDIDATES IN EXISTENCE
----------------------------------------
This script was written before any candidate existed, deliberately. Writing the
consumer first means the generation step can be configured to produce exactly
what the consumer needs, rather than discovering after a paid run that a
measurement was never saved.

So the test cases are built by hand and run every time, not only when real output
is around. Each is a small synthetic complex written in the same mmCIF format the
pipeline emits, constructed to exercise one branch: four correct pairs, a binder
histidine facing a listed target histidine, a binder histidine facing a target
histidine we never listed, a candidate contacting 442, a candidate with no correct
pairs, and one whose side chain is unresolved. The expected verdict for each is
asserted.

Outputs:
  results/findings/10-charge-pair-filter.md
  data/derived/10-candidate-pairs.csv     one row per contact pair
  data/derived/10-candidate-summary.csv   one row per candidate
  results/candidates/shortlist.csv        the ranked survivors, the one file
                                          under results/candidates/ that
                                          .gitignore admits

Run standalone:  python analysis/10_charge_pair_filter.py
                 python analysis/10_charge_pair_filter.py --candidates DIR
                 python analysis/10_charge_pair_filter.py --candidates DIR --parse-only

`--candidates` points at a campaign folder from the design run. `--parse-only` is
for the smoke run against the pipeline's own PD-L1 example: the verdicts there are
meaningless, because PD-L1 is not our target, but it proves the script can read
what the pipeline actually emits, which is the part most likely to be wrong.
"""

import argparse
import csv
import sys
import tempfile
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.mmcifio import MMCIFIO

import egfr_common as common
import target_numbering

ROOT = Path(__file__).resolve().parents[1]
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS
CANDIDATES = ROOT / "results" / "candidates"

# Positions no candidate may contact. 442 is the single human/mouse difference
# inside the original epitope and is in cetuximab's measured contact set, so a
# binder touching it risks behaving differently in the two species at the one
# position where they differ. Standing rule in docs/decisions-log.md.
FORBIDDEN_POSITIONS = {442}

# Positions whose pairs are ranked as less supported, without rejecting the
# candidate. E424 is 6.6 angstroms from the C-terminal cut of the trimmed target
# (results/findings/09-trimmed-target.md); every other anchor is 12.3 or further.
# A trimmed model is least reliable at its ends, so a design that leans on E424 is
# leaning on the part of the target most likely to be wrong. A pair there still
# counts as correct for the verdict, but a design that reaches the same count
# without it ranks above one that needs it. A note in a findings file does not
# survive contact with a shortlist, which is why this is a term in the ranking.
EDGE_RELIANT_POSITIONS = {424}

# Positions whose pairs are ranked as less supported for a reason distinct from
# EDGE_RELIANT_POSITIONS above: mouse cross-reactivity rather than reliability of
# the trimmed model. E400, E421 and E424 are the three anchors nearest N361, a
# sugar-chain attachment point (a sequon) that exists in the human sequence and
# not in the mouse one -- mouse has tyrosine there instead
# (results/findings/14-glycan-sequons.md, docs/explainers/08-sugar-chains-near-
# the-anchors.md). The other five anchors -- E344, H358, D368, H370, E391 -- sit
# nearest N352, a sequon present in both species. A sugar chain at N352 costs
# absolute affinity in both species alike and leaves the mouse-to-human K_D ratio
# untouched, which is what mouse cross-reactivity is scored on; a chain at N361
# does not, because whatever shielding it causes happens on the human target and
# not the mouse one, moving that ratio directly.
#
# This is a tie-break built on a plausible asymmetry, not a measured effect, and
# it should be described that way wherever it is cited. Distance to an
# attachment point is a necessary condition for shielding and not evidence of
# it: nothing measures whether a chain at 16-17 A actually reaches these three
# anchors, and the 20-30 A reach the distance bands themselves rest on is marked
# Unverified in `docs/decisions-log.md`. Both target histidines, H358 and H370,
# are on the N352 side, so the half of the pairing rule that needs an acidic
# residue on the binder is unaffected. A pair here still counts as correct for
# the verdict; a design that reaches the same count without it ranks above one
# that needs it, in the same shape as EDGE_RELIANT_POSITIONS. E424 sits in both
# sets, so a design leaning on it is demoted on two independent grounds rather
# than double-counted within one.
N361_SEQUON_RELIANT_POSITIONS = {400, 421, 424}

# How many correct pairs a candidate needs. The switch is partial rather than
# all-or-nothing at pH 6.5 -- only a fraction of a histidine's copies carry the
# extra charge at any moment -- so one pair produces a weak effect and three or
# four have to be stacked. Three is the floor, four is the design target.
MIN_CORRECT_PAIRS = 3
PAIR_TARGET = 4

# The atoms that actually carry the charge. For aspartic and glutamic acid that
# is the two oxygens of the carboxylate group at the end of the side chain; for
# histidine the two nitrogens of the five-membered imidazole ring. A pair is only
# counted as correct when both residues have these atoms resolved, because
# without them there is no evidence the charges reach each other.
CHARGE_ATOMS = {"D": ("OD1", "OD2"), "E": ("OE1", "OE2"),
                "H": ("ND1", "NE2")}
ACIDIC = ("D", "E")

# Two charged groups further apart than this are reported separately. A salt
# bridge -- two opposite charges holding each other -- is usually taken to need
# about 4 angstroms between the charged groups. 6 is deliberately loose, because
# a predicted structure places side chains less reliably than a measured one.
# This is a reported second opinion, not the headline test, and the findings file
# gives both numbers so neither is hidden.
CHARGE_GROUP_REACH = 6.0

# The design pipeline's own interface cutoff, read from its source on 1 October
# 2026: `cutoff: float = 4.0` in bindcraft/filters.py. Ours is 4.5, the usual
# choice in the literature and the one steps 04 and 06 already use. Ours is the
# looser of the two, so our contact set should contain the pipeline's rather than
# equal it, and the cross-check below is written that way.
PIPELINE_CUTOFF = 4.0

# Set by --break-rule, which disables one branch of the rule on purpose so the
# self-test can confirm the test cases above actually catch it. A check that has
# never been seen to fail is not known to work, and that applies to a test case
# as much as to an assertion: a test that would pass with the rule removed is
# testing nothing.
BROKEN_RULE = None

# Each entry is a branch of the rule, the test case that is supposed to catch it
# being removed, and what goes wrong if nothing does.
BREAKABLE_RULES = {
    "his-his": (
        "binder-his-faces-listed-target-his",
        "a binder histidine facing a target histidine would be scored as "
        "harmless, and candidates whose switch cancels itself would reach the "
        "shortlist"),
    "forbidden-position": (
        "contacts-442",
        "a candidate touching position 442 would be kept, risking "
        "species-specific behaviour at the one position where human and mouse "
        "differ inside the epitope"),
    "unresolved-counts": (
        "unresolved-side-chain",
        "a charge pair with no evidence its charged groups exist would be "
        "counted as correct, inflating every pair count"),
    "numbering-reconcile": (
        "renumbered-target",
        "a candidate whose target was renumbered by the design run would be "
        "scored on the returned numbers, giving confident verdicts about the "
        "wrong residues"),
    "edge-reliance": (
        "e424-pair-among-four",
        "a design leaning on E424, 6.6 angstroms from the cut in a trimmed "
        "model, would rank level with an equivalent design that does not"),
    "sequon-reliance": (
        "e421-pair-among-four",
        "a design leaning on E400, E421 or E424 -- the three anchors nearest "
        "the human-only N361 sequon -- would rank level with an equivalent "
        "design that leans on the five anchors nearest the shared N352 "
        "sequon instead, despite the asymmetric risk to the mouse "
        "cross-reactivity ratio"),
    "pdae-tiebreak": (
        "pdae tie-break",
        "candidates with identical pair counts would be ordered by name, so the "
        "more confident interface would not be preferred among equals"),
}

VERDICT_REJECTED = "rejected"
VERDICT_MEETS = "meets the pair target"
VERDICT_BELOW = "below the pair target"
VERDICT_UNSCORABLE = "not scored"

# The structure file the design run was handed, which is what a returned target
# chain is compared against. Written by analysis/09; not tracked, because it is
# regenerable.
INPUT_TARGET = common.STRUCT_DIR / "6aru_domain3.pdb"


# ---------------------------------------------------------------------------
# the rule itself
# ---------------------------------------------------------------------------

def classify(target_aa, binder_aa):
    """One contact pair, classified by the rule in CLAUDE.md.

    Returns one of "correct", "his-his", "neutral".
    """
    if target_aa == "H" and binder_aa == "H":
        if BROKEN_RULE == "his-his":
            return "neutral"
        return "his-his"
    if target_aa in ACIDIC and binder_aa == "H":
        return "correct"
    if target_aa == "H" and binder_aa in ACIDIC:
        return "correct"
    return "neutral"


def charge_group_distance(residue_a, aa_a, residue_b, aa_b):
    """Closest approach between the two residues' charge-carrying atoms.

    Returns (distance, both resolved). When either side's charge atoms are
    missing from the structure the distance is None and the flag is False, and
    the pair is never counted as correct on that basis -- an unresolved side
    chain is a thing we do not know, not a thing we know to be fine.
    """
    wanted_a = CHARGE_ATOMS.get(aa_a)
    wanted_b = CHARGE_ATOMS.get(aa_b)
    if not wanted_a or not wanted_b:
        return None, False
    atoms_a = [residue_a[name] for name in wanted_a if name in residue_a]
    atoms_b = [residue_b[name] for name in wanted_b if name in residue_b]
    if len(atoms_a) != len(wanted_a) or len(atoms_b) != len(wanted_b):
        return None, False
    best = min(float(np.linalg.norm(a.coord - b.coord))
               for a in atoms_a for b in atoms_b)
    return best, True


# ---------------------------------------------------------------------------
# reading one candidate complex
# ---------------------------------------------------------------------------

def load_complex(path):
    """Parse a candidate structure, via the one shared parser.

    Moved to `egfr_common` when step 13 came to need the same thing: step 13 reads
    the very same candidate files, and two parsers could disagree about one file.
    """
    return common.load_complex(path)


def identify_chains(model, numbering, human, emit):
    """Which chain is the target and which is the designed binder, by measurement.

    The rule and the reasoning now live in `egfr_common.identify_chains`, because
    step 13 has to reach the same verdict about the same file as this step does.
    """
    return common.identify_chains(model, numbering, human, emit)


def unscored(name, path, reason, numbering_status, numbering_detail):
    """The result for a candidate this script refuses to score. Every field a
    scored result has is present and zero, so nothing downstream has to special
    case it, and the reason is carried where a reader will see it."""
    return dict(
        design=name, path=path, pairs=[], verdict=VERDICT_UNSCORABLE,
        reasons=[reason], correct_pairs=0, correct_pairs_tight=0,
        edge_reliant_pairs=0, supported_pairs=0,
        n361_reliant_pairs=0, anchor_group_supported_pairs=0,
        unresolved_pairs=0,
        his_his_pairs=0, forbidden_contacts=0, neutral_pairs=0, total_pairs=0,
        target_positions=[], anchors_paired=[],
        numbering_status=numbering_status, numbering_detail=numbering_detail)


def analyse_candidate(name, path, numbering, human, target_chain=None,
                      binder_chains=None, translate=True, input_target=None):
    """Every contact pair in one candidate, classified, and the verdict.

    The contact calculation is not implemented here. It is
    `egfr_common.contact_pairs`, which steps 04 and 06 also reach, at the same
    4.5 angstrom cutoff. Steps 04 and 06 once each kept their own copy of that
    calculation and disagreed for months, which is why there is now one copy.
    """
    model = load_complex(path)

    # Before anything is scored: is the target in this file numbered the way the
    # structure the run was handed is numbered? Everything below translates the
    # returned numbers into ours by an offset measured for the input, so a target
    # the run renumbered would be read as different residues, and the verdicts
    # would come out confident and wrong. A candidate that cannot be reconciled is
    # refused, not scored on the assumption that nothing moved.
    numbering_status, numbering_detail = "not checked", ""
    if input_target is not None and translate and BROKEN_RULE != "numbering-reconcile":
        check = target_numbering.check_complex(input_target, path)
        numbering_status = check["status"]
        numbering_detail = check["detail"]
        if not check["ok"]:
            return unscored(
                name, path,
                "target numbering cannot be reconciled with the input: "
                + target_numbering.describe(check),
                numbering_status, numbering_detail)

    pairs = common.contact_pairs(
        model, target_chain, binder_chains,
        numbering=numbering if translate else None)

    target_residues = {r.id[1]: r for r in common.protein_residues(
        model[target_chain])}
    binder_residues = {}
    for cid in binder_chains:
        for residue in common.protein_residues(model[cid]):
            binder_residues[(cid, residue.id[1])] = residue

    rows = []
    for pair in pairs:
        target_aa, binder_aa = pair["receptor_aa"], pair["partner_aa"]
        verdict = classify(target_aa, binder_aa)
        distance, resolved = charge_group_distance(
            target_residues[pair["receptor_resnum"]], target_aa,
            binder_residues[(pair["partner_chain"], pair["partner_resnum"])],
            binder_aa)
        rows.append(dict(
            design=name,
            target_pos=pair["receptor_pos"],
            target_resnum=pair["receptor_resnum"],
            target_aa=target_aa,
            binder_chain=pair["partner_chain"],
            binder_pos=pair["partner_resnum"],
            binder_aa=binder_aa,
            min_dist=pair["min_dist"],
            atom_pairs=pair["atom_pairs"],
            classification=verdict,
            charge_group_dist=distance,
            charge_atoms_resolved=resolved,
            charge_groups_reach=(resolved and distance <= CHARGE_GROUP_REACH),
            forbidden=(translate and BROKEN_RULE != "forbidden-position"
                       and pair["receptor_pos"] in FORBIDDEN_POSITIONS),
        ))

    if BROKEN_RULE == "unresolved-counts":
        # The broken version counts a pair whose charge atoms are missing as if
        # they had been seen. This is what the unresolved-side-chain test case
        # exists to catch.
        correct = [r for r in rows if r["classification"] == "correct"]
        unresolved = []
    else:
        correct = [r for r in rows
                   if r["classification"] == "correct"
                   and r["charge_atoms_resolved"]]
        unresolved = [r for r in rows
                      if r["classification"] == "correct"
                      and not r["charge_atoms_resolved"]]
    his_his = [r for r in rows if r["classification"] == "his-his"]
    forbidden = [r for r in rows if r["forbidden"]]
    tight = [r for r in correct if r["charge_groups_reach"]]
    edge_reliant = ([] if BROKEN_RULE == "edge-reliance"
                    else [r for r in correct
                          if r["target_pos"] in EDGE_RELIANT_POSITIONS])
    n361_reliant = ([] if BROKEN_RULE == "sequon-reliance"
                    else [r for r in correct
                          if r["target_pos"] in N361_SEQUON_RELIANT_POSITIONS])

    reasons = []
    if his_his:
        reasons.append(
            "binder histidine faces target histidine "
            + ", ".join(f"H{r['target_pos']}" for r in sorted(
                his_his, key=lambda r: r["target_pos"])))
    if forbidden:
        reasons.append(
            "contacts "
            + ", ".join(str(r["target_pos"]) for r in sorted(
                forbidden, key=lambda r: r["target_pos"])))

    if reasons:
        verdict = VERDICT_REJECTED
    elif len(correct) >= MIN_CORRECT_PAIRS:
        verdict = VERDICT_MEETS
    else:
        verdict = VERDICT_BELOW

    return dict(
        design=name, path=path, pairs=rows, verdict=verdict,
        reasons=reasons, correct_pairs=len(correct),
        correct_pairs_tight=len(tight),
        edge_reliant_pairs=len(edge_reliant),
        supported_pairs=len(correct) - len(edge_reliant),
        n361_reliant_pairs=len(n361_reliant),
        anchor_group_supported_pairs=len(correct) - len(n361_reliant),
        unresolved_pairs=len(unresolved), his_his_pairs=len(his_his),
        forbidden_contacts=len(forbidden),
        neutral_pairs=sum(1 for r in rows if r["classification"] == "neutral"),
        total_pairs=len(rows),
        target_positions=sorted({r["target_pos"] for r in rows}),
        anchors_paired=sorted({r["target_pos"] for r in correct}),
        numbering_status=numbering_status, numbering_detail=numbering_detail,
    )


# ---------------------------------------------------------------------------
# building the test cases
# ---------------------------------------------------------------------------
#
# Each synthetic complex is two chains written in the same mmCIF format the
# design pipeline emits. The target chain holds real EGFR residues at real
# positions, numbered the way the trimmed target is numbered, so the test also
# exercises the translation back into our numbering rather than only the
# classification. The binder chain holds whatever residue the case needs,
# numbered from 1.
#
# Geometry is built rather than taken from anywhere: each target residue sits on
# a line, its side chain pointing one way, and its partner sits facing it with
# the gap set so the two side-chain tips are 3 angstroms apart. Target residues
# are spaced far enough apart that no binder residue reaches a neighbour it was
# not meant to face, which is checked by the pair counts the tests assert.

# (atom name, element, how far out along the side chain, sideways offset)
BACKBONE = [("N", "N", -1.2, -1.2), ("CA", "C", 0.0, 0.0),
            ("C", "C", -1.2, 1.2), ("O", "O", -1.6, 1.8)]
SIDE_CHAINS = {
    "D": [("CB", "C", 1.0, 0.0), ("CG", "C", 2.0, 0.0),
          ("OD1", "O", 2.8, 0.6), ("OD2", "O", 2.8, -0.6)],
    "E": [("CB", "C", 1.0, 0.0), ("CG", "C", 2.0, 0.0), ("CD", "C", 2.8, 0.0),
          ("OE1", "O", 3.6, 0.6), ("OE2", "O", 3.6, -0.6)],
    "H": [("CB", "C", 1.0, 0.0), ("CG", "C", 2.0, 0.0),
          ("ND1", "N", 2.8, 0.7), ("CD2", "C", 2.8, -0.7),
          ("CE1", "C", 3.6, 0.5), ("NE2", "N", 3.6, -0.5)],
    "S": [("CB", "C", 1.0, 0.0), ("OG", "O", 2.0, 0.0)],
    "T": [("CB", "C", 1.0, 0.0), ("OG1", "O", 2.0, 0.5),
          ("CG2", "C", 2.0, -0.5)],
    "L": [("CB", "C", 1.0, 0.0), ("CG", "C", 2.0, 0.0),
          ("CD1", "C", 2.8, 0.7), ("CD2", "C", 2.8, -0.7)],
    "I": [("CB", "C", 1.0, 0.0), ("CG1", "C", 2.0, 0.6),
          ("CG2", "C", 2.0, -0.6), ("CD1", "C", 2.8, 0.9)],
    "A": [("CB", "C", 1.0, 0.0)],
    "V": [("CB", "C", 1.0, 0.0), ("CG1", "C", 2.0, 0.6),
          ("CG2", "C", 2.0, -0.6)],
}
ONE_TO_THREE = {"A": "ALA", "D": "ASP", "E": "GLU", "H": "HIS", "I": "ILE",
                "L": "LEU", "S": "SER", "T": "THR", "V": "VAL"}

TARGET_SPACING = 14.0      # angstroms between neighbouring target residues
TIP_SEPARATION = 3.0       # angstroms between the two facing side-chain tips


def side_chain_for(aa, with_side_chain=True):
    if not with_side_chain:
        return [("CB", "C", 1.0, 0.0)]
    return SIDE_CHAINS.get(aa, [("CB", "C", 1.0, 0.0)])


def reach_of(atoms):
    return max(depth for _name, _el, depth, _lat in atoms)


def build_residue(resnum, aa, atoms, x, y_anchor, direction, serial):
    residue = Residue((" ", resnum, " "), ONE_TO_THREE[aa], " ")
    for name, element, depth, lateral in atoms:
        coord = np.array([x + lateral, y_anchor + direction * depth, 0.0],
                         dtype=float)
        residue.add(Atom(name, coord, 30.0, 1.0, " ", name, serial, element))
        serial += 1
    return residue, serial


def build_test_complex(path, contacts, human, target_shift=0):
    """Write one synthetic candidate complex.

    `contacts` is a list of (target uniprot position, binder amino acid,
    binder has its side chain). The target residue's own amino acid is read out
    of the human sequence rather than supplied, so a test cannot assert a verdict
    against a residue identity that is not real.
    """
    structure = Structure("test")
    model = Model(0)
    structure.add(model)
    target_chain, binder_chain = Chain("A"), Chain("B")
    model.add(target_chain)
    model.add(binder_chain)

    serial = 1
    for index, (pos, binder_aa, binder_side_chain) in enumerate(contacts):
        target_aa = human[pos - 1]
        x = index * TARGET_SPACING
        target_atoms = BACKBONE + side_chain_for(target_aa)
        binder_atoms = BACKBONE + side_chain_for(binder_aa, binder_side_chain)
        gap = reach_of(target_atoms) + reach_of(binder_atoms) + TIP_SEPARATION

        # The target keeps the numbering the trimmed target file has, which is
        # ours minus the offset analysis/02 measured, so the translation back is
        # exercised rather than bypassed.
        # `target_shift` renumbers the target on purpose, which is how a design
        # run that renumbered its target would look.
        target_resnum = pos - common.SIGNAL_PEPTIDE_LEN + target_shift
        residue, serial = build_residue(target_resnum, target_aa, target_atoms,
                                        x, 0.0, +1.0, serial)
        target_chain.add(residue)
        residue, serial = build_residue(index + 1, binder_aa, binder_atoms,
                                        x, gap, -1.0, serial)
        binder_chain.add(residue)

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    return path


def test_cases(human):
    """Every branch of the rule, one case each, with the verdict asserted.

    `contacts` entries are (target position, binder amino acid, binder side chain
    resolved). Expected values are written out rather than computed, so a change
    in the rule shows up as a failing test rather than as a quietly different
    answer.
    """
    return [
        dict(
            name="four-correct-pairs",
            why="The design target: four pairs that all switch the right way. "
                "Two acidic target residues faced by binder histidines, and two "
                "target histidines faced by binder acidic residues.",
            contacts=[(344, "H", True), (368, "H", True),
                      (370, "E", True), (358, "D", True)],
            expect=dict(verdict=VERDICT_MEETS, correct_pairs=4,
                        his_his_pairs=0, forbidden_contacts=0,
                        unresolved_pairs=0),
        ),
        dict(
            name="binder-his-faces-listed-target-his",
            why="A binder histidine facing H370, one of the two target "
                "histidines in our anchor set. Both turn positive at pH 6.5 and "
                "push apart, which can cancel a correct pair elsewhere on the "
                "same face. Rejected rather than scored, even though two correct "
                "pairs are present.",
            contacts=[(344, "H", True), (368, "H", True), (370, "H", True)],
            expect=dict(verdict=VERDICT_REJECTED, correct_pairs=2,
                        his_his_pairs=1, forbidden_contacts=0,
                        unresolved_pairs=0),
        ),
        dict(
            name="binder-his-faces-unlisted-target-his",
            why="The same fault against H418, a target histidine that is not in "
                "our anchor set at all. The rejection has to apply to any target "
                "histidine, not only the ones we listed, because the physical "
                "problem is identical.",
            contacts=[(344, "H", True), (418, "H", True)],
            expect=dict(verdict=VERDICT_REJECTED, correct_pairs=1,
                        his_his_pairs=1, forbidden_contacts=0,
                        unresolved_pairs=0),
        ),
        dict(
            name="contacts-442",
            why="Four correct pairs, and one contact at position 442. That is "
                "the single human/mouse difference inside the original epitope "
                "and sits in cetuximab's contact set, so the candidate is "
                "rejected despite having everything else right. This is the case "
                "that proves a good pair count cannot buy its way past a hard "
                "rule.",
            contacts=[(344, "H", True), (368, "H", True), (391, "H", True),
                      (400, "H", True), (442, "A", True)],
            expect=dict(verdict=VERDICT_REJECTED, correct_pairs=4,
                        his_his_pairs=0, forbidden_contacts=1,
                        unresolved_pairs=0),
        ),
        dict(
            name="no-correct-pairs",
            why="A well-formed interface with no charge pair anywhere in it. "
                "Not rejected, because nothing forbidden happens: it is kept, "
                "reported, and ranked last. A candidate is discarded for "
                "breaking a rule, never for being weak.",
            contacts=[(365, "A", True), (369, "L", True), (371, "V", True),
                      (372, "A", True)],
            expect=dict(verdict=VERDICT_BELOW, correct_pairs=0,
                        his_his_pairs=0, forbidden_contacts=0,
                        unresolved_pairs=0),
        ),
        dict(
            name="unresolved-side-chain",
            why="Three pairs that would all be correct, but one binder "
                "histidine has no side chain in the structure, so there is no "
                "evidence its charge reaches anything. It is counted separately "
                "and not as correct, which drops the candidate from three "
                "correct pairs to two and below the target. If an unresolved "
                "side chain were quietly counted, this case would pass.",
            contacts=[(344, "H", True), (391, "H", True), (368, "H", False)],
            expect=dict(verdict=VERDICT_BELOW, correct_pairs=2,
                        his_his_pairs=0, forbidden_contacts=0,
                        unresolved_pairs=1),
        ),
        dict(
            name="e424-pair-among-four",
            why="Four correct pairs, one of them on E424, which sits 6.6 "
                "angstroms from the cut in the trimmed target. Still meets the "
                "pair target, because the pair is correct, but only three of "
                "the four count as supported, so it ranks below the "
                "four-correct-pairs case that reaches the same count without "
                "E424. E424 is also one of the three anchors nearest N361, so "
                "it is demoted on two independent grounds at once, not "
                "double-counted within either field.",
            contacts=[(344, "H", True), (368, "H", True),
                      (370, "E", True), (424, "H", True)],
            expect=dict(verdict=VERDICT_MEETS, correct_pairs=4,
                        edge_reliant_pairs=1, supported_pairs=3,
                        n361_reliant_pairs=1, anchor_group_supported_pairs=3,
                        his_his_pairs=0, forbidden_contacts=0,
                        unresolved_pairs=0),
        ),
        dict(
            name="e421-pair-among-four",
            why="Four correct pairs, one of them on E421, one of the three "
                "anchors nearest N361, a sugar-chain attachment point that "
                "exists in human and not in mouse. Still meets the pair "
                "target, because the pair is correct, but only three of the "
                "four count toward the anchor-group-supported total, so it "
                "ranks below the four-correct-pairs case that reaches the "
                "same count using only anchors nearest the shared N352 "
                "sequon. E421 is not near the trimmed target's cut edge, so "
                "this isolates the new demotion from the E424 edge-reliance "
                "one above. Named to sort alphabetically before "
                "four-correct-pairs, so the ranking assertion below is a "
                "real test of the demotion term rather than a pass that "
                "would happen anyway from the name-based final tie-break.",
            contacts=[(344, "H", True), (368, "H", True),
                      (370, "E", True), (421, "H", True)],
            expect=dict(verdict=VERDICT_MEETS, correct_pairs=4,
                        edge_reliant_pairs=0, supported_pairs=4,
                        n361_reliant_pairs=1, anchor_group_supported_pairs=3,
                        his_his_pairs=0, forbidden_contacts=0,
                        unresolved_pairs=0),
        ),
        dict(
            name="renumbered-target",
            why="Four correct pairs, but the target in the returned file is "
                "numbered one place off from the structure the run was handed, "
                "which is what a design run that cropped and renumbered its "
                "target would return. Every pair would still look right, read "
                "against the wrong residues, so the candidate is refused rather "
                "than scored.",
            contacts=[(344, "H", True), (368, "H", True),
                      (370, "E", True), (358, "D", True)],
            target_shift=1,
            expect=dict(verdict=VERDICT_UNSCORABLE, correct_pairs=0,
                        numbering_status=target_numbering.MISMATCH),
        ),
    ]


def input_target_path():
    """The structure the design run is handed, which returned targets are checked
    against. Stops with the way to make it if it is missing, because a check that
    quietly skips itself when its reference is absent protects nothing."""
    if not INPUT_TARGET.exists():
        raise SystemExit(f"{INPUT_TARGET} is missing. Run "
                         f"analysis/09_trim_target.py to write it.")
    return INPUT_TARGET


def run_tests(numbering, human, emit):
    """Build every test case, classify it, and compare against the expectation."""
    emit("1. The test cases, built by hand and run every time")
    emit()
    emit("   No candidate exists yet. These are synthetic complexes written in")
    emit("   the same mmCIF format the design pipeline emits, each built to")
    emit("   exercise one branch of the rule. The target residues are real EGFR")
    emit("   residues at real positions, numbered the way the trimmed target is")
    emit("   numbered, so the translation back into our numbering is exercised")
    emit("   too rather than bypassed.")
    emit()
    failures = []
    results = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for case in test_cases(human):
            path = build_test_complex(
                Path(tmpdir) / f"{case['name']}.cif", case["contacts"], human,
                target_shift=case.get("target_shift", 0))
            result = analyse_candidate(case["name"], path, numbering, human,
                                       target_chain="A", binder_chains=["B"],
                                       input_target=input_target_path())
            results.append((case, result))
            emit(f"   {case['name']}")
            emit(f"     {case['why']}")
            emit("     contacts: " + ", ".join(
                f"{human[p - 1]}{p} faced by binder {b}"
                + ("" if side else " (side chain unresolved)")
                for p, b, side in case["contacts"]))
            for field, want in sorted(case["expect"].items()):
                got = result[field]
                ok = got == want
                emit(f"     [{'PASS' if ok else 'FAIL'}] {field}: {got}"
                     + ("" if ok else f", expected {want}"))
                if not ok:
                    failures.append(f"{case['name']}: {field} was {got}, "
                                    f"expected {want}")
            if result["reasons"]:
                emit("     reasons given: " + "; ".join(result["reasons"]))
            emit()
    by_name = {case["name"]: result for case, result in results}
    ordered = [s["design"] for s in rank(
        [by_name["e424-pair-among-four"], by_name["four-correct-pairs"]])]
    ok = ordered[0] == "four-correct-pairs"
    emit("   Ranking: four correct pairs without E424 against four with it")
    emit(f"     [{'PASS' if ok else 'FAIL'}] order: {' then '.join(ordered)}")
    if not ok:
        failures.append("e424 ranking: the design leaning on E424 ranked "
                        "level with or above the equivalent one that does not")
    emit()

    ordered = [s["design"] for s in rank(
        [by_name["e421-pair-among-four"], by_name["four-correct-pairs"]])]
    ok = ordered[0] == "four-correct-pairs"
    emit("   Ranking: four correct pairs on the N352 side against four with "
         "one on the N361 side")
    emit(f"     [{'PASS' if ok else 'FAIL'}] order: {' then '.join(ordered)}")
    if not ok:
        failures.append("n361 ranking: the design leaning on E400/E421/E424 "
                        "ranked level with or above the equivalent one "
                        "leaning on the N352-side anchors instead")
    emit()

    # The tie-break on the pipeline's own i_pDAE. The names are chosen so that
    # ordering by name alone would put them in the wrong order, which is what makes
    # the case fail when the term is switched off.
    base = by_name["four-correct-pairs"]
    level = [dict(base, design="tie-a-doubtful", metrics={"i_pDAE": "0.41"}),
             dict(base, design="tie-m-no-metric"),
             dict(base, design="tie-z-confident", metrics={"i_pDAE": "0.82"})]
    ordered = [s["design"] for s in rank(level)]
    want = ["tie-z-confident", "tie-a-doubtful", "tie-m-no-metric"]
    ok = ordered == want
    emit("   Ranking: three candidates level on every pair term, ordered by i_pDAE")
    emit(f"     [{'PASS' if ok else 'FAIL'}] order: {' then '.join(ordered)}")
    if not ok:
        failures.append("pdae tie-break: level candidates were not ordered with "
                        "the more confident interface first and the one with no "
                        f"value last (got {ordered})")
    emit()

    fewer = dict(base, design="pairs-fewer-confident", supported_pairs=3,
                 correct_pairs=3, correct_pairs_tight=3,
                 metrics={"i_pDAE": "0.95"})
    more = dict(base, design="pairs-more-doubtful", metrics={"i_pDAE": "0.30"})
    ordered = [s["design"] for s in rank([fewer, more])]
    ok = ordered[0] == "pairs-more-doubtful"
    emit("   Ranking: more pairs but a worse i_pDAE against fewer pairs and a better one")
    emit(f"     [{'PASS' if ok else 'FAIL'}] order: {' then '.join(ordered)}")
    if not ok:
        failures.append("pair count leads: a more confident interface outranked "
                        "a candidate with more correct pairs")
    emit()
    emit(f"   {len(test_cases(human))} cases, "
         f"{'all passed' if not failures else f'{len(failures)} assertion(s) failed'}.")
    emit()
    if failures:
        emit("   The rule and the tests disagree. Do not run this script against")
        emit("   real candidates until they agree: a filter that is wrong about")
        emit("   a synthetic case it was handed is wrong about a real one it was")
        emit("   not.")
        emit()
    return failures, results


# ---------------------------------------------------------------------------
# reading a real campaign folder
# ---------------------------------------------------------------------------

def find_candidate_structures(folder):
    """The accepted-design complexes in a campaign folder.

    The design pipeline writes them to `3_Ranked/<design>_seq<n>[_<target>].cif`,
    alongside `<design>_seq<n>_monomer.cif`, which is the binder on its own and
    has no interface to measure. Read from its own output documentation on
    1 October 2026.

    If that folder is not there, every mmCIF below the path is taken instead and
    the run says so, because an output layout that has moved should produce a
    visible fallback rather than an empty result.
    """
    ranked = folder / "3_Ranked"
    if ranked.is_dir():
        found = sorted(p for p in ranked.glob("*.cif")
                       if not p.stem.endswith("_monomer"))
        if found:
            return found, "3_Ranked/"
    found = sorted(p for p in folder.rglob("*.cif")
                   if not p.stem.endswith("_monomer"))
    return found, "every mmCIF below the folder (3_Ranked/ not found)"


def read_metrics_table(folder, emit):
    """The pipeline's own per-design metrics, carried through untouched.

    Looks for `3_Ranked/!_Ranked.csv`, which the pipeline ranks best-first by
    `i_pDAE`. Every column is carried through to our output as-is. We do not
    reorder by any of them and we do not recompute any of them: they are the
    pipeline's answers, kept so a reader can see them beside the pair counts.
    The one exception to not reordering is `i_pDAE`, which `rank` uses as its
    final tie-break.

    Returns ({design name: {column: value}}, column names, where it came from).
    """
    for candidate in (folder / "3_Ranked" / "!_Ranked.csv",
                      folder / "3_Ranked" / "!_Refolded.csv",
                      folder / "summary.csv"):
        if candidate.is_file():
            break
    else:
        matches = sorted(folder.rglob("*.csv"))
        candidate = matches[0] if matches else None
    if candidate is None or not candidate.is_file():
        emit("   No metrics table found. The pipeline's own numbers will be")
        emit("   absent from the output, so the i_pDAE tie-break has nothing to")
        emit("   work with and the ranking rests on the pair terms alone.")
        return {}, [], None

    with candidate.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {}, [], candidate
    columns = list(rows[0].keys())
    key = next((c for c in ("design", "Design", "name") if c in columns), None)
    if key is None:
        emit(f"   {candidate.name} has no design-name column "
             f"({columns[:6]}...), so its numbers cannot be joined to the")
        emit("   structures. Carried through as nothing rather than guessed.")
        return {}, columns, candidate
    return ({row[key]: row for row in rows if row.get(key)},
            columns, candidate)


def join_metrics(name, metrics):
    """Find the metrics row for a structure file.

    The structure is named `<design>_seq<n>[_<target>].cif` and the table's design
    column may hold either the whole thing or the part before the target suffix,
    so an exact match is tried first and then the longest key that the file name
    starts with. A wrong join would attach one candidate's numbers to another, so
    an ambiguous one returns nothing rather than a guess.
    """
    if name in metrics:
        return metrics[name]
    prefixes = sorted((k for k in metrics if name.startswith(k)), key=len,
                      reverse=True)
    return metrics[prefixes[0]] if prefixes else {}


def cross_check_interface(result, row, emit):
    """Our contact set against the pipeline's own, for the same candidate.

    The pipeline reports `Interface_Target_Residues` as a comma-separated list of
    one-letter-plus-number tokens, measured at a 4.0 angstrom cutoff. Ours is
    measured at 4.5, the cutoff steps 04 and 06 use. Ours is the looser of the
    two, so our set should contain the pipeline's rather than equal it, and a
    residue the pipeline reports that we do not is a real disagreement worth
    stopping on.

    Returns (verdict word, detail), or None when the column is absent.
    """
    raw = row.get("Interface_Target_Residues")
    if not raw:
        return None
    theirs = set()
    for token in raw.replace("/", ",").split(","):
        token = token.strip()
        if len(token) > 1 and token[1:].lstrip("-").isdigit():
            theirs.add(int(token[1:]))
    if not theirs:
        return None
    # Their numbers are the target structure's own, so translate ours back.
    ours = {r["target_resnum"] for r in result["pairs"]}
    missing = sorted(theirs - ours)
    extra = sorted(ours - theirs)
    if not missing:
        return ("contains", f"{len(ours)} residues at 4.5 A, containing all "
                            f"{len(theirs)} the pipeline reports at "
                            f"{PIPELINE_CUTOFF} A"
                            + (f", plus {len(extra)} the looser cutoff adds"
                               if extra else ""))
    return ("DISAGREES", f"the pipeline reports {missing} as interface residues "
                         f"and we do not, although our cutoff is the looser one")


# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Score design candidates by pH-switching charge pairs.")
    parser.add_argument("--candidates", type=Path, default=None, metavar="DIR",
                        help="A campaign folder from the design run.")
    parser.add_argument("--parse-only", action="store_true",
                        help="Report contacts without translating numbering or "
                             "issuing verdicts. For the smoke run against the "
                             "pipeline's own example target, where the verdicts "
                             "would be meaningless but the parsing is the point.")
    parser.add_argument("--break-rule", choices=list(BREAKABLE_RULES),
                        default=None,
                        help="Disable one branch of the rule on purpose, so the "
                             "test case covering it is seen to fail. Writes no "
                             "output files.")
    parser.add_argument("--self-test", action="store_true",
                        help="Break each branch of the rule in turn and confirm "
                             "the test cases catch it. A test that would pass "
                             "with the rule removed is testing nothing.")
    args = parser.parse_args(argv)

    if args.self_test:
        return run_self_test()

    global BROKEN_RULE
    BROKEN_RULE = args.break_rule

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    emit("=" * 72)
    emit("CHARGE-PAIR FILTER")
    emit("=" * 72)
    emit()
    emit("Scores design candidates by how many of their contacts form a charge")
    emit("pair that switches on as the surroundings turn acidic, and rejects the")
    emit("ones that break a hard rule. Residue numbers are positions in the")
    emit("human record P00533 in UniProt, the public sequence archive, unless a")
    emit("line says they are the structure file's own.")
    emit()
    emit("Ranked by correct pairs first, because the pairs are what produce the pH")
    emit("switch. Two demotions never discard a candidate but rank it below an")
    emit("equivalent one that does not need the same anchor: leaning on E424, near")
    emit("the trimmed target's cut edge, and leaning on E400, E421 or E424, the")
    emit("three anchors nearest the human-only N361 sugar-chain attachment point.")
    emit("The pipeline's own i_pDAE (higher is better) is carried through and")
    emit("breaks ties among candidates whose pair terms are all equal. Neither")
    emit("demotion nor i_pDAE overrides a pair term and neither discards a")
    emit("candidate.")
    emit()

    human = common.human_sequence()
    numbering = common.load_numbering()

    failures, _test_results = run_tests(numbering, human, emit)
    failures.extend(run_campaign_test(numbering, human, emit))

    summaries = []
    all_pairs = []
    if args.candidates is None:
        emit("3. Real candidates")
        emit()
        emit("   None given. Run again with --candidates pointing at a campaign")
        emit("   folder once the design run has produced one. The tests above")
        emit("   are what can be checked before then, and they check the rule")
        emit("   rather than the pipeline's output format.")
        emit()
    else:
        folder = args.candidates
        emit("3. Real candidates")
        emit()
        if not folder.is_dir():
            emit(f"   {folder} is not a directory.")
            failures.append(f"{folder} is not a directory")
        else:
            summaries, all_pairs, found_failures = process_campaign(
                folder, numbering, human, args.parse_only, emit)
            failures.extend(found_failures)

    if BROKEN_RULE is None:
        write_outputs(summaries, all_pairs, emit,
                      parse_only=args.parse_only, had_candidates=bool(summaries))
    else:
        emit(f"   Running with --break-rule {BROKEN_RULE}, so no file is")
        emit("   written. The point of this run is the failures above.")
        emit()

    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} FAILURE(S).")
        for failure in failures:
            emit(f"  - {failure}")
    elif summaries:
        meets = sum(1 for s in summaries if s["verdict"] == VERDICT_MEETS)
        rejected = sum(1 for s in summaries if s["verdict"] == VERDICT_REJECTED)
        refused = sum(1 for s in summaries if s["verdict"] == VERDICT_UNSCORABLE)
        emit(f"RESULT: PASSED. {len(summaries)} candidates: {meets} meet the "
             f"pair target, {rejected} rejected, {refused} not scored "
             f"(numbering), {len(summaries) - meets - rejected - refused} below "
             f"the target and kept.")
    else:
        emit("RESULT: PASSED. The rule is exercised and correct on every test "
             "case; no real candidates were supplied.")
    emit("=" * 72)

    if BROKEN_RULE is not None:
        # A deliberately broken run must not replace a real findings file.
        return 1 if failures else 0

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "10-charge-pair-filter.md").write_text(
        "# The charge-pair filter\n\n"
        "Computed output of `analysis/10_charge_pair_filter.py`. "
        "Do not hand-edit.\n\n"
        "Takes the candidate complexes a binder-design run produces, works out\n"
        "which target residue each binder contact position faces, and scores each\n"
        "candidate by how many of those pairs switch on as the surroundings turn\n"
        "acidic. The design pipeline has no pH term in its objective, so this is\n"
        "the step that makes the submission pH-conditional rather than generic.\n\n"
        "EGFR is the epidermal growth factor receptor, the protein being designed\n"
        "against. UniProt is the public sequence archive whose numbering this\n"
        "project uses. mmCIF is the structure file format the design pipeline\n"
        "writes.\n\n"
        "**Ranked by correct pairs first.** The pairs are what produce the pH\n"
        "switch, and the design pipeline has no interest in pH. It ranks its own\n"
        "output by `i_pDAE`, a measure of how confident it is in the interface;\n"
        "here that value is carried through and used only as the last tie-break,\n"
        "among candidates whose pair terms are all equal, with the more confident\n"
        "interface first. It never overrides a pair term.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 1 if failures else 0


def build_test_campaign(folder, human):
    """A synthetic campaign folder laid out the way the design pipeline lays one.

    Covers the parts most likely to break when real output first arrives, none of
    which the classification tests touch: finding the accepted complexes among
    the binder-only ones, reading the metrics table, joining a metrics row to a
    structure by name, identifying which chain is the target by measurement, and
    comparing our contact set against the pipeline's own.

    The layout and the column names were read out of the pipeline's own source
    and output documentation on 1 October 2026, not recalled. If the real output
    differs from this, that is the thing to fix, and this is where it shows.
    """
    ranked = folder / "3_Ranked"
    ranked.mkdir(parents=True, exist_ok=True)

    designs = [
        ("design_a1b2c3_seq1", [(344, "H", True), (368, "H", True),
                                (370, "E", True), (358, "D", True)]),
        ("design_d4e5f6_seq1", [(344, "H", True), (370, "H", True)]),
    ]
    rows = []
    for name, contacts in designs:
        build_test_complex(ranked / f"{name}_EGFR_domain3.cif", contacts, human)
        # The binder on its own, which the pipeline also writes and which has no
        # interface to measure. It must be skipped rather than scored.
        build_test_complex(ranked / f"{name}_monomer.cif", contacts[:1], human)
        interface = ",".join(
            f"{human[pos - 1]}{pos - common.SIGNAL_PEPTIDE_LEN}"
            for pos, _b, _s in contacts)
        rows.append({
            "rank": len(rows) + 1, "trajectory": name.split("_seq")[0],
            "design": name, "length": 72, "outcome": "accepted",
            "i_pDAE": 0.21, "i_pTM": 0.83, "pLDDT": 0.91, "pTM": 0.78,
            "i_pAE": 0.24, "Unbound_Binder_pLDDT": 0.88, "Target_pLDDT": 0.94,
            "Binder_Sequence": "H" * 72,
            "Interface_Binder_Residues": ",".join(
                f"{b}{i + 1}" for i, (_p, b, _s) in enumerate(contacts)),
            "Interface_Target_Residues": interface,
            "Interface_BuriedArea": 812.0, "Interface_Residues": len(contacts),
            "Hotspot_Contact_Fraction": 0.75, "Surface_Hydrophobicity": 0.31,
            "Binder_Net_Charge": 4.0, "Binder_pI": 9.2,
            "bindcraft_version": "1.0.1",
        })
    with (ranked / "!_Ranked.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return folder


def run_campaign_test(numbering, human, emit):
    """Read a synthetic campaign folder end to end and assert what comes back."""
    emit("2. Reading a campaign folder, end to end")
    emit()
    emit("   The classification tests above say nothing about whether this")
    emit("   script can read what the design pipeline actually emits, which is")
    emit("   the part most likely to be wrong. So a synthetic campaign folder is")
    emit("   built in the layout the pipeline documents -- accepted complexes and")
    emit("   binder-only files together in 3_Ranked/, alongside a !_Ranked.csv")
    emit("   with its real column names -- and read back through the same code")
    emit("   path a real folder would take.")
    emit()
    failures = []
    quiet = []
    with tempfile.TemporaryDirectory() as tmpdir:
        folder = build_test_campaign(Path(tmpdir) / "campaign", human)
        # The run's own narration is collected rather than printed: what matters
        # here is the assertions below, not a second copy of a report about a
        # folder that does not exist outside this test.
        summaries, pairs, found = process_campaign(
            folder, numbering, human, False, lambda text="": quiet.append(text))
        failures.extend(found)

        checks = [
            ("the two accepted complexes are found and the binder-only files "
             "are skipped", len(summaries), 2),
            ("contact pairs are read from them",
             len(pairs) >= 6, True),
            ("the metrics table is joined to every structure",
             sum(1 for s in summaries if s.get("metrics")), 2),
            ("the pipeline's own i_pDAE is carried through",
             all(s["metrics"].get("i_pDAE") for s in summaries), True),
            ("our contact set contains the pipeline's",
             all((s.get("interface_cross_check") or (None,))[0] == "contains"
                 for s in summaries), True),
            ("the candidate with a binder histidine facing H370 is rejected",
             next(s["verdict"] for s in summaries if "d4e5f6" in s["design"]),
             VERDICT_REJECTED),
            ("the candidate with four correct pairs meets the target",
             next(s["verdict"] for s in summaries if "a1b2c3" in s["design"]),
             VERDICT_MEETS),
            ("ranking puts the four-pair candidate first",
             rank(summaries)[0]["design"].startswith("design_a1b2c3"), True),
        ]
        for label, got, want in checks:
            ok = got == want
            emit(f"   [{'PASS' if ok else 'FAIL'}] {label}: {got}"
                 + ("" if ok else f", expected {want}"))
            if not ok:
                failures.append(f"campaign folder: {label} was {got}, "
                                f"expected {want}")
    emit()
    emit("   Which chain is which was decided by measurement inside that run,")
    emit("   not by the letter. The pipeline puts the target on chain A and the")
    emit("   binder on B, the reverse of BindCraft version 1, so a parser that")
    emit("   trusted the letter would read the wrong molecule and report a full")
    emit("   set of plausible nonsense.")
    emit()
    return failures


def run_self_test():
    """Remove each branch of the rule in turn and confirm a test case notices.

    The test cases above all pass. That on its own does not establish that they
    test anything: a case that would still pass with the rule deleted is
    decoration. So each branch is disabled in a separate process and the exit
    code is genuinely observed, and the run is only sound if every branch has at
    least one case that fails without it.
    """
    import subprocess

    print("=" * 72)
    print("SELF-TEST: is each branch of the rule actually covered by a case?")
    print("=" * 72)
    print()
    print("Each run below disables one branch of the rule and runs the same test")
    print("cases. Each should fail, and should name the case that caught it. A")
    print("branch that can be removed with every case still passing is a branch")
    print("nothing is testing.")
    print()
    failures = []
    for rule, (expected_case, consequence) in BREAKABLE_RULES.items():
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()),
             "--break-rule", rule],
            capture_output=True, text=True,
            cwd=str(Path(__file__).resolve().parent))
        caught = [line.strip() for line in proc.stdout.splitlines()
                  if line.strip().startswith("- ")]
        named = any(expected_case in line for line in caught)
        ok = proc.returncode != 0 and named
        print(f"   [{'PASS' if ok else 'FAIL'}] {rule}")
        print(f"     if unnoticed: {consequence}")
        print(f"     exit code {proc.returncode}, "
              f"{len(caught)} assertion(s) failed")
        for line in caught[:4]:
            print(f"       {line}")
        if not ok:
            if proc.returncode == 0:
                failures.append(
                    f"{rule}: removing this branch broke no test case")
            else:
                failures.append(
                    f"{rule}: tests failed but not the case meant to cover it "
                    f"({expected_case})")
        print()

    print("=" * 72)
    if failures:
        print(f"SELF-TEST FAILED: {len(failures)} branch(es) are not covered.")
        for failure in failures:
            print(f"  - {failure}")
        print("Add or fix a test case before trusting this filter.")
    else:
        print(f"SELF-TEST PASSED. All {len(BREAKABLE_RULES)} branches of the "
              "rule are covered by a")
        print("test case that fails when the branch is removed.")
    print("=" * 72)
    return 1 if failures else 0


def process_campaign(folder, numbering, human, parse_only, emit,
                     input_target=None):
    """Every accepted design in a campaign folder, classified."""
    failures = []
    structures, where = find_candidate_structures(folder)
    emit(f"   {folder}")
    emit(f"   {len(structures)} candidate complex(es), from {where}")
    if not structures:
        emit("   Nothing to score.")
        return [], [], failures
    emit()

    metrics, columns, metrics_path = read_metrics_table(folder, emit)
    if metrics_path is not None:
        emit(f"   Metrics table: {metrics_path.relative_to(folder)}, "
             f"{len(metrics)} rows, {len(columns)} columns.")
        emit("   Every column is carried through to our output unchanged. None")
        emit("   of them is used to rank anything here.")
        strengthish = [c for c in columns
                       if any(term in c for term in
                              ("i_pTM", "i_pAE", "i_pDAE", "pLDDT", "pTM",
                               "BuriedArea", "Interface_Residues",
                               "Hydrophobicity", "Contact_Fraction",
                               "Net_Charge", "pI"))]
        emit("   Of those, the ones describing interface quality, confidence or")
        emit(f"   charge: {', '.join(strengthish) if strengthish else 'none found'}")
        emit()

    emit("   Which chain is the target and which the binder, by measurement:")
    emit()
    first = load_complex(structures[0])
    target_chain, binder_chains, _evidence = identify_chains(
        first, numbering, human, emit)
    emit()
    if target_chain is None:
        if not parse_only:
            emit("   No chain reads as human EGFR. Stopping.")
            emit()
            emit("   If this is the smoke run against the pipeline's own example")
            emit("   target, that is the expected answer and --parse-only is the")
            emit("   way to run it: the verdicts would be meaningless because the")
            emit("   target is not ours, but the parsing is what is being tested.")
            failures.append("no chain in the candidate structures reads as "
                            "human EGFR; use --parse-only if this is the smoke run")
            return [], [], failures
        chain_ids = [c.id for c in first if common.protein_residues(c)]
        target_chain, binder_chains = chain_ids[0], chain_ids[1:]
        emit("   No chain reads as human EGFR, which is expected for the smoke")
        emit("   run. Falling back to the pipeline's documented convention, that")
        emit("   the first chain is the target and the rest are the binder:")
        emit(f"   target {target_chain}, binder {', '.join(binder_chains)}.")
        emit()
    else:
        emit(f"   Target is chain {target_chain}; binder is "
             f"{', '.join(binder_chains)}.")
        emit("   Taken from the measurement above rather than from the letter.")
        emit("   The pipeline puts the target on A and the binder on B, which is")
        emit("   the reverse of BindCraft version 1, so a parser that assumed a")
        emit("   letter would read the wrong molecule and report a full set of")
        emit("   plausible nonsense.")
        emit()

    if parse_only:
        emit("   Running with --parse-only: contacts are reported, verdicts are")
        emit("   not, and nothing is written to the shortlist. The point is to")
        emit("   prove this script can read what the pipeline emits.")
        emit()

    if not parse_only:
        input_target = input_target or input_target_path()
        first_check = target_numbering.check_complex(input_target, structures[0])
        emit("   Target numbering, first candidate against the structure the run")
        emit(f"   was handed ({Path(input_target).name}):")
        emit(f"   {target_numbering.describe(first_check)}.")
        emit("   Every candidate is checked the same way, and one that cannot be")
        emit("   reconciled is refused rather than scored.")
        emit()
    else:
        input_target = None

    summaries, all_pairs = [], []
    for path in structures:
        result = analyse_candidate(
            path.stem, path, numbering, human,
            target_chain=target_chain, binder_chains=binder_chains,
            translate=not parse_only, input_target=input_target)
        row = join_metrics(path.stem, metrics)
        result["metrics"] = row
        checked = cross_check_interface(result, row, emit) if row else None
        result["interface_cross_check"] = checked
        if checked and checked[0] == "DISAGREES":
            failures.append(f"{path.stem}: {checked[1]}")
        summaries.append(result)
        all_pairs.extend(result["pairs"])

    emit(f"   Scored {len(summaries)} candidates, {len(all_pairs)} contact pairs.")
    emit()

    if any(s.get("interface_cross_check") for s in summaries):
        agreeing = sum(1 for s in summaries
                       if (s.get("interface_cross_check") or (None,))[0]
                       == "contains")
        emit("   Cross-check against the pipeline's own interface list:")
        emit(f"     {agreeing} of {len(summaries)} candidates have our 4.5 A")
        emit(f"     contact set containing the pipeline's {PIPELINE_CUTOFF} A one,")
        emit("     which is what a looser cutoff should give. A residue the")
        emit("     pipeline reports that we do not would be a real disagreement")
        emit("     and stops the run.")
        emit()

    if not parse_only:
        report_results(summaries, emit)
    return summaries, all_pairs, failures


def report_results(summaries, emit):
    """The ranked table, and what it says about the design."""
    ranked = rank(summaries)
    emit("3b. The candidates, ranked by correct pairs")
    emit()
    emit("   Ties are broken first by how many correct pairs rest on E424, near the")
    emit("   trimmed target's cut edge, then by how many rest on E400, E421 or E424,")
    emit("   the three anchors nearest the human-only N361 sugar-chain attachment")
    emit("   point (docs/explainers/08-sugar-chains-near-the-anchors.md) -- both")
    emit("   demotions, never exclusions. Remaining ties are broken by how many of")
    emit("   those pairs also have their charged groups within reach of each other,")
    emit("   then by how many distinct target positions are paired, then by the")
    emit("   pipeline's own i_pDAE with the more confident interface first.")
    emit()
    emit("   | rank | design | correct | N361-side | of those, in reach | "
         "unresolved | neutral | verdict |")
    emit("   |---|---|---|---|---|---|---|---|")
    for index, summary in enumerate(ranked, start=1):
        emit(f"   | {index} | {summary['design']} | {summary['correct_pairs']} | "
             f"{summary['n361_reliant_pairs']} | "
             f"{summary['correct_pairs_tight']} | "
             f"{summary['unresolved_pairs']} | {summary['neutral_pairs']} | "
             f"{summary['verdict']}"
             + (f" ({'; '.join(summary['reasons'])})"
                if summary["reasons"] else "") + " |")
    emit()

    meets = [s for s in ranked if s["verdict"] == VERDICT_MEETS]
    rejected = [s for s in ranked if s["verdict"] == VERDICT_REJECTED]
    below = [s for s in ranked if s["verdict"] == VERDICT_BELOW]
    unscored_list = [s for s in ranked if s["verdict"] == VERDICT_UNSCORABLE]
    at_target = [s for s in meets if s["correct_pairs"] >= PAIR_TARGET]
    emit(f"   {len(meets)} candidates reach {MIN_CORRECT_PAIRS} correct pairs, "
         f"of which {len(at_target)} reach {PAIR_TARGET}.")
    emit(f"   {len(rejected)} rejected for breaking a hard rule.")
    emit(f"   {len(unscored_list)} not scored, because the target numbering "
         f"could not be reconciled with the input.")
    for s in unscored_list:
        emit(f"     {s['design']}: {s['numbering_detail']}")
    emit(f"   {len(below)} below the pair target, kept and ranked last rather")
    emit("   than discarded.")
    emit()
    emit("   What this changes about the design. The shortlist is the candidates")
    emit("   that were not rejected, in this order. Choosing among them is a")
    emit("   judgement to make by hand using the pipeline's own numbers carried")
    emit("   through alongside. The aim is the widest gap between pH 6.5 and pH 7.4")
    emit("   with affinity at pH 6.5 as high as the switch allows.")
    emit()


def pipeline_confidence(summary):
    """A sort key from the pipeline's own `i_pDAE`, so the more confident interface
    sorts first.

    BindCraft2 defines `i_pDAE` as a distance-masked interface TM confidence between
    0 and 1 and ranks it higher-is-better: its own `rank.py` lists the metrics where
    lower is better, and this is not among them. The key is the negated value so
    that an ascending sort puts the higher reading first. The first version of this
    function sorted the other way round, on a glossary entry that said lower was
    better, and so preferred the less confident interface.

    A candidate with no value, or one that is not a number, gets infinity and so
    sorts after every candidate that has one. With `--break-rule pdae-tiebreak`
    every candidate gets the same value, which removes the term.
    """
    if BROKEN_RULE == "pdae-tiebreak":
        return 0.0
    value = (summary.get("metrics") or {}).get("i_pDAE")
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("inf")
    return float("inf") if number != number else -number


def rank(summaries):
    """Order the candidates. The one place the ranking rule lives.

    By supported pairs, meaning correct pairs not resting on a position near the
    trimmed target's cut edge (E424), then by anchor-group-supported pairs,
    meaning correct pairs not resting on the three anchors nearest the
    human-only N361 sequon (E400, E421, E424), then by all correct pairs, then
    by how many of those have their charged groups within reach, then by how
    many distinct target positions are involved. Candidates still level after
    all of those are ordered by the pipeline's own `i_pDAE`, the higher
    reading first, and last by name so the order is stable between runs.

    The two demotion terms are independent and both never-reject: each only
    distinguishes among candidates that are otherwise equal on the term before
    it, and raw `correct_pairs` always gets a say once both are tied. Neither
    can make a design with strictly more correct pairs lose to one with fewer.
    E424 sits in both flagged sets, so a design leaning on it is demoted on
    both grounds rather than having one demotion double-counted.

    The pair terms come first because the charge pairs are what produce the switch.
    `i_pDAE` is deliberately the last term before the name: it breaks ties among
    equals and never overrides a pair term.
    """
    return sorted(summaries,
                  key=lambda s: (-s["supported_pairs"],
                                 -s["anchor_group_supported_pairs"],
                                 -s["correct_pairs"],
                                 -s["correct_pairs_tight"],
                                 -len(s["anchors_paired"]),
                                 pipeline_confidence(s),
                                 s["design"]))


def write_outputs(summaries, all_pairs, emit, parse_only, had_candidates):
    DERIVED.mkdir(parents=True, exist_ok=True)
    CANDIDATES.mkdir(parents=True, exist_ok=True)

    pair_columns = ["design", "target_pos", "target_resnum", "target_aa",
                    "binder_chain", "binder_pos", "binder_aa", "min_dist",
                    "atom_pairs", "classification", "charge_group_dist",
                    "charge_atoms_resolved", "charge_groups_reach", "forbidden"]
    with (DERIVED / "10-candidate-pairs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_columns)
        writer.writeheader()
        for row in all_pairs:
            out = dict(row)
            if out["min_dist"] is not None:
                out["min_dist"] = f"{out['min_dist']:.3f}"
            if out["charge_group_dist"] is not None:
                out["charge_group_dist"] = f"{out['charge_group_dist']:.3f}"
            writer.writerow(out)

    metric_columns = []
    for summary in summaries:
        for column in summary.get("metrics", {}):
            if column not in metric_columns:
                metric_columns.append(column)
    summary_columns = ["design", "verdict", "reasons", "correct_pairs",
                       "supported_pairs", "edge_reliant_pairs",
                       "anchor_group_supported_pairs", "n361_reliant_pairs",
                       "correct_pairs_in_reach", "unresolved_pairs",
                       "his_his_pairs", "forbidden_contacts", "neutral_pairs",
                       "total_pairs", "target_numbering",
                       "target_positions_paired",
                       "target_positions_contacted"]
    with (DERIVED / "10-candidate-summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_columns + metric_columns)
        writer.writeheader()
        for summary in rank(summaries):
            writer.writerow(summary_row(summary))

    # The shortlist holds every candidate that broke no hard rule, in order. A
    # candidate is dropped here for a histidine facing a histidine or a contact
    # at 442, never for the strength of its binding in either direction.
    survivors = [s for s in rank(summaries)
                 if s["verdict"] not in (VERDICT_REJECTED, VERDICT_UNSCORABLE)]
    with (CANDIDATES / "shortlist.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["rank", "meets_pair_target"] + summary_columns
            + metric_columns)
        writer.writeheader()
        for index, summary in enumerate(survivors, start=1):
            row = summary_row(summary)
            row["rank"] = index
            row["meets_pair_target"] = summary["verdict"] == VERDICT_MEETS
            writer.writerow(row)

    emit("4. Files written")
    emit()
    emit(f"   data/derived/10-candidate-pairs.csv      "
         f"{len(all_pairs)} contact pairs")
    emit(f"   data/derived/10-candidate-summary.csv    "
         f"{len(summaries)} candidates")
    emit(f"   results/candidates/shortlist.csv         "
         f"{len(survivors)} not rejected")
    if not had_candidates:
        emit()
        emit("   All three are empty apart from their headers, because no real")
        emit("   candidate exists yet. They are written anyway so the columns")
        emit("   the design run has to fill are visible before it runs rather")
        emit("   than after.")
    if parse_only:
        emit()
        emit("   Written from a --parse-only run, so the verdict columns carry")
        emit("   no meaning. Do not read anything into them.")
    emit()


def summary_row(summary):
    return dict(
        design=summary["design"],
        verdict=summary["verdict"],
        reasons="; ".join(summary["reasons"]),
        correct_pairs=summary["correct_pairs"],
        supported_pairs=summary["supported_pairs"],
        edge_reliant_pairs=summary["edge_reliant_pairs"],
        anchor_group_supported_pairs=summary["anchor_group_supported_pairs"],
        n361_reliant_pairs=summary["n361_reliant_pairs"],
        correct_pairs_in_reach=summary["correct_pairs_tight"],
        unresolved_pairs=summary["unresolved_pairs"],
        his_his_pairs=summary["his_his_pairs"],
        forbidden_contacts=summary["forbidden_contacts"],
        neutral_pairs=summary["neutral_pairs"],
        total_pairs=summary["total_pairs"],
        target_numbering=summary.get("numbering_status", "not checked"),
        target_positions_paired=" ".join(str(p)
                                         for p in summary["anchors_paired"]),
        target_positions_contacted=" ".join(str(p) for p
                                            in summary["target_positions"]),
        **summary.get("metrics", {}))


if __name__ == "__main__":
    sys.exit(main())
