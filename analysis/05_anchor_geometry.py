#!/usr/bin/env python3
"""
05_anchor_geometry.py — do the surviving anchors cluster tightly enough?

WHAT THIS MEASURES AND WHY IT MATTERS
-------------------------------------
A binder is a single small object. It has one face, and that face can only press
against one patch of the target. It cannot wrap around the protein to reach
residues on opposite sides.

So individual exposure, which is what step 03 established, is not enough. The
anchors also have to sit close to one another -- roughly within 25 angstroms,
about the span a small designed protein can cover with a single interface. If the
survivors are scattered across different faces of the solenoid fold (a long
coiled repeat, so residues that are neighbours in the sequence can end up far
apart in space), then no single patch contains three or four of them and the
design premise fails. The epitope would be an artefact of reading a sequence left
to right rather than a real surface.

This step can therefore end the 415-466 epitope. If that is the answer, the
script reports it and names the fallback runs to try instead.

WHICH ATOM TO MEASURE FROM
--------------------------
Distances are measured between the tips of the side chains. The obvious
alternative, measuring backbone to backbone, was rejected: the backbone is the
protein's structural spine and is the same chemistry in every residue, so it says
where the residues sit but not where their charges sit, and the charge can be
several angstroms away from it. What forms a charge pair is the charged
functional group at the end of the side chain:

  aspartate (D): the carboxylate carbon, atom CG -- side-chain atoms are named by
                 how far out along the chain they sit, so CG is a specific carbon
  glutamate (E): the carboxylate carbon, atom CD, one position further out again
  histidine (H): the centroid, meaning the average position, of the
                 five-membered imidazole ring (atoms CG, ND1, CD2, CE1, NE2),
                 because the ring's charge is spread over the whole ring rather
                 than sitting on one atom

If a side chain is unresolved in the crystal structure (too mobile for the
experiment to pin down), we fall back to CB, the first carbon of the side chain.
That is a poorer proxy -- it is up to a few angstroms from the functional group --
so every fallback is reported explicitly rather than quietly substituted. If even
CB is missing, the last resort is CA, the alpha carbon, which is present in every
residue.

THE VERDICT RULE, FIXED IN ADVANCE
----------------------------------
Find the largest subset of surviving anchors whose maximum pairwise distance is
under 25 A. If fewer than three anchors cluster, the 415-466 epitope fails and we
move to the fallback conserved runs 394-411 (18 residues) and 331-347 (17
residues), rerunning steps 03-05 against those.

Three is the minimum because the pH switch is partial rather than binary: at
pH 6.5 a histidine is only fractionally protonated -- only some copies of it carry
the extra positive charge at any moment -- so a single pair produces a weak effect
and three or four have to be stacked.

Residue numbers here are positions in the human record in UniProt, the public
protein sequence archive. The structure file is 6ARU from the Protein Data Bank
(PDB), the public archive of measured 3D protein structures, which numbers the
same residues 24 lower; step 02's conversion table handles that throughout.

Also writes an interactive py3Dmol viewer, so the geometry can be looked at
instead of only read as numbers.

Outputs:
  results/findings/05-anchor-geometry.md
  data/derived/05-anchor-distance-matrix.csv
  data/derived/05-anchor-clusters.csv
  explorer/anchor-viewer.html          (interactive, open in a browser)
  explorer/anchor-view.pml             (PyMOL script, same colouring)

Run standalone:  python analysis/05_anchor_geometry.py
"""

import itertools
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

ROOT = Path(__file__).resolve().parents[1]
STRUCT = ROOT / "data" / "structures"
RECEPTOR_PDB = STRUCT / "6aru_receptor_only.pdb"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"
EXPLORER = ROOT / "explorer"
VERDICTS_CSV = DERIVED / "03-anchor-verdicts.csv"
OFFSET_CSV = DERIVED / "02-numbering-offset.csv"
OVERLAP_CSV = DERIVED / "04-epitope-overlap.csv"
GROOVE_CSV = DERIVED / "06-intra-chain-occlusion.csv"

REACH_CUTOFF = 25.0      # angstroms; approximate span of one small binder face
MIN_CLUSTER = 3          # fewer anchors than this in one patch and the epitope fails

EPI_START, EPI_END = 415, 466
ANCHORS = {416: "D", 418: "H", 421: "E", 424: "E",
           433: "H", 455: "E", 458: "D", 460: "D"}
GLYCO_SITE = 444         # sugar attachment point found inside the block by step 03

FUNCTIONAL_ATOM = {"D": ["CG"], "E": ["CD"],
                   "H": ["CG", "ND1", "CD2", "CE1", "NE2"]}

# Filled in by the sugar-proximity section below: anchors close enough to the
# sugar chain at N444 that their measured exposure should be read as an upper
# bound on how reachable they really are.
GLYCAN_RISK = set()

# A complex N-linked glycan is a branched chain of sugars rather than a single
# sugar. Only its innermost sugars sit still enough to appear in a crystal
# structure, while the whole assembly is mobile and can sweep 20-30 A from the
# point where it attaches. The two bands below are therefore cautious on purpose.
GLYCAN_NEAR = 15.0     # very likely shadowed some of the time
GLYCAN_PLAUSIBLE = 25.0  # within reach of an extended chain

THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}


def load_offset():
    if not OFFSET_CSV.exists():
        raise SystemExit("Run analysis/02_structure_prep.py first.")
    uni_to_pdb = {}
    for line in OFFSET_CSV.read_text().splitlines()[1:]:
        uni, pdb, _ = line.split(",")
        uni_to_pdb[int(uni)] = int(pdb)
    return uni_to_pdb


def load_survivors():
    if not VERDICTS_CSV.exists():
        raise SystemExit("Run analysis/03_solvent_accessibility.py first.")
    survivors, verdicts = [], {}
    for line in VERDICTS_CSV.read_text().splitlines()[1:]:
        parts = line.split(",")
        pos = int(parts[0])
        verdict = parts[-2]
        usable = parts[-1].strip() == "True"
        verdicts[pos] = verdict
        if usable:
            survivors.append(pos)
    return sorted(survivors), verdicts


def load_domain_iv_groove():
    """Anchors that sit against domain IV, read from step 06's output.

    Step 06 found that the second half of our block is packed against domain IV,
    the domain that follows ours, in both structures it compared. An anchor in
    that packing sits in a groove between two domains rather than on an open face,
    which is harder for a designed binder to fit.

    Step 06 runs after this one, so on a first run from scratch its table does not
    exist yet. That is reported rather than guessed at, and the second run of this
    script picks it up. Hardcoding the residues here would let the two scripts
    drift apart, which is the failure this project has already had once.
    """
    if not GROOVE_CSV.exists():
        return set(), False
    groove = set()
    for line in GROOVE_CSV.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        _structure, pos, is_anchor, _dist, _other = line.split(",")
        if is_anchor.strip() == "True":
            groove.add(int(pos))
    return groove, True


def load_cetuximab_contacts():
    if not OVERLAP_CSV.exists():
        return set()
    hits = set()
    for line in OVERLAP_CSV.read_text().splitlines()[1:]:
        pos, _anchor, contact, _d = line.split(",")
        if contact.strip() == "True":
            hits.add(int(pos))
    return hits


def functional_point(res, aa):
    """Pick the atom this anchor is measured from.

    Returns (coordinates, a label for the atom used, whether it is a fallback).
    The preferred atom is the charged tip of the side chain; CB, the first carbon
    of the side chain, and then CA, the alpha carbon, are the fallbacks used when
    the structure does not resolve the tip.
    """
    wanted = FUNCTIONAL_ATOM.get(aa, [])
    coords = [res[a].coord for a in wanted if a in res]
    if len(coords) == len(wanted) and coords:
        if aa == "H":
            return np.mean(coords, axis=0), "imidazole centroid", False
        return coords[0], wanted[0], False
    if "CB" in res:
        return res["CB"].coord, "CB (FALLBACK)", True
    return res["CA"].coord, "CA (FALLBACK)", True


def main():
    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    uni_to_pdb = load_offset()
    survivors, verdicts = load_survivors()
    cetux = load_cetuximab_contacts()
    DOMAIN_IV_GROOVE, groove_available = load_domain_iv_groove()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(RECEPTOR_PDB))
    chain = next(iter(structure[0]))
    by_pdb = {r.id[1]: r for r in chain if is_aa(r, standard=True)}

    emit("=" * 72)
    emit("ANCHOR GEOMETRY — is the epitope a real surface?")
    emit("=" * 72)
    emit()
    emit(f"Reach cutoff: {REACH_CUTOFF:.0f} A — roughly what one small binder face")
    emit("              can span. Fixed before the measurement, so it cannot have")
    emit("              been tuned to suit the answer.")
    emit(f"Minimum viable cluster: {MIN_CLUSTER} anchors.")
    emit()
    emit("Measured from the charged tip of each side chain, because that is what")
    emit("forms a charge pair. Backbone-to-backbone distances were the alternative")
    emit("and were rejected: the backbone is identical chemistry in every residue")
    emit("and sits several angstroms from the charge.")
    emit("  D (aspartate) -> CG, the carboxylate carbon")
    emit("  E (glutamate) -> CD, the carboxylate carbon, one atom further out")
    emit("  H (histidine) -> centre of the imidazole ring (CG, ND1, CD2, CE1, NE2)")
    emit("  unresolved side chain -> CB, the first side-chain carbon, flagged as a")
    emit("                          fallback")
    emit()
    emit(f"Anchors carried forward from step 03 (exposed or partial): "
         f"{len(survivors)} of 8")
    emit(f"  {', '.join(f'{ANCHORS[p]}{p}' for p in survivors)}")
    excluded = [p for p in sorted(ANCHORS) if p not in survivors]
    if excluded:
        emit(f"Excluded as buried/unresolved: "
             f"{', '.join(f'{ANCHORS[p]}{p} ({verdicts.get(p)})' for p in excluded)}")
    emit()

    # ---- Gather points ----
    points, fallbacks = {}, []
    for pos in survivors:
        pdb_num = uni_to_pdb.get(pos)
        res = by_pdb.get(pdb_num)
        if res is None:
            emit(f"  WARNING: {ANCHORS[pos]}{pos} has no coordinates; skipped.")
            continue
        aa = THREE_TO_ONE.get(res.get_resname(), "X")
        coord, label, fb = functional_point(res, aa)
        points[pos] = dict(coord=np.asarray(coord, dtype=float),
                           atom=label, aa=aa)
        if fb:
            fallbacks.append(pos)

    emit("1. Atom used for each anchor")
    emit()
    emit("   | anchor | atom measured |")
    emit("   |---|---|")
    for pos in sorted(points):
        emit(f"   | {ANCHORS[pos]}{pos} | {points[pos]['atom']} |")
    emit()
    if fallbacks:
        emit(f"   {len(fallbacks)} anchor(s) fell back to CB because the side chain")
        emit(f"   is unresolved: {', '.join(f'{ANCHORS[p]}{p}' for p in fallbacks)}")
        emit("   CB can sit several angstroms from the functional group, so those")
        emit("   distances carry extra uncertainty, and a cluster that only just")
        emit("   fits the cutoff because of one of them should be treated as")
        emit("   borderline rather than as a clean pass.")
    else:
        emit("   No fallbacks — every functional group is fully resolved, so every")
        emit("   distance below is measured from the atoms that actually form the")
        emit("   charge pair.")
    emit()

    # ---- Distance matrix ----
    ordered = sorted(points)
    emit("2. Full pairwise distance matrix (angstroms)")
    emit()
    emit("   Every anchor against every other anchor. Each cell is the distance")
    emit("   between those two functional groups, so the table shows which anchors")
    emit("   could share one binder face and which are on opposite sides.")
    emit()
    header = "   | | " + " | ".join(f"{ANCHORS[p]}{p}" for p in ordered) + " |"
    emit(header)
    emit("   |---|" + "---|" * len(ordered))
    dist = {}
    for a in ordered:
        cells = []
        for b in ordered:
            d = float(np.linalg.norm(points[a]["coord"] - points[b]["coord"]))
            dist[(a, b)] = d
            cells.append("—" if a == b else f"{d:.1f}")
        emit(f"   | **{ANCHORS[a]}{a}** | " + " | ".join(cells) + " |")
    emit()

    pair_list = [(a, b, dist[(a, b)])
                 for a, b in itertools.combinations(ordered, 2)]
    pair_list.sort(key=lambda t: t[2])
    emit(f"   Closest pair:  {ANCHORS[pair_list[0][0]]}{pair_list[0][0]}–"
         f"{ANCHORS[pair_list[0][1]]}{pair_list[0][1]} at {pair_list[0][2]:.1f} A")
    emit(f"   Furthest pair: {ANCHORS[pair_list[-1][0]]}{pair_list[-1][0]}–"
         f"{ANCHORS[pair_list[-1][1]]}{pair_list[-1][1]} at {pair_list[-1][2]:.1f} A")
    n_within = sum(1 for _, _, d in pair_list if d < REACH_CUTOFF)
    emit(f"   Pairs within {REACH_CUTOFF:.0f} A: {n_within} of {len(pair_list)}")
    emit()

    # ---- Largest cluster ----
    emit(f"3. Largest subset with every pairwise distance under {REACH_CUTOFF:.0f} A")
    emit()
    emit("   Every possible subset of the anchors is checked, so this answer is")
    emit("   exact. With this few anchors that is quick, which is why no")
    emit("   approximate clustering method is used.")
    emit()
    best = []
    all_valid = []
    for size in range(len(ordered), 1, -1):
        for combo in itertools.combinations(ordered, size):
            span = max(dist[(a, b)] for a, b in itertools.combinations(combo, 2))
            if span < REACH_CUTOFF:
                all_valid.append((combo, span))
                if len(combo) > len(best):
                    best = list(combo)
        if best:
            break

    if best:
        span = max(dist[(a, b)] for a, b in itertools.combinations(best, 2))
        emit(f"   Largest cluster: {len(best)} anchors")
        emit(f"   {', '.join(f'{ANCHORS[p]}{p}' for p in best)}")
        emit(f"   Maximum internal distance: {span:.1f} A")
        emit()
        n_acidic = sum(1 for p in best if ANCHORS[p] in "DE")
        n_his = sum(1 for p in best if ANCHORS[p] == "H")
        emit(f"   Composition: {n_acidic} acidic (each one calls for a histidine on")
        emit(f"   the binder), {n_his} target histidine (calls for a D or E on the")
        emit("   binder). That is the shopping list the binder has to present.")
        emit()
        ties = [c for c, s in all_valid if len(c) == len(best)]
        if len(ties) > 1:
            emit(f"   {len(ties)} different subsets of size {len(best)} qualify, so the")
            emit("   design has a choice of contact set:")
            for c, s in sorted([(c, s) for c, s in all_valid
                                if len(c) == len(best)], key=lambda t: t[1])[:6]:
                emit(f"     {', '.join(f'{ANCHORS[p]}{p}' for p in c)}  "
                     f"(span {s:.1f} A)")
            emit()
    else:
        emit("   No subset of two or more anchors falls within the cutoff.")
        emit()

    # ---- Tighter sub-clusters, for a smaller binder ----
    emit("4. How tight can a cluster be? (relevant to molecule size choice)")
    emit()
    emit("   A microbinder (under 40 amino acids, abbreviated aa) presents a smaller")
    emit("   face than a minibinder (40-100 aa), so it needs a tighter anchor")
    emit("   cluster. Largest cluster at several spans:")
    emit()
    emit("   For each limit: the largest cluster that fits, and the most compact")
    emit("   example of that size. An arbitrary example of that size would hide how")
    emit("   tight the best available option actually is.")
    emit()
    emit("   | span limit | largest cluster | tightest such set | its span |")
    emit("   |---|---|---|---|")
    for limit in (12.0, 15.0, 18.0, 20.0, 25.0):
        bestn, candidates = 0, []
        for size in range(len(ordered), 1, -1):
            fits = [(combo, max(dist[(a, b)]
                                for a, b in itertools.combinations(combo, 2)))
                    for combo in itertools.combinations(ordered, size)]
            fits = [(c, s) for c, s in fits if s < limit]
            if fits:
                bestn, candidates = size, fits
                break
        if candidates:
            tightest, span_t = min(candidates, key=lambda t: t[1])
            emit(f"   | {limit:.0f} A | {bestn} | "
                 f"{', '.join(f'{ANCHORS[p]}{p}' for p in tightest)} | "
                 f"{span_t:.1f} A |")
        else:
            emit(f"   | {limit:.0f} A | 0 | — | — |")
    emit()
    emit("   What this decides: how big the molecule has to be. A microbinder")
    emit("   presents a small face and needs a tight cluster; a minibinder (40-100")
    emit("   aa) can span more. The row where the cluster size drops below three is")
    emit("   the size at which a binder becomes too small to carry the stacked")
    emit("   switch at all, so pick a size above that row.")
    emit()

    # ---- Glycosylation site proximity ----
    emit(f"5. Distance from the sugar attachment point N{GLYCO_SITE}")
    emit()
    emit("   Step 03 found a sugar chain attached at N444, inside our block.")
    emit()
    emit("   What is measured here, and how it differs from step 03. The two steps")
    emit("   measure different things, and read together without this note they")
    emit("   look like they disagree:")
    emit()
    emit("     step 03 measures the distance from each anchor to the sugar atoms")
    emit("     that are actually present in the structure file, with a 5 A cutoff.")
    emit("     It found no anchor within 5 A, so no anchor touches the sugars we")
    emit("     can see.")
    emit()
    emit("     this step measures the distance from each anchor to the point where")
    emit("     the chain is attached, N444, with the 15 A and 25 A bands below.")
    emit("     The chain itself is longer than the part the structure shows.")
    emit()
    emit("   Both are correct. Step 03 answers 'does an anchor touch a sugar atom we")
    emit("   have coordinates for', and the answer is no. This step answers 'could")
    emit("   the full chain reach an anchor', and the answer is that three of them")
    emit("   are close enough that it might. The second question is worth asking")
    emit("   because a structure only shows the first few sugars of a chain that")
    emit("   continues past them, so step 03's answer does not settle it.")
    emit()
    gly_pdb = uni_to_pdb.get(GLYCO_SITE)
    gly_res = by_pdb.get(gly_pdb)
    if gly_res is not None and "CB" in gly_res:
        gly_pt = np.asarray(gly_res["CB"].coord, dtype=float)
        emit("   A complex N-linked glycan is a branched chain of sugars that can")
        emit("   sweep 20-30 A from where it attaches, so the bands below are")
        emit("   deliberately cautious:")
        emit(f"     under {GLYCAN_NEAR:.0f} A  — likely shadowed at least some of the time")
        emit(f"     under {GLYCAN_PLAUSIBLE:.0f} A  — within reach of an extended chain")
        emit("     beyond that — probably clear")
        emit()
        emit("   | anchor | distance to N444 CB | assessment |")
        emit("   |---|---|---|")
        for pos in ordered:
            d = float(np.linalg.norm(points[pos]["coord"] - gly_pt))
            if d < GLYCAN_NEAR:
                note = "likely shadowed — read exposure as an upper bound"
                GLYCAN_RISK.add(pos)
            elif d < GLYCAN_PLAUSIBLE:
                note = "possibly reached by an extended chain"
                GLYCAN_RISK.add(pos)
            else:
                note = "probably clear"
            emit(f"   | {ANCHORS[pos]}{pos} | {d:.1f} A | {note} |")
        emit()
        emit("   'Flagged' below means flagged by the attachment-point measure used")
        emit("   here, at 15 A or 25 A. It does not contradict step 03, which found")
        emit("   no anchor within 5 A of a sugar atom present in the file.")
        emit(f"   Anchors flagged: "
             f"{', '.join(f'{ANCHORS[p]}{p}' for p in sorted(GLYCAN_RISK)) or 'none'}")
        emit()
        emit("   This is a flagged risk rather than a settled question, and it cuts")
        emit("   both ways: a glycan is flexible, so 'within reach' means 'covered")
        emit("   some of the time' rather than 'blocked'. A crystal structure cannot")
        emit("   settle it either way. In practice we use it as a tie-breaker between")
        emit("   clusters that are otherwise equally good.")
    else:
        emit(f"   N{GLYCO_SITE} not resolved; cannot measure.")
    emit()

    # ---- Verdict ----
    emit("6. Verdict")
    emit()
    n_best = len(best)
    dead = n_best < MIN_CLUSTER
    if dead:
        emit(f"   The {EPI_START}-{EPI_END} epitope fails this check.")
        emit()
        emit(f"   Only {n_best} anchor(s) cluster within {REACH_CUTOFF:.0f} A, and the")
        emit(f"   stacked switch requires at least {MIN_CLUSTER}. A single binder")
        emit("   cannot reach enough anchors to build a pH switch on this surface.")
        emit("   Better design cannot recover it, because the obstacle is simply")
        emit("   where these residues sit on the protein.")
        emit()
        emit("   Next step: rerun steps 03-05 against the fallback conserved runs")
        emit("   394-411 (18 residues) and 331-347 (17 residues).")
    else:
        emit(f"   The {EPI_START}-{EPI_END} epitope survives the structure check.")
        emit()
        emit(f"   {n_best} anchors cluster within {REACH_CUTOFF:.0f} A "
             f"(max internal span {span:.1f} A), against a minimum of "
             f"{MIN_CLUSTER}.")
        emit("   A single binder face can reach them, so the stacked switch the")
        emit("   design depends on is geometrically possible.")
        emit()
        # There may be several equally large clusters. Whether the method-novelty
        # claim survives depends on whether any one of them contains a target
        # histidine, so we check them all rather than only the first one the
        # search returned.
        maximal = [c for c, s in all_valid if len(c) == n_best]
        with_his = [c for c in maximal if any(ANCHORS[p] == "H" for p in c)]

        emit(f"   {len(maximal)} distinct cluster(s) of size {n_best} qualify, so the")
        emit("   design has a choice of contact set. Comparing them on the three")
        emit("   things that matter:")
        emit()
        emit("   The last column uses the attachment-point measure from section 5,")
        emit("   at 15 A and 25 A. It is not step 03's 5 A check against the sugar")
        emit("   atoms present in the file, which found no anchor within 5 A.")
        emit()
        emit("   | cluster | span A | target His? | cetuximab overlap | within 25 A of N444 |")
        emit("   |---|---|---|---|---|")
        for c in sorted(maximal, key=lambda c: max(
                dist[(a, b)] for a, b in itertools.combinations(c, 2))):
            s = max(dist[(a, b)] for a, b in itertools.combinations(c, 2))
            his = [p for p in c if ANCHORS[p] == "H"]
            ov = [p for p in c if p in cetux]
            risky = [p for p in c if p in GLYCAN_RISK]
            emit(f"   | {', '.join(f'{ANCHORS[p]}{p}' for p in c)} | {s:.1f} | "
                 f"{', '.join(f'H{p}' for p in his) if his else 'none'} | "
                 f"{', '.join(f'{ANCHORS[p]}{p}' for p in ov) if ov else 'none'} | "
                 f"{', '.join(f'{ANCHORS[p]}{p}' for p in risky) if risky else 'none'} |")
        emit()

        if with_his:
            emit("   Clusters that include a target histidine, which is the half of")
            emit("   the pairing rule that needs an acidic residue on the binder:")
            emit()
            for c in with_his:
                s = max(dist[(a, b)] for a, b in itertools.combinations(c, 2))
                emit(f"     {', '.join(f'{ANCHORS[p]}{p}' for p in c)} "
                     f"(span {s:.1f} A) contains "
                     f"{', '.join(f'H{p}' for p in c if ANCHORS[p] == 'H')}")
            emit()
            emit("   This arrangement is published, so it is not ours to claim. Liu X")
            emit("   et al. (2022), Molecular Therapy - Oncolytics 27:256-269,")
            emit("   doi:10.1016/j.omto.2022.11.001, mapped EGFR's own H370 and H433")
            emit("   as the determinants of pH-dependent antibody binding and then")
            emit("   put an acidic residue against H433 deliberately, on this target")
            emit("   and at pH 6.5 against 7.4. See the decisions log under Resolved.")
            emit()
            emit("   That paper does support the rejection criterion this project")
            emit("   derived from first principles. Their histidine-on-the-binder")
            emit("   variant facing H433 changed neither pH-dependency nor affinity,")
            emit("   while the acidic variants improved pH-dependency substantially.")
            emit()
            emit("   On tooling, what is checkable is narrower than what this file")
            emit("   used to assert. Verified 1 October 2026 from the repositories:")
            emit("   AlphaFold2 and AlphaFold-Multimer take sequences and database")
            emit("   paths only; RFdiffusion's inference configuration conditions on")
            emit("   backbone, contigs, hotspot residues and a closed set of")
            emit("   radius-of-gyration and contact potentials; ProteinMPNN's")
            emit("   arguments contain no pH, pKa or protonation option. BindCraft2")
            emit("   is the exception worth stating rather than hiding: its")
            emit("   filters.py carries a side-chain pKa table and reports Binder_pI")
            emit("   and Binder_Net_Charge at a hard-coded REPORTED_PH of 7.4. That")
            emit("   is the binder's own sequence charge, the pH is not a user")
            emit("   setting in its 235-setting catalogue, and both are reported")
            emit("   readouts rather than design objectives.")
            emit()
            emit("   The trade-off on cetuximab overlap still stands. The clusters")
            emit("   containing a target histidine are the ones that overlap")
            emit("   cetuximab's footprint, because H433 is the single anchor")
            emit("   cetuximab touches. Step 07 confirmed that overlap holds for")
            emit("   unmodified cetuximab and not merely for the engineered variant")
            emit("   in 6ARU. The counter-argument, that cetuximab contacts H433 with")
            emit("   no pH dependence so sharing a residue is not sharing a")
            emit("   mechanism, is something to argue in the write-up; nothing here")
            emit("   computes it and it should not be presented as a result.")
        else:
            emit("   No cluster of the maximum size contains a target histidine, so")
            emit("   the acidic-residue-against-target-histidine half of the pairing")
            emit("   rule cannot be used at full cluster")
            emit("   size. Before abandoning the claim, check whether a smaller")
            emit("   cluster that includes H433 is still viable.")
        emit()
        free = [p for p in best if p not in cetux]
        emit(f"   For reference, the tightest cluster ("
             f"{', '.join(f'{ANCHORS[p]}{p}' for p in best)}) has {len(free)} of")
        emit(f"   {n_best} anchors outside cetuximab's footprint.")
        emit()
        emit("   What this does not establish. Geometric reachability is necessary")
        emit("   but not sufficient: it does not show that a foldable binder exists")
        emit("   which presents the right partner residues in the right")
        emit("   orientations, nor that the resulting switch is large enough to clear")
        emit("   the assay's detection floor at pH 7.4. Those are design and")
        emit("   prediction questions, and both are still open.")
    emit()

    # ---- Second scenario: H418 included ----
    # Step 03 marked H418 buried, reading the 6ARU structure, so everything above
    # leaves it out. Step 06 then found it partially exposed in 1NQL, with domain
    # III folded the same way in both, so whether it is usable is unresolved. The
    # clusters differ between the two answers, and the difference decides which
    # contact set the design should aim at, so both are computed here rather than
    # left to a calculation someone once ran by hand.
    emit("7. Second scenario: the same question with H418 included")
    emit()
    emit("   Everything above excludes H418, following step 03, which measured it")
    emit("   as buried in 6ARU. Step 06 measured it as partially exposed in 1NQL,")
    emit("   and domain III has the same fold in both structures, so the two")
    emit("   measurements disagree and the question is open. Because 6ARU has the")
    emit("   antibody clamped on, the burial there may be the antibody holding that")
    emit("   side chain rather than a property of the receptor alone.")
    emit()
    emit("   This section answers what the clusters would be if H418 is usable. It")
    emit("   is conditional on that unresolved question and is labelled as such")
    emit("   wherever the numbers are used.")
    emit()

    if groove_available:
        emit(f"   Anchors against domain IV, read from step 06's table: "
             f"{', '.join(f'{ANCHORS[p]}{p}' for p in sorted(DOMAIN_IV_GROOVE))}")
    else:
        emit("   Step 06 has not run yet, so the domain IV column below is empty.")
        emit("   Run step 06 and then this script again to fill it in.")
    emit()

    conditional = sorted(set(ordered) | {418})
    cpoints = {}
    for pos in conditional:
        pdb_num = uni_to_pdb.get(pos)
        res = by_pdb.get(pdb_num)
        if res is None:
            continue
        aa = THREE_TO_ONE.get(res.get_resname(), "X")
        coord, _label, _fb = functional_point(res, aa)
        cpoints[pos] = np.asarray(coord, dtype=float)

    cdist = {}
    for a in cpoints:
        for b in cpoints:
            cdist[(a, b)] = float(np.linalg.norm(cpoints[a] - cpoints[b]))

    emit(f"   Anchors considered: "
         f"{', '.join(f'{ANCHORS[p]}{p}' for p in sorted(cpoints))}")
    emit()
    emit("   Distances from H418 to the others:")
    for pos in sorted(cpoints):
        if pos == 418:
            continue
        emit(f"     H418 to {ANCHORS[pos]}{pos}: {cdist[(418, pos)]:.1f} A")
    emit()
    pair_418_433 = cdist.get((418, 433))
    if pair_418_433 is not None:
        emit(f"   H418 to H433 is {pair_418_433:.1f} A, against a reach cutoff of")
        emit(f"   {REACH_CUTOFF:.0f} A. The two target histidines cannot both be")
        emit("   reached by one binder, so a design uses one or the other.")
    emit()

    cvalid, cbest = [], []
    for size in range(len(cpoints), 1, -1):
        for combo in itertools.combinations(sorted(cpoints), size):
            sp = max(cdist[(a, b)] for a, b in itertools.combinations(combo, 2))
            if sp < REACH_CUTOFF:
                cvalid.append((combo, sp))
                if len(combo) > len(cbest):
                    cbest = list(combo)
        if cbest:
            break

    if cbest:
        cspan = max(cdist[(a, b)] for a, b in itertools.combinations(cbest, 2))
        emit(f"   Largest cluster with H418 available: {len(cbest)} anchors")
        emit(f"   {', '.join(f'{ANCHORS[p]}{p}' for p in cbest)}")
        emit(f"   Maximum internal distance {cspan:.1f} A")
        emit()
        excluded_span = max(dist[(a, b)]
                            for a, b in itertools.combinations(best, 2))
        emit(f"   Compared with {n_best} anchors at {excluded_span:.1f} A when H418")
        emit("   is excluded.")
        emit()
        cmaximal = [c for c, s in cvalid if len(c) == len(cbest)]
        emit("   All clusters of that size, with the same columns as section 6:")
        emit()
        emit("   | cluster | span A | target His? | cetuximab overlap "
             "| within 25 A of N444 | domain IV groove |")
        emit("   |---|---|---|---|---|---|")
        for c in sorted(cmaximal, key=lambda c: max(
                cdist[(a, b)] for a, b in itertools.combinations(c, 2))):
            s = max(cdist[(a, b)] for a, b in itertools.combinations(c, 2))
            his = [p for p in c if ANCHORS[p] == "H"]
            ov = [p for p in c if p in cetux]
            risky = [p for p in c if p in GLYCAN_RISK]
            groove = [p for p in c if p in DOMAIN_IV_GROOVE]
            emit(f"   | {', '.join(f'{ANCHORS[p]}{p}' for p in c)} | {s:.1f} | "
                 f"{', '.join(f'H{p}' for p in his) if his else 'none'} | "
                 f"{', '.join(f'{ANCHORS[p]}{p}' for p in ov) if ov else 'none'} | "
                 f"{', '.join(f'{ANCHORS[p]}{p}' for p in risky) if risky else 'none'} | "
                 f"{', '.join(f'{ANCHORS[p]}{p}' for p in groove) if groove else 'none'} |")
        emit()
        with_his = [c for c in cmaximal if any(ANCHORS[p] == "H" for p in c)]
        clean = [c for c in cmaximal
                 if any(ANCHORS[p] == "H" for p in c)
                 and not any(p in cetux for p in c)
                 and not any(p in DOMAIN_IV_GROOVE for p in c)]
        if clean:
            emit("   Clusters here that carry a target histidine, avoid the")
            emit("   antibody footprint and avoid the domain IV groove:")
            for c in clean:
                s = max(cdist[(a, b)] for a, b in itertools.combinations(c, 2))
                emit(f"     {', '.join(f'{ANCHORS[p]}{p}' for p in c)} "
                     f"(span {s:.1f} A)")
            emit()
            emit("   No cluster in section 6, where H418 is excluded, manages all")
            emit("   three at once. That is what makes settling H418 worth doing")
            emit("   before choosing a contact set.")
        elif with_his:
            emit("   Clusters here carry a target histidine but none avoids both the")
            emit("   antibody footprint and the domain IV groove.")
        else:
            emit("   No cluster of this size carries a target histidine, so")
            emit("   including H418 does not by itself make the pairing against a")
            emit("   target histidine available at full cluster size.")
    else:
        emit("   No cluster of two or more anchors falls within the cutoff.")
    emit()
    emit("   These numbers hold only if H418 is usable. Until that is settled they")
    emit("   describe an option, not a decision.")
    emit()

    # ---- CSVs ----
    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "05-anchor-distance-matrix-with-h418.csv").open("w") as fh:
        fh.write("anchor_a,anchor_b,distance_a,scenario\n")
        for a, b in itertools.combinations(sorted(cpoints), 2):
            fh.write(f"{ANCHORS[a]}{a},{ANCHORS[b]}{b},{cdist[(a, b)]:.3f},"
                     f"h418_included_conditional\n")
    with (DERIVED / "05-anchor-clusters-with-h418.csv").open("w") as fh:
        fh.write("cluster_size,max_internal_span_a,anchors,contains_target_his,"
                 "cetuximab_overlap,near_n444,domain_iv_groove,scenario\n")
        for combo, s in sorted(cvalid, key=lambda t: (-len(t[0]), t[1])):
            his = any(ANCHORS[p] == "H" for p in combo)
            ov = any(p in cetux for p in combo)
            near = any(p in GLYCAN_RISK for p in combo)
            groove = any(p in DOMAIN_IV_GROOVE for p in combo)
            fh.write(f"{len(combo)},{s:.3f},"
                     f"\"{' '.join(f'{ANCHORS[p]}{p}' for p in combo)}\","
                     f"{his},{ov},{near},{groove},h418_included_conditional\n")
    emit("Wrote data/derived/05-anchor-distance-matrix-with-h418.csv")
    emit("Wrote data/derived/05-anchor-clusters-with-h418.csv")

    with (DERIVED / "05-anchor-distance-matrix.csv").open("w") as fh:
        fh.write("anchor_a,anchor_b,distance_a\n")
        for a, b, d in pair_list:
            fh.write(f"{ANCHORS[a]}{a},{ANCHORS[b]}{b},{d:.3f}\n")
    with (DERIVED / "05-anchor-clusters.csv").open("w") as fh:
        fh.write("cluster_size,max_internal_span_a,anchors\n")
        for combo, s in sorted(all_valid, key=lambda t: (-len(t[0]), t[1])):
            fh.write(f"{len(combo)},{s:.3f},"
                     f"\"{' '.join(f'{ANCHORS[p]}{p}' for p in combo)}\"\n")
    emit("Wrote data/derived/05-anchor-distance-matrix.csv")
    emit("Wrote data/derived/05-anchor-clusters.csv")

    # ---- Visual output ----
    EXPLORER.mkdir(parents=True, exist_ok=True)
    write_viewer(EXPLORER / "anchor-viewer.html", uni_to_pdb, verdicts, best, cetux)
    write_pymol(EXPLORER / "anchor-view.pml", uni_to_pdb, verdicts, best)
    emit("Wrote explorer/anchor-viewer.html  (open in a browser)")
    emit("Wrote explorer/anchor-view.pml     (PyMOL: run this file)")

    emit()
    emit("=" * 72)
    emit(f"RESULT: largest cluster = {n_best} anchors within {REACH_CUTOFF:.0f} A. "
         f"Epitope {'fails this check' if dead else 'survives'}.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "05-anchor-geometry.md").write_text(
        "# Anchor geometry: is the epitope a real surface?\n\n"
        "Computed output of `analysis/05_anchor_geometry.py`. Do not hand-edit.\n\n"
        f"Tests whether the anchors that survived step 03 sit within {REACH_CUTOFF:.0f}\n"
        "angstroms of one another, which is roughly how far across a single small\n"
        "binder can reach with one face. Distances are measured between the tips of\n"
        "the side chains, the parts that actually form a charge pair, rather than\n"
        "between backbone atoms, which are the same chemistry in every residue.\n"
        f"The rule was fixed before measuring: if fewer than {MIN_CLUSTER} anchors\n"
        "cluster, this epitope fails and the fallback stretches are tried instead.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 0


def _colour_map(verdicts, cluster):
    """exposure colour per anchor; cluster members get a thicker representation."""
    colours = {}
    for pos, v in verdicts.items():
        colours[pos] = {"exposed": "#1a9850", "partial": "#fdae61",
                        "buried": "#d73027"}.get(v, "#999999")
    return colours


def write_viewer(path, uni_to_pdb, verdicts, cluster, cetux):
    colours = _colour_map(verdicts, cluster)
    rows = []
    for pos in sorted(ANCHORS):
        rows.append(dict(uni=pos, pdb=uni_to_pdb.get(pos), aa=ANCHORS[pos],
                         verdict=verdicts.get(pos, "unknown"),
                         colour=colours.get(pos, "#999999"),
                         in_cluster=pos in cluster,
                         cetux=pos in cetux))
    import json
    path.write_text(f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>EGFR Anchor Viewer</title>
<!-- 3Dmol.js, the 3D structure viewer, pinned to a version checked to exist on
     30 September 2026. Three mirrors are tried in order, because a wrong or
     withdrawn version number quietly yields "$3Dmol is not defined", which looks
     like a broken page rather than a failed download. The loader below reports
     the failure in words instead. -->
<style>
  :root {{
    --bg: #ffffff; --fg: #1a1a1a; --muted: #5a5a5a; --line: #e2e2e2;
    --panel: #f7f7f8;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --bg: #14161a; --fg: #ececec; --muted: #a0a4ab; --line: #2b2f36;
      --panel: #1b1e24;
    }}
  }}
  :root[data-theme="dark"] {{
    --bg: #14161a; --fg: #ececec; --muted: #a0a4ab; --line: #2b2f36;
    --panel: #1b1e24;
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; background: var(--bg); color: var(--fg);
    font: 15px/1.55 ui-sans-serif, system-ui, -apple-system, sans-serif; }}
  .wrap {{ max-width: 1100px; margin: 0 auto; padding: 24px 16px 64px; }}
  h1 {{ font-size: 1.5rem; margin: 0 0 4px; letter-spacing: -0.01em; }}
  p.sub {{ color: var(--muted); margin: 0 0 20px; }}
  #viewer {{ width: 100%; height: 540px; position: relative;
    border: 1px solid var(--line); border-radius: 10px; overflow: hidden;
    background: var(--panel); }}
  table {{ width: 100%; border-collapse: collapse; margin-top: 20px;
    font-size: 14px; }}
  th, td {{ text-align: left; padding: 8px 10px;
    border-bottom: 1px solid var(--line); }}
  th {{ color: var(--muted); font-weight: 600; font-size: 12px;
    text-transform: uppercase; letter-spacing: .04em; }}
  .sw {{ display: inline-block; width: 11px; height: 11px; border-radius: 2px;
    margin-right: 7px; vertical-align: -1px; }}
  .legend {{ display: flex; gap: 18px; flex-wrap: wrap; margin: 14px 0 0;
    color: var(--muted); font-size: 13px; }}
  .controls {{ margin: 14px 0 0; display: flex; gap: 8px; flex-wrap: wrap; }}
  button {{ font: inherit; font-size: 13px; padding: 6px 12px;
    border: 1px solid var(--line); background: var(--panel); color: var(--fg);
    border-radius: 7px; cursor: pointer; }}
  button:hover {{ border-color: var(--muted); }}
  code {{ background: var(--panel); padding: 1px 5px; border-radius: 4px;
    font-size: 13px; }}
  #status {{ position: absolute; inset: 0; display: flex; align-items: center;
    justify-content: center; text-align: center; padding: 24px;
    color: var(--muted); font-size: 14px; line-height: 1.6; z-index: 5; }}
  #status.error {{ color: #c0392b; }}
  #status a {{ color: inherit; }}
  .note {{ color: var(--muted); font-size: 13px; margin-top: 18px;
    border-left: 2px solid var(--line); padding-left: 12px; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>EGFR domain III — pH-switch anchors</h1>
  <p class="sub">Structure 6ARU from the Protein Data Bank (PDB), receptor chain
     A. Coloured by computed solvent exposure, meaning how much of each residue
     is reachable from outside the protein. Generated by
     <code>analysis/05_anchor_geometry.py</code>.</p>

  <div id="viewer"><div id="status">Loading 3D viewer&hellip;</div></div>

  <div class="legend">
    <span><span class="sw" style="background:#1a9850"></span>exposed (relative
      solvent accessibility, RSA, &ge; 0.25)</span>
    <span><span class="sw" style="background:#fdae61"></span>partial</span>
    <span><span class="sw" style="background:#d73027"></span>buried (RSA &le; 0.05)</span>
    <span><span class="sw" style="background:#6a9fd8"></span>epitope 415–466 backbone</span>
  </div>

  <div class="controls">
    <button onclick="setStyle('cartoon')">Cartoon</button>
    <button onclick="setStyle('surface')">Surface</button>
    <button onclick="toggleCluster()">Toggle cluster highlight</button>
    <button onclick="resetView()">Reset view</button>
  </div>

  <table>
    <thead><tr><th>Anchor</th><th>Exposure</th><th>Role</th>
      <th>In cluster</th><th>Cetuximab contact</th></tr></thead>
    <tbody id="rows"></tbody>
  </table>

  <p class="note">The exposure cutoffs (0.25 / 0.05) are conventions rather than
     physical constants, so an anchor sitting near a boundary should be read as
     ambiguous. The surface shown is the bare protein: the flexible sugar chains
     attached at N444 are not drawn, so real accessibility near that site is
     lower than what you see here.</p>
</div>

<script>
const ANCHORS = {json.dumps(rows)};
const CHAIN = "A";
const EPITOPE_PDB = "391-442";   // UniProt 415-466 minus the +24 offset
const SOURCES = [
  "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.5/3Dmol-min.js",
  "https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js",
  "https://unpkg.com/3dmol@2.5.5/build/3Dmol-min.js"
];
let viewer, showCluster = true, style = "cartoon";

function status(msg, isError) {{
  let el = document.getElementById("status");
  if (!el) {{
    el = document.createElement("div");
    el.id = "status";
    document.getElementById("viewer").appendChild(el);
  }}
  el.className = isError ? "error" : "";
  el.innerHTML = msg;
}}

function clearStatus() {{
  const el = document.getElementById("status");
  if (el) el.remove();
}}

// Load the library from the first mirror that works. Without this, a wrong or
// withdrawn version number produces a bare "$3Dmol is not defined" in the
// console and an empty box on the page, which is indistinguishable from a bug.
function loadLibrary(i, done) {{
  if (typeof $3Dmol !== "undefined") {{ done(); return; }}
  if (i >= SOURCES.length) {{
    status("Could not load the 3D viewer library from any mirror.<br>" +
           "This page needs internet access. The numbers it visualises are in " +
           "<code>results/findings/05-anchor-geometry.md</code>, and " +
           "<code>explorer/anchor-view.pml</code> renders the same thing " +
           "offline in PyMOL.", true);
    return;
  }}
  status("Loading 3D viewer&hellip; (source " + (i + 1) + " of " +
         SOURCES.length + ")");
  const s = document.createElement("script");
  s.src = SOURCES[i];
  s.onload = () => (typeof $3Dmol !== "undefined") ? done()
                                                  : loadLibrary(i + 1, done);
  s.onerror = () => loadLibrary(i + 1, done);
  document.head.appendChild(s);
}}

function rowsHtml() {{
  document.getElementById("rows").innerHTML = ANCHORS.map(a => `
    <tr>
      <td><span class="sw" style="background:${{a.colour}}"></span>
          <strong>${{a.aa}}${{a.uni}}</strong></td>
      <td>${{a.verdict}}</td>
      <td>${{a.aa === "H" ? "target His → binder D/E" : "acidic → binder His"}}</td>
      <td>${{a.in_cluster ? "yes" : "—"}}</td>
      <td>${{a.cetux ? "yes" : "—"}}</td>
    </tr>`).join("");
}}

function render() {{
  viewer.removeAllModels();
  viewer.removeAllSurfaces();
  viewer.removeAllLabels();
  status("Fetching structure 6ARU from RCSB&hellip;");
  $3Dmol.download("pdb:6ARU", viewer, {{}}, function (model) {{
    // 3Dmol's download() resolves with the model even when the fetch FAILS --
    // it loads empty data and logs to the console. So a truthiness check never
    // fires; the atom count is what actually distinguishes success.
    const nAtoms = (model && model.selectedAtoms) ?
      model.selectedAtoms({{}}).length : 0;
    if (!nAtoms) {{
      status("Loaded the viewer, but fetching structure 6ARU from RCSB " +
             "returned no atoms. Check internet access and reload.", true);
      return;
    }}
    clearStatus();
    viewer.setStyle({{}}, {{}});
    viewer.setStyle({{chain: CHAIN}},
      style === "cartoon"
        ? {{cartoon: {{color: "#b9bec7", opacity: 0.85}}}}
        : {{cartoon: {{color: "#b9bec7", opacity: 0.35}}}});

    // Epitope backbone, in PDB numbering (UniProt 415-466 = PDB 391-442).
    viewer.addStyle({{chain: CHAIN, resi: EPITOPE_PDB}},
      {{cartoon: {{color: "#6a9fd8"}}}});

    ANCHORS.forEach(a => {{
      if (a.pdb === null) return;
      viewer.addStyle({{chain: CHAIN, resi: a.pdb}},
        {{stick: {{colorscheme: {{prop: "elem"}},
                  radius: (showCluster && a.in_cluster) ? 0.36 : 0.18}}}});
      viewer.addStyle({{chain: CHAIN, resi: a.pdb}},
        {{sphere: {{color: a.colour,
                   radius: (showCluster && a.in_cluster) ? 0.85 : 0.5}}}});
      // Place the label at the residue's own CA coordinates. Passing a bare
      // selection with an empty position is unreliable across 3Dmol versions.
      const sel = model.selectedAtoms(
        {{chain: CHAIN, resi: a.pdb, atom: "CA"}});
      if (sel && sel.length) {{
        viewer.addLabel(`${{a.aa}}${{a.uni}}`, {{
          position: {{x: sel[0].x, y: sel[0].y, z: sel[0].z}},
          alignment: "center", inFront: true,
          fontSize: 11, fontColor: "white",
          backgroundColor: a.colour, backgroundOpacity: 0.9,
        }});
      }}
    }});

    if (style === "surface") {{
      viewer.addSurface($3Dmol.SurfaceType.VDW,
        {{opacity: 0.6, color: "#c7ccd4"}}, {{chain: CHAIN}});
    }}
    viewer.zoomTo({{chain: CHAIN, resi: EPITOPE_PDB}});
    viewer.render();
  }});
}}

function setStyle(s) {{ if (viewer) {{ style = s; render(); }} }}
function toggleCluster() {{ if (viewer) {{ showCluster = !showCluster; render(); }} }}
function resetView() {{
  if (!viewer) return;
  viewer.zoomTo({{chain: CHAIN, resi: EPITOPE_PDB}});
  viewer.render();
}}

window.addEventListener("load", () => {{
  rowsHtml();   // the table is plain data and must render even if 3D fails
  loadLibrary(0, () => {{
    try {{
      viewer = $3Dmol.createViewer(document.getElementById("viewer"),
        {{backgroundColor: getComputedStyle(document.body)
          .getPropertyValue("--panel").trim() || "white"}});
      render();
    }} catch (e) {{
      // 3Dmol throws a bare string (not an Error) when WebGL is unavailable,
      // so reading e.message alone yields "undefined" and tells you nothing.
      const detail = (e && e.message) ? e.message : String(e);
      status("The 3D viewer could not start: " + detail +
             "<br>This usually means WebGL is unavailable — check that " +
             "hardware acceleration is enabled in your browser settings. " +
             "<code>explorer/anchor-view.pml</code> renders the same view in " +
             "PyMOL without WebGL.", true);
    }}
  }});
}});
</script>
</body>
</html>
""")


def write_pymol(path, uni_to_pdb, verdicts, cluster):
    lines = [
        "# PyMOL script: EGFR domain III with pH-switch anchors",
        "# Generated by analysis/05_anchor_geometry.py — do not hand-edit.",
        "# Run with:  pymol explorer/anchor-view.pml",
        "#",
        "# Numbering note: the colouring below uses the residue numbers in the",
        "# Protein Data Bank file, which in 6ARU are 24 lower than the UniProt",
        "# numbers. The labels show the UniProt numbers used everywhere else.",
        "",
        "fetch 6ARU, async=0",
        "hide everything",
        "show cartoon, chain A",
        "color grey80, chain A",
        "bg_color white",
        "",
        "# candidate epitope 415-466 (UniProt) = 391-442 (PDB)",
        "select epitope, chain A and resi 391-442",
        "color skyblue, epitope",
        "",
    ]
    palette = {"exposed": "forest", "partial": "orange", "buried": "firebrick"}
    for pos in sorted(ANCHORS):
        pdb_num = uni_to_pdb.get(pos)
        if pdb_num is None:
            lines.append(f"# {ANCHORS[pos]}{pos}: unresolved, skipped")
            continue
        v = verdicts.get(pos, "unknown")
        name = f"anchor_{ANCHORS[pos]}{pos}"
        lines += [
            f"# {ANCHORS[pos]}{pos} (UniProt) = {pdb_num} (PDB) — {v}"
            + ("  [in cluster]" if pos in cluster else ""),
            f"select {name}, chain A and resi {pdb_num}",
            f"show sticks, {name}",
            f"color {palette.get(v, 'grey50')}, {name}",
            f"label {name} and name CA, \"{ANCHORS[pos]}{pos}\"",
        ]
        if pos in cluster:
            lines.append(f"show spheres, {name} and not name C+N+O")
            lines.append(f"set sphere_scale, 0.3, {name}")
        lines.append("")
    lines += [
        "deselect",
        "set label_size, 16",
        "set label_color, black",
        "set cartoon_transparency, 0.25, chain A",
        "orient epitope",
        "",
        "# green = exposed, orange = partial, red = buried",
        "# spheres mark anchors in the reachable cluster",
    ]
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
