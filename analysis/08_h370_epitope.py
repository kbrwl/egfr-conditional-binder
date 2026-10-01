#!/usr/bin/env python3
"""
08_h370_epitope.py — is an epitope built around H370 better than 415-466?

WHY THIS EXISTS
---------------
Liu et al. 2022 (Molecular Therapy - Oncolytics 27:256-269) mapped which
histidines on EGFR actually carry its pH-dependent antibody binding, by changing
each one to alanine and seeing which change removed the effect. Two came out:
H370 and H433.

Our epitope, 415-466, contains H433 and not H370. H433 is also the histidine that
paper used for its own design, so building on it puts us on top of published work.
H370 is the other experimentally implicated histidine, it sits outside our epitope,
and nobody has built a binder against it as far as we found. That makes it worth a
measurement.

Abbreviations used here, expanded on first use because this file gets read alone:
  PDB      Protein Data Bank, the public archive of measured 3D structures.
  UniProt  the public protein sequence archive; our residue numbers are positions
           in its records.
  SASA     solvent-accessible surface area, how much of a residue's surface water
           can reach, in square angstroms.
  RSA      relative solvent accessibility, that area divided by the largest value
           for that residue type, which makes residue types comparable.
  CB       the first carbon of a side chain.
  Fab      the gripping arm of an antibody, separated from the rest of it.
  aa       amino acids.

HOW THE CANDIDATE EPITOPE IS DEFINED, AND WHY NOT BY A SEQUENCE RANGE
---------------------------------------------------------------------
The obvious approach is to pick a window of the sequence around 370 and rerun the
earlier steps against it. We do report that window, because the conserved run
365-376 is real and step 01 already found it. But the window is not what decides
anything, for the reason this whole project exists: residues next to each other in
the sequence can point in opposite directions, so a window is not a patch.

So the candidate set is built from the structure instead. Take H370's imidazole
ring centre, collect every residue whose side-chain tip lies within the working
reach of a small binder, and then filter on the three things that actually matter:
identical in human and mouse, reachable by water, and of a type that can carry a
charge pair. A residue 40 positions away in the sequence that happens to sit beside
H370 in space is as useful as its immediate neighbour, and a sequence window would
miss it.

WHAT THIS ANSWERS
-----------------
  (a) which candidate anchors fall inside cetuximab's measured contact set
  (b) every human/mouse difference inside the candidate region
  (c) whether enough acidic anchors cluster within reach alongside H370 to support
      the three or four pairs the switch needs

Then it compares the result against 415-466 on the same terms.

Outputs:
  results/findings/08-h370-epitope.md
  data/derived/08-h370-neighbourhood.csv
  data/derived/08-h370-clusters.csv

Run standalone:  python analysis/08_h370_epitope.py
"""

import itertools
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

import egfr_common as common

ROOT = Path(__file__).resolve().parents[1]
RECEPTOR_PDB = ROOT / "data" / "structures" / "6aru_receptor_only.pdb"
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS

CENTRE = 370              # the histidine this epitope is built around
NEIGHBOURHOOD = 25.0      # angstroms from the centre, the working reach from step 05
REACH_CUTOFF = 25.0       # same working reach, used for clustering
MIN_CLUSTER = 3           # fewer than this and a stacked switch cannot be built

# The conserved sequence run containing H370, as step 01 found it.
RUN_START, RUN_END = 365, 376

FUNCTIONAL_ATOM = {"D": ["CG"], "E": ["CD"],
                   "H": ["CG", "ND1", "CD2", "CE1", "NE2"]}


def functional_point(residue, aa):
    """The position of the part of the side chain that forms a charge pair.

    For aspartate and glutamate that is the carboxylate carbon; for histidine the
    centre of the imidazole ring, because the charge is spread across the ring
    rather than sitting on one atom. Falls back to CB, the first side-chain carbon,
    when the tip is not resolved, and says so, because CB can be several angstroms
    from the group that actually does the work.
    """
    wanted = FUNCTIONAL_ATOM.get(aa, [])
    coords = [residue[a].coord for a in wanted if a in residue]
    if coords and len(coords) == len(wanted):
        if aa == "H":
            return np.mean(coords, axis=0), False
        return np.asarray(coords[0]), False
    if "CB" in residue:
        return np.asarray(residue["CB"].coord), True
    return np.asarray(residue["CA"].coord), True


def mouse_sequence_by_uniprot():
    """Mouse residue at each human UniProt position, via the challenge construct.

    The mouse construct starts at UniProt 25, so its position 1 is our 25. That
    +24 relationship holds upstream of a two-residue insertion near position 638,
    which is far downstream of anything here.
    """
    mouse = None
    for header, seq in common.read_fasta(common.CHALLENGE_FASTA).items():
        if "MOUSE" in header:
            mouse = seq
    if mouse is None:
        raise KeyError("no mouse challenge construct found")
    return lambda pos: mouse[pos - common.SIGNAL_PEPTIDE_LEN - 1]


def load_cetuximab_contacts():
    """EGFR positions the antibody touches, from step 04's computed table."""
    path = DERIVED / "04-cetuximab-contacts.csv"
    if not path.exists():
        raise SystemExit("Run analysis/04_cetuximab_contacts.py first.")
    contacts = {}
    for line in path.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        parts = line.split(",")
        contacts[int(parts[0])] = float(parts[3])
    return contacts


def main():
    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    human = common.human_sequence()
    mouse_at = mouse_sequence_by_uniprot()
    numbering = common.load_numbering()
    cetux = load_cetuximab_contacts()

    emit("=" * 72)
    emit(f"AN EPITOPE BUILT AROUND H{CENTRE}, COMPARED WITH 415-466")
    emit("=" * 72)
    emit()
    emit(f"H{CENTRE} is one of the two histidines Liu et al. 2022 showed experimentally")
    emit("to carry EGFR's pH-dependent antibody binding, by mutating each histidine")
    emit("to alanine. The other is H433, which sits in our current epitope and which")
    emit("that paper built its own design against.")
    emit()

    structure = PDBParser(QUIET=True).get_structure("receptor", str(RECEPTOR_PDB))
    model = structure[0]
    ShrakeRupley().compute(model, level="R")
    chain = next(iter(model))
    by_uniprot = {}
    for residue in common.protein_residues(chain):
        pos = numbering.uniprot_of(residue.id[1])
        if pos is not None:
            by_uniprot[pos] = residue

    centre_residue = by_uniprot.get(CENTRE)
    if centre_residue is None:
        emit(f"H{CENTRE} has no coordinates in this structure. Cannot proceed.")
        return 1
    centre_aa = common.THREE_TO_ONE.get(centre_residue.get_resname(), "X")
    if centre_aa != "H":
        emit(f"Position {CENTRE} reads {centre_aa}, not histidine. Stopping.")
        return 1
    centre_point, _ = functional_point(centre_residue, "H")

    # ---- (b) conservation, both for the sequence run and the spatial region ----
    emit(f"1. Species conservation (question b)")
    emit()
    emit(f"   The conserved sequence run containing H{CENTRE}, as step 01 found it:")
    emit(f"   {RUN_START}-{RUN_END}, {RUN_END - RUN_START + 1} residues, "
         f"{human[RUN_START - 1:RUN_END]}")
    run_diffs = [(p, human[p - 1], mouse_at(p))
                 for p in range(RUN_START, RUN_END + 1)
                 if human[p - 1] != mouse_at(p)]
    emit(f"   Human/mouse differences inside it: "
         f"{[f'{h}{p}{m}' for p, h, m in run_diffs] if run_diffs else 'none'}")
    emit()
    emit("   The residues immediately flanking that run, which is where it ends and")
    emit("   why:")
    for p in (RUN_START - 1, RUN_END + 1):
        emit(f"     {p}: human {human[p - 1]}, mouse {mouse_at(p)}"
             f"{'  <- differs' if human[p - 1] != mouse_at(p) else ''}")
    emit()

    # ---- the spatial neighbourhood ----
    emit(f"2. Every residue within {NEIGHBOURHOOD:.0f} A of H{CENTRE}")
    emit()
    emit("   Measured between side-chain tips, the parts that form a charge pair.")
    emit("   'conserved' compares human against the official mouse construct.")
    emit()
    rows = []
    for pos, residue in sorted(by_uniprot.items()):
        aa = common.THREE_TO_ONE.get(residue.get_resname(), "X")
        point, fallback = functional_point(residue, aa)
        distance = float(np.linalg.norm(point - centre_point))
        if distance > NEIGHBOURHOOD:
            continue
        max_asa = common.MAX_ASA_THEORETICAL.get(aa, 200.0)
        rsa = float(residue.sasa) / max_asa
        mouse_aa = mouse_at(pos)
        rows.append(dict(pos=pos, aa=aa, mouse=mouse_aa, conserved=(aa == mouse_aa),
                         distance=distance, rsa=rsa,
                         exposure=common.classify_exposure(rsa),
                         cetux=pos in cetux, fallback=fallback,
                         point=point))

    emit(f"   {len(rows)} residues within reach. Of those, the ones that could carry")
    emit("   a charge pair — acidic (D or E) or histidine:")
    emit()
    emit("   | pos | aa | mouse | conserved | dist from H370 | RSA | exposure | cetuximab |")
    emit("   |---|---|---|---|---|---|---|---|")
    switchable = [r for r in rows if r["aa"] in "DEH"]
    for r in sorted(switchable, key=lambda r: r["distance"]):
        emit(f"   | {r['pos']} | {r['aa']} | {r['mouse']} | "
             f"{'yes' if r['conserved'] else 'NO'} | {r['distance']:.1f} A | "
             f"{r['rsa']:.3f} | {r['exposure']} | "
             f"{'touched' if r['cetux'] else '-'} |")
    emit()

    # ---- candidate anchors ----
    candidates = [r for r in switchable
                  if r["conserved"] and r["exposure"] in ("exposed", "partial")]
    emit(f"3. Candidate anchors: conserved, reachable, and able to carry a pair")
    emit()
    if not candidates:
        emit("   None.")
    for r in sorted(candidates, key=lambda r: r["distance"]):
        role = common.anchor_role(r["aa"])
        emit(f"     {r['aa']}{r['pos']}: {r['distance']:.1f} A from H{CENTRE}, "
             f"RSA {r['rsa']:.3f} ({r['exposure']}), {role}")
    emit()
    acidic = [r for r in candidates if r["aa"] in "DE"]
    histidines = [r for r in candidates if r["aa"] == "H"]
    emit(f"   {len(acidic)} acidic, {len(histidines)} histidine.")
    emit()

    # ---- (a) cetuximab overlap ----
    emit("4. Cetuximab overlap (question a)")
    emit()
    touched = [r for r in candidates if r["cetux"]]
    emit(f"   Candidate anchors inside cetuximab's measured contact set: "
         f"{len(touched)} of {len(candidates)}")
    for r in touched:
        emit(f"     {r['aa']}{r['pos']} at {cetux[r['pos']]:.2f} A from the antibody")
    if not touched:
        emit("     none")
    emit()
    near_run = sorted(p for p in cetux if RUN_START - 6 <= p <= RUN_END + 6)
    emit(f"   For context, the contacts cetuximab makes near this region: {near_run}")
    emit(f"   H{CENTRE} itself: "
         f"{'touched' if CENTRE in cetux else 'not touched'} by cetuximab.")
    emit()

    # ---- (c) clustering ----
    emit(f"5. Do they cluster within {REACH_CUTOFF:.0f} A of each other? (question c)")
    emit()
    if len(candidates) < 2:
        emit("   Too few candidates to cluster.")
        best, best_span, valid = [], None, []
    else:
        points = {r["pos"]: r["point"] for r in candidates}
        aa_of = {r["pos"]: r["aa"] for r in candidates}
        dist = {}
        for a in points:
            for b in points:
                dist[(a, b)] = float(np.linalg.norm(points[a] - points[b]))
        emit("   Pairwise distances (angstroms):")
        ordered = sorted(points)
        emit("   | | " + " | ".join(f"{aa_of[p]}{p}" for p in ordered) + " |")
        emit("   |---|" + "---|" * len(ordered))
        for a in ordered:
            cells = ["—" if a == b else f"{dist[(a, b)]:.1f}" for b in ordered]
            emit(f"   | **{aa_of[a]}{a}** | " + " | ".join(cells) + " |")
        emit()
        valid, best = [], []
        for size in range(len(ordered), 1, -1):
            for combo in itertools.combinations(ordered, size):
                span = max(dist[(x, y)]
                           for x, y in itertools.combinations(combo, 2))
                if span < REACH_CUTOFF:
                    valid.append((combo, span))
                    if len(combo) > len(best):
                        best = list(combo)
            if best:
                break
        if best:
            best_span = max(dist[(x, y)]
                            for x, y in itertools.combinations(best, 2))
            emit(f"   Largest cluster: {len(best)} anchors, "
                 f"{', '.join(f'{aa_of[p]}{p}' for p in best)}, "
                 f"maximum internal distance {best_span:.1f} A")
            n_his = sum(1 for p in best if aa_of[p] == "H")
            emit(f"   Of those, {len(best) - n_his} acidic and {n_his} histidine.")
            has_centre = CENTRE in best
            emit(f"   Includes H{CENTRE}: {'yes' if has_centre else 'NO'}")
        else:
            best_span = None
            emit("   No two candidates fall within the cutoff.")
    emit()

    # ---- are they on one face, or wrapped around the protein? ----
    # The distance test above asks whether the anchors fit inside a ball roughly
    # 25 angstroms across. A binder presents a face, closer to a disc than a ball,
    # so two anchors on opposite sides of that ball are within the distance limit
    # and still unreachable together. With four anchors that gap between the test
    # and the thing being tested is small. With eight it is not, so it is measured
    # rather than assumed.
    emit("5b. Are they on one face, or wrapped around the protein?")
    emit()
    emit("   The distance test asks whether the anchors fit inside a ball. A binder")
    emit("   presents something closer to a flat face, so anchors on opposite sides")
    emit("   of that ball pass the distance test and are still unreachable together.")
    emit()
    emit("   Method: take the direction each anchor points away from the protein's")
    emit("   centre, then measure the angle between those directions. Anchors on one")
    emit("   face point roughly the same way. A wide spread means the set wraps")
    emit("   around the protein, however close together the distances look.")
    emit()
    all_ca = [r["CA"].coord for r in common.protein_residues(chain) if "CA" in r]
    centroid = np.mean(np.asarray(all_ca, dtype=float), axis=0)

    def outward(pos_point):
        v = np.asarray(pos_point, dtype=float) - centroid
        norm = np.linalg.norm(v)
        return v / norm if norm else v

    def angle_between(u, v):
        cosine = float(np.clip(np.dot(u, v), -1.0, 1.0))
        return float(np.degrees(np.arccos(cosine)))

    FACE_LIMIT = 90.0   # widest angular spread we treat as reachable by one face
    if best:
        normals = {p: outward(points[p]) for p in best}
        worst = max(
            (angle_between(normals[a], normals[b]), a, b)
            for a, b in itertools.combinations(best, 2))
        emit(f"   Widest angle within the {len(best)}-anchor cluster: "
             f"{worst[0]:.0f} degrees, between {aa_of[worst[1]]}{worst[1]} and "
             f"{aa_of[worst[2]]}{worst[2]}")
        emit()
        # Largest subset that passes distance AND angular spread.
        face_best, face_valid = [], []
        for size in range(len(best), 1, -1):
            for combo in itertools.combinations(sorted(best), size):
                span = max(dist[(x, y)]
                           for x, y in itertools.combinations(combo, 2))
                spread = max(angle_between(normals[x], normals[y])
                             for x, y in itertools.combinations(combo, 2))
                if span < REACH_CUTOFF and spread < FACE_LIMIT:
                    face_valid.append((combo, span, spread))
                    if len(combo) > len(face_best):
                        face_best = list(combo)
            if face_best:
                break
        if face_best:
            fspan = max(dist[(x, y)]
                        for x, y in itertools.combinations(face_best, 2))
            fspread = max(angle_between(normals[x], normals[y])
                          for x, y in itertools.combinations(face_best, 2))
            emit(f"   Largest subset that also fits one face, under "
                 f"{FACE_LIMIT:.0f} degrees of spread:")
            emit(f"     {len(face_best)} anchors, "
                 f"{', '.join(f'{aa_of[p]}{p}' for p in face_best)}")
            emit(f"     span {fspan:.1f} A, angular spread {fspread:.0f} degrees")
            n_his_face = sum(1 for p in face_best if aa_of[p] == "H")
            emit(f"     {len(face_best) - n_his_face} acidic, "
                 f"{n_his_face} histidine; includes H{CENTRE}: "
                 f"{'yes' if CENTRE in face_best else 'no'}")
            if len(face_best) < len(best):
                emit()
                emit(f"   So the honest figure is {len(face_best)} rather than "
                     f"{len(best)}. The distance test alone overstates it, because")
                emit("   some of those anchors sit round the curve of the protein.")
        else:
            emit(f"   No subset of two or more passes both tests. The cluster is")
            emit("   spread around the protein rather than sitting on one face.")
        emit()
        emit(f"   The {FACE_LIMIT:.0f} degree limit is a working convention chosen")
        emit("   here, not a measured property of binders, and it is cruder than")
        emit("   docking a real backbone against the surface. Treat it as a filter")
        emit("   that removes clearly unreachable sets rather than as a guarantee")
        emit("   about the ones that remain.")
    else:
        face_best, face_valid = [], []
        emit("   No cluster to test.")
    emit()

    # ---- verdict and comparison ----
    emit("6. Verdict, and how it compares with 415-466")
    emit()
    n_best = len(best)
    viable = n_best >= MIN_CLUSTER and CENTRE in best
    # Apply the same face test to the 415-466 clusters, so the two epitopes are
    # judged on identical terms. Without this we would be holding the new
    # candidate to a stricter standard than the incumbent.
    emit("   The same face test applied to the 415-466 clusters, so both epitopes")
    emit("   are judged on identical terms:")
    emit()
    incumbent_sets = {
        "D416 E421 E424 E455": [416, 421, 424, 455],
        "E424 E455 D458 D460": [424, 455, 458, 460],
        "H433 E455 D458 D460": [433, 455, 458, 460],
        "D416 E424 E455 D460": [416, 424, 455, 460],
        "D416 H418 E421 E424 E455 (conditional on H418)":
            [416, 418, 421, 424, 455],
    }
    for label, members in incumbent_sets.items():
        pts, ok = {}, True
        for pos in members:
            residue = by_uniprot.get(pos)
            if residue is None:
                ok = False
                break
            aa = common.THREE_TO_ONE.get(residue.get_resname(), "X")
            pt, _fb = functional_point(residue, aa)
            pts[pos] = pt
        if not ok:
            emit(f"     {label}: not all residues resolved")
            continue
        norms = {p: outward(pts[p]) for p in pts}
        span = max(float(np.linalg.norm(pts[a] - pts[b]))
                   for a, b in itertools.combinations(members, 2))
        spread = max(angle_between(norms[a], norms[b])
                     for a, b in itertools.combinations(members, 2))
        verdict = "one face" if spread < FACE_LIMIT else "WRAPS AROUND"
        emit(f"     {label}: span {span:.1f} A, spread {spread:.0f} deg — {verdict}")
    emit()

    emit("   415-466, from steps 03 and 05:")
    emit("     52 residues, one human/mouse difference (S442G)")
    emit("     7 of 8 anchors reachable, largest cluster 4 anchors within 22.6 A")
    emit("     contains H433, which cetuximab touches and which Liu et al. used")
    emit("     two anchors, D458 and D460, sit in the groove against domain IV")
    emit("     a sugar chain is attached at N444, inside the block")
    emit()
    emit(f"   An epitope around H{CENTRE}, from this script:")
    emit(f"     conserved run {RUN_START}-{RUN_END}, "
         f"{len(run_diffs)} human/mouse difference(s) inside it")
    emit(f"     {len(candidates)} candidate anchors conserved and reachable")
    emit(f"     largest cluster {n_best} anchors"
         + (f" within {best_span:.1f} A" if best_span is not None else ""))
    emit(f"     contains H{CENTRE}: {'yes' if CENTRE in best else 'no'}")
    emit(f"     candidate anchors cetuximab touches: {len(touched)}")
    emit()
    if viable:
        emit(f"   This epitope is viable: {n_best} anchors cluster within reach and")
        emit(f"   H{CENTRE} is among them, against a minimum of {MIN_CLUSTER}.")
    elif n_best >= MIN_CLUSTER:
        emit(f"   {n_best} anchors cluster, but H{CENTRE} is not among them, so the")
        emit("   experimentally implicated histidine cannot be paired against from")
        emit("   this cluster. That removes the reason for preferring this epitope.")
    else:
        emit(f"   Not viable: only {n_best} anchor(s) cluster within reach, against a")
        emit(f"   minimum of {MIN_CLUSTER}. A stacked switch cannot be built here.")
    emit()
    emit("   The argument for this epitope, if it is viable: H433 is the histidine")
    emit("   Liu et al. designed against, so a binder built on it sits on top of")
    emit(f"   published work, while H{CENTRE} is the other histidine their experiment")
    emit("   implicated and no published binder uses it. Both are experimentally")
    emit("   supported; only one is taken.")
    emit()
    emit("   What this does not settle. Cluster geometry is a necessary condition and")
    emit("   not a sufficient one, exactly as in step 05. It says nothing about")
    emit("   whether a foldable binder can present the right partners in the right")
    emit("   orientations, and nothing about whether the switch clears the assay's")
    emit("   detection floor.")
    emit()

    # ---- CSVs ----
    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "08-h370-neighbourhood.csv").open("w") as fh:
        fh.write("uniprot_pos,human_aa,mouse_aa,conserved,distance_from_h370_a,"
                 "sasa_a2,rsa,exposure,cetuximab_contact,switchable,"
                 "candidate_anchor,cb_fallback\n")
        cand_pos = {r["pos"] for r in candidates}
        for r in sorted(rows, key=lambda r: r["pos"]):
            residue = by_uniprot[r["pos"]]
            fh.write(f"{r['pos']},{r['aa']},{r['mouse']},{r['conserved']},"
                     f"{r['distance']:.3f},{float(residue.sasa):.2f},"
                     f"{r['rsa']:.4f},{r['exposure']},{r['cetux']},"
                     f"{r['aa'] in 'DEH'},{r['pos'] in cand_pos},"
                     f"{r['fallback']}\n")
    with (DERIVED / "08-h370-clusters.csv").open("w") as fh:
        fh.write("cluster_size,max_internal_span_a,anchors,contains_h370\n")
        for combo, span in sorted(valid, key=lambda t: (-len(t[0]), t[1])):
            fh.write(f"{len(combo)},{span:.3f},"
                     f"\"{' '.join(str(p) for p in combo)}\","
                     f"{CENTRE in combo}\n")
    emit("Wrote data/derived/08-h370-neighbourhood.csv")
    emit("Wrote data/derived/08-h370-clusters.csv")

    emit()
    emit("=" * 72)
    emit(f"RESULT: {len(candidates)} candidate anchors around H{CENTRE}; "
         f"largest cluster {n_best}"
         + (f" within {best_span:.1f} A" if best_span is not None else "")
         + f"; includes H{CENTRE}: {CENTRE in best}.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "08-h370-epitope.md").write_text(
        f"# An epitope built around H{CENTRE}\n\n"
        "Computed output of `analysis/08_h370_epitope.py`. Do not hand-edit.\n\n"
        f"H{CENTRE} is one of two histidines that Liu et al. 2022 (Molecular Therapy\n"
        "- Oncolytics 27:256-269, doi:10.1016/j.omto.2022.11.001) showed carry EGFR's\n"
        "pH-dependent antibody binding, found by mutating each histidine to alanine.\n"
        "The other is H433, which sits in our current epitope and which that paper\n"
        "built its own design against.\n\n"
        "The candidate set is defined from the structure rather than from a sequence\n"
        "window, because residues adjacent in sequence can point in opposite\n"
        "directions and a window is therefore not a patch.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
