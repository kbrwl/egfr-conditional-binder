#!/usr/bin/env python3
"""Settle where the C-terminal trim boundary goes, by prediction rather than argument.

Step 09 cut the target down to domain III, our residues 310-480, and measured two
costs of that cut. One disulfide bond -- a chemical staple between the sulfur atoms
of two cysteine residues -- is severed, because C470 is inside the fragment and its
partner C499 is outside it. And E424, one of the eight anchors the binder is aimed
at, sits 6.6 angstroms from the cut at 480, which is the part of a trimmed model
least likely to be right.

Extending the boundary past 499 would address both at once: the staple becomes whole
and E424 gains margin from the new end. It costs about 20 more residues of target to
model, and it pulls in the start of domain IV, which may not fold properly on its
own and could make things worse rather than better.

That is an argument with two sides and no measurement in it, which is what this
script replaces. The method is the one already recorded against the Unverified entry
for the trim in docs/decisions-log.md: predict each candidate fragment on its own,
then compare it against the same residues in the intact measured structure. A
fragment that reproduces the intact structure is one a design made against it can be
trusted on. A fragment that does not is one where the designed binder was shaped
against a surface the real protein does not have.

Abbreviations, expanded here because this file gets read on its own. EGFR is the
epidermal growth factor receptor, the protein being designed against. UniProt is the
public archive of protein sequences, and every residue number here is a position in
its human record P00533 unless the line says otherwise. PDB is the Protein Data
Bank, the public archive of measured three-dimensional structures; 6ARU is the entry
this project uses as its structure of record. RMSD, root-mean-square deviation, is
the average distance between two structures' matching atoms once one has been moved
on top of the other as well as it can be -- so a small number means the same shape.

THE THREE CANDIDATES

  310-480   as step 09 cut it, with the severed staple and E424 near the end
  310-480   with C470 mutated to serine, which removes the unpaired cysteine at
            the cost of one deliberate difference from the real human sequence
  310-499   whole staple, both cysteines inside, and more margin for E424

The third boundary is 499 because that is where the staple's far end sits. The next
disulfide bond in the chain is C506-C515, so a boundary anywhere from 499 to 505
severs nothing; 499 is the cheapest of those and is the one measured here.

WHAT PREDICTS THE FRAGMENTS, AND WHAT THAT IS WORTH

ESMFold, reached through its public programming interface at api.esmatlas.com. It
predicts a structure from a single sequence, with no multiple sequence alignment --
no column of related proteins from other species to read conservation out of -- which
makes it meaningfully less accurate than AlphaFold2 in absolute terms. It is used
here for two reasons: it needs no graphics card and no account, and the project's
graphics-card access is still blocked, so the alternative is no measurement at all.

That limit governs how the numbers below may be read. An absolute RMSD from this
predictor is not evidence about how well the real fragment folds. What is usable is
the comparison between the three candidates, because all three go through the same
instrument with the same weaknesses, so a difference between them is attributable to
the boundary rather than to the predictor. Any conclusion here that rests on an
absolute value rather than on a difference is marked as such and should be treated as
unsettled.

HOW THE COMPARISON IS KEPT HONEST

The three fragments are not the same length, so an RMSD over each fragment's own
residues compares a 171-residue fit against a 190-residue fit, and the longer one is
penalised for residues the shorter one never had to get right. Every candidate is
therefore also measured over the common core, our 310-480, which all three contain.
The common-core number is the like-for-like comparison and is the one that decides
the question. Both are reported so neither is hidden.

Three superpositions are made per candidate, because they answer different questions
and a single number hides the one that matters:

  whole fragment    fitted on its own residues: does the fragment as a whole come
                    back the right shape
  common core       fitted on 310-480 only: the comparable number across candidates
  anchor face       fitted on the eight anchor residues alone: is the surface the
                    binder is actually designed against the right shape, which can
                    be true even when the fragment as a whole wobbles

Written by analysis/11_trim_boundary.py. Findings go to
results/findings/11-trim-boundary.md and tables to data/derived/.

    ./.venv/bin/python analysis/11_trim_boundary.py
    ./.venv/bin/python analysis/11_trim_boundary.py --self-test
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer

sys.path.insert(0, str(Path(__file__).resolve().parent))
import egfr_common as E  # noqa: E402

# The eight residues on one face of EGFR that step 08 chose and step 09 checked.
# Read back from step 09's committed table rather than retyped, so that if step 08's
# answer ever changes this script stops instead of quietly measuring the old set.
ANCHOR_EXPECTATION = [344, 358, 368, 370, 391, 400, 421, 424]

# The anchor the whole question is about: 6.6 angstroms from the cut at 480.
EDGE_ANCHOR = 424

# The severed staple. C470 is inside the 310-480 fragment, C499 outside it.
STAPLE = (470, 499)

# A disulfide bond measures about 2.05 angstroms between the two sulfur (SG) atoms.
# Anything under this counts as formed; the measured bonds in 6ARU are all 2.03-2.05.
SG_BOND_MAX = 2.5

# How many residues at each end of a fragment count as "the cut edge" for reporting
# the local damage. Ten is a judgement, chosen because it is about three turns of
# chain, not a measured threshold.
EDGE_WINDOW = 10

# Residues within this distance of E424 in the intact structure form its
# neighbourhood: the local surface a binder reaching for E424 would actually touch.
NEIGHBOURHOOD_RADIUS = 10.0

# ESMFold's public interface. Predicts from one sequence, no alignment, no account.
ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_RETRIES = 4
ESMFOLD_BACKOFF = 20          # seconds, doubled each retry

# Predictions are cached on disk so a rerun reproduces the same numbers without
# calling the interface again, and so the numbers do not silently change under a
# model update between runs. data/structures/ is not tracked; these are regenerable.
PREDICTION_DIR = E.STRUCT_DIR / "predictions"

# Set by --self-test to shift the predicted-to-our-numbering mapping by one, which
# must make the residue identity check below fail. A check never seen to fail is not
# known to work.
BREAK_MAPPING = False


def log(line=""):
    print(line)


# ----------------------------------------------------------------------------
# The intact structure, which is what everything is compared against
# ----------------------------------------------------------------------------

def load_intact():
    """6ARU's receptor chain, keyed by our residue numbers.

    Uses the offset analysis/02 measured and committed rather than an offset written
    in here, because an offset error raises nothing and returns a plausible residue
    in the wrong place.
    """
    path = E.STRUCT_DIR / "6aru_receptor_only.pdb"
    if not path.exists():
        raise SystemExit(
            f"{path} is missing. Run analysis/02_structure_prep.py first.")
    model = PDBParser(QUIET=True).get_structure("intact", str(path))[0]
    numbering = E.load_numbering()
    chain = next(iter(model))
    residues = {}
    for residue in E.protein_residues(chain):
        our = numbering.uniprot_of(residue.id[1])
        if our is not None:
            residues[our] = residue
    return model, residues, numbering


def fragment_sequence(intact, first, last, substitutions=None):
    """The one-letter sequence of our residues `first` to `last`, read out of the
    measured structure rather than out of the reference sequence, so that what is
    predicted is exactly the stretch of chain the design run would be handed.

    `substitutions` is {position: one-letter} for a deliberate mutation.
    """
    substitutions = substitutions or {}
    letters = []
    for position in range(first, last + 1):
        if position not in intact:
            raise SystemExit(
                f"residue {position} has no coordinates in 6ARU, so the fragment "
                f"{first}-{last} cannot be built from it")
        letter = E.THREE_TO_ONE.get(intact[position].get_resname(), "X")
        letters.append(substitutions.get(position, letter))
    return "".join(letters)


# ----------------------------------------------------------------------------
# Prediction
# ----------------------------------------------------------------------------

def predict(sequence, label, emit):
    """Predict one fragment's structure, caching the answer on disk.

    Returns the path to a PDB-format file whose B-factor column holds pLDDT, the
    predictor's own confidence at each residue, which ESMFold's interface reports on
    a 0 to 1 scale. It is rescaled to the conventional 0 to 100 elsewhere in this
    script and labelled where it is printed.
    """
    PREDICTION_DIR.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(sequence.encode()).hexdigest()[:16]
    cached = PREDICTION_DIR / f"{label}-{digest}.pdb"
    if cached.exists() and cached.stat().st_size > 0:
        emit(f"     {label}: cached prediction reused ({cached.name})")
        return cached

    emit(f"     {label}: predicting {len(sequence)} residues via ESMFold ...")
    body = sequence.encode()
    delay = ESMFOLD_BACKOFF
    last_error = None
    for attempt in range(1, ESMFOLD_RETRIES + 1):
        request = urllib.request.Request(
            ESMFOLD_URL, data=body,
            headers={"Content-Type": "text/plain"}, method="POST")
        try:
            with urllib.request.urlopen(request, timeout=600) as response:
                text = response.read().decode()
            if "ATOM" not in text:
                raise ValueError(f"no coordinates in reply: {text[:200]!r}")
            cached.write_text(text)
            emit(f"     {label}: predicted, written to {cached.name}")
            return cached
        except (urllib.error.URLError, ValueError, TimeoutError) as exc:
            last_error = exc
            if attempt < ESMFOLD_RETRIES:
                emit(f"     {label}: attempt {attempt} failed ({exc}); "
                     f"retrying in {delay}s")
                time.sleep(delay)
                delay *= 2
    raise SystemExit(
        f"ESMFold could not be reached after {ESMFOLD_RETRIES} attempts: "
        f"{last_error}.\nNothing is reported rather than reporting a boundary "
        f"choice made without the measurement it depends on.")


def load_prediction(path, first, intact, emit, label, substitutions=None):
    """Read a prediction and map its residues onto our numbering.

    The predictor is handed a bare sequence and numbers its answer from 1, so our
    position is `first + i - 1` for the i-th residue. That is arithmetic, and
    arithmetic is exactly what went wrong between steps 04 and 06, so every mapped
    position is checked by asking what amino acid it actually turns out to be and
    comparing that against the intact structure. A shift of one would leave the
    numbers looking reasonable and every residue identity wrong.
    """
    substitutions = substitutions or {}
    model = PDBParser(QUIET=True).get_structure("pred", str(path))[0]
    residues = E.protein_residues(next(iter(model)))
    offset = first + (1 if BREAK_MAPPING else 0)
    mapped, mismatches = {}, []
    for index, residue in enumerate(residues, start=1):
        our = offset + index - 1
        mapped[our] = residue
        predicted_aa = E.THREE_TO_ONE.get(residue.get_resname(), "X")
        if our in substitutions:
            expected = substitutions[our]
        elif our in intact:
            expected = E.THREE_TO_ONE.get(intact[our].get_resname(), "X")
        else:
            continue
        if predicted_aa != expected:
            mismatches.append((our, expected, predicted_aa))
    if mismatches:
        shown = ", ".join(f"{p}: expected {e}, got {g}"
                          for p, e, g in mismatches[:6])
        emit(f"   [FAIL] {label}: {len(mismatches)} residue(s) of the prediction "
             f"do not read as the amino acid at that position ({shown})")
        raise SystemExit(
            f"{label}: the prediction's residues do not line up with 6ARU's. "
            f"{len(mismatches)} mismatch(es), first few: {shown}. Refusing to "
            f"report distances measured between residues that are not the same "
            f"residue.")
    emit(f"   [PASS] {label}: all {len(mapped)} predicted residues read as the "
         f"expected amino acid, so the mapping is right by identity and not only "
         f"by arithmetic")
    return model, mapped


def plddt(residue):
    """The predictor's confidence at one residue, rescaled to 0-100.

    ESMFold's interface writes pLDDT into the B-factor column on a 0 to 1 scale.
    Above about 70 is usually read as a reliable backbone, below about 50 as
    essentially no information. Those thresholds are the field's convention rather
    than anything measured here.
    """
    values = [atom.get_bfactor() for atom in residue]
    return 100.0 * float(np.mean(values)) if values else float("nan")


# ----------------------------------------------------------------------------
# Comparison
# ----------------------------------------------------------------------------

def ca_pairs(mapped, intact, positions):
    """Matching alpha-carbon atoms for the positions present in both structures."""
    fixed, moving, used = [], [], []
    for position in positions:
        a, b = intact.get(position), mapped.get(position)
        if a is not None and b is not None and "CA" in a and "CA" in b:
            fixed.append(a["CA"])
            moving.append(b["CA"])
            used.append(position)
    return fixed, moving, used


def fit(model, mapped, intact, positions):
    """Move the prediction onto the intact structure using `positions` only, and
    return the RMSD of that fit plus the per-residue deviation everywhere.

    The deviations are measured after the move, so a residue's number says how far
    it sits from where the measured structure puts it once the two have been lined
    up as well as the chosen positions allow.
    """
    fixed, moving, used = ca_pairs(mapped, intact, positions)
    if len(fixed) < 3:
        return None, {}, []
    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    superimposer.apply(list(model.get_atoms()))
    deviations = {}
    for position, residue in mapped.items():
        other = intact.get(position)
        if other is not None and "CA" in residue and "CA" in other:
            deviations[position] = float(residue["CA"] - other["CA"])
    return float(superimposer.rms), deviations, used


def rmsd_over(deviations, positions):
    values = [deviations[p] for p in positions if p in deviations]
    if not values:
        return float("nan")
    return float(np.sqrt(np.mean(np.square(values))))


def charge_group_centre(residue):
    """The middle of a residue's charged group, which is the part a charge pair
    actually depends on. For glutamic and aspartic acid that is the pair of
    carboxyl oxygens; a contact measured between backbones can be right while the
    charged groups point away from each other.
    """
    names = {"GLU": ("OE1", "OE2"), "ASP": ("OD1", "OD2"),
             "HIS": ("ND1", "NE2")}.get(residue.get_resname())
    if not names:
        return None
    coords = [residue[n].coord for n in names if n in residue]
    return np.mean(coords, axis=0) if len(coords) == len(names) else None


def neighbourhood(intact, centre_position, radius):
    """Positions whose non-hydrogen atoms come within `radius` of the centre
    residue in the intact structure: the local surface around it."""
    centre = intact.get(centre_position)
    if centre is None:
        return []
    centre_atoms = [a for a in centre if a.element != "H"]
    near = []
    for position, residue in intact.items():
        if position == centre_position:
            continue
        for atom in residue:
            if atom.element == "H":
                continue
            if min(float(atom - c) for c in centre_atoms) <= radius:
                near.append(position)
                break
    return sorted(near)


def staple_distance(mapped):
    """Sulfur-to-sulfur distance of the C470-C499 staple in a prediction, or None
    when the fragment does not contain both cysteines."""
    a, b = mapped.get(STAPLE[0]), mapped.get(STAPLE[1])
    if a is None or b is None or "SG" not in a or "SG" not in b:
        return None
    return float(a["SG"] - b["SG"])


# ----------------------------------------------------------------------------
# One candidate, measured
# ----------------------------------------------------------------------------

def measure(candidate, intact, core, anchors, e424_neighbours, emit):
    label = candidate["label"]
    first, last = candidate["first"], candidate["last"]
    substitutions = candidate.get("substitutions") or {}

    sequence = fragment_sequence(intact, first, last, substitutions)
    path = predict(sequence, label, emit)
    model, mapped = load_prediction(path, first, intact, emit, label,
                                    substitutions)

    own = [p for p in range(first, last + 1)]
    result = {
        "label": label, "description": candidate["description"],
        "first": first, "last": last, "length": len(sequence),
        "sequence": sequence, "substitutions": substitutions,
        "prediction_file": path.name,
    }

    # 1. Fitted on the fragment's own residues: the fragment as a whole.
    rms_own, dev_own, _ = fit(model, mapped, intact, own)
    result["rmsd_whole_fragment"] = rms_own

    # 2. Fitted on the common core only: the comparable number across candidates.
    rms_core, dev_core, used_core = fit(model, mapped, intact, core)
    result["rmsd_common_core"] = rms_core
    result["core_residues_compared"] = len(used_core)
    result["anchor_rmsd_core_fit"] = rmsd_over(dev_core, anchors)
    result["anchor_deviations"] = {p: dev_core.get(p) for p in anchors}
    result["e424_deviation"] = dev_core.get(EDGE_ANCHOR)
    result["e424_neighbourhood_rmsd"] = rmsd_over(dev_core, e424_neighbours)

    # The last residues of this fragment, and the stretch 471-480 that every
    # candidate shares, so the cut's local damage is comparable as well as reported.
    own_edge = [p for p in range(max(first, last - EDGE_WINDOW + 1), last + 1)]
    shared_edge = [p for p in range(480 - EDGE_WINDOW + 1, 481)]
    result["own_c_edge_rmsd"] = rmsd_over(dev_core, own_edge)
    result["own_c_edge_range"] = (own_edge[0], own_edge[-1])
    result["shared_edge_rmsd"] = rmsd_over(dev_core, shared_edge)

    # 3. Fitted on the eight anchors alone: the surface designs are made against.
    rms_face, dev_face, used_face = fit(model, mapped, intact, anchors)
    result["anchor_face_rmsd"] = rms_face
    result["anchors_fitted"] = len(used_face)

    # Charge-group placement for E424, measured under the anchor-face fit, because
    # that is the frame a binder designed against this face would be built in.
    predicted_centre = charge_group_centre(mapped.get(EDGE_ANCHOR)) \
        if mapped.get(EDGE_ANCHOR) else None
    intact_centre = charge_group_centre(intact[EDGE_ANCHOR])
    if predicted_centre is not None and intact_centre is not None:
        result["e424_charge_group_shift"] = float(
            np.linalg.norm(predicted_centre - intact_centre))
    else:
        result["e424_charge_group_shift"] = None

    # Confidence, which is the predictor's own opinion rather than a comparison.
    result["plddt_mean"] = float(np.mean([plddt(r) for r in mapped.values()]))
    result["plddt_anchors"] = float(np.mean(
        [plddt(mapped[p]) for p in anchors if p in mapped]))
    result["plddt_e424"] = plddt(mapped[EDGE_ANCHOR]) \
        if EDGE_ANCHOR in mapped else None
    result["plddt_own_c_edge"] = float(np.mean(
        [plddt(mapped[p]) for p in own_edge if p in mapped]))

    # The staple itself.
    result["staple_sg_distance"] = staple_distance(mapped)
    cys470 = mapped.get(STAPLE[0])
    result["has_c470"] = cys470 is not None and cys470.get_resname() == "CYS"
    result["has_c499"] = STAPLE[1] in mapped

    return result


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------

def number(value, places=2, suffix=" A"):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "n/a"
    return f"{value:.{places}f}{suffix}"


def report(results, anchors, core, e424_neighbours, intact_e424_distance, emit):
    emit()
    emit("=" * 72)
    emit("WHERE THE C-TERMINAL TRIM BOUNDARY GOES")
    emit("=" * 72)
    emit()
    emit("Residue numbers are positions in the human record P00533 in UniProt,")
    emit("the public protein sequence archive. Each candidate fragment is")
    emit("predicted from its own sequence alone and compared against the same")
    emit("residues in 6ARU, the measured structure from the Protein Data Bank.")
    emit()
    emit("RMSD, root-mean-square deviation, is the average distance between two")
    emit("structures' matching atoms once one has been placed on top of the other")
    emit("as well as it can be. Smaller means the same shape. It is measured here")
    emit("between alpha carbons, one backbone atom per residue, so it reports the")
    emit("fold rather than the side chains.")
    emit()

    emit("1. What is being compared, and what the instrument is worth")
    emit()
    emit("   The predictor is ESMFold, reached through its public interface. It")
    emit("   predicts from a single sequence with no multiple sequence alignment,")
    emit("   which makes it less accurate than AlphaFold2 in absolute terms. It is")
    emit("   used because it needs no graphics card and no account, and this")
    emit("   project's graphics-card access is still blocked, so the alternative")
    emit("   is no measurement at all.")
    emit()
    emit("   So the absolute numbers below are NOT evidence about how well the real")
    emit("   fragment folds. The differences between the three candidates are the")
    emit("   usable result, because all three pass through the same instrument.")
    emit()
    for result in results:
        emit(f"     {result['label']}: our {result['first']}-{result['last']}, "
             f"{result['length']} residues — {result['description']}")
    emit()

    emit("2. Does the fragment come back the right shape?")
    emit()
    emit("   Two numbers per candidate, because the fragments are different")
    emit("   lengths and one number would not be a fair comparison.")
    emit()
    emit("   'whole fragment' fits each candidate on its own residues, so the")
    emit("   171-residue candidates are judged on 171 residues and the")
    emit("   190-residue one on 190. That is each fragment's own quality but it")
    emit("   is not like-for-like: the longer fragment is being asked to get more")
    emit("   residues right.")
    emit()
    emit("   'common core' fits every candidate on our 310-480 only, which all")
    emit("   three contain. This is the like-for-like comparison and it is the")
    emit("   number that decides the question.")
    emit()
    emit("   | candidate | whole fragment | common core 310-480 | residues in core |")
    emit("   |---|---|---|---|")
    for result in results:
        emit(f"   | {result['label']} | "
             f"{number(result['rmsd_whole_fragment'])} | "
             f"{number(result['rmsd_common_core'])} | "
             f"{result['core_residues_compared']} |")
    emit()

    emit("3. How well does it do across the eight anchors specifically?")
    emit()
    emit("   The anchors are the eight residues the binder is aimed at, so their")
    emit("   geometry matters more than the fragment's overall shape. Two")
    emit("   measurements, which answer different questions:")
    emit()
    emit("   'under the core fit' is how far the anchors sit from where 6ARU puts")
    emit("   them once the whole core has been lined up. It mixes the anchor")
    emit("   geometry together with any overall drift.")
    emit()
    emit("   'anchor face fit' lines the structures up on the eight anchors alone")
    emit("   and reports how well they can be made to agree. That is the number a")
    emit("   designed binder depends on: the binder is built against this face, so")
    emit("   what matters is whether the face has the right shape, which can hold")
    emit("   even when the rest of the fragment drifts.")
    emit()
    emit("   | candidate | anchors under the core fit | anchor face fit |")
    emit("   |---|---|---|")
    for result in results:
        emit(f"   | {result['label']} | "
             f"{number(result['anchor_rmsd_core_fit'])} | "
             f"{number(result['anchor_face_rmsd'])} |")
    emit()
    emit("   Per-anchor distance from 6ARU under the common-core fit:")
    emit()
    header = "   | candidate | " + " | ".join(f"{p}" for p in anchors) + " |"
    emit(header)
    emit("   |---" * (len(anchors) + 1) + "|")
    for result in results:
        cells = " | ".join(number(result["anchor_deviations"].get(p), 1, "")
                           for p in anchors)
        emit(f"   | {result['label']} | {cells} |")
    emit()

    emit("4. What happens around E424 and the C-terminal cut")
    emit()
    emit(f"   E424 is the anchor the question is about: step 09 measured it "
         f"{intact_e424_distance} from")
    emit("   the cut at 480, and every other anchor at 12.3 angstroms or more.")
    emit()
    emit("   Four measurements. E424's own distance from where 6ARU puts it. The")
    emit("   shift in its charged group, the pair of carboxyl oxygens that a")
    emit("   charge pair actually depends on, measured under the anchor-face fit")
    emit("   because that is the frame a binder would be built in. The RMSD of")
    emit(f"   E424's neighbourhood, the {len(e424_neighbours)} residues within "
         f"{NEIGHBOURHOOD_RADIUS:.0f} angstroms of it in")
    emit("   6ARU, which is the local surface a binder reaching for E424 touches.")
    emit("   And the last ten residues of each fragment, which is the cut edge")
    emit("   itself.")
    emit()
    emit("   | candidate | E424 | E424 charge group | E424 neighbourhood | "
         "own last 10 | shared 471-480 |")
    emit("   |---|---|---|---|---|---|")
    for result in results:
        emit(f"   | {result['label']} | "
             f"{number(result['e424_deviation'])} | "
             f"{number(result['e424_charge_group_shift'])} | "
             f"{number(result['e424_neighbourhood_rmsd'])} | "
             f"{number(result['own_c_edge_rmsd'])} "
             f"({result['own_c_edge_range'][0]}-{result['own_c_edge_range'][1]}) | "
             f"{number(result['shared_edge_rmsd'])} |")
    emit()
    emit("   The 'shared 471-480' column is the same ten residues for every")
    emit("   candidate, so it is the comparable one. For 310-480 those ten")
    emit("   residues are the fragment's own end; for 310-499 they sit 19")
    emit("   residues inside it, which is the whole point of extending.")
    emit()

    emit("5. The staple, and whether the prediction forms it")
    emit()
    emit("   A disulfide bond measures about 2.05 angstroms between the two sulfur")
    emit("   atoms, written SG. Below "
         f"{SG_BOND_MAX} angstroms counts as formed here.")
    emit()
    for result in results:
        if not result["has_c499"]:
            state = ("C499 is outside this fragment, so the staple cannot form "
                     "and C470 is unpaired")
        elif result["staple_sg_distance"] is None:
            state = "both cysteines present but no SG atoms to measure"
        else:
            distance = result["staple_sg_distance"]
            formed = "formed" if distance <= SG_BOND_MAX else "NOT formed"
            state = (f"C470 SG to C499 SG is {distance:.2f} angstroms, "
                     f"so the staple is {formed}")
        if not result["has_c470"]:
            state = ("C470 is mutated to serine, so there is no unpaired "
                     "cysteine and no staple to form")
        emit(f"     {result['label']}: {state}")
    emit()

    emit("6. The predictor's own confidence")
    emit()
    emit("   pLDDT is ESMFold's confidence at each residue, on a 0 to 100 scale")
    emit("   after rescaling from the 0 to 1 its interface reports. Above about 70")
    emit("   is usually read as a reliable backbone and below about 50 as")
    emit("   essentially no information; those are the field's conventions rather")
    emit("   than anything measured here. It is the model's opinion of itself, so")
    emit("   it is reported beside the comparison against 6ARU and not instead of")
    emit("   it: a prediction can be confident and wrong.")
    emit()
    emit("   | candidate | whole fragment | anchors | E424 | own last 10 |")
    emit("   |---|---|---|---|---|")
    for result in results:
        emit(f"   | {result['label']} | "
             f"{number(result['plddt_mean'], 1, '')} | "
             f"{number(result['plddt_anchors'], 1, '')} | "
             f"{number(result['plddt_e424'], 1, '')} | "
             f"{number(result['plddt_own_c_edge'], 1, '')} |")
    emit()


def write_tables(results, anchors):
    E.DERIVED.mkdir(parents=True, exist_ok=True)

    summary = E.DERIVED / "11-boundary-candidates.csv"
    columns = ["candidate", "first", "last", "length", "substitutions",
               "rmsd_whole_fragment_a", "rmsd_common_core_a",
               "core_residues_compared", "anchor_rmsd_core_fit_a",
               "anchor_face_rmsd_a", "e424_deviation_a",
               "e424_charge_group_shift_a", "e424_neighbourhood_rmsd_a",
               "own_c_edge_rmsd_a", "own_c_edge_first", "own_c_edge_last",
               "shared_edge_471_480_rmsd_a", "staple_sg_distance_a",
               "plddt_mean", "plddt_anchors", "plddt_e424", "plddt_own_c_edge",
               "prediction_file"]
    lines = [",".join(columns)]
    for r in results:
        subs = ";".join(f"{p}{a}" for p, a in sorted(r["substitutions"].items()))
        row = [r["label"], r["first"], r["last"], r["length"], subs or "none",
               fmt(r["rmsd_whole_fragment"]), fmt(r["rmsd_common_core"]),
               r["core_residues_compared"], fmt(r["anchor_rmsd_core_fit"]),
               fmt(r["anchor_face_rmsd"]), fmt(r["e424_deviation"]),
               fmt(r["e424_charge_group_shift"]),
               fmt(r["e424_neighbourhood_rmsd"]), fmt(r["own_c_edge_rmsd"]),
               r["own_c_edge_range"][0], r["own_c_edge_range"][1],
               fmt(r["shared_edge_rmsd"]), fmt(r["staple_sg_distance"]),
               fmt(r["plddt_mean"], 1), fmt(r["plddt_anchors"], 1),
               fmt(r["plddt_e424"], 1), fmt(r["plddt_own_c_edge"], 1),
               r["prediction_file"]]
        lines.append(",".join(str(c) for c in row))
    summary.write_text("\n".join(lines) + "\n")

    per_anchor = E.DERIVED / "11-boundary-anchor-deviations.csv"
    lines = ["candidate,uniprot_pos,amino_acid,ca_deviation_a,plddt"]
    for r in results:
        for position in anchors:
            lines.append(",".join(str(c) for c in [
                r["label"], position, E.ANCHORS.get(position, ""),
                fmt(r["anchor_deviations"].get(position)), ""]))
    per_anchor.write_text("\n".join(lines) + "\n")
    return summary, per_anchor


def fmt(value, places=3):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return f"{value:.{places}f}"


# ----------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Settle the C-terminal trim boundary by predicting each "
                    "candidate fragment and comparing it against 6ARU.")
    parser.add_argument("--self-test", action="store_true",
                        help="shift the predicted-to-our numbering mapping by "
                             "one and confirm the residue identity check stops "
                             "the run, because a check never seen to fail is "
                             "not known to work")
    args = parser.parse_args(argv)

    global BREAK_MAPPING

    lines = []

    def emit(line=""):
        lines.append(line)
        log(line)

    intact, intact_residues, _ = load_intact()
    del intact

    # The anchor set, read back from step 09's committed table rather than retyped.
    anchor_csv = E.DERIVED / "09-hotspot-numbering.csv"
    if not anchor_csv.exists():
        raise SystemExit(f"{anchor_csv} is missing. Run analysis/09 first.")
    anchors = []
    for line in anchor_csv.read_text().splitlines()[1:]:
        if line.strip():
            anchors.append(int(line.split(",")[0]))
    anchors.sort()
    if anchors != ANCHOR_EXPECTATION:
        raise SystemExit(
            f"the anchor set changed: step 09 gives {anchors}, this script "
            f"expected {ANCHOR_EXPECTATION}. Stopping rather than measuring a "
            f"boundary against the wrong residues.")

    if args.self_test:
        log()
        log("=" * 72)
        log("SELF-TEST: breaking the numbering mapping on purpose")
        log("=" * 72)
        log()
        log("The prediction is numbered from 1 and mapped onto our numbering by")
        log("arithmetic. If that arithmetic were off by one, every distance below")
        log("would be measured between two different residues and would still")
        log("look like a plausible number. The check that catches it compares the")
        log("amino acid at each mapped position against 6ARU. Here the mapping is")
        log("shifted by one on purpose; the check must stop the run.")
        log()
        BREAK_MAPPING = True
        candidate = {"label": "310-480", "first": E.D3_START, "last": E.D3_END,
                     "description": "as step 09 cut it"}
        try:
            measure(candidate, intact_residues, [], anchors, [], log)
        except SystemExit as exc:
            log()
            log("=" * 72)
            log("SELF-TEST PASSED. The shifted mapping was caught and the run")
            log("stopped. Reported reason:")
            log(f"  {str(exc).splitlines()[0]}")
            log("=" * 72)
            return 0
        log()
        log("=" * 72)
        log("SELF-TEST FAILED. A mapping shifted by one was not caught, so the")
        log("identity check is not protecting anything.")
        log("=" * 72)
        return 1

    core = list(range(E.D3_START, E.D3_END + 1))
    e424_neighbours = neighbourhood(intact_residues, EDGE_ANCHOR,
                                   NEIGHBOURHOOD_RADIUS)

    candidates = [
        {"label": "310-480", "first": 310, "last": 480,
         "description": "as step 09 cut it, severed staple, E424 near the end"},
        {"label": "310-480-C470S", "first": 310, "last": 480,
         "substitutions": {470: "S"},
         "description": "the unpaired cysteine mutated to serine, one "
                        "deliberate difference from the real sequence"},
        {"label": "310-499", "first": 310, "last": 499,
         "description": "extended past the staple's far end, whole staple, "
                        "more margin for E424"},
    ]

    log()
    log("=" * 72)
    log("PREDICTING THE THREE CANDIDATE FRAGMENTS")
    log("=" * 72)
    log()
    log("   Each is predicted from its sequence alone, with the answer cached on")
    log("   disk so a rerun reproduces these numbers without calling the")
    log("   interface again.")
    log()
    results = []
    for candidate in candidates:
        results.append(measure(candidate, intact_residues, core, anchors,
                               e424_neighbours, log))

    report(results, anchors, core, e424_neighbours, "6.6 angstroms", emit)
    summary, per_anchor = write_tables(results, anchors)

    emit("7. Files written")
    emit()
    emit(f"   {summary.relative_to(E.ROOT)}")
    emit(f"   {per_anchor.relative_to(E.ROOT)}")
    emit()

    best_core = min(results, key=lambda r: r["rmsd_common_core"])
    best_face = min(results, key=lambda r: r["anchor_face_rmsd"])
    emit("=" * 72)
    emit(f"RESULT: best common-core fit is {best_core['label']} at "
         f"{number(best_core['rmsd_common_core'])}; best anchor-face fit is "
         f"{best_face['label']} at {number(best_face['anchor_face_rmsd'])}.")
    emit("Whether those differences are large enough to decide the boundary is")
    emit("judged in results/findings/11-trim-boundary.md, not here.")
    emit("=" * 72)

    E.FINDINGS.mkdir(parents=True, exist_ok=True)
    out = E.FINDINGS / "11-trim-boundary.md"
    head = [
        "# Where the C-terminal trim boundary goes",
        "",
        "Computed output of `analysis/11_trim_boundary.py`. Do not hand-edit.",
        "",
        "Predicts each candidate target fragment from its sequence alone and",
        "compares it against the same residues in 6ARU, the measured structure,",
        "to decide whether the trim should stop at 480 as step 09 cut it or",
        "extend past the severed disulfide staple at 499.",
        "",
        "The predictor is ESMFold via its public interface, which works from a",
        "single sequence with no multiple sequence alignment and is therefore",
        "less accurate than AlphaFold2 in absolute terms. **The absolute numbers",
        "here are not evidence about how well the real fragment folds.** The",
        "differences between the three candidates are the usable result, because",
        "all three pass through the same instrument. This was the available",
        "instrument because the project's graphics-card access is still blocked.",
        "",
        "```",
    ]
    out.write_text("\n".join(head + lines + ["```", ""]))
    log()
    log(f"Wrote {out.relative_to(E.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
