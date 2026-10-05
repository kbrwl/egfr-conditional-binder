"""23_signal_peptide_distance.py -- how far is the signal peptide from our anchors?

WHY THIS EXISTS
---------------
On 4 October 2026 the organisers said the test article is the full ectodomain
starting at methionine 1, so residues 1-24 -- the signal peptide, the short
leader that is normally cut off when the protein is made -- are present in the
protein that goes on the instrument rather than removed from it. They added that
designs aimed near residue 25 should model 1-24.

That is their general statement about all entrants' designs. It is not a
measurement of ours. This step replaces it with one: how far is the start of the
mature region from each of the four anchors our designs actually pair against?

If the answer is tens of angstroms, the signal peptide is irrelevant to us
whatever it does, and the methods document can say so with a number instead of
repeating the organisers' caveat. If it is close, the design is exposed to a
stretch of protein nobody modelled.

Abbreviations, expanded because this file gets read on its own:
  EGFR  epidermal growth factor receptor, the protein we design against.
  PDB   the Protein Data Bank, the public archive of measured three-dimensional
        structures, and the older of its two file formats.
  aa    amino acids, the unit a protein's length is counted in.

WHICH STRUCTURE, AND WHAT IT CANNOT SHOW
-----------------------------------------
Measured on 1NQL, the tethered conformation, because that is the form the
organisers say the assay uses.

**No structure resolves the signal peptide.** 1-24 is absent from 1NQL and from
every other EGFR structure, because it is cleaved before the protein is
crystallised and is disordered in any case. So this cannot measure the distance
to residues 1-24 themselves. What it measures is the distance to the first
residue of the mature region that the structure does resolve, which is the point
the unmodelled leader would be attached to. The true distance from an anchor to
any part of the signal peptide is that distance plus however far the leader
extends, which is unknown and depends on a conformation nobody has measured.

This is therefore a lower bound, and it is reported as one. A lower bound is
enough to answer the question that matters: if the attachment point is already
far away, the leader hanging off it cannot reach.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import egfr_common as common

ROOT = Path(__file__).resolve().parents[1]
STRUCT = common.STRUCT_DIR
FINDINGS = common.FINDINGS
DERIVED = common.DERIVED

TETHERED = STRUCT / "1nql.pdb"
MATURE_START = 25           # first residue of the mature extracellular region
TIGHT4 = {344: "E", 358: "H", 368: "D", 370: "H"}

# The bar an anchor is judged against is how far the leader could physically
# reach, not a round number. The signal peptide is 24 residues and a polypeptide
# chain stretched out straight spans about 3.5 A a residue, so a fully extended
# leader reaches roughly 84 A from where it attaches. That is the worst case and
# a deliberately unkind one: a disordered chain samples compact shapes far more
# often than extended ones, and the fraction of time it spends anywhere near
# full extension is small.
#
# Using the worst case as the bar means a "cannot reach" verdict is strong and a
# "could reach" verdict is weak -- it says only that distance alone does not rule
# the contact out. That asymmetry is deliberate: the purpose here is to find out
# whether the organisers' caveat applies to our designs, and the honest way to
# retire a caveat is against the worst case rather than the likely one.
SIGNAL_PEPTIDE_RESIDUES = 24
EXTENDED_RESIDUE_SPAN_A = 3.5
CLEARLY_DISTANT_A = SIGNAL_PEPTIDE_RESIDUES * EXTENDED_RESIDUE_SPAN_A


def main():
    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    if not TETHERED.exists():
        sys.exit(f"{TETHERED} is missing. Run analysis/02 to fetch structures.")

    model = common.load_complex(TETHERED)
    human = common.human_sequence()
    chain_id, numbering, identity, _others = common.receptor_chain_of(model, human)
    chain = model[chain_id]
    residues = common.protein_residues(chain)

    resolved = sorted(p for p in (numbering.uniprot_of(r.id[1]) for r in residues)
                      if p is not None)
    first_resolved = resolved[0]

    emit("=" * 72)
    emit("SIGNAL PEPTIDE DISTANCE — how far is residue 25 from our anchors?")
    emit("=" * 72)
    emit()
    emit("The organisers confirmed on 4 October 2026 that the test article "
         "starts at")
    emit("methionine 1, so the signal peptide (residues 1-24) is present rather "
         "than")
    emit("cleaved, and said designs aimed near residue 25 should model it. This "
         "step")
    emit("measures where our own anchors sit relative to that end.")
    emit()
    emit(f"1. The structure: {TETHERED.name}, the tethered conformation, "
         f"chain {chain_id}")
    emit(f"   identified as EGFR by sequence identity, not by the chain letter: "
         f"{identity:.1f}%")
    emit()
    emit(f"   first residue resolved      {first_resolved} "
         f"(UniProt numbering)")
    emit(f"   mature region begins at     {MATURE_START}")
    emit(f"   residues {MATURE_START}-{first_resolved - 1} are not in the file"
         if first_resolved > MATURE_START else
         "   the mature start is resolved")
    emit()
    emit("   The signal peptide itself is resolved in no EGFR structure. The "
         "distances")
    emit("   below are to the first resolved residue, which is where the "
         "unmodelled")
    emit("   leader attaches, and are therefore LOWER BOUNDS on the distance to "
         "any")
    emit("   part of it.")
    emit()

    anchor_point = {}
    start_res = None
    for res in residues:
        pos = numbering.uniprot_of(res.id[1])
        if pos == first_resolved:
            start_res = res
        if pos in TIGHT4:
            # functional_point returns (coords, atom label, was it a fallback).
            # The label is carried through so the table says when a charged tip
            # was unresolved and a backbone atom stood in for it.
            anchor_point[pos] = common.functional_point(res, TIGHT4[pos])

    if start_res is None:
        sys.exit("could not locate the first resolved residue")

    start_atom = start_res["CA"] if "CA" in start_res else next(iter(start_res))

    emit("2. Distance from the mature N-terminus to each tight4 anchor")
    emit()
    emit(f"   Bar: {CLEARLY_DISTANT_A:.0f} A, a {SIGNAL_PEPTIDE_RESIDUES}-residue "
         f"leader at {EXTENDED_RESIDUE_SPAN_A} A a residue, fully extended.")
    emit()
    emit("   | anchor | distance (A) | reading |")
    emit("   |---|---|---|")
    rows = []
    for pos in sorted(TIGHT4):
        entry = anchor_point.get(pos)
        if entry is None:
            emit(f"   | {TIGHT4[pos]}{pos} | not resolved | cannot measure |")
            rows.append((pos, None, ""))
            continue
        point, atom_label, _fallback = entry
        d = float((((point - start_atom.get_coord()) ** 2).sum()) ** 0.5)
        reading = ("beyond even a fully extended leader" if d >= CLEARLY_DISTANT_A
                   else "within a fully extended leader's span")
        emit(f"   | {TIGHT4[pos]}{pos} ({atom_label}) | {d:.1f} | {reading} |")
        rows.append((pos, d, atom_label))
    emit()

    measured = [d for _, d, _ in rows if d is not None]
    nearest = min(measured) if measured else None

    emit("3. What this changes about the design")
    emit()
    if nearest is not None and nearest >= CLEARLY_DISTANT_A:
        emit(f"   The nearest of the four anchors is {nearest:.1f} A from the "
             f"point the signal")
        emit(f"   peptide attaches to, beyond the {CLEARLY_DISTANT_A:.0f} A a "
             "fully extended leader")
        emit("   could span. The organisers' caveat about modelling residues "
             "1-24 does not")
        emit("   apply to these designs, and the methods document states this "
             "distance")
        emit("   instead of repeating the caveat.")
    elif nearest is not None:
        emit(f"   The nearest anchor is {nearest:.1f} A from the signal "
             "peptide's attachment")
        emit(f"   point, inside the {CLEARLY_DISTANT_A:.0f} A a fully extended "
             f"{SIGNAL_PEPTIDE_RESIDUES}-residue leader could")
        emit("   span, so distance alone does not rule out contact with the "
             "anchor face.")
        emit("   It is a weak result in the unfavourable direction: full "
             "extension is the")
        emit("   worst case and a disordered chain rarely approaches it, and "
             "nothing here")
        emit("   measures what the leader actually does. The design does not "
             "account for")
        emit("   it and the write-up says so rather than treating the distance "
             "as clearance.")
    emit()
    emit("   This is a lower bound measured on one structure, and the leader's "
         "own")
    emit("   extent is unknown. It bounds the question rather than closing it.")
    emit()
    emit("=" * 72)
    emit(f"RESULT: nearest tight4 anchor is {nearest:.1f} A from the mature "
         f"N-terminus (lower bound)." if nearest else "RESULT: nothing measurable.")
    emit("=" * 72)

    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "23-signal-peptide-distance.csv").open("w") as fh:
        fh.write("anchor_uniprot_pos,anchor_aa,atom_measured_from,"
                 "distance_to_mature_n_terminus_a,"
                 "first_resolved_residue,structure\n")
        for pos, d, atom_label in rows:
            fh.write(f"{pos},{TIGHT4[pos]},{atom_label},"
                     f"{'' if d is None else f'{d:.3f}'},"
                     f"{first_resolved},{TETHERED.name}\n")
    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "23-signal-peptide-distance.md").write_text(
        "# How far is the signal peptide from our anchors?\n\n"
        "Computed output of `analysis/23_signal_peptide_distance.py`. "
        "Do not hand-edit.\n\n```\n" + "\n".join(out) + "\n```\n")
    print(f"\nWrote results/findings/23-signal-peptide-distance.md")
    print(f"Wrote data/derived/23-signal-peptide-distance.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
