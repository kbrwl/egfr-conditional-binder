#!/usr/bin/env python3
"""
03_solvent_accessibility.py — which anchors are actually on the surface?

WHAT THIS MEASURES AND WHY IT MATTERS
-------------------------------------
Solvent accessibility is how much of a residue is exposed to the surrounding
water. It is computed by rolling a water-sized probe ball over the protein's
surface and measuring the area that probe can touch for each residue. A residue
with a large accessible area sticks out and a binder can reach it. A residue with
an area near zero is buried inside the protein's core: it exists, it is
conserved, and it is completely useless to us, because nothing can touch it.

This is the step that decides how much material we actually have. Our eight
anchors were chosen by sequence analysis, which cannot see direction. Domain III
is a solenoid -- a spiral staircase -- so some of those eight point outward and
some point into the core. Only the outward ones are real.

WHY THE RECEPTOR ALONE, NOT THE COMPLEX
---------------------------------------
This runs on the receptor-only structure written by step 02, with the cetuximab
Fab removed. That is essential, not tidiness. Cetuximab is sitting directly on
the surface we care about. Computing on the complex would measure our epitope as
buried -- but it is not buried, it is merely covered by an antibody that will not
be anywhere near our assay. We would reject a perfectly good epitope for the
wrong reason.

RAW AREA IS NOT COMPARABLE BETWEEN RESIDUES
-------------------------------------------
Amino acids differ enormously in size: tryptophan is far larger than glycine. So
50 square angstroms of exposed tryptophan is mostly buried, while 50 square
angstroms of exposed glycine is wide open. The raw number therefore cannot be
compared across residue types.

The fix is RELATIVE solvent accessibility (RSA): divide the measured area by the
largest area that amino acid type could possibly expose. The result runs roughly
0 to 1 and IS comparable. RSA is the number to read; raw area is not.

REFERENCE MAXIMA — SOURCE AND STATUS
------------------------------------
Maxima are from Table 1 of:

  Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013)
  "Maximum Allowed Solvent Accessibilites of Residues in Proteins."
  PLOS ONE 8(11): e80635. doi:10.1371/journal.pone.0080635

Both published columns are transcribed below, fetched from the journal's
manuscript XML on 30 September 2026 rather than recalled from memory. The
THEORETICAL column is used for the headline numbers. The EMPIRICAL column is then
used as a sensitivity check: if a residue's classification changes depending on
which column you normalise by, that residue's verdict is not robust and the
script says so rather than presenting one of the two as the answer.

Several such tables exist in the literature and they disagree by up to 20%, so
which one was used has to be stated for any RSA number to be reproducible.

CLASSIFICATION CUTOFFS ARE CONVENTIONS, NOT PHYSICS
---------------------------------------------------
  RSA >= 0.25  exposed
  RSA <= 0.05  buried
  in between   partially exposed

Nothing physical changes in the molecule at 0.25. These are conventional
thresholds that the field has settled on for convenience. A residue at 0.24 and
one at 0.26 are essentially the same thing. This script therefore flags any
residue within 0.05 of a cutoff as BORDERLINE instead of letting the label imply
more precision than exists.

Outputs:
  results/findings/03-solvent-accessibility.md
  data/derived/03-epitope-rsa.csv
  data/derived/03-anchor-verdicts.csv

Run standalone:  python analysis/03_solvent_accessibility.py
"""

import sys
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

ROOT = Path(__file__).resolve().parents[1]
STRUCT_DIR = ROOT / "data" / "structures"
RECEPTOR_PDB = STRUCT_DIR / "6aru_receptor_only.pdb"
COMPLEX_PDB = STRUCT_DIR / "6aru.pdb"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"

# Offset established empirically by step 02: UniProt = PDB residue number + 24.
# Not assumed -- read back from step 02's output so the two cannot drift apart.
OFFSET_CSV = DERIVED / "02-numbering-offset.csv"

EPI_START, EPI_END = 415, 466
ANCHORS = {416: "D", 418: "H", 421: "E", 424: "E",
           433: "H", 455: "E", 458: "D", 460: "D"}

EXPOSED_CUT, BURIED_CUT = 0.25, 0.05
BORDERLINE_MARGIN = 0.05

THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}

# Tien et al. 2013, PLOS ONE 8(11): e80635, Table 1.
# Verified against the published manuscript XML, 30 September 2026.
MAX_ASA_THEORETICAL = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}
MAX_ASA_EMPIRICAL = {
    "A": 121.0, "R": 265.0, "N": 187.0, "D": 187.0, "C": 148.0,
    "E": 214.0, "Q": 214.0, "G": 97.0, "H": 216.0, "I": 195.0,
    "L": 191.0, "K": 230.0, "M": 203.0, "F": 228.0, "P": 154.0,
    "S": 143.0, "T": 163.0, "W": 264.0, "Y": 255.0, "V": 165.0,
}

# Sugar residue names to look for when checking glycan occlusion.
GLYCANS = {"NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "SIA", "BGC", "GLC"}


def classify(rsa):
    if rsa >= EXPOSED_CUT:
        return "exposed"
    if rsa <= BURIED_CUT:
        return "buried"
    return "partial"


def is_borderline(rsa):
    return (abs(rsa - EXPOSED_CUT) < BORDERLINE_MARGIN
            or abs(rsa - BURIED_CUT) < BORDERLINE_MARGIN)


def load_offset():
    """Read the offset from step 02's computed output rather than hardcoding it."""
    if not OFFSET_CSV.exists():
        raise SystemExit("Run analysis/02_structure_prep.py first — "
                         f"{OFFSET_CSV.relative_to(ROOT)} is missing.")
    offsets = {}
    for line in OFFSET_CSV.read_text().splitlines()[1:]:
        uni, pdb, off = line.split(",")
        offsets[int(pdb)] = int(uni)
    return offsets


def main():
    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    if not RECEPTOR_PDB.exists():
        raise SystemExit("Run analysis/02_structure_prep.py first — "
                         "receptor-only structure is missing.")

    pdb_to_uniprot = load_offset()

    emit("=" * 72)
    emit("SOLVENT ACCESSIBILITY OF THE CANDIDATE EPITOPE")
    emit("=" * 72)
    emit()
    emit("Input:  data/structures/6aru_receptor_only.pdb  (Fab REMOVED)")
    emit("Method: Shrake-Rupley probe rolling, Bio.PDB.SASA.ShrakeRupley")
    emit("Norm:   Tien et al. 2013, PLOS ONE 8(11):e80635, Table 1")
    emit("        theoretical column for headline values,")
    emit("        empirical column as a sensitivity check")
    emit(f"Numbering: UniProt = PDB + 24 (established empirically by step 02)")
    emit()
    emit("The Fab is removed deliberately. Cetuximab sits on the surface we are")
    emit("measuring, so computing on the complex would report our epitope as")
    emit("buried when it is only covered by an antibody absent from our assay.")
    emit()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(RECEPTOR_PDB))
    sr = ShrakeRupley()
    sr.compute(structure[0], level="R")

    chain = next(iter(structure[0]))
    records = {}
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        pdb_num = res.id[1]
        uni = pdb_to_uniprot.get(pdb_num)
        if uni is None:
            continue
        aa = THREE_TO_ONE.get(res.get_resname(), "X")
        sasa = float(res.sasa)
        rsa_t = sasa / MAX_ASA_THEORETICAL[aa]
        rsa_e = sasa / MAX_ASA_EMPIRICAL[aa]
        records[uni] = dict(pdb_num=pdb_num, aa=aa, sasa=sasa,
                            rsa_t=rsa_t, rsa_e=rsa_e,
                            cls_t=classify(rsa_t), cls_e=classify(rsa_e))

    emit(f"Computed SASA for {len(records)} residues in the receptor chain.")
    emit()

    # ---- Full table for the epitope ----
    emit(f"1. Every residue in {EPI_START}-{EPI_END}")
    emit()
    emit("   RSA is the number to read. 'anchor' marks our eight candidates.")
    emit("   'BORDERLINE' means within 0.05 of a classification cutoff, so the")
    emit("   label should not be read as precise.")
    emit()
    emit("   | UniProt | aa | PDB# | SASA A^2 | RSA | class | flags |")
    emit("   |---|---|---|---|---|---|---|")
    missing = []
    for pos in range(EPI_START, EPI_END + 1):
        r = records.get(pos)
        if r is None:
            missing.append(pos)
            emit(f"   | {pos} | ? | — | — | — | UNRESOLVED | no coordinates |")
            continue
        flags = []
        if pos in ANCHORS:
            flags.append("ANCHOR")
        if is_borderline(r["rsa_t"]):
            flags.append("BORDERLINE")
        if r["cls_t"] != r["cls_e"]:
            flags.append(f"NORM-SENSITIVE({r['cls_e']} if empirical)")
        emit(f"   | {pos} | {r['aa']} | {r['pdb_num']} | {r['sasa']:7.1f} | "
             f"{r['rsa_t']:.3f} | {r['cls_t']} | {', '.join(flags)} |")
    emit()

    if missing:
        emit(f"   {len(missing)} residue(s) unresolved: {missing}")
        emit("   No coordinates means no answer. Not 'buried' -- unknown.")
        emit()

    # Epitope-level summary
    present = [records[p] for p in range(EPI_START, EPI_END + 1) if p in records]
    n_exp = sum(1 for r in present if r["cls_t"] == "exposed")
    n_par = sum(1 for r in present if r["cls_t"] == "partial")
    n_bur = sum(1 for r in present if r["cls_t"] == "buried")
    emit(f"   Epitope summary ({len(present)} resolved): "
         f"{n_exp} exposed, {n_par} partial, {n_bur} buried")
    emit(f"   Mean RSA across the block: "
         f"{sum(r['rsa_t'] for r in present) / len(present):.3f}")
    emit()

    # ---- Anchor verdicts ----
    emit("2. VERDICT ON EACH OF THE EIGHT ANCHORS")
    emit()
    emit("   This determines how much material we have left to work with.")
    emit()
    emit("   | anchor | role | SASA A^2 | RSA (theor.) | RSA (emp.) | VERDICT | usable? |")
    emit("   |---|---|---|---|---|---|---|")
    anchor_rows, usable = [], []
    for pos in sorted(ANCHORS):
        want = ANCHORS[pos]
        role = ("acidic -> binder HIS" if want in "DE"
                else "target HIS -> binder D/E")
        r = records.get(pos)
        if r is None:
            emit(f"   | {want}{pos} | {role} | — | — | — | UNRESOLVED | no |")
            anchor_rows.append((pos, want, role, "", "", "", "unresolved", False))
            continue
        assert r["aa"] == want, f"anchor identity mismatch at {pos}"
        verdict = r["cls_t"].upper()
        note = []
        if is_borderline(r["rsa_t"]):
            note.append("BORDERLINE")
        if r["cls_t"] != r["cls_e"]:
            note.append("NORM-SENSITIVE")
        # Usable = exposed or partial. Buried anchors are dead.
        ok = r["cls_t"] in ("exposed", "partial")
        if ok:
            usable.append(pos)
        emit(f"   | {want}{pos} | {role} | {r['sasa']:7.1f} | {r['rsa_t']:.3f} | "
             f"{r['rsa_e']:.3f} | {verdict}{' *' + ','.join(note) if note else ''} "
             f"| {'YES' if ok else 'no'} |")
        anchor_rows.append((pos, want, role, f"{r['sasa']:.1f}",
                            f"{r['rsa_t']:.4f}", f"{r['rsa_e']:.4f}",
                            r["cls_t"], ok))
    emit()

    n_usable = len(usable)
    emit(f"   Anchors surviving as exposed or partially exposed: {n_usable} of 8")
    emit(f"   Survivors: {', '.join(f'{ANCHORS[p]}{p}' for p in usable)}")
    dead = [p for p in sorted(ANCHORS) if p not in usable]
    if dead:
        emit(f"   Lost to burial or absence: "
             f"{', '.join(f'{ANCHORS[p]}{p}' for p in dead)}")
    emit()

    # ---- Sensitivity summary ----
    norm_sensitive = [p for p in sorted(ANCHORS)
                      if p in records and records[p]["cls_t"] != records[p]["cls_e"]]
    borderline = [p for p in sorted(ANCHORS)
                  if p in records and is_borderline(records[p]["rsa_t"])]
    emit("3. How solid are these verdicts?")
    emit()
    emit("   Two ways a verdict could be an artefact of an arbitrary choice:")
    emit()
    emit(f"   a) Which reference table we normalise by.")
    if norm_sensitive:
        emit(f"      {len(norm_sensitive)} anchor(s) change class between the")
        emit(f"      theoretical and empirical columns: "
             f"{', '.join(f'{ANCHORS[p]}{p}' for p in norm_sensitive)}")
        emit("      Their verdicts are NOT robust and should be treated as")
        emit("      ambiguous rather than settled.")
    else:
        emit("      No anchor changes class between the two columns. Verdicts are")
        emit("      robust to this choice.")
    emit()
    emit(f"   b) Where the cutoffs sit (0.25 / 0.05 are conventions, not physics).")
    if borderline:
        emit(f"      {len(borderline)} anchor(s) sit within {BORDERLINE_MARGIN} of a")
        emit(f"      cutoff: {', '.join(f'{ANCHORS[p]}{p}' for p in borderline)}")
        emit("      Read these as 'somewhere near the boundary', not as the label.")
    else:
        emit(f"      No anchor sits within {BORDERLINE_MARGIN} of a cutoff.")
    emit()

    # ---- Glycan caveat: a real occlusion risk the sequence analysis cannot see ----
    emit("4. CAVEAT NOT IN THE ORIGINAL PLAN: glycan occlusion")
    emit()
    emit("   EGFR is a glycoprotein -- sugar chains are attached to it at specific")
    emit("   points. Those chains are large, and the receptor-only file used above")
    emit("   contains protein atoms ONLY, because step 02 stripped everything that")
    emit("   was not a standard amino acid. So the surface measured above is the")
    emit("   BARE protein. If a sugar chain sits over our epitope, the real")
    emit("   accessible surface is smaller than computed here.")
    emit()
    emit("   Checking the original complex for sugar atoms near the epitope:")
    complex_struct = parser.get_structure("complex", str(COMPLEX_PDB))
    sugars = [r for r in complex_struct[0].get_residues()
              if r.get_resname().strip() in GLYCANS]
    emit(f"   Sugar residues present in 6ARU: {len(sugars)}")
    if sugars:
        epi_atoms = []
        for pos in range(EPI_START, EPI_END + 1):
            if pos in records:
                pdb_num = records[pos]["pdb_num"]
                try:
                    res = complex_struct[0]["A"][pdb_num]
                except KeyError:
                    continue
                for atom in res:
                    if atom.element != "H":
                        epi_atoms.append((pos, atom))
        close = {}
        for sug in sugars:
            for satom in sug:
                for pos, eatom in epi_atoms:
                    d = satom - eatom
                    if d < 5.0 and (pos not in close or d < close[pos][0]):
                        close[pos] = (d, sug.get_resname().strip(),
                                      sug.get_parent().id)
        if close:
            emit(f"   Epitope residues within 5 A of a sugar atom: {len(close)}")
            for pos in sorted(close):
                d, name, ch = close[pos]
                mark = " <-- ANCHOR" if pos in ANCHORS else ""
                emit(f"     {records[pos]['aa']}{pos}: {d:.2f} A from "
                     f"{name} (chain {ch}){mark}")
            hit_anchors = [p for p in close if p in ANCHORS]
            if hit_anchors:
                emit()
                emit("   CONSEQUENCE: at least one anchor has a sugar chain nearby.")
                emit("   Its usable surface is smaller than the bare-protein number")
                emit("   above, and glycans are flexible so the real extent is not")
                emit("   fixed by this structure. Treat those anchors' RSA as an")
                emit("   UPPER BOUND, not a measurement.")
            else:
                emit()
                emit("   No anchor is within 5 A of a sugar. Glycan occlusion is")
                emit("   not a concern for the anchors specifically.")
        else:
            emit("   No epitope residue is within 5 A of any sugar atom in this")
            emit("   structure. Glycan occlusion is not a concern here.")
        emit()
        emit("   LIMIT: a crystal structure resolves only the innermost, most")
        emit("   ordered sugars. Real glycan chains extend considerably further")
        emit("   and are mobile. Absence of a modelled sugar is weak evidence of")
        emit("   absence. Flagged as a residual risk, not resolved.")
    emit()

    # ---- What this means ----
    emit("5. WHAT THIS MEANS FOR THE DESIGN")
    emit()
    if n_usable >= 4:
        emit(f"   {n_usable} of 8 anchors are reachable. The pairing rule needs")
        emit("   three or four pairs stacked, because the pH switch is partial")
        emit("   rather than binary, so this is enough material to work with --")
        emit("   PROVIDED they cluster within reach of a single binder, which is")
        emit("   step 05's question and is not answered here.")
    elif n_usable == 3:
        emit("   Exactly 3 anchors reachable. That is the bare minimum for a")
        emit("   stacked switch and leaves no redundancy: if step 05 finds any of")
        emit("   the three sits apart from the others, the epitope fails.")
    else:
        emit(f"   Only {n_usable} anchors reachable. Fewer than three means the")
        emit("   stacked switch the design depends on cannot be built here.")
        emit("   Consider the fallback blocks 394-411 and 331-347.")
    emit()
    emit("   NOTE ON WHAT EXPOSURE DOES AND DOES NOT TELL US. An exposed anchor")
    emit("   is reachable. It is not therefore a good contact point: it still has")
    emit("   to point in a compatible direction, sit in a pocket a designed")
    emit("   backbone can present a partner to, and cluster with the others.")
    emit("   Exposure is a filter that removes impossible options, not evidence")
    emit("   that the remaining ones work.")
    emit()

    # ---- CSVs ----
    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "03-epitope-rsa.csv").open("w") as fh:
        fh.write("uniprot_pos,pdb_resnum,aa,sasa_a2,rsa_theoretical,"
                 "rsa_empirical,class_theoretical,class_empirical,is_anchor,"
                 "borderline,norm_sensitive\n")
        for pos in range(EPI_START, EPI_END + 1):
            r = records.get(pos)
            if r is None:
                fh.write(f"{pos},,,,,,unresolved,unresolved,"
                         f"{pos in ANCHORS},,\n")
                continue
            fh.write(f"{pos},{r['pdb_num']},{r['aa']},{r['sasa']:.2f},"
                     f"{r['rsa_t']:.4f},{r['rsa_e']:.4f},{r['cls_t']},"
                     f"{r['cls_e']},{pos in ANCHORS},"
                     f"{is_borderline(r['rsa_t'])},"
                     f"{r['cls_t'] != r['cls_e']}\n")
    with (DERIVED / "03-anchor-verdicts.csv").open("w") as fh:
        fh.write("uniprot_pos,aa,role,sasa_a2,rsa_theoretical,rsa_empirical,"
                 "verdict,usable\n")
        for r in anchor_rows:
            fh.write(f"{r[0]},{r[1]},\"{r[2]}\",{r[3]},{r[4]},{r[5]},"
                     f"{r[6]},{r[7]}\n")
    emit("Wrote data/derived/03-epitope-rsa.csv")
    emit("Wrote data/derived/03-anchor-verdicts.csv")

    emit()
    emit("=" * 72)
    emit(f"RESULT: {n_usable} of 8 anchors usable "
         f"({', '.join(f'{ANCHORS[p]}{p}' for p in usable) or 'none'}).")
    emit("Clustering is step 05's question. This step only removed the buried.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "03-solvent-accessibility.md").write_text(
        "# Solvent accessibility of the candidate epitope\n\n"
        "Computed output of `analysis/03_solvent_accessibility.py`. Do not hand-edit.\n\n"
        "Normalisation maxima are from Table 1 of Tien MZ, Meyer AG, Sydykova DK,\n"
        "Spielman SJ, Wilke CO (2013), \"Maximum Allowed Solvent Accessibilites of\n"
        "Residues in Proteins\", PLOS ONE 8(11): e80635,\n"
        "doi:10.1371/journal.pone.0080635 — transcribed from the journal's\n"
        "manuscript XML on 30 September 2026.\n\n"
        "The classification cutoffs (RSA >= 0.25 exposed, <= 0.05 buried) are\n"
        "**conventional, not physical constants**. Nothing changes in the molecule\n"
        "at 0.25. Residues near a cutoff are flagged BORDERLINE.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
