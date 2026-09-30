#!/usr/bin/env python3
"""
03_solvent_accessibility.py — which anchors are actually on the surface?

WHAT THIS MEASURES AND WHY IT MATTERS
-------------------------------------
Solvent accessibility is how much of a residue is exposed to the surrounding
water. It is computed by rolling a water-sized probe ball over the protein's
surface and measuring the area that probe can touch for each residue -- the
solvent-accessible surface area, or SASA, reported in square angstroms. A residue
with a large accessible area sticks out, and a binder can reach it. A residue
with an area near zero is buried inside the protein's core: it exists and it is
conserved, and it is still useless to us, because nothing can touch it.

This step decides how much material we actually have. Our eight anchors were
chosen by sequence analysis, which cannot see direction. Domain III is a solenoid
-- a spiral staircase -- so some of those eight point outward and some point into
the core. Only the outward ones are real candidates.

WHY THE RECEPTOR ALONE RATHER THAN THE COMPLEX
------------------------------------------------
This runs on the receptor-only structure written by step 02, with the cetuximab
Fab removed. A Fab is the gripping arm of an antibody, cut free of the rest. The
alternative, running on the complex as deposited, was rejected because cetuximab
sits directly on the surface we care about: the complex would measure our epitope
as buried, when it is only covered by an antibody that will be nowhere near our
assay, and we would drop a workable epitope on the strength of that.

RAW AREA CANNOT BE COMPARED BETWEEN RESIDUE TYPES
-------------------------------------------------
Amino acids (aa) differ enormously in size: tryptophan is far larger than
glycine. So 50 square angstroms of exposed tryptophan is mostly buried, while 50
square angstroms of exposed glycine is wide open. The raw area therefore says
nothing on its own about how exposed a residue is for its own size.

The fix is relative solvent accessibility (RSA): the measured area divided by the
largest area that amino acid type could possibly expose. That puts every residue
type on one scale, running roughly 0 to 1, and RSA is the number to read in the
tables below.

REFERENCE MAXIMA — SOURCE AND STATUS
------------------------------------
The largest-possible areas used for that division are from Table 1 of:

  Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013)
  "Maximum Allowed Solvent Accessibilites of Residues in Proteins."
  PLOS ONE 8(11): e80635. doi:10.1371/journal.pone.0080635

Both published columns are transcribed below, fetched from the journal's
manuscript XML on 30 September 2026 rather than recalled from memory. The
theoretical column gives the headline numbers. The empirical column is then used
as a sensitivity check: if a residue's classification changes depending on which
column you divide by, that residue's verdict is not robust, and the script says
so rather than presenting one of the two columns as the answer.

Several such tables exist in the literature and they disagree by up to 20%, so
which one was used has to be stated for any RSA number to be reproducible.

WHERE THE CLASSIFICATION CUTOFFS COME FROM
------------------------------------------
  RSA >= 0.25  exposed
  RSA <= 0.05  buried
  in between   partially exposed

These are conventional thresholds that the field has settled on for convenience,
and nothing physical changes in the molecule at 0.25. A residue at 0.24 and one
at 0.26 are essentially the same thing. This script therefore flags any residue
within 0.05 of a cutoff as BORDERLINE, so the label is not read as more precise
than it is.

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

# Numbering offset measured by step 02: a position in UniProt, the public protein
# sequence archive, equals the residue number in the Protein Data Bank (PDB) file --
# the public archive of measured 3D structures -- plus 24. Read back from step 02's
# output here rather than typed in again, so the two scripts cannot drift apart.
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

# Sugar residue names to look for when checking whether an attached sugar chain
# covers the epitope. NAG is N-acetylglucosamine, the first sugar of such a chain;
# the others are further sugars that appear in these chains.
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
    """Read the offset from step 02's computed output instead of hardcoding it here.

    Hardcoding would let this script and step 02 drift apart without any error.
    """
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
    emit("Input:  data/structures/6aru_receptor_only.pdb  (Fab removed)")
    emit("Method: Shrake-Rupley probe rolling, Bio.PDB.SASA.ShrakeRupley")
    emit("Norm:   Tien et al. 2013, PLOS ONE 8(11):e80635, Table 1")
    emit("        theoretical column for headline values,")
    emit("        empirical column as a sensitivity check")
    emit(f"Numbering: UniProt = PDB + 24 (measured by step 02, not assumed)")
    emit()
    emit("The Fab -- the gripping arm of the cetuximab antibody -- was removed on")
    emit("purpose. It sits on the surface we are measuring, so running this on the")
    emit("whole complex would report our epitope as buried when it is only covered")
    emit("by an antibody that will not be present in our assay.")
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
    emit("   'aa' is the amino acid at that position, 'SASA A^2' its accessible area")
    emit("   in square angstroms, and RSA that area corrected for residue size, which")
    emit("   is the number to read. 'ANCHOR' marks our eight candidates. 'BORDERLINE'")
    emit("   means the value sits within 0.05 of a classification cutoff, so the label")
    emit("   on that row is less precise than it looks.")
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
        emit("   No coordinates means no answer for those positions: unknown, and not")
        emit("   the same thing as buried.")
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
    emit("2. Verdict on each of the eight anchors")
    emit()
    emit("   This is what decides how much material the design has left to work with.")
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
        # Usable means exposed or partially exposed. A buried anchor cannot be
        # reached by a binder, so it drops out of the design here.
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
    emit("   Two ways a verdict could be an artefact of a choice we made arbitrarily:")
    emit()
    emit(f"   a) Which column of maximum areas we divide by.")
    if norm_sensitive:
        emit(f"      {len(norm_sensitive)} anchor(s) change class between the")
        emit(f"      theoretical and empirical columns: "
             f"{', '.join(f'{ANCHORS[p]}{p}' for p in norm_sensitive)}")
        emit("      Their verdicts are not robust, so treat those anchors as")
        emit("      ambiguous rather than settled.")
    else:
        emit("      No anchor changes class between the two columns. Verdicts are")
        emit("      robust to this choice.")
    emit()
    emit(f"   b) Where the cutoffs sit. 0.25 and 0.05 are conventions the field uses.")
    if borderline:
        emit(f"      {len(borderline)} anchor(s) sit within {BORDERLINE_MARGIN} of a")
        emit(f"      cutoff: {', '.join(f'{ANCHORS[p]}{p}' for p in borderline)}")
        emit("      Read these as sitting near the boundary rather than as whatever")
        emit("      label the row happens to carry.")
    else:
        emit(f"      No anchor sits within {BORDERLINE_MARGIN} of a cutoff.")
    emit()

    # ---- Glycans: a way the surface can be blocked that sequence analysis misses ----
    emit("4. A caveat not in the original plan: sugar chains in the way")
    emit()
    emit("   EGFR is a glycoprotein: sugar chains, called glycans, are attached to it")
    emit("   at specific points, and those chains are large. The receptor-only file")
    emit("   used above holds protein atoms only, because step 02 dropped everything")
    emit("   that was not a standard amino acid, so the areas above describe the bare")
    emit("   protein. If a sugar chain sits over our epitope, the surface a binder")
    emit("   could really reach is smaller than the numbers above.")
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
            emit()
            emit("   What this 5 A test settles, and what it does not. It asks whether")
            emit("   an anchor touches a sugar atom that is actually present in the")
            emit("   file. A structure shows only the first few sugars of a chain that")
            emit("   continues past them, so a 'no' here does not mean the full chain")
            emit("   cannot reach that far. Step 05 asks the wider question, measuring")
            emit("   each anchor's distance to the attachment point N444 -- the")
            emit("   asparagine the chain hangs off -- in 15 A and 25 A bands, and it")
            emit("   flags three anchors on that basis. The two tests answer different")
            emit("   questions, so the results do not contradict each other.")
            hit_anchors = [p for p in close if p in ANCHORS]
            if hit_anchors:
                emit()
                emit("   What follows from this: at least one anchor has a sugar chain")
                emit("   nearby, so its usable surface is smaller than the bare-protein")
                emit("   number above. Glycans are flexible, so this structure does not")
                emit("   fix how far the chain actually extends. Treat those anchors'")
                emit("   RSA as an upper bound on what a binder could reach.")
            else:
                emit()
                emit("   No anchor is within 5 A of a sugar. Glycan occlusion is")
                emit("   not a concern for the anchors specifically.")
        else:
            emit("   No epitope residue is within 5 A of any sugar atom in this")
            emit("   structure. Glycan occlusion is not a concern here.")
        emit()
        emit("   One limit on all of the above: a crystal structure resolves only the")
        emit("   innermost, most ordered sugars. Real glycan chains extend considerably")
        emit("   further and they move about, so a sugar missing from the file is weak")
        emit("   evidence that nothing is there. This stays on the list as a residual")
        emit("   risk rather than something we have settled.")
    emit()

    # ---- What this means ----
    emit("5. What this means for the design")
    emit()
    if n_usable >= 4:
        emit(f"   {n_usable} of 8 anchors are reachable. The pairing rule needs three")
        emit("   or four pairs stacked, because the pH switch is partial rather than")
        emit("   all-or-nothing, so this is enough material to work with, provided the")
        emit("   anchors cluster within reach of a single binder. That is step 05's")
        emit("   question and is not answered here.")
    elif n_usable == 3:
        emit("   Exactly 3 anchors reachable. That is the bare minimum for a")
        emit("   stacked switch and leaves no redundancy: if step 05 finds any of")
        emit("   the three sits apart from the others, the epitope fails.")
    else:
        emit(f"   Only {n_usable} anchors reachable. Fewer than three means the")
        emit("   stacked switch the design depends on cannot be built here.")
        emit("   Consider the fallback blocks 394-411 and 331-347.")
    emit()
    emit("   What exposure does and does not tell us: an exposed anchor is reachable,")
    emit("   and that is all it says. To be a good contact point it still has to point")
    emit("   in a compatible direction, sit in a pocket a designed backbone can present")
    emit("   a partner to, and cluster with the other anchors. Exposure removes the")
    emit("   impossible options; it is no evidence that the ones left over will work.")
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
    emit("Whether they cluster is step 05's question. This step removed the buried.")
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
        "conventions the field has settled on, rather than physical constants.\n"
        "Nothing changes in the molecule at 0.25, so residues near a cutoff are\n"
        "flagged BORDERLINE.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
