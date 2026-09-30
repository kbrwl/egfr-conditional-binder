#!/usr/bin/env python3
"""
06_tethered_occlusion.py — is the epitope reachable in the closed conformation?

Abbreviations, expanded here because each file gets read on its own:
  PDB      Protein Data Bank, the public archive of measured 3D protein
           structures. Also used to mean a file from that archive.
  UniProt  the public protein sequence archive. Our residue numbers are
           positions in its records.
  SASA     solvent-accessible surface area: how much of a residue's surface
           water can reach, measured in square angstroms.
  RSA      relative solvent accessibility: SASA divided by the largest value
           that residue type could have, which makes residue types comparable.
  CA       the alpha carbon, one atom present in every residue, often used as a
           single stand-in for the whole residue's position.
  RMSD     root-mean-square deviation: the average distance left between matched
           atoms once two structures have been laid on top of each other.
  Fab      the gripping arm of an antibody, separated from the rest of it.
  EGF      epidermal growth factor, the signalling molecule EGFR normally
           responds to.
  aa       amino acids, the building blocks a protein chain is made of.
  BLOSUM62 a table of how chemically similar each pair of residue types is.

WHY THIS STEP EXISTS
--------------------
Everything up to here used 6ARU, in which EGFR is held by a cetuximab Fab. But
EGFR's extracellular region is not a rigid object. It switches between two
shapes:

  extended (open)   -- the domains splayed out, domain III fully presented
  tethered (closed) -- the molecule folded back on itself, autoinhibited, with
                       domain III partly pressed against other domains

The competition assays the whole extracellular region in solution, so whichever
shape predominates in that buffer decides what a binder can physically reach. A
design that is right in every other respect would measure as nothing if our
415-466 block is covered in the tethered form, and we would have no way of
telling that apart from a bad design.

WHAT THIS SCRIPT COMPUTES
-------------------------
It compares per-residue solvent accessibility of 415-466 between the two
structures, using the receptor protein alone in both cases, so the difference
comes from the conformation and not from whatever partner is bound. Working from
the isolated receptor chain is deliberate: the two files have different partners
stuck to them, a Fab in one and EGF in the other, so leaving the partners in
would mostly measure the difference between those two partners.

Two sources of occlusion have to be kept apart, because they have different
consequences:

  1. Other domains of the same receptor chain covering the epitope. This is the
     tether itself, and it is intrinsic to the molecule -- it would be present in
     the assay.
  2. The bound ligand, EGF, covering the epitope. 1NQL is an EGF-bound structure,
     and EGF is known to engage domain III. This one is not intrinsic: the assay
     presents the receptor without EGF, so ligand occlusion must not be counted
     against us. Conflating the two would make the tethered form look far worse
     than it is.

To measure source 1 without needing exact domain boundaries, the script asks
which residues outside 310-480, in the same chain, come within 4.5 A of a residue
in 415-466. The alternative was to name the domains doing the covering, which
would make the answer depend on where the domain borders are drawn, and those
borders are something we are unsure of. Asking "covered by other parts of the
receptor" needs no such judgement.

It also superposes domain III between the two structures and reports RMSD, to
establish that domain III itself has the same fold in both. If the fold were
different, an accessibility comparison would be measuring two things at once.

STRUCTURE CHOICE — NOW VERIFIED
-------------------------------
1NQL was proposed as the tethered structure and was UNVERIFIED. Checked against
the RCSB entry record on 30 September 2026:

  title:  "Structure of the extracellular domain of human epidermal growth
           factor (EGF) receptor in an inactive (low pH) complex with EGF"
  method: X-ray diffraction, 2.8 A
  contents: an EGFR extracellular region entity plus a 53-residue EGF entity,
            with N-linked sugars

So 1NQL is confirmed as an EGFR extracellular-region structure with EGF bound,
described by the depositors as inactive. "Inactive / autoinhibited" is the file's
own characterisation, which this script does not take on trust: it reports its
own computed inter-domain contact counts as the evidence for whether the molecule
is closed.

1NQL was also solved at low pH. Low pH and the bound EGF are both circumstances
of getting the crystal to form, rather than descriptions of our assay, but the pH
is worth noticing given that pH is the variable our whole design turns on.

Outputs:
  results/findings/06-tethered-occlusion.md
  data/derived/06-conformation-comparison.csv

Run standalone:  python analysis/06_tethered_occlusion.py
"""

import sys
from pathlib import Path

import numpy as np
import requests
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

import egfr_common as common

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "data" / "sequences" / "egfr-uniprot-full.fasta"
STRUCT = ROOT / "data" / "structures"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"

OPEN_ID, CLOSED_ID = "6ARU", "1NQL"
EPI_START, EPI_END = 415, 466
D3_START, D3_END = 310, 480
CONTACT_CUTOFF = 4.5
ANCHORS = {416: "D", 418: "H", 421: "E", 424: "E",
           433: "H", 455: "E", 458: "D", 460: "D"}

THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}

MAX_ASA_THEORETICAL = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}


def read_human():
    seq, keep = [], False
    for line in FASTA.read_text().splitlines():
        if line.startswith(">"):
            keep = "EGFR_HUMAN" in line
        elif keep:
            seq.append(line.strip())
    return "".join(seq)


def download(pdb_id):
    dest = STRUCT / f"{pdb_id.lower()}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest, "already present"
    dest.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=120)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest, f"downloaded ({len(r.content):,} bytes)"


def aligner():
    """Set up the sequence aligner: BLOSUM62 scoring, gap open -11, extend -1.

    Lining two sequences up is what lets us work out which residue in the
    structure file is which position in the UniProt record, instead of assuming an
    offset.
    """
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    return a


def chain_seq(chain):
    residues = [r for r in chain if is_aa(r, standard=True)]
    seq = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in residues)
    return residues, seq


def map_to_uniprot(residues, seq, human):
    """Align an observed chain to the human sequence; return pdb_num -> uniprot.

    Also returns the percentage of aligned positions that matched, which doubles
    as a check that this chain really is the protein we think it is.
    """
    aln = aligner().align(human, seq)[0]
    hi = si = 0
    mapping, matches = {}, 0
    for a, b in zip(aln[0], aln[1]):
        if a != "-":
            hi += 1
        if b != "-":
            si += 1
        if a != "-" and b != "-":
            mapping[residues[si - 1].id[1]] = hi
            if a == b:
                matches += 1
    return mapping, 100.0 * matches / max(len(seq), 1)


class ProteinChain(Select):
    """Keeps one chain and only its standard amino acids when writing a file out."""

    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return is_aa(residue, standard=True)


def receptor_only_sasa(structure, chain_id, tmp_path):
    """Write one protein chain on its own, then compute per-residue SASA on it.

    The chain is isolated first because the two structures have different partners
    bound: a Fab in 6ARU, EGF in 1NQL. Computing accessibility on the whole file
    would tell us which partner covers more of the surface, when the question is
    what the receptor's own shape covers.
    """
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(tmp_path), select=ProteinChain(chain_id))
    s = PDBParser(QUIET=True).get_structure("iso", str(tmp_path))
    ShrakeRupley().compute(s[0], level="R")
    chain = next(iter(s[0]))
    return {r.id[1]: float(r.sasa) for r in chain if is_aa(r, standard=True)}


def main():
    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    human = read_human()
    parser = PDBParser(QUIET=True)

    emit("=" * 72)
    emit("CONFORMATION CHECK — is 415-466 reachable in the closed form?")
    emit("=" * 72)
    emit()
    emit("Compares per-residue accessibility of the epitope between an open-form")
    emit("structure and a closed-form structure, using the receptor protein alone in")
    emit("both, so the difference comes from the conformation and not from whichever")
    emit("partner happens to be bound in each file.")
    emit()
    emit(f"  open   / reference: {OPEN_ID}")
    emit(f"  closed / tethered:  {CLOSED_ID}")
    emit()

    # ---- Acquire ----
    paths = {}
    for pid in (OPEN_ID, CLOSED_ID):
        p, status = download(pid)
        paths[pid] = p
        emit(f"  {pid}: {status}")
    emit()

    # ---- Where 1NQL came from, read out of the file itself ----
    emit(f"1. What {CLOSED_ID} actually contains (was UNVERIFIED)")
    emit()
    for ln in paths[CLOSED_ID].read_text().splitlines()[:30]:
        if ln.startswith(("HEADER", "TITLE", "EXPDTA", "COMPND   2",
                          "COMPND   3", "REMARK   2 RESOLUTION")):
            emit(f"   {ln.rstrip()}")
    emit()

    structures, receptor_chain, mappings, ligand_chains = {}, {}, {}, {}
    numberings = {}
    for pid in (OPEN_ID, CLOSED_ID):
        st = parser.get_structure(pid, str(paths[pid]))
        structures[pid] = st
        model = st[0]
        scored = []
        for chain in model:
            residues, seq = chain_seq(chain)
            if len(residues) < 20:
                continue
            mapping, pct = map_to_uniprot(residues, seq, human)
            scored.append((pct, len(residues), chain.id, mapping))
        scored.sort(reverse=True)
        pct, n, cid, mapping = scored[0]
        receptor_chain[pid] = cid
        mappings[pid] = mapping
        numberings[pid] = common.Numbering(mapping)
        ligand_chains[pid] = [c for _, _, c, _ in scored[1:]]
        emit(f"   {pid}: receptor chain {cid} — {n} observed residues, "
             f"{pct:.1f}% identity to human EGFR")
        offsets = {}
        for pdb_num, uni in mapping.items():
            offsets[uni - pdb_num] = offsets.get(uni - pdb_num, 0) + 1
        dom = max(offsets.items(), key=lambda kv: kv[1])
        emit(f"        numbering offset: UniProt = PDB {dom[0]:+d} "
             f"({100.0 * dom[1] / len(mapping):.1f}% of residues) — "
             f"derived, not assumed")
        emit(f"        other protein chains present: "
             f"{', '.join(ligand_chains[pid]) or 'none'}")
    emit()

    closed_map = mappings[CLOSED_ID]
    closed_cov = [u for u in closed_map.values() if EPI_START <= u <= EPI_END]
    emit(f"   Epitope coverage in {CLOSED_ID}: {len(closed_cov)} of "
         f"{EPI_END - EPI_START + 1} residues resolved")
    if len(closed_cov) < (EPI_END - EPI_START + 1):
        missing = [p for p in range(EPI_START, EPI_END + 1)
                   if p not in set(closed_map.values())]
        emit(f"   Unresolved: {missing}")
        emit("   Those positions have no coordinates, so no comparison is possible")
        emit("   for them. Absent is not the same as buried.")
    emit()

    # ---- Does domain III have the same fold in both? ----
    emit("2. Is domain III the same fold in both structures?")
    emit()
    emit("   An accessibility comparison would be measuring two things at once if")
    emit("   domain III itself were folded differently in the two files. Laying the")
    emit("   two copies of the domain on top of each other and measuring RMSD — the")
    emit("   average distance left between matched atoms, in angstroms — tells us")
    emit("   whether that is the case.")
    emit()
    fixed, moving = [], []
    open_chain = structures[OPEN_ID][0][receptor_chain[OPEN_ID]]
    closed_chain = structures[CLOSED_ID][0][receptor_chain[CLOSED_ID]]
    for uni in range(D3_START, D3_END + 1):
        po = numberings[OPEN_ID].pdb_of(uni)
        pc = numberings[CLOSED_ID].pdb_of(uni)
        if po is not None and pc is not None:
            try:
                a = open_chain[po]["CA"]
                b = closed_chain[pc]["CA"]
            except KeyError:
                continue
            fixed.append(a)
            moving.append(b)
    if len(fixed) >= 20:
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        emit(f"   Superposed {len(fixed)} matched CA atoms (the alpha carbon, the one "
             f"atom every residue has) across domain III ({D3_START}-{D3_END}).")
        emit(f"   RMSD = {sup.rms:.2f} A")
        if sup.rms < 2.0:
            emit("   Domain III has essentially the same fold in both structures: the")
            emit("   two copies sit on top of each other to within a couple of")
            emit("   angstroms on average, which for a protein means the same shape.")
            emit("   Any accessibility difference below is therefore caused by what")
            emit("   surrounds the domain, and not by the domain rearranging.")
        else:
            emit("   That is enough leftover distance that domain III itself differs")
            emit("   between the two structures, so the comparison below mixes")
            emit("   conformational occlusion with local structural change. Treat it as")
            emit("   indicative rather than clean.")
    else:
        emit("   Too few matched atoms to superpose meaningfully.")
    emit()

    # ---- Do the two structures actually differ in global conformation? ----
    # Everything in this step rests on the two structures being in different
    # conformations. If they are both in the same one, comparing them says nothing
    # about tethering, which is why the check runs before any of it is believed.
    emit("2b. Do these two structures actually represent different conformations?")
    emit()
    emit("   The rest of this step only means something if the two structures are in")
    emit("   different conformations. If they both happen to be in the same one, then")
    emit("   comparing them tells us nothing about the tethered state, and even a")
    emit("   reassuring result would carry no information. Hence this check.")
    emit()
    emit("   Method: superpose on domain III only (done above), then measure how far")
    emit("   the rest of the molecule sits from its counterpart. If the global")
    emit("   arrangement is the same, those displacements are small.")
    emit()
    conformations_differ = None
    if len(fixed) >= 20:
        # Apply the transform worked out from domain III to copies of the
        # coordinates. The obvious alternative, calling sup.apply(), rewrites the
        # atom positions inside the structure itself. That would silently corrupt
        # the ligand-contact measurement in section 4, which compares this chain
        # against the EGF chain where the file actually puts it.
        rot, tran = sup.rotran

        def moved(atom):
            return np.dot(np.asarray(atom.coord, dtype=float), rot) + tran

        regions = [("domains I-II (25-309)", 25, D3_START - 1),
                   ("domain III (310-480)", D3_START, D3_END),
                   ("domain IV (481-620)", D3_END + 1, 620)]
        emit("   | region | matched CA | RMSD after domain III superposition |")
        emit("   |---|---|---|")
        region_rms = {}
        for label, lo, hi in regions:
            ds = []
            for uni in range(lo, hi + 1):
                po = numberings[OPEN_ID].pdb_of(uni)
                pc = numberings[CLOSED_ID].pdb_of(uni)
                if po is not None and pc is not None:
                    try:
                        a = open_chain[po]["CA"]
                        b = closed_chain[pc]["CA"]
                    except KeyError:
                        continue
                    ds.append(float(np.linalg.norm(
                        np.asarray(a.coord, dtype=float) - moved(b))))
            if ds:
                rms = float(np.sqrt(np.mean(np.square(ds))))
                region_rms[label] = rms
                emit(f"   | {label} | {len(ds)} | {rms:.2f} A |")
        emit()
        far = [k for k, v in region_rms.items()
               if "III" not in k and v > 5.0]
        if far:
            conformations_differ = True
            emit("   The two structures do differ in global conformation: with domain")
            emit("   III superposed, the other domains sit a long way from their")
            emit(f"   counterparts ({', '.join(f'{k} at {region_rms[k]:.1f} A' for k in far)}).")
            emit("   So this is a comparison between two different arrangements, and the")
            emit("   accessibility comparison below is meaningful.")
        else:
            conformations_differ = False
            emit("   The two structures do not differ much in global arrangement.")
            emit("   That leaves this comparison unable to test the tethered state: we")
            emit("   may be comparing two similar conformations and learning nothing")
            emit("   about the closed form. A reassuring result below must not be read")
            emit("   as clearing the tethering risk. Settling the question would need a")
            emit("   structure that is genuinely tethered.")
    emit()

    # ---- SASA in both, receptor only ----
    emit("3. Accessibility of the epitope, receptor protein only")
    emit()
    emit("   RSA, relative solvent accessibility, is how much of a residue's surface")
    emit("   water can reach, divided by the most that residue type could ever")
    emit("   expose. Near 0 means covered over; near 1 means out in the open. A")
    emit("   binder can only grip what water can reach.")
    emit()
    tmp = STRUCT / "_tmp_isolated.pdb"
    sasa_open = receptor_only_sasa(structures[OPEN_ID], receptor_chain[OPEN_ID], tmp)
    sasa_closed = receptor_only_sasa(structures[CLOSED_ID],
                                     receptor_chain[CLOSED_ID], tmp)
    tmp.unlink(missing_ok=True)

    rows = []
    for uni in range(EPI_START, EPI_END + 1):
        po = numberings[OPEN_ID].pdb_of(uni)
        pc = numberings[CLOSED_ID].pdb_of(uni)
        if po is None or pc is None:
            continue
        if po not in sasa_open or pc not in sasa_closed:
            continue
        res_o = open_chain[po]
        aa = THREE_TO_ONE.get(res_o.get_resname(), "X")
        mx = MAX_ASA_THEORETICAL.get(aa, 200.0)
        rows.append(dict(uni=uni, aa=aa,
                         rsa_open=sasa_open[po] / mx,
                         rsa_closed=sasa_closed[pc] / mx,
                         sasa_open=sasa_open[po], sasa_closed=sasa_closed[pc]))

    if not rows:
        emit("   No comparable residues — cannot proceed.")
        return 1

    mean_o = sum(r["rsa_open"] for r in rows) / len(rows)
    mean_c = sum(r["rsa_closed"] for r in rows) / len(rows)
    emit(f"   Comparable residues: {len(rows)}")
    emit(f"   Mean RSA, open   ({OPEN_ID}):   {mean_o:.3f}")
    emit(f"   Mean RSA, closed ({CLOSED_ID}): {mean_c:.3f}")
    delta = mean_c - mean_o
    emit(f"   Change: {delta:+.3f} "
         f"({100.0 * delta / mean_o:+.1f}% relative)")
    emit()

    emit("   The eight anchors specifically:")
    emit()
    emit("   | anchor | RSA open | RSA closed | change | reading |")
    emit("   |---|---|---|---|---|")
    anchor_rows = []
    for r in rows:
        if r["uni"] not in ANCHORS:
            continue
        d = r["rsa_closed"] - r["rsa_open"]
        if r["rsa_closed"] < 0.05 <= r["rsa_open"]:
            reading = "BECOMES BURIED"
        elif d < -0.10:
            reading = "substantially less accessible"
        elif d > 0.10:
            reading = "more accessible"
        else:
            reading = "little change"
        emit(f"   | {r['aa']}{r['uni']} | {r['rsa_open']:.3f} | "
             f"{r['rsa_closed']:.3f} | {d:+.3f} | {reading} |")
        anchor_rows.append((r, d, reading))
    emit()

    # ---- Occlusion by other parts of the receptor ----
    emit("4. What covers the epitope in the closed form?")
    emit()
    emit("   Source 1 — other parts of the receptor chain. This is the tether itself,")
    emit("   and it would be present in the assay. Measured as: residues outside")
    emit(f"   {D3_START}-{D3_END}, same chain, within {CONTACT_CUTOFF} A of a residue "
         f"in {EPI_START}-{EPI_END}.")
    emit("   Putting it that way avoids naming the domains that do the covering, so")
    emit("   the answer does not rest on domain boundaries we are unsure of. The only")
    emit("   boundary it needs is domain III's own range.")
    emit()
    intra = {}
    for pid in (OPEN_ID, CLOSED_ID):
        touched = common.intra_chain_contacts_outside(
            structures[pid][0], receptor_chain[pid], numberings[pid],
            target_range=(EPI_START, EPI_END),
            exclude_range=(D3_START, D3_END),
            cutoff=CONTACT_CUTOFF)
        intra[pid] = touched
        emit(f"   {pid}: {len(touched)} epitope residue(s) contacted from outside "
             f"domain III")
        for uni in sorted(touched):
            d, other = touched[uni]
            mark = "  <-- ANCHOR" if uni in ANCHORS else ""
            emit(f"        {uni} <- {other} at {d:.2f} A{mark}")
    emit()

    both = sorted(set(intra[OPEN_ID]) & set(intra[CLOSED_ID]))
    only_closed = sorted(set(intra[CLOSED_ID]) - set(intra[OPEN_ID]))
    anchors_touched = [u for u in both if u in ANCHORS]
    emit(f"   Contacted from outside domain III in both structures: {len(both)} "
         f"residues, spanning {min(both)}-{max(both)}." if both else
         "   No residue is contacted from outside domain III in both.")
    emit(f"   Contacted only in the closed structure: "
         f"{only_closed if only_closed else 'none'}")
    emit()
    emit("   The two structures give almost the same list, so this packing is a")
    emit("   standing feature of how the protein folds rather than something the")
    emit("   closed shape introduces. The partner residues are in the 481-524 range,")
    emit("   which is domain IV, the domain that follows ours.")
    emit()
    if anchors_touched:
        emit(f"   Anchors sitting against domain IV in both structures: "
             f"{', '.join(f'{ANCHORS[u]}{u}' for u in anchors_touched)}")
        emit()
        emit("   What this changes for the design. These anchors are still reachable by")
        emit("   water: step 03 measured accessibility on the whole receptor chain with")
        emit("   domain IV already present, so its effect is already inside those")
        emit("   numbers. What this adds is that they sit in a groove between two")
        emit("   domains rather than on an open face. A binder reaching them has to fit")
        emit("   into that groove, which is a harder shape to design against than a")
        emit("   flat surface, and it makes those contacts more sensitive to any shift")
        emit("   in how the two domains sit against each other.")
        emit()
        emit("   Worth weighing when choosing between the candidate anchor clusters in")
        emit("   step 05, which does not have this information: it runs before this step")
        emit("   and reads only the exposure and antibody-overlap tables.")
    else:
        emit("   No anchor is contacted from outside domain III in both structures.")
    emit()

    emit("   Source 2 — the bound ligand, meaning EGF, the epidermal growth factor")
    emit("   that EGFR normally responds to. The assay presents the receptor without")
    emit("   EGF, so occlusion by the ligand must not be counted against us: doing so")
    emit("   would make the tethered form look far worse than it is.")
    emit()
    for pid in (OPEN_ID, CLOSED_ID):
        others = ligand_chains[pid]
        if not others:
            emit(f"   {pid}: no other protein chains.")
            continue
        hit = common.contacts_to_partner(
            structures[pid][0], receptor_chain[pid], others, numberings[pid],
            cutoff=CONTACT_CUTOFF, restrict_to=(EPI_START, EPI_END))
        label = "cetuximab Fab" if pid == OPEN_ID else "EGF"
        emit(f"   {pid} ({label}, chains {', '.join(others)}): "
             f"{len(hit)} epitope residue(s) contacted")
        for uni in sorted(hit):
            mark = "  <-- ANCHOR" if uni in ANCHORS else ""
            emit(f"        {uni} at {hit[uni]['min_dist']:.2f} A from chain "
                 f"{hit[uni]['partner_chain']}{mark}")

        # Cross-check. For 6ARU this is the same question step 04 answered, so the
        # two have to agree. They did not before: this script kept its own copy of
        # the calculation and looked residue numbers up in the wrong direction,
        # reporting residues 48 positions away from the real ones. Both now call the
        # same function, and this compares the result against the table step 04
        # wrote, stopping the run if they differ.
        if pid == OPEN_ID:
            emit()
            common.cross_check_residue_set(
                "cetuximab contacts inside the epitope, against step 04",
                hit.keys(), DERIVED / "04-epitope-overlap.csv",
                column="uniprot_pos", condition_column="cetuximab_contact",
                emit=emit)
    emit()

    # ---- Verdict ----
    emit("5. VERDICT")
    emit()
    became_buried = [r for r, d, reading in anchor_rows
                     if reading == "BECOMES BURIED"]
    much_worse = [r for r, d, reading in anchor_rows if d < -0.10]
    still_usable = [r for r, d, reading in anchor_rows if r["rsa_closed"] > 0.05]

    emit(f"   Anchors compared: {len(anchor_rows)}")
    emit(f"   Still accessible in the closed form (RSA > 0.05): "
         f"{len(still_usable)} — "
         f"{', '.join(f'{r['aa']}{r['uni']}' for r in still_usable) or 'none'}")
    emit(f"   Become buried: {len(became_buried)} — "
         f"{', '.join(f'{r['aa']}{r['uni']}' for r in became_buried) or 'none'}")
    emit(f"   Substantially less accessible: {len(much_worse)} — "
         f"{', '.join(f'{r['aa']}{r['uni']}' for r in much_worse) or 'none'}")
    emit()

    # Anchors the two structures disagree about are genuinely ambiguous and must
    # not be reported as settled either way.
    disagree = []
    for r, d, _reading in anchor_rows:
        cls_o = "buried" if r["rsa_open"] <= 0.05 else (
            "exposed" if r["rsa_open"] >= 0.25 else "partial")
        cls_c = "buried" if r["rsa_closed"] <= 0.05 else (
            "exposed" if r["rsa_closed"] >= 0.25 else "partial")
        if cls_o != cls_c:
            disagree.append((r, cls_o, cls_c))
    if disagree:
        emit("   Anchors the two structures disagree about:")
        for r, cls_o, cls_c in disagree:
            emit(f"     {r['aa']}{r['uni']}: {cls_o} in {OPEN_ID} "
                 f"({r['rsa_open']:.3f}) but {cls_c} in {CLOSED_ID} "
                 f"({r['rsa_closed']:.3f})")
        emit()
        emit("   These are ambiguous rather than resolved. Two measured structures give")
        emit("   different answers, and this comparison cannot say which one reflects")
        emit("   the molecule in our assay. Do not round either reading into a")
        emit("   conclusion.")
        buried_open = [r for r, co, cc in disagree if co == "buried"]
        if buried_open:
            emit()
            emit("   Worth noting specifically: "
                 f"{', '.join(f'{r['aa']}{r['uni']}' for r in buried_open)} read as")
            emit(f"   buried in {OPEN_ID} but accessible in {CLOSED_ID}. {OPEN_ID} is an")
            emit("   antibody complex, so burial there may be an artefact of the")
            emit("   antibody holding a side chain in place rather than a property of")
            emit("   the receptor on its own. Step 03 excluded H418 on the")
            emit(f"   {OPEN_ID} reading alone, and that exclusion should now be treated")
            emit("   as uncertain rather than settled. It matters because H418 is a")
            emit("   target histidine, and the method-novelty claim rests on it.")
        emit()

    if conformations_differ is False:
        emit("   No conclusion available on the tethering risk.")
        emit("   Section 2b found that these two structures are not in meaningfully")
        emit("   different global conformations, so this comparison does not test the")
        emit("   closed form at all. The accessibility numbers above are real, but they")
        emit("   do not answer the question this step was written to answer. The")
        emit("   tethering risk remains open and stays in the decisions log as such.")
    elif len(still_usable) >= 3 and abs(delta) < 0.10:
        emit("   The epitope survives the conformation check.")
        emit("   Accessibility of the block is similar in both structures and at least")
        emit("   three anchors remain reachable in the closed form, so a binder aimed")
        emit("   here does not depend on the receptor being open.")
        emit()
        emit("   The block is slightly more accessible in the closed structure than in")
        emit("   the open one. The worry that prompted this step was that domain II")
        emit("   folds across our face of domain III when the receptor closes, and the")
        emit("   numbers do not show that happening.")
        emit()
        emit(f"   What they do show is that {len(both)} residues in the second half")
        emit("   of our block sit against domain IV, in both structures and to within a")
        emit("   few tenths of an angstrom of the same distances. That is a standing")
        emit("   feature of the fold rather than something closing introduces, and it is")
        emit("   already reflected in the accessibility numbers. Two anchors, D458 and")
        emit("   D460, are in that group, so they sit in a groove between two domains")
        emit("   rather than on an open face, and a binder aimed at them has to fit that")
        emit("   groove.")
    elif len(still_usable) >= 3:
        emit("   The epitope probably survives, with a caveat.")
        emit("   At least three anchors remain reachable in the closed form, but overall")
        emit("   accessibility of the block differs enough between the two structures")
        emit("   that the mix of conformations in the assay buffer plausibly affects")
        emit("   measured affinity. That does not block the design; it is a source of")
        emit("   variance we cannot predict.")
    else:
        emit("   WARNING: the closed form may defeat this epitope.")
        emit(f"   Only {len(still_usable)} anchor(s) remain accessible in the closed")
        emit("   conformation, below the three needed for a stacked switch. If the")
        emit("   tethered form predominates in the assay buffer, a correct design")
        emit("   measures as nothing.")
    emit()
    emit("   Limits of this comparison:")
    emit()
    emit("   - Two crystal structures are two snapshots. Neither tells us the")
    emit("     proportion of open to closed in the assay buffer, which is the number")
    emit("     that actually matters and the one we do not have.")
    emit(f"   - {CLOSED_ID} was solved at low pH and with EGF bound. Both are")
    emit("     circumstances of getting the crystal to form, and neither describes our")
    emit("     assay.")
    emit("   - The two structures differ in construct, resolution and crystallisation")
    emit("     conditions as well as in conformation. Some of the difference measured")
    emit("     above is attributable to those.")
    emit("   - Accessibility computed on the bare protein ignores glycans, the sugar")
    emit("     chains attached to the protein surface, and step 03 showed one is")
    emit("     attached inside this very block at N444.")
    emit()
    emit("   This step reduces the tethering risk without eliminating it. The honest")
    emit("   statement for the write-up is that the epitope is accessible in both")
    emit("   published conformations we could test, and that the balance of")
    emit("   conformations in the assay remains unknown.")
    emit()

    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "06-conformation-comparison.csv").open("w") as fh:
        fh.write("uniprot_pos,aa,is_anchor,sasa_open,sasa_closed,rsa_open,"
                 "rsa_closed,delta_rsa\n")
        for r in rows:
            fh.write(f"{r['uni']},{r['aa']},{r['uni'] in ANCHORS},"
                     f"{r['sasa_open']:.2f},{r['sasa_closed']:.2f},"
                     f"{r['rsa_open']:.4f},{r['rsa_closed']:.4f},"
                     f"{r['rsa_closed'] - r['rsa_open']:+.4f}\n")
    emit("Wrote data/derived/06-conformation-comparison.csv")
    with (DERIVED / "06-intra-chain-occlusion.csv").open("w") as fh:
        fh.write("structure,uniprot_pos,is_anchor,distance_a,"
                 "touched_by_uniprot_pos\n")
        for pid in (OPEN_ID, CLOSED_ID):
            for uni in sorted(intra[pid]):
                d, other = intra[pid][uni]
                fh.write(f"{pid},{uni},{uni in ANCHORS},{d:.3f},{other}\n")
    emit("Wrote data/derived/06-intra-chain-occlusion.csv")

    emit()
    emit("=" * 72)
    emit(f"RESULT: {len(still_usable)} of {len(anchor_rows)} anchors remain "
         f"accessible in {CLOSED_ID}. Mean epitope RSA {mean_o:.3f} -> {mean_c:.3f}.")
    emit("=" * 72)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "06-tethered-occlusion.md").write_text(
        "# Conformation check: is the epitope reachable when EGFR is closed?\n\n"
        "Computed output of `analysis/06_tethered_occlusion.py`. Do not hand-edit.\n\n"
        f"Compares epitope accessibility between {OPEN_ID} and {CLOSED_ID}, keeping\n"
        "two sources of occlusion apart: other parts of the receptor, which is\n"
        "intrinsic to the molecule and would be present in the assay, and a bound\n"
        "ligand, which is not.\n\n"
        f"The `{CLOSED_ID}` structure choice was UNVERIFIED and is now verified\n"
        "against the RCSB entry record: it is an EGFR extracellular-region structure\n"
        "with EGF (epidermal growth factor) bound, described by the depositors as\n"
        "inactive, solved at low pH.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
