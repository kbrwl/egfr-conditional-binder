"""
Shared code for the EGFR analysis scripts.

Anything computed in more than one place lives here, so two scripts cannot give
two answers to the same question. This file exists because that happened: step 04
and step 06 both worked out which residues the antibody touches, from the same
file with the same settings, and got different answers. One of them had looked a
residue number up in the wrong direction.

Abbreviations, expanded here because each file gets read on its own:
  PDB      Protein Data Bank, the public archive of measured 3D protein
           structures. Also used to mean a file from that archive.
  UniProt  the public protein sequence archive. Our residue numbers are
           positions in its records.
  SASA     solvent-accessible surface area: how much of a residue's surface
           water can reach, measured in square angstroms.
  RSA      relative solvent accessibility: SASA divided by the largest value
           that residue type could have. Comparable between residue types.
  CA       the alpha carbon, one atom present in every residue. Often used as a
           single stand-in for the whole residue's position.
  Fab      the gripping arm of an antibody, separated from the rest of it.
  aa       amino acids, the building blocks a protein chain is made of.
  BLOSUM62 a table of how chemically similar each pair of residue types is.
"""

from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

ROOT = Path(__file__).resolve().parents[1]
SEQ_DIR = ROOT / "data" / "sequences"
STRUCT_DIR = ROOT / "data" / "structures"
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"

UNIPROT_FASTA = SEQ_DIR / "egfr-uniprot-full.fasta"
CHALLENGE_FASTA = SEQ_DIR / "egfr-challenge-constructs.fasta"

# Residue ranges, as positions in the full human UniProt record P00533.
SIGNAL_PEPTIDE_LEN = 24     # residues 1-24 are cut off when the protein is made
ECD_START, ECD_END = 25, 645    # the part of EGFR outside the cell
D3_START, D3_END = 310, 480      # domain III, our working definition
EPI_START, EPI_END = 415, 466    # the candidate epitope, the patch we aim at

# The eight positions we could build a pH switch against, with the amino acid each
# one is. D and E are acidic and always carry a negative charge, so we put a
# histidine opposite them. H is histidine, which changes charge with pH, so we put
# an acidic residue opposite it.
ANCHORS = {416: "D", 418: "H", 421: "E", 424: "E",
           433: "H", 455: "E", 458: "D", 460: "D"}

# The atom each anchor is measured from: the charged tip of its side chain, because
# that is what forms a charge pair. D uses the carboxylate carbon CG, E the
# carboxylate carbon CD one atom further out, and H the centre of the imidazole ring.
# Shared by steps 05 and 14 so the two measure from the same point.
FUNCTIONAL_ATOM = {"D": ["CG"], "E": ["CD"],
                   "H": ["CG", "ND1", "CD2", "CE1", "NE2"]}

# A complex N-linked glycan is a branched chain of sugars rather than a single
# sugar. Only its innermost sugars sit still enough to appear in a crystal
# structure, while the whole assembly is mobile and can sweep 20-30 angstroms from
# the point where it attaches. The two bands are therefore cautious on purpose.
# UNVERIFIED: the 20-30 angstrom reach is from memory and is not measured here.
GLYCAN_NEAR = 15.0       # very likely shadowed some of the time
GLYCAN_PLAUSIBLE = 25.0  # within reach of an extended chain

# The tag on the target in the assay. Both the human and the mouse protein carry a
# C-terminal His tag that the organisers expect to leave on (docs/competition-qa-log.md).
# A binder that grips the tag looks pH-selective and binds anything with a His tag,
# so the design run is told to avoid it as an off-target.
#
# The sequence is the six histidines with a serine-glycine flank on each side. The
# flank is there because BindCraft2 samples a window of at least 10 residues from a
# sequence target by default, and a bare six-residue peptide is shorter than that.
# It also resembles a tag at the end of a chain more than a free peptide does.
# UNVERIFIED: that BindCraft2 accepts this as an off-target has not been run.
HIS_TAG_NAME = "HisTag"
HIS_TAG_SEQUENCE = "GSHHHHHHGS"

# A negative weight tells BindCraft2 to push the binder away from a target, and the
# magnitude sets how hard relative to the binding objective. Half strength is the
# value in BindCraft2's own cross-reactive example. It is a starting value and not
# a tuned one.
HIS_TAG_WEIGHT = -0.5

# BindCraft2's own default ceiling on an off-target's interface confidence: a
# design is rejected if the binder is predicted to bind an off-target at or above
# this (max_detarget_iptm, docs/source/design-guide/02-setting-up-a-design.md).
# It is the tool's default and has not been calibrated for this tag.
DETARGET_IPTM_CEILING = 0.4

# Two residues are counted as touching if any pair of their atoms is this close,
# measured in angstroms (an angstrom is a ten-billionth of a metre). 4.5 is the
# usual choice in the literature for protein-protein contact.
CONTACT_CUTOFF = 4.5

# The site inside our epitope where a sugar chain is attached to the protein.
GLYCOSYLATION_SITE = 444

# Converts three-letter residue names as they appear in structure files (ASP)
# into the single letters used in sequences (D).
THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}

# Largest surface area each residue type can expose, in square angstroms, used to
# turn a raw area into a comparable fraction. From Table 1 of: Tien MZ, Meyer AG,
# Sydykova DK, Spielman SJ, Wilke CO (2013), "Maximum Allowed Solvent
# Accessibilites of Residues in Proteins", PLOS ONE 8(11): e80635. Read off the
# journal's own manuscript file on 30 September 2026 rather than from memory.
#
# The paper gives two columns. We use "theoretical" for our headline numbers and
# "empirical" to check whether a result depends on that choice. Several such
# tables exist and they differ by as much as 20%, so a result is only reproducible
# if it says which one it used.
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
ASA_REFERENCE = ("Tien et al. 2013, PLOS ONE 8(11):e80635, Table 1 "
                 "(doi:10.1371/journal.pone.0080635)")

# Bands used when judging how exposed a residue is. These are the values the field
# has settled on for convenience. Nothing in the molecule changes at 0.25, so a
# residue just either side of a line is really the same thing as one just the
# other side, and results near a line are reported as uncertain.
EXPOSED_CUT, BURIED_CUT = 0.25, 0.05
BORDERLINE_MARGIN = 0.05

# Sugar residue names, so they can be told apart from amino acids in a file.
GLYCAN_NAMES = {"NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "SIA", "BGC", "GLC"}


def read_fasta(path):
    """Read a sequence file into {header line without '>': sequence}."""
    records, header, chunks = {}, None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks)
    return records


def human_sequence():
    """The full human EGFR sequence, all 1210 residues including the part that
    gets cut off. Every residue number in this project is a position in this."""
    for header, seq in read_fasta(UNIPROT_FASTA).items():
        if "EGFR_HUMAN" in header:
            return seq
    raise KeyError("no human record in " + str(UNIPROT_FASTA))


def mouse_sequence():
    for header, seq in read_fasta(UNIPROT_FASTA).items():
        if "EGFR_MOUSE" in header:
            return seq
    raise KeyError("no mouse record in " + str(UNIPROT_FASTA))


def make_aligner(mode="global"):
    """An aligner set up the same way everywhere: BLOSUM62 scoring, gap open -11,
    gap extend -1.

    Lining two sequences up is necessary before comparing them, because a single
    extra residue early on shifts everything after it. BLOSUM62 is a table of how
    chemically similar each pair of residues is, so the alignment prefers swaps
    that are plausible in real proteins. The gap costs are the standard defaults:
    starting a gap is expensive and continuing one is cheap, which matches the
    fact that real insertions tend to happen in one go rather than singly.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    return aligner


class Numbering:
    """Translates between the residue numbers in a structure file and ours.

    Structure files carry their own residue numbers, and no rule says which
    convention they follow. Files of secreted proteins often count from the mature
    protein, whose residue 1 is UniProt residue 25, because that is the molecule
    that was actually crystallised. Nothing in the file records the choice. Both
    structures used in this project turn out to be off by 24.

    The class exists because the two directions used to be passed around as plain
    dictionaries, and step 06 indexed a UniProt-to-structure dictionary with a
    structure number. The two numbering ranges overlap, so that lookup succeeded
    and handed back a number 48 positions away from the right one, raising no
    error. Asking for `uniprot_of(...)` or `pdb_of(...)` by name cannot go wrong
    the same way: give either method the wrong kind of number and it returns
    nothing, instead of a plausible wrong answer that survives into a result.
    """

    def __init__(self, pdb_to_uniprot):
        self._pdb_to_uniprot = dict(pdb_to_uniprot)
        self._uniprot_to_pdb = {u: p for p, u in pdb_to_uniprot.items()}

    def uniprot_of(self, pdb_resnum):
        """Our residue number, given the structure file's number."""
        return self._pdb_to_uniprot.get(pdb_resnum)

    def pdb_of(self, uniprot_pos):
        """The structure file's number, given ours."""
        return self._uniprot_to_pdb.get(uniprot_pos)

    def uniprot_positions(self):
        return sorted(self._pdb_to_uniprot.values())

    def dominant_offset(self):
        """The single value of (our number minus the file's number) that covers
        most of the chain, with the fraction of residues that agree on it. One
        value covering nearly everything means one consistent convention."""
        counts = {}
        for pdb_num, uni in self._pdb_to_uniprot.items():
            counts[uni - pdb_num] = counts.get(uni - pdb_num, 0) + 1
        offset, n = max(counts.items(), key=lambda kv: kv[1])
        return offset, n / max(len(self._pdb_to_uniprot), 1), counts

    def __len__(self):
        return len(self._pdb_to_uniprot)


def numbering_from_alignment(residues, human=None):
    """Work out a structure's numbering by lining its own sequence up against the
    human sequence, rather than assuming an offset.

    `residues` is the list of amino-acid residues read from one chain of the file.
    Returns a Numbering, plus what percentage of the aligned positions matched,
    which doubles as a check that this chain really is the protein we think it is.
    """
    human = human or human_sequence()
    seq = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in residues)
    alignment = make_aligner().align(human, seq)[0]
    human_i = struct_i = matched = 0
    pdb_to_uniprot = {}
    for human_char, struct_char in zip(alignment[0], alignment[1]):
        if human_char != "-":
            human_i += 1
        if struct_char != "-":
            struct_i += 1
        if human_char != "-" and struct_char != "-":
            pdb_to_uniprot[residues[struct_i - 1].id[1]] = human_i
            if human_char == struct_char:
                matched += 1
    pct = 100.0 * matched / max(len(seq), 1)
    return Numbering(pdb_to_uniprot), pct


def load_numbering(csv_path=None):
    """Read a Numbering back from the table step 02 wrote, so later steps use the
    offset that was actually measured rather than one written in by hand."""
    csv_path = Path(csv_path or (DERIVED / "02-numbering-offset.csv"))
    if not csv_path.exists():
        raise SystemExit("Run analysis/02_structure_prep.py first: "
                         f"{csv_path} is missing.")
    pdb_to_uniprot = {}
    for line in csv_path.read_text().splitlines()[1:]:
        uni, pdb, _offset = line.split(",")
        pdb_to_uniprot[int(pdb)] = int(uni)
    return Numbering(pdb_to_uniprot)


def functional_point(res, aa):
    """Pick the atom an anchor is measured from.

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


def protein_residues(chain):
    return [r for r in chain if is_aa(r, standard=True)]


def heavy_atoms(residues):
    """Every atom except hydrogen.

    Hydrogen atoms are too light to show up in most measured structures and are
    simply absent from these files, so "heavy atom" means the atoms we actually
    have positions for.
    """
    return [a for r in residues for a in r if a.element != "H"]


def contact_pairs(model, receptor_chain_id, partner_chain_ids, numbering=None,
                  cutoff=CONTACT_CUTOFF, restrict_to=None):
    """Every pair of residues in contact across an interface, one row per pair.

    This is the one place the contact calculation lives. Steps 04, 06 and 10 all
    reach it, so they cannot produce different answers for the same question.
    `contacts_to_partner` below is a summary of this function's output rather
    than a second calculation.

    Two residues are counted as in contact when any pair of their non-hydrogen
    atoms is within `cutoff` angstroms. Hydrogens are absent from measured
    structures, so "non-hydrogen" means the atoms we have positions for.

    Returns a list of dicts, one per contacting pair of residues:

        receptor_pos      our residue number, or the file's own number when no
                          numbering is supplied
        receptor_resnum   the number in the structure file
        receptor_aa       single-letter amino acid
        partner_chain     chain the partner residue is in
        partner_resnum    the partner's number in the structure file
        partner_aa        the partner's single-letter amino acid
        min_dist          closest approach between the two residues, angstroms
        receptor_atom     the receptor atom at that closest approach
        partner_atom      the partner atom at that closest approach
        atom_pairs        how many atom pairs fall inside the cutoff, which says
                          whether the two residues brush past each other or sit
                          against each other

    Distances are found with Bio.PDB.NeighborSearch, which sorts the atoms into a
    spatial index first, rather than by comparing every atom against every other
    atom. With a few thousand atoms on each side the direct approach means tens of
    millions of comparisons; the index gives the same answer in a fraction of the
    time.

    `numbering` is an optional `Numbering`. With one, receptor residues it cannot
    translate are skipped, which is how steps 04 and 06 restrict themselves to
    residues the alignment accounted for. Without one, the file's own numbers are
    reported unchanged, which is what a designed binder needs: its residues are
    positions in a sequence that does not exist in any archive.

    `restrict_to` optionally limits the search to a range of receptor numbers,
    given as (first, last), in whichever numbering is being reported.
    """
    partner = heavy_atoms(
        [r for cid in partner_chain_ids for r in protein_residues(model[cid])])
    if not partner:
        return []
    search = NeighborSearch(partner)

    pairs = []
    for residue in protein_residues(model[receptor_chain_id]):
        if numbering is None:
            receptor_pos = residue.id[1]
        else:
            receptor_pos = numbering.uniprot_of(residue.id[1])
            if receptor_pos is None:
                continue
        if restrict_to and not (restrict_to[0] <= receptor_pos <= restrict_to[1]):
            continue
        # Collect every atom-level hit first, then reduce to one row per pair of
        # residues, so that a pair's reported distance is its closest approach
        # rather than whichever atom happened to be looked at last.
        per_partner = {}
        for atom in residue:
            if atom.element == "H":
                continue
            for near in search.search(atom.coord, cutoff):
                distance = atom - near
                partner_residue = near.get_parent()
                key = (partner_residue.get_parent().id, partner_residue.id[1])
                entry = per_partner.get(key)
                if entry is None or distance < entry[0]:
                    per_partner[key] = (distance, atom.get_id(), near.get_id(),
                                        partner_residue,
                                        (entry[4] + 1) if entry else 1)
                else:
                    per_partner[key] = entry[:4] + (entry[4] + 1,)
        receptor_aa = THREE_TO_ONE.get(residue.get_resname(), "X")
        for (chain_id, partner_resnum), entry in per_partner.items():
            distance, receptor_atom, partner_atom, partner_residue, n = entry
            pairs.append(dict(
                receptor_pos=receptor_pos, receptor_resnum=residue.id[1],
                receptor_aa=receptor_aa, partner_chain=chain_id,
                partner_resnum=partner_resnum,
                partner_aa=THREE_TO_ONE.get(partner_residue.get_resname(), "X"),
                min_dist=distance, receptor_atom=receptor_atom,
                partner_atom=partner_atom, atom_pairs=n))
    pairs.sort(key=lambda p: (p["receptor_pos"], p["partner_chain"],
                             p["partner_resnum"]))
    return pairs


def contacts_to_partner(model, receptor_chain_id, partner_chain_ids, numbering,
                        cutoff=CONTACT_CUTOFF, restrict_to=None):
    """Which receptor residues touch the given partner chains, one row each.

    Returns {our residue number: {aa, pdb_resnum, min_dist, partner_chain,
    receptor_atom, partner_atom}}, where min_dist is the closest approach in
    angstroms to any of the partner chains.

    A summary of `contact_pairs` above rather than a second calculation: it keeps
    only the closest partner for each receptor residue. Steps 04 and 06 use this
    because the question they ask is "which of our residues does the antibody
    touch". Step 10 uses `contact_pairs` because the question it asks is which
    residue faces which, and a summary cannot answer that.
    """
    found = {}
    for pair in contact_pairs(model, receptor_chain_id, partner_chain_ids,
                              numbering, cutoff, restrict_to):
        existing = found.get(pair["receptor_pos"])
        if existing is not None and existing["min_dist"] <= pair["min_dist"]:
            continue
        found[pair["receptor_pos"]] = dict(
            aa=pair["receptor_aa"], pdb_resnum=pair["receptor_resnum"],
            min_dist=pair["min_dist"], partner_chain=pair["partner_chain"],
            receptor_atom=pair["receptor_atom"],
            partner_atom=pair["partner_atom"])
    return found


def intra_chain_contacts_outside(model, chain_id, numbering,
                                 target_range, exclude_range,
                                 cutoff=CONTACT_CUTOFF):
    """Which residues in `target_range` are touched by residues of the same chain
    from outside `exclude_range`.

    Used to ask whether other parts of the receptor fold across our epitope. The
    alternative was to name the domains doing the covering, which would make the
    answer depend on where each domain is taken to begin and end, and we are not
    confident about those boundaries. Asking "touched from outside domain III"
    needs no such judgement.

    Returns {our residue number: (distance, our number of the residue touching it)}.
    """
    residues = protein_residues(model[chain_id])
    outside, target = [], []
    for residue in residues:
        uniprot_pos = numbering.uniprot_of(residue.id[1])
        if uniprot_pos is None:
            continue
        if not (exclude_range[0] <= uniprot_pos <= exclude_range[1]):
            outside.extend(a for a in residue if a.element != "H")
        if target_range[0] <= uniprot_pos <= target_range[1]:
            target.extend((uniprot_pos, a) for a in residue
                          if a.element != "H")
    if not outside or not target:
        return {}
    search = NeighborSearch(outside)
    found = {}
    for uniprot_pos, atom in target:
        for near in search.search(atom.coord, cutoff):
            distance = atom - near
            other = numbering.uniprot_of(near.get_parent().id[1])
            if uniprot_pos not in found or distance < found[uniprot_pos][0]:
                found[uniprot_pos] = (distance, other)
    return found


def classify_exposure(rsa):
    """Turn a relative solvent accessibility value into a word."""
    if rsa >= EXPOSED_CUT:
        return "exposed"
    if rsa <= BURIED_CUT:
        return "buried"
    return "partial"


def is_borderline(rsa):
    """True if the value sits close enough to a dividing line that the word
    attached to it should not be read as precise."""
    return (abs(rsa - EXPOSED_CUT) < BORDERLINE_MARGIN
            or abs(rsa - BURIED_CUT) < BORDERLINE_MARGIN)


def anchor_role(aa):
    if aa in "DE":
        return "acidic, so the binder gets a histidine here"
    return "target histidine, so the binder gets an acidic residue here"


class CrossCheckError(AssertionError):
    """Raised when the same thing, computed twice, comes out differently."""


def cross_check_residue_set(label, computed, reference_csv, column="uniprot_pos",
                            condition_column=None, emit=print):
    """Compare a set of residue numbers against a set another script wrote out,
    and stop with an error if they differ.

    A disagreement between two scripts stops the run here, instead of sitting in
    two findings files until somebody happens to notice. That is exactly how step
    04 and step 06 came to disagree about the antibody's contacts: one of them
    looked a residue number up in the wrong direction, and because the two
    numbering ranges overlap, the wrong lookup returned a plausible answer 48
    positions away rather than failing.

    `condition_column` names a true/false column; when given, only rows where it
    reads True are counted.
    """
    path = Path(reference_csv)
    if not path.exists():
        emit(f"   Cross-check skipped: {path.name} not present yet.")
        return
    lines = path.read_text().splitlines()
    header = lines[0].split(",")
    idx = header.index(column)
    cond = header.index(condition_column) if condition_column else None
    reference = set()
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split(",")
        if cond is not None and parts[cond].strip() != "True":
            continue
        reference.add(int(parts[idx]))
    computed = set(computed)
    if computed == reference:
        emit(f"   [PASS] cross-check {label}: {len(computed)} residues, "
             f"identical to {path.name}")
        return
    only_here = sorted(computed - reference)
    only_there = sorted(reference - computed)
    message = (f"cross-check FAILED for {label}: this script and {path.name} "
               f"disagree. Only here: {only_here}. Only in {path.name}: "
               f"{only_there}.")
    emit(f"   [FAIL] {message}")
    raise CrossCheckError(message)


def load_primary_cluster(csv_path=None):
    """The anchor positions, read back from step 08's committed output.

    Step 08 decided which residues cluster on one face around H370 and wrote the
    answer to `data/derived/08-h370-clusters.csv`. That answer is computed in one
    place and read everywhere else, which is the rule in `CLAUDE.md` that came out
    of steps 04 and 06 disagreeing about the same quantity.

    This lives here rather than in step 09 because steps 09, 12 and 13 all need
    the same eight positions, and a second copy of the parsing is a second chance
    to read a different answer out of the same file.

    Returns the sorted positions of the largest cluster that contains H370.
    """
    path = Path(csv_path or (DERIVED / "08-h370-clusters.csv"))
    if not path.exists():
        raise SystemExit(
            f"{path} is missing. Run analysis/08_h370_epitope.py first; later "
            "steps take the anchor set from its output rather than keeping their "
            "own copy of it.")
    best = None
    for line in path.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        # cluster_size,max_internal_span_a,"344 358 ...",contains_h370
        _before, _, rest = line.partition('"')
        anchors, _, after = rest.partition('"')
        members = [int(tok) for tok in anchors.split()]
        if after.strip(", ").strip() != "True":
            continue
        if best is None or len(members) > len(best):
            best = members
    if best is None:
        raise SystemExit(
            f"{path} contains no cluster that includes H370. Rerun "
            "analysis/08_h370_epitope.py.")
    return sorted(best)


def load_complex(path):
    """Parse a structure, whichever of the two formats it is in.

    The design pipeline writes mmCIF (macromolecular crystallographic information
    file, the newer standard format). The structures we hold locally are PDB
    format. Synthetic test cases are written as mmCIF so that the tests exercise
    the same parser real design output will go through.

    Shared by steps 10 and 13 so the two cannot read one file differently.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    parser = MMCIFParser(QUIET=True) if suffix == ".cif" else PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))[0]


def identify_chains(model, numbering, human, emit=None):
    """Which chain is the target and which is the designed binder, by measurement.

    The design pipeline sorts its output chains so that the target takes A and the
    binder takes the next letter. That is the opposite of BindCraft version 1,
    where the designed binder was chain A, so a parser that assumes a letter will
    read the wrong molecule and report a full set of plausible nonsense.

    So the letter is not trusted. Each chain is scored by how many of its residues
    translate through the numbering step 02 measured and then read as the amino
    acid the human EGFR sequence has at that position. The target chain scores
    near one; a designed binder, whose residue numbers are positions in a sequence
    that exists nowhere, scores near zero.

    Shared by steps 10 and 13, which must agree about which chain is which.

    Returns (target chain id, [binder chain ids], evidence rows), or
    (None, [], evidence) when no chain looks like EGFR.
    """
    evidence = []
    for chain in model:
        residues = protein_residues(chain)
        if not residues:
            continue
        translatable = matching = 0
        for residue in residues:
            pos = numbering.uniprot_of(residue.id[1])
            if pos is None or pos > len(human):
                continue
            translatable += 1
            if human[pos - 1] == THREE_TO_ONE.get(residue.get_resname()):
                matching += 1
        share = matching / translatable if translatable else 0.0
        evidence.append(dict(chain=chain.id, residues=len(residues),
                             translatable=translatable, matching=matching,
                             share=share))

    if emit is not None:
        emit("   | chain | residues | translate | read as human EGFR | share |")
        emit("   |---|---|---|---|---|")
        for row in evidence:
            emit(f"   | {row['chain']} | {row['residues']} | "
                 f"{row['translatable']} | {row['matching']} | "
                 f"{row['share']:.2f} |")

    ranked = sorted(evidence, key=lambda r: -r["share"])
    if not ranked or ranked[0]["share"] < 0.8 or ranked[0]["matching"] < 3:
        return None, [], evidence
    target = ranked[0]["chain"]
    binders = [row["chain"] for row in evidence if row["chain"] != target]
    return target, binders, evidence


def receptor_chain_of(model, human=None, min_residues=20):
    """Find the chain that is EGFR in a measured structure, and its numbering.

    A structure file may hold the receptor plus whatever was bound to it, and
    which chain is which is not recorded anywhere a program can read. Each protein
    chain is therefore aligned against the human sequence and the best match wins.
    The numbering comes out of that same alignment rather than from an assumed
    offset, because a wrong offset shifts every result and raises no error.

    Returns (chain id, Numbering, percent identity, [other protein chain ids]).
    """
    human = human or human_sequence()
    scored = []
    for chain in model:
        residues = protein_residues(chain)
        if len(residues) < min_residues:
            continue
        numbering, pct = numbering_from_alignment(residues, human)
        scored.append((pct, len(residues), chain.id, numbering))
    if not scored:
        raise SystemExit("No protein chain long enough to identify in this model.")
    scored.sort(key=lambda row: (-row[0], -row[1], row[2]))
    pct, _n, chain_id, numbering = scored[0]
    others = [row[2] for row in scored[1:]]
    return chain_id, numbering, pct, others


def shell_counts(model, chain_id, numbering, centre_positions, exclude_range,
                 shells=(8.0, 12.0)):
    """How crowded each centre residue is by the same chain from outside a range.

    `intra_chain_contacts_outside` answers whether something touches a residue at
    4.5 angstroms. That is the right question for "is this covered over", and the
    wrong one for "can a protein get here": a binder is a body roughly 20
    angstroms across, so a residue can be free of contacts and still sit at the
    bottom of a cleft too tight to reach into. These wider shells measure that.

    Returns {our position: {"nearest": (distance, our position of the other
    residue), shell: count of residues inside it}}.
    """
    residues = protein_residues(model[chain_id])
    outside, centres = [], {}
    for residue in residues:
        pos = numbering.uniprot_of(residue.id[1])
        if pos is None:
            continue
        if not (exclude_range[0] <= pos <= exclude_range[1]):
            outside.append(residue)
        if pos in centre_positions:
            centres[pos] = residue
    result = {}
    for pos, residue in centres.items():
        atoms = [a for a in residue if a.element != "H"]
        row = {shell: 0 for shell in shells}
        nearest = (float("inf"), None)
        for other in outside:
            best = min((a - b) for a in atoms
                       for b in other if b.element != "H")
            other_pos = numbering.uniprot_of(other.id[1])
            if best < nearest[0]:
                nearest = (best, other_pos)
            for shell in shells:
                if best < shell:
                    row[shell] += 1
        row["nearest"] = nearest
        result[pos] = row
    return result
