"""
Shared code for the EGFR analysis scripts.

Anything computed in more than one place lives here, so the scripts cannot
disagree with each other. This file exists because two of them did disagree: step
04 and step 06 both worked out which residues the antibody touches, from the same
file with the same settings, and got different answers, because one of them looked
up a residue number in the wrong direction.

Abbreviations used in this file, expanded on first use below:
  PDB      Protein Data Bank, the public archive of measured 3D protein shapes.
           Also used to mean a file from that archive.
  UniProt  the public archive of protein sequences. Our residue numbers are
           positions in its records.
  SASA     solvent-accessible surface area: how much of a residue's surface
           water can reach, measured in square angstroms.
  RSA      relative solvent accessibility: SASA divided by the largest value
           that residue type could have. Comparable between residue types.
  CA       the "alpha carbon", one specific atom present in every residue. Often
           used as a single stand-in for the whole residue's position.
  Fab      the gripping arm of an antibody, separated from the rest of it.
"""

from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import NeighborSearch
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

# The eight positions we could build a pH switch against, with the amino acid
# each one is. D and E are acidic and always carry a negative charge, so we put a
# histidine opposite them. H is histidine, which changes charge with pH, so we put
# an acidic residue opposite it.
ANCHORS = {416: "D", 418: "H", 421: "E", 424: "E",
           433: "H", 455: "E", 458: "D", 460: "D"}

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

    Structure files carry their own residue numbers and there is no rule about
    which convention they follow. Files of secreted proteins often count from the
    mature protein, whose residue 1 is UniProt residue 25, because that is the
    molecule that was actually crystallised. Nothing in the file records the
    choice. Both structures used in this project turn out to be off by 24.

    This class exists because the two directions were previously passed around as
    plain dictionaries, and step 06 indexed a UniProt-to-structure dictionary with
    a structure number. Since the two ranges overlap, that lookup succeeded and
    returned a number 48 positions away from the right one, with no error raised.
    Asking for `uniprot_of(...)` or `pdb_of(...)` by name cannot go wrong in the
    same way, and passing the wrong kind of number now returns nothing rather than
    a plausible wrong answer.
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


def protein_residues(chain):
    return [r for r in chain if is_aa(r, standard=True)]


def heavy_atoms(residues):
    """Every atom except hydrogen.

    Hydrogen atoms are too light to show up in most measured structures and are
    simply absent from these files, so "heavy atom" means the atoms we actually
    have positions for.
    """
    return [a for r in residues for a in r if a.element != "H"]


def contacts_to_partner(model, receptor_chain_id, partner_chain_ids, numbering,
                        cutoff=CONTACT_CUTOFF, restrict_to=None):
    """Which receptor residues touch the given partner chains.

    Returns {our residue number: {aa, pdb_resnum, min_dist, partner_chain,
    receptor_atom, partner_atom}}, where min_dist is the closest approach in
    angstroms.

    This is the one place the contact calculation lives. Step 04 and step 06 both
    call it, so they cannot produce different answers for the same question.

    Distances are found with Bio.PDB.NeighborSearch, which sorts the atoms into a
    spatial index first, rather than by comparing every atom against every other
    atom. With a few thousand atoms on each side the direct approach means tens of
    millions of comparisons; the index gives the same answer in a fraction of the
    time.

    `restrict_to` optionally limits the search to a range of our residue numbers,
    given as (first, last).
    """
    partner = heavy_atoms(
        [r for cid in partner_chain_ids for r in protein_residues(model[cid])])
    if not partner:
        return {}
    search = NeighborSearch(partner)

    found = {}
    for residue in protein_residues(model[receptor_chain_id]):
        uniprot_pos = numbering.uniprot_of(residue.id[1])
        if uniprot_pos is None:
            continue
        if restrict_to and not (restrict_to[0] <= uniprot_pos <= restrict_to[1]):
            continue
        best = None
        for atom in residue:
            if atom.element == "H":
                continue
            for near in search.search(atom.coord, cutoff):
                distance = atom - near
                if best is None or distance < best[0]:
                    best = (distance, near.get_parent().get_parent().id,
                            atom.get_id(), near.get_id())
        if best is not None:
            found[uniprot_pos] = dict(
                aa=THREE_TO_ONE.get(residue.get_resname(), "X"),
                pdb_resnum=residue.id[1], min_dist=best[0],
                partner_chain=best[1], receptor_atom=best[2],
                partner_atom=best[3])
    return found


def intra_chain_contacts_outside(model, chain_id, numbering,
                                 target_range, exclude_range,
                                 cutoff=CONTACT_CUTOFF):
    """Which residues in `target_range` are touched by residues of the same chain
    from outside `exclude_range`.

    Used to ask whether other parts of the receptor fold across our epitope.
    Phrasing it as "touched from outside domain III" avoids having to decide where
    each domain begins and ends, so the answer does not depend on boundaries we
    are not certain about.

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

    The point is that a disagreement between two scripts should stop the run
    rather than sit in two findings files waiting for someone to notice. Step 04
    and step 06 disagreed about the antibody's contacts for exactly that reason.

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
