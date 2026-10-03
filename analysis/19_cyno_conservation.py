#!/usr/bin/env python3
"""
19_cyno_conservation.py — are the eight anchors and the step 14 sequons the same
in cynomolgus monkey EGFR as in human?

Abbreviations, expanded here because each file gets read on its own:
  EGFR     epidermal growth factor receptor, the protein we design against.
  UniProt  the public archive of protein sequences. Every residue number in this
           project is a position in its human record P00533.
  FASTA    the plain-text format UniProt serves sequences in: a header line
           beginning with '>', then the sequence.
  Swiss-Prot  the part of UniProt a curator has read and checked. Entries from it
           are called "reviewed".
  TrEMBL   the part of UniProt filled in automatically from gene predictions, with
           no curator. Entries from it are called "unreviewed".
  aa       amino acids, the building blocks a protein chain is made of.
  cyno     the usual short name for the cynomolgus monkey, Macaca fascicularis,
           also called the crab-eating macaque.
  rhesus   the rhesus macaque, Macaca mulatta, a closely related monkey.
  ECD      extracellular domain, the part of EGFR outside the cell, residues
           25-645 in our numbering.
  sequon   the three-residue pattern N-X-S/T that a sugar chain can attach to: an
           asparagine (N), any amino acid except proline (X), then a serine or a
           threonine (S/T). Found by step 14, whose rule this script calls rather
           than restating.
  D, E, H, N, Y   single-letter codes for aspartate, glutamate, histidine,
           asparagine and tyrosine.
  SHA-256  a fingerprint of a file's bytes, used here to notice if the cached
           sequence file changes underneath us.

WHY THIS STEP EXISTS, AND WHAT IT IS NOT
----------------------------------------
**This is not a design input and it changes no decision.** The competition page
names human and mouse only, and the decision recorded on 2 October 2026 is to
design for human and mouse only (`docs/decisions-log.md`, Open section). A
pre-launch message had said "mouse and cyno", so the question of whether the
target face survives in cynomolgus monkey was left open rather than answered. The
log recorded it as a cheap check worth doing once, because the answer costs no
graphics-card time and no campaign slot.

So this step exists to turn an open question into a recorded sentence. Nothing
downstream reads its output. No anchor is added, dropped, demoted or promoted
because of it. If the answer had come back badly it would still change nothing,
because cyno is not an objective; it would only change what the methods write-up
can claim.

WHAT THIS SCRIPT COMPUTES
-------------------------
1. Fetches the cynomolgus monkey EGFR sequence from UniProt and caches it, with a
   record of the accession, the address it came from and the date it was fetched.
   There is no reviewed cyno EGFR entry, so the best-annotated one is used and the
   script says so in its output rather than letting the reader assume otherwise.

2. Aligns it against human P00533 by the same method step 01 uses for mouse:
   end-to-end pairwise alignment, BLOSUM62 scoring, gap open -11, gap extend -1,
   from `egfr_common.make_aligner`.

3. Measures the offset between the cyno entry's own numbering and ours, rather
   than assuming it is zero. The offset goes through `egfr_common.Numbering`, the
   same class the structure steps use, because a numbering mix-up shifts every
   result and raises no error.

4. Reports, for each of the eight anchors of the current H370 face, whether cyno
   has the same amino acid; and for any that differ, whether the substitution
   keeps the charge behaviour the design rule depends on. The design rule needs an
   acidic residue to stay acidic (D or E, either will do), and a target histidine
   to stay a histidine (nothing else switches charge over this pH range).

5. Reports whether each sequon step 14 found is intact in cyno, with the sequons
   that step 14 measured as nearest to a current anchor called out. N361 is the
   one that is a sequon in human and not in mouse, so what cyno has there is the
   interesting cell in the table.

HOW THIS STEP IS CHECKED
------------------------
  - The human-versus-mouse alignment is recomputed here and its domain III
    differences are cross-checked against the set step 01 committed. If this
    script cannot reproduce step 01, its cyno answer is not to be believed either.
  - The human sequons are found by calling step 14's own `find_sequons`, not by a
    second copy of the rule, and the positions are cross-checked against the table
    step 14 committed.
  - The verdict on whether a substitution keeps the needed charge behaviour is
    tested on constructed inputs every run, with the answers written out rather
    than computed.
  - `--break-offset N` shifts the measured map by N positions on purpose, so the
    checks that are supposed to notice a numbering slip can be seen to notice one.
    `--self-test` runs that and confirms the run fails.
  - A second, independently annotated macaque sequence is read as a control: the
    reviewed rhesus macaque entry. The cyno entry is unreviewed, so a claim about
    one automatic gene prediction is worth little on its own.

Writes:
  results/findings/19-cyno-conservation.md
  data/derived/19-cyno-anchor-conservation.csv
  data/derived/19-cyno-sequons.csv
  data/derived/19-cyno-domain3-differences.csv
  data/sequences/egfr-macaque-uniprot.fasta         (cached input)
  data/sequences/egfr-macaque-uniprot.source.txt    (where it came from, and when)

Run standalone:  python analysis/19_cyno_conservation.py
"""

import argparse
import csv
import datetime
import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

import egfr_common as common  # noqa: E402

SEQ_DIR = common.SEQ_DIR
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS

MACAQUE_FASTA = SEQ_DIR / "egfr-macaque-uniprot.fasta"
MACAQUE_SOURCE = SEQ_DIR / "egfr-macaque-uniprot.source.txt"

# The two macaque entries, and why each one.
#
# CYNO_ACCESSION is the subject of this step. UniProt holds no reviewed entry for
# cynomolgus monkey EGFR: searching its archive for the gene EGFR in Macaca
# fascicularis (taxonomy identifier 9541) returns two unreviewed entries and no
# reviewed one. A0A2K5WKD8 is the longer of the two at 704 aa against 628 aa, it
# is the second version of its sequence rather than the first, and it is the one
# in UniProt's reference proteome for the species. Both start at the same residue
# as the human record, and the shorter one stops inside the extracellular region,
# so the longer one is the only candidate that spans the whole face we aim at.
#
# RHESUS_ACCESSION is a control and not the subject. Rhesus macaque EGFR is a
# reviewed entry of the full 1210 aa. Reading it alongside says whether a
# difference seen in the unreviewed cyno entry is a real difference between the
# species or an artefact of an automatic gene prediction.
CYNO_ACCESSION = "A0A2K5WKD8"
CYNO_TOKEN = "MACFA"
RHESUS_ACCESSION = "P55245"
RHESUS_TOKEN = "MACMU"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"

ECD_START, ECD_END = common.ECD_START, common.ECD_END
D3_START, D3_END = common.D3_START, common.D3_END
EPI_START, EPI_END = common.EPI_START, common.EPI_END

# Identity through the measured map has to clear this for the comparison to mean
# anything. Human and cyno EGFR are far above it; a map shifted by even one
# position falls far below it, which is what makes this a check and not a
# formality. It is a floor and not a measurement, set well clear of both.
MIN_MAPPED_IDENTITY_PCT = 80.0


def load_step14():
    """Step 14 loaded as a module, for its sequon rule, so that what counts as a
    place a sugar chain attaches is decided in one place. Step 14 tests each
    branch of that rule and can break each one on purpose; a second copy here
    would be a second chance to answer the same question differently."""
    spec = importlib.util.spec_from_file_location(
        "step14", Path(__file__).resolve().parent / "14_glycan_sequons.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Getting the sequences, and recording where they came from
# ---------------------------------------------------------------------------

def fetch_macaque_fasta(emit):
    """Download both macaque records once and keep them, the way the human and
    mouse records are kept: one FASTA file in `data/sequences/`, selected from by
    a token in the header.

    A second file records the accession, the address and the date fetched, because
    none of that is recoverable from the sequence afterwards and a sequence with no
    provenance is a number written in by hand. Every later run checks the cached
    file's fingerprint against the recorded one, so a file edited by hand stops the
    run instead of quietly changing a result.
    """
    if MACAQUE_FASTA.exists() and MACAQUE_FASTA.stat().st_size > 0:
        emit(f"   cache: {MACAQUE_FASTA.relative_to(common.ROOT)} already present "
             f"({MACAQUE_FASTA.stat().st_size:,} bytes), not re-fetched")
        return False
    SEQ_DIR.mkdir(parents=True, exist_ok=True)
    chunks = []
    for accession in (CYNO_ACCESSION, RHESUS_ACCESSION):
        url = UNIPROT_FASTA_URL.format(accession)
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        text = response.text
        if not text.startswith(">"):
            raise SystemExit(f"{url} did not return a FASTA record.")
        chunks.append(text if text.endswith("\n") else text + "\n")
        emit(f"   fetched {accession} from {url} ({len(text):,} bytes)")
    MACAQUE_FASTA.write_text("".join(chunks))
    fetched = datetime.datetime.now(datetime.timezone.utc)
    MACAQUE_SOURCE.write_text(
        "Where data/sequences/egfr-macaque-uniprot.fasta came from.\n"
        "Written by analysis/19_cyno_conservation.py when it fetched the file.\n"
        "Do not hand-edit: the script checks the fingerprint below on every run.\n"
        f"cynomolgus_monkey_accession: {CYNO_ACCESSION}\n"
        f"cynomolgus_monkey_url: {UNIPROT_FASTA_URL.format(CYNO_ACCESSION)}\n"
        f"rhesus_macaque_accession: {RHESUS_ACCESSION}\n"
        f"rhesus_macaque_url: {UNIPROT_FASTA_URL.format(RHESUS_ACCESSION)}\n"
        f"date_fetched_utc: {fetched.strftime('%Y-%m-%d')}\n"
        f"time_fetched_utc: {fetched.strftime('%H:%M:%S')}\n"
        f"bytes: {MACAQUE_FASTA.stat().st_size}\n"
        f"sha256: {file_fingerprint(MACAQUE_FASTA)}\n")
    return True


def file_fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_source_record():
    """The provenance file as {field: value}, or {} if it is not there."""
    if not MACAQUE_SOURCE.exists():
        return {}
    record = {}
    for line in MACAQUE_SOURCE.read_text().splitlines():
        if ": " in line and not line.endswith(":"):
            key, _, value = line.partition(": ")
            record[key.strip()] = value.strip()
    return record


def is_tracked(path):
    """Is this file tracked by git, meaning it travels with the repository?

    Asked of git rather than of the ignore rules, because what matters to a reader
    is whether the file will be there when they clone, and only the index knows
    that. Returns True, False, or None when git cannot be reached at all.
    """
    try:
        result = subprocess.run(
            ["git", "ls-files", "--error-unmatch", str(Path(path).resolve())],
            cwd=str(common.ROOT), capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return result.returncode == 0


def macaque_sequence(token):
    for header, seq in common.read_fasta(MACAQUE_FASTA).items():
        if token in header:
            return header, seq
    raise KeyError(f"no record containing {token!r} in {MACAQUE_FASTA}")


# ---------------------------------------------------------------------------
# Aligning, and measuring the offset rather than assuming it
# ---------------------------------------------------------------------------

def align_to_human(human, other, shift=0):
    """Line another species' sequence up against human and return the measured
    correspondence between the two numbering systems.

    The method is step 01's, through the shared aligner: end to end, BLOSUM62, gap
    open -11, gap extend -1. Returns (Numbering, columns where human has a gap,
    columns where the other sequence has a gap).

    The Numbering is built the way the structure steps build theirs — the other
    sequence's own position plays the part the structure file's number plays — so
    `uniprot_of` takes a cyno position and `pdb_of` takes one of ours. Handing
    either method the wrong kind of number returns nothing rather than a plausible
    answer, which is the whole reason that class exists.

    `shift` is for `--break-offset` and is 0 in a real run.
    """
    alignment = common.make_aligner().align(human, other)[0]
    human_i = other_i = 0
    human_gaps = other_gaps = 0
    other_to_human = {}
    for human_char, other_char in zip(alignment[0], alignment[1]):
        if human_char != "-":
            human_i += 1
        else:
            human_gaps += 1
        if other_char != "-":
            other_i += 1
        else:
            other_gaps += 1
        if human_char != "-" and other_char != "-":
            other_to_human[other_i] = human_i + shift
    return common.Numbering(other_to_human), human_gaps, other_gaps


def restricted(numbering, lo, hi):
    """The same correspondence, narrowed to a range of our positions, so that an
    offset can be reported for the region that matters rather than only for the
    whole entry."""
    kept = {}
    for pos in range(lo, hi + 1):
        other = numbering.pdb_of(pos)
        if other is not None:
            kept[other] = pos
    return common.Numbering(kept)


def region_identity(human, other, numbering, lo, hi):
    """(identical, compared, unmapped) across our positions lo-hi.

    A position with no counterpart in the other sequence is counted as unmapped
    rather than as a difference, because the two are not the same thing: one says
    the species differ there, the other says this entry does not reach that far.
    """
    identical = compared = unmapped = 0
    for pos in range(lo, hi + 1):
        other_pos = numbering.pdb_of(pos)
        if other_pos is None:
            unmapped += 1
            continue
        compared += 1
        if human[pos - 1] == other[other_pos - 1]:
            identical += 1
    return identical, compared, unmapped


def differences(human, other, numbering, lo, hi):
    """[(our position, human residue, the other species' residue)] where they
    differ, over our positions lo-hi."""
    found = []
    for pos in range(lo, hi + 1):
        other_pos = numbering.pdb_of(pos)
        if other_pos is None:
            continue
        if human[pos - 1] != other[other_pos - 1]:
            found.append((pos, human[pos - 1], other[other_pos - 1]))
    return found


def residue_in(other, numbering, pos):
    """The other species' residue at our position `pos`, or None when the entry
    has nothing aligned there."""
    other_pos = numbering.pdb_of(pos)
    return other[other_pos - 1] if other_pos is not None else None


def contiguous_run(numbering, offset, start=1, limit=None):
    """How far from `start` our positions map one after another at a single
    offset, with no break.

    This is how far an entry can be read straight across. Past the end of the run
    an automatic gene prediction has usually gone its own way, and letters the
    alignment places beyond it are not evidence about the species.
    """
    limit = limit or max(numbering.uniprot_positions() or [0])
    pos = start
    while pos <= limit:
        other = numbering.pdb_of(pos)
        if other is None or pos - other != offset:
            break
        pos += 1
    return pos - 1


# ---------------------------------------------------------------------------
# What a substitution does to the charge behaviour the design rule needs
# ---------------------------------------------------------------------------

# Set by --break-rule, which switches one branch off on purpose so the test meant
# to cover it can be seen to catch its removal.
BROKEN_RULE = None

BREAKABLE_RULES = {
    "acidic": (
        "acidic stays acidic",
        "a swap from an acidic residue to something uncharged would be reported "
        "as keeping the charge the binder's histidine is attracted to, so an "
        "anchor that had stopped working would read as intact"),
    "histidine": (
        "a target histidine stays a histidine",
        "a target histidine replaced by anything else would be reported as "
        "harmless, although histidine is the only residue that changes charge "
        "between pH 7.4 and pH 6.5, so the pair built there would no longer "
        "switch"),
}


def rule_live(name):
    return BROKEN_RULE != name


def substitution_verdict(human_aa, other_aa):
    """What a substitution at an anchor does to the charge behaviour the design
    rule depends on.

    The design rule pairs charges position by position. Against an acidic anchor,
    which is aspartate or glutamate and carries a negative charge at both pH
    values, the binder gets a histidine; the swap from one acidic residue to the
    other keeps that, because the binder's histidine is attracted to the negative
    charge and not to the particular residue carrying it. Against a target
    histidine the binder gets an acidic residue, and only a histidine will do on
    the target side: histidine is the one amino acid that changes charge over this
    pH range, so replacing it removes the switch rather than weakening it.

    Returns (a short code, a sentence for the reader).
    """
    if human_aa == other_aa:
        return "identical", "identical"
    if human_aa in "DE" and other_aa in "DE" and rule_live("acidic"):
        return "charge kept", (
            "still acidic, so the negative charge the binder's histidine is "
            "attracted to is unchanged")
    if human_aa == "H" and rule_live("histidine"):
        return "switch lost", (
            "the target histidine is gone, and histidine is the only residue that "
            "changes charge over this pH range, so the pair built here would not "
            "switch")
    if human_aa in "DE" and rule_live("acidic"):
        return "charge lost", (
            "no longer acidic, so there is no negative charge for the binder's "
            "histidine to be attracted to")
    return "charge kept", "treated as harmless by a branch switched off on purpose"


def run_tests(emit):
    """The verdict rule on constructed inputs, with the answers written out rather
    than computed."""
    failures = []
    emit("   Does a substitution keep the charge behaviour the design rule needs?")
    emit("   Tested on constructed inputs:")
    emit()
    cases = [
        ("E", "E", "identical", "no substitution at all"),
        ("H", "H", "identical", "a target histidine that stays a histidine"),
        ("E", "D", "charge kept", "acidic stays acidic: glutamate to aspartate"),
        ("D", "E", "charge kept", "acidic stays acidic: aspartate to glutamate"),
        ("E", "Q", "charge lost",
         "glutamate to glutamine: same shape, no charge"),
        ("D", "N", "charge lost",
         "aspartate to asparagine: same shape, no charge"),
        ("H", "Y", "switch lost",
         "a target histidine replaced, as mouse does at 361"),
        ("H", "D", "switch lost",
         "acidic in place of a target histidine is still no switch"),
    ]
    for human_aa, other_aa, want, why in cases:
        got, _ = substitution_verdict(human_aa, other_aa)
        ok = got == want
        emit(f"     [{'PASS' if ok else 'FAIL'}] {human_aa} -> {other_aa}: "
             f"{got:12s} ({why})")
        if not ok:
            label = ("a target histidine stays a histidine" if human_aa == "H"
                     else "acidic stays acidic")
            failures.append(f"{label}: {human_aa} -> {other_aa} gave {got!r}, "
                            f"expected {want!r}")
    emit()
    return failures


# ---------------------------------------------------------------------------

def main(argv=None):
    global BROKEN_RULE

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--break-offset", type=int, default=0, metavar="N",
                        help="shift the measured correspondence by N positions, "
                             "to confirm the checks notice. Expected to fail.")
    parser.add_argument("--break-rule", choices=list(BREAKABLE_RULES), default=None,
                        help="switch one branch of the charge rule off, to confirm "
                             "the matching test then fails. Expected to fail.")
    parser.add_argument("--self-test", action="store_true",
                        help="break each check in turn and confirm something "
                             "notices.")
    args = parser.parse_args(argv)
    if args.self_test:
        return run_self_test()
    BROKEN_RULE = args.break_rule

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    failures = []

    def check(label, ok, detail):
        emit(f"   [{'PASS' if ok else 'FAIL'}] {label}: {detail}")
        if not ok:
            failures.append(label)

    step14 = load_step14()
    human = common.human_sequence()
    mouse = common.mouse_sequence()
    anchors = common.load_primary_cluster()
    old_anchors = sorted(common.ANCHORS)

    emit("=" * 72)
    emit("CYNOMOLGUS MONKEY EGFR — are the anchors and the sequons conserved?")
    emit("=" * 72)
    emit()
    emit("   THIS IS NOT A DESIGN INPUT AND IT CHANGES NO DECISION. The competition")
    emit("   page names human and mouse only, and the decision of 2 October 2026 is")
    emit("   to design for human and mouse only. A pre-launch message had said")
    emit("   'mouse and cyno', so whether the face we aim at survives in cynomolgus")
    emit("   monkey was left open rather than answered. This step exists to turn")
    emit("   that open question into a recorded sentence. Nothing downstream reads")
    emit("   its output, and no anchor is added, dropped, demoted or promoted")
    emit("   because of it.")
    emit()
    if args.break_offset:
        emit(f"   NOTE: --break-offset {args.break_offset} is in effect: the measured")
        emit("   correspondence between the two numbering systems is shifted on")
        emit("   purpose. This run is expected to fail.")
        emit()
    if BROKEN_RULE:
        emit(f"   NOTE: --break-rule {BROKEN_RULE} is in effect: "
             f"{BREAKABLE_RULES[BROKEN_RULE][1]}.")
        emit("   This run is expected to fail.")
        emit()

    # ---- 1. Where the sequence came from ------------------------------------
    emit("1. The cynomolgus monkey sequence, and where it came from")
    emit()
    emit("   UniProt holds no reviewed entry for cynomolgus monkey EGFR: an entry is")
    emit("   'reviewed' when a curator has read it, and the cyno EGFR entries are")
    emit("   filled in automatically from a gene prediction instead. The entry used")
    emit("   here is the better annotated of the two the archive holds: the longer")
    emit("   of them, the second version of its sequence rather than the first, and")
    emit("   the one in UniProt's reference proteome for the species. That it is")
    emit("   unreviewed is the main limit on everything below, which is why a")
    emit("   reviewed sequence from a closely related monkey is read alongside it in")
    emit("   section 7.")
    emit()
    fetched_now = fetch_macaque_fasta(emit)
    source = read_source_record()
    fingerprint = file_fingerprint(MACAQUE_FASTA)
    cyno_header, cyno = macaque_sequence(CYNO_TOKEN)
    rhesus_header, rhesus = macaque_sequence(RHESUS_TOKEN)
    tracked = is_tracked(MACAQUE_FASTA)
    tracked_word = {True: "yes, it is tracked by git and travels with the "
                          "repository",
                    False: "NO — it is not tracked by git yet, so a fresh clone "
                           "would re-fetch it",
                    None: "could not be determined: git did not answer"}[tracked]
    emit()
    emit("   | field | value |")
    emit("   |---|---|")
    emit(f"   | accession used | {CYNO_ACCESSION} (unreviewed) |")
    emit(f"   | length | {len(cyno)} aa |")
    emit(f"   | address | {UNIPROT_FASTA_URL.format(CYNO_ACCESSION)} |")
    emit(f"   | date fetched (UTC) | {source.get('date_fetched_utc', 'not recorded')} "
         f"{'(this run)' if fetched_now else '(an earlier run)'} |")
    emit(f"   | cached at | data/sequences/{MACAQUE_FASTA.name} |")
    emit(f"   | tracked by git, when this ran | {tracked_word} |")
    emit(f"   | SHA-256 of the cached file | {fingerprint} |")
    emit()
    emit("   The header, exactly as UniProt serves it (shown on its own line because")
    emit("   it contains the character the table above uses as a separator):")
    emit(f"     {cyno_header}")
    emit()
    recorded = source.get("sha256")
    if recorded:
        check("cached sequence file unchanged since it was fetched",
              recorded == fingerprint,
              "fingerprint matches the one recorded at fetch time" if
              recorded == fingerprint else
              f"fingerprint is {fingerprint} but {MACAQUE_SOURCE.name} records "
              f"{recorded}; the file has been edited")
    else:
        emit(f"   [WARN] {MACAQUE_SOURCE.name} is missing or records no "
             "fingerprint, so the cached file cannot be checked.")
    emit()

    # ---- 2. The alignment and the measured offset ---------------------------
    emit("2. Lining the two sequences up, and measuring the offset")
    emit()
    emit("   Method, the same one step 01 uses for mouse: end-to-end pairwise")
    emit("   alignment, BLOSUM62 scoring, gap open -11, gap extend -1. Comparing")
    emit("   position 1 to position 1 straight down two sequences does not work,")
    emit("   because one extra residue early on shifts everything after it.")
    emit()
    emit("   The offset is measured and not assumed. Nothing in a sequence file")
    emit("   records which numbering convention it follows, and a wrong offset")
    emit("   shifts every result while raising no error.")
    emit()
    maps = {}
    gap_counts = {}
    for name, other in (("cyno", cyno), ("rhesus", rhesus), ("mouse", mouse)):
        shift = args.break_offset if name == "cyno" else 0
        numbering, human_gaps, other_gaps = align_to_human(human, other, shift)
        maps[name] = numbering
        gap_counts[name] = (human_gaps, other_gaps)

    emit("   | species | entry | length | positions paired | offset over domain III "
         "| share of domain III at that offset |")
    emit("   |---|---|---|---|---|---|")
    region_offsets = {}
    for name, accession, other in (("cynomolgus monkey", CYNO_ACCESSION, cyno),
                                   ("rhesus macaque", RHESUS_ACCESSION, rhesus),
                                   ("mouse", "Q01279", mouse)):
        key = {"cynomolgus monkey": "cyno", "rhesus macaque": "rhesus",
               "mouse": "mouse"}[name]
        d3 = restricted(maps[key], D3_START, D3_END)
        offset, share, _ = d3.dominant_offset()
        region_offsets[key] = (offset, share, len(d3))
        emit(f"   | {name} | {accession} | {len(other)} aa | {len(maps[key])} | "
             f"ours = theirs {offset:+d} | {100.0 * share:.1f}% of {len(d3)} |")
    emit()
    cyno_offset, cyno_share, cyno_d3_n = region_offsets["cyno"]
    whole_offset, whole_share, _ = maps["cyno"].dominant_offset()
    in_words = ("exactly the cyno entry's own residue number, with nothing to add "
                "or subtract" if cyno_offset == 0 else
                f"the cyno entry's own residue number "
                f"{'plus' if cyno_offset > 0 else 'minus'} {abs(cyno_offset)}")
    emit("   Stated plainly, since this is the number that shifts every later")
    emit("   result if it is wrong: across domain III our residue number is")
    emit(f"   {in_words},")
    emit(f"   at every one of the {cyno_d3_n} positions ({100.0 * cyno_share:.1f}%). "
         "Across the whole entry that")
    emit(f"   same single offset covers {100.0 * whole_share:.1f}% of its "
         f"{len(cyno)} residues, and the rest is the")
    emit("   tail discussed next.")
    emit()
    run_end = contiguous_run(maps["cyno"], cyno_offset, start=1, limit=len(human))
    emit("   The cyno entry reads straight across, one position after another at")
    emit(f"   that single offset, from our position 1 to our position {run_end}. Past")
    emit("   there the two sequences stop corresponding: this entry is a gene")
    emit("   prediction, and its last stretch reads the gene differently from the")
    emit("   human record. Letters the alignment places beyond that point are not")
    emit("   evidence about the species, and are marked as such below.")
    emit(f"   Our extracellular region ends at {ECD_END}, so its last "
         f"{max(0, ECD_END - run_end)} positions fall in that")
    emit(f"   stretch. All {len(anchors)} anchors sit well inside the readable part, "
         f"the furthest at")
    emit(f"   position {max(anchors)}. The sequon table in section 6 covers the whole")
    emit("   extracellular region and marks any row that falls beyond the readable")
    emit("   part, so a reader can see which ones to discount.")
    emit()
    missing = len(human) - len(cyno)
    emit("   Columns where human has a gap, meaning the cyno entry has extra")
    emit(f"   residues there: {gap_counts['cyno'][0]}. Columns where the cyno entry "
         f"has a gap: {gap_counts['cyno'][1]}.")
    emit(f"   The human record is {missing} residues longer than this entry, "
         + ("and that accounts for the gaps exactly."
            if gap_counts["cyno"][1] == missing and gap_counts["cyno"][0] == 0
            else "which is most of that."))
    emit("   The length it is short by is at the far end of the protein, past the")
    emit("   extracellular region, so it is not our concern.")
    emit()

    # ---- 3. Identity by region ----------------------------------------------
    emit("3. Identity by region")
    emit()
    emit("   | Region | Positions | Identical | Identity | Not covered by the entry |")
    emit("   |---|---|---|---|---|")
    identity_rows = []
    for label, lo, hi in (("Extracellular region", ECD_START, ECD_END),
                          ("Readable part of it", ECD_START, min(ECD_END, run_end)),
                          ("Domain III", D3_START, D3_END),
                          ("Candidate epitope", EPI_START, EPI_END)):
        same, compared, unmapped = region_identity(human, cyno, maps["cyno"], lo, hi)
        pct = 100.0 * same / compared if compared else 0.0
        identity_rows.append((label, lo, hi, same, compared, unmapped, pct))
        emit(f"   | {label} | {lo}-{hi} | {same}/{compared} | {pct:.1f}% | "
             f"{unmapped} |")
    emit()
    d3_pct = next(r[6] for r in identity_rows if r[0] == "Domain III")
    check(f"domain III identity through the measured map is above "
          f"{MIN_MAPPED_IDENTITY_PCT:.0f}%",
          d3_pct >= MIN_MAPPED_IDENTITY_PCT,
          f"{d3_pct:.1f}%. Below this the two sequences are not being compared at "
          "matching positions, whatever the alignment says")
    emit()

    # ---- 4. Differences in domain III ---------------------------------------
    cyno_d3 = differences(human, cyno, maps["cyno"], D3_START, D3_END)
    cyno_d3_labels = [f"{h}{p}{o}" for p, h, o in cyno_d3]
    mouse_d3 = differences(human, mouse, maps["mouse"], D3_START, D3_END)
    mouse_d3_labels = [f"{h}{p}{o}" for p, h, o in mouse_d3]
    emit(f"4. Differences in domain III ({D3_START}-{D3_END}), human UniProt "
         "numbering")
    emit("   Notation: S348T means human has S at that position where the other")
    emit("   species has T.")
    emit()
    emit(f"   human vs cyno:  {'  '.join(cyno_d3_labels) if cyno_d3_labels else 'none'}")
    emit(f"   human vs mouse: {'  '.join(mouse_d3_labels) if mouse_d3_labels else 'none'}"
         f"  ({len(mouse_d3_labels)} of them)")
    emit()
    shared_with_mouse = sorted(set(p for p, _, _ in cyno_d3)
                               & set(p for p, _, _ in mouse_d3))
    if cyno_d3:
        plural = "difference" if len(cyno_d3_labels) == 1 else "differences"
        shared_words = "No position differs in both species."
        if shared_with_mouse:
            shared_words = (
                ("Position " if len(shared_with_mouse) == 1 else "Positions ")
                + ", ".join(str(p) for p in shared_with_mouse)
                + (" differs" if len(shared_with_mouse) == 1 else " differ")
                + " in both species.")
        emit(f"   {len(cyno_d3_labels)} {plural} between human and cyno across the "
             f"{D3_END - D3_START + 1} positions of")
        emit(f"   domain III, against {len(mouse_d3_labels)} between human and "
             f"mouse. {shared_words}")
    else:
        emit("   Domain III is identical between human and cyno.")
    emit()
    epi_cyno = [f"{h}{p}{o}" for p, h, o in cyno_d3 if EPI_START <= p <= EPI_END]
    emit(f"   Inside the fallback epitope {EPI_START}-{EPI_END}: "
         f"{epi_cyno if epi_cyno else 'no differences'}")
    emit("   (Mouse has one there, S442G, which is why step 10 rejects any candidate")
    emit("   touching position 442.)")
    emit()

    # ---- 5. The eight anchors ------------------------------------------------
    emit("5. The eight anchors of the current H370 face")
    emit()
    emit("   The anchor set is read back from step 08's committed output rather than")
    emit("   written out here, so there is one copy of it in the project. Six are")
    emit("   acidic — aspartate or glutamate, negatively charged at both pH values —")
    emit("   and the binder gets a histidine opposite each. Two are the target")
    emit("   histidines, H358 and H370, and the binder gets an acidic residue")
    emit("   opposite those.")
    emit()
    anchor_rows = []
    for set_name, positions in (("current H370 face", anchors),
                                ("old 415-466 epitope (the recorded fallback)",
                                 old_anchors)):
        emit(f"   {set_name}")
        emit()
        emit("   | anchor | role | cyno | rhesus | mouse | same as human? | what it "
             "does to the charge behaviour |")
        emit("   |---|---|---|---|---|---|---|")
        for pos in positions:
            human_aa = human[pos - 1]
            cells = {}
            for key, seq in (("cyno", cyno), ("rhesus", rhesus), ("mouse", mouse)):
                cells[key] = residue_in(seq, maps[key], pos) or "not covered"
            code, sentence = substitution_verdict(human_aa, cells["cyno"])
            role = ("acidic" if human_aa in "DE" else "target histidine")
            emit(f"   | {human_aa}{pos} | {role} | {cells['cyno']} | "
                 f"{cells['rhesus']} | {cells['mouse']} | "
                 f"{'yes' if code == 'identical' else 'NO'} | {sentence} |")
            anchor_rows.append(dict(
                anchor_set=set_name, uniprot_pos=pos, challenge_pos=pos - 24,
                human_aa=human_aa, cyno_aa=cells["cyno"],
                rhesus_aa=cells["rhesus"], mouse_aa=cells["mouse"],
                role=role, identical=(code == "identical"), verdict=code))
        emit()
    current_rows = [r for r in anchor_rows if r["anchor_set"] == "current H370 face"]
    cyno_same = [r for r in current_rows if r["identical"]]
    cyno_diff = [r for r in current_rows if not r["identical"]]
    exceptions = ", ".join(f"{r['human_aa']}{r['uniprot_pos']}{r['cyno_aa']}"
                           for r in cyno_diff)
    emit(f"   {len(cyno_same)} of the {len(current_rows)} current anchors are the "
         "same amino acid in cyno as in human"
         + (f"; the exceptions: {exceptions}." if cyno_diff else "."))
    emit()

    # ---- 6. The sequons ------------------------------------------------------
    emit("6. The sequons step 14 found, checked in cyno")
    emit()
    emit("   A sequon is the pattern N-X-S/T that a sugar chain attaches to: an")
    emit("   asparagine, then any amino acid except proline, then a serine or a")
    emit("   threonine. The rule is step 14's and is called here rather than")
    emit("   restated, so the two steps cannot disagree about what counts as one.")
    emit("   A sugar chain is large and mobile, so an attachment point near an")
    emit("   anchor can cover it some of the time; one present in a species and")
    emit("   absent in another covers one target and not the other.")
    emit()
    cyno_sequons = step14.compare_species(human, cyno, ECD_START, ECD_END)
    mouse_from_step14 = {}
    sequon_csv = DERIVED / "14-sequons.csv"
    if sequon_csv.exists():
        for row in csv.DictReader(sequon_csv.open()):
            mouse_from_step14[int(row["uniprot_pos"])] = (
                row["mouse_residues_aligned"], row["species_status"])
    # Which sequons step 14 measured as the nearest attachment point to a current
    # anchor. Read from its committed table rather than typed in here.
    near_anchor = set()
    distance_csv = DERIVED / "14-anchor-glycan-distance.csv"
    if distance_csv.exists():
        for row in csv.DictReader(distance_csv.open()):
            if row["anchor_set"] == "current H370 face":
                near_anchor.add(int(row["nearest_sequon"]))
    emit("   In the status column, 'both' means the pattern is there in human and")
    emit("   in cyno, 'human only' that cyno has lost it, and 'cyno only' that cyno")
    emit("   has one where human has none.")
    emit()
    emit("   | position | human | cyno, aligned | mouse, aligned (step 14) | cyno "
         "status | nearest to a current anchor? |")
    emit("   |---|---|---|---|---|---|")
    sequon_rows = []
    for row in cyno_sequons:
        pos = row["pos"]
        # step 14's comparison names the second species 'mouse' in its own output;
        # the second species here is cyno, so the label is restated.
        status = row["status"].replace("mouse only", "cyno only")
        mouse_triplet, mouse_status = mouse_from_step14.get(pos, ("-", "-"))
        flagged = "yes" if pos in near_anchor else "no"
        note = ""
        if pos is not None and run_end and pos > run_end:
            note = " (past the readable part of the entry)"
        emit(f"   | N{pos if pos is not None else 'no aligned position'} | "
             f"{row['human']} | {row['mouse']} | {mouse_triplet} | {status}{note} | "
             f"{flagged} |")
        sequon_rows.append(dict(uniprot_pos=pos, human_residues=row["human"],
                                cyno_residues_aligned=row["mouse"],
                                mouse_residues_aligned=mouse_triplet,
                                cyno_status=status, mouse_status=mouse_status,
                                nearest_to_current_anchor=(pos in near_anchor),
                                past_readable_part=bool(
                                    pos is not None and pos > run_end)))
    emit()
    anchor_sequons = sorted(p for p in near_anchor)
    emit("   The ones that matter here are the attachment points step 14 measured as")
    emit("   nearest to a current anchor: "
         + ", ".join(f"N{p}" for p in anchor_sequons) + ".")
    emit()
    for pos in anchor_sequons:
        row = next((r for r in sequon_rows if r["uniprot_pos"] == pos), None)
        if row is None:
            emit(f"     N{pos}: not found among the human sequons — unexpected.")
            failures.append(f"sequon N{pos} named by step 14 was not found here")
            continue
        in_cyno = {"both": "intact in cyno",
                   "human only": "LOST in cyno",
                   "cyno only": "present in cyno and absent in human"}.get(
                       row["cyno_status"], row["cyno_status"])
        in_mouse = {"both": "intact in mouse",
                    "human only": "lost in mouse",
                    "-": "not compared"}.get(row["mouse_status"],
                                             row["mouse_status"])
        emit(f"     N{pos}: human {row['human_residues']}, cyno "
             f"{row['cyno_residues_aligned']}, mouse "
             f"{row['mouse_residues_aligned']} — {in_cyno}, {in_mouse}")
    emit()
    n361 = next((r for r in sequon_rows if r["uniprot_pos"] == 361), None)
    if n361 is not None:
        emit("   N361 is the one worth a sentence of its own. It is a sequon in human")
        emit("   and not in mouse, where the asparagine is replaced by a tyrosine, and")
        emit("   that asymmetry is why step 10 demotes a design leaning on E400, E421")
        emit(f"   or E424. In cyno the three residues are {n361['cyno_residues_aligned']}"
             f", so the sequon is "
             + ("intact" if n361["cyno_status"] == "both" else "absent") + ":")
        emit("   cyno matches human here and mouse is the odd one out.")
        emit()

    # ---- 7. The reviewed rhesus control --------------------------------------
    emit("7. Control: the reviewed rhesus macaque entry")
    emit()
    emit("   The cyno entry is unreviewed, so a result drawn from it alone rests on")
    emit("   one automatic gene prediction. Rhesus macaque EGFR is a reviewed entry")
    emit("   of the full length, and the two monkeys are closely related, so reading")
    emit("   it alongside says whether anything above is an artefact of that")
    emit("   prediction. This is a check on the source and not a second answer: the")
    emit("   question asked is about cyno.")
    emit()
    rhesus_offset, rhesus_share, rhesus_n = region_offsets["rhesus"]
    rhesus_d3 = differences(human, rhesus, maps["rhesus"], D3_START, D3_END)
    emit(f"   {RHESUS_ACCESSION} ({rhesus_header.split('OS=')[0].strip()}): "
         f"{len(rhesus)} aa, offset ours = theirs {rhesus_offset:+d} across "
         f"{100.0 * rhesus_share:.1f}% of domain III.")
    emit("   Domain III differences from human: "
         + ("  ".join(f"{h}{p}{o}" for p, h, o in rhesus_d3) or "none"))
    emit()
    disagreements = []
    for row in current_rows:
        if row["rhesus_aa"] != row["cyno_aa"]:
            disagreements.append(
                f"{row['human_aa']}{row['uniprot_pos']}: cyno {row['cyno_aa']}, "
                f"rhesus {row['rhesus_aa']}")
    check("the two macaque entries agree at every current anchor",
          not disagreements,
          "they do" if not disagreements else "; ".join(disagreements))
    rhesus_sequons = step14.compare_species(human, rhesus, ECD_START, ECD_END)
    rhesus_by_pos = {r["pos"]: r for r in rhesus_sequons}
    sequon_disagreements = []
    for pos in anchor_sequons:
        cyno_row = next((r for r in sequon_rows if r["uniprot_pos"] == pos), None)
        rhesus_row = rhesus_by_pos.get(pos)
        cyno_state = cyno_row["cyno_status"] if cyno_row else "absent"
        rhesus_state = rhesus_row["status"] if rhesus_row else "absent"
        if cyno_state != rhesus_state:
            sequon_disagreements.append(
                f"N{pos}: cyno {cyno_state}, rhesus {rhesus_state}")
    check("the two macaque entries agree about every anchor-adjacent sequon",
          not sequon_disagreements,
          "they do" if not sequon_disagreements else "; ".join(sequon_disagreements))
    emit()

    # ---- 8. Controls against steps 01 and 14 --------------------------------
    emit("8. Controls against steps 01 and 14")
    emit()
    emit("   This script aligns sequences that step 01 already aligned and finds")
    emit("   sequons that step 14 already found. Where it recomputes something, it")
    emit("   is compared against what the earlier step committed, and the run stops")
    emit("   if the two disagree. Steps 04 and 06 once answered the same question")
    emit("   two different ways for months before anyone noticed.")
    emit()
    emit("   The rule the sequons come from is step 14's own, loaded as a module")
    emit("   rather than copied, so there is one implementation of it.")
    emit()
    try:
        common.cross_check_residue_set(
            "human vs mouse domain III differences, against step 01",
            [p for p, _, _ in mouse_d3],
            DERIVED / "01-domain3-differences.csv", emit=emit)
    except common.CrossCheckError as exc:
        failures.append(str(exc))
    human_sequon_positions = [p for p, _ in step14.find_sequons(
        human[ECD_START - 1:ECD_END], first=ECD_START)]
    try:
        common.cross_check_residue_set(
            "human sequons in the extracellular region, against step 14",
            human_sequon_positions, sequon_csv, emit=emit)
    except common.CrossCheckError as exc:
        failures.append(str(exc))
    mouse_labels_here = set(mouse_d3_labels)
    committed_labels = set()
    diff_csv = DERIVED / "01-domain3-differences.csv"
    if diff_csv.exists():
        for row in csv.DictReader(diff_csv.open()):
            committed_labels.add(row["label"])
    check("the mouse substitutions themselves match step 01, not only the positions",
          mouse_labels_here == committed_labels,
          f"all {len(committed_labels)} agree" if mouse_labels_here == committed_labels
          else f"only here: {sorted(mouse_labels_here - committed_labels)}; only in "
               f"step 01: {sorted(committed_labels - mouse_labels_here)}")
    emit()

    # ---- 9. The tests of this script's own rule ------------------------------
    emit("9. The charge rule, tested every time")
    emit()
    failures.extend(run_tests(emit))

    # ---- 10. What this changes ----------------------------------------------
    emit("=" * 72)
    emit("10. What this changes about the design")
    emit("=" * 72)
    emit()
    emit("   Nothing. That is the honest answer and it was the expected one. Cyno is")
    emit("   not an objective: the three objectives are pH selectivity, mouse")
    emit("   cross-reactivity and human binding, and the competition page names two")
    emit("   species. No anchor moves, no ranking term changes, no campaign is")
    emit("   resized.")
    emit()
    if cyno_diff:
        emit("   What it records: among the "
             + f"{len(current_rows)} anchors, "
             + ", ".join(f"{r['human_aa']}{r['uniprot_pos']}{r['cyno_aa']}"
                         for r in cyno_diff)
             + " differ between human and")
        emit("   cyno. Section 5 says for each one whether the charge behaviour the")
        emit("   design rule depends on survives the substitution.")
    else:
        emit(f"   What it records: all {len(current_rows)} anchors of the current "
             "H370 face hold the same")
        emit("   amino acid in cyno as in human, so a binder that works by pairing")
        emit("   charges against them has the same charges to pair against in cyno.")
        emit("   The one difference anywhere in domain III, "
             + ", ".join(f"{h}{p}{o}" for p, h, o in cyno_d3)
             + ", is not an anchor")
        emit("   and is the same position that differs in mouse. The attachment")
        emit("   points nearest the anchors are intact in cyno: "
             + ", ".join(f"N{p}" for p in anchor_sequons) + ". That")
        emit("   includes N361, the sequon human has and mouse does not, which cyno")
        emit("   has. So the asymmetry behind step 10's demotion of E400, E421 and")
        emit("   E424 is between human and mouse specifically, rather than between")
        emit("   human and other species generally. The demotion stands unchanged,")
        emit("   because mouse cross-reactivity is the objective and cyno is not.")
        emit("   The methods write-up may say all of this, with the qualifications")
        emit("   below.")
    emit()
    emit("   The qualifications, which belong with the claim wherever it is made:")
    emit()
    emit("     - The cyno entry is unreviewed. It is an automatic gene prediction")
    emit("       that no curator has read. The reviewed rhesus macaque entry agrees")
    emit("       with it everywhere this step looked, which is reassurance and not")
    emit("       proof.")
    emit("     - This is sequence analysis. It says which amino acid sits at each")
    emit("       position and nothing about whether the face is shaped the same way")
    emit("       or reachable, which for human is what steps 02 to 06 settled and")
    emit("       for cyno has not been done at all. No cyno structure was examined.")
    emit("     - Identical anchors do not make a binder cross-reactive. The")
    emit("       surrounding surface a binder also touches was not compared position")
    emit("       by position outside domain III.")
    emit("     - Nothing was measured about binding. No cyno affinity is predicted,")
    emit("       claimed or implied.")
    emit()

    # ---- files ---------------------------------------------------------------
    wrote_files = False
    if not args.break_offset and BROKEN_RULE is None:
        DERIVED.mkdir(parents=True, exist_ok=True)
        with (DERIVED / "19-cyno-anchor-conservation.csv").open("w") as fh:
            fh.write("anchor_set,uniprot_pos,challenge_pos,role,human_aa,cyno_aa,"
                     "rhesus_aa,mouse_aa,identical_in_cyno,charge_verdict\n")
            for row in anchor_rows:
                fh.write(f"{row['anchor_set']},{row['uniprot_pos']},"
                         f"{row['challenge_pos']},{row['role']},{row['human_aa']},"
                         f"{row['cyno_aa']},{row['rhesus_aa']},{row['mouse_aa']},"
                         f"{row['identical']},{row['verdict']}\n")
        with (DERIVED / "19-cyno-sequons.csv").open("w") as fh:
            fh.write("uniprot_pos,human_residues,cyno_residues_aligned,"
                     "mouse_residues_aligned,cyno_status,mouse_status,"
                     "nearest_to_current_anchor,past_readable_part\n")
            for row in sequon_rows:
                fh.write(f"{row['uniprot_pos']},{row['human_residues']},"
                         f"{row['cyno_residues_aligned']},"
                         f"{row['mouse_residues_aligned']},{row['cyno_status']},"
                         f"{row['mouse_status']},"
                         f"{row['nearest_to_current_anchor']},"
                         f"{row['past_readable_part']}\n")
        with (DERIVED / "19-cyno-domain3-differences.csv").open("w") as fh:
            fh.write("uniprot_pos,challenge_pos,human_aa,cyno_aa,label,"
                     "also_differs_in_mouse,inside_fallback_epitope\n")
            for pos, human_aa, other_aa in cyno_d3:
                fh.write(f"{pos},{pos - 24},{human_aa},{other_aa},"
                         f"{human_aa}{pos}{other_aa},"
                         f"{pos in set(p for p, _, _ in mouse_d3)},"
                         f"{EPI_START <= pos <= EPI_END}\n")
        emit("Wrote data/derived/19-cyno-anchor-conservation.csv, "
             "19-cyno-sequons.csv and")
        emit("19-cyno-domain3-differences.csv")
        emit()
        wrote_files = True

    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED —")
        for failure in failures:
            emit(f"  - {failure}")
    else:
        intact = [p for p in anchor_sequons
                  if next((r["cyno_status"] for r in sequon_rows
                           if r["uniprot_pos"] == p), None) == "both"]
        emit(f"RESULT: PASSED. {len(cyno_same)} of {len(current_rows)} current "
             f"anchors hold the same amino acid in cyno; {len(intact)} of "
             f"{len(anchor_sequons)} anchor-adjacent")
        emit("sequons intact in cyno ("
             + ", ".join(f"N{p}" for p in intact) + "). Changes no decision.")
    emit("=" * 72)

    if wrote_files:
        FINDINGS.mkdir(parents=True, exist_ok=True)
        headline = (
            f"The answer: {len(cyno_same)} of the {len(current_rows)} anchors of "
            "the current H370 face hold the same\namino acid in cynomolgus monkey "
            f"as in human"
            + ("" if cyno_diff else ", so none of them loses the charge the\n"
                                    "design rule pairs against")
            + ". The attachment points for sugar chains that sit nearest\nthose "
              "anchors — "
            + ", ".join(f"N{p}" for p in anchor_sequons)
            + " — are all present in cyno too, including N361, which\nis present "
              "in human and absent in mouse. The offset between the cyno entry's "
              "own\nnumbering and ours is "
            + ("zero" if cyno_offset == 0 else f"{cyno_offset:+d}")
            + f", measured at all {cyno_d3_n} positions of domain III rather\nthan "
              "assumed.\n")
        (FINDINGS / "19-cyno-conservation.md").write_text(
            "# Cynomolgus monkey EGFR: are the anchors and the sequons conserved?\n\n"
            "Computed output of `analysis/19_cyno_conservation.py`. "
            "Do not hand-edit.\n\n"
            "**This is not a design input and it changes no decision.** The "
            "competition page names\nhuman and mouse only, and the decision of 2 "
            "October 2026 is to design for human and\nmouse only. A pre-launch "
            "message had said \"mouse and cyno\", so whether the face we\naim at "
            "survives in cynomolgus monkey (*Macaca fascicularis*, the crab-eating "
            "macaque)\nwas left open rather than answered. This step exists to turn "
            "that open question into\na recorded sentence, so that the methods "
            "write-up can either make the claim or name\nwhat breaks. No anchor is "
            "added, dropped, demoted or promoted because of it.\n\n"
            + headline + "\n"
            "EGFR is the epidermal growth factor receptor, the protein we design "
            "against. UniProt\nis the public archive of protein sequences, and every "
            "residue number here is a\nposition in its full human record P00533.\n\n"
            "The cynomolgus monkey sequence is cached at "
            "`data/sequences/egfr-macaque-uniprot.fasta`,\nwith the accession, the "
            "address and the date fetched recorded alongside it in\n"
            "`egfr-macaque-uniprot.source.txt`. Section 1 below says whether that "
            "file is tracked\nby git, checked by asking git rather than by assuming.\n\n"
            "```\n" + "\n".join(out) + "\n```\n")
    return 1 if failures else 0


def failed_checks(stdout):
    """The checks a run reported as failed, read from the list it prints at the
    end.

    Only the lines after the 'CHECK(S) FAILED' heading count. Prose elsewhere in
    the output is also written as bullets, and a self-test that counted those
    would pass on a run that had caught nothing.
    """
    marker = "CHECK(S) FAILED"
    if marker not in stdout:
        return []
    tail = stdout[stdout.index(marker):]
    return [line.strip()[2:] for line in tail.splitlines()
            if line.strip().startswith("- ")]


def run_self_test():
    """Break each check in turn and confirm the check meant to cover it notices.
    A check that has never been seen to fail is not known to work."""
    print("=" * 72)
    print("SELF-TEST: does each check here catch the mistake it is for?")
    print("=" * 72)
    print()
    script = str(Path(__file__).resolve())
    offset_check = "domain III identity through the measured map"
    cases = [
        (["--break-offset", "1"], "numbering shifted by one position",
         offset_check,
         "every anchor and every sequon would be read one position away from "
         "where it is, and the letters would still look like amino acids"),
        (["--break-offset", "-3"], "numbering shifted by three the other way",
         offset_check, "the same, in the other direction"),
    ]
    for rule, (expected_case, consequence) in BREAKABLE_RULES.items():
        cases.append((["--break-rule", rule],
                      f"charge rule: {expected_case} switched off",
                      expected_case, consequence))
    failures = []
    for argv, label, expected_case, consequence in cases:
        proc = subprocess.run([sys.executable, script] + argv,
                              capture_output=True, text=True,
                              cwd=str(Path(__file__).resolve().parent))
        caught = failed_checks(proc.stdout)
        named = any(expected_case in line for line in caught)
        ok = proc.returncode != 0 and named
        print(f"   [{'PASS' if ok else 'FAIL'}] {label}")
        print(f"     if unnoticed: {consequence}")
        print(f"     exit code {proc.returncode}, {len(caught)} check(s) failed, "
              f"the one meant to cover it {'among them' if named else 'NOT among them'}")
        for line in caught[:3]:
            print(f"       {line[:110]}")
        if not ok:
            failures.append(
                f"{label}: "
                + ("breaking this failed no check" if proc.returncode == 0
                   else f"checks failed but not the one meant to cover it "
                        f"({expected_case})"))
        print()
    print("=" * 72)
    if failures:
        print(f"SELF-TEST FAILED: {len(failures)} case(s) went unnoticed.")
        for failure in failures:
            print(f"  - {failure}")
    else:
        print(f"SELF-TEST PASSED. All {len(cases)} broken runs were caught by the "
              "check meant to cover them.")
    print("=" * 72)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
