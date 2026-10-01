"""Did the design run leave the target's residue numbering alone?

BindCraft2 registers a filter metric called Target_Crop_Length, which suggests it
may crop the target. Its source says the metric counts residues that are not
padding, and padding is a batching device, so the evidence read so far points at no
cropping and no renumbering. That is a reading of source code, not a measurement of
what comes out, and the cost of being wrong is high and silent: analysis/10 turns
the numbers in a returned structure back into ours, and a renumbered target would
produce confident verdicts about the wrong residues.

This module answers the question by comparing the target chain in a returned
complex against the structure file that went in. It has no dependencies beyond the
standard library, because it runs in two places: on the rented machine during the
smoke run, where Biopython is not guaranteed to be installed, and locally inside
analysis/10, which refuses to score a candidate it cannot reconcile. One
implementation serves both, so the two cannot disagree.

Abbreviations: mmCIF is the structure file format BindCraft2 writes; PDB is the
older format the trimmed target is supplied in (and also the Protein Data Bank, the
public archive of measured structures). `auth_seq_id` is the mmCIF field holding the
residue number the author of the file gave, which is the field Biopython reads as
the residue number, and `label_seq_id` is the field holding a renumbering from 1 that
the format requires. Which of the two BindCraft2 fills with the input's numbers is
exactly what this checks rather than assumes.
"""

from __future__ import annotations

import shlex
from pathlib import Path

STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}

# Outcomes of a comparison. The first two are reconciled: every residue in the
# returned target sits at a number where the input has the same amino acid, so a
# position read from the returned file means the same thing in the input's
# numbering. The last two are not.
IDENTICAL = "identical"
CROPPED = "cropped, numbers preserved"
RENUMBERED = "renumbered"
MISMATCH = "does not match the input"
RECONCILED = (IDENTICAL, CROPPED)


def _read_pdb(path):
    chains = {}
    last = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        name = line[17:20].strip()
        if name not in STANDARD:
            continue
        chain = line[21].strip() or "?"
        try:
            number = int(line[22:26])
        except ValueError:
            continue
        key = (number, line[26].strip())
        if last.get(chain) == key:
            continue
        last[chain] = key
        chains.setdefault(chain, []).append(
            dict(number=number, insertion=key[1], name=name, label_seq=None))
    return chains


def _cif_loop_rows(text):
    """Rows of the `_atom_site` loop as dicts keyed by column name."""
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() == "loop_":
            columns = []
            cursor = index + 1
            while cursor < len(lines) and lines[cursor].startswith("_"):
                columns.append(lines[cursor].strip())
                cursor += 1
            if columns and columns[0].startswith("_atom_site."):
                tokens = []
                while cursor < len(lines):
                    row = lines[cursor]
                    if (row.startswith("_") or row.startswith("loop_")
                            or row.startswith("#") or row.startswith("data_")):
                        break
                    if row.strip():
                        tokens.extend(shlex.split(row, posix=True))
                    cursor += 1
                width = len(columns)
                for start in range(0, len(tokens) - width + 1, width):
                    yield dict(zip(columns, tokens[start:start + width]))
                return
            index = cursor
        else:
            index += 1


def _read_cif(path):
    chains = {}
    last = {}
    for row in _cif_loop_rows(Path(path).read_text()):
        if row.get("_atom_site.group_PDB", "ATOM") != "ATOM":
            continue
        name = row.get("_atom_site.auth_comp_id",
                       row.get("_atom_site.label_comp_id", ""))
        if name not in STANDARD:
            continue
        chain = row.get("_atom_site.auth_asym_id",
                        row.get("_atom_site.label_asym_id", "?"))
        number_text = row.get("_atom_site.auth_seq_id",
                              row.get("_atom_site.label_seq_id"))
        try:
            number = int(number_text)
        except (TypeError, ValueError):
            continue
        insertion = row.get("_atom_site.pdbx_PDB_ins_code", "?")
        insertion = "" if insertion in ("?", ".") else insertion
        key = (number, insertion)
        if last.get(chain) == key:
            continue
        last[chain] = key
        label = row.get("_atom_site.label_seq_id")
        try:
            label = int(label)
        except (TypeError, ValueError):
            label = None
        chains.setdefault(chain, []).append(
            dict(number=number, insertion=insertion, name=name, label_seq=label))
    return chains


def read_residues(path):
    """{chain id: [residue dicts in file order]} for the standard amino acids.

    Each residue carries `number` (the number a reader such as Biopython reports),
    `insertion`, `name` (three letters) and `label_seq` (mmCIF only, else None).
    """
    path = Path(path)
    if path.suffix.lower() == ".cif":
        return _read_cif(path)
    return _read_pdb(path)


def complex_size(chains):
    """Total residues across every chain: the figure the run time scales with."""
    return sum(len(residues) for residues in chains.values())


def _pairs(residues):
    return [(r["number"], r["name"]) for r in residues]


def _agreement(output_residues, input_lookup):
    if not output_residues:
        return 0.0
    hits = sum(1 for r in output_residues
               if input_lookup.get(r["number"]) == r["name"])
    return hits / len(output_residues)


def find_target(output_chains, input_chains):
    """Which returned chain is the target, and which input chain it came from.

    By measurement rather than by letter: the pair (output chain, input chain)
    whose residues agree best by number and amino acid. A designed binder agrees
    with nothing, because its sequence did not exist before the run.
    Returns (output chain, input chain, share), or (None, None, 0.0).
    """
    best = (None, None, 0.0)
    for out_id, out_residues in output_chains.items():
        for in_id, in_residues in input_chains.items():
            lookup = dict(_pairs(in_residues))
            share = _agreement(out_residues, lookup)
            if share > best[2]:
                best = (out_id, in_id, share)
    return best


def _ranges(numbers):
    numbers = sorted(numbers)
    if not numbers:
        return "none"
    runs, start, previous = [], numbers[0], numbers[0]
    for n in numbers[1:]:
        if n == previous + 1:
            previous = n
            continue
        runs.append((start, previous))
        start = previous = n
    runs.append((start, previous))
    return ", ".join(f"{a}" if a == b else f"{a}-{b}" for a, b in runs)


def compare(input_residues, output_residues):
    """Compare one returned target chain against the input chain it came from.

    Returns a dict. `status` is one of IDENTICAL, CROPPED, RENUMBERED, MISMATCH;
    `ok` is whether it is reconciled; `detail` says exactly how it differs.

    RENUMBERED is the case that must not be scored: the amino acid sequence is
    present in the input but at different numbers, so a verdict computed on the
    returned numbers would be about the wrong residues. `recoverable_from` says
    whether the returned file carries anything that would let the numbers be mapped
    back, and is reported rather than acted on.
    """
    in_pairs, out_pairs = _pairs(input_residues), _pairs(output_residues)
    in_lookup = dict(in_pairs)
    in_numbers, out_numbers = set(in_lookup), {n for n, _ in out_pairs}

    same_count = len(in_pairs) == len(out_pairs)
    missing = sorted(in_numbers - out_numbers)
    extra = sorted(out_numbers - in_numbers)
    wrong_identity = [n for n, name in out_pairs
                      if n in in_lookup and in_lookup[n] != name]

    result = dict(
        input_count=len(in_pairs), output_count=len(out_pairs),
        same_count=same_count,
        same_numbering=(out_pairs == in_pairs),
        input_range=(min(in_numbers), max(in_numbers)) if in_numbers else None,
        output_range=(min(out_numbers), max(out_numbers)) if out_numbers else None,
        missing_from_output=missing, not_in_input=extra,
        wrong_identity=wrong_identity, recoverable_from="nothing found",
    )

    if out_pairs == in_pairs:
        result.update(status=IDENTICAL, ok=True,
                      detail=f"all {len(in_pairs)} residues present at the "
                             f"input's own numbers, same amino acids")
        return result

    if not extra and not wrong_identity:
        result.update(
            status=CROPPED, ok=True,
            detail=f"{len(out_pairs)} of the input's {len(in_pairs)} residues "
                   f"returned, every one at the input's own number with the same "
                   f"amino acid; absent from the output: {_ranges(missing)}")
        return result

    # Not reconciled. Say whether the sequence survived under different numbers.
    in_seq = "".join(ONE[name] for _, name in in_pairs)
    out_seq = "".join(ONE[name] for _, name in out_pairs)
    where = in_seq.find(out_seq)
    unique = where >= 0 and in_seq.count(out_seq) == 1
    if unique:
        implied = [in_pairs[where + i][0] for i in range(len(out_pairs))]
        offsets = {a - b for a, b in zip(implied, [n for n, _ in out_pairs])}
        label_ok = all(r["label_seq"] is not None for r in output_residues)
        label_matches = label_ok and [
            r["label_seq"] for r in output_residues] == implied
        how = []
        if len(offsets) == 1:
            how.append(f"a constant offset of {next(iter(offsets)):+d}")
        else:
            how.append("a non-constant mapping")
        if label_matches:
            how.append("and label_seq_id equals the input numbers")
        elif label_ok:
            labels = [r["label_seq"] for r in output_residues]
            how.append(f"label_seq_id runs {labels[0]}-{labels[-1]}, which "
                       f"is not the input numbering")
        result["recoverable_from"] = (
            "the amino acid sequence matches input residues "
            f"{implied[0]}-{implied[-1]} exactly, so the mapping can be "
            f"recovered by sequence match ({'; '.join(how)}), but nothing the "
            f"file states says so")
        result.update(
            status=RENUMBERED, ok=False,
            detail=f"the returned target's {len(out_pairs)} residues are numbered "
                   f"{result['output_range'][0]}-{result['output_range'][1]} but "
                   f"read as input residues {implied[0]}-{implied[-1]}, so a number "
                   f"read from the returned file is not the input's number for "
                   f"that residue")
        return result

    result.update(
        status=MISMATCH, ok=False,
        detail=f"{len(wrong_identity)} returned residue(s) sit at a number where "
               f"the input has a different amino acid"
               + (f", {len(extra)} number(s) not in the input at all"
                  if extra else "")
               + (f", counts {len(in_pairs)} in against {len(out_pairs)} out"
                  if not same_count else ""))
    return result


def check_complex(input_path, output_path):
    """Compare a returned complex's target chain against the input structure.

    Returns the dict from `compare`, plus `target_chain`, `input_chain`,
    `chain_share` and `complex_residues`; or a dict with status MISMATCH and
    ok False when no returned chain resembles the input at all.
    """
    input_chains = read_residues(input_path)
    output_chains = read_residues(output_path)
    out_id, in_id, share = find_target(output_chains, input_chains)
    size = complex_size(output_chains)
    if out_id is None or share < 0.5:
        # Nothing agrees by number. Try the same match ignoring numbers, so the
        # report can say a renumbered target was found rather than no target.
        in_seqs = {c: "".join(ONE[r["name"]] for r in rs)
                   for c, rs in input_chains.items()}
        for cid, residues in output_chains.items():
            seq = "".join(ONE[r["name"]] for r in residues)
            for icid, in_seq in in_seqs.items():
                if len(seq) >= 20 and seq in in_seq:
                    result = compare(input_chains[icid], residues)
                    result.update(target_chain=cid, input_chain=icid,
                                  chain_share=share, complex_residues=size)
                    return result
        return dict(status=MISMATCH, ok=False, target_chain=None,
                    input_chain=None, chain_share=share, complex_residues=size,
                    detail="no returned chain resembles the input target by "
                           "number or by sequence")
    result = compare(input_chains[in_id], output_chains[out_id])
    result.update(target_chain=out_id, input_chain=in_id, chain_share=share,
                  complex_residues=size)
    return result


def describe(result):
    """One plain-language sentence for a comparison result."""
    text = f"target numbering {result['status']}: {result['detail']}"
    if not result["ok"] and result.get("recoverable_from"):
        text += f". Recoverable? {result['recoverable_from']}"
    return text


def _self_test(input_path):
    """Exercise every outcome against a real structure, in both file formats.

    Writes altered copies of `input_path` and checks each is classified as it
    should be. The mmCIF copies are written by Biopython, so this also confirms
    that this module's reader and Biopython's agree on which number is the
    residue's number, which is the property analysis/10 depends on.
    """
    import tempfile
    from Bio.PDB import PDBParser
    from Bio.PDB.mmcifio import MMCIFIO

    structure = PDBParser(QUIET=True).get_structure("t", str(input_path))
    residues = [r for r in structure.get_residues()]
    failures = []

    def expect(label, got, want):
        ok = got == want
        print(f"   [{'PASS' if ok else 'FAIL'}] {label}: {got}"
              + ("" if ok else f", expected {want}"))
        if not ok:
            failures.append(label)

    with tempfile.TemporaryDirectory() as tmp:
        def write(name, mutate):
            fresh = PDBParser(QUIET=True).get_structure("t", str(input_path))
            chain = next(iter(next(iter(fresh))))
            chain = next(iter(fresh[0]))
            mutate(chain)
            io = MMCIFIO()
            io.set_structure(fresh)
            path = Path(tmp) / name
            io.save(str(path))
            return path

        def untouched(chain):
            pass

        def crop(chain):
            for r in list(chain)[:20]:
                chain.detach_child(r.id)

        def shift(chain):
            # Renumber from 1, which is what a run that cropped and renumbered
            # would write. Two passes so no number collides mid-way.
            items = list(chain)
            for r in items:
                chain.detach_child(r.id)
            for i, r in enumerate(items, start=1):
                r.id = (" ", i, " ")
                chain.add(r)

        def swap_identity(chain):
            list(chain)[5].resname = "TRP"

        print("\n   target_numbering self-test, against", Path(input_path).name)
        for name, mutate, want in [
                ("same.cif", untouched, IDENTICAL),
                ("cropped.cif", crop, CROPPED),
                ("renumbered.cif", shift, RENUMBERED),
                ("changed.cif", swap_identity, MISMATCH)]:
            out = write(name, mutate)
            result = check_complex(input_path, out)
            expect(name, result["status"], want)
            if want == RENUMBERED:
                expect("renumbered: not ok", result["ok"], False)
                print("        ", describe(result))
        expect("a PDB input compared with itself",
               check_complex(input_path, input_path)["status"], IDENTICAL)
    print(f"\n   {len(residues)} residues in the reference. "
          f"{'ALL PASSED' if not failures else str(len(failures)) + ' FAILED'}")
    return 1 if failures else 0


if __name__ == "__main__":
    import sys
    raise SystemExit(_self_test(
        sys.argv[1] if len(sys.argv) > 1 else
        Path(__file__).resolve().parents[1] / "data" / "structures"
        / "6aru_domain3.pdb"))
