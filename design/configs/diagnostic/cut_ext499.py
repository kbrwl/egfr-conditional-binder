"""Cut a diagnostic 310-499 fragment from 6ARU, the same way analysis/09 cuts 310-480.

Why this is a scratchpad script and not a numbered analysis step: it produces one
diagnostic input for one test run, to answer whether the extended fragment clears
BindCraft2's target-confidence gate on its own. If the answer says the boundary
should move, that is a decision for the owner and the move belongs in analysis/09
behind its own constant, as docs/decisions-log.md already specifies. This script
does not change any committed boundary.

The method is copied from analysis/09_trim_target.py rather than reinvented: same
source file (the receptor chain alone), same RangeSelect by structure-file residue
number, same writer. The check at the end is the point -- over the residues the two
fragments share, this cut must reproduce 6aru_domain3.pdb exactly. If it does not,
the two were not made the same way and the comparison between them would be
measuring the method rather than the boundary.
"""

from pathlib import Path

from Bio.PDB import PDBParser, PDBIO, Select

# Resolved from this file's own location rather than hard-coded, so the script
# runs wherever the repository is checked out.
ROOT = Path(__file__).resolve().parents[3]
RECEPTOR = ROOT / "data" / "structures" / "6aru_receptor_only.pdb"
EXISTING = ROOT / "data" / "structures" / "6aru_domain3.pdb"
OUT = ROOT / "data" / "structures" / "6aru_domain3_ext499.pdb"

# Our numbering is the full UniProt record; 6ARU numbers by the mature protein, so
# UniProt = PDB + 24. Verified empirically by analysis/02 across 100% of aligned
# residues; never assumed, in either direction.
OFFSET = 24
TRIM_START, TRIM_END = 310, 499          # the extended boundary under test
EXISTING_END = 480                        # the committed boundary, for the overlap check


class RangeSelect(Select):
    """Keep exactly the residues whose structure-file numbers are listed.

    Copied from analysis/09 so the two fragments are cut identically.
    """

    def __init__(self, chain_id, keep_resnums):
        self.chain_id = chain_id
        self.keep = set(keep_resnums)

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[1] in self.keep


def main():
    structure = PDBParser(QUIET=True).get_structure("receptor", str(RECEPTOR))
    model = structure[0]
    chain = next(iter(model))
    keep_pdb = [p - OFFSET for p in range(TRIM_START, TRIM_END + 1)]

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(OUT), select=RangeSelect(chain.id, keep_pdb))

    written = PDBParser(QUIET=True).get_structure("ext", str(OUT))[0]
    written_chain = next(iter(written))
    residues = [r for r in written_chain if r.id[0] == " "]
    numbers = sorted(r.id[1] for r in residues)
    print(f"wrote {OUT.name}: {len(residues)} residues, "
          f"file numbers {numbers[0]}-{numbers[-1]}, "
          f"our numbering {numbers[0] + OFFSET}-{numbers[-1] + OFFSET}")

    # The overlap check. Over the residues both fragments contain, the atom lines
    # must be byte-identical, or the two cuts differ by method and not only by
    # boundary.
    def atom_lines(path, upper_pdb):
        out = []
        for line in Path(path).read_text().splitlines():
            if line.startswith("ATOM"):
                try:
                    number = int(line[22:26])
                except ValueError:
                    continue
                if number <= upper_pdb:
                    out.append(line)
        return out

    mine = atom_lines(OUT, EXISTING_END - OFFSET)
    theirs = atom_lines(EXISTING, EXISTING_END - OFFSET)
    if mine == theirs:
        print(f"overlap check PASSED: {len(mine)} atom lines over our "
              f"{TRIM_START}-{EXISTING_END} are byte-identical to "
              f"{EXISTING.name}, so the two cuts differ only in where they stop.")
    else:
        print(f"overlap check FAILED: {len(mine)} lines here against "
              f"{len(theirs)} in {EXISTING.name}. The two fragments were not cut "
              f"the same way, so comparing them would measure the method.")
        for a, b in zip(mine, theirs):
            if a != b:
                print(f"  first difference:\n    here  {a}\n    there {b}")
                break
        return 1

    # Does the extension actually include the disulfide partner the cut severed?
    # C470 is inside both fragments; C499 is the partner, outside 310-480.
    present = {r.id[1] + OFFSET for r in residues}
    for position in (470, 499):
        print(f"  residue {position} present: {position in present}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
