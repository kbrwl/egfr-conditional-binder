"""24_structural_novelty.py -- score each submitted design on the organisers' novelty scale.

The repository's earlier novelty check was sequence-only, against Swiss-Prot, and
the methods document named the structural half as unrun. This runs it: each
submitted design's predicted structure searched against the Protein Data Bank
with FoldSeek, aligned with TM-align, and scored against the scale the
organisers publish at adaptyvbio.com/blog/novelty.

Abbreviations, expanded because this file gets read on its own:
  PDB       the Protein Data Bank, the public archive of measured
            three-dimensional protein structures.
  FoldSeek  a search tool that finds structures resembling a query structure,
            the structural counterpart of a sequence search.
  TM-score  a measure of how alike two protein folds are, from 0 to 1, where
            about 0.5 is the conventional line above which two structures are
            usually taken to share a fold.
  aa        amino acids, the unit a protein's length is counted in.

THE SCALE, AS THE ORGANISERS PUBLISH IT
----------------------------------------
For non-antibodies, which is what we submit:

  level 4   sequence similarity 30% or below AND less than moderate structural
  level 3   exactly one of: sequence above 30%, or moderate structural
  level 2   sequence above 70%; or high structural; or sequence above 30% AND
            moderate structural
  level 1   sequence above 70% AND moderate structural

"Moderate structural" means more than 70% of the sequence covered by domains
matching a known structure at a TM-score of 0.5 or better. "High" means the same
coverage at 0.8 or better. **The gate is level 3**, so 3 and 4 both pass and 2
does not.

HOW THE SEQUENCE HALF IS READ, WHICH IS THE CONSEQUENTIAL JUDGEMENT HERE
-------------------------------------------------------------------------
Our MMseqs2 search against Swiss-Prot returned best hits at 31-42% raw identity,
which looks like it crosses the 30% line. It does not, and the distinction
decides whether these designs read as level 3 or level 2.

Those identities are over short partial alignments and every one of them is
statistically insignificant: the best e-value across 2,239 hits was 0.30, and
none reached 0.001. An e-value of 0.30 means a hit that good is expected by
chance in a database that size. A 48-residue query routinely finds 35-40%
identity over 40-odd residues of an unrelated protein, and that is not homology.
The organisers' own wording is that a design is de novo if it hits nothing "with
any homology", which is a question about significance rather than about a raw
percentage.

FoldSeek's own sequence identities over these structural alignments agree and
are reported below: they run far under 30%.

So the sequence half is read as **at or below 30%**. This is a judgement, it is
the single judgement this score rests on, and it is recorded here rather than
buried. If it were read the other way, every design would fall to level 2.

**The platform scores novelty automatically on upload and its score is the one
that counts.** This is our own reading of a published rule, computed so that an
unpleasant surprise at upload is less likely, not a substitute for it.
"""
import argparse
import collections
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HITS = ROOT / "data" / "derived" / "24-foldseek-pdb-hits.m8"
SUMMARY = ROOT / "data" / "derived" / "24-structural-novelty.csv"
FINDINGS = ROOT / "results" / "findings" / "24-structural-novelty.md"

# The command that produced the hits file, recorded because the PDB database is
# 2.2 GB and is not kept in the repository, so this script scores a search it
# does not itself run.
SEARCH_COMMAND = (
    "foldseek easy-search <binder monomers, PDB format> <PDB database> "
    "hits.m8 tmp --alignment-type 1 "
    "--format-output query,target,fident,alnlen,qcov,tcov,evalue,bits,"
    "alntmscore,qtmscore,ttmscore -e 10 --max-seqs 2000"
)
FIELDS = ["query", "target", "fident", "alnlen", "qcov", "tcov", "evalue",
          "bits", "alntmscore", "qtmscore", "ttmscore"]

MODERATE_TM = 0.5
HIGH_TM = 0.8
COVERAGE = 0.70
# Read within this of a threshold, a value is reported as near the boundary
# rather than as the label it was given.
BOUNDARY = 0.05

# Our sequence half, established in docs/decisions-log.md on 3 October 2026 and
# argued in this file's header.
SEQUENCE_OVER_30 = False


def level(seq_over_30, structural):
    """The organisers' scale. Returns (level, the reason in words)."""
    if structural == "high":
        return 2, "high structural similarity"
    if seq_over_30 and structural == "moderate":
        return 2, "sequence above 30% and moderate structural"
    if seq_over_30 or structural == "moderate":
        which = "sequence above 30%" if seq_over_30 else "moderate structural"
        return 3, f"exactly one condition met: {which}"
    return 4, "sequence at or below 30% and less than moderate structural"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--hits", type=Path, default=HITS)
    args = parser.parse_args()
    if not args.hits.exists():
        sys.exit(f"missing {args.hits}; run the FoldSeek search first:\n"
                 f"  {SEARCH_COMMAND}")

    rows = []
    for line in args.hits.read_text().splitlines():
        if not line.strip():
            continue
        rows.append(dict(zip(FIELDS, line.split("\t"))))

    by_query = collections.defaultdict(list)
    for r in rows:
        by_query[r["query"]].append(r)

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    emit("=" * 72)
    emit("STRUCTURAL NOVELTY — each submitted design against the PDB")
    emit("=" * 72)
    emit()
    emit(f"{len(rows)} hits across {len(by_query)} designs, from:")
    emit(f"  {SEARCH_COMMAND}")
    emit()
    emit("The structure searched is the binder alone, not the complex. A "
         "complex would")
    emit("match EGFR trivially and tell us nothing about whether our own "
         "design")
    emit("resembles something already known.")
    emit()
    emit("1. Best structural match per design, by TM-score over the query")
    emit()
    emit("   | design | closest PDB entry | query coverage | TM-score | "
         "seq identity | structural verdict |")
    emit("   |---|---|---|---|---|---|")

    results = []
    for query in sorted(by_query):
        hits = by_query[query]
        # The rule is about coverage AND fold similarity together, so the hit
        # that matters is the best TM-score among hits that actually cover the
        # design, not the best TM-score overall.
        covering = [h for h in hits if float(h["qcov"]) >= COVERAGE]
        pool = covering or hits
        best = max(pool, key=lambda h: float(h["qtmscore"]))
        tm = float(best["qtmscore"])
        qcov = float(best["qcov"])
        if qcov >= COVERAGE and tm >= HIGH_TM:
            structural = "high"
        elif qcov >= COVERAGE and tm >= MODERATE_TM:
            structural = "moderate"
        else:
            structural = "below moderate"
        lvl, why = level(SEQUENCE_OVER_30, structural)
        near = ""
        if abs(tm - HIGH_TM) <= BOUNDARY:
            near = " (near the 0.8 line)"
        elif abs(tm - MODERATE_TM) <= BOUNDARY:
            near = " (near the 0.5 line)"
        emit(f"   | {query} | {best['target']} | {qcov:.2f} | {tm:.3f}{near} "
             f"| {float(best['fident']) * 100:.1f}% | {structural} |")
        results.append(dict(design=query, target=best["target"], qcov=qcov,
                            tm=tm, fident=float(best["fident"]),
                            structural=structural, level=lvl, reason=why,
                            hits=len(hits), near_boundary=bool(near)))
    emit()

    emit("2. Novelty level")
    emit()
    emit("   Sequence half read as at or below 30%: no significant Swiss-Prot "
         "homology")
    emit("   (best e-value 0.30 over 2,239 hits, none below 0.001). The "
         "reasoning and")
    emit("   what turns on it are in this script's header.")
    emit()
    emit("   | design | level | why |")
    emit("   |---|---|---|")
    for r in results:
        emit(f"   | {r['design']} | **{r['level']}** | {r['reason']} |")
    emit()

    levels = collections.Counter(r["level"] for r in results)
    failing = [r for r in results if r["level"] < 3]
    emit("3. What this changes about the submission")
    emit()
    emit("   " + ", ".join(f"level {k}: {v} design(s)"
                           for k, v in sorted(levels.items())))
    emit()
    if failing:
        emit("   **These designs do not clear the gate and must be removed:**")
        for r in failing:
            emit(f"     {r['design']} — {r['reason']}")
    else:
        emit("   Every design clears the level 3 gate. None is removed on "
             "novelty.")
    emit()
    boundary = [r for r in results if r["near_boundary"]]
    if boundary:
        emit("   Reported as near a classification boundary rather than as the "
             "label given:")
        for r in boundary:
            emit(f"     {r['design']} at TM {r['tm']:.3f}. A small change in "
                 "the alignment would")
            emit("     move it across, and the platform's own run may land the "
                 "other side.")
        emit()
    emit("   This is our reading of a published rule. The platform scores "
         "novelty")
    emit("   automatically on upload and that score is the one that counts.")
    emit()
    emit("=" * 72)
    emit(f"RESULT: {len(results)} designs scored; "
         + ", ".join(f"{v} at level {k}" for k, v in sorted(levels.items()))
         + ("; ALL CLEAR THE GATE." if not failing
            else f"; {len(failing)} BELOW THE GATE."))
    emit("=" * 72)

    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["design", "closest_pdb_entry", "query_coverage",
                    "tm_score_over_query", "sequence_identity",
                    "structural_verdict", "novelty_level", "reason",
                    "hits_considered", "near_boundary"])
        for r in results:
            w.writerow([r["design"], r["target"], f"{r['qcov']:.3f}",
                        f"{r['tm']:.4f}", f"{r['fident']:.4f}",
                        r["structural"], r["level"], r["reason"], r["hits"],
                        "yes" if r["near_boundary"] else "no"])
    FINDINGS.parent.mkdir(parents=True, exist_ok=True)
    FINDINGS.write_text(
        "# Structural novelty of the submitted designs\n\n"
        "Computed output of `analysis/24_structural_novelty.py`. "
        "Do not hand-edit.\n\n```\n" + "\n".join(out) + "\n```\n")
    print(f"\nWrote {SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {FINDINGS.relative_to(ROOT)}")
    return 1 if failing else 0


if __name__ == "__main__":
    sys.exit(main())
