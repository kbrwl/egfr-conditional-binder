#!/usr/bin/env python3
"""
17_occupancy_table.py — how much of the binder on the chip is occupied, at the
assay's top concentration, for a given affinity.

Abbreviations, expanded here because this file gets read on its own:
  K_D   the dissociation constant, an affinity measured as a concentration.
        Smaller means a tighter grip. It is defined as the concentration of
        the other molecule at which half of the binding sites are occupied.
  nM    nanomolar, a billionth of a mole per litre: a unit of concentration.
  µM    micromolar, a millionth of a mole per litre: 1000 times looser than
        1 nM for the same number.

WHY THIS STEP EXISTS
---------------------
`docs/explainers/07-the-assay-and-what-it-changes.md`, `docs/decisions-log.md`
and `CLAUDE.md` each state, in prose, what share of the binder is occupied at
the assay's top concentration for a few illustrative K_D values -- for example
that a K_D of 10 micromolar leaves only about 9% of the binder occupied at the
assay's top concentration of 1000 nanomolar. Those figures are arithmetic, not
measurement, but CLAUDE.md's rule that a number which matters should be
regenerable by a script applies to arithmetic as much as to anything measured:
three documents were quoting the same small calculation typed out by hand.

THE FORMULA
-----------
K_D is, by definition, the analyte concentration at which half the binding
sites are occupied. The fraction occupied at any analyte concentration follows
from that definition alone:

    occupied fraction = [analyte] / ([analyte] + K_D)

Set [analyte] = K_D and the fraction is exactly 0.5, which is what makes K_D
the number it is. The function living in `egfr_common.occupancy_fraction` is
this formula and nothing more.

WHAT THIS CHANGES ABOUT THE DESIGN
-----------------------------------
The organisers flow the target at a top concentration of 1000 nM
(`docs/competition-qa-log.md`). A design in the high nanomolar range for human
binding at pH 6.5 sits near full occupancy and gives a strong, easily fitted
signal. A design deliberately weakened toward the far end of the instrument's
reportable range (about 0.1 nM to 10 micromolar) gives a progressively smaller
one, and far enough along it reads as no detectable binding at pH 6.5 as well
as at pH 7.4 -- which fails human binding rather than demonstrating pH
selectivity. This is the arithmetic behind the rule in `CLAUDE.md` to maximise
the pH gap rather than hold affinity down deliberately (`docs/decisions-
log.md`, Ruled out).

What this does not establish: where exactly the instrument stops returning a
fittable curve. That depends on the size of the design and the density of
binder on the chip as well as on occupancy, and the organisers have not said.
The table below gives the occupancy the formula implies; it does not give the
detection floor.

HOW THIS IS CHECKED
--------------------
Three properties any correct occupancy formula must have, asserted every run:
occupancy at K_D = 0 is 1 (an infinitely tight binder occupies everything),
occupancy at K_D equal to the analyte concentration is exactly 0.5 (the
definition of K_D), and occupancy strictly decreases as K_D grows (a weaker
binder occupies less, never more). `--break-rule flip` swaps the formula for
its reciprocal-shaped inverse, K_D / (K_D + analyte), to confirm the checks
actually fail when the formula is wrong rather than passing by construction.

Writes:
  results/findings/17-occupancy.md
  data/derived/17-occupancy.csv

Run standalone:  python analysis/17_occupancy_table.py
                 python analysis/17_occupancy_table.py --break-rule flip
                 python analysis/17_occupancy_table.py --self-test
"""

import argparse
import subprocess
import sys
from pathlib import Path

import egfr_common as common

DERIVED = common.DERIVED
FINDINGS = common.FINDINGS

# Set by --break-rule, which swaps in a wrong formula on purpose so the checks
# below can be confirmed to catch it. A check that has never been seen to fail
# is not known to work.
BROKEN_RULE = None


def occupancy(kd_nm, analyte_nm=common.ASSAY_TOP_CONCENTRATION_NM):
    """The formula under test, via the shared function -- or its deliberately
    wrong inverse when `--break-rule flip` is in effect."""
    if BROKEN_RULE == "flip":
        return kd_nm / (kd_nm + analyte_nm)
    return common.occupancy_fraction(kd_nm, analyte_nm)


# The K_D values the explainer and the decisions log illustrate with, in
# nanomolar: 100 nM, 1 micromolar, 3 micromolar, 10 micromolar. Chosen because
# they are the values already quoted in prose and span the high-nanomolar
# range down to the far end of the instrument's reportable range (about
# 0.1 nM to 10 micromolar, docs/competition-qa-log.md).
ILLUSTRATIVE_KD_NM = [100.0, 1_000.0, 3_000.0, 10_000.0]


def format_kd(kd_nm):
    """A K_D in nanomolar, written the way the prose documents write it."""
    if kd_nm < 1000:
        return f"{kd_nm:.0f} nM"
    return f"{kd_nm / 1000.0:.0f} µM"


def run_checks(emit):
    """The three properties any correct occupancy formula must have."""
    failures = []
    emit("   Checks on the formula itself, every run:")
    emit()

    zero = occupancy(0.0)
    ok = abs(zero - 1.0) < 1e-9
    emit(f"   [{'PASS' if ok else 'FAIL'}] K_D = 0 (an infinitely tight "
         f"binder): occupancy = {zero:.6f}, expected 1.0")
    if not ok:
        failures.append(f"occupancy(K_D=0) = {zero}, expected 1.0")

    at_kd = occupancy(common.ASSAY_TOP_CONCENTRATION_NM)
    ok = abs(at_kd - 0.5) < 1e-9
    emit(f"   [{'PASS' if ok else 'FAIL'}] K_D equal to the analyte "
         f"concentration (the definition of K_D): occupancy = {at_kd:.6f}, "
         "expected 0.5")
    if not ok:
        failures.append(f"occupancy(K_D=analyte) = {at_kd}, expected 0.5")

    sample = [1.0, 10.0, 100.0, 1_000.0, 10_000.0, 100_000.0]
    values = [occupancy(kd) for kd in sample]
    ok = all(a > b for a, b in zip(values, values[1:]))
    emit(f"   [{'PASS' if ok else 'FAIL'}] occupancy strictly decreases as "
         f"K_D grows, over {sample[0]:.0f} to {sample[-1]:.0f} nM: "
         + " > ".join(f"{v:.3f}" for v in values))
    if not ok:
        failures.append("occupancy is not strictly decreasing in K_D: "
                         f"{values}")
    emit()
    return failures


def main(argv=None):
    global BROKEN_RULE
    parser = argparse.ArgumentParser(
        description="The fraction of binder occupied at the assay's top "
                    "concentration, as a function of K_D.")
    parser.add_argument("--break-rule", choices=["flip"], default=None,
                        help="Swap in a deliberately wrong formula, to "
                             "confirm the checks catch it. Writes no output "
                             "files.")
    parser.add_argument("--self-test", action="store_true",
                        help="Run with --break-rule flip in a subprocess and "
                             "confirm it fails.")
    args = parser.parse_args(argv)

    if args.self_test:
        return run_self_test()

    BROKEN_RULE = args.break_rule

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    emit("=" * 72)
    emit("OCCUPANCY — how much of the binder is occupied at the assay's top")
    emit("concentration, for a given K_D")
    emit("=" * 72)
    emit()
    emit("K_D, the dissociation constant, is the concentration of the other")
    emit("molecule at which half the binding sites are occupied. The fraction")
    emit("occupied at any concentration follows from that definition alone:")
    emit("[analyte] / ([analyte] + K_D). Nothing here is measured; it is the")
    emit("arithmetic consequence of what K_D means, made regenerable rather")
    emit("than typed into three separate documents by hand.")
    emit()
    if BROKEN_RULE:
        emit(f"   NOTE: --break-rule {BROKEN_RULE} is in effect: the formula "
             "is deliberately wrong.")
        emit("   This run is expected to fail its own checks.")
        emit()

    failures = run_checks(emit)

    conc = common.ASSAY_TOP_CONCENTRATION_NM
    emit(f"1. Occupancy at the assay's top concentration, {conc:.0f} nM")
    emit()
    emit("   | K_D | occupied fraction |")
    emit("   |---|---|")
    rows = []
    for kd in ILLUSTRATIVE_KD_NM:
        frac = occupancy(kd, conc)
        rows.append((kd, frac))
        emit(f"   | {format_kd(kd)} | {100 * frac:.0f}% |")
    emit()
    emit("   A design in the high nanomolar range for human binding at pH 6.5")
    emit("   sits near full occupancy and gives a strong, easily fitted")
    emit("   signal. A design deliberately weakened toward the far end of the")
    emit("   instrument's reportable range gives a progressively smaller")
    emit("   one, and far enough along it reads as no detectable binding at")
    emit("   pH 6.5 as well as at pH 7.4, which fails human binding rather")
    emit("   than demonstrating pH selectivity. Where exactly the instrument")
    emit("   stops returning a fittable curve is not established here: that")
    emit("   depends on the size of the design and the density on the chip")
    emit("   as well as on occupancy, and the organisers have not said.")
    emit()

    if BROKEN_RULE is not None:
        emit("   Running with --break-rule, so no file is written. The point")
        emit("   of this run is the failures above.")
        emit()

    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} FAILURE(S).")
        for failure in failures:
            emit(f"  - {failure}")
    else:
        emit(f"RESULT: PASSED. Formula checks hold; occupancy tabulated for "
             f"{len(ILLUSTRATIVE_KD_NM)} illustrative K_D values at "
             f"{conc:.0f} nM.")
    emit("=" * 72)

    if BROKEN_RULE is not None:
        return 1 if failures else 0

    DERIVED.mkdir(parents=True, exist_ok=True)
    with (DERIVED / "17-occupancy.csv").open("w") as fh:
        fh.write("kd_nm,analyte_nm,occupied_fraction\n")
        for kd, frac in rows:
            fh.write(f"{kd:.1f},{conc:.1f},{frac:.6f}\n")

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / "17-occupancy.md").write_text(
        "# Occupancy at the assay's top concentration\n\n"
        "Computed output of `analysis/17_occupancy_table.py`. "
        "Do not hand-edit.\n\n"
        "K_D, the dissociation constant, is the concentration of the other\n"
        "molecule at which half the binding sites are occupied; it is a\n"
        "measure of affinity, with a smaller number meaning a tighter grip.\n"
        "nM is nanomolar, a billionth of a mole per litre; µM is\n"
        "micromolar, a thousand times looser for the same number. The\n"
        "fraction of binder occupied at any concentration follows from the\n"
        "definition of K_D alone: [analyte] / ([analyte] + K_D). Nothing\n"
        "here is measured; it is the arithmetic consequence of what K_D\n"
        "means, computed once rather than typed by hand into\n"
        "`docs/explainers/07-the-assay-and-what-it-changes.md`,\n"
        "`docs/decisions-log.md` and `CLAUDE.md`, which all previously\n"
        "quoted it.\n\n"
        "```\n" + "\n".join(out) + "\n```\n"
    )
    return 1 if failures else 0


def run_self_test():
    """Confirm the checks actually fail when the formula is wrong."""
    print("=" * 72)
    print("SELF-TEST: do the checks catch a deliberately wrong formula?")
    print("=" * 72)
    print()
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--break-rule", "flip"],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent))
    caught = [line.strip() for line in proc.stdout.splitlines()
              if line.strip().startswith("- ")]
    ok = proc.returncode != 0 and bool(caught)
    print(f"   [{'PASS' if ok else 'FAIL'}] flip")
    print(f"     exit code {proc.returncode}, {len(caught)} check(s) failed")
    for line in caught:
        print(f"       {line}")
    print()
    print("=" * 72)
    if ok:
        print("SELF-TEST PASSED. The checks fail when the formula is wrong.")
    else:
        print("SELF-TEST FAILED. The checks did not notice a wrong formula.")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
