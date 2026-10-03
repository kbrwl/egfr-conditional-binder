#!/usr/bin/env python3
"""
15_histag_counterscreen.py — does any candidate grip the His tag on the target?

Abbreviations, expanded here because each file gets read on its own:
  EGFR     epidermal growth factor receptor, the protein we design against.
  His tag  a short run of histidine residues, usually six, added to the end of a
           protein so it can be caught on a metal column.
  pH       how acidic something is. Tumour tissue sits around 6.5 and blood
           around 7.4, which is the difference the design has to detect.
  i_pTM    BindCraft2's confidence in the interface of a predicted complex, 0 to 1.
  FASTA    a plain text format for a protein sequence.
  SPR      surface plasmon resonance, the instrument the competition measures
           binding with.

WHY THIS STEP EXISTS
--------------------
Both the human and the mouse target carry a C-terminal His tag, and the organisers
expect to screen with it left on (docs/competition-qa-log.md). Histidine is the one
amino acid that gains a charge between pH 7.4 and pH 6.5, so a binder that grips
the tag looks pH-selective and binds anything carrying a His tag. Amir at Anthropic
said such a binder should not be selected, and Tudor at Adaptyv said the organisers
will weight in silico scoring more heavily where a design might be hitting a tag.

This project is more exposed to that than most, for a structural reason. The
pairing rule puts acidic residues on the binder facing the target's own histidines,
H358 and H370. A pocket shaped to hold a protonated histidine holds the histidines
of a tag as readily. The tag is also a floppy tail free to thread into a pocket,
where H358 and H370 are held in place by the fold.

This step is written before any candidate exists, for the same reason step 10 was:
a screen written afterwards turns out to need a measurement the design run did not
save.

DECISIONS, RECORDED BEFORE THE SCRIPT WAS WRITTEN
-------------------------------------------------
1. The predictor is BindCraft2's own AlphaFold2, used through its detargeting
   objective, with the tag as an off-target. Read from BindCraft2's documentation
   and source on 2 October 2026. A negative target weight tells BindCraft2 to push
   the binder away from that target during design, and the tag is given as a
   sequence, which BindCraft2 co-folds with the binder as a disordered region. Its
   output tables then record `i_pTM_detarget`, `i_pAE_detarget` and
   `Interface_Residues_detarget`.

   Rejected: `bindcraft score`. Its help says it scores one design structure on
   checks taken off coordinates and reports the prediction readings from the design's
   own stamp, so it cannot predict a binder against a peptide. Rejected: placing a
   peptide by geometry. A six-residue tail has no fixed pose to place.
   Rejected: a separate prediction run after the fact. It would need a prediction
   runner that does not exist yet, and the in-loop version costs nothing extra and
   steers designs away from the tag instead of only finding it afterwards.

2. 'Accepts it' means: BindCraft2 reads the binder's interface with the tag at an
   `i_pTM_detarget` of at least 0.4. That figure is BindCraft2's own default ceiling
   for rejecting an off-target (`max_detarget_iptm`). It is the tool's default and
   NOT a calibrated threshold for this tag: calibrating one needs known tag binders
   and known non-binders, and neither exists yet. The threshold therefore rests on
   the tool's choice and on nothing measured here.

   The screen asks whether the binder is predicted to bind the tag at all, and not
   only through its acidic pocket. A binder that grips the tag by any residues
   confounds the assay in the same way.

3. A blank is never a pass. BindCraft2's documentation says not to treat a blank
   as zero. A candidate with no off-target reading is reported 'not recorded',
   which means no evidence either way, and is not counted as clear. That is the
   case for every design generated without the tag in the campaign file.

4. A pass is weak evidence. A prediction that finds no interface does not show
   there is none; AlphaFold2 can miss a real one, especially for a floppy tail. The
   screen can show that a risk is present and cannot show that it is absent. The
   methods write-up should say so.

HOW THIS STEP IS CHECKED
------------------------
  - The campaign files step 09 wrote are checked for the settings this depends on:
    the tag as a sequence target with a negative weight, spelling six histidines,
    the sequence file present, and the termini property on.
  - The rule is tested on constructed rows run every time: a low reading, a reading
    just under 0.4, a reading at exactly 0.4, a reading well above, a blank, and a
    value that is not a number.
  - `--break-rule ceiling` removes the threshold and `--break-rule blank-is-zero`
    treats a missing reading as zero. The matching test must then fail, and
    `--self-test` runs both.

Writes:
  results/findings/15-histag-counterscreen.md
  data/derived/15-tag-screen.csv

Run standalone:  python analysis/15_histag_counterscreen.py
Against a campaign:  python analysis/15_histag_counterscreen.py --candidates DIR
"""

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import egfr_common as common  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DERIVED = common.DERIVED
FINDINGS = common.FINDINGS
CONFIG_DIR = ROOT / "design" / "configs"
BASE_CONFIG = CONFIG_DIR / "egfr-domain3-h370.json"
NOTAG_CONFIG = CONFIG_DIR / "egfr-domain3-h370-notag.json"
TAG_FASTA = common.SEQ_DIR / "his-tag-offtarget.fasta"

# The column BindCraft2 writes for an off-target's interface confidence.
READING = "i_pTM_detarget"

# What it actually wrote, on the first real output, 3 October 2026. There is no
# `i_pTM_detarget` column in `2_Refolded/!_Refolded.csv`. Instead `i_pTM` carries
# one value per target, separated by a semicolon and in the order the campaign
# file lists them: `0.37;0.08` means 0.37 against EGFR and 0.08 against the tag.
# Without this the screen read every candidate as 'not recorded' while the
# reading sat in the table, which is the quietest way for a counter-screen to
# fail — it reports nothing rather than reporting wrong.
PACKED_READING = "i_pTM"


def detarget_reading(row):
    """The off-target interface confidence for one candidate, or None.

    Prefers an explicit `i_pTM_detarget` column. Falls back to the second
    semicolon-separated field of `i_pTM`, which is where the pipeline actually
    puts it when the campaign declares one off-target. A row carrying a trailing
    semicolon and nothing after it has no off-target reading, and returns None
    rather than zero: absent is not the same as zero, and zero would read as a
    clean pass.
    """
    direct = row.get(READING)
    if direct not in (None, ""):
        return direct
    packed = row.get(PACKED_READING)
    if not packed or ";" not in str(packed):
        return None
    fields = str(packed).split(";")
    if len(fields) < 2 or not fields[1].strip():
        return None
    return fields[1].strip()

VERDICT_ACCEPTS = "accepts the tag"
VERDICT_CLEAR = "no tag binding predicted"
VERDICT_UNRECORDED = "not recorded"

BROKEN_RULE = None
BREAKABLE_RULES = {
    "ceiling": (
        "ceiling",
        "a design predicted to bind the tag at high confidence would be called "
        "clear, and a tag binder dressed up as a pH switch would reach the "
        "shortlist"),
    "blank-is-zero": (
        "blank is not a pass",
        "a design generated without the tag in its campaign file would be called "
        "clear, when nothing was measured and nothing can be said"),
}


def verdict_for(raw):
    """Classify one candidate from its raw off-target reading.

    `raw` is whatever the campaign table holds: a number as text, an empty string,
    None, or something that is not a number.
    """
    try:
        value = float(raw)
        if value != value:
            raise ValueError
    except (TypeError, ValueError):
        if BROKEN_RULE == "blank-is-zero":
            value = 0.0
        else:
            return VERDICT_UNRECORDED, None
    if BROKEN_RULE != "ceiling" and value >= common.DETARGET_IPTM_CEILING:
        return VERDICT_ACCEPTS, value
    return VERDICT_CLEAR, value


def run_tests(emit):
    """The rule on constructed readings, with the answer written out."""
    emit(f"   The rule, on constructed {READING} readings (ceiling "
         f"{common.DETARGET_IPTM_CEILING}):")
    emit()
    cases = [
        ("0.12", VERDICT_CLEAR, "ceiling", "well below the ceiling"),
        ("0.39", VERDICT_CLEAR, "ceiling", "just below the ceiling"),
        ("0.40", VERDICT_ACCEPTS, "ceiling", "exactly at the ceiling"),
        ("0.71", VERDICT_ACCEPTS, "ceiling", "well above the ceiling"),
        ("", VERDICT_UNRECORDED, "blank is not a pass", "an empty cell"),
        (None, VERDICT_UNRECORDED, "blank is not a pass", "no column at all"),
        ("nan", VERDICT_UNRECORDED, "blank is not a pass", "not a number"),
    ]
    failures = []
    for raw, want, label, why in cases:
        got, _ = verdict_for(raw)
        ok = got == want
        shown = "none" if raw is None else repr(raw)
        emit(f"     [{'PASS' if ok else 'FAIL'}] {shown:7s} -> {got}   ({why})")
        if not ok:
            failures.append(f"{label}: reading {shown} gave '{got}', expected "
                            f"'{want}'")
    emit()
    return failures


def run_campaign_test(emit):
    """Read a constructed campaign folder end to end, through the same reader step
    10 uses, so the path real output will take is exercised and not only the rule.
    Layout and column names are BindCraft2's: `3_Ranked/!_Ranked.csv`, with the
    off-target columns it documents."""
    import tempfile

    emit("   Reading a constructed campaign folder, end to end:")
    emit()
    failures = []
    header = ("rank,design,i_pDAE,i_pTM,i_pTM_detarget,i_pAE_detarget,"
              "Interface_Residues_detarget\n")
    rows = [("design_a1b2c3_seq1", "0.12", VERDICT_CLEAR),
            ("design_d4e5f6_seq1", "0.46", VERDICT_ACCEPTS),
            ("design_g7h8i9_seq1", "", VERDICT_UNRECORDED),
            ("design_j0k1l2_seq1", "0.40", VERDICT_ACCEPTS)]
    with tempfile.TemporaryDirectory() as tmp:
        ranked = Path(tmp) / "campaign" / "3_Ranked"
        ranked.mkdir(parents=True)
        text = header + "".join(
            f"{i},{name},0.8,0.8,{value},0.4,5\n"
            for i, (name, value, _v) in enumerate(rows, start=1))
        (ranked / "!_Ranked.csv").write_text(text)
        quiet = []
        step10 = load_step10()
        metrics, _cols, _src = step10.read_metrics_table(
            Path(tmp) / "campaign", lambda t="": quiet.append(t))
        for name, _value, want in rows:
            got, _ = verdict_for(metrics.get(name, {}).get(READING))
            ok = got == want
            emit(f"     [{'PASS' if ok else 'FAIL'}] {name}: {got}")
            if not ok:
                failures.append(f"campaign reading: {name} gave '{got}', "
                                f"expected '{want}'")
    emit()
    return failures


def check_configs(emit):
    """The campaign files must carry what this screen depends on. This is the one
    part of the step that exercises real files before any candidate exists."""
    problems = []
    for path in (BASE_CONFIG, NOTAG_CONFIG):
        if not path.exists():
            problems.append(f"{path.name} is missing; run analysis/09_trim_target.py")
    if problems:
        return problems

    base = json.loads(BASE_CONFIG.read_text())
    notag = json.loads(NOTAG_CONFIG.read_text())

    def has(cond, ok_text, bad_text):
        emit(f"   [{'PASS' if cond else 'FAIL'}] {ok_text if cond else bad_text}")
        if not cond:
            problems.append(bad_text)

    has(base.get("termini_accessible") is True,
        f"{BASE_CONFIG.name}: termini_accessible is on",
        f"{BASE_CONFIG.name}: termini_accessible is not on, so nothing keeps the "
        "C-terminal end out of the interface")
    has(notag.get("termini_accessible") is True,
        f"{NOTAG_CONFIG.name}: termini_accessible is on",
        f"{NOTAG_CONFIG.name}: termini_accessible is not on")
    tags = [t for t in notag.get("targets", []) if t.get("name") == common.HIS_TAG_NAME]
    has(len(tags) == 1, f"{NOTAG_CONFIG.name}: one target named {common.HIS_TAG_NAME}",
        f"{NOTAG_CONFIG.name}: expected exactly one target named "
        f"{common.HIS_TAG_NAME}, found {len(tags)}")
    if tags:
        has(tags[0].get("weight", 0) < 0,
            f"the tag target has a negative weight ({tags[0].get('weight')}), so "
            "it is avoided",
            "the tag target does not have a negative weight, so BindCraft2 would "
            "try to BIND it")
        fasta = (CONFIG_DIR / tags[0].get("target_path", "")).resolve()
        has(fasta.exists(), f"the tag sequence file exists ({fasta.name})",
            f"the tag sequence file {tags[0].get('target_path')} does not exist")
        if fasta.exists():
            seq = "".join(l.strip() for l in fasta.read_text().splitlines()
                          if not l.startswith(">"))
            has("HHHHHH" in seq,
                f"the tag sequence {seq} contains six consecutive histidines",
                f"the tag sequence {seq!r} does not contain six consecutive "
                "histidines")
    has(len([t for t in base.get("targets", []) if t.get("weight", 1) < 0]) == 0,
        f"{BASE_CONFIG.name}: no off-targets, so it is the plain campaign",
        f"{BASE_CONFIG.name} carries an off-target but is the plain campaign")
    has(notag.get("targets", [{}])[0] == base.get("targets", [{}])[0],
        "both files give the EGFR target and hotspots identically",
        "the two campaign files disagree about the EGFR target or its hotspots")
    emit()
    return problems


def load_step10():
    """Step 10 loaded as a module, for its reader of the campaign's metrics table,
    so the two steps read the same file the same way."""
    spec = importlib.util.spec_from_file_location(
        "step10", Path(__file__).resolve().parent / "10_charge_pair_filter.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv=None):
    global BROKEN_RULE

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--candidates", type=Path, default=None, metavar="DIR",
                        help="a campaign folder from the design run.")
    parser.add_argument("--break-rule", choices=list(BREAKABLE_RULES), default=None,
                        help="switch one rule off, to confirm its test then fails. "
                             "Expected to fail.")
    parser.add_argument("--self-test", action="store_true",
                        help="switch each rule off in turn and confirm a test "
                             "notices.")
    args = parser.parse_args(argv)
    if args.self_test:
        return run_self_test()
    BROKEN_RULE = args.break_rule
    write_files = BROKEN_RULE is None

    out = []

    def emit(text=""):
        print(text)
        out.append(text)

    failures = []

    emit("=" * 72)
    emit("HIS-TAG COUNTER-SCREEN — does any candidate grip the tag on the target?")
    emit("=" * 72)
    emit()
    emit("Both targets carry a C-terminal His tag the organisers expect to leave on.")
    emit("A binder that grips it looks pH-selective and binds anything with a His")
    emit("tag. Our pairing rule builds acidic pockets to grip the target's own")
    emit("histidines, which a tag's histidines fit as readily.")
    emit()
    if BROKEN_RULE:
        emit(f"   NOTE: --break-rule {BROKEN_RULE} is in effect: "
             f"{BREAKABLE_RULES[BROKEN_RULE][1]}.")
        emit("   This run is expected to fail and writes no files.")
        emit()

    emit("1. What the screen rests on")
    emit()
    emit("   Predictor: BindCraft2's own AlphaFold2, through its detargeting")
    emit("   objective, with the tag given as a sequence target at a negative")
    emit("   weight. Not `bindcraft score`, which reads existing structures and")
    emit("   cannot predict a complex.")
    emit(f"   Accepts it: {READING} of at least {common.DETARGET_IPTM_CEILING}, "
         "BindCraft2's own default ceiling")
    emit("   for rejecting an off-target. That is the tool's default and not a")
    emit("   threshold calibrated for this tag. Calibrating one needs known tag")
    emit("   binders and known non-binders, and neither exists yet.")
    emit(f"   The tag: {common.HIS_TAG_SEQUENCE}, weight {common.HIS_TAG_WEIGHT}. "
         "The serine-glycine flank is there")
    emit("   because BindCraft2 samples a window of at least 10 residues from a")
    emit("   sequence target by default, and six residues is shorter than that.")
    emit()

    emit("2. The rule, tested every time")
    emit()
    failures.extend(run_tests(emit))

    failures.extend(run_campaign_test(emit))

    emit("3. The campaign files this depends on")
    emit()
    problems = check_configs(emit)
    failures.extend(problems)

    emit("4. Real candidates")
    emit()
    rows = []
    if args.candidates is None:
        emit("   None given. Run again with --candidates pointing at a campaign")
        emit("   folder once the design run has produced one. Until then this step")
        emit("   has checked its rule and its configuration files and nothing else.")
    elif not args.candidates.is_dir():
        emit(f"   {args.candidates} is not a directory.")
        failures.append(f"{args.candidates} is not a directory")
    else:
        step10 = load_step10()
        metrics, _columns, source = step10.read_metrics_table(args.candidates, emit)
        emit(f"   Metrics table: {source or 'none found'}")
        emit()
        emit(f"   | design | {READING} | verdict |")
        emit("   |---|---|---|")
        for name, row in sorted(metrics.items()):
            verdict, value = verdict_for(detarget_reading(row))
            rows.append((name, value, verdict))
            shown = "—" if value is None else f"{value:.2f}"
            emit(f"   | {name} | {shown} | {verdict} |")
        emit()
        counts = {}
        for _n, _v, verdict in rows:
            counts[verdict] = counts.get(verdict, 0) + 1
        for verdict in (VERDICT_ACCEPTS, VERDICT_CLEAR, VERDICT_UNRECORDED):
            emit(f"   {verdict}: {counts.get(verdict, 0)}")
        if counts.get(VERDICT_UNRECORDED):
            emit()
            emit(f"   {counts[VERDICT_UNRECORDED]} candidate(s) have no off-target "
                 "reading. That is what a campaign run without the tag in its")
            emit("   file produces. They are NOT cleared: nothing was measured.")
    emit()

    emit("=" * 72)
    emit("5. What this changes about the design")
    emit("=" * 72)
    emit()
    emit(f"   Use `{NOTAG_CONFIG.name}` for the campaign, not the plain file, if")
    emit("   the tag-avoiding variant starts. Designs generated without the tag")
    emit("   cannot be screened afterwards without a prediction run that does not")
    emit("   yet exist, so the choice has to be made before the campaign and not")
    emit("   after it.")
    emit()
    emit("   Not run: detargeting against a sequence target has not been tried on")
    emit("   real hardware, and the weight and the ceiling are starting values.")
    emit("   A clear verdict is weak evidence, since a prediction that finds no")
    emit("   interface does not show there is none. It can show a risk is present")
    emit("   and cannot show it is absent.")
    emit()

    if write_files:
        DERIVED.mkdir(parents=True, exist_ok=True)
        with (DERIVED / "15-tag-screen.csv").open("w") as fh:
            fh.write(f"design,{READING},verdict\n")
            for name, value, verdict in rows:
                fh.write(f"{name},{'' if value is None else f'{value:.4f}'},"
                         f"{verdict}\n")
        emit("Wrote data/derived/15-tag-screen.csv"
             + ("" if rows else " (header only; no real candidates yet)"))
        emit()

    emit("=" * 72)
    if failures:
        emit(f"RESULT: {len(failures)} CHECK(S) FAILED —")
        for failure in failures:
            emit(f"  - {failure}")
    elif rows:
        emit(f"RESULT: PASSED. {len(rows)} candidates screened.")
    else:
        emit("RESULT: PASSED. The rule is exercised and the campaign files carry "
             "what it depends on; no real candidates were supplied.")
    emit("=" * 72)

    if write_files:
        FINDINGS.mkdir(parents=True, exist_ok=True)
        (FINDINGS / "15-histag-counterscreen.md").write_text(
            "# The His-tag counter-screen\n\n"
            "Computed output of `analysis/15_histag_counterscreen.py`. "
            "Do not hand-edit.\n\n"
            "Both targets carry a C-terminal His tag the organisers expect to leave "
            "on, and a\nbinder that grips it looks pH-selective and binds anything "
            "with a His tag. Written\nbefore any candidate exists. No candidate has "
            "been screened: the evidence below is\nthe rule's own tests and a check "
            "of the campaign files.\n\n"
            "```\n" + "\n".join(out) + "\n```\n")
    return 1 if failures else 0


def run_self_test():
    """Remove each rule in turn and confirm a test notices."""
    print("=" * 72)
    print("SELF-TEST: is each rule of the tag screen covered by a test?")
    print("=" * 72)
    print()
    failures = []
    for rule, (expected_case, consequence) in BREAKABLE_RULES.items():
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--break-rule", rule],
            capture_output=True, text=True,
            cwd=str(Path(__file__).resolve().parent))
        caught = [ln.strip() for ln in proc.stdout.splitlines()
                  if ln.strip().startswith("- ")]
        named = any(expected_case in ln for ln in caught)
        ok = proc.returncode != 0 and named
        print(f"   [{'PASS' if ok else 'FAIL'}] {rule}")
        print(f"     if unnoticed: {consequence}")
        print(f"     exit code {proc.returncode}, {len(caught)} assertion(s) failed")
        for ln in caught[:3]:
            print(f"       {ln}")
        if not ok:
            failures.append(f"{rule}: "
                            + ("removing this rule broke no test"
                               if proc.returncode == 0
                               else f"tests failed but not the one meant to cover "
                                    f"it ({expected_case})"))
        print()
    print("=" * 72)
    if failures:
        print(f"SELF-TEST FAILED: {len(failures)} rule(s) not covered.")
        for failure in failures:
            print(f"  - {failure}")
    else:
        print(f"SELF-TEST PASSED. All {len(BREAKABLE_RULES)} rules are covered by a "
              "test that")
        print("fails when the rule is removed.")
    print("=" * 72)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
