"""Step 18 — the complete record of a design campaign: every attempt, kept and judged.

WHAT THIS IS FOR

A design campaign spends money to make attempts, and throws most of them away. The
tooling we have records the survivors well and the failures barely at all:

  - BindCraft2 writes its accepted designs into `3_Ranked/` with a ranked table.
    What it tried and discarded is left in `1_Trajectories/` with no summary that
    says which is which or why.
  - `analysis/10_charge_pair_filter.py` reads `3_Ranked/` only. So every attempt
    BindCraft2 discarded is retained on disk and never tested by us at all. There
    is no per-attempt record of what was tried and why it did not survive.

That gap is the reason this step exists. A record containing only successes makes
the pipeline look better than it is, and it makes the rates that size the main run
impossible to check afterwards, because the denominator is missing.

This step produces one row per attempt, covering every attempt, with the outcome at
each gate and the reason. It adds no rule and changes no verdict: the charge-pair
rules live in step 10 and are imported from it rather than restated here, so the two
cannot drift apart. Step 15 reuses step 10 the same way, and this follows that.

Abbreviations, expanded here because this file gets read on its own:
  EGFR   epidermal growth factor receptor, the protein we design against.
  mmCIF  the structure file format BindCraft2 writes its complexes in.
  PDB    the Protein Data Bank, the public archive of measured structures, and
         also the older structure file format.
  aa     amino acids, the unit protein length is counted in.
  SHA-256  a fingerprint of a file's contents. Two files with the same fingerprint
         hold the same bytes. Recorded so that a file can later be shown to be the
         one this inventory described, rather than assumed to be.

THE TWO GATES AN ATTEMPT PASSES THROUGH, WHICH MEASURE DIFFERENT THINGS

  Gate 1, BindCraft2's own thresholds. Did the attempt fold into a confident
  interface? BindCraft2 has no pH term anywhere in what it optimises -- it reports
  binder charge at a hard-coded pH 7.4 as a readout it then ignores -- so passing
  this gate says nothing about the pH switch.

  Gate 2, our charge-pair floor. Does the design carry at least three correct
  charge pairs, place no histidine on the binder facing a histidine on the target,
  and avoid contacting position 442? This is step 10's rule, imported.

An attempt can fail gate 1 and still be interesting: if designs BindCraft2 discards
would have passed our gate, then its thresholds are throwing away the thing we are
paying for, and the main run should be configured differently. Nothing could notice
that before this step, because nothing looked.

WHAT COUNTS AS "KEPT"

Two different things, easy to confuse, so both are reported separately:

  retained      the file exists, on the results volume and copied back here. The
                campaign runner copies the whole project folder, so this should be
                everything the run wrote.
  committed     the file is tracked in version control. Bulk design output
                deliberately is not, because it is large; `.gitignore` admits only
                the final shortlist. That is a decision about repository size and
                says nothing about whether the data is retained.

This step writes its ledger into `data/derived/`, which is tracked, so the record
of what happened survives even where the structures themselves are not committed.

HOW IT HANDLES NOT KNOWING BINDCRAFT2'S TABLE SCHEMA

No EGFR campaign has ever been run, so the exact columns BindCraft2 writes for our
target have never been seen. This step therefore discovers columns rather than
assuming them, and reports what it found instead of failing when a name is absent.
A missing column becomes a blank cell and a line in the findings saying which
column was missing; it never becomes a silently wrong verdict.

Outputs:
  data/derived/18-trajectory-ledger-<campaign>.csv   one row per attempt, every gate
  data/derived/18-file-inventory-<campaign>.csv      one row per retained file, with its size
                                          and fingerprint
  results/findings/18-campaign-inventory-<campaign>.md

Run standalone:  python analysis/18_campaign_inventory.py --campaign DIR
                 python analysis/18_campaign_inventory.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
FINDINGS = ROOT / "results" / "findings"

# BindCraft2's own output subdirectories, in the order an attempt passes through
# them. Used to say which stage a retained file came from, and to keep the
# template structures BindCraft2 ships for other binder formats out of the record
# -- those live in `scaffolds/` and are not output of any run. A smoke run on
# 2 October mistook four of them for designed complexes and reported a false
# numbering failure as a result; `docs/decisions-log.md` has the account.
STAGE_DIRS = ("1_Trajectories", "2_Refolded", "3_Ranked")

# Column names that may carry BindCraft2's own accept/reject decision. Several are
# listed because the schema differs between versions. The first one present is used
# and the findings say which.
#
# `terminated` is the one BindCraft 2 v1.0.1 actually writes, discovered by running
# the validation campaign on 2 October 2026 and reading what came back rather than
# by assuming. It works the opposite way round from the others: it is blank when an
# attempt survived, and when an attempt was stopped it holds the name of the stage
# that stopped it. That makes it better than a plain yes-or-no column, because the
# stage is itself the reason, which is what this step exists to record.
ACCEPT_COLUMNS = ("Accepted", "accepted", "Passed", "passed", "Success",
                  "success", "PassedFilters", "passed_filters",
                  # What BindCraft2 v1.0.1 actually writes, found on the first
                  # real output 3 October 2026: `2_Refolded/!_Refolded.csv` has
                  # an `outcome` column reading `accepted` or `rejected`. Without
                  # this name the verdict sat unread and the step guessed from
                  # which folder a file was in instead.
                  "outcome")

# Values in an accept column that mean the design was accepted. `accepted` is
# here for the `outcome` column above; the rest predate it.
ACCEPTED_VALUES = ("1", "true", "yes", "y", "accepted", "pass", "passed")
TERMINATED_COLUMN = "terminated"

# The stages an attempt passes through, in order, as BindCraft2 names them in its
# own output. Used to turn a `terminated` value into a plain-language reason.
STAGE_MEANING = {
    "screen": "stopped at the screen stage, the first and cheapest check, before "
              "any structure was written",
    "anneal": "stopped at the anneal stage, while the design was still being "
              "relaxed into shape",
    "harden": "stopped at the harden stage",
    "mutate": "stopped at the mutate stage",
    "refine": "stopped at the refine stage",
    "final": "stopped at the final check",
}

# Column names that may carry the design's identifier, matched the same way.
NAME_COLUMNS = ("Design", "design", "Name", "name", "DesignName",
                "design_name", "Trajectory", "trajectory")


def load_step10():
    """Import step 10 so its rules are used rather than restated.

    `10_charge_pair_filter.py` cannot be imported by name because a module name
    may not begin with a digit, so it is loaded from its path. This is the same
    mechanism `analysis/15_histag_counterscreen.py` uses, and it is used for the
    same reason: the charge-pair rule is one rule, and a second copy of it here
    would be a second thing to keep in step with CLAUDE.md's design rule.
    """
    path = Path(__file__).resolve().parent / "10_charge_pair_filter.py"
    spec = importlib.util.spec_from_file_location("step10", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["step10"] = module
    spec.loader.exec_module(module)
    return module


def fingerprint(path, chunk=1 << 20):
    """SHA-256 of a file, read in chunks so a large structure file does not have to
    fit in memory at once. Returns None if the file cannot be read, rather than
    stopping the inventory, because a file that cannot be read is itself a finding.
    """
    digest = hashlib.sha256()
    try:
        with open(path, "rb") as handle:
            while True:
                block = handle.read(chunk)
                if not block:
                    break
                digest.update(block)
    except OSError:
        return None
    return digest.hexdigest()


def stage_of(path, campaign):
    """Which of BindCraft2's stages a file belongs to, or where else it came from.

    Returns one of the names in STAGE_DIRS, or "scaffold" for the template
    structures BindCraft2 ships for other binder formats, or "other" for anything
    else the run wrote, such as its metadata and caches. Every retained file gets a
    label; nothing is dropped for being unrecognised, because an unrecognised file
    is exactly the kind of thing a complete record should show.
    """
    try:
        parts = path.relative_to(campaign).parts
    except ValueError:
        return "other"
    if "scaffolds" in parts:
        return "scaffold"
    for stage in STAGE_DIRS:
        if stage in parts:
            return stage
    return "other"


def inventory_files(campaign):
    """Every file under the campaign folder, with its stage, size and fingerprint.

    This is the "nothing was quietly dropped" record. It covers all files and not
    only structures, because the tables, the campaign metadata and the model-weight
    checkpoints are the provenance of the structures and are worth as much.
    """
    rows = []
    for path in sorted(campaign.rglob("*")):
        if not path.is_file():
            continue
        try:
            size = path.stat().st_size
        except OSError:
            size = None
        rows.append(dict(
            path=str(path.relative_to(campaign)),
            stage=stage_of(path, campaign),
            suffix=path.suffix.lower(),
            bytes=size,
            sha256=fingerprint(path),
        ))
    return rows


def first_present(columns, candidates):
    """The first candidate column name that exists, or None.

    Used instead of a fixed column name because BindCraft2's schema for our target
    has never been seen. Returning None is a reportable outcome and not an error.
    """
    for name in candidates:
        if name in columns:
            return name
    return None


def read_tables(campaign):
    """Every table the run wrote, as (relative path, header, rows).

    Read with the standard library's CSV reader rather than a dataframe library, so
    that a malformed table degrades into fewer rows instead of raising, and so this
    step keeps the same dependency footprint as the rest of `analysis/`.
    """
    tables = {}
    for path in sorted(campaign.rglob("*.csv")):
        try:
            with path.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle)
                rows = list(reader)
                tables[str(path.relative_to(campaign))] = (
                    reader.fieldnames or [], rows)
        except OSError:
            continue
    return tables


def bindcraft_verdicts(tables):
    """What BindCraft2 decided about each attempt, from its own tables.

    Returns (verdicts, notes). `verdicts` maps a design name to a dict holding the
    decision, the table it came from, and every metric that table carried for it,
    so nothing in the table is discarded on the way into the ledger.

    `notes` records what could not be determined and why -- a missing accept column,
    a table with no recognisable name column -- so the findings can say which part of
    the record is incomplete rather than presenting a partial record as a whole one.
    """
    verdicts, notes = {}, []
    for table_name, (columns, rows) in sorted(tables.items()):
        if not rows:
            continue
        if table_name.endswith("_losses.csv"):
            # Handled by `trajectory_traces`, which joins it by the folder it sits
            # in. It has no design-name column by design, so reporting it as
            # unjoinable here would be noise that hides the real gaps.
            continue
        name_column = first_present(columns, NAME_COLUMNS)
        if name_column is None:
            notes.append(f"{table_name}: no recognisable design-name column among "
                         f"{columns}, so its {len(rows)} rows could not be joined "
                         f"to an attempt")
            continue
        accept_column = first_present(columns, ACCEPT_COLUMNS)
        has_terminated = TERMINATED_COLUMN in columns
        if accept_column is None and not has_terminated:
            notes.append(f"{table_name}: no recognisable accept/reject column, so "
                         f"BindCraft2's own decision is taken from which folder a "
                         f"design's files are in instead")
        for row in rows:
            name = (row.get(name_column) or "").strip()
            if not name:
                continue
            record = verdicts.setdefault(name, dict(design=name, tables=[],
                                                    metrics={}))
            record["tables"].append(table_name)
            if accept_column is not None:
                record["bindcraft_accepted"] = row.get(accept_column)
                record["bindcraft_accept_column"] = accept_column
            if has_terminated:
                # Blank means the attempt was not stopped early. A value names the
                # stage that stopped it, which is recorded as the reason.
                stage = (row.get(TERMINATED_COLUMN) or "").strip()
                record["terminated_at"] = stage
                if stage:
                    record["bindcraft_accepted"] = "false"
                    record["bindcraft_accept_column"] = TERMINATED_COLUMN
                    record["termination_reason"] = STAGE_MEANING.get(
                        stage, f"stopped at the {stage} stage")
                else:
                    # A blank `terminated` means the trajectory ran to the end of
                    # the stage pipeline. **That is not the same as a design being
                    # accepted** and must not be recorded as one. An earlier
                    # version asserted acceptance here, and on the round-one short
                    # band it reported one accepted design where BindCraft2
                    # accepted none: the one trajectory that completed went on to
                    # produce ten candidates and every one was then rejected.
                    # Acceptance is a per-design verdict and comes from the accept
                    # column or from a structure in 3_Ranked/; this table is
                    # per-trajectory and cannot carry it.
                    record["completed_pipeline"] = True
            # Keep every metric the table carried, prefixed by its table, so two
            # tables reporting the same metric name do not overwrite each other.
            for key, value in row.items():
                if key and key != name_column:
                    record["metrics"][f"{table_name}:{key}"] = value
    return verdicts, notes


def trajectory_traces(campaign):
    """The per-attempt optimisation trace, joined by the folder it sits in.

    BindCraft2 writes each attempt a `<name>_losses.csv` inside a folder named
    after the attempt, holding one row per optimisation round -- 100 rows for the
    validation run. The table has no column naming the design, because the folder
    already does, so it cannot be joined the way the campaign-wide tables are.

    This is the finest-grained record of what happened inside a single attempt, and
    it exists even for attempts that were stopped before any structure was written.
    For those, it is the only record there is, which is why it is picked up
    separately rather than left in the "could not be joined" pile.

    Returns a map from design name to (row count, final row), the final row being
    the attempt's state when it stopped.
    """
    traces = {}
    for path in sorted(campaign.rglob("*_losses.csv")):
        name = path.parent.name
        try:
            with path.open(newline="", errors="replace") as handle:
                rows = list(csv.DictReader(handle))
        except OSError:
            continue
        if rows:
            traces[name] = (len(rows), rows[-1])
    return traces


def designs_by_stage(campaign):
    """Design names seen at each stage, from the structure files themselves.

    The folder a design's structure is in is BindCraft2's decision made visible: a
    design with a file in `3_Ranked/` was accepted. This is used when no table
    carries an explicit accept column, and as a cross-check when one does. Binder-only
    files, whose names end `_monomer`, are excluded because they hold no interface
    and step 10 cannot score them.
    """
    seen = {stage: set() for stage in STAGE_DIRS}
    for path in campaign.rglob("*.cif"):
        if path.stem.endswith("_monomer"):
            continue
        stage = stage_of(path, campaign)
        if stage in seen:
            seen[stage].add(path.stem)
    return seen


def build_ledger(campaign, scored_by_design=None):
    """One row per attempt, joining every gate, with gaps named rather than hidden.

    `scored_by_design` maps a design name to step 10's verdict for it, where one
    was computed. An attempt with no entry gets `our_verdict` of "not scored" and a
    reason, which is the field this whole step exists to stop being blank.
    """
    scored_by_design = scored_by_design or {}
    tables = read_tables(campaign)
    verdicts, notes = bindcraft_verdicts(tables)
    stages = designs_by_stage(campaign)
    traces = trajectory_traces(campaign)

    names = set(verdicts) | set(traces)
    if stages:
        names |= set().union(*stages.values())
    rows = []
    for name in sorted(names):
        record = verdicts.get(name, {})
        in_ranked = any(name.startswith(s) or name == s for s in stages["3_Ranked"]) \
            or name in stages["3_Ranked"]
        declared = record.get("bindcraft_accepted")
        if declared is not None:
            accepted = str(declared).strip().lower() in ACCEPTED_VALUES
            accepted_from = f"table column {record.get('bindcraft_accept_column')}"
        else:
            accepted = in_ranked
            accepted_from = "presence of a structure in 3_Ranked/"

        scored = scored_by_design.get(name)
        if scored is not None:
            our_verdict = scored.get("verdict")
            our_reasons = "; ".join(scored.get("reasons") or []) or ""
            correct_pairs = scored.get("correct_pairs")
        else:
            our_verdict = "not scored"
            correct_pairs = None
            our_reasons = (
                "no complex structure was available to score"
                if not accepted else
                "accepted by BindCraft2 but absent from the scored set, which "
                "should not happen and is a defect to chase")

        trace_rounds, final_row = traces.get(name, (None, {}))

        rows.append(dict(
            design=name,
            bindcraft_accepted=accepted,
            bindcraft_decision_from=accepted_from,
            bindcraft_terminated_at=record.get("terminated_at", ""),
            bindcraft_reason=record.get("termination_reason", ""),
            in_1_trajectories=name in stages["1_Trajectories"],
            in_2_refolded=name in stages["2_Refolded"],
            in_3_ranked=name in stages["3_Ranked"],
            our_verdict=our_verdict,
            our_correct_pairs=correct_pairs,
            our_reasons=our_reasons,
            trace_rounds=trace_rounds,
            final_state=";".join(
                f"{k}={v}" for k, v in sorted(final_row.items())
                if k in ("phase", "round") or "plddt" in k.lower()
                or "iptm" in k.lower()),
            tables=";".join(record.get("tables", [])),
            metric_count=len(record.get("metrics", {})),
        ))
    return rows, notes, tables


def coverage(rows):
    """How much of the campaign the record actually covers, as counts.

    Reported because a ledger that silently omitted half the attempts would look
    exactly like one that covered everything. The numbers here are what make the
    completeness checkable instead of asserted.
    """
    total = len(rows)
    accepted = sum(1 for r in rows if r["bindcraft_accepted"])
    scored = sum(1 for r in rows if r["our_verdict"] != "not scored")
    return dict(
        attempts=total,
        accepted_by_bindcraft=accepted,
        rejected_by_bindcraft=total - accepted,
        scored_by_us=scored,
        not_scored_by_us=total - scored,
        accepted_but_unscored=sum(
            1 for r in rows
            if r["bindcraft_accepted"] and r["our_verdict"] == "not scored"),
    )


def write_outputs(rows, files, notes, tables, campaign, emit=print):
    """Write the two tables and the findings file."""
    DERIVED.mkdir(parents=True, exist_ok=True)
    FINDINGS.mkdir(parents=True, exist_ok=True)

    ledger_columns = ["design", "bindcraft_accepted", "bindcraft_decision_from",
                      "bindcraft_terminated_at", "bindcraft_reason",
                      "in_1_trajectories", "in_2_refolded", "in_3_ranked",
                      "our_verdict", "our_correct_pairs", "our_reasons",
                      "trace_rounds", "final_state",
                      "tables", "metric_count"]
    # Named after the campaign, so running this on a second campaign does not
    # overwrite the first. The fixed names this used to write meant the long
    # band's ledger replaced the short band's on 3 October 2026, and the two
    # campaigns being separable is the whole point of recording per
    # configuration rather than pooling. `stem` is the campaign folder's name.
    stem = Path(campaign).name
    with (DERIVED / f"18-trajectory-ledger-{stem}.csv").open(
            "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ledger_columns)
        writer.writeheader()
        writer.writerows(rows)

    with (DERIVED / f"18-file-inventory-{stem}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["path", "stage", "suffix", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(files)

    counts = coverage(rows)
    by_stage = {}
    for row in files:
        entry = by_stage.setdefault(row["stage"], [0, 0])
        entry[0] += 1
        entry[1] += row["bytes"] or 0

    lines = []
    lines.append("# Finding 18 — the complete record of the campaign")
    lines.append("")
    lines.append("Computed by `analysis/18_campaign_inventory.py`. "
                 "Do not hand-edit; rerun the script.")
    lines.append("")
    lines.append(f"Campaign folder: `{campaign}`")
    lines.append("")
    lines.append("## 1. Coverage")
    lines.append("")
    lines.append("Every attempt the campaign made, and what is known about each. "
                 "A record holding only the survivors would make the pipeline look "
                 "better than it is and would leave the acceptance rate with no "
                 "denominator, so the counts below are the point of this file.")
    lines.append("")
    lines.append("| | count |")
    lines.append("|---|---|")
    for key, value in counts.items():
        lines.append(f"| {key.replace('_', ' ')} | {value} |")
    lines.append("")
    if counts["accepted_but_unscored"]:
        lines.append(f"**{counts['accepted_but_unscored']} designs were accepted by "
                     "BindCraft2 and carry no verdict from us.** That combination "
                     "should not occur and is a defect to chase rather than a "
                     "property of the run.")
        lines.append("")
    lines.append("## 2. Files retained")
    lines.append("")
    lines.append("Retained means the file exists and was copied back. It is a "
                 "separate question from whether it is committed to version "
                 "control: bulk design output deliberately is not, because it is "
                 "large, and `.gitignore` admits only the final shortlist. The "
                 "fingerprints are in `data/derived/18-file-inventory-" + stem + ".csv`, which "
                 "is tracked, so the record survives even where the structures "
                 "themselves are not committed.")
    lines.append("")
    lines.append("| stage | files | bytes |")
    lines.append("|---|---|---|")
    for stage, (count, size) in sorted(by_stage.items()):
        lines.append(f"| {stage} | {count} | {size:,} |")
    lines.append("")
    lines.append("## 3. Tables the run wrote")
    lines.append("")
    for name, (columns, table_rows) in sorted(tables.items()):
        lines.append(f"- `{name}` — {len(table_rows)} rows, "
                     f"{len(columns)} columns")
    if not tables:
        lines.append("- none")
    lines.append("")
    lines.append("## 4. What could not be determined, and why")
    lines.append("")
    if notes:
        for note in notes:
            lines.append(f"- {note}")
    else:
        lines.append("- Nothing. Every table joined to an attempt and every "
                     "attempt carried a decision.")
    lines.append("")
    lines.append("## 5. What this does not establish")
    lines.append("")
    lines.append("- Nothing here says a design works. It says what was tried, what "
                 "survived each gate and why, which is the record needed to check "
                 "the rates afterwards.")
    lines.append("- The charge-pair verdicts are step 10's and carry step 10's "
                 "limits: a predicted structure scored against a predicted "
                 "interface, with the pH switch inferred from geometry and never "
                 "measured.")
    (FINDINGS / f"18-campaign-inventory-{stem}.md").write_text("\n".join(lines) + "\n")

    emit(f"   data/derived/18-trajectory-ledger-{stem}.csv   {len(rows)} attempts")
    emit(f"   data/derived/18-file-inventory-{stem}.csv      {len(files)} files")
    emit(f"   results/findings/18-campaign-inventory-{stem}.md")
    return counts


def _self_test():
    """Exercise the inventory on a constructed campaign folder, including the cases
    that have gone wrong before.

    CLAUDE.md's rule is that a check never seen to fail is not known to work, so
    each case below is paired with a deliberately broken input whose outcome must
    differ.
    """
    import tempfile

    failures = []

    def check(label, got, want):
        if got != want:
            failures.append(f"{label}: got {got!r}, wanted {want!r}")

    with tempfile.TemporaryDirectory() as temp:
        campaign = Path(temp) / "campaign"
        for relative in ("1_Trajectories/d1.cif", "1_Trajectories/d2.cif",
                         "1_Trajectories/d3.cif",
                         "2_Refolded/d1.cif",
                         "3_Ranked/d1.cif", "3_Ranked/d1_monomer.cif",
                         "scaffolds/VHH.cif"):
            path = campaign / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("x")
        (campaign / "1_Trajectories" / "trajectory_stats.csv").write_text(
            "Design,Accepted,i_pTM\nd1,True,0.8\nd2,False,0.3\nd3,False,0.2\n")

        stages = designs_by_stage(campaign)
        check("monomer files are not counted as designs",
              sorted(stages["3_Ranked"]), ["d1"])
        check("the shipped scaffold is not counted as a design",
              "VHH" in stages["1_Trajectories"], False)

        rows, notes, tables = build_ledger(
            campaign, scored_by_design={"d1": dict(verdict="meets the pair target",
                                                   reasons=[], correct_pairs=4)})
        by_name = {r["design"]: r for r in rows}

        check("every attempt appears, not only the accepted one",
              sorted(by_name), ["d1", "d2", "d3"])
        check("BindCraft2's own accept column is used when present",
              by_name["d1"]["bindcraft_decision_from"], "table column Accepted")
        check("a rejected attempt is recorded as rejected",
              by_name["d2"]["bindcraft_accepted"], False)
        check("a rejected attempt still gets a row with a stated reason",
              by_name["d2"]["our_verdict"], "not scored")
        check("the reason for not scoring a rejected attempt is recorded",
              by_name["d2"]["our_reasons"],
              "no complex structure was available to score")
        check("our verdict is carried through for a scored design",
              (by_name["d1"]["our_verdict"], by_name["d1"]["our_correct_pairs"]),
              ("meets the pair target", 4))
        check("metrics from the table are counted, not dropped",
              by_name["d2"]["metric_count"] > 0, True)

        counts = coverage(rows)
        check("coverage counts all three attempts", counts["attempts"], 3)
        check("coverage separates accepted from rejected",
              (counts["accepted_by_bindcraft"], counts["rejected_by_bindcraft"]),
              (1, 2))
        check("coverage reports what we did not score",
              counts["not_scored_by_us"], 2)
        check("an accepted design with no verdict is flagged as a defect",
              counts["accepted_but_unscored"], 0)

        files = inventory_files(campaign)
        check("every file is inventoried, scaffolds included and labelled",
              len(files), 8)
        check("the scaffold is labelled as a scaffold and not as output",
              [f["stage"] for f in files if f["path"].endswith("VHH.cif")],
              ["scaffold"])
        check("files carry a fingerprint",
              all(f["sha256"] for f in files if f["suffix"] == ".cif"), True)
        # Identical contents must fingerprint identically, or the record cannot be
        # used to show a file is the one described.
        digests = {f["sha256"] for f in files if f["suffix"] == ".cif"}
        check("identical files share a fingerprint", len(digests), 1)

        # Broken on purpose: with the accept column removed, the decision must fall
        # back to which folder the structure is in and say so, rather than silently
        # marking everything rejected.
        (campaign / "1_Trajectories" / "trajectory_stats.csv").write_text(
            "Design,i_pTM\nd1,0.8\nd2,0.3\nd3,0.2\n")
        rows2, notes2, _ = build_ledger(campaign)
        by_name2 = {r["design"]: r for r in rows2}
        check("without an accept column the folder decides",
              by_name2["d1"]["bindcraft_accepted"], True)
        check("and it says where the decision came from",
              by_name2["d1"]["bindcraft_decision_from"],
              "presence of a structure in 3_Ranked/")
        if not any("accept/reject" in n for n in notes2):
            failures.append("a missing accept column should be reported in the "
                            "notes and was not")

        # Broken on purpose: a table with no recognisable name column must be
        # reported rather than quietly contributing nothing.
        (campaign / "stray.csv").write_text("foo,bar\n1,2\n")
        _, notes3, _ = build_ledger(campaign)
        if not any("stray.csv" in n for n in notes3):
            failures.append("a table with no design-name column should be reported "
                            "in the notes and was not")

    # BindCraft 2 v1.0.1's real schema, as the validation campaign wrote it on
    # 2 October 2026: a `terminated` column that is blank for a surviving attempt
    # and names the stage that stopped a failed one, plus a per-attempt losses
    # table joined by its folder name rather than by a column.
    with tempfile.TemporaryDirectory() as temp:
        campaign = Path(temp) / "real"
        (campaign / "1_Trajectories").mkdir(parents=True)
        (campaign / "1_Trajectories" / "!_Trajectories.csv").write_text(
            "trajectory,design,length,terminated\n"
            "1,run_l79_aaa,79,screen\n"
            "2,run_l62_bbb,62,\n")
        for name, rounds in (("run_l79_aaa", 100), ("run_l62_bbb", 40)):
            folder = campaign / "1_Trajectories" / name
            folder.mkdir(parents=True)
            body = "".join(
                f"screen,{i},0.35,0.26\n" for i in range(rounds))
            (folder / f"{name}_losses.csv").write_text(
                "phase,round,EGFR_domain3.target_plddt,EGFR_domain3.iptm\n" + body)
        (campaign / "3_Ranked").mkdir(parents=True)
        (campaign / "3_Ranked" / "run_l62_bbb.cif").write_text("x")

        rows, notes, _ = build_ledger(campaign)
        by_name = {r["design"]: r for r in rows}
        check("a blank terminated column means the attempt survived",
              by_name["run_l62_bbb"]["bindcraft_accepted"], True)
        check("a filled terminated column means the attempt was stopped",
              by_name["run_l79_aaa"]["bindcraft_accepted"], False)
        check("the stage that stopped it is recorded",
              by_name["run_l79_aaa"]["bindcraft_terminated_at"], "screen")
        check("and is turned into a plain-language reason",
              by_name["run_l79_aaa"]["bindcraft_reason"],
              STAGE_MEANING["screen"])
        check("the decision is attributed to the terminated column",
              by_name["run_l79_aaa"]["bindcraft_decision_from"],
              "table column terminated")
        check("the per-attempt trace is joined by its folder name",
              (by_name["run_l79_aaa"]["trace_rounds"],
               by_name["run_l62_bbb"]["trace_rounds"]), (100, 40))
        check("the attempt's final state is carried, including why it failed",
              "EGFR_domain3.target_plddt=0.35"
              in by_name["run_l79_aaa"]["final_state"], True)
        # An attempt stopped before any structure was written still gets a full
        # row. This is the case the whole step exists for: nothing else in the
        # pipeline records it at all.
        check("an attempt with no structure is still recorded",
              by_name["run_l79_aaa"]["in_3_ranked"], False)
        if not by_name["run_l79_aaa"]["our_reasons"]:
            failures.append("an attempt with no structure must still carry a "
                            "stated reason for not being scored")
        # Broken on purpose: the losses table has no design-name column, so if the
        # folder-name join were removed it would land in the notes instead of
        # becoming a trace. Confirm it did not.
        if any("_losses.csv" in n for n in notes):
            failures.append("the per-attempt losses table should be joined by its "
                            "folder name, not reported as unjoinable")

    if failures:
        print("\n".join(f"FAIL  {f}" for f in failures))
        return False
    print("self-test PASSED: every attempt is recorded with its outcome at both "
          "gates, rejected attempts carry a stated reason, scaffolds and "
          "binder-only files are excluded from the design count, files are "
          "fingerprinted, and each broken-on-purpose case was reported rather "
          "than hidden.")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--campaign", type=Path, default=None, metavar="DIR",
                        help="a campaign folder returned by the design run")
    parser.add_argument("--scored", type=Path, default=None, metavar="CSV",
                        help="step 10's candidate summary, to join our verdicts "
                             "from; defaults to data/derived/10-candidate-summary.csv")
    parser.add_argument("--self-test", action="store_true",
                        help="run the checks and exit")
    args = parser.parse_args(argv)

    if args.self_test:
        return 0 if _self_test() else 1

    if args.campaign is None:
        print("No campaign given. Run again with --campaign pointing at a folder "
              "the design run returned. Nothing was written.")
        return 1
    if not args.campaign.is_dir():
        print(f"Not a directory: {args.campaign}")
        return 1

    scored_path = args.scored or (DERIVED / "10-candidate-summary.csv")
    scored_by_design = {}
    if scored_path.exists():
        with scored_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                name = (row.get("design") or "").strip()
                if not name:
                    continue
                pairs = row.get("correct_pairs")
                scored_by_design[name] = dict(
                    verdict=row.get("verdict"),
                    reasons=[row["reasons"]] if row.get("reasons") else [],
                    correct_pairs=int(pairs) if (pairs or "").isdigit() else None)
        print(f"Joined {len(scored_by_design)} verdicts from {scored_path}")
    else:
        print(f"No scored summary at {scored_path}. Every attempt will be "
              f"recorded as not scored, with that stated as the reason.")

    rows, notes, tables = build_ledger(args.campaign, scored_by_design)
    files = inventory_files(args.campaign)
    print("\nFiles written")
    counts = write_outputs(rows, files, notes, tables, args.campaign)
    print("\nCoverage")
    for key, value in counts.items():
        print(f"   {key.replace('_', ' '):28s} {value}")
    if notes:
        print("\nWhat could not be determined")
        for note in notes:
            print(f"   - {note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
