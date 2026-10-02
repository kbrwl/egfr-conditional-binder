"""
bindcraft2_smoke.py — get BindCraft2 running end to end on a rented GPU.

PURPOSE: a working pipeline, not a good design. This runs BindCraft2 against the
target that ships with it (PD-L1), to prove the whole chain works: image builds,
GPU is visible to JAX, model weights download, a campaign runs, and a sequence
comes out. Only once that works do we point it at EGFR.

Abbreviations, expanded on first use because this file gets read alone:
  GPU   graphics processing unit, the hardware the folding models need.
  VRAM  the memory on that GPU.
  JAX   the numerical library BindCraft2's folding models run on. If it cannot
        see the GPU it silently falls back to the processor and runs far slower,
        which is the first thing to check when a run seems slow.
  CUDA  NVIDIA's GPU programming layer. Version mismatches between CUDA, the
        driver and the JAX wheels are the usual cause of install failures.
  PD-L1 programmed death-ligand 1, an unrelated protein used here only because
        BindCraft2 ships a worked example against it.

WHY MODAL RATHER THAN RUNPOD OR COLAB
-------------------------------------
Checked 1 October 2026 against modal.com/pricing:
  - billing is per second, and idle time is not billed, which suits an exercise
    that is mostly failed attempts
  - the Starter plan includes $30 per month of free compute, which may cover
    this entire exercise
  - it takes a container image definition directly, so BindCraft2's own install
    script runs once at image build time and is then cached. Image builds are
    billed as processor time rather than GPU time.
Colab was rejected because sessions time out and disconnect and the GPU type is
not guaranteed, both of which are bad for a multi-hour job. RunPod is a
reasonable second choice but bills by the hour on a running pod, so idle time
during debugging costs money.

COST, AT THE PRICES READ ON 1 OCTOBER 2026
------------------------------------------
  Nvidia L4        $0.000222/s  = $0.80/hour
  Nvidia A10       $0.000306/s  = $1.10/hour
  Nvidia A100 40GB $0.000583/s  = $2.10/hour
  Nvidia A100 80GB $0.000694/s  = $2.50/hour
The $30 monthly free allowance is about 37 hours on an L4 or 14 on an A100 40GB.

BindCraft2's documentation does not state a minimum VRAM, so the choice of card
was a judgement rather than a requirement when nothing had run yet. It no longer
is: the PD-L1 smoke run (2 October 2026) measured peak card memory at 17,950 MiB
on an A100 40GB (40,960 MiB), sampled once a second across two workers sharing
the card. That is an upper bound -- JAX can reserve more than it needs at
start-up -- and it fits inside an L4's 24,564 MiB with about 6 GiB of headroom
even without correcting for the over-estimate. So the default below is now L4,
which is the larger lever: $30 of free credit buys roughly 37 hours on an L4
against 14 on an A100 40GB (see the costs below), for a card that has not been
shown to be too small. No speed ratio between the two cards has been measured,
so the trajectory-rate figures from the A100 smoke run do not carry over
directly; the pilot campaign (docs/decisions-log.md) is what measures that on
this card. There is no runtime flag to pick the card -- `gpu=DEFAULT_GPU` is
fixed at decoration time -- so switching back to A100-40GB if an out-of-memory
failure appears means changing the constant below and redeploying; the earlier
docstring's "drop to L4 with --gpu l4" described a flag that was never
implemented, and is corrected here rather than repeated.

HOW TO RUN
----------
One-time, and only you can do it because it opens a browser:
    ./.venv/bin/modal setup

Then:
    ./.venv/bin/modal run design/modal/bindcraft2_smoke.py::check
        Builds the image and prints what GPU JAX can see. Cheap, a minute or two
        of GPU time. Run this first: it separates "the environment is broken"
        from "the design run is broken".

    ./.venv/bin/modal run design/modal/bindcraft2_smoke.py::smoke
        The actual end-to-end run against PD-L1. Writes results into a Modal
        volume and prints any sequences it produced.

    ./.venv/bin/modal volume get bindcraft2-results / ./results/candidates/
        Copies whatever it produced back to this machine.
"""

import os

import modal

APP_NAME = "bindcraft2-smoke"
VOLUME_NAME = "bindcraft2-results"
WEIGHTS_VOLUME = "bindcraft2-weights"

# Default card. See the note above: switched from A100-40GB to L4 on 2 October
# 2026 once the smoke run showed peak memory fits inside it.
DEFAULT_GPU = "L4"

# BindCraft2's own installer is run at image build time, so it is cached and not
# repeated on every run. The repository is cloned to a fixed path because its
# documentation warns that moving or renaming the directory breaks the
# `bindcraft` command, the install being editable.
INSTALL_DIR = "/opt/BindCraft2"

image = (
    # A CUDA development image rather than a plain Python one, because the
    # installer picks accelerator wheels by looking at the machine and needs the
    # CUDA toolchain present to choose the GPU ones.
    modal.Image.from_registry(
        "nvidia/cuda:12.8.1-devel-ubuntu24.04", add_python="3.12"
    )
    .apt_install("git", "wget", "curl", "build-essential", "ca-certificates")
    .run_commands(
        f"git clone --depth 1 https://github.com/PacesaLab/BindCraft2.git {INSTALL_DIR}",
        # install.sh creates its own virtual environment inside the clone.
        f"cd {INSTALL_DIR} && bash install.sh",
    )
    # Model weights are fetched into a volume rather than baked into the image,
    # so a rebuild does not re-download roughly 20 GB.
    .env({"BINDCRAFT_HOME": INSTALL_DIR})
    # The numbering check shared with analysis/10, copied in last so a change to it
    # does not rebuild the image above. Plain string paths rather than pathlib
    # arithmetic, because this file is imported again inside the container, where
    # the directory layout is different and parent-directory arithmetic can fail.
    .add_local_file(
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "..", "analysis", "target_numbering.py"),
        "/opt/target_numbering.py")
)

# The image is passed to the app so that every function runs inside it. The first
# version of this file defined `image` above and never used it, so both functions
# ran in Modal's default container, where /opt/BindCraft2 does not exist, and the
# check below exited successfully having tested nothing.
app = modal.App(APP_NAME, image=image)
results = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)
weights = modal.Volume.from_name(WEIGHTS_VOLUME, create_if_missing=True)


def _run(cmd, cwd=INSTALL_DIR, venv=True, want_output=False):
    """Run a shell command inside the BindCraft2 environment and stream output.

    `venv=False` skips activating BindCraft2's virtual environment, for commands
    that are not part of it such as nvidia-smi: the activation path is relative to
    `cwd`, so running those from "/" with it switched on fails before the command
    starts. `want_output=True` returns (exit code, stdout) instead of the exit
    code alone.
    """
    import subprocess

    activate = ". .venv/bin/activate && " if venv else ""
    full = f"cd {cwd} && {activate}{cmd}"
    print(f"\n$ {cmd}\n", flush=True)
    proc = subprocess.run(["bash", "-lc", full], capture_output=True, text=True)
    if proc.stdout:
        print(proc.stdout, flush=True)
    if proc.stderr:
        print("--- stderr ---", flush=True)
        print(proc.stderr, flush=True)
    if want_output:
        return proc.returncode, proc.stdout
    return proc.returncode


# ---------------------------------------------------------------------------
# What the smoke run measures besides "did it work"
# ---------------------------------------------------------------------------
#
# Explainer 04 has no figure for time per trajectory, peak graphics-card memory,
# or the size of the complex, because the only ones available were quoted from
# documentation for a different version of the tool. These come from our own
# hardware. The size of the complex matters as much as the other two: run time and
# memory scale with the total residue count, so a time measured on PD-L1 only
# transfers to EGFR through that number.
#
# And the question flagged in analysis/09: BindCraft2 registers a filter metric
# called Target_Crop_Length, which suggests it may crop the target. Its source says
# the metric counts non-padding residues, and padding is a batching device, so the
# reading so far is that nothing is cropped. That is a reading and not a
# measurement. If the returned target were cropped and renumbered, analysis/10
# would translate its numbers into ours wrongly and still print confident
# verdicts, so the check is built into the run: the target chain of every returned
# complex is compared against the structure file that went in, and the answer is
# written down whichever way it falls.

def _start_gpu_poll(path):
    """Sample the card once a second into a file while the design runs.

    Returns the process, or None if nvidia-smi is not there, in which case the
    report says memory was not measured rather than omitting the line.
    """
    import subprocess

    try:
        return subprocess.Popen(
            ["nvidia-smi", "--query-gpu=name,memory.used,memory.total",
             "--format=csv,noheader,nounits", "-l", "1"],
            stdout=open(path, "w"), stderr=subprocess.DEVNULL)
    except OSError:
        return None


def _read_gpu_poll(path):
    """(card name, peak memory used in MiB, total memory in MiB) from the samples.

    Peak is the largest figure seen at one-second sampling, so a short spike
    between samples is missed and the true peak is at least this. The figure is
    the card's total used memory, which JAX may inflate by reserving most of the
    card at start-up whether or not it needs it, so it bounds the requirement
    from above and does not state it.
    """
    name, peak, total = None, None, None
    try:
        for line in open(path).read().splitlines():
            parts = [x.strip() for x in line.split(",")]
            if len(parts) != 3:
                continue
            name, used, total = parts[0], int(parts[1]), int(parts[2])
            peak = used if peak is None else max(peak, used)
    except (OSError, ValueError):
        pass
    return name, peak, total


def _find_input_structure(campaign_json):
    """The structure file the campaign was given, or every candidate for it.

    The PD-L1 example names a preset target ("hPDL1") rather than a file path, so
    the file is looked for under the install. Returns (path or None, candidates).
    Not finding exactly one is reported as such: a guess here would make the
    numbering comparison a comparison against the wrong file.
    """
    import json
    import pathlib

    root = pathlib.Path(INSTALL_DIR)
    explicit = None
    try:
        spec = json.loads((root / campaign_json).read_text())
        for entry in spec.get("targets", []) or []:
            if entry.get("target_path"):
                explicit = (root / campaign_json).parent / entry["target_path"]
        key = str(spec.get("target", "")).lower()
    except (OSError, ValueError):
        key = ""
    if explicit is not None and explicit.exists():
        return explicit, [explicit]
    found = []
    for suffix in ("*.pdb", "*.cif"):
        for path in root.rglob(suffix):
            text = str(path).lower()
            if ".venv" in text or "/results/" in text or "ranked" in text:
                continue
            if key and key in path.name.lower().replace("_", ""):
                found.append(path)
    return (found[0] if len(found) == 1 else None), found


def _count_trajectories(project):
    """Trajectory count read from whatever the run wrote, and where it came from.

    BindCraft2 attempts trajectories until enough designs pass, so the number
    attempted is not the number requested. It is taken as the row count of a table
    in the project folder whose path mentions trajectories, and None if there is
    no such table, in which case the report prints every table's row count so the
    right one can be chosen by hand rather than guessed.
    """
    tables = {}
    for path in project.rglob("*.csv"):
        try:
            tables[str(path.relative_to(project))] = max(
                0, len(path.read_text(errors="replace").splitlines()) - 1)
        except OSError:
            continue
    for name, rows in sorted(tables.items()):
        if "traj" in name.lower() and rows > 0:
            return rows, name, tables
    return None, None, tables


@app.function(gpu=DEFAULT_GPU, timeout=60 * 30,
              volumes={"/weights": weights})
def check():
    """Cheapest possible test: does the GPU exist and can JAX see it?

    Run this before the design job. If JAX reports a processor device rather than
    a GPU one, the design run will work but take so long it looks broken, and the
    BindCraft2 troubleshooting notes call that out as the usual first fault.
    """
    problems = []

    rc, _ = _run("nvidia-smi", cwd="/", venv=False, want_output=True)
    if rc != 0:
        problems.append("nvidia-smi failed: no GPU is visible to the container")

    rc, out = _run('python -c "'
                   "import jax; print('jax', jax.__version__); "
                   "print('devices', jax.devices()); "
                   "print('default backend', jax.default_backend())"
                   '"', want_output=True)
    if rc != 0:
        problems.append("JAX did not import inside BindCraft2's environment")
    elif "default backend gpu" not in out:
        problems.append("JAX is not using the GPU: its default backend is not "
                        "'gpu', so a design run would fall back to the processor")

    # The subcommands and listing flags below are BindCraft2's own, read from its
    # cli.py on 1 October 2026. `--list-targets` is a flag of `design` rather
    # than of `bindcraft` itself, which the earlier version of this file had
    # wrong and which would have burned a GPU minute finding out.
    if _run("bindcraft design --help") != 0:
        problems.append("`bindcraft design --help` failed: the command is not "
                        "installed or not on the path")
    _run("bindcraft design --list-targets || true")
    _run("bindcraft design --list-settings | head -40 || true")

    # A check that exits successfully whatever happens inside it is not a check.
    # The first version of this function did exactly that.
    if problems:
        raise RuntimeError("check FAILED:\n  - " + "\n  - ".join(problems))
    print("\ncheck PASSED: GPU visible, JAX on the GPU backend, bindcraft installed.")


@app.function(gpu=DEFAULT_GPU, timeout=60 * 60 * 4,
              volumes={"/results": results, "/weights": weights})
def smoke(target: str = "examples/pdl1.json"):
    """End-to-end run against BindCraft2's own example target.

    Deliberately not our epitope. The point is to find out whether the pipeline
    works at all, and using the shipped example removes our own configuration as
    a possible cause of failure.
    """
    import pathlib
    import shutil

    rc = _run("bindcraft fetch-weights || true")
    print(f"fetch-weights exit code {rc}")

    import json
    import sys
    import time

    poll_path = "/tmp/gpu-poll.csv"
    poller = _start_gpu_poll(poll_path)
    started = time.monotonic()
    rc = _run(f"bindcraft design {target}")
    wall_seconds = time.monotonic() - started
    if poller is not None:
        poller.terminate()
    print(f"\ndesign exit code {rc}, wall-clock {wall_seconds:.0f} s")

    # Copy anything that looks like output back into the volume, then report any
    # sequences found, because a sequence is the thing that proves the chain.
    #
    # TWO THINGS HERE WERE WRONG AND ARE THE REASON THIS BLOCK IS COMMENTED.
    #
    # The patterns did not include *.cif. BindCraft2 writes its designed
    # complexes as mmCIF, not PDB -- `3_Ranked/<design>_seq<n>[_<target>].cif` --
    # so the old list copied back every summary table and none of the structures.
    # That matters more than it looks: BindCraft2 does not record which binder
    # residue faces which target residue anywhere in its CSVs, under any setting.
    # It writes two separate residue lists. The pairing the charge-pair filter in
    # analysis/10 runs on has to be recomputed from the complex itself, so losing
    # the .cif files means losing the run and paying to generate it again.
    #
    # The copy also flattened everything into one directory. analysis/10 finds
    # candidates by looking in 3_Ranked/ and skipping the *_monomer.cif files
    # that hold the binder alone, and flattening both destroys that and lets two
    # files with the same base name overwrite each other.
    out = pathlib.Path("/results/pdl1-smoke")
    out.mkdir(parents=True, exist_ok=True)
    produced = []
    source = pathlib.Path(INSTALL_DIR)
    for pattern in ("*.csv", "*.cif", "*.fasta", "*.fa", "*.pdb", "*.json",
                    "*.log", "*.txt"):
        for path in source.rglob(pattern):
            if ".venv" in str(path) or "/examples/" in str(path):
                continue
            relative = path.relative_to(source)
            destination = out / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            try:
                shutil.copy2(path, destination)
                produced.append(str(relative))
            except Exception as exc:                     # noqa: BLE001
                print(f"could not copy {path}: {exc}")
    results.commit()

    print(f"\nFiles copied to the volume: {len(produced)}")
    for name in sorted(set(produced)):
        print(f"  {name}")

    # Restricted to BindCraft2's own output subdirectories. The broad copy above
    # also pulls in template scaffold structures BindCraft2 ships for other
    # binder modalities (scaffolds/ARP.cif, Fab.cif, VHH.cif, scFv.cif, found by
    # the 2 October smoke run) -- unrelated structures that correctly fail the
    # numbering check below and would otherwise make its verdict say the target
    # "did not come back" when every real designed complex did.
    OUTPUT_DIRS = ("1_Trajectories", "2_Refolded", "3_Ranked")
    complexes = [p for p in out.rglob("*.cif")
                 if not p.stem.endswith("_monomer")
                 and any(d in p.parts for d in OUTPUT_DIRS)]
    print(f"\nDesigned complexes (what analysis/10 reads): {len(complexes)}")
    if not complexes:
        print("  NONE. Without these the charge-pair filter has nothing to run")
        print("  on, because the pairing it needs is in the structures and not")
        print("  in any table BindCraft2 writes. Find out where they went before")
        print("  spending anything on a real campaign.")

    for path in sorted(out.rglob("*.csv")):
        text = path.read_text(errors="replace").splitlines()
        print(f"\n--- {path.relative_to(out)}, first 5 lines ---")
        for line in text[:5]:
            print(line)

    # ---- the measurements this run exists to produce -----------------------
    sys.path.insert(0, "/opt")
    import target_numbering

    card, peak_mib, total_mib = _read_gpu_poll(poll_path)
    try:
        project = pathlib.Path(INSTALL_DIR) / json.loads(
            (pathlib.Path(INSTALL_DIR) / target).read_text())["project_folder"]
    except (OSError, ValueError, KeyError):
        project = pathlib.Path(INSTALL_DIR) / "results"
    trajectories, trajectory_table, tables = _count_trajectories(project)
    input_path, input_candidates = _find_input_structure(target)

    numbering_results, sizes = [], []
    if input_path is not None:
        for path in sorted(complexes):
            result = target_numbering.check_complex(input_path, path)
            result["file"] = str(path.relative_to(out))
            numbering_results.append(result)
            sizes.append(result["complex_residues"])

    report = dict(
        design_exit_code=rc, wall_clock_seconds=round(wall_seconds, 1),
        card=card, card_memory_mib=total_mib, peak_memory_used_mib=peak_mib,
        trajectories=trajectories, trajectory_table=trajectory_table,
        seconds_per_trajectory=(round(wall_seconds / trajectories, 1)
                                if trajectories else None),
        table_row_counts=tables,
        input_structure=str(input_path) if input_path else None,
        input_candidates=[str(c) for c in input_candidates],
        complexes_found=len(complexes),
        complex_residues=sorted(set(sizes)),
        numbering=[{k: v for k, v in r.items()
                    if k not in ("missing_from_output", "not_in_input",
                                 "wrong_identity")}
                   | {"missing_count": len(r.get("missing_from_output", [])),
                      "not_in_input_count": len(r.get("not_in_input", [])),
                      "wrong_identity_count": len(r.get("wrong_identity", []))}
                   for r in numbering_results],
    )
    (out / "smoke-report.json").write_text(json.dumps(report, indent=2))
    results.commit()

    print("\n" + "=" * 70)
    print("WHAT THE SMOKE RUN MEASURED")
    print("=" * 70)
    print(f"card                     {card or 'NOT MEASURED (no nvidia-smi)'}"
          f"  ({total_mib} MiB)")
    print(f"peak card memory used    "
          f"{peak_mib if peak_mib is not None else 'NOT MEASURED'} MiB, sampled "
          f"once a second. An upper bound on the need: JAX can reserve most of "
          f"the card at start-up whether or not it uses it.")
    print(f"wall-clock, whole run    {wall_seconds:.0f} s")
    if trajectories:
        print(f"trajectories attempted   {trajectories} (rows of {trajectory_table})")
        print(f"wall-clock per trajectory {wall_seconds / trajectories:.1f} s")
    else:
        print("wall-clock per trajectory NOT COMPUTED: no table with 'traj' in "
              "its path. Row counts of every table, to pick the right one:")
        for name, rows in sorted(tables.items()):
            print(f"    {rows:6d}  {name}")
    print(f"complex size             "
          f"{', '.join(str(s) for s in sorted(set(sizes))) or 'NOT MEASURED'} "
          f"residues (target plus binder). Run time and memory scale with this; "
          f"it is what lets a PD-L1 figure be carried to EGFR.")

    print("\nCROPPING AND NUMBERING OF THE TARGET")
    if input_path is None:
        print("  NOT ANSWERED. The input structure could not be identified, so "
              "nothing was compared. Candidates found:")
        for c in input_candidates:
            print(f"    {c}")
        print("  Fix this before reading anything into the run: an unanswered "
              "question is not an answer of 'no change'.")
    elif not numbering_results:
        print("  NOT ANSWERED. No designed complex came back to compare.")
    else:
        statuses = sorted({r["status"] for r in numbering_results})
        for r in numbering_results:
            print(f"  {r['file']}: {target_numbering.describe(r)}")
        print(f"  Outcomes across {len(numbering_results)} complex(es): "
              f"{', '.join(statuses)}")
        if all(r["ok"] for r in numbering_results):
            print("  The target came back at the input's own numbers. "
                  "analysis/10 can translate it.")
        else:
            print("  The target did NOT come back at the input's numbers. "
                  "analysis/10 will refuse these. Do not start an EGFR campaign "
                  "before working out how to map them back.")
    return sorted(set(produced))


@app.function(timeout=60 * 10)
def inspect_cli():
    """Print BindCraft2's own help for each subcommand. No GPU, so it costs almost
    nothing. Written 2 October 2026 so that a design decision about scoring a
    candidate against something other than the main target rests on what the
    command says it does, not on a guess about its interface."""
    for sub in ("score", "rank", "filter", "fetch-weights"):
        _run(f"bindcraft {sub} --help")
    _run("bindcraft design --list-modalities || true")
    _run("bindcraft design --list-properties || true")


@app.function(timeout=60 * 10)
def shell(cmd: str):
    """Run one read-only shell command inside the BindCraft2 image, on the
    processor with no GPU. For reading BindCraft2's own documentation and source
    before configuring a paid run from it: modal run ...::shell --cmd "ls docs"."""
    _run(cmd)


@app.local_entrypoint()
def main():
    print("Run one of these explicitly rather than this entrypoint:")
    print("  modal run design/modal/bindcraft2_smoke.py::check")
    print("  modal run design/modal/bindcraft2_smoke.py::smoke")
