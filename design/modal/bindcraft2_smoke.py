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
here is a judgement rather than a requirement. UNVERIFIED either way. A100 40GB
is the default below because the folding step is the memory-hungry part and a
failure from running out of memory costs more time than the price difference
costs money. Drop to L4 with --gpu l4 once something is known to work.

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

import modal

APP_NAME = "bindcraft2-smoke"
VOLUME_NAME = "bindcraft2-results"
WEIGHTS_VOLUME = "bindcraft2-weights"

# Default card. See the note above on why this rather than something cheaper.
DEFAULT_GPU = "A100-40GB"

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
)

app = modal.App(APP_NAME)
results = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)
weights = modal.Volume.from_name(WEIGHTS_VOLUME, create_if_missing=True)


def _run(cmd, cwd=INSTALL_DIR):
    """Run a shell command inside the BindCraft2 environment and stream output."""
    import subprocess

    full = f"cd {cwd} && . .venv/bin/activate && {cmd}"
    print(f"\n$ {cmd}\n", flush=True)
    proc = subprocess.run(["bash", "-lc", full], capture_output=True, text=True)
    if proc.stdout:
        print(proc.stdout, flush=True)
    if proc.stderr:
        print("--- stderr ---", flush=True)
        print(proc.stderr, flush=True)
    return proc.returncode


@app.function(gpu=DEFAULT_GPU, timeout=60 * 30,
              volumes={"/weights": weights})
def check():
    """Cheapest possible test: does the GPU exist and can JAX see it?

    Run this before the design job. If JAX reports a processor device rather than
    a GPU one, the design run will work but take so long it looks broken, and the
    BindCraft2 troubleshooting notes call that out as the usual first fault.
    """
    _run("nvidia-smi", cwd="/")
    _run('python -c "'
         "import jax; print('jax', jax.__version__); "
         "print('devices', jax.devices()); "
         "print('default backend', jax.default_backend())"
         '"')
    # The subcommands and listing flags below are BindCraft2's own, read from its
    # cli.py on 1 October 2026. `--list-targets` is a flag of `design` rather
    # than of `bindcraft` itself, which the earlier version of this file had
    # wrong and which would have burned a GPU minute finding out.
    _run("bindcraft design --help")
    _run("bindcraft design --list-targets || true")
    _run("bindcraft design --list-settings | head -40 || true")
    print("\nIf 'default backend' above is not 'gpu', stop and fix that first.")


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

    rc = _run(f"bindcraft design {target}")
    print(f"\ndesign exit code {rc}")

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

    complexes = [p for p in out.rglob("*.cif")
                 if not p.stem.endswith("_monomer")]
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
    return sorted(set(produced))


@app.local_entrypoint()
def main():
    print("Run one of these explicitly rather than this entrypoint:")
    print("  modal run design/modal/bindcraft2_smoke.py::check")
    print("  modal run design/modal/bindcraft2_smoke.py::smoke")
