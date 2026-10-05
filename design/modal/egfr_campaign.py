"""egfr_campaign.py — run OUR EGFR campaign on a rented GPU, and measure what it costs.

PURPOSE. `bindcraft2_smoke.py` proved the pipeline works by running BindCraft2's
own shipped PD-L1 example, deliberately not our configuration, so that our
configuration could not be the cause of a first failure. That worked. This file is
the next step: the same pipeline pointed at our own trimmed EGFR fragment and our
own campaign file. It is a separate file rather than an edit to the smoke script,
so the run that is known to work stays available to go back to.

It does two things, in this order, and the order is not optional:

  validate  one trajectory, to find out whether the campaign file starts up at all
            against the real fragment. Costs cents. `docs/decisions-log.md`,
            Next actions item 2, names this as a prerequisite for item 4.
  pilot     150-200 trajectories, to measure the two rates that size the main run.

Abbreviations, expanded here because this file gets read on its own:
  EGFR  epidermal growth factor receptor, the protein we design against.
  GPU   graphics processing unit, the hardware the folding models need. Also
        called the card.
  JAX   the numerical library BindCraft2's folding models run on.
  PDB   the Protein Data Bank, the public archive of measured three-dimensional
        structures, and also the older of the two file formats it serves them in.
  mmCIF the newer structure file format, which is what BindCraft2 writes.
  FASTA the plain-text format for a sequence with no structure attached.
  aa    amino acids, the unit protein length is counted in.
  His   histidine, the one amino acid that changes charge between pH 7.4 and
        pH 6.5. A His tag is a short run of them added to a protein so it can be
        purified; the assay's target carries one and is expected to keep it.

WHAT THE TWO RATES ARE, AND WHY THEY ARE THE THING WORTH PAYING FOR
-------------------------------------------------------------------
Sizing the main run needs two numbers that can only come from a real run against
our own target. Neither is known. The campaign file currently asks for up to 2,000
trajectories and 200 accepted designs, and `design/configs/egfr-domain3-h370.md`
says plainly that those counts were never derived from a measured runtime.

  candidates per card-hour
      How many designs clear BindCraft2's own filters per hour of card time. A
      trajectory is one attempt at a binder; most attempts are discarded. The
      PD-L1 smoke run accepted 10 of 30 attempts (33%) at 184.7 seconds an
      attempt, but that was a 115-residue target on an A100 40GB card, and ours
      is a 171-residue fragment on an L4. Both differences push the number and
      neither has been measured.

  charge-pair survival rate
      Of the designs BindCraft2 accepts, the share that reach the floor of three
      correct charge pairs that `analysis/10_charge_pair_filter.py` requires.
      BindCraft2 has no interest in pH -- it reports binder charge at a hard-coded
      pH 7.4 as a suppressed readout and has no pH term in its objective -- so it
      is not selecting for the thing we need, and this rate is the cost of that
      indifference. It has never been measured on real output because no real
      output has existed.

Multiplied together they give designs-that-pass-our-filter per card-hour, which is
what converts a budget into a campaign size. The main run is not launched from
here: that decision goes back to the owner, per `CLAUDE.md`.

THE FILES THIS NEEDS INSIDE THE CONTAINER, AND HOW THEY GET THERE
-----------------------------------------------------------------
The smoke run needed nothing of ours, because it ran BindCraft2's own example.
This run needs four of our files, none of which is in the image:

    data/structures/6aru_domain3.pdb          the trimmed target, 171 residues
    data/sequences/his-tag-offtarget.fasta    the His tag, as an off-target
    design/configs/egfr-domain3-h370.json     the campaign
    design/configs/egfr-domain3-h370-notag.json  the same, plus the tag off-target

The mechanism is `modal.Image.add_local_file(local, remote)`, checked against the
installed client (Modal 1.6.0) rather than recalled: `copy_local_file` and
`copy_local_dir`, which earlier versions had, are gone from this one. The default
`copy=False` adds the files when a container starts instead of baking them into an
image layer, which is what we want -- editing a campaign file then costs nothing,
where `copy=True` would rebuild the roughly 20 GB image behind it.

**The remote paths mirror the repository layout, and that is load-bearing.**
BindCraft2 resolves `target_path` relative to the campaign file's own directory
(its `docs/source/reference.md`, line 41, read 2 October 2026). Our campaign files
carry `"target_path": "../../data/structures/6aru_domain3.pdb"`. So the config has
to sit two directories below a root that also holds `data/`, exactly as it does on
disk. Flattening the four files into one directory would break the lookup. The
alternative -- rewriting the paths to absolute ones inside the container -- was
rejected because it would mean the file BindCraft2 reads is not the file in the
repository, and the committed config is generated by `analysis/09_trim_target.py`
and documented as not to be hand-edited.

WHAT WAS CHECKED BEFORE SPENDING ANYTHING, AND WHAT CAME BACK
-------------------------------------------------------------
Read from BindCraft2's own `docs/source/reference.md` and `README.md` inside the
image on 2 October 2026, using the no-GPU `shell` helper in `bindcraft2_smoke.py`:

  - `targets[].target_path` takes "PDB, mmCIF or FASTA -- supply a structure or
    sequence target" (reference.md line 68). So handing it a FASTA file for the
    His tag is a supported shape and not a hopeful guess. `docs/decisions-log.md`
    lists detargeting against a sequence target under Unverified; this moves the
    question from "is the shape even accepted" to "what does it do to the run",
    which only running it answers.
  - `crop_fasta_sequence` takes a window such as `[10,40]`, and "`false` uses the
    full sequence" (line 80). The notag config sets `false`, which is right for a
    ten-residue tag there is no reason to crop.
  - `termini_accessible` is **not** among the 134 names `--list-settings` prints,
    which was alarming, because `docs/decisions-log.md` records it as set in both
    campaign files and an invented key could make BindCraft2 refuse the file.
    It is valid. BindCraft2's own `README.md` tabulates `"termini_accessible":
    true` as the campaign-file spelling of the `--termini-accessible` flag, and
    `settings/core/reference.json` carries the key. `--list-settings` lists what
    `--set KEY=VALUE` accepts, which is a narrower set than the keys a campaign
    file may hold. Recorded because the next person to check will have the same
    scare.
  - `--set max_trajectories=N` and `--set project_folder=PATH` are both real, from
    `bindcraft design --help`. The pilot is bounded with the first of these rather
    than by editing the committed JSON, whose counts describe the eventual main
    run and are still unconfirmed.

HOW THE CARD TIME IS BOUNDED
----------------------------
Two independent limits, because either alone has a failure mode.

`--set max_trajectories=N` stops the run when it has spent N attempts. That is the
limit that should fire, and the run ends by itself.

The Modal function timeout is the backstop for the case where a trajectory hangs
or the rate is far worse than measured. It is set from a spend ceiling in dollars
at the L4 price, so the ceiling is stated in the unit the budget is in rather than
in seconds. At $0.80 an hour a dollar is 4,500 seconds of card time.

A timeout that killed the container would lose every file the run had written, so
the project folder is copied into the results volume and committed every five
minutes by a background thread. A run cut short then still yields both rates from
however many trajectories finished, because both are rates. Writing the project
folder straight onto the volume was the alternative and was rejected: a Modal
volume is not a complete POSIX filesystem, and BindCraft2 writes and rewrites many
small files and tables as it goes.

THE OUTPUT-COPY BUG THIS FILE DOES NOT REPEAT
---------------------------------------------
`smoke()` globbed every `*.cif` under the whole BindCraft2 install directory, which
swept up four template scaffolds BindCraft2 ships for other binder formats
(`scaffolds/ARP.cif`, `Fab.cif`, `VHH.cif`, `scFv.cif`). The numbering check then
correctly reported that a nanobody scaffold does not resemble PD-L1, and the run's
printed verdict said the target "did NOT come back at the input's numbers" and
warned against starting an EGFR campaign. That was wrong: all 91 real designed
complexes read identical. Fixed 2 October 2026 and recorded in
`docs/decisions-log.md`.

This file copies the project folder alone, and the numbering check looks only at
files under `1_Trajectories/`, `2_Refolded/` and `3_Ranked/`. Both restrictions,
rather than one, because they fail differently: a wrong copy loses files we paid
for, and a wrong check list produces a confident false verdict.

HOW TO RUN
----------
    ./.venv/bin/modal run design/modal/egfr_campaign.py::validate
        One trajectory against the notag campaign. Prints whether BindCraft2
        accepted the campaign file, whether it accepted the His tag as a sequence
        off-target, and the seconds one trajectory took on this card for our
        fragment -- which is the figure the pilot is sized from.

    ./.venv/bin/modal run design/modal/egfr_campaign.py::validate --config plain
        The same against the campaign without the tag off-target, as the fallback
        if the notag one does not start.

    ./.venv/bin/modal run --detach design/modal/egfr_campaign.py::pilot \
        --trajectories 150 --spend-ceiling-usd 8.0
        The pilot. Detached, because it outlasts a terminal. `--detach` is what
        the smoke run used.

    ./.venv/bin/modal volume get bindcraft2-results /egfr-pilot ./results/candidates/
        Bring the output back for analysis/10.

    ./.venv/bin/python design/modal/egfr_campaign.py
        Self-test of the arithmetic and the file-layout helpers, locally, with no
        Modal account and no card. See `_self_test`.
"""

import os

import modal

APP_NAME = "egfr-campaign"
VOLUME_NAME = "bindcraft2-results"
WEIGHTS_VOLUME = "bindcraft2-weights"

# Same card as the smoke script's default, and for the same measured reason: that
# run peaked at 18,218 MiB on an A100 40GB, which fits an L4's 24,564 MiB. The
# figure is an upper bound, since JAX can reserve more of the card than it uses.
# No speed ratio between the two cards has ever been measured, and measuring it is
# part of what the pilot is for, so the per-trajectory times from the A100 smoke
# run are not assumed to carry over. There is no runtime flag for the card:
# `gpu=` is fixed when the function is decorated, so switching back to an A100
# means editing this constant.
DEFAULT_GPU = "L4"

# BindCraft2 is installed here by the image build. The path is fixed because its
# own documentation warns that moving the directory breaks the `bindcraft`
# command, the install being editable.
INSTALL_DIR = "/opt/BindCraft2"

# Our files are mirrored under here, keeping the repository's own layout so that
# the relative `target_path` in a campaign file resolves. See the module docstring.
WORK_DIR = "/work"

# L4 price per second, modal.com/pricing read 1 October 2026: $0.80 an hour.
# Used only to turn a spend ceiling in dollars into a timeout in seconds, so the
# ceiling can be stated in the unit the $30 monthly allowance is in.
L4_USD_PER_SECOND = 0.80 / 3600.0

# The two campaign files, by the short name the command line takes. The notag one
# is tried first: a design generated without the His tag as an off-target cannot
# be screened for the tag afterwards, so if it runs, it is the one to use.
CONFIGS = {
    "notag": "design/configs/egfr-domain3-h370-notag.json",
    "plain": "design/configs/egfr-domain3-h370.json",
    # A diagnostic, not a campaign. Same settings as `notag` with one difference,
    # the target fragment extended from 310-480 to 310-499 so the disulfide the cut
    # severed has its partner back. `design/configs/diagnostic/README.md` says what
    # question it answers and why it is kept apart from the generated files.
    "ext499": "design/configs/diagnostic/ext499-notag.json",
    # The committed fragment and settings, with the hotspot list cut from eight
    # anchors to the tight four on the N352 side. The one variable changed, so an
    # improvement in interface confidence can be attributed to it.
    "tight4": "design/configs/diagnostic/tight4-notag.json",
    # Round one, 3 October 2026: the tight4 hotspots (E344, H358, D368, H370; config
    # numbers A320, A334, A344, A346), same committed 310-480 fragment, same His-tag
    # off-target, with the one variable between the two being the binder length band.
    # `batch_r1_short` and `batch_r1_long` below say why.
    "r1-short": "design/configs/diagnostic/r1-short-notag.json",
    "r1-long": "design/configs/diagnostic/r1-long-notag.json",
}

# The interface-confidence floors, lowered. BindCraft2's defaults are
# `min_iptm_<stage>` of 0.5 at anneal, harden and mutate and 0.7 at the final
# check. Eleven trajectories against our target produced `i_pTM` readings from
# 0.20 to 0.62, so the default floors stop everything before a structure is
# written and our own charge-pair filter never sees a candidate.
#
# Lowering them is a change to a rejection threshold and was authorised by the
# owner on 2 October 2026, on a stated ground: BindCraft2's 0.7 is its own
# conservative default and not the competition's bar. The organisers said on
# 2 October that in silico scoring is "only a smaller component in the overall
# selection criteria" and that worse confidence on a rigid or obstructed epitope
# will not be penalised much (`docs/competition-qa-log.md`), and this project's own
# ranking leads on charge-pair count with `i_pDAE` only as a tie-break. So a design
# at 0.5 is one we can still judge on the criterion that matters to us.
#
# What this gives up, stated rather than buried: these are designs the predictor
# does not confidently believe bind, and a low `i_pTM` can mean exactly that. The
# floors are a diagnostic setting, not a proposal for the main run.
LOWERED_IPTM_FLOORS = {
    "min_iptm_anneal": 0.35,
    "min_iptm_harden": 0.35,
    "min_iptm_mutate": 0.35,
    "min_iptm_final": 0.5,
}

# Round one's settings: the lowered floors above, plus `save_design_sequences`
# turned on. Added 3 October 2026, after `egfr-tight4-tight4` reached the final
# stage on two of its six trajectories (i_pTM 0.51 and 0.61 in the optimisation
# losses, our highest readings against this target to date) and still left no
# recoverable sequence or structure on disk or on the Modal volume -- every
# `1_Trajectories/<design>/` folder held only `<design>_losses.csv`, a table of
# scalar metrics, because `save_design_sequences` defaults to false. Without it,
# a trajectory that clears every gate still writes nothing `analysis/10` or
# `analysis/18` can read. `save_failed_trajectories` and `save_failed_refolds`
# were already true and are left as they are.
ROUND1_SETTINGS = dict(LOWERED_IPTM_FLOORS, save_design_sequences=True)

# ---------------------------------------------------------------------------
# Overriding BindCraft2's worker-memory estimate, from outside the package
# ---------------------------------------------------------------------------
#
# BindCraft2 decides how many design workers share a card in
# `bindcraft/design_workers.py`, in `design_workers_per_gpu` (lines 61-70). Lines
# 68-69 clamp the count by an estimate of how much card memory a worker needs, and
# that clamp is unconditional -- it overrides an explicitly requested worker count
# just as it overrides the automatic one. For our short band the estimate is
# 11.78 GB a worker, which on the L4's 22.49 GB leaves room for exactly one once
# 4 GB of headroom is held back. No campaign setting and no environment variable
# gets past it: `max_workers_per_gpu` and `design_workers` can only lower the
# count, and `BINDCRAFT_WORKERS_PER_GPU` feeds the same clamped variable.
#
# The estimate carries a safety multiplier of 2.0. Round one then measured a real
# peak of 8,894 MiB, which is 8.69 GB, so two real workers need about 17.4 GB of
# the 22.49 GB card and leave roughly 5 GB spare. The card fits two workers and the
# estimate does not. The full derivation, with source line numbers, is at the end
# of Pipeline status in `docs/decisions-log.md`.
#
# So the multiplier is recalibrated against that measured peak, at run time, from
# here. Three constraints shaped how:
#
#   - The vendored package is not edited and the image is not rebuilt. The patch is
#     a `sitecustomize.py` written into the container at run time and put on
#     PYTHONPATH, which Python imports at interpreter start-up. BindCraft2's own
#     virtual environment has no `sitecustomize` of its own, so nothing is shadowed
#     (checked 5 October 2026).
#   - `bindcraft design` runs as a subprocess, so patching this process would do
#     nothing. PYTHONPATH reaches the subprocess and the design workers it spawns.
#   - The patch is applied when `bindcraft.design_workers` is first imported, by a
#     finder on `sys.meta_path`, rather than by importing that module at start-up.
#     Importing it eagerly would pull JAX into every Python process the container
#     starts, including ones that never touch a card.
#
# 1.5 is chosen so that two workers are planned and three are not: at 256 residues
# it gives 8.84 GB, and 18.49 // 8.84 = 2, where three would need 6.16 GB or below.
# It sits just above the 1.475 the measured peak implies, so it stays marginally
# conservative. This is a prediction from one band's measured peak, not a measured
# result, which is why the patch reports what it computed and the caller checks the
# run's own peak memory against a 20,000 MiB ceiling.
CALIBRATED_MEMORY_SAFETY_FACTOR = 1.5

WORKER_MEMORY_PATCH = '''"""Written into the container at run time by design/modal/egfr_campaign.py.

Recalibrates BindCraft2's design-worker memory estimate against a measured peak, so
that the worker count it plans reflects what the card holds rather than a 2.0 safety
multiplier. Explained where it is generated; derived in docs/decisions-log.md.

Nothing here edits BindCraft2. The module is imported by Python at start-up because
it is named sitecustomize and sits on PYTHONPATH, and it patches the target module
only when that module is first imported.
"""
import os
import sys
from importlib.machinery import PathFinder

TARGET = "bindcraft.design_workers"
FACTOR = float(os.environ.get("BC2_CALIBRATED_SAFETY_FACTOR", "1.5"))


def _apply(module):
    original_estimate = module.estimate_design_memory_gb
    resident_gb = module.DESIGN_MODEL_RESIDENT_GB
    bytes_per_pair = module.DESIGN_ACTIVATION_BYTES_PER_RESIDUE_PAIR

    def calibrated(residue_count):
        return FACTOR * (resident_gb + bytes_per_pair * int(residue_count) ** 2 / 1e9)

    module.estimate_design_memory_gb = calibrated

    original_per_gpu = module.design_workers_per_gpu

    def reported(settings, free_gb, residue_count):
        workers = original_per_gpu(settings, free_gb, residue_count)
        count = int(residue_count or 0)
        print("[worker-memory-patch] free %.2f GB, headroom %.1f GB, residues %d, "
              "BindCraft2 estimate %.2f GB, calibrated %.2f GB, workers per card %d"
              % (free_gb, module.GPU_MEMORY_HEADROOM_GB, count,
                 original_estimate(count), calibrated(count), workers), flush=True)
        return workers

    module.design_workers_per_gpu = reported

    original_plan = module.plan_design_workers

    def planned(settings, residue_count=None, trajectory_budget=None):
        plan = original_plan(settings, residue_count, trajectory_budget)
        print("[worker-memory-patch] planned %d worker(s); memory fractions %s"
              % (len(plan), [worker.get("memory_fraction") for worker in plan]),
              flush=True)
        return plan

    module.plan_design_workers = planned


class _PatchingLoader:
    def __init__(self, inner):
        self._inner = inner

    def create_module(self, spec):
        return self._inner.create_module(spec)

    def exec_module(self, module):
        self._inner.exec_module(module)
        try:
            _apply(module)
            print("[worker-memory-patch] applied, safety multiplier %.3f" % FACTOR,
                  flush=True)
        except Exception as exc:
            print("[worker-memory-patch] FAILED to apply: %r" % (exc,), flush=True)


class _PatchFinder(PathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname != TARGET:
            return None
        spec = super().find_spec(fullname, path, target)
        if spec is None or spec.loader is None:
            return None
        spec.loader = _PatchingLoader(spec.loader)
        return spec


try:
    if not any(isinstance(finder, _PatchFinder) for finder in sys.meta_path):
        sys.meta_path.insert(0, _PatchFinder())
except Exception as exc:
    print("[worker-memory-patch] FAILED to install: %r" % (exc,), flush=True)
'''

PATCH_DIR = "/tmp/bc2-worker-memory-patch"


def write_worker_memory_patch(factor=CALIBRATED_MEMORY_SAFETY_FACTOR):
    """Write the sitecustomize patch into the container and return its directory.

    Returns the shell prefix that puts it on PYTHONPATH for one command, rather than
    exporting it for the whole container, so that only `bindcraft design` and the
    workers it spawns are affected and the read-only helpers are not.
    """
    import pathlib

    directory = pathlib.Path(PATCH_DIR)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "sitecustomize.py").write_text(WORKER_MEMORY_PATCH)
    print(f"worker-memory patch written to {directory}/sitecustomize.py, "
          f"safety multiplier {factor}")
    return f"PYTHONPATH={directory} BC2_CALIBRATED_SAFETY_FACTOR={factor} "


# BindCraft2's own output subdirectories, inside the project folder. The numbering
# check looks nowhere else; see the module docstring on the bug that caused.
OUTPUT_DIRS = ("1_Trajectories", "2_Refolded", "3_Ranked")

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.join(_HERE, "..", "..")


def _local(*parts):
    """A path in the repository, from this file's location.

    Plain string joining rather than pathlib arithmetic, because this module is
    imported again inside the container where the layout differs and
    parent-directory arithmetic on a different tree can silently produce a path
    that exists but is the wrong one.
    """
    return os.path.join(_REPO, *parts)


image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.8.1-devel-ubuntu24.04", add_python="3.12"
    )
    .apt_install("git", "wget", "curl", "build-essential", "ca-certificates")
    .run_commands(
        f"git clone --depth 1 https://github.com/PacesaLab/BindCraft2.git {INSTALL_DIR}",
        f"cd {INSTALL_DIR} && bash install.sh",
    )
    .env({"BINDCRAFT_HOME": INSTALL_DIR})
    # Our own files, added last so that editing one does not rebuild the install
    # above. `copy=False` is the default and is what we want: these arrive when a
    # container starts rather than becoming an image layer.
    .add_local_file(_local("data", "structures", "6aru_domain3.pdb"),
                    f"{WORK_DIR}/data/structures/6aru_domain3.pdb")
    .add_local_file(_local("data", "sequences", "his-tag-offtarget.fasta"),
                    f"{WORK_DIR}/data/sequences/his-tag-offtarget.fasta")
    .add_local_file(_local("design", "configs", "egfr-domain3-h370.json"),
                    f"{WORK_DIR}/design/configs/egfr-domain3-h370.json")
    .add_local_file(_local("design", "configs", "egfr-domain3-h370-notag.json"),
                    f"{WORK_DIR}/design/configs/egfr-domain3-h370-notag.json")
    # The diagnostic fragment and its config. Mounted alongside rather than in
    # place of the committed pair, so a run can be pointed at either without a
    # rebuild and the two can be compared in the same session.
    .add_local_file(_local("data", "structures", "6aru_domain3_ext499.pdb"),
                    f"{WORK_DIR}/data/structures/6aru_domain3_ext499.pdb")
    .add_local_file(_local("design", "configs", "diagnostic", "ext499-notag.json"),
                    f"{WORK_DIR}/design/configs/diagnostic/ext499-notag.json")
    .add_local_file(_local("design", "configs", "diagnostic", "tight4-notag.json"),
                    f"{WORK_DIR}/design/configs/diagnostic/tight4-notag.json")
    # Round one's two length-band configs. Same target structure and His-tag
    # off-target as tight4, already mounted above; no new data file needed.
    .add_local_file(_local("design", "configs", "diagnostic", "r1-short-notag.json"),
                    f"{WORK_DIR}/design/configs/diagnostic/r1-short-notag.json")
    .add_local_file(_local("design", "configs", "diagnostic", "r1-long-notag.json"),
                    f"{WORK_DIR}/design/configs/diagnostic/r1-long-notag.json")
    .add_local_file(_local("analysis", "target_numbering.py"),
                    "/opt/target_numbering.py")
)

app = modal.App(APP_NAME, image=image)
results = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)
weights = modal.Volume.from_name(WEIGHTS_VOLUME, create_if_missing=True)


# ---------------------------------------------------------------------------
# Helpers that do not need Modal, a card or a network, so `_self_test` can
# exercise them. Everything with arithmetic in it lives here for that reason.
# ---------------------------------------------------------------------------

def timeout_for_spend(ceiling_usd, usd_per_second=L4_USD_PER_SECOND):
    """Seconds of card time a dollar ceiling buys, as a whole number.

    Stated in dollars at the call site because the budget is $30 a month and
    about $3.23 of it went on the PD-L1 smoke run, so what is left is tracked in
    dollars and not in seconds.
    """
    return int(ceiling_usd / usd_per_second)


def cost_usd(seconds, usd_per_second=L4_USD_PER_SECOND):
    """What a run of this many seconds cost on this card.

    Card time only. Modal bills image builds as processor time, which is cheaper
    and is not counted here, so this understates the bill slightly rather than
    overstating it.
    """
    return seconds * usd_per_second


def rates(wall_seconds, trajectories, accepted):
    """The two throughput figures, plus the acceptance share, or None where unknown.

    `accepted` is designs that cleared BindCraft2's own filters. The charge-pair
    survival rate is not here: it comes from `analysis/10_charge_pair_filter.py`
    reading the structures afterwards, and inventing a placeholder for it would
    put a number in a report that nothing computed.

    Returns a dict. A field is None when its input is missing or zero, rather than
    0.0, because "no trajectory table was written" and "no trajectory was
    accepted" are different outcomes and a zero would merge them.
    """
    hours = wall_seconds / 3600.0 if wall_seconds else None
    return dict(
        card_hours=round(hours, 3) if hours else None,
        seconds_per_trajectory=(round(wall_seconds / trajectories, 1)
                                if wall_seconds and trajectories else None),
        acceptance_rate=(round(accepted / trajectories, 4)
                         if trajectories and accepted is not None else None),
        candidates_per_card_hour=(round(accepted / hours, 2)
                                  if hours and accepted else None),
    )


def structure_target(spec):
    """The campaign's structure target, and the off-targets, kept apart.

    A campaign may hold several targets and ours holds two: the EGFR fragment as a
    structure, and the His tag as a sequence in a FASTA file. The numbering check
    has to be run against the structure one. The smoke script's equivalent took
    whichever entry came last, which here would pick the tag and compare every
    returned complex against a ten-residue sequence with no coordinates.

    Returns (structure entry or None, list of the others). A structure entry is one
    whose `target_path` ends in a structure suffix; `.fasta` and `.fa` are
    sequence targets.
    """
    structure, others = None, []
    for entry in spec.get("targets") or []:
        path = str(entry.get("target_path") or "")
        if path.lower().endswith((".pdb", ".cif", ".mmcif")):
            if structure is None:
                structure = entry
            else:
                others.append(entry)
        else:
            others.append(entry)
    return structure, others


def count_accepted(project_dir):
    """Accepted designs, counted two ways, with where each count came from.

    Two ways because they can disagree and the disagreement is informative. The
    ranked table `3_Ranked/!_Ranked.csv` has a row per accepted design; the
    structures in `3_Ranked/` are the files `analysis/10` actually reads, ignoring
    the `*_monomer.cif` files that hold the binder with no target. If the table
    says more than the directory holds, something was not written and the
    candidates we paid for are not all there.

    Returns a dict, with None where the thing being counted is absent rather than
    empty. `project_dir` is a pathlib path.
    """
    ranked_rows, ranked_table = None, None
    for path in project_dir.rglob("*.csv"):
        if path.name.startswith("!_Ranked") or "ranked" in path.name.lower():
            try:
                lines = path.read_text(errors="replace").splitlines()
            except OSError:
                continue
            ranked_rows = max(0, len(lines) - 1)
            ranked_table = str(path.relative_to(project_dir))
            break
    complexes = [p for p in project_dir.rglob("*.cif")
                 if not p.stem.endswith("_monomer")
                 and "3_Ranked" in p.parts]
    return dict(ranked_table_rows=ranked_rows, ranked_table=ranked_table,
                ranked_complexes=len(complexes) if complexes else 0)


def numbering_candidates(project_dir):
    """Complexes the numbering check should read: the run's own output, nothing else.

    Restricted to `OUTPUT_DIRS` and excluding `*_monomer.cif`. This is the
    restriction the smoke script lacked, which let four of BindCraft2's shipped
    template scaffolds into the check and made its verdict say the target had not
    come back when every real complex said it had.
    """
    return sorted(p for p in project_dir.rglob("*.cif")
                  if not p.stem.endswith("_monomer")
                  and any(d in p.parts for d in OUTPUT_DIRS))


# ---------------------------------------------------------------------------
# Instrumentation carried over from bindcraft2_smoke.py. Unchanged in behaviour,
# because it is the part that is already known to work on real output.
# ---------------------------------------------------------------------------

def _run(cmd, cwd=INSTALL_DIR, venv=True, want_output=False, stream=False):
    """Run a shell command inside BindCraft2's environment and print what it said.

    `venv=False` skips activating BindCraft2's virtual environment, whose
    activation path is relative to `cwd`, for commands that are not part of it
    such as nvidia-smi. `want_output=True` returns (exit code, stdout).

    `stream=True` prints each line as it arrives instead of holding everything
    until the command returns, and still returns the output when asked for it.

    The streaming option exists because of what the buffered version cost on the
    PD-L1 smoke run. `docs/decisions-log.md` records it: "its console output only
    appears when the command finishes, because `_run` captures it, so progress is
    not visible", and that run was then left unattended and its final report never
    collected. For a one-trajectory check the difference does not matter. For a
    pilot of 150 trajectories it is the difference between noticing a stall in ten
    minutes and noticing it in six hours, and the gap is paid for in card time.

    A hang and healthy progress look identical in silence, which is the specific
    reason this is not merely a convenience: silence is not evidence that the run
    is working.

    Output is merged into one stream rather than kept as separate standard output
    and standard error, because the interleaving is what makes a progress line
    readable next to the warning that preceded it. The cost is that the two cannot
    be told apart afterwards, which is why the buffered path is kept for the short
    commands where separating them is the more useful property.
    """
    import subprocess

    activate = ". .venv/bin/activate && " if venv else ""
    full = f"cd {cwd} && {activate}{cmd}"
    print(f"\n$ {cmd}\n", flush=True)

    if not stream:
        proc = subprocess.run(["bash", "-lc", full], capture_output=True, text=True)
        if proc.stdout:
            print(proc.stdout, flush=True)
        if proc.stderr:
            print("--- stderr ---", flush=True)
            print(proc.stderr, flush=True)
        if want_output:
            return proc.returncode, proc.stdout
        return proc.returncode

    collected = []
    proc = subprocess.Popen(["bash", "-lc", full], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in proc.stdout:
        collected.append(line)
        print(line.rstrip("\n"), flush=True)
    proc.wait()
    if want_output:
        return proc.returncode, "".join(collected)
    return proc.returncode


def _start_gpu_poll(path):
    """Sample the card once a second into a file while the design runs.

    Returns the process, or None if nvidia-smi is missing, in which case the
    report says memory was not measured rather than leaving the line out.
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

    Peak is the largest figure seen at one-second sampling, so a spike between two
    samples is missed and the true peak is at least this. The figure is the card's
    total used memory, which JAX may inflate by reserving most of the card at
    start-up whether it needs it or not, so it bounds the requirement from above
    and does not state it.
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


def _count_trajectories(project):
    """Trajectories attempted, from whatever table the run wrote, and which table.

    BindCraft2 keeps attempting until enough designs pass, so attempts are not
    what was requested. Taken as the row count of a table whose path mentions
    trajectories; None when there is no such table, in which case every table's
    row count is printed so the right one can be chosen by hand rather than
    guessed.
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


def _snapshot(project_dir, destination, every_seconds=300):
    """Copy the project folder into the results volume every few minutes, forever.

    Runs on a daemon thread. The point is that the Modal function timeout is a
    backstop that kills the container, and a killed container loses everything on
    its local disk including a multi-hour run's output. With this, a run cut short
    still has its trajectories, its tables and its structures on the volume, and
    both rates can be computed from however much finished.

    Errors are printed and swallowed. A failed snapshot must not end a run that is
    otherwise working, and the final copy in `_campaign` is the one that matters.
    """
    import shutil
    import threading
    import time

    def loop():
        while True:
            time.sleep(every_seconds)
            try:
                if project_dir.exists():
                    shutil.copytree(project_dir, destination, dirs_exist_ok=True)
                    results.commit()
                    print(f"[snapshot] {destination} committed", flush=True)
            except Exception as exc:                          # noqa: BLE001
                print(f"[snapshot] failed, continuing: {exc}", flush=True)

    thread = threading.Thread(target=loop, daemon=True)
    thread.start()
    return thread


# ---------------------------------------------------------------------------
# The campaign itself
# ---------------------------------------------------------------------------

def _campaign(config_key, trajectories, run_name, final_designs=None,
              extra_settings=None, workers=None):
    """Run one campaign, measure it, and write a report beside its output.

    Shared by `validate` and `pilot`, which differ only in how many trajectories
    they spend and how long they are allowed to take. One implementation so that
    the thing validated is the thing the pilot then runs.
    """
    import json
    import pathlib
    import shutil
    import sys
    import time

    config_relative = CONFIGS[config_key]
    config_path = pathlib.Path(WORK_DIR) / config_relative

    print("=" * 72)
    print(f"CAMPAIGN  {run_name}")
    print("=" * 72)
    print(f"campaign file   {config_path}")
    print(f"trajectories    {trajectories} (via --set max_trajectories)")

    # Does the campaign file exist, and does it name files that exist? Checked
    # before any card time is spent, because a missing mount is the most likely
    # way this run differs from the smoke run, and finding out from BindCraft2's
    # own error message costs a container start-up.
    if not config_path.exists():
        raise RuntimeError(
            f"the campaign file is not in the container at {config_path}. "
            "The image's add_local_file calls and CONFIGS disagree.")
    spec = json.loads(config_path.read_text())
    structure, off_targets = structure_target(spec)
    if structure is None:
        raise RuntimeError(f"{config_relative} names no structure target")

    input_path = (config_path.parent / structure["target_path"]).resolve()
    print(f"structure target {structure.get('name')}  ->  {input_path}")
    if not input_path.exists():
        raise RuntimeError(
            f"the target structure is not in the container at {input_path}. "
            "target_path resolves from the campaign file's own directory, so "
            "the mounted layout must mirror the repository's.")
    for entry in off_targets:
        other = (config_path.parent / entry["target_path"]).resolve()
        print(f"off-target       {entry.get('name')} weight "
              f"{entry.get('weight')}  ->  {other}"
              f"{'' if other.exists() else '   MISSING'}")
        if not other.exists():
            raise RuntimeError(f"off-target file missing in container: {other}")
    print(f"hotspots         {structure.get('hotspots')}")
    print("                 (these are the trimmed file's own numbers; add 24 "
          "for our UniProt numbering)")

    rc = _run("bindcraft fetch-weights || true")
    print(f"fetch-weights exit code {rc}")

    # The project folder goes on the container's local disk and is snapshotted to
    # the volume, rather than being written onto the volume directly. See the
    # module docstring on why.
    project = pathlib.Path("/tmp") / run_name
    volume_out = pathlib.Path("/results") / run_name
    volume_out.mkdir(parents=True, exist_ok=True)

    overrides = [f"--set 'max_trajectories={trajectories}'",
                 f"--set 'project_folder={project}'"]
    # `workers_per_gpu` alone does not get two workers onto an L4: BindCraft2 clamps
    # whatever is asked for by its own memory estimate. The patch recalibrates that
    # estimate against round one's measured peak; without it this setting is ignored.
    # See WORKER_MEMORY_PATCH above and docs/decisions-log.md, end of Pipeline status.
    patch_prefix = ""
    if workers:
        overrides.append(f"--set 'workers_per_gpu={workers}'")
        print(f"workers requested  {workers}")
        if workers > 1:
            patch_prefix = write_worker_memory_patch()
    if final_designs is not None:
        overrides.append(f"--set 'number_of_final_designs={final_designs}'")
    for key, value in sorted((extra_settings or {}).items()):
        overrides.append(f"--set '{key}={value}'")
        print(f"override         {key}={value}")

    poll_path = "/tmp/gpu-poll.csv"
    poller = _start_gpu_poll(poll_path)
    _snapshot(project, volume_out)

    started = time.monotonic()
    # Streamed, so a stalled run is visible while it is still cheap to stop. See
    # `_run`.
    rc = _run(patch_prefix + f"bindcraft design {config_path} " + " ".join(overrides),
              stream=True)
    wall_seconds = time.monotonic() - started
    if poller is not None:
        poller.terminate()

    print(f"\ndesign exit code {rc}, wall-clock {wall_seconds:.0f} s, "
          f"card time about ${cost_usd(wall_seconds):.2f}")

    # Final copy of the project folder, and nothing else. The smoke script copied
    # the whole install tree; see the module docstring.
    if project.exists():
        shutil.copytree(project, volume_out, dirs_exist_ok=True)
    results.commit()

    card, peak_mib, total_mib = _read_gpu_poll(poll_path)
    attempted, trajectory_table, tables = _count_trajectories(volume_out)
    accepted = count_accepted(volume_out)
    complexes = numbering_candidates(volume_out)

    sys.path.insert(0, "/opt")
    import target_numbering

    numbering, sizes = [], []
    for path in complexes:
        result = target_numbering.check_complex(input_path, path)
        result["file"] = str(path.relative_to(volume_out))
        numbering.append(result)
        sizes.append(result["complex_residues"])

    accepted_count = accepted["ranked_complexes"] or accepted["ranked_table_rows"]
    throughput = rates(wall_seconds, attempted, accepted_count)

    report = dict(
        run_name=run_name, config=config_relative,
        requested_trajectories=trajectories,
        extra_settings=dict(extra_settings or {}),
        design_exit_code=rc,
        wall_clock_seconds=round(wall_seconds, 1),
        estimated_cost_usd=round(cost_usd(wall_seconds), 2),
        card=card, card_memory_mib=total_mib, peak_memory_used_mib=peak_mib,
        trajectories_attempted=attempted, trajectory_table=trajectory_table,
        table_row_counts=tables,
        accepted=accepted,
        input_structure=str(input_path),
        off_targets=[{"name": e.get("name"), "path": e.get("target_path"),
                      "weight": e.get("weight")} for e in off_targets],
        complexes_checked=len(complexes),
        complex_residues=sorted(set(sizes)),
        **throughput,
        numbering=[{k: v for k, v in r.items()
                    if k not in ("missing_from_output", "not_in_input",
                                 "wrong_identity")}
                   | {"missing_count": len(r.get("missing_from_output", [])),
                      "not_in_input_count": len(r.get("not_in_input", [])),
                      "wrong_identity_count": len(r.get("wrong_identity", []))}
                   for r in numbering],
    )
    (volume_out / "campaign-report.json").write_text(json.dumps(report, indent=2))
    results.commit()

    # ---- the report ------------------------------------------------------
    print("\n" + "=" * 72)
    print(f"WHAT {run_name} MEASURED")
    print("=" * 72)
    print(f"startup                  {'OK' if rc == 0 else f'FAILED, exit {rc}'}")
    print(f"card                     {card or 'NOT MEASURED (no nvidia-smi)'} "
          f"({total_mib} MiB)")
    print(f"peak card memory used    "
          f"{peak_mib if peak_mib is not None else 'NOT MEASURED'} MiB, sampled "
          f"once a second. An upper bound on what is needed: JAX can reserve more "
          f"of the card than it uses.")
    print(f"wall-clock, whole run    {wall_seconds:.0f} s "
          f"(about ${cost_usd(wall_seconds):.2f} of card time)")
    if attempted:
        print(f"trajectories attempted   {attempted} (rows of {trajectory_table})")
        print(f"seconds per trajectory   {throughput['seconds_per_trajectory']}")
    else:
        print("trajectories attempted   NOT COUNTED: no table with 'traj' in its "
              "path. Row counts of every table, to pick the right one by hand:")
        for name, rows in sorted(tables.items()):
            print(f"    {rows:6d}  {name}")
    print(f"accepted designs         "
          f"{accepted['ranked_complexes']} structures in 3_Ranked/, "
          f"{accepted['ranked_table_rows']} rows in "
          f"{accepted['ranked_table'] or 'no ranked table'}")
    if (accepted["ranked_table_rows"] is not None
            and accepted["ranked_table_rows"] != accepted["ranked_complexes"]):
        print("    THESE DISAGREE. The table and the structures should match; "
              "fewer structures than rows means output we paid for is missing.")
    print(f"acceptance rate          {throughput['acceptance_rate']} "
          f"(accepted / attempted)")
    print(f"candidates per card-hour {throughput['candidates_per_card_hour']}  "
          f"<-- the first of the two rates the pilot exists to measure")
    print(f"complex size             "
          f"{', '.join(str(s) for s in sorted(set(sizes))) or 'NOT MEASURED'} "
          f"residues (target plus binder). Run time and memory scale with this.")

    print("\nCROPPING AND NUMBERING OF THE TARGET")
    if not numbering:
        print("  NOT ANSWERED. No designed complex came back to compare. That is "
              "expected if no trajectory was accepted, and is not a pass.")
    else:
        statuses = sorted({r["status"] for r in numbering})
        for r in numbering[:10]:
            print(f"  {r['file']}: {target_numbering.describe(r)}")
        if len(numbering) > 10:
            print(f"  ... and {len(numbering) - 10} more")
        print(f"  Outcomes across {len(numbering)} complex(es): "
              f"{', '.join(statuses)}")
        if all(r["ok"] for r in numbering):
            print("  The target came back at the input's own numbers, for every "
                  "file checked. analysis/10 can translate it.")
        else:
            print("  The target did NOT come back at the input's numbers for at "
                  "least one file. analysis/10 will refuse those. The files "
                  "checked are the run's own output only, so BindCraft2's shipped "
                  "scaffolds are not the explanation this time.")

    print("\nWHAT IS NOT MEASURED HERE")
    print("  The charge-pair survival rate. That needs "
          "analysis/10_charge_pair_filter.py to read these structures and count "
          "correct pairs against the floor of three, which happens off the card.")
    return report


@app.function(gpu=DEFAULT_GPU, timeout=60 * 45,
              volumes={"/results": results, "/weights": weights})
def validate(config: str = "notag"):
    """One trajectory, to find out whether the campaign file starts up at all.

    `docs/decisions-log.md` Next actions item 2: neither campaign file has ever
    been run, `termini_accessible` and the His-tag off-target have never been
    exercised, and this is named there as a prerequisite before the EGFR campaign
    spends real trajectories.

    What a pass means: BindCraft2 accepted the campaign file, found the trimmed
    fragment and the tag, and carried one trajectory far enough to stop by itself.
    What it does not mean: that the designs are any good, or that one trajectory's
    timing predicts a long run's, since the first trajectory carries the cost of
    compiling the folding models and is therefore the slowest. The smoke run saw
    that directly -- its two trajectories that included compilation took 207 and
    634 seconds against 185 to 574 for the rest.

    `number_of_final_designs` is forced to 1 as well. Left at the committed 200,
    the run would stop on trajectories rather than on designs, which is the same
    outcome here but only by accident.
    """
    return _campaign(config, trajectories=1,
                     run_name=f"egfr-validate-{config}", final_designs=1)


@app.function(gpu=DEFAULT_GPU, timeout=timeout_for_spend(8.0),
              volumes={"/results": results, "/weights": weights})
def pilot(trajectories: int = 150, config: str = "notag",
          spend_ceiling_usd: float = 8.0):
    """The pilot: enough trajectories to measure both rates, and no more.

    150 to 200 trajectories is the band `design/configs/egfr-domain3-h370.md`
    names as what the main run should be sized from. The default is the bottom of
    it, because the figure that matters is a rate and a rate does not improve much
    between 150 and 200 attempts while the cost rises by a third.

    `spend_ceiling_usd` is recorded in the report and printed, but the timeout is
    fixed when this function is decorated, so passing a different value here does
    not move the timeout. Changing the ceiling means editing the decorator. That
    is stated rather than hidden, because a parameter that looks like a limit and
    is not one is worse than no parameter: the earlier smoke script documented a
    `--gpu l4` flag that was never implemented.

    `number_of_final_designs` is left as the committed 200 on purpose. The pilot
    should stop because it ran out of trajectories, not because it collected
    enough designs, or the acceptance rate would be measured against a run that
    stopped early for the other reason.
    """
    print(f"spend ceiling stated as ${spend_ceiling_usd:.2f}; the timeout this "
          f"function was deployed with is {timeout_for_spend(8.0)} s, which is "
          f"${8.0:.2f} of L4 card time")
    return _campaign(config, trajectories=trajectories,
                     run_name=f"egfr-pilot-{config}")


@app.function(gpu=DEFAULT_GPU, timeout=60 * 45,
              volumes={"/results": results, "/weights": weights})
def probe_fragment(config: str = "ext499", trajectories: int = 1):
    """Test A: does a different target fragment clear the confidence gate on its own?

    The validation run was stopped at its first stage on `pLDDT.EGFR_domain3 = 0.35`
    against BindCraft2's default floor of 0.6. That number is the predictor's
    confidence in our target, and it barely moves with the binder -- across the
    PD-L1 smoke run's ten accepted designs it varied by 0.01 -- so the same gate
    would be expected to stop essentially every attempt.

    This runs the same campaign against the fragment extended to 310-499, which
    gives back the disulfide partner the 310-480 cut severed, and **changes no
    filter**. That is the point: if the extended fragment clears 0.6 by itself, the
    fragment was the problem and no rejection rule has to be touched, which is a
    better outcome than relaxing a threshold.

    One trajectory is enough, because the figure being read is a property of the
    target. Whether the attempt is ultimately accepted does not matter here and
    should not be read as the result.
    """
    return _campaign(config, trajectories=trajectories,
                     run_name=f"egfr-probe-{config}", final_designs=1)


@app.function(gpu=DEFAULT_GPU, timeout=timeout_for_spend(2.5),
              volumes={"/results": results, "/weights": weights})
def diagnose(trajectories: int = 10, config: str = "notag",
             min_target_plddt: float = 0.25):
    """Test B: with the confidence gate lowered, do designs clear everything else?

    Keeps the committed 310-480 fragment exactly as decided and lowers only the
    target-confidence floor, through BindCraft2's own `min_target_plddt_final`
    setting. The question is what happens at the gates after it: interface
    confidence, interface error, clashes, and the binder's own fold.

    **What lowering this threshold gives up, stated plainly.** The filter exists to
    catch a target the predictor cannot place confidently inside the complex, which
    can mean the binder is distorting it. Lowering the floor accepts designs where
    that has not been ruled out. The argument for doing it anyway is that our target
    is supplied as a structure whose coordinates are held, that an independent
    predictor rates the same fragment 80.6, and that the default being overridden is
    named for disordered sequence targets rather than structured ones. That is an
    argument, not a measurement, and the number it licenses is a diagnostic setting
    and not a campaign setting.

    0.25 is chosen to sit below the 0.35 the fragment actually scored, so the gate
    stops rejecting everything while still catching an attempt that collapses well
    below the observed baseline. It is a floor relative to a measured value rather
    than a round number, and it is not a proposal for the main run.

    Ten trajectories, which is enough to see whether designs flow through at all and
    to get a first rate, and small enough that being wrong costs about a dollar.
    """
    return _campaign(config, trajectories=trajectories,
                     run_name=f"egfr-diagnose-{config}",
                     extra_settings={"min_target_plddt_final": min_target_plddt})


@app.function(gpu=DEFAULT_GPU, timeout=timeout_for_spend(2.5),
              volumes={"/results": results, "/weights": weights})
def batch_tight4(trajectories: int = 8, config: str = "tight4"):
    """Arm 1: fewer hotspots, everything else held.

    Eleven trajectories against the eight-anchor patch produced no accepted design,
    and `i_pTM` never exceeded 0.62 against a final gate of 0.7. The hypothesis with
    a mechanism behind it is that eight hotspots spread over 24.3 angstroms, with 58
    degrees of angular spread, ask for more interface than a 30-100 residue binder
    can make well, producing a diffuse contact set that satisfies none of them.
    BindCraft2's own PD-L1 example, which accepted 33%, names four hotspots.

    So this cuts the list to four — E344, H358, D368, H370, config numbers A320,
    A334, A344, A346 — and changes **nothing else**. The committed 310-480 fragment,
    the His-tag off-target, `termini_accessible`, the binder length band and every
    filter stay as they are.

    **Why the fragment is deliberately not changed at the same time.** An earlier
    version of this arm also moved to the extended 310-499 fragment, on the grounds
    that it scored best on complex confidence. That was wrong twice over: complex
    confidence is not the gate that is failing, since most trajectories cleared the
    0.6 `pLDDT` floor and then died on interface confidence instead; and moving two
    variables at once would have made any improvement unattributable. Eight
    trajectories buy one signal, so the signal is spent on the one hypothesis with a
    mechanism.

    **`forced_targeting` was considered and cannot be used.** It would address the
    same diagnosis by concentrating contact on the declared hotspots, at no extra
    cost. BindCraft2 refuses it here: the His-tag off-target is a FASTA sequence,
    which triggers its `disordered target` feature, and that feature declares a
    conflict with forced targeting because "forced targeting rebuilds around a
    resolved backbone, which a sequence target does not carry"
    (`bindcraft/settings.py`). Conflicting combinations are refused before a
    campaign starts. Using it would mean dropping the His-tag detargeting, which is
    the one part of the design rule currently working — 32 readings from 0.01 to
    0.08 against a 0.4 ceiling. Not worth trading.
    """
    return _campaign(config, trajectories=trajectories,
                     run_name=f"egfr-tight4-{config}")


@app.function(gpu=DEFAULT_GPU, timeout=timeout_for_spend(2.5),
              volumes={"/results": results, "/weights": weights})
def batch_floors(trajectories: int = 8, config: str = "notag"):
    """Arm 2: the eight-anchor patch, with the interface floors lowered.

    Keeps the committed fragment and all eight hotspots, and lowers only the
    interface-confidence floors (`LOWERED_IPTM_FLOORS`, documented there with what
    it gives up and on whose authority).

    **This is the arm that is meant to produce something.** Nothing in this project
    has ever written a designed complex against EGFR, so
    `analysis/10_charge_pair_filter.py` has never scored a real candidate and the
    charge-pair survival rate has been missing since step 10 was written. Observed
    `i_pTM` runs 0.20 to 0.62, so floors at 0.35 and a final gate at 0.5 should let
    most trajectories through to a structure.

    That rate may matter more than the designs this arm produces, because it answers
    a question about the premise rather than about one campaign: whether the
    charge-pair mechanism this whole project is built on is reachable on this
    epitope at all. A design that reaches a structure and then carries zero correct
    charge pairs would be a far more serious finding than a design BindCraft2
    declined to finish.
    """
    return _campaign(config, trajectories=trajectories,
                     run_name=f"egfr-floors-{config}",
                     extra_settings=LOWERED_IPTM_FLOORS)


@app.function(gpu=DEFAULT_GPU, timeout=timeout_for_spend(3.5),
              volumes={"/results": results, "/weights": weights})
def batch_r1_short(trajectories: int = 10, config: str = "r1-short"):
    """Round one, short band: the tight four hotspots, binder length narrowed to 30-60 aa.

    `egfr-tight4-tight4` (`batch_tight4`) reached BindCraft2's final stage on two of
    its six trajectories, at i_pTM 0.51 and 0.61 -- our highest interface-confidence
    readings against this target to date -- while running the committed 30-100 aa
    binder length band whole. That single wide band cannot say how much of the cost
    of a trajectory comes from the binder's length, because any one trajectory could
    have drawn a length from anywhere in it. Round one splits the band in two to
    answer that: this function runs the bottom half, 30-60 aa; `batch_r1_long` runs
    the top half, 60-100 aa.

    Everything else is held exactly as `batch_tight4` set it: the tight four hotspots
    (E344, H358, D368, H370; config numbers A320, A334, A344, A346), the committed
    310-480 fragment, the His-tag off-target, and `termini_accessible`. The one
    deliberate difference between this function and `batch_r1_long` is the length
    band named in each one's name. Binder length is also the one variable on that
    list that changes compute cost per trajectory -- a longer chain is a larger
    folding problem for BindCraft2's structure predictor to solve at every stage --
    so it is the variable round two's size will be read off of, once both bands have
    a measured card-time-per-trajectory figure to compare.

    `extra_settings` is `ROUND1_SETTINGS`, not `LOWERED_IPTM_FLOORS` alone, because
    of what the tight4 run above did not leave behind: its two trajectories that
    reached the final stage wrote a `<design>_losses.csv` table of scores and
    nothing else, no sequence and no structure, because `save_design_sequences`
    defaults to false and nothing in that run turned it on. `ROUND1_SETTINGS` is the
    same lowered floors with that setting turned on, so a trajectory that again
    clears every gate this time leaves something `analysis/10` or `analysis/18` can
    actually read.

    Ten trajectories by default, inside this function's $3.50 card-time ceiling --
    sized to get a per-trajectory cost reading for the short band without spending
    the budget meant to also cover the long band in `batch_r1_long`.
    """
    return _campaign(config, trajectories=trajectories,
                     run_name=f"egfr-r1-short-{config}", extra_settings=ROUND1_SETTINGS)


@app.function(gpu=DEFAULT_GPU, timeout=timeout_for_spend(3.5),
              volumes={"/results": results, "/weights": weights})
def batch_r1_long(trajectories: int = 10, config: str = "r1-long"):
    """Round one, long band: the tight four hotspots, binder length narrowed to 60-100 aa.

    The counterpart to `batch_r1_short`, which that function's docstring explains in
    full: both are round one of splitting the committed 30-100 aa binder length band
    that `egfr-tight4-tight4` (`batch_tight4`) ran whole, in order to measure how
    much of a trajectory's cost tracks the binder's length rather than the hotspot
    set or the fragment, neither of which changes between the two functions.

    This one runs the top half of the band, 60-100 aa, against the same tight four
    hotspots (E344, H358, D368, H370; config numbers A320, A334, A344, A346), the
    same committed 310-480 fragment, and the same His-tag off-target as
    `batch_r1_short`. The length band is the only deliberate difference between the
    two functions.

    `extra_settings` is `ROUND1_SETTINGS` for the same reason it is in
    `batch_r1_short`: `egfr-tight4-tight4` reached the final stage on two of its six
    trajectories and still left no recoverable sequence or structure, because
    `save_design_sequences` defaults to false and was never turned on there.
    `ROUND1_SETTINGS` carries the same lowered interface floors with that setting
    turned on.

    Ten trajectories by default, inside this function's own $3.50 card-time
    ceiling, run alongside `batch_r1_short` rather than instead of it -- the
    comparison that sizes round two needs a cost-per-trajectory figure for both
    bands, not just one.
    """
    return _campaign(config, trajectories=trajectories,
                     run_name=f"egfr-r1-long-{config}", extra_settings=ROUND1_SETTINGS)


@app.function(gpu=DEFAULT_GPU, timeout=timeout_for_spend(7.75),
              volumes={"/results": results, "/weights": weights})
def batch_r2_short(trajectories: int = 40, config: str = "r1-short",
                   workers: int = 2, length_band: str = "33,60"):
    """Round two, short band: round one's configuration unchanged, on two workers.

    Everything about the design is held exactly as `batch_r1_short` ran it -- the
    tight four hotspots (E344, H358, D368, H370; config numbers A320, A334, A344,
    A346), the committed 310-480 fragment, the 30-60 aa binder length band, the
    His-tag off-target, `termini_accessible`, and `ROUND1_SETTINGS` with its lowered
    interface floors and `save_design_sequences` turned on. Nothing about what is
    being designed changes between round one and round two. Two settings change:
    how many attempts are paid for, and how many of them run at once.

    WHY TWO WORKERS, AND WHY IT TAKES A PATCH TO GET THEM
    -----------------------------------------------------
    Round one ran one worker on the L4 and peaked at 8,894 MiB of the card's 23,034,
    which is 39% of it. Two workers would therefore use about 17.4 GB and leave
    roughly 5 GB spare. The decisions log recorded for two campaigns that the second
    worker was out of reach because BindCraft2 holds fan-out to the trajectory
    budget. That reason was wrong, and the real one is in `design_workers_per_gpu`
    (`bindcraft/design_workers.py`, lines 61-70): an unconditional clamp by a memory
    estimate that budgets 11.78 GB for a worker on a complex of our size and so
    allows exactly one. Asking for two with `workers_per_gpu` does not get two; the
    clamp overrides the request. `WORKER_MEMORY_PATCH` recalibrates the estimate at
    run time against round one's measured peak, from outside the package. The
    derivation, with source line numbers, is at the end of Pipeline status in
    `docs/decisions-log.md`.

    The second worker is a throughput lever and not a design change. It does not
    alter what is designed, only how many attempts fit in an hour of card time, so
    nothing it produces is comparable to round one on any axis except cost and rate.

    WHY THE BAND IS 33-60 AND NOT 30-60
    -----------------------------------
    The first attempt at this run used round one's 30-60 band and exposed a
    partition this project had never seen, because it had never run two workers.
    BindCraft2 groups binder lengths into padding buckets of 32
    (`padded_prediction_length`): 30, 31 and 32 pad to 32, and 33 to 60 pad to 64.
    That is two buckets. `assign_worker_length_buckets`
    (`bindcraft/design_workers.py`, lines 111-118) then pins each worker to one
    bucket when there are at least as many workers as buckets, so worker 1 was
    locked to lengths 30-32 and worker 0 to 33-60.

    The two workers draw from one shared trajectory counter rather than from a
    per-bucket share -- observed directly, worker 0 on trajectory 2 while worker 1
    was on trajectory 1 -- so the split is about half each by count. Half the
    budget would therefore have gone to a three-length window at the bottom of the
    band. No trajectory at any length below 33 has ever completed in this project,
    and round one's single completion was at length 48.

    Dropping to 33-60 leaves one bucket, and line 112 returns the plan unchanged
    when there is only one, so no worker gets a length restriction and both draw
    from the whole band. Authorised by the owner on 5 October 2026 as a
    three-length trim, recorded in `docs/decisions-log.md` so it does not read as
    configuration drift. It is set with `--set binder_lengths`, leaving the
    committed campaign file as round one ran it.

    SIZING
    ------
    Forty trajectories, against a $8.00 card-time ceiling expressed as this
    function's timeout -- 10 hours at the L4's $0.80 an hour. Round one's billed
    cost was $0.316 a trajectory at one worker, so 40 attempts at one worker would
    be $12.64, beyond both the ceiling and the $10.87 of credit left. At two workers
    with no contention that halves to about $6.32; with contention at 1.6 times
    throughput rather than 2, about $7.90. The ceiling is what bounds the bad case,
    so the run may stop on its timeout with fewer than 40 attempts rather than on
    its trajectory budget. That is the intended behaviour and not a failure.

    WHAT STOPS IT
    -------------
    Two conditions, both decided before launch, both reported by this function:
    the run dies with an out-of-memory error, or peak card memory exceeds 20,000
    MiB. The caller also compares the first completed trajectory's wall-clock
    against round one's 1,040 s mean design time; two workers that are thrashing
    rather than sharing show up there before most of the money is spent.
    """
    low, high = (int(part) for part in length_band.split(","))
    settings = dict(ROUND1_SETTINGS, binder_lengths=f"[{low},{high}]")
    return _campaign(config, trajectories=trajectories,
                     run_name=f"egfr-r2-short-{config}",
                     extra_settings=settings, workers=workers)


@app.function(timeout=60 * 10, volumes={"/results": results})
def inspect_volume(path: str = "/results"):
    """List what is on the results volume. No card, so it costs almost nothing.

    For checking on a detached run without paying for a GPU container, and for
    finding out what came back before copying it down.
    """
    _run(f"find {path} -maxdepth 4 | head -200", cwd="/", venv=False)
    _run(f"du -sh {path}/* 2>/dev/null || true", cwd="/", venv=False)


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _self_test():
    """Check the arithmetic and the layout helpers, and that each fails when broken.

    Run with `./.venv/bin/python design/modal/egfr_campaign.py`. No Modal account,
    no card, no network. CLAUDE.md's rule is that a check never seen to fail is
    not known to work, so each case below is paired with a deliberately broken
    input that must come out differently.
    """
    import json
    import pathlib
    import tempfile

    failures = []

    def check(label, got, want):
        if got != want:
            failures.append(f"{label}: got {got!r}, wanted {want!r}")

    # Cost arithmetic. $0.80 an hour means an hour costs $0.80 and a dollar buys
    # 4,500 seconds.
    check("an hour on an L4 costs $0.80", round(cost_usd(3600), 4), 0.8)
    check("$8 buys 36,000 s", timeout_for_spend(8.0), 36000)
    check("$1 buys 4,500 s", timeout_for_spend(1.0), 4500)

    # The rates. The PD-L1 smoke run's real figures: 5540 s, 30 trajectories, 10
    # accepted. 10 accepted in 5540 s is 10 / 1.5389 hours = 6.5 an hour.
    smoke = rates(5540, 30, 10)
    check("smoke run card hours", smoke["card_hours"], 1.539)
    check("smoke run s/trajectory", smoke["seconds_per_trajectory"], 184.7)
    check("smoke run acceptance", smoke["acceptance_rate"], 0.3333)
    check("smoke run candidates/card-hour", smoke["candidates_per_card_hour"], 6.5)

    # A run that accepted nothing must not report a rate of zero as though it had
    # measured one, and a run with no trajectory table must not divide by it.
    check("nothing accepted gives no rate",
          rates(5540, 30, 0)["candidates_per_card_hour"], None)
    check("no trajectory table gives no rate",
          rates(5540, None, 10)["seconds_per_trajectory"], None)

    # Separating the structure target from the sequence off-target. This is the
    # case the smoke script's version got wrong.
    notag = json.loads(
        pathlib.Path(_local("design", "configs",
                            "egfr-domain3-h370-notag.json")).read_text())
    structure, off = structure_target(notag)
    check("the structure target is the EGFR fragment",
          structure["name"], "EGFR_domain3")
    check("the His tag is an off-target and not the structure",
          [e["name"] for e in off], ["HisTag"])
    check("the off-target keeps its negative weight", off[0]["weight"], -0.5)

    plain = json.loads(
        pathlib.Path(_local("design", "configs",
                            "egfr-domain3-h370.json")).read_text())
    structure, off = structure_target(plain)
    check("the plain config has one target and no off-target",
          (structure["name"], off), ("EGFR_domain3", []))

    # The diagnostic config must differ from the committed one in the target file
    # and in nothing else, or a comparison between them measures more than the
    # boundary. Asserted rather than trusted, because it is hand-written.
    ext = json.loads(
        pathlib.Path(_local("design", "configs", "diagnostic",
                            "ext499-notag.json")).read_text())
    differing = {k for k in set(notag) | set(ext)
                 if notag.get(k) != ext.get(k)}
    check("the diagnostic config differs only in its name, folder and targets",
          differing, {"campaign_name", "project_folder", "targets"})
    ext_structure, ext_off = structure_target(ext)
    check("the diagnostic points at the extended fragment",
          ext_structure["target_path"].endswith("6aru_domain3_ext499.pdb"), True)
    check("the diagnostic keeps the same hotspots, so the aim is unchanged",
          ext_structure["hotspots"], structure_target(notag)[0]["hotspots"])
    check("the diagnostic keeps the His tag off-target at the same weight",
          [(e["name"], e["weight"]) for e in ext_off], [("HisTag", -0.5)])

    # Arm 1 must differ from the committed campaign in the hotspot list and in
    # nothing else, because that is the whole claim the arm rests on. Asserted
    # because an earlier draft of it also changed the target fragment, which would
    # have made any improvement unattributable.
    tight = json.loads(
        pathlib.Path(_local("design", "configs", "diagnostic",
                            "tight4-notag.json")).read_text())
    differing = {k for k in set(notag) | set(tight) if notag.get(k) != tight.get(k)}
    check("the tight-hotspot config differs only in name, folder and targets",
          differing, {"campaign_name", "project_folder", "targets"})
    tight_structure, tight_off = structure_target(tight)
    check("arm 1 keeps the committed fragment, not the extended one",
          tight_structure["target_path"].endswith("6aru_domain3.pdb"), True)
    check("arm 1 narrows the hotspots to the tight four",
          tight_structure["hotspots"], "A320,A334,A344,A346")
    check("arm 1 keeps the His tag off-target at the same weight",
          [(e["name"], e["weight"]) for e in tight_off], [("HisTag", -0.5)])
    # The four must be a subset of the committed eight, or this is a different
    # epitope rather than a narrower aim at the same one.
    eight = set(structure_target(notag)[0]["hotspots"].split(","))
    four = set(tight_structure["hotspots"].split(","))
    check("the tight four are a subset of the committed eight",
          four <= eight, True)
    check("the tight four keep both target histidines, H358 and H370",
          {"A334", "A346"} <= four, True)

    # The lowered floors must only ever lower. A typo raising one would tighten the
    # gate the arm exists to open, and would look like a failed hypothesis.
    defaults = {"min_iptm_anneal": 0.5, "min_iptm_harden": 0.5,
                "min_iptm_mutate": 0.5, "min_iptm_final": 0.7}
    check("every lowered floor is below BindCraft2's default",
          all(LOWERED_IPTM_FLOORS[k] < v for k, v in defaults.items()), True)
    check("the lowered floors touch only the interface metric",
          all("iptm" in k for k in LOWERED_IPTM_FLOORS), True)

    # Broken on purpose: a campaign whose only target is a FASTA has no structure
    # to compare numbering against, and must say so rather than returning the
    # FASTA as though it were one.
    check("a FASTA-only campaign has no structure target",
          structure_target({"targets": [{"name": "HisTag",
                                         "target_path": "x.fasta"}]})[0], None)
    # And if the order were reversed in the file, the structure must still win.
    reversed_spec = {"targets": [{"name": "HisTag", "target_path": "x.fasta"},
                                 {"name": "EGFR", "target_path": "y.pdb"}]}
    check("target order in the file does not decide which is the structure",
          structure_target(reversed_spec)[0]["name"], "EGFR")

    # The mirrored layout. This is what makes the relative target_path resolve,
    # so it is worth asserting rather than trusting.
    config_in_container = pathlib.Path(WORK_DIR) / CONFIGS["notag"]
    resolved = (config_in_container.parent / notag["targets"][0]["target_path"])
    check("the mounted layout resolves target_path",
          str(pathlib.Path(os.path.normpath(resolved))),
          f"{WORK_DIR}/data/structures/6aru_domain3.pdb")
    # Broken on purpose: flattening the files into one directory breaks it, which
    # is the mistake the mirrored layout exists to prevent.
    flat = pathlib.Path("/work") / "egfr-domain3-h370-notag.json"
    flat_resolved = os.path.normpath(
        flat.parent / notag["targets"][0]["target_path"])
    if flat_resolved == f"{WORK_DIR}/data/structures/6aru_domain3.pdb":
        failures.append("flattening the mount should break target_path and did not")

    # Every local file the image claims to add must actually be there, or the run
    # fails after a container has started rather than before one is asked for.
    for relative in list(CONFIGS.values()) + [
            "data/structures/6aru_domain3.pdb",
            "data/structures/6aru_domain3_ext499.pdb",
            "data/sequences/his-tag-offtarget.fasta",
            "analysis/target_numbering.py"]:
        if not os.path.exists(_local(*relative.split("/"))):
            failures.append(f"missing locally, so it cannot be mounted: {relative}")

    # The streamed and buffered paths through `_run` must agree on what a command
    # said and on whether it succeeded, or the pilot would be watched through a
    # code path the validation run never exercised. Run locally with venv off and
    # cwd at the root, because BindCraft2's install is not on this machine.
    rc_buffered, out_buffered = _run(
        "printf 'alpha\\nbeta\\n'", cwd="/", venv=False, want_output=True)
    rc_streamed, out_streamed = _run(
        "printf 'alpha\\nbeta\\n'", cwd="/", venv=False, want_output=True,
        stream=True)
    check("streaming and buffering agree on the output",
          out_streamed, out_buffered)
    check("streaming and buffering agree on success",
          (rc_streamed, rc_buffered), (0, 0))

    # A failing command must report its exit code through the streamed path too.
    # Without this, a campaign that died would be read as one that finished.
    rc_fail, out_fail = _run("printf 'partial\\n'; exit 3", cwd="/", venv=False,
                             want_output=True, stream=True)
    check("a failed command's exit code survives streaming", rc_fail, 3)
    check("output written before a failure is still returned",
          out_fail, "partial\n")

    # Streaming merges standard error into the stream, which is the property the
    # progress watch depends on: BindCraft2's progress lines and its warnings have
    # to arrive in one ordered stream to be readable together.
    _, out_merged = _run("printf 'to-stdout\\n'; printf 'to-stderr\\n' >&2",
                         cwd="/", venv=False, want_output=True, stream=True)
    check("streaming merges standard error into the stream",
          sorted(out_merged.split()), ["to-stderr", "to-stdout"])
    # Broken on purpose: the buffered path returns standard output alone, so a
    # watch built on it would miss anything written to standard error. That
    # asymmetry is deliberate and is asserted so it stays deliberate.
    _, out_split = _run("printf 'to-stdout\\n'; printf 'to-stderr\\n' >&2",
                        cwd="/", venv=False, want_output=True)
    if "to-stderr" in out_split:
        failures.append("the buffered path returned standard error, so the "
                        "documented difference between the two paths is wrong")

    # Counting accepted designs, and the scaffold-file trap. A tree holding both
    # real output and BindCraft2's shipped template scaffolds must count only the
    # real output, and must not offer a scaffold to the numbering check.
    with tempfile.TemporaryDirectory() as temp:
        root = pathlib.Path(temp)
        for relative in ("3_Ranked/design_1_seq1_EGFR.cif",
                         "3_Ranked/design_1_seq1_monomer.cif",
                         "3_Ranked/design_2_seq1_EGFR.cif",
                         "1_Trajectories/traj_1.cif",
                         "2_Refolded/refold_1.cif",
                         "scaffolds/VHH.cif",
                         "scaffolds/Fab.cif"):
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("")
        (root / "3_Ranked").mkdir(exist_ok=True)
        (root / "3_Ranked" / "!_Ranked.csv").write_text(
            "name,i_pDAE\nd1,0.22\nd2,0.19\n")
        (root / "1_Trajectories" / "trajectory_stats.csv").write_text(
            "name\na\nb\nc\n")

        counts = count_accepted(root)
        check("ranked structures counted, monomers excluded",
              counts["ranked_complexes"], 2)
        check("ranked table rows counted", counts["ranked_table_rows"], 2)
        check("trajectories counted from the trajectory table",
              _count_trajectories(root)[0], 3)
        # Sorted by full path, so the files group by stage -- trajectories,
        # then refolded, then ranked -- which is the order they were produced in
        # and the order a reader of the report wants them in.
        names = [str(p.relative_to(root)) for p in numbering_candidates(root)]
        check("the numbering check sees only the run's own output, by stage",
              names, ["1_Trajectories/traj_1.cif", "2_Refolded/refold_1.cif",
                      "3_Ranked/design_1_seq1_EGFR.cif",
                      "3_Ranked/design_2_seq1_EGFR.cif"])
        if any(n in names for n in ("VHH.cif", "Fab.cif")):
            failures.append("a shipped scaffold reached the numbering check, "
                            "which is the 2 October bug")

    print("\n".join(f"FAIL  {f}" for f in failures) if failures
          else "self-test PASSED: arithmetic, target separation, mounted layout, "
               "accepted-design counting and the scaffold restriction all hold, "
               "and each broken-on-purpose case failed as it should.")
    return not failures


@app.local_entrypoint()
def main():
    print(__doc__.split("HOW TO RUN")[-1])


if __name__ == "__main__":
    raise SystemExit(0 if _self_test() else 1)
