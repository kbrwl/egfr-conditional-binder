#!/usr/bin/env bash
#
# collect_campaigns.sh -- bring finished campaigns home and score them, in one command.
#
# WHY THIS EXISTS. A campaign runs detached on rented hardware and writes its output
# to a Modal volume. Getting from there to an answer takes several steps in a fixed
# order, and the PD-L1 smoke run on 2 October 2026 was left unattended and its final
# report never collected, because the steps were in somebody's head rather than in a
# file. This is that sequence, written down.
#
# It is safe to run at any time, including while a campaign is still going: the
# snapshots mean the volume always holds whatever has finished so far, so an early
# run gives a partial answer rather than an error.
#
# WHAT IT DOES, per campaign named on the command line or found on the volume:
#   1. copies the campaign folder down into results/candidates/
#   2. if any designed complex came back, scores it with analysis/10 (charge pairs)
#      and preserves that campaign's copy of the three output files, because
#      analysis/10 writes to fixed paths and a second campaign would overwrite them
#   3. builds the per-attempt ledger with analysis/18, which records every attempt
#      including the ones that failed, with the stage and the reason
#   4. prints the headline numbers
#
# USAGE
#   bash design/modal/collect_campaigns.sh                     # everything on the volume
#   bash design/modal/collect_campaigns.sh egfr-tight4-tight4   # one campaign
#
set -uo pipefail

cd "$(dirname "$0")/../.." || exit 1
PY=./.venv/bin/python
MODAL=./.venv/bin/modal
VOLUME=bindcraft2-results
OUT=results/candidates

if [ "$#" -gt 0 ]; then
    CAMPAIGNS=("$@")
else
    # Everything on the volume. `modal volume ls` prints one name a line at the top
    # level, which is what the campaign folders are.
    mapfile -t CAMPAIGNS < <($MODAL volume ls "$VOLUME" 2>/dev/null | tr -d ' ')
fi

if [ "${#CAMPAIGNS[@]}" -eq 0 ]; then
    echo "No campaigns found on the volume $VOLUME. Nothing to collect."
    exit 1
fi

echo "============================================================"
echo "Collecting ${#CAMPAIGNS[@]} campaign(s): ${CAMPAIGNS[*]}"
echo "============================================================"

mkdir -p "$OUT"

for name in "${CAMPAIGNS[@]}"; do
    [ -z "$name" ] && continue
    echo ""
    echo "------------------------------------------------------------"
    echo "$name"
    echo "------------------------------------------------------------"

    $MODAL volume get "$VOLUME" "$name" "$OUT/" --force >/dev/null 2>&1 \
        || { echo "  could not copy $name down; skipping"; continue; }
    folder="$OUT/$name"
    [ -d "$folder" ] || { echo "  $folder is not there after copying; skipping"; continue; }

    # Designed complexes, excluding the *_monomer.cif files that hold the binder
    # with no target and so carry no interface to measure.
    complexes=$(find "$folder" -name '*.cif' ! -name '*_monomer.cif' \
                     -path '*3_Ranked*' 2>/dev/null | wc -l | tr -d ' ')
    echo "  designed complexes in 3_Ranked/: $complexes"

    if [ "$complexes" -gt 0 ]; then
        echo "  scoring charge pairs with analysis/10 ..."
        $PY analysis/10_charge_pair_filter.py --candidates "$folder" 2>&1 | tail -20
        # analysis/10 writes to fixed paths, so each campaign's copy is preserved
        # before the next campaign overwrites them.
        for f in data/derived/10-candidate-pairs.csv \
                 data/derived/10-candidate-summary.csv \
                 results/candidates/shortlist.csv; do
            [ -f "$f" ] && cp "$f" "${f%.*}--$name.${f##*.}"
        done
        echo "  preserved this campaign's copies as *--$name.csv"
    else
        echo "  no complexes, so there is nothing for analysis/10 to score."
        echo "  That is a real outcome and not an error: a trajectory stopped"
        echo "  before the stage that writes a structure leaves none behind."
    fi

    echo "  building the per-attempt ledger with analysis/18 ..."
    $PY analysis/18_campaign_inventory.py --campaign "$folder" 2>&1 \
        | grep -E 'attempts|accepted|rejected|scored|Files written|ledger|inventory|could not' \
        | sed 's/^/    /'

    if [ -f "$folder/campaign-report.json" ]; then
        echo "  headline numbers from the run's own report:"
        $PY -c "
import json, sys
r = json.load(open('$folder/campaign-report.json'))
for k in ('card', 'peak_memory_used_mib', 'wall_clock_seconds',
          'estimated_cost_usd', 'trajectories_attempted',
          'seconds_per_trajectory', 'acceptance_rate',
          'candidates_per_card_hour', 'complexes_checked'):
    print(f'    {k:26s} {r.get(k)}')
print(f\"    {'extra_settings':26s} {r.get('extra_settings')}\")
" 2>/dev/null || echo "    could not read it"
    else
        echo "  no campaign-report.json yet, so the run had not finished when the"
        echo "  last snapshot was taken. The trajectory table below is still valid."
    fi

    table="$folder/1_Trajectories/!_Trajectories.csv"
    if [ -f "$table" ]; then
        echo "  every attempt, with where it stopped and what it cost:"
        $PY -c "
import csv
rows = list(csv.DictReader(open('$table')))
total = 0.0
for r in rows:
    t = dict(p.split('=', 1) for p in r.get('Timing', '').split(';') if '=' in p)
    d = float(t.get('design', 0) or 0); total += d
    stopped = r.get('terminated') or 'PASSED'
    print(f\"    {r.get('trajectory', '?'):>3}  len={r.get('length', '?'):>3}  \"
          f\"{stopped:<8}  {d:7.1f}s\")
if rows:
    print(f'    {len(rows)} attempts, {total:.0f} s of design, '
          f'mean {total/len(rows):.0f} s, about \${total*0.80/3600:.2f}')
" 2>/dev/null || echo "    could not read the trajectory table"
    fi
done

echo ""
echo "============================================================"
echo "Done. The ledger is data/derived/18-trajectory-ledger.csv and the"
echo "findings are results/findings/18-campaign-inventory.md."
echo "Both are tracked, so the record survives even though the bulk"
echo "structures under results/candidates/ are not committed."
echo "============================================================"
