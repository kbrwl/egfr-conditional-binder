# Finding 18 — the complete record of the campaign

Computed by `analysis/18_campaign_inventory.py`. Do not hand-edit; rerun the script.

Campaign folder: `results/candidates/egfr-tight4-tight4`

## 1. Coverage

Every attempt the campaign made, and what is known about each. A record holding only the survivors would make the pipeline look better than it is and would leave the acceptance rate with no denominator, so the counts below are the point of this file.

| | count |
|---|---|
| attempts | 6 |
| accepted by bindcraft | 0 |
| rejected by bindcraft | 6 |
| scored by us | 0 |
| not scored by us | 6 |
| accepted but unscored | 0 |

## 2. Files retained

Retained means the file exists and was copied back. It is a separate question from whether it is committed to version control: bulk design output deliberately is not, because it is large, and `.gitignore` admits only the final shortlist. The fingerprints are in `data/derived/18-file-inventory.csv`, which is tracked, so the record survives even where the structures themselves are not committed.

| stage | files | bytes |
|---|---|---|
| 1_Trajectories | 7 | 134,391 |
| other | 45 | 54,998,527 |

## 3. Tables the run wrote

- `1_Trajectories/!_Trajectories.csv` — 6 rows, 14 columns
- `1_Trajectories/egfr-domain3-h370-tight4_detarget_l31_37c084aa3b34c399/egfr-domain3-h370-tight4_detarget_l31_37c084aa3b34c399_losses.csv` — 240 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-tight4_detarget_l41_8d70e8e0582565b7/egfr-domain3-h370-tight4_detarget_l41_8d70e8e0582565b7_losses.csv` — 280 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-tight4_detarget_l60_f130f4efac9a6117/egfr-domain3-h370-tight4_detarget_l60_f130f4efac9a6117_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-tight4_detarget_l68_341198730c9639d7/egfr-domain3-h370-tight4_detarget_l68_341198730c9639d7_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-tight4_detarget_l79_ca6b16aa615baea5/egfr-domain3-h370-tight4_detarget_l79_ca6b16aa615baea5_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-tight4_detarget_l88_520cd9438896db72/egfr-domain3-h370-tight4_detarget_l88_520cd9438896db72_losses.csv` — 280 rows, 30 columns

## 4. What could not be determined, and why

- Nothing. Every table joined to an attempt and every attempt carried a decision.

## 5. What this does not establish

- Nothing here says a design works. It says what was tried, what survived each gate and why, which is the record needed to check the rates afterwards.
- The charge-pair verdicts are step 10's and carry step 10's limits: a predicted structure scored against a predicted interface, with the pH switch inferred from geometry and never measured.
