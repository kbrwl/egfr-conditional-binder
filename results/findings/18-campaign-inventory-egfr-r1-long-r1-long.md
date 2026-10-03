# Finding 18 — the complete record of the campaign

Computed by `analysis/18_campaign_inventory.py`. Do not hand-edit; rerun the script.

Campaign folder: `results/candidates/egfr-r1-long-r1-long`

## 1. Coverage

Every attempt the campaign made, and what is known about each. A record holding only the survivors would make the pipeline look better than it is and would leave the acceptance rate with no denominator, so the counts below are the point of this file.

| | count |
|---|---|
| attempts | 10 |
| accepted by bindcraft | 0 |
| rejected by bindcraft | 10 |
| scored by us | 0 |
| not scored by us | 10 |
| accepted but unscored | 0 |

## 2. Files retained

Retained means the file exists and was copied back. It is a separate question from whether it is committed to version control: bulk design output deliberately is not, because it is large, and `.gitignore` admits only the final shortlist. The fingerprints are in `data/derived/18-file-inventory-egfr-r1-long-r1-long.csv`, which is tracked, so the record survives even where the structures themselves are not committed.

| stage | files | bytes |
|---|---|---|
| 1_Trajectories | 21 | 10,541,032 |
| other | 1 | 2,676 |

## 3. Tables the run wrote

- `1_Trajectories/!_Trajectories.csv` — 10 rows, 14 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l100_6a853d0897bc270f/egfr-domain3-h370-r1-long_detarget_l100_6a853d0897bc270f_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l100_cf5bf8f31b8a7748/egfr-domain3-h370-r1-long_detarget_l100_cf5bf8f31b8a7748_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l63_86868b8ece95bfa0/egfr-domain3-h370-r1-long_detarget_l63_86868b8ece95bfa0_losses.csv` — 240 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l67_9d835c085ebf8e22/egfr-domain3-h370-r1-long_detarget_l67_9d835c085ebf8e22_losses.csv` — 150 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l70_47af44f4006bcf71/egfr-domain3-h370-r1-long_detarget_l70_47af44f4006bcf71_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l74_6128581b478859ae/egfr-domain3-h370-r1-long_detarget_l74_6128581b478859ae_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l76_c05b4747060af5a1/egfr-domain3-h370-r1-long_detarget_l76_c05b4747060af5a1_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l81_5342e17194d92345/egfr-domain3-h370-r1-long_detarget_l81_5342e17194d92345_losses.csv` — 150 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l87_4d2f3c6736bfb423/egfr-domain3-h370-r1-long_detarget_l87_4d2f3c6736bfb423_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-long_detarget_l94_bc82ebce2c3cf818/egfr-domain3-h370-r1-long_detarget_l94_bc82ebce2c3cf818_losses.csv` — 100 rows, 30 columns

## 4. What could not be determined, and why

- Nothing. Every table joined to an attempt and every attempt carried a decision.

## 5. What this does not establish

- Nothing here says a design works. It says what was tried, what survived each gate and why, which is the record needed to check the rates afterwards.
- The charge-pair verdicts are step 10's and carry step 10's limits: a predicted structure scored against a predicted interface, with the pH switch inferred from geometry and never measured.
