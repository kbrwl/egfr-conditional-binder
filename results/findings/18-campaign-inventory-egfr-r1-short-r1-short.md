# Finding 18 — the complete record of the campaign

Computed by `analysis/18_campaign_inventory.py`. Do not hand-edit; rerun the script.

Campaign folder: `results/candidates/egfr-r1-short-r1-short`

## 1. Coverage

Every attempt the campaign made, and what is known about each. A record holding only the survivors would make the pipeline look better than it is and would leave the acceptance rate with no denominator, so the counts below are the point of this file.

| | count |
|---|---|
| attempts | 40 |
| accepted by bindcraft | 0 |
| rejected by bindcraft | 40 |
| scored by us | 10 |
| not scored by us | 30 |
| accepted but unscored | 0 |

## 2. Files retained

Retained means the file exists and was copied back. It is a separate question from whether it is committed to version control: bulk design output deliberately is not, because it is large, and `.gitignore` admits only the final shortlist. The fingerprints are in `data/derived/18-file-inventory-egfr-r1-short-r1-short.csv`, which is tracked, so the record survives even where the structures themselves are not committed.

| stage | files | bytes |
|---|---|---|
| 1_Trajectories | 21 | 7,078,612 |
| 2_Refolded | 31 | 2,527,754 |
| other | 1 | 14,956 |

## 3. Tables the run wrote

- `1_Trajectories/!_Trajectories.csv` — 10 rows, 51 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l32_9e1c275495126373/egfr-domain3-h370-r1-short_detarget_l32_9e1c275495126373_losses.csv` — 150 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l36_f69e1780445bdc40/egfr-domain3-h370-r1-short_detarget_l36_f69e1780445bdc40_losses.csv` — 240 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l38_b07af46868aa55f7/egfr-domain3-h370-r1-short_detarget_l38_b07af46868aa55f7_losses.csv` — 240 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l38_eb38ef4b3f33db61/egfr-domain3-h370-r1-short_detarget_l38_eb38ef4b3f33db61_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l44_caf1da1213a96faa/egfr-domain3-h370-r1-short_detarget_l44_caf1da1213a96faa_losses.csv` — 240 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4/egfr-domain3-h370-r1-short_detarget_l48_b495b644c4981ad4_losses.csv` — 280 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l54_431efedac9a0f2c0/egfr-domain3-h370-r1-short_detarget_l54_431efedac9a0f2c0_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l56_43567bfc1793b68d/egfr-domain3-h370-r1-short_detarget_l56_43567bfc1793b68d_losses.csv` — 280 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l57_0b53e72e622621d8/egfr-domain3-h370-r1-short_detarget_l57_0b53e72e622621d8_losses.csv` — 240 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l60_076e72ac0632960e/egfr-domain3-h370-r1-short_detarget_l60_076e72ac0632960e_losses.csv` — 250 rows, 30 columns
- `2_Refolded/!_Refolded.csv` — 10 rows, 50 columns

## 4. What could not be determined, and why

- Nothing. Every table joined to an attempt and every attempt carried a decision.

## 5. What this does not establish

- Nothing here says a design works. It says what was tried, what survived each gate and why, which is the record needed to check the rates afterwards.
- The charge-pair verdicts are step 10's and carry step 10's limits: a predicted structure scored against a predicted interface, with the pH switch inferred from geometry and never measured.
