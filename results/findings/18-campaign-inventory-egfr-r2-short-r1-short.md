# Finding 18 — the complete record of the campaign

Computed by `analysis/18_campaign_inventory.py`. Do not hand-edit; rerun the script.

Campaign folder: `results/candidates/egfr-r2-short-r1-short`

## 1. Coverage

Every attempt the campaign made, and what is known about each. A record holding only the survivors would make the pipeline look better than it is and would leave the acceptance rate with no denominator, so the counts below are the point of this file.

| | count |
|---|---|
| attempts | 23 |
| accepted by bindcraft | 0 |
| rejected by bindcraft | 23 |
| scored by us | 0 |
| not scored by us | 23 |
| accepted but unscored | 0 |

## 2. Files retained

Retained means the file exists and was copied back. It is a separate question from whether it is committed to version control: bulk design output deliberately is not, because it is large, and `.gitignore` admits only the final shortlist. The fingerprints are in `data/derived/18-file-inventory-egfr-r2-short-r1-short.csv`, which is tracked, so the record survives even where the structures themselves are not committed.

| stage | files | bytes |
|---|---|---|
| 1_Trajectories | 47 | 13,451,774 |
| other | 32 | 24,884,131 |

## 3. Tables the run wrote

- `1_Trajectories/!_Trajectories.csv` — 23 rows, 13 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l33_f5d4ccd08dae55f0/egfr-domain3-h370-r1-short_detarget_l33_f5d4ccd08dae55f0_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l34_01def70b4402b69a/egfr-domain3-h370-r1-short_detarget_l34_01def70b4402b69a_losses.csv` — 150 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l35_598960479a168cf6/egfr-domain3-h370-r1-short_detarget_l35_598960479a168cf6_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l37_d0f5ea2e88e6bf01/egfr-domain3-h370-r1-short_detarget_l37_d0f5ea2e88e6bf01_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l39_422c5de1d503f531/egfr-domain3-h370-r1-short_detarget_l39_422c5de1d503f531_losses.csv` — 280 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l40_33993c77fdc787dd/egfr-domain3-h370-r1-short_detarget_l40_33993c77fdc787dd_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l40_acbe7da37b5a948e/egfr-domain3-h370-r1-short_detarget_l40_acbe7da37b5a948e_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l41_b5edac2b613523c8/egfr-domain3-h370-r1-short_detarget_l41_b5edac2b613523c8_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l42_00624845fd22da85/egfr-domain3-h370-r1-short_detarget_l42_00624845fd22da85_losses.csv` — 150 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l43_82ad44be8ae3d7bc/egfr-domain3-h370-r1-short_detarget_l43_82ad44be8ae3d7bc_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l43_b94f35b52a82abc3/egfr-domain3-h370-r1-short_detarget_l43_b94f35b52a82abc3_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l45_25919fefd2f09b9f/egfr-domain3-h370-r1-short_detarget_l45_25919fefd2f09b9f_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l47_609d2e68c78d107c/egfr-domain3-h370-r1-short_detarget_l47_609d2e68c78d107c_losses.csv` — 280 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l51_0c907544060bd0a7/egfr-domain3-h370-r1-short_detarget_l51_0c907544060bd0a7_losses.csv` — 240 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l52_968c4b6a602256e3/egfr-domain3-h370-r1-short_detarget_l52_968c4b6a602256e3_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l53_43903366ec75868e/egfr-domain3-h370-r1-short_detarget_l53_43903366ec75868e_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l56_40964173b4fd7b26/egfr-domain3-h370-r1-short_detarget_l56_40964173b4fd7b26_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l56_96ebc1997a70648e/egfr-domain3-h370-r1-short_detarget_l56_96ebc1997a70648e_losses.csv` — 280 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l56_b28fbc1e8a391e8e/egfr-domain3-h370-r1-short_detarget_l56_b28fbc1e8a391e8e_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l58_59d74b7ee3dd3a65/egfr-domain3-h370-r1-short_detarget_l58_59d74b7ee3dd3a65_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l59_401074e6866d9f99/egfr-domain3-h370-r1-short_detarget_l59_401074e6866d9f99_losses.csv` — 250 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l59_b9a9d1b79eb863f9/egfr-domain3-h370-r1-short_detarget_l59_b9a9d1b79eb863f9_losses.csv` — 100 rows, 30 columns
- `1_Trajectories/egfr-domain3-h370-r1-short_detarget_l60_bd77320f6718bc34/egfr-domain3-h370-r1-short_detarget_l60_bd77320f6718bc34_losses.csv` — 250 rows, 30 columns

## 4. What could not be determined, and why

- Nothing. Every table joined to an attempt and every attempt carried a decision.

## 5. What this does not establish

- Nothing here says a design works. It says what was tried, what survived each gate and why, which is the record needed to check the rates afterwards.
- The charge-pair verdicts are step 10's and carry step 10's limits: a predicted structure scored against a predicted interface, with the pH switch inferred from geometry and never measured.
