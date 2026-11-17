# Geographic Prediction of *Salmonella enterica* Serotype Typhimurium from Pork
Random Forest Model for Geographic Prediction of *Salmonella enterica* Serotype Typhimurium from Pork in China and the United States Trained Based on Gene Presence/Absence

This repository contains the R code for predicting the country of origin (*China* vs. *United States*) of pork-associated *Salmonella* Typhimurium isolates with a random forest trained on pan-genome gene presence/absence data, and for testing whether the prediction generalizes beyond clonal population structure.

## 📂 Scripts
- `Main Workflow.R` — **Main script** (feature filtering, random forest, random hold-out and lineage-aware validation, feature-importance stability, Figure 10). All settings are in the `cfg` list at the top.

## 🧠 What `Main Workflow.R` does
1. **Loads the data**: isolate table (`Isolates.xlsx`, sheet *Isolates*), the pan-genome gene presence/absence matrix (`Gene Presence Absence Matrix.xlsx`), and the 70 phylogenetic clusters (`Clusters.xlsx`, sheet *Clusters*). Place the three files in `cfg$data_dir`.
2. **Filters features**: near-zero-variance genes are removed with the default `caret::nearZeroVar` criteria (`freqCut = 95/5`, `uniqueCut = 10`), implemented as a vectorised function and **recomputed on every training set**. A gene is removed if the ratio of its most common to second most common state exceeds 19 **and** the percentage of unique values (100 × distinct states / number of isolates) is below 10, or if it is constant. For binary data this removes genes present in fewer than 5% or more than 95% of isolates (for 589 isolates, genes present in 30–559 isolates are retained). On the full data the script keeps **1,986 of 15,477** gene clusters (reported in `run_log.txt`).
3. **Fits random forests** with `ranger` (probability forests, 750 trees). **No hyperparameters are tuned**: `mtry = floor(sqrt(p))`, minimum node size 10 (the `ranger` default for probability forests) and the Gini splitting rule are the software defaults. The main model is a **balanced random forest**: every tree is grown on equal numbers of isolates from each country, drawn without replacement. Inverse-frequency class-weighted and unweighted forests are also fitted. Isolates with a predicted probability ≥ 0.5 for the United States are classified as United States.
4. **Validates under several designs**: stratified random hold-out (55% of each country for training, 45% for testing; 10 repeats), unstratified random hold-out, grouped 10-fold cross-validation holding out whole clusters (5 repeats), leave-one-cluster-out, and the same designs within ST19; plus lineage-only baselines (ST, tree cluster).
5. **Feature importance**: mean decrease in impurity within folds (stability across cluster-blocked folds), and permutation importance for the final model trained on all isolates.
6. **Produces a figure** (SVG and PDF).

## 📊 Outputs (`Results/`)
`run_log.txt`, `sessionInfo.txt`, `ml_summary.csv`, `ml_metrics_per_repeat.csv`, `ml_predictions_first_repeat.csv`, `ml_importance_stability_random_split.csv`, `ml_importance_stability_tree_cluster_blocked.csv`, `ml_importance_random_vs_tree_cluster_blocked.csv`, `ml_importance_full_model_permutation.csv`, `lineage_ST_by_country.csv`, `lineage_treecluster_by_country.csv`, `Figure_*.svg/pdf`, and a cached `pangenome_matrix.rds`.

## ▶️ Running it
```r
source("Main Workflow.R")   # set cfg$fast <- TRUE first for a quick smoke test
```
Expected values in `run_log.txt`: nearZeroVar keeps 1,986 of 15,477 clusters on the full data; Clusters: 70 clusters, 52 singletons, largest 313. Leave-one-cluster-out is slow (`cfg$do_loco`).

## 🧩 Dependencies
R packages: `readxl`, `dplyr`, `tidyr`, `readr`, `ranger`, `pROC`, `ggplot2`, `svglite`, and (optional) `patchwork`. Versions: R 4.5.3, ranger 0.18.0. A fixed random seed of 123 is used throughout (`cfg$seed`).

---

**Citation:**
This workflow was developed for binary classification of *Salmonella* Typhimurium isolates from pork based on pan-genomic features.
If used in a publication, please cite the corresponding methods description.
