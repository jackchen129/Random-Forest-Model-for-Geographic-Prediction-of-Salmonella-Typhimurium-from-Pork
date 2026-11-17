# =====================================================================
# Random forest geographic classifier: near-zero-variance filtering and lineage-aware
# validation.
# Pork-associated Salmonella Typhimurium, China (n = 491) vs United States (n = 98)
#
# What it does
#   1. Loads the isolate table (Isolates.xlsx), the pan-genome matrix and the 70 clusters (Clusters.xlsx)
#   2. Random forest (balanced, 750 trees) under: random 0.55:0.45 splits (stratified and unstratified);
#      grouped 10-fold CV and leave-one-cluster-out using the 70 clusters (Clusters.xlsx);
#      the same inside ST19; class-weighted / unweighted forests;
#      lineage-only baselines (ST, tree cluster)
#   3. Feature-importance stability across blocked folds
#   4. Consolidated a figure (SVG + PDF)
#
# Settings (no hyperparameter tuning is performed in any design)
#   ranger probability forests, 750 trees, default mtry = floor(sqrt(p)), default minimum node size
#   (10 for probability forests), Gini splitting rule (default). Main model = balanced forest: every
#   tree is grown on equal numbers of isolates from each country, drawn without replacement.
#   Feature importance: mean decrease in impurity within folds; permutation importance for the
#   final model trained on all isolates.
# Feature filtering: caret::nearZeroVar defaults (freqCut = 95/5, uniqueCut = 10), vectorised and
#   recomputed on every training set. On the full data it should keep 1,986 of 15,477 gene clusters
#   (reported in run_log.txt).
#
# Usage: put the three input files in cfg$data_dir, then source() this script or run
#        `Rscript Main Workflow.R`. Outputs go to cfg$out_dir
#        (including run_log.txt). Set cfg$fast <- TRUE for a quick smoke test.
# =====================================================================

## 0. CONFIG ----------------------------------------------------------
cfg <- list(
  data_dir = ".",                    # folder holding the input files
  # input file names (spaces or underscores both work; any numeric upload prefix is ignored).
  f_s1     = "Isolates.xlsx",
  f_gpa    = "Gene Presence Absence Matrix.xlsx",
  f_s6     = "Clusters.xlsx",        # the 70 clusters (sheet "Clusters")
  out_dir  = "Random Forest",
  seed     = 123,
  threads  = max(1, parallel::detectCores() - 1),
  fast     = FALSE,                  # TRUE = tiny runs (smoke test only)
  num_trees       = 750,
  nzv_freq_ratio  = 95/5,            # caret::nearZeroVar default freqCut
  nzv_unique_cut  = 10,              # caret::nearZeroVar default uniqueCut
  k_folds         = 10,              # grouped k-fold (lineage-blocked)
  cv_repeats      = 5,
  primary_blocked = "tree_cluster_blocked",   # design shown in Figure / importance stability (the 70 clusters, Clusters.xlsx)
  holdout_repeats = 10,              # repeats of the random hold-out designs
  do_loco         = TRUE             # leave-one-cluster-out (slow)
)
run <- list(ml = TRUE)

if (cfg$fast) {
  cfg$k_folds <- 5; cfg$cv_repeats <- 1; cfg$holdout_repeats <- 2
  cfg$do_loco <- FALSE; cfg$num_trees <- 150
}

need <- c("readxl", "dplyr", "tidyr", "readr", "ranger", "pROC", "ggplot2", "svglite")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(cfg$seed)

log_file <- file.path(cfg$out_dir, "run_log.txt"); if (file.exists(log_file)) file.remove(log_file)
logmsg   <- function(...) {                       # prints to the console AND appends to revision_results/run_log.txt
  m <- paste0(format(Sys.time(), "%H:%M:%S "), ...)
  message(m); cat(m, "\n", file = log_file, append = TRUE, sep = "")
}
save_csv <- function(df, name) {
  readr::write_csv(as.data.frame(df), file.path(cfg$out_dir, name)); invisible(df)
}
save_plot <- function(p, name, w = 7, h = 5) {       # name without extension -> .svg + .pdf (both vector)
  base <- sub("\\.(png|pdf|svg)$", "", name)
  svg_dev <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else grDevices::svg   # svglite preferred; base svg() needs cairo
  ggplot2::ggsave(file.path(cfg$out_dir, paste0(base, ".svg")), p, width = w, height = h, device = svg_dev)
  ggplot2::ggsave(file.path(cfg$out_dir, paste0(base, ".pdf")), p, width = w, height = h, device = "pdf")
}
find_file <- function(name, rx = NA, required = TRUE) {
  # exact file name from cfg, ignoring case, a numeric upload prefix, and spaces vs underscores
  f <- list.files(cfg$data_dir, full.names = TRUE); f <- f[!grepl("^~\\$", basename(f))]
  norm <- function(x) tolower(gsub("[ _]+", "_", sub("^[0-9]+_", "", x)))
  hit <- f[norm(basename(f)) == norm(name)]
  if (!length(hit) && !is.na(rx)) hit <- f[grepl(rx, basename(f), ignore.case = TRUE)]
  if (!length(hit)) {
    if (required) stop("Input file '", name, "' not found in ", normalizePath(cfg$data_dir))
    return(NA_character_)
  }
  hit[1]
}


## 1. DATA LOADING ----------------------------------------------------
logmsg("1. Loading data")
f_s1    <- find_file(cfg$f_s1)
f_gpa   <- find_file(cfg$f_gpa)
f_s6    <- find_file(cfg$f_s6)

# --- Isolates.xlsx: isolates (title row 1, header row 2) ---
iso <- readxl::read_excel(f_s1, sheet = "Isolates", skip = 1) |> filter(!is.na(Country))
iso <- iso |> transmute(
  Isolate = Isolate,
  country = factor(ifelse(Country == "United States", "US", Country), levels = c("China", "US")),
  ST      = ifelse(grepl("^NT", `Sequence type`), "NT", as.character(`Sequence type`)))
ids <- iso$Isolate

# --- pan-genome matrix (cached after first read; the xlsx is slow) ---
cache <- file.path(cfg$out_dir, "pangenome_matrix.rds")
if (file.exists(cache)) {
  X <- readRDS(cache)
} else {
  hdr <- readxl::read_excel(f_gpa, n_max = 1)
  g   <- readxl::read_excel(f_gpa, col_types = c("text", rep("numeric", ncol(hdr) - 1)))
  X   <- t(as.matrix(g[, -1])); colnames(X) <- g[[1]]
  X[is.na(X)] <- 0; X[X > 0] <- 1; storage.mode(X) <- "integer"
  saveRDS(X, cache)
}
stopifnot(all(ids %in% rownames(X)))
X <- X[ids, , drop = FALSE]

# --- the 70 clusters (Clusters.xlsx, sheet "Clusters": Isolate, Cluster, Cluster type) ---
sh <- readxl::excel_sheets(f_s6); sh <- if ("Clusters" %in% sh) "Clusters" else sh[length(sh)]
raw <- readxl::read_excel(f_s6, sheet = sh, col_names = FALSE)
hr  <- which(apply(raw, 1, function(r) any(tolower(as.character(r)) == "isolate")))[1]
tcs <- raw[(hr + 1):nrow(raw), ]; names(tcs) <- as.character(unlist(raw[hr, ]))
tc6 <- as.character(tcs$Cluster[match(iso$Isolate, tcs$Isolate)])
stopifnot(!anyNA(tc6))
iso$tree_cluster <- paste0("C", tc6)
logmsg("   Clusters.xlsx clusters: ", n_distinct(iso$tree_cluster), " clusters, ",
       sum(table(iso$tree_cluster) == 1), " singletons, largest ", max(table(iso$tree_cluster)))

# --- descriptive lineage tables (support the choice of blocking) ---
save_csv(iso |> count(ST, country) |> pivot_wider(names_from = country, values_from = n, values_fill = 0), "lineage_ST_by_country.csv")
save_csv(iso |> count(tree_cluster, country) |> pivot_wider(names_from = country, values_from = n, values_fill = 0) |>
           arrange(desc(China + US)), "lineage_treecluster_by_country.csv")
logmsg("   ", nrow(iso), " isolates: China ", sum(iso$country == "China"), " / US ", sum(iso$country == "US"), "; ",
       n_distinct(iso$ST), " STs; ", n_distinct(iso$tree_cluster), " clusters")

## 2. LINEAGE-BLOCKED RANDOM FOREST -------------------------
# Designs compared: (i) random hold-out (0.55:0.45, stratified by
# country; plus an unstratified 0.55:0.45 split), (ii) grouped 10-fold CV holding out whole clusters (the 70 clusters of Clusters.xlsx),
# (iii) leave-one-cluster-out, (iv) the same inside ST19, the only ST with many isolates from
# both countries. Leave-one-ST-out is NOT run: 98.3% of isolates are ST34 or ST19, so holding
# out an ST removes one country's main lineage.
if (run$ml) {
  logmsg("2. Random forest: lineage-blocked validation")
  y <- iso$country
  # caret::nearZeroVar defaults (freqCut 95/5, uniqueCut 10), vectorised for speed
  keep_cols <- function(M) {
    n <- nrow(M); n1 <- colSums(M == 1L); n0 <- n - n1
    second <- pmin(n1, n0); fr <- pmax(n1, n0) / second
    pu <- 100 * ((n1 > 0) + (n0 > 0)) / n
    second > 0 & !(fr > cfg$nzv_freq_ratio & pu < cfg$nzv_unique_cut)
  }
  logmsg("   nearZeroVar (freqCut=", round(cfg$nzv_freq_ratio, 2), ", uniqueCut=", cfg$nzv_unique_cut,
         ") keeps ", sum(keep_cols(X)), " of ", ncol(X), " clusters on the full data (1,986)")
  
  # mode: "balanced"   = balanced random forest (each tree sees equal numbers of both classes; MAIN model)
  #       "classweight"= ranger class.weights (inverse frequency) - inverse-frequency weights; in
  #                      probability forests this changes predictions very little, so it is reported as a check
  #       "unweighted" = plain forest
  fit_rf <- function(tr, mode = "balanced", importance = "impurity", seed = 1) {
    keep <- keep_cols(X[tr, , drop = FALSE]); ytr <- y[tr]
    a <- list(x = X[tr, keep, drop = FALSE], y = ytr, num.trees = cfg$num_trees, probability = TRUE,
              importance = importance, seed = seed, num.threads = cfg$threads)
    if (mode == "classweight") a$class.weights <- as.numeric(length(ytr) / (2 * table(ytr)))
    if (mode == "balanced") { a$sample.fraction <- rep(min(table(ytr)) - 0.5, 2) / length(ytr); a$replace <- FALSE }  # -0.5 avoids float rounding above class size
    list(rf = do.call(ranger::ranger, a), keep = keep)
  }
  pred_rf <- function(fit, te) predict(fit$rf, X[te, fit$keep, drop = FALSE],
                                       num.threads = cfg$threads)$predictions[, "US"]
  lookup_p <- function(grp, tr, te) {            # lineage-only baseline: US rate of the group in ing
    prior <- mean(y[tr] == "US")
    tab <- tapply(y[tr] == "US", grp[tr], function(v) (sum(v) + prior) / (length(v) + 1))
    p <- as.numeric(tab[as.character(grp[te])]); p[is.na(p)] <- prior; p
  }
  metrics <- function(truth, p, thr = 0.5) {
    pred <- ifelse(p >= thr, "US", "China")
    tp <- sum(pred == "US" & truth == "US"); tn <- sum(pred == "China" & truth == "China")
    fp <- sum(pred == "US" & truth == "China"); fn <- sum(pred == "China" & truth == "US")
    n <- length(truth); acc <- (tp + tn) / n
    r_us <- tp / (tp + fn); r_cn <- tn / (tn + fp)
    pe <- ((tp + fp) * (tp + fn) + (tn + fn) * (tn + fp)) / n^2
    mcc <- (tp * tn - fp * fn) / sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    auc <- if (length(unique(truth)) == 2)
      as.numeric(pROC::auc(pROC::roc(truth, p, levels = c("China", "US"), direction = "<", quiet = TRUE))) else NA
    data.frame(n = n, accuracy = acc, balanced_accuracy = mean(c(r_us, r_cn), na.rm = TRUE),
               recall_China = r_cn, recall_US = r_us, kappa = (acc - pe) / (1 - pe),
               MCC = mcc, AUC = auc, majority_baseline = max(table(truth)) / n)
  }
  make_splits <- function(sc, idx, rep_i) {
    set.seed(cfg$seed + 1000 * rep_i)
    if (sc$type == "holdout") {
      tr <- if (sc$strat) unlist(lapply(split(idx, y[idx]), function(v) sample(v, round(sc$frac * length(v)))))
      else sample(idx, round(sc$frac * length(idx)))
      return(list(list(tr = tr, te = setdiff(idx, tr))))
    }
    g <- as.character(iso[[sc$block]][idx])
    if (sc$type == "loco") return(lapply(unique(g), function(u) list(tr = idx[g != u], te = idx[g == u])))
    gs <- data.frame(g = unique(g)); gs$n <- as.vector(table(g)[gs$g])
    gs$us <- as.vector(tapply(y[idx] == "US", g, sum)[gs$g])
    gs <- gs[order(-gs$n * runif(nrow(gs), 0.7, 1.3)), ]            # big clusters first, jittered per repeat
    fs <- rep(0, cfg$k_folds); fu <- fs; asg <- integer(nrow(gs))
    for (i in seq_len(nrow(gs))) {                                     # greedy: keep US count and size balanced
      f <- which.min(fs + 4 * fu); asg[i] <- f; fs[f] <- fs[f] + gs$n[i]; fu[f] <- fu[f] + gs$us[i]
    }
    fold <- setNames(asg, gs$g)[g]
    Filter(function(z) length(z$te) > 0, lapply(seq_len(cfg$k_folds), function(f) list(tr = idx[fold != f], te = idx[fold == f])))
  }
  
  schemes <- list(
    random_stratified_55   = list(type = "holdout", strat = TRUE,  frac = 0.55, reps = cfg$holdout_repeats, unweighted = TRUE),   
    random_unstratified_55 = list(type = "holdout", strat = FALSE, frac = 0.55, reps = cfg$holdout_repeats, unweighted = FALSE),  # unstratified random split
    tree_cluster_blocked   = list(type = "kfold", block = "tree_cluster", reps = cfg$cv_repeats, unweighted = TRUE))
  if (cfg$do_loco) {
    schemes$loco_tree_cluster <- list(type = "loco", block = "tree_cluster", reps = 1, unweighted = FALSE)
  }
  
  imp_store <- list()
  run_scheme <- function(nm, sc, idx) {
    mets <- list(); preds <- list()
    for (r in seq_len(sc$reps)) {
      pool <- list(te = integer(0), RF_balanced = numeric(0), RF_classweight = numeric(0),
                   RF_unweighted = numeric(0), lookup_ST = numeric(0), lookup_tree_cluster = numeric(0))
      for (s in make_splits(sc, idx, r)) {
        f <- fit_rf(s$tr, "balanced", seed = r)
        imp_store[[nm]][[length(imp_store[[nm]]) + 1]] <<- f$rf$variable.importance
        pool$te <- c(pool$te, s$te); pool$RF_balanced <- c(pool$RF_balanced, pred_rf(f, s$te))
        if (sc$unweighted) {
          pool$RF_classweight <- c(pool$RF_classweight, pred_rf(fit_rf(s$tr, "classweight", seed = r), s$te))
          pool$RF_unweighted  <- c(pool$RF_unweighted,  pred_rf(fit_rf(s$tr, "unweighted",  seed = r), s$te))
        }
        pool$lookup_ST <- c(pool$lookup_ST, lookup_p(iso$ST, s$tr, s$te))
        pool$lookup_tree_cluster <- c(pool$lookup_tree_cluster, lookup_p(iso$tree_cluster, s$tr, s$te))
      }
      for (m in c("RF_balanced", "RF_classweight", "RF_unweighted", "lookup_ST", "lookup_tree_cluster")) {
        if (!length(pool[[m]])) next
        mets[[length(mets) + 1]] <- cbind(scheme = nm, model = m, rep = r, metrics(y[pool$te], pool[[m]]))
        if (r == 1) preds[[length(preds) + 1]] <- data.frame(scheme = nm, model = m, Isolate = ids[pool$te],
                                                             truth = y[pool$te], p_US = pool[[m]])
      }
    }
    list(metrics = bind_rows(mets), preds = bind_rows(preds))
  }
  
  res <- list()
  all_idx <- seq_len(nrow(iso))
  for (nm in names(schemes)) { logmsg("   scheme: ", nm); res[[nm]] <- run_scheme(nm, schemes[[nm]], all_idx) }
  # same designs restricted to ST19 (only ST with many isolates from both countries)
  idx19 <- which(iso$ST == "19")
  for (nm in c("random_stratified_55", "tree_cluster_blocked")) {
    nn <- paste0("ST19only_", nm); logmsg("   scheme: ", nn)
    sc <- schemes[[nm]]; sc$unweighted <- FALSE
    res[[nn]] <- run_scheme(nn, sc, idx19)
  }
  
  ml_metrics <- bind_rows(lapply(res, `[[`, "metrics"))
  ml_preds   <- bind_rows(lapply(res, `[[`, "preds"))
  ml_summary <- ml_metrics |> group_by(scheme, model) |>
    summarise(across(c(accuracy, balanced_accuracy, recall_China, recall_US, kappa, MCC, AUC, majority_baseline),
                     list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)), .names = "{.col}_{.fn}"),
              n_repeats = n(), .groups = "drop")
  save_csv(ml_metrics, "ml_metrics_per_repeat.csv"); save_csv(ml_summary, "ml_summary.csv")
  save_csv(ml_preds, "ml_predictions_first_repeat.csv")
  print(ml_summary |> transmute(scheme, model, acc = round(accuracy_mean, 3), bal_acc = round(balanced_accuracy_mean, 3),
                                AUC = round(AUC_mean, 3), recall_US = round(recall_US_mean, 3),
                                baseline = round(majority_baseline_mean, 3)), n = 50)
  
  # --- feature importance: stability across lineage-blocked folds, and vs random split ---
  imp_mat <- function(lst) { f <- unique(unlist(lapply(lst, names)))
  sapply(lst, function(v) { o <- setNames(rep(NA_real_, length(f)), f); o[names(v)] <- v; o }) }
  summarise_imp <- function(lst, top = 20) {
    M <- imp_mat(lst); rk <- apply(M, 2, function(v) rank(-v, na.last = "keep", ties.method = "min"))
    data.frame(feature = rownames(M), mean_importance = rowMeans(M, na.rm = TRUE),
               folds_in_top = rowSums(rk <= top, na.rm = TRUE), folds_present = rowSums(!is.na(M)),
               stringsAsFactors = FALSE) |> arrange(desc(folds_in_top), desc(mean_importance))
  }
  imp_random <- summarise_imp(imp_store$random_stratified_55); imp_random$annotated <- !grepl("^group_", imp_random$feature)
  save_csv(imp_random, "ml_importance_stability_random_split.csv")
  imp_by <- list()
  for (nm in "tree_cluster_blocked") {
    ib <- summarise_imp(imp_store[[nm]]); ib$annotated <- !grepl("^group_", ib$feature); imp_by[[nm]] <- ib
    save_csv(ib, paste0("ml_importance_stability_", nm, ".csv"))
    common <- intersect(ib$feature, imp_random$feature)
    top_b <- head(arrange(ib, desc(mean_importance))$feature, 20); top_r <- head(arrange(imp_random, desc(mean_importance))$feature, 20)
    ic <- data.frame(design = nm, spearman_mean_importance = cor(ib$mean_importance[match(common, ib$feature)],
                                                                 imp_random$mean_importance[match(common, imp_random$feature)], method = "spearman"),
                     top20_overlap = length(intersect(top_b, top_r)))
    save_csv(ic, paste0("ml_importance_random_vs_", nm, ".csv")); print(ic)
  }
  imp_blocked <- imp_by[[cfg$primary_blocked]]
  
  # --- final model on all isolates, permutation importance ---
  logmsg("   final model (all isolates, permutation importance)")
  ff <- fit_rf(all_idx, "balanced", importance = "permutation", seed = cfg$seed)
  imp_full <- data.frame(feature = names(ff$rf$variable.importance), importance = ff$rf$variable.importance) |>
    arrange(desc(importance)) |> mutate(annotated = !grepl("^group_", feature))
  save_csv(imp_full, "ml_importance_full_model_permutation.csv")
  top_unannotated <- head(imp_full$feature[!imp_full$annotated], 50)
  
  # --- consolidated figure: confusion matrix, ROC, probabilities, importance ---
  # Every panel and the combined figure are saved as SVG and PDF (vector).
  lab_primary <- "Cluster-blocked CV"
  pb <- ml_preds |> filter(scheme == cfg$primary_blocked, model == "RF_balanced")
  pr <- ml_preds |> filter(scheme == "random_stratified_55", model == "RF_balanced")   # Training/Testing split design (0.55:0.45, stratified by country)
  cm <- pb |> mutate(pred = factor(ifelse(p_US >= 0.5, "US", "China"), levels = c("China", "US"))) |>
    count(truth, pred) |> group_by(truth) |> mutate(pct_of_actual = round(100 * n / sum(n), 1)) |> ungroup()
  roc_df <- function(d, lab) { r <- pROC::roc(d$truth, d$p_US, levels = c("China", "US"), direction = "<", quiet = TRUE)
  data.frame(design = sprintf("%s (AUC %.3f)", lab, as.numeric(pROC::auc(r))),
             one_minus_specificity = 1 - r$specificities, sensitivity_US = r$sensitivities) }
  roc_all <- bind_rows(roc_df(pr, "Random split 55:45"), roc_df(pb, lab_primary))
  imp_top <- imp_blocked |> arrange(desc(mean_importance)) |> head(15)
  
  pA <- ggplot(cm, aes(pred, truth, fill = n)) + geom_tile(colour = "white") +
    geom_text(aes(label = paste0(n, "\n(", pct_of_actual, "%)")), size = 4.5) +
    scale_fill_gradient(low = "white", high = "#2171B5", guide = "none") +
    labs(title = paste0("A  Confusion matrix (", lab_primary, ")"), x = "Predicted", y = "Actual") + theme_minimal()
  pB <- ggplot(roc_all, aes(one_minus_specificity, sensitivity_US, colour = design)) +
    geom_line(linewidth = 1) + geom_abline(lty = 2, colour = "grey50") +
    labs(title = "B  ROC", x = "1 - specificity", y = "Sensitivity (US)", colour = NULL) + theme_minimal() +
    theme(legend.position = "bottom")
  pC <- ggplot(pb, aes(p_US, fill = truth)) + geom_density(alpha = 0.5) +
    labs(title = "C  Predicted P(US), blocked CV", x = "Probability of US", y = "Density", fill = "Actual") + theme_minimal()
  pD <- ggplot(imp_top, aes(reorder(feature, mean_importance), mean_importance, fill = annotated)) + geom_col() + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#984EA3", `FALSE` = "grey60"), labels = c(`FALSE` = "unannotated", `TRUE` = "annotated")) +
    labs(title = "D  Top 15 features (blocked folds)", x = NULL, y = "Mean impurity importance", fill = NULL) + theme_minimal()
  panels <- list(A_confusion_matrix = pA, B_ROC_curve = pB, C_predicted_probability = pC, D_top15_importance = pD)
  for (nm in names(panels)) save_plot(panels[[nm]], paste0("Figure_panel", nm), 6, 5)
  if (has_patchwork) save_plot((pA | pB) / (pC | pD), "Figure_consolidated", 11, 9)
}


writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "sessionInfo.txt"))
logmsg("Done. Results in ", normalizePath(cfg$out_dir))
