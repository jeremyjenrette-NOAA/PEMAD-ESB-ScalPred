library(dplyr)
library(ggplot2)
source("./gamfunc.R")   # also sources cvfunc.R
source("./procfunc.R")

# ======================================================================
# 0. Configuration
# ======================================================================
N_FOLDS          <- 5       # must match makedat.R
DET_ENGINE       <- "gam"   # "gam" (REML) recommended; "bam" only for image/block REs

# Fixed-effect formula selection:
#   "cv"  = every candidate cross-fitted on the block folds, OOF log loss,
#           one-SE rule toward lower edf (respects within-image/block clustering)
#   "aic" = previous behaviour (AIC on all of test_GAM_train; optional nested)
FORMULA_SELECT   <- "cv"
NESTED_SELECTION <- FALSE   # only used when FORMULA_SELECT = "aic"

# Random-effect structures to compare (names from RE_SPECS in cvfunc.R), chosen
# by the same OOF one-SE rule. Image/block REs gave no OOF gain in the
# 2026-10-05 run and are slow, so the default compares cruise-level terms only.
RE_COMPARE <- TRUE
RE_SPECS_TO_COMPARE <- c("none", "cruise", "cruise_fs")
RE_SPEC <- "auto"   # "auto" or a fixed name, e.g. "cruise_fs"

# Inference-only report of image/block/cruise SDs (one bam fit per region).
# Uses an additive fixed part: bam(discrete) diverges on tensor products.
RE_VARCOMP_REPORT <- TRUE
VARCOMP_FORMULA   <- "M05_Boxsize"

# Holdout predictions: "full" = full-data GAM (the deployed model);
# "fold_avg" = mean of the K fold models (same models as the OOF training inputs).
TEST_PRED_SOURCE <- "full"

# ======================================================================
# 1. Load Data
# ======================================================================
models <- readRDS("../data/processed/eval_2226.rds")
if (!"fold" %in% names(models[[1]]$det)) {
  stop("eval_2226.rds has no `fold` column. Re-run makedat.R first.")
}

# Standardize names if not already done in the RDS builder
names(models)[names(models) == "YOLO"] <- "YOLOv12"
names(models)[names(models) == "Cascade"] <- "Cascade R-CNN"

# ======================================================================
# 2. Define Candidate Detection-Level GAM Structures (Systematic Stepwise)
# ======================================================================
candidate_forms <- list(
  # STAGE 1: Additive Isotropic Splines (Proving independent covariate importance)
  M01_ConfOnly = truedetect ~ s(conf, k = 5),
  M02_Spatial  = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7),
  M03_Depth    = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5),
  M04_Optical  = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(field_of_view_sq_meter, k = 5),
  M05_Boxsize  = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(boxsize, k = 5),
  M06_Environ  = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(fluorometer_backscatter_ntu, k = 5),
  
  # STAGE 2: Tensor Products (Proving Confidence interacts with physical scale)
  M07_teDepth   = truedetect ~ te(conf, bottom_depth, k = 5) + s(latitude, longitude, k = 7),
  M08_teAlt     = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5),
  M09_teOptics  = truedetect ~ te(conf, altitude, k = 5) + te(conf, field_of_view_sq_meter, k = 5) + s(latitude, longitude, k = 7),
  M10_teBoxsize = truedetect ~ te(conf, boxsize, k = 5) + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7),
  
  # STAGE 3: The "Previous Best" & Fine-Tuning
  M11_PrevBest  = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7),
  M12_PB_Depth  = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(bottom_depth, k = 5),
  M13_PB_Backsc = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(fluorometer_backscatter_ntu, k = 5),
  M14_PB_Temp   = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(ctd_temperature_celsius, k = 5),
  M15_All_te    = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + te(conf, bottom_depth, k = 5)
)

model_names_to_run <- c("YOLOv12", "Cascade R-CNN")
re_summary <- list()

# ======================================================================
# 3. Iterative Pipeline
# ======================================================================
for (mod_name in model_names_to_run) {
  cat("\n========================================\n")
  cat("Running Pipeline for:", mod_name, "\n")
  cat("========================================\n")
  
  mod_data <- models[[mod_name]]
  safe_mod_name <- gsub(" ", "", mod_name)   # "CascadeR-CNN" (keeps existing file names)
  mod_data$det <- add_re_factors(mod_data$det)
  train_det <- det_fit_rows(mod_data$det %>% filter(dataset == "test_GAM_train"))
  
  # A. Optimal F1 thresholds per region, tuned on test_GAM_train ONLY
  cat("Evaluating PR on test_GAM_train to find optimal F1 thresholds...\n")
  temp_model_list <- list()
  temp_model_list[[mod_name]] <- list(
    img = mod_data$img %>% filter(dataset == "test_GAM_train"),
    det = mod_data$det %>% filter(dataset == "test_GAM_train")
  )
  out_pr <- evaluate_pr_models(temp_model_list, stratify_region = TRUE)
  best_pts <- out_pr$pr_all %>%
    group_by(model) %>%
    filter(f1 == max(f1, na.rm = TRUE)) %>%
    slice_max(conf, n = 1) %>% 
    ungroup()
  thresh_gb  <- best_pts %>% filter(grepl(" GB$", model)) %>% pull(conf)
  thresh_mab <- best_pts %>% filter(grepl(" MAB$", model)) %>% pull(conf)
  cat("  -> GB F1 Threshold:", thresh_gb, "\n")
  cat("  -> MAB F1 Threshold:", thresh_mab, "\n")
  f1_list <- list(GB = thresh_gb, MAB = thresh_mab)
  
  # B. Fixed-effect formula (no REs)
  if (FORMULA_SELECT == "cv") {
    cat("Selecting fixed-effect formula by block-fold CV (no REs)...\n")
    # Cached: the ~20-min CV is reused on re-runs. Delete the .rds to refit
    # (required after changing candidate_forms, folds, data or engine).
    cv_cache <- paste0("../data/processed/cache_formula_cv_", safe_mod_name, "_2226.rds")
    fcv <- if (file.exists(cv_cache)) readRDS(cv_cache) else NULL
    data_stamp <- as.numeric(file.mtime("../data/processed/eval_2226.rds"))
    if (!is.null(fcv) && (is.null(fcv$fold_diag) || nrow(fcv$fold_diag) == 0 || any(is.na(fcv$scores$complexity)))) {
      cat("  (cached formula CV is incomplete - refitting)\n"); fcv <- NULL
    }
    if (!is.null(fcv) && !identical(fcv$data_stamp, data_stamp)) {
      cat("  (eval_2226.rds changed since the formula CV was cached - refitting)\n"); fcv <- NULL
    }
    if (!is.null(fcv)) {
      cat("  (loading cached formula CV from", cv_cache, ")\n")
    } else {
      fcv <- cv_select_detection_formula(mod_data$det, candidate_forms, K = N_FOLDS, engine = DET_ENGINE)
      fcv$data_stamp <- data_stamp
      saveRDS(fcv, cv_cache)
    }
    write.csv(fcv$fold_diag, paste0("../data/processed/", safe_mod_name, "_formula_cv_folds_2226.csv"), row.names = FALSE)
    n_bad <- sum(!fcv$fold_diag$usable)
    if (n_bad > 0) {
      cat("  !!", n_bad, "unusable fold fits:\n")
      print(fcv$fold_diag %>% filter(!usable) %>% select(candidate, region, fold, converged, dev_expl, notes), n = Inf, width = Inf)
    }
    fsel <- choose_detection_formula(fcv)
    print(fsel$one_se %>% arrange(region, log_loss), n = Inf, width = Inf)
    write.csv(fsel$one_se, paste0("../data/processed/", safe_mod_name, "_formula_cv_2226.csv"), row.names = FALSE)
    base_ids <- fsel$chosen
  } else {
    cat("Selecting fixed-effect formula by AIC (no REs)...\n")
    base_sel <- lapply(c(GB = "GB", MAB = "MAB"), function(r)
      fit_best_detection_gam(train_det %>% filter(region == r), candidate_forms, engine = DET_ENGINE))
    for (r in names(base_sel)) { cat(" ", r, "\n"); print(base_sel[[r]]$table) }
    base_ids <- lapply(base_sel, `[[`, "id")
    rm(base_sel)
  }
  base_formulas <- lapply(base_ids, function(id) candidate_forms[[id]])
  cat("  -> Formula: GB =", base_ids$GB, "| MAB =", base_ids$MAB, "\n")
  
  # C. Random-effect structure: compare on OOF scores (one-SE rule)
  if (RE_COMPARE) {
    cat("Comparing random-effect structures (OOF, K =", N_FOLDS, ")...\n")
    re_cmp <- compare_re_specs(mod_data$det, base_formulas,
                               re_specs = RE_SPECS[RE_SPECS_TO_COMPARE],
                               K = N_FOLDS, engine = DET_ENGINE)
    print(re_cmp$one_se %>% arrange(region, log_loss), n = Inf, width = Inf)
    write.csv(re_cmp$one_se, paste0("../data/processed/", safe_mod_name, "_re_compare_2226.csv"), row.names = FALSE)
    write.csv(re_cmp$sds,    paste0("../data/processed/", safe_mod_name, "_re_sd_2226.csv"), row.names = FALSE)
    re_summary[[mod_name]] <- re_cmp$one_se %>% mutate(model = mod_name, .before = 1)
  }
  re_choice <- if (RE_SPEC == "auto") {
    if (!RE_COMPARE) stop("RE_SPEC = 'auto' requires RE_COMPARE = TRUE.")
    re_cmp$chosen
  } else list(GB = RE_SPEC, MAB = RE_SPEC)
  re_terms <- lapply(re_choice, function(nm) RE_SPECS[[nm]])
  cat("  -> RE structure: GB =", re_choice$GB, "| MAB =", re_choice$MAB, "\n")
  
  # D. Final detection GAMs on all of test_GAM_train
  cat("Fitting final detection GAMs...\n")
  if (FORMULA_SELECT == "cv") {
    final_ids <- base_ids
    gams <- lapply(c(GB = "GB", MAB = "MAB"), function(r) {
      fit <- try_fit_det_gam(with_re(base_formulas[[r]], re_terms[[r]]),
                             train_det %>% filter(region == r), engine = DET_ENGINE, label = paste("final", r))
      if (!fit_ok(fit)) stop("Final detection GAM unusable for ", r)
      fit
    })
  } else {
    # formula re-selected by AIC under the chosen RE structure
    final_sel <- lapply(c(GB = "GB", MAB = "MAB"), function(r)
      fit_best_detection_gam(train_det %>% filter(region == r), candidate_forms,
                             re_terms = re_terms[[r]], engine = DET_ENGINE))
    gams <- lapply(final_sel, `[[`, "fit")
    final_ids <- lapply(final_sel, `[[`, "id")
    rm(final_sel)
  }
  cat("  -> Final formulas: GB =", final_ids$GB, "| MAB =", final_ids$MAB, "\n")
  saveRDS(gams, file = paste0("../data/processed/", safe_mod_name, "_detgams_2226.rds"))
  
  # E. Cross-fit: OOF pred_p for test_GAM_train (+ fold-averaged test preds)
  cat("Cross-fitting detection GAMs (K =", N_FOLDS, ", nested selection =", NESTED_SELECTION, ")...\n")
  cf <- crossfit_detection_gams(
    det_df             = mod_data$det,
    candidate_formulas = candidate_forms,
    fixed_formulas     = list(GB = candidate_forms[[final_ids$GB]],
                              MAB = candidate_forms[[final_ids$MAB]]),
    K                  = N_FOLDS,
    nested_selection   = FORMULA_SELECT == "aic" && NESTED_SELECTION,
    re_terms           = re_terms,
    engine             = DET_ENGINE
  )
  mod_data$det <- cf$det
  print(cf$selection, width = Inf)
  write.csv(cf$selection, paste0("../data/processed/", safe_mod_name, "_cv_selection_2226.csv"), row.names = FALSE)
  
  # F. Test on holdout data
  cat("Testing GAMs on holdout data...\n")
  eval_results <- test_calibration_gams(
    gams             = gams, 
    det_df           = mod_data$det, 
    img_df           = mod_data$img, 
    f1_thresholds    = f1_list,
    test_pred_source = TEST_PRED_SOURCE
  )
  eval_results$re_choice <- re_choice
  eval_results$formulas  <- final_ids
  saveRDS(eval_results, file = paste0("../data/processed/", safe_mod_name, "_deteval_2226.rds"))
  
  print(eval_results$metrics, width = Inf)
  cat("Leakage check on test_GAM_train (in-sample vs out-of-fold sums):\n")
  print(eval_results$leakage_check, width = Inf)
  
  # G. Generate and Save Diagnostics
  cat("Generating and saving plots...\n")
  plots <- generate_gam_plots(eval_res = eval_results, model_name = mod_name, dataset_label = "2022-2024, 2026")
  
  # Figure file names drop hyphens as before
  safe_mod_name <- gsub(" |-", "", mod_name)
  
  ggsave(paste0("../figures/diag3/2226_", safe_mod_name, "_count_strat.png"), plot = plots$p_zoomed, width = 9, height = 5, bg = "white")
  ggsave(paste0("../figures/diag3/2226_", safe_mod_name, "_resid_strat.png"), plot = plots$p_resid, width = 5, height = 4.5, bg = "white")
  ggsave(paste0("../figures/diag3/2226_", safe_mod_name, "_fn_strat.png"), plot = plots$p_fn, width = 5, height = 4.5, bg = "white")
  ggsave(paste0("../figures/ms_figures/2226_", safe_mod_name, "_fn_strat.png"), plot = plots$p_fn, width = 5, height = 4.5, bg = "white")
  ggsave(paste0("../figures/diag3/2226_", safe_mod_name, "_sum_strat.png"), plot = plots$p_sum, width = 6, height = 5, bg = "white")
  
  # H. Variance components (inference only; does not feed the predictions)
  if (RE_VARCOMP_REPORT) {
    cat("Variance-component report (image/block/cruise SDs)...\n")
    vc_tab <- report_variance_components(train_det, list(GB = candidate_forms[[VARCOMP_FORMULA]],
                                                         MAB = candidate_forms[[VARCOMP_FORMULA]]))
    print(vc_tab, width = Inf)
    write.csv(vc_tab, paste0("../data/processed/", gsub(" ", "", mod_name), "_varcomp_2226.csv"), row.names = FALSE)
  }
}

if (length(re_summary)) {
  write.csv(bind_rows(re_summary), "../data/processed/re_compare_all_2226.csv", row.names = FALSE)
}
cat("\nPipeline complete! All plots saved to ../figures/diag3/\n")
