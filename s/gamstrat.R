library(dplyr)
library(ggplot2)
source("./gamfunc.R")
source("./procfunc.R")
source("./cvfunc.R")

N_FOLDS          <- 5     # must match makedat.R
NESTED_SELECTION <- TRUE  # re-select the detection-GAM formula inside each fold

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

# ======================================================================
# 3. Model Selection (Example: Testing structures on YOLO for Georges Bank)
# ======================================================================
cat("--- Running Model Selection Phase ---\n")
selection_results_yolo_GB <- compare_detection_gams(
  det_df = models$YOLOv12$det, 
  candidate_formulas = candidate_forms, 
  region_focus = "GB"
)

selection_results_yolo_MAB <- compare_detection_gams(
  det_df = models$YOLOv12$det, 
  candidate_formulas = candidate_forms, 
  region_focus = "MAB"
)

selection_results_cascade_GB <- compare_detection_gams(
  det_df = models$`Cascade R-CNN`$det, 
  candidate_formulas = candidate_forms, 
  region_focus = "GB"
)

selection_results_cascade_MAB <- compare_detection_gams(
  det_df = models$`Cascade R-CNN`$det, 
  candidate_formulas = candidate_forms, 
  region_focus = "MAB"
)

print(selection_results_yolo_GB)
print(selection_results_yolo_MAB)
print(selection_results_cascade_GB)
print(selection_results_cascade_MAB)

# Final (full test_GAM_train) formulas: lowest AIC, the same rule applied inside
# each cross-fitting fold, so the final model and the OOF models are selected
# consistently. Inspect the printed tables; override below if needed.
best_formula_yolo_gb     <- pick_best_formula(selection_results_yolo_GB, candidate_forms)
best_formula_yolo_mab    <- pick_best_formula(selection_results_yolo_MAB, candidate_forms)
best_formula_cascade_gb  <- pick_best_formula(selection_results_cascade_GB, candidate_forms)
best_formula_cascade_mab <- pick_best_formula(selection_results_cascade_MAB, candidate_forms)
# Previous manual choices:
# best_formula_yolo_gb  <- candidate_forms$M15_All_te 
# best_formula_yolo_mab <- candidate_forms$M15_All_te
# best_formula_cascade_gb  <- candidate_forms$M14_PB_Temp 
# best_formula_cascade_mab <- candidate_forms$M14_PB_Temp

model_names_to_run <- c("YOLOv12", "Cascade R-CNN")

# ======================================================================
# 3. Iterative Pipeline
# ======================================================================
for (mod_name in model_names_to_run) {
  cat("\n========================================\n")
  cat("Running Pipeline for:", mod_name, "\n")
  cat("========================================\n")
  
  mod_data <- models[[mod_name]]
  
  # A. Dynamically Calculate Optimal F1 Thresholds per Region
  #    Tuned on test_GAM_train ONLY so the F1 baseline never sees holdout images.
  cat("Evaluating PR on test_GAM_train to find optimal F1 thresholds...\n")
  temp_model_list <- list()
  temp_model_list[[mod_name]] <- list(
    img = mod_data$img %>% filter(dataset == "test_GAM_train"),
    det = mod_data$det %>% filter(dataset == "test_GAM_train")
  )
  
  # We use evaluate_pr_models to get the curves, then filter for the max F1
  out_pr <- evaluate_pr_models(temp_model_list, stratify_region = TRUE)
  best_pts <- out_pr$pr_all %>%
    group_by(model) %>%
    filter(f1 == max(f1, na.rm = TRUE)) %>%
    slice_max(conf, n = 1) %>% 
    ungroup()
  
  # Extract specific thresholds safely
  thresh_gb  <- best_pts %>% filter(grepl(" GB$", model)) %>% pull(conf)
  thresh_mab <- best_pts %>% filter(grepl(" MAB$", model)) %>% pull(conf)
  
  cat("  -> GB F1 Threshold:", thresh_gb, "\n")
  cat("  -> MAB F1 Threshold:", thresh_mab, "\n")
  
  f1_list <- list(GB = thresh_gb, MAB = thresh_mab)
  
  if (mod_name == "YOLOv12") {
    formula_gb <- best_formula_yolo_gb;    formula_mab <- best_formula_yolo_mab
  } else {
    formula_gb <- best_formula_cascade_gb; formula_mab <- best_formula_cascade_mab
  }
  
  # B. Train Calibration GAMs on all of test_GAM_train (used for test_GAM_test)
  gams <- fit_calibration_gams(
    det_df = mod_data$det,
    formula_gb = formula_gb,
    formula_mab = formula_mab
  )
  
  # B2. Cross-fit: out-of-fold pred_p for every test_GAM_train detection.
  #     These OOF sums (not in-sample fitted values) feed the synergy GAM.
  cat("Cross-fitting detection GAMs (K =", N_FOLDS, ", nested selection =", NESTED_SELECTION, ")...\n")
  cf <- crossfit_detection_gams(
    det_df             = mod_data$det,
    candidate_formulas = candidate_forms,
    fixed_formulas     = list(GB = formula_gb, MAB = formula_mab),
    K                  = N_FOLDS,
    nested_selection   = NESTED_SELECTION
  )
  mod_data$det <- cf$det
  print(cf$selection)
  write.csv(cf$selection, paste0("../data/processed/", gsub(" ", "", mod_name),
                                 "_cv_selection_2226.csv"), row.names = FALSE)
  
  if (mod_name == "Cascade R-CNN") mod_name <- gsub("\\s", "", mod_name)
  saveRDS(gams, file = paste0("../data/processed/", mod_name, 
                                       "_detgams_2226.rds"))
  
  # C. Test GAMs on Holdout Data
  cat("Testing GAMs on holdout data...\n")
  eval_results <- test_calibration_gams(
    gams          = gams, 
    det_df        = mod_data$det, 
    img_df        = mod_data$img, 
    f1_thresholds = f1_list
  )
  
  if (mod_name == "Cascade R-CNN") mod_name <- gsub("\\s", "", mod_name)
  saveRDS(eval_results, file = paste0("../data/processed/", mod_name, 
                              "_deteval_2226.rds"))
  
  print(eval_results$metrics)
  cat("Leakage check on test_GAM_train (in-sample vs out-of-fold sums):\n")
  print(eval_results$leakage_check)
  
  # D. Generate and Save Diagnostics
  cat("Generating and saving plots...\n")
  plots <- generate_gam_plots(eval_res = eval_results, model_name = mod_name, dataset_label = "2022-2024, 2026")
  
  # Safe filename string (removes spaces/special characters)
  safe_mod_name <- gsub(" |-", "", mod_name)
  
  ggsave(paste0("../figures/diag3/2226_", safe_mod_name, "_count_strat.png"), plot = plots$p_zoomed, width = 9, height = 5, bg = "white")
  ggsave(paste0("../figures/diag3/2226_", safe_mod_name, "_resid_strat.png"), plot = plots$p_resid, width = 5, height = 4.5, bg = "white")
  ggsave(paste0("../figures/diag3/2226_", safe_mod_name, "_fn_strat.png"), plot = plots$p_fn, width = 5, height = 4.5, bg = "white")
  ggsave(paste0("../figures/ms_figures/2226_", safe_mod_name, "_fn_strat.png"), plot = plots$p_fn, width = 5, height = 4.5, bg = "white")
  ggsave(paste0("../figures/diag3/2226_", safe_mod_name, "_sum_strat.png"), plot = plots$p_sum, width = 6, height = 5, bg = "white")
}

cat("\nPipeline complete! All plots saved to ../figures/diag3/\n")
