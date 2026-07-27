library(dplyr)
library(ggplot2)
source("./gamfunc.R")
source("./procfunc.R")

# ======================================================================
# 1. Load Data
# ======================================================================
models <- readRDS("../data/processed/crab_evalmulti_2426.rds")

# Standardize names if not already done in the RDS builder
names(models)[names(models) == "YOLO"] <- "YOLOv12"
names(models)[names(models) == "Cascade"] <- "Cascade R-CNN"

# ======================================================================
# 2. Define Candidate Detection-Level GAM Structures (Systematic Stepwise)
# ======================================================================
candidate_forms <- list(
  # STAGE 1: Additive Isotropic Splines (Proving independent covariate importance)
  M01_ConfOnly = truedetect ~ spname + s(conf, k = 5),
  M02_Spatial  = truedetect ~ spname + s(conf, k = 5) + s(latitude, longitude, k = 7),
  M03_Depth    = truedetect ~ spname + s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5),
  M04_Optical  = truedetect ~ spname + s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(field_of_view_sq_meter, k = 5),
  M05_Boxsize  = truedetect ~ spname + s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(boxsize, k = 5),
  M06_Environ  = truedetect ~ spname + s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(fluorometer_backscatter_ntu, k = 5),
  
  # STAGE 2: Tensor Products (Proving Confidence interacts with physical scale)
  M07_teDepth   = truedetect ~ spname + te(conf, bottom_depth, k = 5) + s(latitude, longitude, k = 7),
  M08_teAlt     = truedetect ~ spname + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5),
  M09_teOptics  = truedetect ~ spname + te(conf, altitude, k = 5) + te(conf, field_of_view_sq_meter, k = 5) + s(latitude, longitude, k = 7),
  M10_teBoxsize = truedetect ~ spname + te(conf, boxsize, k = 5) + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7),
  
  # STAGE 3: The "Previous Best" & Fine-Tuning
  M11_PrevBest  = truedetect ~ spname + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7),
  M12_PB_Depth  = truedetect ~ spname + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(bottom_depth, k = 5),
  M13_PB_Backsc = truedetect ~ spname + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(fluorometer_backscatter_ntu, k = 5),
  M14_PB_Temp   = truedetect ~ spname + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(ctd_temperature_celsius, k = 5),
  M15_All_te    = truedetect ~ spname + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + te(conf, bottom_depth, k = 5)
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

print(selection_results_yolo_GB, n = Inf)
print(selection_results_yolo_MAB, n = Inf)
print(selection_results_cascade_GB, n = Inf)
print(selection_results_cascade_MAB, n = Inf)

# Pick the best formulas based on the selection results! 
# (You could automate this, but it is safer to manually inspect the AIC table and define them here)
best_formula_yolo_gb  <- candidate_forms$M15_All_te
best_formula_yolo_mab <- candidate_forms$M15_All_te
best_formula_cascade_gb  <- candidate_forms$M14_PB_Temp
best_formula_cascade_mab <- candidate_forms$M12_PB_Depth

# ======================================================================
# 3. Multiclass Iterative Pipeline (Line 81 Onward)
# ======================================================================
model_names_to_run <- c("YOLOv12", "Cascade R-CNN")
target_classes     <- c("jonah_crab", "rock_crab", "cancer_sp")

for (mod_name in model_names_to_run) {
  cat("\n========================================\n")
  cat("Running Multiclass Pipeline for:", mod_name, "\n")
  cat("========================================\n")
  
  mod_data <- models[[mod_name]]
  
  # --------------------------------------------------------------------
  # A. Dynamically Calculate Optimal Species & Region F1 Thresholds
  # --------------------------------------------------------------------
  cat("Evaluating PR curves for optimal species-level F1 thresholds...\n")
  temp_model_list <- list()
  temp_model_list[[mod_name]] <- mod_data
  
  # Calculate species-stratified PR curves per region
  out_pr <- evaluate_pr_models(temp_model_list, stratify_region = TRUE)
  
  # Parse out best threshold per (Region, Species) combination
  f1_thresh_df <- out_pr$pr_all %>%
    mutate(
      region = case_when(
        grepl(" GB$", model)  ~ "GB",
        grepl(" MAB$", model) ~ "MAB",
        TRUE ~ "ALL"
      )
    ) %>%
    group_by(region, species) %>%
    filter(f1 == max(f1, na.rm = TRUE)) %>%
    slice_max(conf, n = 1) %>%
    ungroup() %>%
    select(region, species, conf, precision, recall, f1)
  
  cat("\nOptimal Species F1 Thresholds:\n")
  print(f1_thresh_df)
  
  # --------------------------------------------------------------------
  # B. Train Calibration GAMs
  # --------------------------------------------------------------------
  if (mod_name == "YOLOv12") {
    gams <- fit_calibration_gams(
      det_df = mod_data$det,
      formula_gb = best_formula_yolo_gb,
      formula_mab = best_formula_yolo_mab
    )
  } else {
    gams <- fit_calibration_gams(
      det_df = mod_data$det,
      formula_gb = best_formula_cascade_gb,
      formula_mab = best_formula_cascade_mab
    )
  }
  
  safe_mod_name <- gsub(" |-", "", mod_name)
  saveRDS(gams, file = paste0("../data/processed/", safe_mod_name, "_crab_detgams_multi_2426.rds"))
  
  # --------------------------------------------------------------------
  # C. Test GAMs on Holdout Data (Multiclass Counts)
  # --------------------------------------------------------------------
  cat("\nTesting GAMs on holdout data across crab species...\n")
  eval_results <- test_calibration_gams(
    gams          = gams, 
    det_df        = mod_data$det, 
    img_df        = mod_data$img, 
    f1_thresh_df  = f1_thresh_df,
    target_classes = target_classes
  )
  
  saveRDS(eval_results, file = paste0("../data/processed/", safe_mod_name, "_crab_deteval_multi_2426.rds"))
  
  cat("\nSpecies & Total Performance Metrics (Holdout Set):\n")
  print(eval_results$metrics, n = Inf)
  
  # --------------------------------------------------------------------
  # D. Generate Multiclass Confusion Matrix
  # --------------------------------------------------------------------
  cat("\nGenerating multiclass confusion matrix...\n")
  # Select average model F1 threshold for confusion matrix display
  mean_thresh <- mean(f1_thresh_df$conf, na.rm = TRUE)
  
  cm_df <- generate_confusion_matrix(
    det_data       = eval_results$det_eval %>% filter(dataset == "test_GAM_test"),
    img_data       = eval_results$img_eval %>% filter(dataset == "test_GAM_test"),
    threshold      = mean_thresh,
    target_classes = target_classes
  )
  
  cat("Confusion Matrix Output (Threshold =", round(mean_thresh, 3), "):\n")
  print(cm_df %>% pivot_wider(names_from = predicted, values_from = n, values_fill = 0))
  
  # E. Generate and Save Diagnostics
  cat("\nGenerating diagnostic plots...\n")
  plots <- generate_gam_plots(
    eval_res = eval_results, 
    model_name = mod_name, 
    dataset_label = "2024, 2026",
    target_classes = target_classes
  )
  
  dir.create("../figures/diag_crab_multi", showWarnings = FALSE, recursive = TRUE)
  
  ggsave(paste0("../figures/diag_crab_multi/2426_", safe_mod_name, "_species_strat.png"), plot = plots$p_species_fit, width = 11, height = 7, bg = "white")
  ggsave(paste0("../figures/diag_crab_multi/2426_", safe_mod_name, "_count_strat.png"), plot = plots$p_zoomed, width = 8, height = 5, bg = "white")
  ggsave(paste0("../figures/diag_crab_multi/2426_", safe_mod_name, "_resid_strat.png"), plot = plots$p_resid, width = 8, height = 4.5, bg = "white")
  ggsave(paste0("../figures/diag_crab_multi/2426_", safe_mod_name, "_sum_strat.png"), plot = plots$p_sum, width = 9, height = 5.5, bg = "white")}

cat("\nPipeline complete! All multiclass outputs saved.\n")
