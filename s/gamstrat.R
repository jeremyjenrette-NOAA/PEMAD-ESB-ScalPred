library(dplyr)
library(ggplot2)
source("./gamfunc.R")
source("./procfunc.R")

# ======================================================================
# 1. Load Data
# ======================================================================
models <- readRDS("../data/processed/eval_2226.rds")

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
  M15_All_te    = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + te(conf, bottom_depth, k = 5),
  
  # STAGE 4: 3D Multi-Way Physical & Optical Scale Tensors
  # Jointly models pixel boxsize, camera height/FOV, and detector confidence simultaneously
  M16_3D_OpticsBox = truedetect ~ te(conf, altitude, boxsize, k = c(5, 5, 5)) + 
    s(latitude, longitude, k = 7) + 
    s(bottom_depth, k = 5),
  M17_3D_FOVBox    = truedetect ~ te(conf, field_of_view_sq_meter, boxsize, k = c(5, 5, 5)) + 
    s(latitude, longitude, k = 7) + 
    s(bottom_depth, k = 5),
  
  # STAGE 5: Combined Environmental Synergy
  # Evaluates water-column clarity (backscatter) and temperature together rather than separately
  M18_FullEnviron  = truedetect ~ te(conf, altitude, k = 5) + 
    te(conf, boxsize, k = 7) + 
    s(latitude, longitude, k = 7) + 
    s(bottom_depth, k = 5) + 
    s(ctd_temperature_celsius, k = 5) + 
    s(fluorometer_backscatter_ntu, k = 5),
  
  # STAGE 6: Tensor Interaction ANOVA Decomposition (`ti`)
  # Explicitly separates univariate smooths from pure interaction terms to prevent main-effect bias
  M19_ANOVA_Optics = truedetect ~ s(conf, k = 7) + 
    s(altitude, k = 5) + 
    s(boxsize, k = 5) + 
    ti(conf, altitude, k = c(5, 5)) + 
    ti(conf, boxsize, k = c(5, 5)) + 
    s(latitude, longitude, k = 7) + 
    s(bottom_depth, k = 5),
  
  # STAGE 7: Fine-Scale Spatial Smooths & Overlap Quality (`iu`)
  # Leverages high-knot spatial splines for localized bed density + bounding box overlap quality
  M20_HighSpatial  = truedetect ~ te(conf, altitude, k = 5) + 
    te(conf, boxsize, k = 7) + 
    s(latitude, longitude, k = 15) + 
    s(bottom_depth, k = 5) + 
    s(ctd_temperature_celsius, k = 5)
)

# Append to candidate_forms in your script
candidate_forms_v2 <- c(candidate_forms, list(
  
  # STAGE 4: 3D Multi-Way Tensor Interactions
  # Captures tri-variate interaction between Confidence, Camera Altitude, and Bounding Box Size
  M16_3D_Optics = truedetect ~ te(conf, altitude, boxsize, k = c(5, 5, 5)) + 
    te(conf, field_of_view_sq_meter, k = 5) + 
    s(latitude, longitude, k = 7) + 
    s(bottom_depth, k = 5),
  
  # STAGE 5: Physical Scale Interaction + Environmental Drivers
  # Interacts Confidence with optical FOV and Box Size simultaneously
  M17_Physical_Scale = truedetect ~ te(conf, boxsize, field_of_view_sq_meter, k = c(5, 5, 5)) + 
    te(conf, altitude, k = 5) + 
    s(latitude, longitude, k = 7) + 
    s(ctd_temperature_celsius, k = 5),
  
  # Higher resolution spatial smoothing to capture localized habitat density effects
  M18_HighK_Spatial = truedetect ~ te(conf, altitude, k = 5) + 
    te(conf, field_of_view_sq_meter, k = 7) + 
    te(conf, boxsize, k = 7) + 
    s(latitude, longitude, k = 15) + 
    s(ctd_temperature_celsius, k = 5),
  
  # Full Environmental & Optical Interaction
  M19_Full_OpticEnv = truedetect ~ te(conf, altitude, k = 5) + 
    te(conf, field_of_view_sq_meter, k = 7) + 
    te(conf, boxsize, k = 7) + 
    te(conf, bottom_depth, k = 5) + 
    s(ctd_temperature_celsius, k = 5) + 
    s(fluorometer_backscatter_ntu, k = 5) + 
    s(latitude, longitude, k = 7)
))
# ======================================================================
# 3. Model Selection (Example: Testing structures on YOLO for Georges Bank)
# ======================================================================
cat("--- Running Model Selection Phase ---\n")
# Extract best formulas programmatically by lowest AIC
get_best_formula <- function(selection_results, candidate_list) {
  top_model_id <- selection_results$Model_ID[1]
  cat("  -> Dynamic Best Selection:", top_model_id, "(AIC =", round(selection_results$AIC[1], 1), ")\n")
  return(candidate_list[[top_model_id]])
}

# Dynamic Assignment
best_formula_yolo_gb     <- get_best_formula(selection_results_yolo_GB, candidate_forms_v2)
best_formula_yolo_mab    <- get_best_formula(selection_results_yolo_MAB, candidate_forms_v2)
best_formula_cascade_gb  <- get_best_formula(selection_results_cascade_GB, candidate_forms_v2)
best_formula_cascade_mab <- get_best_formula(selection_results_cascade_MAB, candidate_forms_v2)

# print(selection_results_yolo_GB, n = Inf)
# print(selection_results_yolo_MAB, n = Inf)
# print(selection_results_cascade_GB, n = Inf)
# print(selection_results_cascade_MAB, n = Inf)

# Pick the best formulas based on the selection results! 
# (You could automate this, but it is safer to manually inspect the AIC table and define them here)
# best_formula_yolo_gb  <- candidate_forms$M19_ANOVA_Optics
# best_formula_yolo_mab <- candidate_forms$M19_ANOVA_Optics
# best_formula_cascade_gb  <- candidate_forms$M19_ANOVA_Optics
# best_formula_cascade_mab <- candidate_forms$M19_ANOVA_Optics

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
  cat("Evaluating PR to find optimal F1 thresholds...\n")
  temp_model_list <- list()
  temp_model_list[[mod_name]] <- mod_data
  
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
  
  # B. Train Calibration GAMs
  gams <- fit_calibration_gams(
    det_df = mod_data$det,
    formula_gb = best_formula_yolo_gb,
    formula_mab = best_formula_yolo_mab
  )
  
  } else{
    
  gams <- fit_calibration_gams(
    det_df = mod_data$det,
    formula_gb = best_formula_cascade_gb,
    formula_mab = best_formula_cascade_mab
  )
    
  }
  
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
