# ======================================================================
# 1. New Function: Evaluate Unified Detection GAMs for Manuscript Table
# ======================================================================
compare_unified_detection_gams <- function(det_df, img_df, candidate_formulas) {
  
  # 1. Prepare unified train and test sets (both GB and MAB together)
  train_data <- det_df %>%
    filter(dataset == "test_GAM_train") %>%
    drop_na(truedetect, conf, bottom_depth, latitude, longitude, altitude, 
            field_of_view_sq_meter, fluorometer_backscatter_ntu, ctd_temperature_celsius,
            boxsize, iu, region)
  
  test_det <- det_df %>% filter(dataset == "test_GAM_test")
  test_img <- img_df %>% filter(dataset == "test_GAM_test")
  
  cat("\nEvaluating", length(candidate_formulas), "unified models (n =", nrow(train_data), ")...\n")
  
  results <- purrr::map_dfr(names(candidate_formulas), function(mod_name) {
    form <- candidate_formulas[[mod_name]]
    cat("  Fitting:", mod_name, "...\n")
    
    # Fit unified binomial GAM
    fit <- gam(form, family = binomial(), data = train_data, method = "REML")
    
    # Extract deviance explained
    dev_expl <- summary(fit)$dev.expl
    
    # Predict on holdout detection set
    test_det$pred_p <- predict(fit, newdata = test_det, type = "response", na.action = na.pass)
    
    # Aggregate predictions to the image level
    img_preds <- test_det %>%
      group_by(image_id) %>%
      summarise(predicted_number = sum(pred_p, na.rm = TRUE), .groups = "drop")
    
    eval_df <- test_img %>%
      left_join(img_preds, by = "image_id") %>%
      mutate(predicted_number = replace_na(predicted_number, 0))
    
    # Calculate Image-Level Test Metrics
    r2_test <- summary(lm(n_annotations ~ predicted_number, data = eval_df))$adj.r.squared
    rmse_test <- sqrt(mean((eval_df$n_annotations - eval_df$predicted_number)^2, na.rm = TRUE))
    
    tibble(
      Model = mod_name,
      AIC = AIC(fit),
      Deviance_Explained = dev_expl,
      Test_R2 = r2_test,
      Test_RMSE = rmse_test
    )
  })
  
  return(results %>% arrange(AIC))
}

# ======================================================================
# 2. Define Unified Candidate Detection-Level GAM Structures
# ======================================================================
unified_candidates <- list(
  M01_ConfOnly = truedetect ~ region + s(conf, k = 5),
  M02_Spatial  = truedetect ~ region + s(conf, k = 5) + s(latitude, longitude, k = 7),
  M03_Depth    = truedetect ~ region + s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5),
  M04_Optical  = truedetect ~ region + s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(field_of_view_sq_meter, k = 5),
  M05_Boxsize  = truedetect ~ region + s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(boxsize, k = 5),
  M06_Environ  = truedetect ~ region + s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(fluorometer_backscatter_ntu, k = 5),
  
  M07_teDepth   = truedetect ~ region + te(conf, bottom_depth, k = 5) + s(latitude, longitude, k = 7),
  M08_teAlt     = truedetect ~ region + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5),
  M09_teOptics  = truedetect ~ region + te(conf, altitude, k = 5) + te(conf, field_of_view_sq_meter, k = 5) + s(latitude, longitude, k = 7),
  M10_teBoxsize = truedetect ~ region + te(conf, boxsize, k = 5) + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7),
  
  M11_PrevBest  = truedetect ~ region + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7),
  M12_PB_Depth  = truedetect ~ region + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(bottom_depth, k = 5),
  M13_PB_Backsc = truedetect ~ region + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(fluorometer_backscatter_ntu, k = 5),
  M14_PB_Temp   = truedetect ~ region + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(ctd_temperature_celsius, k = 5),
  M15_All_te    = truedetect ~ region + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + te(conf, bottom_depth, k = 5)
)
# ======================================================================
# 3. Run the Evaluation for YOLOv12 (or Cascade)
# ======================================================================
cat("--- Generating Unified Manuscript Table Metrics ---\n")
unified_results <- compare_unified_detection_gams(
  det_df = models$YOLOv12$det, 
  img_df = models$YOLOv12$img,
  candidate_formulas = unified_candidates
)

print(unified_results)


unified_results_cascade <- compare_unified_detection_gams(
  det_df = models$`Cascade R-CNN`$det, 
  img_df = models$`Cascade R-CNN`$img,
  candidate_formulas = unified_candidates
)

print(unified_results_cascade)

########################################################################

# ======================================================================
# 1. New Function: Evaluate Unified Synergistic GAMs for Manuscript Table
# ======================================================================
compare_unified_image_gams <- function(img_combined, candidate_formulas) {
  
  # Filter train and test sets (Both GB and MAB together)
  train_data <- img_combined %>% filter(dataset == "test_GAM_train") %>% drop_na(region)
  test_data  <- img_combined %>% filter(dataset == "test_GAM_test")
  
  cat("\nEvaluating", length(candidate_formulas), "unified synergistic models (n =", nrow(train_data), ")...\n")
  
  results <- purrr::map_dfr(names(candidate_formulas), function(mod_name) {
    form <- candidate_formulas[[mod_name]]
    cat("  Fitting:", mod_name, "...\n")
    
    # Fit unified Negative Binomial GAM
    fit <- gam(form, family = nb(), data = train_data, method = "REML")
    
    # Extract deviance explained
    dev_expl <- summary(fit)$dev.expl
    
    # Predict on holdout test set (type = "response" gives estimated counts)
    preds <- predict(fit, newdata = test_data, type = "response", na.action = na.pass)
    
    # Calculate Image-Level Test Metrics
    r2_test <- summary(lm(test_data$n_annotations ~ preds))$adj.r.squared
    rmse_test <- sqrt(mean((test_data$n_annotations - preds)^2, na.rm = TRUE))
    
    tibble(
      Model = mod_name,
      AIC = AIC(fit),
      Deviance_Explained = dev_expl,
      Test_R2 = r2_test,
      Test_RMSE = rmse_test
    )
  })
  
  return(results %>% arrange(AIC))
}

# ======================================================================
# 2. Define Systematic Synergistic Candidate Formulas
# ======================================================================
unified_image_candidates <- list(
  # STAGE 1: Single Model Raw Baselines
  M01_YoloRaw     = n_annotations ~ region + s(pred_yolo),
  M02_CascadeRaw  = n_annotations ~ region + s(pred_cascade),
  M03_MeanRaw     = n_annotations ~ region + s(pred_mean),
  
  # STAGE 2: Log-Transformed Single Predictors
  M04_YoloLog     = n_annotations ~ region + s(log_yolo_pred),
  M05_CascadeLog  = n_annotations ~ region + s(log_cascade_pred),
  M06_MeanLog     = n_annotations ~ region + s(log_pred_mean),
  
  # STAGE 3: Adding the Disagreement/Difference Covariate
  M07_YoloDiff    = n_annotations ~ region + s(log_yolo_pred, pred_diff),
  M08_CascadeDiff = n_annotations ~ region + s(log_cascade_pred, pred_diff),
  M09_MeanDiff    = n_annotations ~ region + s(log_pred_mean, pred_diff),
  M10_YoloLogDiff = n_annotations ~ region + s(log_yolo_pred, log_pred_diff),
  
  # STAGE 4: Multi-Model Synergy & Your Champion
  M11_Additive    = n_annotations ~ region + s(log_yolo_pred) + s(log_cascade_pred) + s(log_pred_diff),
  M12_Tensor      = n_annotations ~ region + te(log_yolo_pred, log_cascade_pred),
  M13_TensorDiff  = n_annotations ~ region + ti(log_yolo_pred, log_cascade_pred) + log_pred_diff,
  M14_Champion    = n_annotations ~ region + s(log_yolo_pred, log_pred_diff) + s(log_cascade_pred, log_pred_diff)
)

# ======================================================================
# 3. Run the Evaluation
# ======================================================================
cat("--- Generating Unified Synergistic Manuscript Table Metrics ---\n")
unified_syn_results <- compare_unified_image_gams(
  img_combined = img_combined, 
  candidate_formulas = unified_image_candidates
)

print(unified_syn_results)
