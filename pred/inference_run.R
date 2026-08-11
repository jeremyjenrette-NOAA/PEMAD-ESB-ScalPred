library(dplyr)
library(mgcv)

# 1. Source the function repository
source("inference_functions.R")

# 2. Load the inference dataset and pre-trained GAMs
models    <- readRDS("../data/processed/eval_inference_2226.rds")
yolo_gams <- readRDS("../data/processed/YOLOv12_detgams_2226.rds")
cas_gams  <- readRDS("../data/processed/CascadeR-CNN_detgams_2226.rds")
syn_gams  <- readRDS("../data/processed/Syngams_2226.rds")

# Define region-specific optimal F1 confidence cutoff vectors
yolo_f1_thresh <- c(GB = 0.454, MAB = 0.390)
cas_f1_thresh  <- c(GB = 0.356, MAB = 0.517)

# ======================================================================
# 3. Process YOLOv12
# ======================================================================
cat("Processing YOLOv12 Inference...\n")
yolo_det_calib <- predict_calibration_gams(models$YOLOv12$det, yolo_gams)
yolo_img_calib <- aggregate_image_predictions(yolo_det_calib, models$YOLOv12$img, yolo_f1_thresh)

models$YOLOv12$det <- yolo_det_calib
models$YOLOv12$img <- yolo_img_calib

# ======================================================================
# 4. Process Cascade R-CNN
# ======================================================================
cat("Processing Cascade R-CNN Inference...\n")
cas_det_calib <- predict_calibration_gams(models$`Cascade R-CNN`$det, cas_gams)
cas_img_calib <- aggregate_image_predictions(cas_det_calib, models$`Cascade R-CNN`$img, cas_f1_thresh)

models$`Cascade R-CNN`$det <- cas_det_calib
models$`Cascade R-CNN`$img <- cas_img_calib

# ======================================================================
# 5. Process Region-Stratified Synergistic GAM Fusion
# ======================================================================
cat("Processing Synergistic GAM Fusion...\n")

syn_img_calib <- yolo_img_calib %>%
  rename(
    yolo_pred    = predicted_number,
    yolo_f1_pred = predicted_f1_number
  ) %>%
  inner_join(
    cas_img_calib %>% select(
      image_id, 
      cascade_pred    = predicted_number,
      cascade_f1_pred = predicted_f1_number
    ),
    by = "image_id"
  ) %>%
  filter(!is.na(region) & region %in% names(syn_gams)) %>%
  mutate(
    # Full covariate set for candidate formula swapping
    pred_mean        = (yolo_pred + cascade_pred) / 2,
    pred_diff        = cascade_pred - yolo_pred,
    pred_sum         = yolo_pred + cascade_pred,
    log_yolo_pred    = log1p(yolo_pred),
    log_cascade_pred = log1p(cascade_pred),
    log_pred_mean    = log1p(pred_mean),
    log_pred_sum     = log1p(pred_sum),
    log_pred_diff    = log1p(abs(pred_diff))
  ) %>%
  group_by(region) %>%
  group_modify(~ {
    # Force factor level to character string to ensure exact list key matching
    reg_str <- as.character(.y$region[[1]])
    
    if (reg_str %in% names(syn_gams)) {
      .x$final_abundance <- as.numeric(predict(syn_gams[[reg_str]], newdata = .x, type = "response"))
    } else {
      warning(paste("Region", reg_str, "not found in syn_gams. Setting final_abundance to NA."))
      .x$final_abundance <- NA_real_
    }
    return(.x)
  }) %>%
  ungroup() %>%
  mutate(
    # Zero-detection face-value override
    final_abundance = if_else(yolo_pred == 0 & cascade_pred == 0, 0, final_abundance)
  ) %>%
  filter(
    altitude < 3,
    bottom_depth > 20
  )

# Update master list structure
models$synergistic <- syn_img_calib

# ======================================================================
# 6. Export Calibrated Dataset
# ======================================================================
saveRDS(models, file = "../data/processed/eval_inference_calibrated_2226.rds")
cat("Inference calibration complete and saved!\n")