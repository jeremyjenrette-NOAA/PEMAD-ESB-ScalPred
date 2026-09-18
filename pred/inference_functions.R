library(dplyr)
library(tidyr)
library(mgcv)

# ======================================================================
# 1. Predict Detection-Level Calibration
# ======================================================================
predict_calibration_gams <- function(det_df, gams) {
  # Filter out detections with missing or unassigned regions
  det_df <- det_df %>%
    filter(!is.na(region) & region %in% names(gams))
  
  det_df$pred_p <- NA_real_
  
  # Predict safely across available regions
  for (reg in names(gams)) {
    reg_idx <- which(det_df$region == reg)
    if (length(reg_idx) > 0) {
      det_df$pred_p[reg_idx] <- predict(
        gams[[reg]], 
        newdata = det_df[reg_idx, ], 
        type = "response", 
        na.action = na.pass
      )
    }
  }
  
  return(det_df)
}

# ======================================================================
# 2. Aggregate to Image Level (with F1 Cutoff Support)
# ======================================================================
aggregate_image_predictions <- function(det_df_pred, img_df, f1_thresholds = NULL) {
  # Filter out images with missing regions
  img_df <- img_df %>%
    filter(!is.na(region) & region != "")
  
  # Assign regional F1 cutoff thresholds if provided
  if (!is.null(f1_thresholds)) {
    det_df_pred$regional_f1_thresh <- ifelse(
      det_df_pred$region == "GB", 
      f1_thresholds[["GB"]], 
      f1_thresholds[["MAB"]]
    )
  }
  
  # Sum probabilities and F1 thresholded detections per image
  img_predictions <- det_df_pred %>%
    group_by(image_id) %>%
    summarise(
      raw_detection_number = n(),
      predicted_number     = sum(pred_p, na.rm = TRUE),
      predicted_f1_number  = if (!is.null(f1_thresholds)) sum(conf >= regional_f1_thresh, na.rm = TRUE) else NA_integer_,
      .groups = "drop"
    )
  
  # Join back to master image list, filling 0 for empty images
  img_eval <- img_df %>%
    left_join(img_predictions, by = "image_id") %>%
    mutate(
      raw_detection_number = replace_na(raw_detection_number, 0),
      predicted_number     = replace_na(predicted_number, 0),
      predicted_f1_number  = replace_na(predicted_f1_number, 0)
    )
  
  return(img_eval)
}