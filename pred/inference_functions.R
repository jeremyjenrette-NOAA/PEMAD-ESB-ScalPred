library(dplyr)
library(tidyr)
library(mgcv)

# ======================================================================
# 1. Predict Detection-Level Calibration
# ======================================================================
predict_calibration_gams <- function(det_df, gams) {
  # Initialize the column to safely hold probabilities
  det_df$pred_p <- NA_real_
  
  # Predict strictly for Georges Bank using the GB GAM
  gb_idx <- which(det_df$region == "GB")
  if (length(gb_idx) > 0) {
    det_df$pred_p[gb_idx] <- predict(gams$GB, newdata = det_df[gb_idx, ], 
                                     type = "response", na.action = na.pass)
  }
  
  # Predict strictly for Mid-Atlantic Bight using the MAB GAM
  mab_idx <- which(det_df$region == "MAB")
  if (length(mab_idx) > 0) {
    det_df$pred_p[mab_idx] <- predict(gams$MAB, newdata = det_df[mab_idx, ], 
                                      type = "response", na.action = na.pass)
  }
  
  return(det_df)
}

# ======================================================================
# 2. Aggregate to Image Level
# ======================================================================
aggregate_image_predictions <- function(det_df_pred, img_df) {
  # Sum the probabilities for images that have detections
  img_predictions <- det_df_pred %>%
    group_by(image_id) %>%
    summarise(
      raw_detection_number = n(),
      predicted_number     = sum(pred_p, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Join back to the master image list, filling 0 for empty images
  img_eval <- img_df %>%
    left_join(img_predictions, by = "image_id") %>%
    mutate(
      raw_detection_number = replace_na(raw_detection_number, 0),
      predicted_number     = replace_na(predicted_number, 0)
    )
  
  return(img_eval)
}