library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# ======================================================================
# 1. Core Precision-Recall Evaluation
# ======================================================================
evaluate_pr_curve <- function(
    calib_df,
    img_df,
    conf_grid = seq(0, 1, by = 0.0005), # Adjusted step slightly for speed, but can be 0.0005
    model_id  = "model"
) {
  
  map_dfr(conf_grid, function(th) {
    
    df_th <- calib_df %>% filter(conf >= th)
    
    TP <- sum(df_th$truedetect, na.rm = TRUE)
    FP <- sum(!df_th$truedetect, na.rm = TRUE)
    
    auto_true_img <- df_th %>%
      filter(truedetect) %>%
      count(image_id, name = "n_auto_true")
    
    # Calculate False Negatives using the new 'n_annotations' column
    FN <- img_df %>%
      left_join(auto_true_img, by = "image_id") %>%
      mutate(
        n_auto_true = replace_na(n_auto_true, 0L),
        FN = pmax(n_annotations - n_auto_true, 0L)
      ) %>%
      summarise(FN = sum(FN, na.rm = TRUE)) %>%
      pull(FN)
    
    precision <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
    recall    <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
    f1        <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_
    else 2 * precision * recall / (precision + recall)
    
    tibble(
      model     = model_id,
      conf      = th,
      TP        = TP,
      FP        = FP,
      FN        = FN,
      precision = precision,
      recall    = recall,
      f1        = f1
    )
  })
}

# ======================================================================
# 2. Average Precision (mAP) Calculator via Trapezoidal AUC
# ======================================================================
calculate_ap <- function(recall, precision) {
  valid <- !is.na(recall) & !is.na(precision)
  r <- recall[valid]
  p <- precision[valid]
  
  if (length(r) < 2) return(NA_real_)
  
  # Order by recall ascending to properly calculate AUC
  ord <- order(r)
  r <- r[ord]
  p <- p[ord]
  
  # Area under curve using Trapezoidal rule
  sum(diff(r) * (p[-1] + p[-length(p)]) / 2)
}

# ======================================================================
# 3. Model Wrapper (Handles unified .rds list and Regional Toggling)
# ======================================================================
evaluate_pr_models <- function(
    models_list, 
    conf_grid = seq(0, 1, by = 0.0005),
    stratify_region = TRUE
) {
  
  # Loop over the dynamically named models (e.g., YOLO, Cascade)
  pr_all <- map_dfr(names(models_list), function(mod_name) {
    
    img_data <- models_list[[mod_name]]$img
    det_data <- models_list[[mod_name]]$det
    
    if (stratify_region) {
      # Georges Bank subset
      gb_img <- img_data %>% filter(region == "GB")
      gb_det <- det_data %>% filter(region == "GB")
      pr_gb <- evaluate_pr_curve(gb_det, gb_img, conf_grid, paste0(mod_name, " GB"))
      
      # Mid-Atlantic Bight subset
      mab_img <- img_data %>% filter(region == "MAB")
      mab_det <- det_data %>% filter(region == "MAB")
      pr_mab <- evaluate_pr_curve(mab_det, mab_img, conf_grid, paste0(mod_name, " MAB"))
      
      return(bind_rows(pr_gb, pr_mab))
      
    } else {
      # Pooled across all regions
      return(evaluate_pr_curve(det_data, img_data, conf_grid, mod_name))
    }
  })
  
  # Generate standard PR Plot
  legend_title <- if (stratify_region) "Model by region" else "Model"
  p_pr <- ggplot(pr_all, aes(x = recall, y = precision, color = model)) +
    geom_path(linewidth = 1.2, na.rm = TRUE) +
    theme_minimal() +
    labs(
      title = "Precision–Recall",
      x = "Recall",
      y = "Precision",
      color = legend_title
    )
  
  # Calculate mAP for each evaluated model
  map_summary <- pr_all %>%
    group_by(model) %>%
    summarize(Average_Precision = calculate_ap(recall, precision), .groups = "drop")
  
  list(pr_all = pr_all, p_pr = p_pr, map_summary = map_summary)
}