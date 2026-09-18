library(dplyr)
library(mgcv)
library(tidyr)
library(ggplot2)
library(tools)
source("./gamfunc.R")

<<<<<<< HEAD
# 1. Load data & build multiclass synergy dataset
meta      <- read.csv("../data/raw/dataset_split_crab.csv") %>% janitor::clean_names() %>% rename(n_annotations = total_annotations) %>% mutate(image_id = stringr::str_remove(imagename, "\\.[A-Za-z0-9]+$")) %>% distinct(image_id, .keep_all = TRUE)
yolo_eval <- readRDS("../data/processed/YOLOv12_crab_deteval_multi_2426.rds")
cas_eval  <- readRDS("../data/processed/CascadeRCNN_crab_deteval_multi_2426.rds")

img_combined <- build_synergy_dataset_multi(yolo_eval, cas_eval, meta)

# 2. Define Systematic Candidate Formula Templates
candidate_templates <- list(
  # STAGE 1: Single-Model Log Baselines
  M01_YoloOnly     = "n_{sp} ~ s(log_yolo_{sp})",
  M02_CascadeOnly  = "n_{sp} ~ s(log_cascade_{sp})",
=======
# 1. Load metadata and past evaluation results
meta <- read.csv("../data/raw/dataset_split_2226.csv") %>% 
  janitor::clean_names() %>%
  mutate(image_id = stringr::str_remove(imagename, "\\.[A-Za-z0-9]+$")) %>%
  distinct(image_id, .keep_all = TRUE)

yolo_eval <- readRDS("../data/processed/YOLOv12_deteval_2226.rds")
cas_eval  <- readRDS("../data/processed/CascadeR-CNN_deteval_2226.rds")

# 2. Build Synergistic Dataset & ensure log_ratio is defined for M16_PropDiff
img_combined <- build_synergy_dataset(yolo_eval, cas_eval, meta) %>%
  mutate(log_ratio = log((pred_yolo + 1) / (pred_cascade + 1)))

# ======================================================================
# 2. Define Systematic Synergistic Candidate Formulas
# ======================================================================
image_candidate_forms <- list(
  # STAGE 1: Single Model Raw Baselines
  M01_YoloRaw     = n_annotations ~ s(pred_yolo),
  M02_CascadeRaw  = n_annotations ~ s(pred_cascade),
  M03_MeanRaw     = n_annotations ~ s(pred_mean),
  M035_SumRaw     = n_annotations ~ s(pred_sum),
>>>>>>> 099e1fe8751b6dd1e92581670bbe7f7234962132
  
  # STAGE 2: Additive & Tensor Synergy
  M03_Additive     = "n_{sp} ~ s(log_yolo_{sp}) + s(log_cascade_{sp})",
  M04_Tensor       = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp})",
  
<<<<<<< HEAD
  # STAGE 3: Adding Disagreement / Difference Covariate
  M05_TensorDiff   = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp}) + s(diff_{sp})",
  M06_TensorTi     = "n_{sp} ~ ti(log_yolo_{sp}, log_cascade_{sp}) + s(diff_{sp})",
  
  # STAGE 4: Environmental Interaction
  M07_EnvOptics    = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp}) + s(diff_{sp}) + s(altitude, k = 5)",
  M08_EnvFull      = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp}) + s(diff_{sp}) + s(altitude, k = 5) + s(bottom_depth, k = 5)"
)

# 3. Run Selection Phase
selection_gb  <- compare_synergy_gams_multi(img_combined, candidate_templates, region_focus = "GB")
selection_mab <- compare_synergy_gams_multi(img_combined, candidate_templates, region_focus = "MAB")

# Print top-ranked models per species
print(selection_gb  %>% select(species, Model_ID, AIC, Deviance_Explained, RMSE_Train), n = Inf)
print(selection_mab %>% select(species, Model_ID, AIC, Deviance_Explained, RMSE_Train), n = Inf)

# 4. Extract Winner Formulas automatically (or manually overwrite)
best_gb  <- extract_best_formulas(selection_gb, metric = "AIC")
best_mab <- extract_best_formulas(selection_mab, metric = "AIC")
=======
  # STAGE 3: Disagreement/Difference Covariate
  M07_YoloDiff    = n_annotations ~ s(log_yolo_pred, pred_diff),
  M08_CascadeDiff = n_annotations ~ s(log_cascade_pred, pred_diff),
  M09_MeanDiff    = n_annotations ~ s(log_pred_mean, pred_diff),
  M09_MeanCasDiff = n_annotations ~ s(log_pred_mean, pred_diff) + s(log_cascade_pred, pred_diff),
  M10_YoloLogDiff = n_annotations ~ s(log_yolo_pred, log_pred_diff),
  
  # STAGE 4: Multi-Model Synergy
  M11_Additive    = n_annotations ~ s(log_yolo_pred) + s(log_cascade_pred) + s(log_pred_diff),
  M12_Tensor      = n_annotations ~ te(log_yolo_pred, log_cascade_pred),
  M13_TensorDiff  = n_annotations ~ ti(log_yolo_pred, log_cascade_pred) + log_pred_diff,
  M14_Champion    = n_annotations ~ s(log_yolo_pred, log_pred_diff) + s(log_cascade_pred, log_pred_diff),
  # M15_TensorDiff2 = n_annotations ~ ti(log_yolo_pred, log_cascade_pred) + 
  #   ti(log_cascade_pred, log_pred_diff, k = c(6,6)) + 
  #   ti(log_yolo_pred, log_pred_diff, k = c(6,6)),
  
  # STAGE 5: Advanced ANOVA & Ratio Models
  M15b_FullDecomp = n_annotations ~ s(log_yolo_pred, k = 8) + 
    s(log_cascade_pred, k = 8) + 
    s(log_pred_diff, k = 8) + 
    ti(log_yolo_pred, log_cascade_pred, k = c(6, 6)) + 
    ti(log_cascade_pred, log_pred_diff, k = c(6, 6)) + 
    ti(log_yolo_pred, log_pred_diff, k = c(6, 6)),
  
  M16_PropDiff    = n_annotations ~ s(log_yolo_pred, k = 8) + 
    s(log_cascade_pred, k = 8) + 
    s(log_ratio, k = 8) + 
    ti(log_yolo_pred, log_cascade_pred, k = c(6, 6)) + 
    ti(log_yolo_pred, log_ratio, k = c(6, 6)) + 
    ti(log_cascade_pred, log_ratio, k = c(6, 6)),
  
  M17_Tensor3D    = n_annotations ~ te(log_yolo_pred, log_cascade_pred, log_pred_diff, k = c(5, 5, 5))
)

# ======================================================================
# Helper Function: Automated Formula Selection
# ======================================================================
select_best_synergy_model <- function(eval_table, candidate_list, metric = "AIC") {
  # Sort table by designated metric
  sorted_table <- eval_table %>% arrange(.data[[metric]])
  
  best_id   <- sorted_table$Model_ID[1]
  best_score <- round(sorted_table[[metric]][1], 2)
  best_r2    <- round(sorted_table$R2_Train[1], 3)
  
  cat(sprintf("  [Selected] %s (by %s = %g | R² = %g)\n", best_id, metric, best_score, best_r2))
  
  return(candidate_list[[best_id]])
}

# ======================================================================
# Run Evaluation and Dynamically Assign Best Models
# ======================================================================
cat("--- Synergy Selection: Georges Bank ---\n")
selection_gb <- compare_image_gams(img_combined, image_candidate_forms, "GB")
print(selection_gb, n = Inf)

best_syn_gb <- select_best_synergy_model(selection_gb, image_candidate_forms, metric = "AIC")

cat("\n--- Synergy Selection: Mid-Atlantic Bight ---\n")
selection_mab <- compare_image_gams(img_combined, image_candidate_forms, "MAB")
print(selection_mab, n = Inf)

best_syn_mab <- select_best_synergy_model(selection_mab, image_candidate_forms, metric = "AIC")

# best_syn_gb <- image_candidate_forms$M15b_FullDecomp
# best_syn_mab <- image_candidate_forms$M15b_FullDecomp
# ======================================================================
# 3. Train and Test Final Synergistic Model
# ======================================================================
syn_gams <- fit_synergy_gams(img_combined, best_syn_gb, best_syn_mab)
syn_eval <- test_synergy_gams(syn_gams, img_combined)

saveRDS(syn_gams, file = paste0("../data/processed/Syngams_2226.rds"))
# ======================================================================
# 4. The 12-Panel Plot Generation
# ======================================================================
>>>>>>> 099e1fe8751b6dd1e92581670bbe7f7234962132

# 5. Fit & Test Winning Synergistic GAMs on Holdout Data
syn_gams <- fit_synergy_gams_multi(img_combined, best_gb, best_mab)
syn_eval <- test_synergy_gams_multi(syn_gams, img_combined)

# Print final holdout performance across species and total
print(syn_eval$metrics)

# ======================================================================
# Macro Abundance Comparison Across All Workflows
# ======================================================================
generate_plot2_macro_abundance <- function(
    yolo_eval, cas_eval, syn_eval, 
    target_classes = c("jonah_crab", "rock_crab", "cancer_sp")
) {
  
  # 1. Align holdout test image IDs
  test_ids <- syn_eval$img_eval %>% filter(dataset == "test_GAM_test") %>% pull(image_id)
  
  y_test <- yolo_eval$img_eval %>% filter(image_id %in% test_ids)
  c_test <- cas_eval$img_eval %>% filter(image_id %in% test_ids)
  s_test <- syn_eval$img_eval %>% filter(image_id %in% test_ids)
  
  # 2. Extract totals per species across all 6 workflows
  macro_df <- purrr::map_dfr(target_classes, function(sp) {
    gt_cnt    <- sum(s_test[[paste0("n_", sp)]], na.rm = TRUE)
    y_f1_cnt  <- sum(y_test[[paste0("pred_f1_", sp)]], na.rm = TRUE)
    y_gam_cnt <- sum(y_test[[paste0("pred_gam_", sp)]], na.rm = TRUE)
    c_f1_cnt  <- sum(c_test[[paste0("pred_f1_", sp)]], na.rm = TRUE)
    c_gam_cnt <- sum(c_test[[paste0("pred_gam_", sp)]], na.rm = TRUE)
    syn_cnt   <- sum(s_test[[paste0("pred_synergy_", sp)]], na.rm = TRUE)
    
    tibble(
      species = gsub("_", " ", tools::toTitleCase(sp)),
      `Ground Truth`    = gt_cnt,
      `YOLO F1`         = y_f1_cnt,
      `YOLO Det-GAM`    = y_gam_cnt,
      `Cascade F1`      = c_f1_cnt,
      `Cascade Det-GAM` = c_gam_cnt,
      `Synergistic GAM` = syn_cnt
    )
  }) %>%
    pivot_longer(-species, names_to = "Workflow", values_to = "Count") %>%
    mutate(
      Workflow = factor(Workflow, levels = c(
        "Ground Truth", "YOLO F1", "YOLO Det-GAM", 
        "Cascade F1", "Cascade Det-GAM", "Synergistic GAM"
      ))
    )
  
  # 3. Calculate % error relative to Ground Truth
  gt_lookup <- macro_df %>% filter(Workflow == "Ground Truth") %>% select(species, GT = Count)
  
  macro_df <- macro_df %>%
    left_join(gt_lookup, by = "species") %>%
    mutate(
      Err_Pct = ifelse(Workflow == "Ground Truth", 0, ((Count - GT) / GT) * 100),
      Label = case_when(
        Workflow == "Ground Truth" ~ sprintf("%d", round(Count)),
        TRUE ~ sprintf("%d\n(%+.1f%%)", round(Count), Err_Pct)
      )
    )
  
  # 4. Generate Plot
  p2 <- ggplot(macro_df, aes(x = species, y = Count, fill = Workflow)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75, color = "black", linewidth = 0.3) +
    geom_text(
      aes(label = Label),
      position = position_dodge(width = 0.85),
      vjust = -0.25, size = 2.9, fontface = "bold"
    ) +
    scale_fill_manual(values = c(
      "Ground Truth"    = "grey40", 
      "YOLO F1"         = "#FB6A4A", 
      "YOLO Det-GAM"    = "#CB181D", 
      "Cascade F1"      = "#6BAED6", 
      "Cascade Det-GAM" = "#2171B5", 
      "Synergistic GAM" = "#7570B3"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    theme_minimal(base_size = 13) +
    labs(
      title = "Multi-Class Abundance Estimation across Workflows",
      subtitle = "Sum count and % error relative to Ground Truth evaluated on holdout test set",
      x = "Crab Species",
      y = "Total Abundance Count",
      fill = "Workflow"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
      axis.text.x = element_text(face = "bold", size = 11),
      panel.grid.major.x = element_blank()
    )
  
  return(p2)
}

# Generate and save Plot 2
p_macro_abun <- generate_plot2_macro_abundance(yolo_eval, cas_eval, syn_eval)
p_macro_abun
ggsave("../figures/diag_crab_multi/2426_workflow_abundance_summary.png", plot = p_macro_abun, width = 11, height = 6, dpi = 300, bg = "white")

# ======================================================================
# 2x3 Grid Regional & Species Synergistic Calibration
# ======================================================================
<<<<<<< HEAD
generate_plot1_synergy_calibration <- function(
    syn_eval, 
    target_classes = c("jonah_crab", "rock_crab", "cancer_sp")
) {
=======

# Create Overarching Headers using Grid graphics wrapped for patchwork
gb_header  <- wrap_elements(textGrob("Georges Bank", gp = gpar(fontsize = 16, fontface = "bold")))
mab_header <- wrap_elements(textGrob("Mid-Atlantic Bight", gp = gpar(fontsize = 16, fontface = "bold")))
header_row <- gb_header | mab_header

# Assemble the rows (4 columns wide)
row1 <- p[[1]] | p[[2]] | p[[3]] | p[[4]]
row2 <- p[[5]] | p[[6]] | p[[7]] | p[[8]]
row3 <- p[[9]] | p[[10]]| p[[11]]| p[[12]]

# Create dynamic subcaption string
caption_string <- sprintf("Datasets: 2022-2024, 2026 | Models Evaluated: YOLOv12, Cascade R-CNN, Synergistic GAM\nHoldout Test Images: Georges Bank (n = %d), Mid-Atlantic Bight (n = %d)\nScallop Abundance: Georges Bank (n = %d), Mid-Atlantic Bight (n = %d)", n_gb, n_mab, n_gb_abundance, n_mab_abundance)

# Stitch it all together: Header row on top, followed by the 3 data rows
final_plot <- (header_row / row1 / row2 / row3) + 
  plot_layout(heights = c(0.1, 1, 1, 1)) + # Gives the header a slim vertical footprint
  plot_annotation(
    title = bquote(bold("Comparative Assessment: Optimal F"[1] ~ "vs. Detection GAMs vs. Synergistic Calibration")),
    caption = caption_string,
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.caption = element_text(size = 12, face = "italic", color = "grey30", hjust = 0.5, margin = margin(t = 15))
    )
  )

print(final_plot)

# Save the high-resolution figure
ggsave(
  filename = "../figures/ms_figures/2226_12panel_calibration_summary.png",
  plot = final_plot,
  width = 13,
  height = 12,
  device = "png",
  bg = "white"
)
ggsave(
  filename = "~/saltnoaa/presentations/figures/2226_12panel_calibration_summary.png",
  plot = final_plot,
  width = 17,
  height = 12,
  device = "png",
  bg = "white"
)


###########################################################################

# ======================================================================
# 6. Final Comprehensive Abundance Summary Plot
# ======================================================================

# Ensure we are strictly using the holdout test set
true_val <- sum(s_test$n_annotations, na.rm = TRUE)

# Build the comparison dataframe
sum_df <- data.frame(
  Method = factor(
    c("True Abundance", 
      "YOLOv12 F1", "YOLOv12 Det-GAM", 
      "Cascade F1", "Cascade Det-GAM", 
      "Synergistic GAM"),
    levels = c("True Abundance", 
               "YOLOv12 F1", "YOLOv12 Det-GAM", 
               "Cascade F1", "Cascade Det-GAM", 
               "Synergistic GAM")
  ),
  Value = c(
    true_val,
    sum(y_test$predicted_f1_number, na.rm = TRUE),
    sum(y_test$predicted_number, na.rm = TRUE),
    sum(c_test$predicted_f1_number, na.rm = TRUE),
    sum(c_test$predicted_number, na.rm = TRUE),
    sum(s_test$pred_synergy, na.rm = TRUE)
  )
)

library(forcats) 
library(latex2exp)
library(ggplot2)
library(dplyr)

# 1. Define the Labels in the EXACT same order as the factor levels
label_map <- c(
  "True Abundance"  = "True~Abundance",
  "YOLOv12 F1"      = "YOLO~F[1]",
  "YOLOv12 Det-GAM" = "YOLO~Sigma*p[det]", # Moved to 3rd to match factor levels
  "Cascade F1"      = "Cascade~F[1]",       # Moved to 4th to match factor levels
  "Cascade Det-GAM" = "Cascade~Sigma*p[det]",
  "Synergistic GAM" = "Synergistic~GAM"
)

sum_df <- sum_df %>%
  mutate(
    Error_Pct = ((Value - true_val) / true_val) * 100,
    Label_Text = case_when(
      Method == "True Abundance" ~ as.character(round(Value, 0)),
      TRUE ~ sprintf("%d\n(%+.1f%%)", round(Value, 0), Error_Pct)
    )
  )

# Update global text geom defaults once before plotting
ggplot2::update_geom_defaults("text", list(family = "sans"))

# 2. Generate the plot with parsed expressions directly
p_sum_final <- ggplot(sum_df, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.7, color = "black") +
  geom_hline(yintercept = true_val, linetype = "dashed", color = "black", linewidth = 1.0) +
  geom_text(aes(label = Label_Text), y = 10000, 
            vjust = 0, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c(
    "True Abundance"  = "grey40", 
    "YOLOv12 F1"      = "#FB6A4A", 
    "Cascade F1"      = "#6BAED6", 
    "YOLOv12 Det-GAM" = "#CB181D", 
    "Cascade Det-GAM" = "#2171B5", 
    "Synergistic GAM" = "#7570B3"
  )) +
>>>>>>> 099e1fe8751b6dd1e92581670bbe7f7234962132
  
  # 1. Isolate holdout test set
  syn_test <- syn_eval$img_eval %>% filter(dataset == "test_GAM_test")
  
  # 2. Reshape cleanly using purrr (avoids pivot mismatch bugs)
  sp_long <- purrr::map_dfr(target_classes, function(sp) {
    gt_col   <- paste0("n_", sp)
    pred_col <- paste0("pred_synergy_", sp)
    
    syn_test %>%
      select(image_id, region, true_count = all_of(gt_col), pred_count = all_of(pred_col)) %>%
      mutate(species_raw = sp)
  })
  
  # 3. Explicitly format labels & factor levels
  sp_long <- sp_long %>%
    mutate(
      species_clean = factor(
        case_when(
          species_raw == "jonah_crab" ~ "Jonah Crab",
          species_raw == "rock_crab"  ~ "Rock Crab",
          species_raw == "cancer_sp"  ~ "Cancer Sp",
          TRUE ~ species_raw
        ),
        levels = c("Jonah Crab", "Rock Crab", "Cancer Sp")
      ),
      region_full = factor(
        case_when(
          region == "GB"  ~ "Georges Bank",
          region == "MAB" ~ "Mid-Atlantic Bight",
          TRUE ~ as.character(region)
        ),
        levels = c("Georges Bank", "Mid-Atlantic Bight")
      )
    )
  
  # 4. Format R² & RMSE metrics per panel (Region x Species)
  r2_annotations <- syn_eval$metrics %>%
    filter(species != "TOTAL_CRABS") %>%
    mutate(
      species_clean = factor(
        case_when(
          species == "jonah_crab" ~ "Jonah Crab",
          species == "rock_crab"  ~ "Rock Crab",
          species == "cancer_sp"  ~ "Cancer Sp",
          TRUE ~ species
        ),
        levels = c("Jonah Crab", "Rock Crab", "Cancer Sp")
      ),
      region_full = factor(
        case_when(
          region == "GB"  ~ "Georges Bank",
          region == "MAB" ~ "Mid-Atlantic Bight",
          TRUE ~ as.character(region)
        ),
        levels = c("Georges Bank", "Mid-Atlantic Bight")
      ),
      label_str = sprintf("R² = %.3f\nRMSE = %.2f", r2_synergy, rmse_synergy)
    )
  
  # 5. Explicit High-Contrast Palette
  species_palette <- c(
    "Jonah Crab" = "green",  # Vibrant Teal
    "Rock Crab"  = "blue",  # Vivid Orange
    "Cancer Sp"  = "red"   # Deep Purple
  )
<<<<<<< HEAD
  
  # 6. Generate 2x3 Grid Plot
  p1 <- ggplot(sp_long, aes(x = pred_count, y = true_count, color = species_clean, fill = species_clean)) +
    # 1:1 Reference Line of Perfect Agreement
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    
    # Species-colored Data Points
    geom_point(alpha = 0.7, size = 2.5) +
    
    # Species-colored Linear Regression Line with 95% Confidence Intervals
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1.1) +
    
    # 2 Rows (Regions) x 3 Columns (Species) Facet Grid
    facet_grid(region_full ~ species_clean, scales = "free") +
    
    # Per-Panel R² & RMSE Text Overlay
    geom_text(
      data = r2_annotations,
      aes(x = -Inf, y = Inf, label = label_str),
      hjust = -0.15, vjust = 1.25, inherit.aes = FALSE,
      size = 3.8, fontface = "bold", family = "sans", color = "black"
    ) +
    
    # Apply Manual Color Palette
    scale_color_manual(values = species_palette) +
    scale_fill_manual(values = species_palette) +
    
    # Clean Theme Styling
    theme_bw(base_size = 13) +
    labs(
      title = "Synergistic GAM Multi-Class Calibration",
      subtitle = "Evaluated on Holdout Test Set across Regions (Rows) and Species (Columns)",
      x = "Synergistic GAM Predicted Count (Σ p)",
      y = "Manual Ground Truth Count"
    ) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "grey50", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",  # Top column headers already identify species
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 12),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black")
    )
  
  return(p1)
}

# Re-generate and save Plot 1
p_synergy_reg <- generate_plot1_synergy_calibration(syn_eval)
p_synergy_reg
ggsave("../figures/diag_crab_multi/2426_synergy_regional_calibration.png", plot = p_synergy_reg, width = 11, height = 7, dpi = 300, bg = "white")
=======
p_sum_final
# Save the final render
ggsave("../figures/ms_figures/2226_final_abundance_summary.png",
       plot = p_sum_final, width = 8, height = 6, dpi = 300)
>>>>>>> 099e1fe8751b6dd1e92581670bbe7f7234962132
