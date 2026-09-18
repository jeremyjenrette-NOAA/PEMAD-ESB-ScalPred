library(dplyr)
library(mgcv)
library(tidyr)
library(ggplot2)
library(tools)
source("./gamfunc.R")

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
  
  # STAGE 2: Additive & Tensor Synergy
  M03_Additive     = "n_{sp} ~ s(log_yolo_{sp}) + s(log_cascade_{sp})",
  M04_Tensor       = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp})",
  
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
generate_plot1_synergy_calibration <- function(
    syn_eval, 
    target_classes = c("jonah_crab", "rock_crab", "cancer_sp")
) {
  
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
