# ============================================================#
# gammod.R — Generalized Multi-Class GAM Evaluation
# ============================================================#
library(dplyr)
library(mgcv)
library(tidyr)
library(ggplot2)
library(tools)
source("./gamfunc.R")

# ------------------------------------------------------------#
# 1. Load Data & Build Synergy Dataset
# ------------------------------------------------------------#
models    <- readRDS("../data/processed/star_eval_24.rds")
meta      <- read.csv("../data/raw/dataset_split_star.csv") %>% 
  janitor::clean_names() %>% 
  rename(n_annotations = total_annotations) %>% 
  mutate(image_id = stringr::str_remove(imagename, "\\.[A-Za-z0-9]+$")) %>%
         # region = if_else(longitude >= -71, "GB", "MAB") %>% factor(levels = c("MAB", "GB"))) %>% 
  distinct(image_id, .keep_all = TRUE)

yolo_eval <- readRDS("../data/processed/YOLOv12_star_deteval_24.rds")
cas_eval  <- readRDS("../data/processed/CascadeRCNN_star_deteval_24.rds")

# Autolocate target classes across detection labels and metadata
target_sp <- get_target_classes(img_data = meta, det_data = yolo_eval$det)

# Build synergy dataset with dynamic species log counts and differences
img_combined <- build_synergy_dataset_multi(
  yolo_eval      = yolo_eval, 
  cas_eval       = cas_eval, 
  meta           = meta, 
  target_classes = target_sp
)

# ------------------------------------------------------------#
# 2. Define Systematic Candidate Formula Templates (with k bounds)
# ------------------------------------------------------------#
candidate_templates <- list(
  # ====================================================================
  # STAGE 1: Single-Model Baselines (Raw vs. Log)
  # ====================================================================
  M01a_Yolo_Raw        = "n_{sp} ~ s(yolo_{sp}, k = 4)",
  M01b_Yolo_Log        = "n_{sp} ~ s(log_yolo_{sp}, k = 4)",
  M02a_Cascade_Raw     = "n_{sp} ~ s(cascade_{sp}, k = 4)",
  M02b_Cascade_Log     = "n_{sp} ~ s(log_cascade_{sp}, k = 4)",
  
  # ====================================================================
  # STAGE 2: Additive & Tensor Synergy (Raw vs. Log)
  # ====================================================================
  M03a_Additive_Raw    = "n_{sp} ~ s(yolo_{sp}, k = 4) + s(cascade_{sp}, k = 4)",
  M03b_Additive_Log    = "n_{sp} ~ s(log_yolo_{sp}, k = 4) + s(log_cascade_{sp}, k = 4)",
  
  M04a_Tensor_Raw      = "n_{sp} ~ te(yolo_{sp}, cascade_{sp}, k = c(3, 3))",
  M04b_Tensor_Log      = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp}, k = c(3, 3))",
  
  # ====================================================================
  # STAGE 3: Disagreement / Difference Covariate (Raw vs. Log)
  # ====================================================================
  M05a_TensorDiff_Raw  = "n_{sp} ~ te(yolo_{sp}, cascade_{sp}, k = c(3, 3)) + s(diff_{sp}, k = 3)",
  M05b_TensorDiff_Log  = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp}, k = c(3, 3)) + s(diff_{sp}, k = 3)",
  
  M06a_TensorTi_Raw    = "n_{sp} ~ ti(yolo_{sp}, cascade_{sp}, k = c(3, 3)) + s(diff_{sp}, k = 3)",
  M06b_TensorTi_Log    = "n_{sp} ~ ti(log_yolo_{sp}, log_cascade_{sp}, k = c(3, 3)) + s(diff_{sp}, k = 3)",
  
  # ====================================================================
  # STAGE 4: Environmental Interaction (Raw vs. Log)
  # ====================================================================
  M07a_EnvOptics_Raw   = "n_{sp} ~ te(yolo_{sp}, cascade_{sp}, k = c(3, 3)) + s(diff_{sp}, k = 3) + s(altitude, k = 4)",
  M07b_EnvOptics_Log   = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp}, k = c(3, 3)) + s(diff_{sp}, k = 3) + s(altitude, k = 4)",
  
  M08a_EnvFull_Raw     = "n_{sp} ~ te(yolo_{sp}, cascade_{sp}, k = c(3, 3)) + s(diff_{sp}, k = 3) + s(altitude, k = 4) + s(bottom_depth, k = 4)",
  M08b_EnvFull_Log     = "n_{sp} ~ te(log_yolo_{sp}, log_cascade_{sp}, k = c(3, 3)) + s(diff_{sp}, k = 3) + s(altitude, k = 4) + s(bottom_depth, k = 4)"
)

# ------------------------------------------------------------#
# 3. Model Selection Phase Across Regions
# ------------------------------------------------------------#
regions_present <- if ("region" %in% names(img_combined)) unique(img_combined$region) else "ALL"
# regions_present <- "ALL"
selection_results <- purrr::map(regions_present, function(reg) {
  compare_synergy_gams_multi(
    img_combined        = img_combined, 
    candidate_templates = candidate_templates, 
    region_focus        = reg, 
    target_classes      = target_sp
  )
}) %>% setNames(regions_present)

# Print top-ranked models per region
purrr::walk(names(selection_results), function(reg) {
  cat(sprintf("\n=== Top Models for Region: %s ===\n", reg))
  print(selection_results[[reg]] %>% select(species, Model_ID, AIC, Deviance_Explained, RMSE_Train), n = Inf)
})

# ------------------------------------------------------------#
# 4. Extract Winning Formulas and Fit Holdout GAMs
# ------------------------------------------------------------#
best_formulas_gb  <- extract_best_formulas(selection_results[[1]], metric = "AIC")
best_formulas_mab <- if (length(selection_results) > 1) extract_best_formulas(selection_results[[2]], metric = "AIC") else NULL
# M05a_TensorDiff_Raw
syn_gams <- fit_synergy_gams_multi(img_combined, best_formulas_gb, best_formulas_mab, target_classes = target_sp)
syn_eval <- test_synergy_gams_multi(syn_gams, img_combined, target_classes = target_sp)

print("=== Holdout Metrics Summary ===")
print(syn_eval$metrics)

# ============================================================#
# Plotting Function 1: Macro Abundance Comparison
# ============================================================#
generate_plot2_macro_abundance <- function(
    yolo_eval, 
    cas_eval, 
    syn_eval, 
    target_classes = NULL
) {
  if (is.null(target_classes)) {
    target_classes <- get_target_classes(img_data = syn_eval$img_eval)
  }
  
  test_ids <- syn_eval$img_eval %>% filter(dataset == "test_GAM_test") %>% pull(image_id)
  
  y_test <- yolo_eval$img_eval %>% filter(image_id %in% test_ids)
  c_test <- cas_eval$img_eval  %>% filter(image_id %in% test_ids)
  s_test <- syn_eval$img_eval  %>% filter(image_id %in% test_ids)
  
  macro_df <- purrr::map_dfr(target_classes, function(sp) {
    gt_cnt    <- sum(s_test[[paste0("n_", sp)]], na.rm = TRUE)
    y_f1_cnt  <- if (paste0("pred_f1_", sp) %in% names(y_test)) sum(y_test[[paste0("pred_f1_", sp)]], na.rm = TRUE) else 0
    y_gam_cnt <- if (paste0("pred_gam_", sp) %in% names(y_test)) sum(y_test[[paste0("pred_gam_", sp)]], na.rm = TRUE) else 0
    c_f1_cnt  <- if (paste0("pred_f1_", sp) %in% names(c_test)) sum(c_test[[paste0("pred_f1_", sp)]], na.rm = TRUE) else 0
    c_gam_cnt <- if (paste0("pred_gam_", sp) %in% names(c_test)) sum(c_test[[paste0("pred_gam_", sp)]], na.rm = TRUE) else 0
    syn_cnt   <- if (paste0("pred_synergy_", sp) %in% names(s_test)) sum(s_test[[paste0("pred_synergy_", sp)]], na.rm = TRUE) else 0
    
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
  
  gt_lookup <- macro_df %>% filter(Workflow == "Ground Truth") %>% select(species, GT = Count)
  
  macro_df <- macro_df %>%
    left_join(gt_lookup, by = "species") %>%
    mutate(
      Err_Pct = ifelse(Workflow == "Ground Truth" | GT == 0, 0, ((Count - GT) / GT) * 100),
      Label = case_when(
        Workflow == "Ground Truth" ~ sprintf("%d", round(Count)),
        TRUE ~ sprintf("%d\n(%+.1f%%)", round(Count), Err_Pct)
      )
    )
  
  ggplot(macro_df, aes(x = species, y = Count, fill = Workflow)) +
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
      x = "Class",
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
}

# Render & Save Abundance Comparison
p_macro_abun <- generate_plot2_macro_abundance(yolo_eval, cas_eval, syn_eval, target_classes = target_sp)
p_macro_abun
ggsave("../figures/diag_star1/24_workflow_abundance_summary.png", plot = p_macro_abun, width = 11, height = 6, dpi = 300, bg = "white")

# ============================================================#
# Plotting Function 2: Synergistic Calibration Grid
# ============================================================#
generate_plot1_synergy_calibration <- function(
    syn_eval, 
    target_classes = NULL
) {
  if (is.null(target_classes)) {
    target_classes <- get_target_classes(img_data = syn_eval$img_eval)
  }
  
  syn_test <- syn_eval$img_eval %>% filter(dataset == "test_GAM_test")
  has_region <- "region" %in% names(syn_test)
  
  sp_long <- purrr::map_dfr(target_classes, function(sp) {
    gt_col   <- paste0("n_", sp)
    pred_col <- paste0("pred_synergy_", sp)
    
    df <- syn_test %>%
      select(image_id, true_count = all_of(gt_col), pred_count = all_of(pred_col), any_of("region")) %>%
      mutate(species_raw = sp)
    
    if (!"region" %in% names(df)) df$region <- "ALL"
    return(df)
  })
  
  sp_long <- sp_long %>%
    mutate(
      species_clean = gsub("_", " ", tools::toTitleCase(species_raw)),
      region_full = case_when(
        region == "GB"  ~ "Georges Bank",
        region == "MAB" ~ "Mid-Atlantic Bight",
        TRUE ~ as.character(region)
      )
    )
  
  r2_annotations <- syn_eval$metrics %>%
    filter(!grepl("^TOTAL", species)) %>%
    mutate(
      species_clean = gsub("_", " ", tools::toTitleCase(species)),
      region_full = case_when(
        region == "GB"  ~ "Georges Bank",
        region == "MAB" ~ "Mid-Atlantic Bight",
        TRUE ~ as.character(region)
      ),
      label_str = sprintf("R² = %.3f\nRMSE = %.2f", r2_synergy, rmse_synergy)
    )
  
  p1 <- ggplot(sp_long, aes(x = pred_count, y = true_count, color = species_clean, fill = species_clean)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1.1) +
    geom_text(
      data = r2_annotations,
      aes(x = -Inf, y = Inf, label = label_str),
      hjust = -0.15, vjust = 1.25, inherit.aes = FALSE,
      size = 3.8, fontface = "bold", family = "sans", color = "black"
    ) +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    theme_bw(base_size = 13) +
    labs(
      title = "Synergistic GAM Multi-Class Calibration",
      subtitle = "Evaluated on Holdout Test Set",
      x = "Synergistic GAM Predicted Count (Σ p)",
      y = "Manual Ground Truth Count"
    ) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "grey50", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 12),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black")
    )
  
  # Adaptive Faceting
  if (has_region && length(unique(sp_long$region)) > 1) {
    p1 <- p1 + facet_grid(region_full ~ species_clean, scales = "free")
  } else {
    p1 <- p1 + facet_wrap(~ species_clean, scales = "free")
  }
  
  return(p1)
}

# Render & Save Regional Calibration
p_synergy_reg <- generate_plot1_synergy_calibration(syn_eval, target_classes = target_sp)
p_synergy_reg
ggsave("../figures/diag_star1/24_synergy_calibration.png", plot = p_synergy_reg, width = 11, height = 7, dpi = 300, bg = "white")
