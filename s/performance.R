# ============================================================#
# performance.R
# ============================================================#
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
source("./procfunc.R")

# ============================================================#
# 1. Load Unified Data and Evaluate
# ============================================================#
models <- readRDS("../data/processed/seal_eval_2026.rds")

# Evaluate class-specific PR curves
out_pr <- evaluate_pr_models(models, stratify_region = FALSE, conf_grid = seq(0, 1, by = 0.005))

pr_all      <- out_pr$pr_all
p_pr        <- out_pr$p_pr
map_summary <- out_pr$map_summary

print("=== Class-Specific Average Precision (AP) Summary ===")
print(map_summary)

# ============================================================#
# 2. Optimal F1 Threshold Calculation per Model/Species
# ============================================================#
best_pts <- pr_all %>%
  group_by(model, species) %>%
  filter(f1 == max(f1, na.rm = TRUE)) %>%
  slice_max(conf, n = 1) %>% # Break ties by confidence
  ungroup()

print("=== Optimal F1 Thresholds ===")
print(best_pts %>% select(model, species, conf, f1, precision, recall))

# ============================================================#
# 3. Plot mAP Comparison (Panel C)
# ============================================================#
p_map <- ggplot(map_summary, aes(x = model, y = Average_Precision, fill = species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.3f", Average_Precision)), 
    position = position_dodge(width = 0.8), 
    vjust = -0.8, 
    size = 3.5, 
    fontface = "bold"
  ) +
  scale_y_continuous(limits = c(0, 1.05), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Class-Specific Average Precision (AP)",
    x = "Model Architecture",
    y = "AP Score",
    fill = "Species"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# ============================================================#
# 4. Multi-Class Confusion Matrix Construction (Panel D)
# ============================================================#
# Evaluate confusion matrices at a standard working threshold (e.g., conf = 0.10)
# to see where classification confusion and misses happen
conf_data <- map_dfr(names(models), function(mod_name) {
  generate_confusion_matrix(
    det_data = models[[mod_name]]$det,
    img_data = models[[mod_name]]$img,
    threshold = 0.10 # Working threshold
  ) %>% mutate(model = mod_name)
})

# Normalize within actual classes to display row-wise percentages
conf_normalized <- conf_data %>%
  group_by(model, actual) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    actual = factor(actual, levels = c("adult", "pup", "background")),
    predicted = factor(predicted, levels = c("adult", "pup", "missed"))
  )

p_conf <- ggplot(conf_normalized, aes(x = predicted, y = actual, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%d\n(%.1f%%)", n, pct)), 
    color = "black", 
    fontface = "bold", 
    size = 3
  ) +
  facet_wrap(~ model) +
  scale_fill_gradient(low = "#F7FCF0", high = "#7BCCC4", name = "% of Class") +
  labs(
    title = "Confusion Matrices (Conf >= 0.10)",
    x = "Predicted Class",
    y = "True Class"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ============================================================#
# 5. Master Multi-Panel Assembly
# ============================================================#
top_row <- p_pr
bottom_row <- (p_map / p_conf) + plot_layout(heights = c(1, 2.2))

combined_layout <- (top_row | bottom_row) + 
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = 'A')
combined_layout
# ============================================================#
# 6. Save Outputs
# ============================================================#
# dir.create("../figures/diag_crab1", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "../figures/diag_seal1/26_seal_performance.png", 
  plot = combined_layout, 
  width = 11, 
  height = 9, 
  device = "png",
  bg = "white"
)
print("Success! Multi-class metrics saved directly to diagnostics directory.")