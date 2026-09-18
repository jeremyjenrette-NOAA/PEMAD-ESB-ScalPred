library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
source("./procfunc.R")

# ============================================================#
# 1. Load Unified Data and Evaluate
# ============================================================#
models <- readRDS("../data/processed/eval_2226.rds")

# Toggle 'stratify_region = FALSE' if you want pooled overall performance
out_pr <- evaluate_pr_models(models, stratify_region = TRUE, conf_grid = seq(0, 1, by = 0.0005))

pr_all <- out_pr$pr_all
p_pr   <- out_pr$p_pr
map_summary <- out_pr$map_summary

print("=== Average Precision (mAP) Summary ===")
print(map_summary) # at mAP@0.1

# ============================================================#
# 2. Optimal F1 Threshold Calculation
# ============================================================#
best_pts <- pr_all %>%
  group_by(model) %>%
  filter(f1 == max(f1, na.rm = TRUE)) %>%
  slice_max(conf, n = 1) %>%   # break ties by confidence
  ungroup()

print("=== Optimal F1 Thresholds ===")
print(best_pts %>% select(model, conf, f1, precision, recall))

# ============================================================#
# 3. F1 Curve Visualization (Panel B)
# ============================================================#
p_f1 <- ggplot(pr_all, aes(x = conf, y = f1, color = model)) +
  geom_line(linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "F1 Score vs Confidence Threshold",
    x = "Confidence threshold",
    y = "F1 score",
    color = "Model by region"
  )

# ============================================================#
# 4. Prepare mAP Data and Plot (Panel C)
# ============================================================#
map_plot_data <- map_summary %>%
  mutate(
    Region = case_when(
      str_detect(model, " GB$") ~ "Georges Bank",
      str_detect(model, " MAB$") ~ "Mid-Atlantic Bight",
      TRUE ~ "Overall"
    ),
    Architecture = str_remove(model, " GB$| MAB$")
  )

p_map <- ggplot(map_plot_data, aes(x = Architecture, y = Average_Precision, fill = Region)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.3f", Average_Precision)), 
    position = position_dodge(width = 0.8), 
    vjust = -0.8, 
    size = 4, 
    fontface = "bold"
  ) +
  scale_y_continuous(limits = c(0, 1.05), labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Georges Bank" = "#2C7FB8", "Mid-Atlantic Bight" = "#D95F0E", "Overall" = "grey50")) +
  labs(
    title = "Average Precision (mAP@10)",
    x = "",
    y = "mAP Score",
    fill = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# ============================================================#
# 5. Build the Error Composition Chart (Panel D)
# ============================================================#
error_plot_data <- best_pts %>%
  mutate(
    Region = case_when(
      str_detect(model, " GB$") ~ "Georges Bank",
      str_detect(model, " MAB$") ~ "Mid-Atlantic Bight",
      TRUE ~ "Overall"
    ),
    Architecture = str_remove(model, " GB$| MAB$")
  ) %>%
  select(Architecture, Region, FP, FN) %>%
  tidyr::pivot_longer(cols = c(FP, FN), names_to = "Error_Type", values_to = "Count") %>%
  mutate(
    Error_Type = recode(Error_Type, 
                        "FP" = "False Positives", 
                        "FN" = "False Negatives")
  )

p_error <- ggplot(error_plot_data, aes(x = Architecture, y = Count, fill = Error_Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
  geom_text(
    aes(label = Count), 
    position = position_dodge(width = 0.8), 
    vjust = -0.5, 
    size = 3.5, 
    fontface = "bold"
  ) +
  facet_wrap(~ Region) +
  # Use distinct colors to separate from the regional palette used in p_map
  scale_fill_manual(values = c("False Positives" = "#E41A1C", "False Negatives" = "#377EB8")) +
  # Expand the top of the Y-axis so text labels don't get cut off
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  labs(
    title = "Model Errors at Optimal F1 Threshold",
    x = "",
    y = "Raw Error Count",
    fill = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )
p_error
# ============================================================#
# 6. Final 4-Panel Master Assembly
# ============================================================#
# By keeping top_row and bottom_row separate before combining, 
# we prevent patchwork from awkwardly mashing all 3 different legends together.
top_row <- (p_pr | p_f1) + plot_layout(guides = "collect")
bottom_row <- (p_map | p_error) 

combined_4panel <- (top_row / bottom_row) + 
  plot_annotation(tag_levels = 'A')

combined_2panel <- (top_row) + 
  plot_annotation(tag_levels = 'A')

# print(combined_4panel)

# ============================================================#
# 7. Save
# ============================================================#
ggsave(
  filename = "../figures/ms_figures/2226_model_performance_summary.pdf", 
  plot = combined_2panel, 
  width = 9, 
  height = 5, 
  device = "pdf",
  bg = "white"
)
# ggsave(
#   filename = "~/saltnoaa/presentations/figures/2226_model_performance_summary.png", 
#   plot = combined_2panel, 
#   width = 11, 
#   height = 5, 
#   bg = "white"
# )
