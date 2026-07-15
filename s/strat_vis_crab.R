library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ======================================================================
# 1. Load and Prepare Plotting Data
# ======================================================================
dat_split <- read.csv("../data/raw/dataset_split_crab_multiclass.csv")

# Robustly filter for training data and sanitize logical/character flags
train_data <- dat_split %>%
  filter(as.logical(as.character(is_train)) == TRUE) %>%
  mutate(
    is_empty = as.logical(as.character(is_empty)),
    year     = factor(year)
  )

# Isolate Empty Background images
empties <- train_data %>%
  filter(is_empty == TRUE) %>%
  mutate(label = "Empty Background")

# Isolate Positive Detections and pivot to capture mixed images accurately
positives <- train_data %>%
  filter(is_empty == FALSE) %>%
  pivot_longer(
    cols = c(n_jonah_crab, n_rock_crab, n_cancer_sp),
    names_to = "species_col",
    values_to = "count"
  ) %>%
  filter(count > 0) %>%
  mutate(label = case_when(
    species_col == "n_jonah_crab" ~ "Jonah Crab",
    species_col == "n_rock_crab"  ~ "Rock Crab",
    species_col == "n_cancer_sp"  ~ "Cancer sp.",
    TRUE ~ "Other"
  ))

# Combine back into a clean, factor-ordered plotting frame
plot_data <- bind_rows(empties, positives) %>%
  mutate(label = factor(label, levels = c("Jonah Crab", "Rock Crab", "Cancer sp.", "Empty Background")))

# Colorblind-friendly, high-contrast palette definition
class_colors <- c(
  "Jonah Crab"       = "#E69F00",  # Orange
  "Rock Crab"        = "#56B4E9",  # Sky Blue
  "Cancer sp."       = "#009E73",  # Bluish Green
  "Empty Background" = "#999999"   # Slate Grey
)

# ======================================================================
# 2. Panel A: Spatiotemporal Distribution Map
# ======================================================================
p_space <- ggplot(plot_data, aes(x = longitude, y = latitude, color = label)) +
  geom_point(alpha = 0.5, size = 1.8, stroke = 0) +
  facet_wrap(~ year, ncol = 1) +
  scale_color_manual(values = class_colors) +
  labs(
    title = "Spatiotemporal Training",
    x = "Longitude (°W)",
    y = "Latitude (°N)",
    color = "Training Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# ======================================================================
# 3. Panel B: Bottom Depth Stratification
# ======================================================================
p_depth <- ggplot(plot_data, aes(x = label, y = bottom_depth, fill = label)) +
  # Thin violins show data density curves
  geom_violin(alpha = 0.3, color = NA, scale = "width") +
  # Jittered background dots reveal raw sample density patterns
  geom_jitter(aes(color = label), width = 0.15, alpha = 0.1, size = 0.6) +
  # Heavy central boxplots capture structural quantiles
  geom_boxplot(width = 0.25, color = "#222222", alpha = 0.8, outlier.shape = NA, linewidth = 0.6) +
  scale_fill_manual(values = class_colors) +
  scale_color_manual(values = class_colors) +
  # Reverse Y-axis to naturally depict down-framer depth metrics
  scale_y_reverse(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "Bottom Depth",
    x = "",
    y = "Bottom Depth (meters)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

# ======================================================================
# 4. Master Composite Layout & Export
# ======================================================================
combined_stratification <- (p_space | p_depth) +
  plot_layout(widths = c(1.4, 1.0), guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = 'bottom')
combined_stratification
# dir.create("../figures/data_diagnostics", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = "../figures/diag_crab_multi/2426_training_data_stratification.png",
  plot = combined_stratification,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

print("Success! Training stratification diagnostic graphic saved.")

# ======================================================================
# Quantitative Dataset Summary Calculation
# ======================================================================
cat("\n==================================================\n")
cat("          GLOBAL DATASET QUANTITATIVES            \n")
cat("==================================================\n")

# 1. Total Image Split Composition
image_summary <- dat_split %>%
  group_by(dataset_level1) %>%
  summarize(
    Total_Images    = n(),
    Positive_Images = sum(as.logical(as.character(is_empty)) == FALSE),
    Empty_Images    = sum(as.logical(as.character(is_empty)) == TRUE),
    .groups = "drop"
  )
print(as.data.frame(image_summary))

# 2. Training Split Species Totals
cat("\n--- Training Split Box Counts ---\n")
train_totals <- dat_split %>%
  filter(as.logical(as.character(is_train)) == TRUE) %>%
  summarize(
    Jonah_Crab_Boxes = sum(n_jonah_crab, na.rm = TRUE),
    Rock_Crab_Boxes  = sum(n_rock_crab, na.rm = TRUE),
    Cancer_sp_Boxes  = sum(n_cancer_sp, na.rm = TRUE),
    Total_Boxes      = Jonah_Crab_Boxes + Rock_Crab_Boxes + Cancer_sp_Boxes
  )
print(as.data.frame(train_totals))

# 3. Validation Split Species Totals (How many validated)
cat("\n--- Validation Split Box Counts ---\n")
val_totals <- dat_split %>%
  filter(as.logical(as.character(is_test)) == TRUE) %>%
  summarize(
    Jonah_Crab_Boxes = sum(n_jonah_crab, na.rm = TRUE),
    Rock_Crab_Boxes  = sum(n_rock_crab, na.rm = TRUE),
    Cancer_sp_Boxes  = sum(n_cancer_sp, na.rm = TRUE),
    Total_Boxes      = Jonah_Crab_Boxes + Rock_Crab_Boxes + Cancer_sp_Boxes
  )
print(as.data.frame(val_totals))
cat("==================================================\n")

