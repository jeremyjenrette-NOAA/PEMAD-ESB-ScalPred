library(dplyr)
library(tidyr)
library(ggplot2)

# ======================================================================
# 1. Load Data & Isolate True Positives
# ======================================================================
models <- readRDS("../data/processed/crab_evalmulti_2426.rds")

# Bind detection dataframes from both architectures
raw_detections <- bind_rows(
  models$YOLOv12$det %>% mutate(model = "YOLOv12"),
  models$`Cascade R-CNN`$det %>% mutate(model = "Cascade R-CNN")
)

# Standardize names and isolate True Positives
tp_data <- raw_detections %>%
  filter(truedetect == TRUE) %>%
  mutate(
    gt_label = tolower(trimws(gt_label)),
    spname   = tolower(trimws(spname))
  )

# ======================================================================
# 2. Define Taxonomic Resolution Shifts
# ======================================================================
resolution_data <- tp_data %>%
  # Determine if the human annotator labeled Genus or Species
  mutate(
    gt_level = if_else(gt_label %in% c("cancer_sp", "cancer_sp."), "Genus", "Species")
  ) %>%
  # Evaluate the classification shift between human and model
  mutate(
    outcome = case_when(
      # CASE 1: Human labeled Genus (cancer_sp)
      gt_level == "Genus" & spname == "cancer_sp" ~ "Conservative (Agreed on Genus)",
      gt_level == "Genus" & spname %in% c("jonah_crab", "rock_crab") ~ "Audited (Upgraded to Species)",
      
      # CASE 2: Human labeled Specific Species (jonah / rock)
      gt_level == "Species" & spname == gt_label ~ "Precise (Agreed on Species)",
      gt_level == "Species" & spname == "cancer_sp" ~ "Generalized (Fallback to Genus)",
      gt_level == "Species" & spname %in% c("jonah_crab", "rock_crab") & spname != gt_label ~ "Species Confusion",
      
      TRUE ~ "Other"
    ),
    # Prettier names for the plot facets
    gt_facet_label = if_else(
      gt_level == "Genus", 
      "Ground Truth: Genus (cancer_sp)", 
      "Ground Truth: Species (jonah / rock)"
    )
  )

# ======================================================================
# 3. Calculate Percentages & Filter Labels
# ======================================================================
resolution_summary <- resolution_data %>%
  group_by(model, gt_facet_label, outcome) %>%
  summarize(count = n(), .groups = "drop_last") %>%
  mutate(
    percentage = count / sum(count) * 100,
    # Clean up labels. Hide text labels if the stack is less than 4% wide
    # to prevent cramped, overlapping text in the final output.
    label_text = if_else(percentage >= 4.0, sprintf("%d\n(%.1f%%)", count, percentage), "")
  ) %>%
  ungroup()

# Set factor levels to control the stack sequence logically
resolution_summary$outcome <- factor(
  resolution_summary$outcome,
  levels = c(
    "Precise (Agreed on Species)",
    "Audited (Upgraded to Species)",
    "Conservative (Agreed on Genus)",
    "Generalized (Fallback to Genus)",
    "Species Confusion"
  )
)

# ======================================================================
# 4. Color Palette Mapping
# ======================================================================
resolution_palette <- c(
  "Precise (Agreed on Species)"     = "#2ca02c",  # Deep Forest Green (Success)
  "Audited (Upgraded to Species)"   = "#98df8a",  # Soft Mint Green (Resolved Upgrade)
  "Conservative (Agreed on Genus)"  = "#bcbd22",  # Muted Yellow-Green (Neutral Genus)
  "Generalized (Fallback to Genus)" = "#999999",  # Slate Grey (Uncertainty Fallback)
  "Species Confusion"               = "#d62728"   # Deep Red (Classification Error)
)

# ======================================================================
# 5. Build and Save the Combined Audit Plot
# ======================================================================
p_taxonomic <- ggplot(resolution_summary, aes(x = model, y = percentage, fill = outcome)) +
  geom_col(position = "stack", width = 0.65, color = "black", linewidth = 0.5) +
  geom_text(
    aes(label = label_text),
    position = position_stack(vjust = 0.5),
    color = "black",
    fontface = "bold",
    size = 3.5
  ) +
  facet_wrap(~ gt_facet_label, scales = "free_x") +
  scale_fill_manual(values = resolution_palette) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Taxonomic Resolution",
    subtitle = "Evaluating taxonomic upgrades (auditing) vs. downgrades (generalization) under true detections",
    x = "Model Architecture",
    y = "Percentage of True Positives",
    fill = "Taxonomic Outcome"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 11, color = "grey30", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
p_taxonomic
# Export visual
# dir.create("../figures/diagnostics", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = "../figures/diag_crab_multi/2426_taxonomic_resolution_audit.png",
  plot = p_taxonomic,
  width = 13,
  height = 8,
  dpi = 300,
  bg = "white"
)

print("Success! Unified taxonomic resolution graphic generated.")