library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(patchwork)
library(maps) # Explicitly loaded for map_data("world")
source("./gamfunc.R")

# Load the newly formatted dataset
dat_split <- read.csv("../data/raw/dataset_split_2226.csv")

# ======================================================================
# 1. Safe Fallbacks for Transect Logic
# Just in case these were excluded from the final CSV export
# ======================================================================
if(!"transect_group" %in% names(dat_split)) {
  dat_split <- dat_split %>% 
    mutate(transect_group = sub("_block_.*", "", transect_block_id))
}

if(!"transect_position" %in% names(dat_split)) {
  dat_split <- dat_split %>%
    group_by(transect_group) %>%
    mutate(transect_position = row_number() - 1) %>%
    ungroup()
}

# ======================================================================
# 2. Prepare Plotting Data
# ======================================================================
dat_plot <- dat_split %>%
  mutate(
    # Directly map our new dataset column to your plot labels
    split_stage = case_when(
      dataset == "train" ~ "Detection train",
      dataset == "test_GAM_train" ~ "Calibration train",
      dataset == "test_GAM_test" ~ "Calibration test",
      TRUE ~ "Other"
    ),
    density_class = case_when(
      n_annotations == 0 ~ "0 scallops",
      n_annotations == 1 ~ "1 scallop",
      n_annotations == 2 ~ "2 scallops",
      n_annotations > 2 ~ ">2 scallops",
      TRUE ~ NA_character_
    ),
    transect_group = as.factor(transect_group),
    transect_block_id = as.factor(transect_block_id)
  ) %>%
  filter(!is.na(transect_group), !is.na(transect_position))

# ======================================================================
# 3. Transect Block Plot (p_blocks)
# ======================================================================
p_blocks <- ggplot(
  dat_plot,
  aes(
    x = transect_position,
    y = fct_reorder(transect_group, as.numeric(transect_group)),
    fill = split_stage
  )
) +
  geom_tile(
    aes(alpha = pmin(n_annotations, 15)),
    height = 0.8,
    width = 0.95
  ) +
  scale_alpha_continuous(
    name = "Manual annotations\n(max = 15+)",
    range = c(0.4, 1)
  ) +
  scale_fill_manual(
    values = c(
      "Detection train" = "black",
      "Calibration train"  = "#2C7FB8",
      "Calibration test"   = "#D95F0E",
      "Other"      = "grey80"
    )
  ) +
  labs(
    title = "Transect-block stratification of annotated HabCam images",
    subtitle = "1 out of every 50 images annotated",
    x = "Image position within transect",
    y = "Transect",
    fill = "Dataset split"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

print(p_blocks)

# ======================================================================
# 4. Geospatial Map Plot (p_map_split)
# ======================================================================
world <- map_data("world")

p_map_split <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey55",
    color = "grey40",
    linewidth = 0.2
  ) +
  geom_path(
    data = dat_plot %>% arrange(transect_group, transect_position),
    aes(
      x = longitude,
      y = latitude,
      group = transect_group
    ),
    color = "grey60",
    linewidth = 0.25,
    alpha = 0.45
  ) +
  geom_point(
    data = dat_plot,
    aes(
      x = longitude,
      y = latitude,
      color = split_stage,
      size = pmin(n_annotations, 15)
    ),
    alpha = 0.75
  ) +
  coord_quickmap(
    xlim = range(dat_plot$longitude, na.rm = TRUE) + c(-0.4, 0.4),
    ylim = range(dat_plot$latitude, na.rm = TRUE) + c(-0.4, 0.4)
  ) +
  scale_size_continuous(
    name = "Manual annotations\n(max = 15+)",
    range = c(0.4, 3.2)
  ) +
  scale_color_manual(
    values = c(
      "Detection train" = "black",
      "Calibration train"  = "#2C7FB8",
      "Calibration test"   = "#D95F0E",
      "Other"      = "grey80"
    )
  ) +
  labs(
    title = "Georges Banks and Mid-Atlantic Bight",
    x = "Longitude",
    y = "Latitude",
    color = "Dataset split"
  ) +
  theme_minimal(base_size = 13)

print(p_map_split)

# ======================================================================
# 5. Density Violin Plot (p_density)
# ======================================================================
p_density <- ggplot(
  dat_plot,
  aes(x = split_stage, y = n_annotations, fill = split_stage)
) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.alpha = 0.15) +
  scale_y_continuous(
    trans = "pseudo_log",
    breaks = c(0, 1, 2, 5, 10, 25, 50, 100, 250)
  ) +
  scale_fill_manual(
    values = c(
      "Detection train" = "#4D4D4D",
      "Calibration train"  = "#2C7FB8",
      "Calibration test"   = "#D95F0E",
      "Other"      = "grey80"
    )
  ) +
  labs(
    title = "Scallop density across datasets",
    x = "",
    y = "Scallop annotations per image",
    fill = "Dataset split"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p_density)

# ======================================================================
# 6. Combined Figure via Patchwork
# ======================================================================
combined_figure <- (p_blocks / p_density) + 
  plot_layout(heights = c(1,1)) +
  plot_annotation(tag_levels = 'A')
combined_figure
# Save outputs
# save_cal(p_map_split, id = paste0("2226_spatial_strat"), format = "jpg", 
#          outdir = "../figures/diag2/", width = 8, height = 6)
# save_cal(p_blocks, id = paste0("2226_spatial_strat_grid"), format = "jpg", 
#          outdir = "../figures/diag2/", width = 11.75, height = 6)
# save_cal(p_density, id = paste0("2226_spatial_strat_density"), format = "jpg", 
#          outdir = "../figures/diag2/", width = 8, height = 6)

save_cal(combined_figure, id = paste0("2226_spatial_strat_comb"), format = "png", 
         outdir = "../figures/diag3/", width = 10, height = 12.5)
save_cal(combined_figure, id = paste0("2226_spatial_strat_comb"), format = "pdf", 
         outdir = "../figures/ms_figures/", width = 13.5, height = 12.5)