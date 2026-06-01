library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(patchwork)
source("./gamfunc.R")

dat_split = read.csv("../data/raw/dataset_split_2224.csv")

dat_split <- dat_split %>%
  rename(
    # Imagename              = imagename,
    altitude               = ALTIMETER_ALTITUDE_METER,
    altitude2              = ALTIMETER_ALTITUDE2_METER,
    stereo_altitude        = STEREO_ALTITUDE_METER,
    backscatter            = FLUOROMETER_BACKSCATTER_NTU,
    bottom_depth           = FATHOMETER_OCEAN_DEPTH_METER,
    cdom                   = FLUOROMETER_CDOM,
    chlorophyll            = FLUOROMETER_CHLOROPHYLL,
    # field_of_view_sq_meter = FIELD_OF_VIEW_SQ_METER,
    field_of_view_source   = FIELD_OF_VIEW_SOURCE,
    heading                = VEHICLE_MAGNETIC_HEADING,
    # latitude               = SHIP_LATITUDE,
    # longitude              = SHIP_LONGITUDE,
    # millimeter_per_pixel   = MILLIMETER_PER_PIXEL,
    o2                     = CTD2_DISSOLVED_OXYGEN,
    pitch                  = VEHICLE_PITCH_ANGLE,
    roll                   = VEHICLE_ROLL_ANGLE,
    s                      = CTD_SALINITY,
    t                      = CTD_TEMPERATURE_CELSIUS,
    fluorometer_signal     = FLUOROMETER_SIGNAL,
    datetime               = image_timestamp,
    v_depth                = CTD_VEHICLE_DEPTH_METER,
    gear                   = GEAR,
    cruise_id              = CRUISE_ID,
    habcam_pk              = HABCAM_PK
  ) %>%
  mutate(
    bottom_depth = altitude + v_depth
  ) %>%
  # mutate(
  #   Imagename = basename(Imagename)
  # ) %>%
  distinct(imagename, .keep_all = TRUE)

dat_split <- dat_split %>%
  mutate(density = n_annotations / field_of_view_sq_meter)

dat_plot <- dat_split %>%
  mutate(
    split_stage = case_when(
      is_train == "True" | is_train == TRUE ~ "Detection train",
      is_test_gam_train == "True" | is_test_gam_train == TRUE ~ "GAM train",
      is_test_gam_test == "True" | is_test_gam_test == TRUE ~ "GAM test",
      TRUE ~ "Other"
    ),
    density = n_annotations / field_of_view_sq_meter,
    density_class = case_when(
      n_annotations == 0 ~ "0 scallops",
      n_annotations == 1 ~ "1 scallop",
      n_annotations == 2 ~ "2 scallops",
      n_annotations > 2 ~ ">2 scallops",
      TRUE ~ NA_character_
    ),
    block_pos = transect_position %% 10,
    transect_group = as.factor(transect_group),
    transect_block_id = as.factor(transect_block_id)
  ) %>%
  filter(!is.na(transect_group), !is.na(transect_position))

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
      "GAM train"  = "#2C7FB8",
      "GAM test"   = "#D95F0E",
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

p_blocks

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
      "GAM train"  = "#2C7FB8",
      "GAM test"   = "#D95F0E",
      "Other"      = "grey80"
    )
  ) +
  labs(
    title = "Georges Banks and Mid-Atlantic Bight",
    # subtitle = "Grey lines represent transects",
    x = "Longitude",
    y = "Latitude",
    color = "Dataset split"
  ) +
  theme_minimal(base_size = 13)

p_map_split

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
      "GAM train"  = "#2C7FB8",
      "GAM test"   = "#D95F0E",
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

p_density

combined_figure <- (p_blocks / (p_density | p_map_split)) + 
  plot_layout(heights = c(1,1), ) +
  plot_annotation(tag_levels = 'A')

# Optional: Save the figure with ggsave
# ggsave("combined_figure.png", combined_figure, width = 12, height = 10, bg = "white")

save_cal(p_map_split, id = paste0("2224_spatial_strat"), format = "jpg", 
         outdir = "../figures/diag2/", width = 8, height = 6)
save_cal(p_blocks, id = paste0("2224_spatial_strat_grid"), format = "jpg", 
         outdir = "../figures/diag2/", width = 11.75, height = 6)
save_cal(p_density, id = paste0("2224_spatial_strat_density"), format = "jpg", 
         outdir = "../figures/diag2/", width = 8, height = 6)

save_cal(combined_figure, id = paste0("2224_spatial_strat_comb"), format = "pdf", 
         outdir = "../figures/ms_figures/", width = 13.5, height = 12.5)
  