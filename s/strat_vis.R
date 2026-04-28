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

world <- map_data("world")

dat_plot <- dat_split %>%
  mutate(
    density = n_annotations / field_of_view_sq_meter,
    split_stage = case_when(
      is_train == "True" | is_train == TRUE ~ "YOLO train (70%)",
      is_test_gam_train == "True" | is_test_gam_train == TRUE ~ "GAM train (65% of test)",
      is_test_gam_test == "True" | is_test_gam_test == TRUE ~ "GAM test (35% of test)",
      TRUE ~ "Other"
    )
  ) %>%
  filter(
    !is.na(latitude),
    !is.na(longitude),
    !is.na(density)
  )

spatial_summary <- dat_plot %>%
  group_by(stratum, split_stage) %>%
  summarise(
    lon = mean(longitude, na.rm = TRUE),
    lat = mean(latitude, na.rm = TRUE),
    mean_density = mean(density, na.rm = TRUE),
    n_images = n(),
    .groups = "drop"
  )

p_strat = ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey75",
    linewidth = 0.2
  ) +
  geom_point(
    data = spatial_summary,
    aes(
      x = lon,
      y = lat,
      size = mean_density,
      color = split_stage
    ),
    alpha = 0.8
  ) +
  coord_quickmap(
    xlim = range(dat_plot$longitude, na.rm = TRUE) + c(-0.5, 0.5),
    ylim = range(dat_plot$latitude, na.rm = TRUE) + c(-0.5, 0.5)
  ) +
  scale_size_continuous(name = "Mean density\n(n/m²)") +
  labs(
    title = "Spatial distribution of stratified scallop image splits",
    x = "Longitude",
    y = "Latitude",
    color = "Dataset split"
  ) +
  theme_minimal()

save_cal(p_strat, id = paste0("2224_spatial_strat"), format = "jpg", 
         outdir = "../figures/diag2/", width = 8, height = 6)


