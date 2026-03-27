# Load GAMs
########################################################################
library(tidyverse)
library(mgcv)
library(dplyr)
library(ggplot2)
library(maps)
library(lubridate)
source("../s/gamfunc.R")
source("../s/fitfunc.R")
source("./predfunc.R")

load("../data/processed/YOLOv122224.RData")

model = YOLOv122224
model_name = deparse(substitute(YOLOv122224))

gams <- fit_calibration_gams(
  model = model,
  model_name = model_name
)
gams$GB
########################################################################
# Load predictions + metadata
pred_detections = rbind(
  read.csv("../data/raw/2023_pred/detections_2023.csv") %>% mutate(year = 2023),
  read.csv("../data/raw/2024_pred/detections_2024.csv") %>% mutate(year = 2024)
)
colnames(pred_detections)
table(pred_detections$year)

pred_allimgs = rbind(
  (read.table("../data/raw/2023_pred/completed_2023.txt") %>%
  rename(Imagename = V1) %>%
  mutate(Imagename = basename(Imagename),
         year = 2023)),
  (read.table("../data/raw/2024_pred/completed_2024.txt") %>%
    rename(Imagename = V1) %>%
    mutate(Imagename = basename(Imagename),
           year = 2024))
)
table(pred_allimgs$year)

meta = rbind(
  read.csv("../data/raw/metapred/processedimages23.csv"),
  read.csv("../data/raw/metapred/processedimages24.csv")
)
colnames(meta)


########################################################################
img_inventory = readRDS("../data/processed/img_inventory_2022_2024_compiled.rds")

img_inventory <- img_inventory %>%
  mutate(
    datetime = ymd_hms(datetime_str, tz = "UTC")
  ) %>%
  arrange(datetime)

img_inventory <- img_inventory %>%
  mutate(
    dt = as.numeric(difftime(datetime, lag(datetime), units = "secs")),
    fps_inst = 1 / dt
  ) 

# view_overlap = subset(img_inventory, date_yyyymmdd=="20230615" & hour == "05" & minute == "37")

df <- img_inventory %>%
  count(year_month) %>%
  mutate(year_month = as.Date(paste0(year_month, "-01")))

df <- img_inventory %>%
  count(year_month) %>%
  arrange(year_month)

p_img_bytime = ggplot(df, aes(x = factor(year_month, levels = unique(year_month)), y = n)) +
  geom_col(fill = "#2C7BB6") +
  labs(
    x = "Year-Month",
    y = "Number of Images",
    title = "Temporal distribution of processed images 2022–2024"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_fps = ggplot(img_inventory, aes(x = factor(year), y = 1/dt, fill = factor(year))) +
  geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
  coord_cartesian(ylim = c(0, 10)) +
  labs(
    x = "Year",
    y = "Frames per second",
    title = "Distribution of Effective FPS by Year",
    fill = "Year"
  ) +
  theme_minimal(base_size = 15)

save_cal(p = p_img_bytime, id = "2224_img_bytime", width = 6, height = 6)
save_cal(p = p_fps, id = "2224_fps", width = 6)

########################################################################

out <- compile_calibrated_image_df(
  gams = gams,
  pred_detections = pred_detections,
  pred_allimgs = pred_allimgs,
  meta = meta
)

img_level_final  <- out$img_level_final
pred_calibrated  <- out$pred_calibrated
pred_master      <- out$pred_master

pred_lengths <- pred_calibrated %>%
  mutate(
    width_px  = brx - tlx,
    height_px = bry - tly,
    length_px = (width_px + height_px) / 2,
    length_mm = length_px * millimeter_per_pixel
  ) %>%
  mutate(
    length_px = (brx - tlx + bry - tly) / 2,
    length_mm = length_px * millimeter_per_pixel
  )

p_totalcal = ggplot(img_level_final, aes(x = raw_detection_number, y = predicted_number)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Raw detections (YOLO)",
    y = "Calibrated count (Σ p(detection))",
    title = "Calibration effect: raw vs calibrated counts"
  ) +
  theme_minimal()
save_cal(p = p_totalcal, outdir = "../figures/", id = "2224_totalcal", width = 6)

lims <- quantile(
  c(img_level_final$raw_detection_number,
    img_level_final$predicted_number),
  probs = c(0.05, 0.975),
  na.rm = TRUE
)

p_totalcal +
  coord_cartesian(
    xlim = lims,
    ylim = lims
  )

world <- map_data("world")

df <- img_level_final %>%
  filter(year == 2024) %>%
  mutate(density_cal_plot = density_cal + 1e-6)

p_mapcal = ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey85",
    color = "grey60",
    linewidth = 0.2
  ) +
  geom_point(
    data = df,
    aes(x = longitude, y = latitude, color = density_cal_plot),
    size = 2,
    alpha = 0.7
  ) +
  scale_color_viridis_c(
    option = "viridis",
    trans = "log10",
    limits = quantile(df$density_cal_plot, c(0, 1), na.rm = TRUE),
    oob = scales::squish,
    breaks = scales::log_breaks(n = 5),
    labels = scales::label_number_auto()
  ) +
  coord_quickmap(
    xlim = range(df$longitude, na.rm = TRUE),
    ylim = range(df$latitude, na.rm = TRUE)
  ) +
  labs(
    title = "2024 - Calibrated scallop density",
    color = "Log Density (n/m²)"
  ) +
  theme_minimal()
save_cal(p = p_mapcal, outdir = "../figures/", id = "2224_mapcal", width = 7.5)

p_densedepth = ggplot(img_level_final, aes(x = bottom_depth, y = density_cal)) +
  geom_point(alpha = 0.2) +
  # geom_smooth(method = "gam", formula = y ~ s(x), color = "blue") +
  # xlim(92, 115) +
  labs(
    x = "Bottom depth",
    y = "Density (n/m²)",
    title = "Depth vs calibrated density"
  ) +
  theme_minimal()
save_cal(p = p_densedepth, outdir = "../figures/", id = "2224_densedepth", width = 6)

p_caldepth = plot_calibration_by_depth_envelopes(
  pred_calibrated,
  depth_envelopes = list(c(10, 60), c(60, 110), c(110, 160))
)
save_cal(p = p_caldepth, outdir = "../figures/", id = "2224_caldepth", width = 7)

p_caldense = ggplot(img_level_final, aes(x = density_raw, y = density_cal)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # scale_x_log10() +
  # scale_y_log10() +
  labs(
    x = "Raw density (n/m²)",
    y = "Calibrated density (n/m²)",
    title = "Density estimate calibration"
  ) +
  theme_minimal()
save_cal(p = p_caldense, outdir = "../figures/", id = "2224_caldense", width = 7)

p_corr = pred_calibrated %>%
  mutate(correction = p_detection - conf) %>%
  ggplot(aes(x = conf, y = correction)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red") +
  labs(
    x = "YOLO confidence",
    y = "P(detection) - confidence",
    title = "Confidence correction"
  ) +
  theme_minimal()
save_cal(p = p_corr, outdir = "../figures/", id = "2224_corr", width = 7)

mean_p <- mean(pred_calibrated$p_detection, na.rm = TRUE)

p_pdist = ggplot(pred_calibrated, aes(x = p_detection)) +
  geom_histogram(
    bins = 60,
    fill = "#2C7BB6",
    color = "white",
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = mean_p,
    linetype = "dashed",
    color = "#333333"
  ) +
  annotate(
    "text",
    x = mean_p,
    y = Inf,
    label = paste0("Mean = ", round(mean_p, 3)),
    vjust = 1.5,
    hjust = -0.1,
    size = 4
  ) +
  labs(
    x = expression(P(detection)),
    y = "Count",
    title = "Distribution of p(detection)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
save_cal(p = p_pdist, outdir = "../figures/", id = "2224_dpist", width = 7)


x_pos <- max(img_level_final$bottom_depth, na.rm = TRUE)

y_raw <- mean(img_level_final$density_raw, na.rm = TRUE)
y_cal <- mean(img_level_final$density_cal, na.rm = TRUE)

p_rawcaldepth = ggplot(img_level_final, aes(x = bottom_depth)) +
  geom_smooth(aes(y = density_raw), color = "red", se = FALSE, linewidth = 1.2) +
  geom_smooth(aes(y = density_cal), color = "blue", se = FALSE, linewidth = 1.2) +
  
  annotate("text",
           x = x_pos,
           y = y_raw,
           label = "Raw density",
           color = "red",
           hjust = 1.1,
           size = 4.5) +
  
  annotate("text",
           x = x_pos,
           y = y_cal,
           label = "Calibrated density",
           color = "blue",
           hjust = 1.1,
           size = 4.5) +
  
  labs(
    x = "Bottom depth (m)",
    y = expression("Density (n / m"^2 * ")"),
    title = "Depth-dependent bias correction"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
save_cal(p = p_rawcaldepth, outdir = "../figures/", id = "2224_rawcaldepth", width = 7)

img_level_final %>%
  arrange(desc(raw_detection_number)) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(rank)) +
  geom_line(aes(y = raw_detection_number), color = "red") +
  geom_line(aes(y = predicted_number), color = "blue") +
  labs(
    x = "Image rank",
    y = "Count",
    title = "Detection yield: raw vs calibrated"
  ) +
  theme_minimal()


p_lengthdist = ggplot(pred_lengths, aes(x = length_mm, weight = p_detection)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "white", alpha = 0.9) +
  xlim(0,150) +
  labs(
    x = "Estimated scallop length (mm)",
    y = "Expected number of individuals",
    title = "Expected scallop size distribution"
  ) +
  theme_minimal()
p_lengthdist
save_cal(p = p_lengthdist, outdir = "../figures/", id = "2224_lengthdist", width = 7)

