# Load GAMs
########################################################################
library(tidyverse)
library(mgcv)
library(dplyr)
library(ggplot2)
library(maps)
library(scales)
library(lubridate)
source("../s/gamfunc.R")
source("../s/fitfunc.R")
source("../s/procfunc.R")
source("./predfunc.R")

load("../data/processed/YOLOv122224.RData")

model = YOLOv122224
model_name = deparse(substitute(YOLOv122224))

gams <- fit_calibration_gams(
  model = model,
  model_name = model_name
)
summary(gams$GB)
########################################################################
# Load predictions + metadata
pred_detections = rbind(
  read.csv("../data/raw/2022_pred/detections_2022.csv") %>% mutate(year = 2022),
  read.csv("../data/raw/2023_pred/detections_2023.csv") %>% mutate(year = 2023),
  read.csv("../data/raw/2024_pred/detections_2024.csv") %>% mutate(year = 2024)
)
colnames(pred_detections)
table(pred_detections$year)

pred_allimgs = rbind(
  (read.table("../data/raw/2022_pred/completed_2022.txt") %>%
     rename(Imagename = V1) %>%
     mutate(Imagename = basename(Imagename),
            year = 2022)),
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
  read.csv("../data/raw/metapred/processedimages22.csv"),
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
  scale_y_continuous(
    labels = label_number(scale = 1e-6, suffix = " mil", accuracy = 0.1)
  ) +
  labs(
    x = "Year-Month",
    y = "Number of Images",
    title = "HabCam survey images 2022–2024"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(face = "bold")
  )
p_img_bytime

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
out_pr <- evaluate_pr_models(list(model), 
                             stratify_region = TRUE)
pr_all <- out_pr$pr_all
best_pts <- pr_all %>%
  group_by(model) %>%
  filter(f1 == max(f1, na.rm = TRUE)) %>%
  slice_max(conf, n = 1) %>%   # break ties by confidence
  ungroup()
f1conf_gb = best_pts$conf[1]
f1conf_mab = best_pts$conf[2]

out <- compile_calibrated_image_df(
  gams = gams,
  pred_detections = pred_detections,
  pred_allimgs = pred_allimgs,
  meta = meta,
  f1conf_gb = f1conf_gb,
  f1conf_mab = f1conf_mab
)

img_level_final  <- out$img_level_final
img_level_final_f1 <- out$img_level_final_f1
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

year_summary <- img_level_final %>%
  group_by(year) %>%
  summarise(
    total_p_detection = sum(predicted_number, na.rm = TRUE),
    total_density_cal = sum(density_cal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    img_level_final_f1 %>%
      group_by(year) %>%
      summarise(
        total_raw_detection_f1 = sum(predicted_number_f1, na.rm = TRUE),
        total_density_f1 = sum(density_f1, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "year"
  )

plot_df <- year_summary %>%
  select(year, total_density_cal, total_density_f1) %>%
  pivot_longer(
    cols = c(total_density_cal, total_density_f1),
    names_to = "method",
    values_to = "density"
  ) %>%
  mutate(
    method = recode(method,
                    total_density_cal = "GAM (calibrated)",
                    total_density_f1  = "F1 threshold"
    )
  )

p_summ = ggplot(plot_df, aes(x = factor(year), y = density, fill = method)) +
  geom_col(position = "dodge") +
  scale_fill_manual(
    values = c(
      "GAM (calibrated)" = "#1B9E77",
      "F1 threshold"     = "#D95F02"
    )
  ) +
  scale_y_continuous(
    labels = scales::label_number(scale = 1e-3, suffix = "k")
  ) +
  labs(
    x = "Year",
    y = "Total scallop density (n/m²)",
    fill = "Method",
    title = "Comparison of abundance estimation methods by year"
  ) +
  theme_minimal()
save_cal(p = p_summ, outdir = "../figures/", id = "2224_summ", width = 6)

year_summary %>%
  mutate(ratio = total_density_f1 / total_density_cal) %>%
  ggplot(aes(x = factor(year), y = ratio)) +
  geom_col(fill = "#D95F02") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Year",
    y = "F1 / GAM ratio",
    title = "Relative difference between thresholded and calibrated estimates"
  ) +
  theme_minimal()
########################################################################

p_totalcal = ggplot(img_level_final, aes(x = raw_detection_number, y = predicted_number)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Detection (YOLO)",
    y = "Predicted count",
    title = "Detection vs calibrated count"
  ) +
  theme_minimal()
save_cal(p = p_totalcal, outdir = "../figures/", id = "2224_totalcal", width = 6)

# 1. Create depth bins and counts
depth_bins <- img_level_final %>%
  filter(!is.na(bottom_depth), !is.na(density_cal)) %>%
  mutate(depth_bin = cut(bottom_depth, breaks = 19)) %>%
  group_by(depth_bin) %>%
  summarise(
    depth_mid = mean(bottom_depth, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(depth_mid)

# 2. Rescale counts onto the density axis
max_density <- max(img_level_final$density_cal, na.rm = TRUE)
max_n <- max(depth_bins$n, na.rm = TRUE)

scale_factor <- max_density / max_n

depth_bins <- depth_bins %>%
  mutate(n_scaled = n * scale_factor)

# 3. Plot
p_densedepth <- ggplot() +
  geom_col(
    data = depth_bins,
    aes(x = depth_mid, y = n_scaled),
    width = diff(range(img_level_final$bottom_depth, na.rm = TRUE)) / 20 * 0.9,
    fill = "grey80",
    color = "#635758",
    alpha = 1
  ) +
  geom_point(
    data = img_level_final,
    aes(x = bottom_depth, y = density_cal),
    alpha = 0.25,
    size = 1.6,
    color = "#1F4E79"
  ) +
  scale_y_continuous(
    name = expression("Density (n / m"^2 * ")"),
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Number of images",
      labels = label_number(scale = 1e-3, suffix = "k", accuracy = 1)
    )
  ) +
  labs(
    x = "Bottom depth (m)",
    title = "Predicted scallop density at depth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.y.left = element_text(color = "#1F4E79"),
    axis.title.y.right = element_text(color = "#635758"),
    panel.grid.minor = element_blank()
  )

p_densedepth
save_cal(p = p_densedepth, outdir = "../figures/", id = "2224_densedepth", width = 7)

world <- map_data("world")

df <- img_level_final %>%
  filter(year %in% c(2022, 2024)) %>%
  mutate(
    density_cal_plot = density_cal + 1e-6,
    year = factor(year, levels = c(2022, 2024))
  )

p_mapcal <- ggplot() +
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
    # limits = quantile(df$density_cal_plot, c(0, 1), na.rm = TRUE),
    oob = scales::squish,
    breaks = scales::log_breaks(n = 5),
    labels = scales::label_number_auto()
  ) +
  coord_quickmap(
    xlim = range(df$longitude, na.rm = TRUE),
    ylim = range(df$latitude, na.rm = TRUE)
  ) +
  facet_wrap(~ year, ncol = 1) +
  labs(
    title = "Predicted scallop density",
    color = "Density (n / m²)"
  ) +
  theme_minimal()
save_cal(p = p_mapcal, outdir = "../figures/", id = "2224_mapcal", width = 7)

# should facet this by region
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
    x = "Detection confidence",
    y = "P(detection) - confidence",
    title = "Confidence correction (YOLO)"
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
    title = "Distribution of p(detection)",
    subtitle = "YOLOv12"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
save_cal(p = p_pdist, outdir = "../figures/", id = "2224_dpist", width = 7)

p_lengthdist <- ggplot(pred_lengths, aes(x = length_mm, weight = p_detection)) +
  geom_histogram(
    bins = 50,
    fill = "darkgreen",
    color = "white",
    alpha = 0.9
  ) +
  facet_wrap(~ year) +
  xlim(0, 150) +
  scale_y_continuous(
    labels = label_number(scale = 1e-3, suffix = "k", accuracy = 1)
  ) +
  labs(
    x = "Estimated scallop length (mm)",
    y = "Expected number of individuals",
    title = "Expected scallop size distribution by year"
  ) +
  theme_minimal()
p_lengthdist
save_cal(p = p_lengthdist, outdir = "../figures/", id = "2224_lengthdist", width = 7)

