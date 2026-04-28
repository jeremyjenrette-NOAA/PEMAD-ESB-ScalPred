library(dplyr)
library(ggplot2)
library(scales)
# ================================
# CONFIGURATION
# ================================
load("../data/processed/Cas2224.RData")
load("../data/processed/YOLOv122224.RData")

CONFIG <- list(
  
  # ---- Model selection ----
  
  model = "cas",   # options: "yolo", "cas"
  
  model_cas = Cas2224, # validated model
  model_name_cas = deparse(substitute(Cas2224)), # returns object model name
  
  model_yolo = YOLOv122224,
  model_name_yolo = deparse(substitute(YOLOv122224)),
  # ---- Year selection ----
  year = 2024,
  
  # ---- Paths ----
  paths = list(
    detections = "../data/processed/detections/",
    completed  = "../data/processed/completed/",
    metadata   = "../data/processed/metapred/meta_all2224.csv",
    inventory  = "../data/processed/img_inventory_2022_2024_compiled.rds"
  ),
  
  # ---- Output ----
  output_dir = "diagnostics/",
  
  # ---- Regional split ----
  region_split_lat = 40,
  
  # ---- Plotting ----
  facet_by_region = TRUE,
  facet_by_depth  = TRUE,
  facet_by_year   = TRUE
)

years <- 2022:2024

metadata = read.csv(CONFIG$paths$metadata)
inventory = readRDS(CONFIG$paths$inventory)

source("../s/procfunc.R")
source("../s/fitfunc.R")
source("../s/datfunc.R")
source("../s/gamfunc.R")
source("./predfunc.R")
# ================================
# RUN PIPELINE - MODEL 1
# ================================
gams_cas <- fit_calibration_gams(
  model = CONFIG$model_cas,
  model_name = CONFIG$model_name_cas
)
summary(gams_cas$GB)

detections_cas <- NULL
completed_cas  <- NULL

for (year in years) {
  CONFIG$year <- year
  data <- load_all_data(CONFIG)
  
  # Bind detections
  detections_cas <- bind_rows(
    detections_cas,
    mutate(data$detections, year = year)
  )
  
  # Bind completed
  completed_cas <- bind_rows(
    completed_cas,
    mutate(data$completed, year = year)
  )
  
  # 🔥 Free memory aggressively
  rm(data)
  gc()
}

out_pr <- evaluate_pr_models(list(CONFIG$model_cas), 
                             stratify_region = TRUE)
pr_all <- out_pr$pr_all
best_pts <- pr_all %>%
  group_by(model) %>%
  filter(f1 == max(f1, na.rm = TRUE)) %>%
  slice_max(conf, n = 1) %>%   # break ties by confidence
  ungroup()
f1conf_gb_cas = best_pts$conf[1]
f1conf_mab_cas = best_pts$conf[2]

out_cas <- compile_calibrated_image_df(
  gams = gams_cas, 
  pred_detections = detections_cas, 
  pred_allimgs = completed_cas, 
  meta   = metadata,
  f1conf_gb = f1conf_gb_cas,
  f1conf_mab = f1conf_mab_cas
) 

img_level_final_cas  <- out_cas$img_level_final 
img_level_final_f1_cas <- out_cas$img_level_final_f1
pred_calibrated_cas  <- out_cas$pred_calibrated
pred_master_cas      <- out_cas$pred_master

pred_calibrated_cas$tlx = as.numeric(pred_calibrated_cas$tlx)
pred_calibrated_cas$tly = as.numeric(pred_calibrated_cas$tly)

pred_lengths_cas <- pred_calibrated_cas %>%
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

# ---- Diagnostics ----
p_length_cas <- plot_length_distribution(pred_lengths_cas, CONFIG)
p_depth_cas  <- plot_depth_effects(pred_calibrated_cas, CONFIG)

# ================================
# RUN PIPELINE - MODEL 2
# ================================
CONFIG$model = "yolo"

gams_yolo <- fit_calibration_gams(
  model = CONFIG$model_yolo,
  model_name = CONFIG$model_name_yolo
)
summary(gams_yolo$GB)

years <- 2022:2024

detections_yolo <- NULL
completed_yolo  <- NULL

for (year in years) {
  CONFIG$year <- year
  data <- load_all_data(CONFIG)
  
  # Bind detections
  detections_yolo <- bind_rows(
    detections_yolo,
    mutate(data$detections, year = year)
  )
  
  # Bind completed
  completed_yolo <- bind_rows(
    completed_yolo,
    mutate(data$completed, year = year)
  )
  
  # 🔥 Free memory aggressively
  rm(data)
  gc()
}

out_pr <- evaluate_pr_models(list(CONFIG$model_yolo), 
                             stratify_region = TRUE)
pr_all <- out_pr$pr_all
best_pts <- pr_all %>%
  group_by(model) %>%
  filter(f1 == max(f1, na.rm = TRUE)) %>%
  slice_max(conf, n = 1) %>%   # break ties by confidence
  ungroup()
f1conf_gb_yolo = best_pts$conf[1]
f1conf_mab_yolo = best_pts$conf[2]

out_yolo <- compile_calibrated_image_df(
  gams = gams_yolo, 
  pred_detections = detections_yolo, 
  pred_allimgs = completed_yolo, 
  meta   = metadata
  # f1conf_gb = f1conf_gb_yolo,
  # f1conf_mab = f1conf_mab_yolo
) 

img_level_final_yolo  <- out_yolo$img_level_final 
# img_level_final_f1_yolo <- out_yolo$img_level_final_f1
pred_calibrated_yolo  <- out_yolo$pred_calibrated
pred_master_yolo      <- out_yolo$pred_master

pred_calibrated_yolo$tlx = as.numeric(pred_calibrated_yolo$tlx)
pred_calibrated_yolo$tly = as.numeric(pred_calibrated_yolo$tly)

pred_lengths_yolo <- pred_calibrated_yolo %>%
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

# ---- Diagnostics ----
p_length_yolo <- plot_length_distribution(pred_lengths_yolo, CONFIG)
p_depth_yolo  <- plot_depth_effects(pred_calibrated_yolo, CONFIG)

data_list = list(img_level_final_yolo=img_level_final_yolo, 
                 img_level_final_cas=img_level_final_cas, 
                 pred_calibrated_yolo=pred_calibrated_yolo, 
                 pred_calibrated_cas=pred_calibrated_cas, 
                 pred_master_yolo=pred_master_yolo, 
                 pred_master_cas=pred_master_cas,
                 img_cas=img_cas,
                 img_yolo=img_yolo,
                 img_combined=img_combined,
                 models=models)

saveRDS(data_list, "../data/processed/2224_allpred.rds")
# ================================
# RUN PIPELINE - IMAGE-LEVEL
# ================================

# get the GT predicted images for both models
img_cas  <- run_model_pipeline(CONFIG$model_cas, CONFIG$model_name_cas)
img_yolo  <- run_model_pipeline(CONFIG$model_yolo, CONFIG$model_name_yolo)

# source("../s/gammod.R")
colnames(img_combined)
disagreement
gammod

# max_gt = max(img_combined$true_number)

gam_img = models$m7
summary(gam_img)

# next, join img_level_final_cas and img_level_final_yolo by imagename

yolo_df <- img_level_final_yolo %>%
  rename(
    pred_yolo = predicted_number,
    raw_yolo  = raw_detection_number,
    conf_yolo = mean_conf,
    density_raw_yolo = density_raw,
    density_cal_yolo = density_cal
  ) %>%
  mutate(log_yolo_pred = log1p(pred_yolo)) %>%
  filter(year == 2022)

cas_df <- img_level_final_cas %>%
  rename(
    pred_cascade = predicted_number,
    raw_cas  = raw_detection_number,
    conf_cas = mean_conf,
    density_raw_cas = density_raw,
    density_cal_cas = density_cal
  )  %>%
  mutate(log_cascade_pred = log1p(pred_cascade)) %>%
  filter(year == 2022)

img_joined <- yolo_df %>%
  left_join(
    cas_df %>%
      select(imagename, pred_cascade, raw_cas, conf_cas,
             density_raw_cas, density_cal_cas, log_cascade_pred),
    by = "imagename"
  ) %>%
  mutate(
    # pred_yolo = pmin(pred_yolo, max_gt),
    # pred_cascade = pmin(pred_cascade, max_gt),
    pred_diff = pred_cascade - pred_yolo,
    log_pred_diff = log1p(abs(pred_diff))
  ) %>%
  mutate(
    pred_yolocas = predict(gam_img, newdata = ., type = "response")
    # pred_yolocas = if_else(pred_cascade > max_gt, pmin(pred_yolocas, pred_cascade), pred_yolocas)
  )

sum_df <- tibble(
  model = c("YOLO", "Cascade", "Img-Lvl Calibrated"),
  total = c(
    sum(img_joined$pred_yolo, na.rm = TRUE),
    sum(img_joined$pred_cascade, na.rm = TRUE),
    sum(img_joined$pred_yolocas, na.rm = TRUE)
  )
)

p_abundance <- ggplot(sum_df, aes(x = model, y = total, fill = model)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(aes(label = comma(round(total))),
            vjust = -0.5, size = 5) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  labs(
    x = "",
    y = "Total predicted scallops",
    title = "Comparison of total abundance estimates"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")
p_abundance

p_yolo_c <- ggplot(img_joined, aes(x = pred_yolo, y = pred_yolocas)) +
  geom_point(alpha = 0.2, size = 1.5, color = "#7570B3") +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              color = "black", linewidth = 1.2) +
  
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  # coord_cartesian(
  #   xlim = c(0, 5),
  #   ylim = c(0, 5)
  # ) +
  
  labs(
    x = "YOLO predicted count",
    y = "Image calibrated count",
    title = "Image Calibration vs. YOLO predictions"
  ) +
  
  theme_minimal(base_size = 15)
p_yolo_c


p_cas_c <- ggplot(img_joined, aes(x = pred_cascade, y = pred_yolocas)) +
  geom_point(alpha = 0.2, size = 1.5, color = "#1B9E77") +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              color = "black", linewidth = 1.2) +
  
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  # coord_cartesian(
  #   xlim = c(0, 5),
  #   ylim = c(0, 5)
  # ) +
  
  labs(
    x = "Cascade predicted count",
    y = "Image calibrated count",
    title = "Image Calibration vs. Cascade predictions"
  ) +
  
  theme_minimal(base_size = 15)
p_cas_c

p_yolocas = ggplot(img_joined, aes(x = pred_yolo, y = pred_cascade)) +
  geom_point(alpha = 0.2, size = 1.5, color = "#7570B3") +
  
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "black") +
  
  geom_smooth(method = "lm", color = "darkred", linewidth = 1) +
  
  labs(
    x = "YOLO predictions",
    y = "Cascade predictions",
    title = "Agreement between detection models"
  ) +
  
  theme_minimal(base_size = 15)

p_diff = ggplot(img_joined, aes(x = pred_diff, y = pred_yolocas)) +
  geom_point(alpha = 0.2, size = 1.5, color = "black") +
  
  labs(
    x = "Prediction difference (YOLO - Cascade)",
    y = "GAM calibrated count",
    title = "Effect of model disagreement on calibrated abundance"
  ) +
  
  theme_minimal(base_size = 15)

p_logdense = ggplot(img_joined, aes(x = pred_yolo, y = pred_yolocas)) +
  geom_point(alpha = 0.2, size = 1.5) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "YOLO (log scale)",
    y = "GAM calibrated (log scale)",
    title = "Calibration across density scales"
  ) +
  theme_minimal(base_size = 15)

save_cal(p_abundance, id = paste("2224_aubndance_est_m7",sep=""),
         outdir = "../figures/diag1", width = 6, height = 7, format = "jpg")
save_cal(p_yolo_c, id = paste("2022_yolo_compare_m7",sep=""),
         outdir = "../figures/diag1", width = 6, height = 7, format = "jpg")
save_cal(p_cas_c, id = paste("2022_cas_compare_m7",sep=""),
         outdir = "../figures/diag1", width = 6, height = 7, format = "jpg")
save_cal(p_diff, id = paste("2224_p_diff",sep=""),
         outdir = "../figures/diag1", width = 6, height = 7, format = "jpg")

sum(img_joined$pred_yolo)
sum(img_joined$pred_cascade)
sum(img_joined$pred_yolocas)

plot(img_joined$pred_yolo, img_joined$pred_yolocas)
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(img_joined$pred_cascade, img_joined$pred_yolocas)
plot(img_joined$pred_yolo, img_joined$pred_cascade)
plot(img_joined$pred_diff, img_joined$pred_yolocas)

sum(img_joined$pred_yolo) / sum(img_joined$pred_yolocas)
cor(img_joined$pred_yolo, img_joined$pred_yolocas)
cor(img_joined$pred_cascade, img_joined$pred_yolocas)
