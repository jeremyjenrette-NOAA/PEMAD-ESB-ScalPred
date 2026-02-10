library(tidyverse)
source("./gamfunc.R")

load("../data/processed/Cas2024v2.RData")
load("../data/processed/YOLO2024.RData")

model = Cas2024v2
model_name = deparse(substitute(Cas2024v2))

gams <- fit_calibration_gams(
  model = model,
  model_name = model_name
)

pred_gb <- predict_calibration_by_depth(
  gam        = gams$GB,
  calib_df   = get_calib_data(model, model_name, "GB"),
  depth_breaks = c(40, 80, 120),
  region     = "GB",
  model_name = model_name
)

pred_mab <- predict_calibration_by_depth(
  gams$MAB,
  calib_df   = get_calib_data(model, model_name, "MAB"),
  depth_breaks = c(40, 50, 60, 70),
  region = "MAB",
  model_name = model_name
)

pred_all <- bind_rows(pred_gb, pred_mab)

p_cal <- plot_calibration_by_depth(pred_all, model_name)
p_cal



res_gb <- compute_image_level_counts(
  calib_df  = get_calib_data(model, model_name, "GB"),
  gam       = gams$GB,
  region    = "GB",
  model_name = model_name
)

res_mab <- compute_image_level_counts(
  calib_df  = get_calib_data(model, model_name, "MAB"),
  gam       = gams$MAB,
  region    = "MAB",
  model_name = model_name
)

str(res_gb)

img_all <- bind_rows(
  res_gb$img,
  res_mab$img
)

metrics_all <- bind_rows(
  res_gb$metrics |> mutate(region = "GB"),
  res_mab$metrics |> mutate(region = "MAB")
)

p_counts <- plot_image_level_fit(
  img_df     = img_all,
  metrics_df = metrics_all
)

p_counts

p_counts_zoomed = plot_image_level_fit_zoom(img_all,
                          metrics_all,
                          model_name,
                          zoom_q = 0.85)
p_counts_zoomed

save_cal(p_cal, id = paste("bydepth_",model_name,sep=""), 
         outdir = "../figures/")

save_cal(p_counts, id = paste("fit_",model_name,sep=""), 
         outdir = "../figures/")

save_cal(p_counts_zoomed, id = paste("fit_zoomed_",model_name,sep=""), 
         outdir = "../figures/")
