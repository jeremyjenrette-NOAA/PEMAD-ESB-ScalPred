library(tidyverse)
library(mgcv)
source("./gamfunc.R")
source("./fitfunc.R")
source("./procfunc.R")

load("../data/processed/Cas2224.RData")
load("../data/processed/YOLOv122224.RData")

model = YOLOv12strat2224
model_name = deparse(substitute(YOLOv12strat2224))

gams <- fit_calibration_gams(
  model = model,
  model_name = model_name
)

pred_gb <- predict_calibration_by_depth(
  gam        = gams$GB,
  calib_df   = get_data(model, model_name, "GB"),
  depth_breaks = c(40, 80, 120),
  region     = "GB",
  model_name = model_name
)

pred_mab <- predict_calibration_by_depth(
  gam = gams$MAB,
  calib_df   = get_data(model, model_name, "MAB"),
  depth_breaks = c(40, 50, 60, 70),
  region = "MAB",
  model_name = model_name
)

pred_all <- bind_rows(pred_gb, pred_mab)

p_cal <- plot_calibration_by_depth(pred_all, model_name)
p_cal
#######################################################
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
#######################################################
res_gb <- compute_image_level_counts(
  calib_df  = get_data(model, model_name, "GB"),
  gam       = gams$GB,
  region    = "GB",
  model_name = model_name,
  f1_thresh = f1conf_gb
)

res_mab <- compute_image_level_counts(
  calib_df  = get_data(model, model_name, "MAB"),
  gam       = gams$MAB,
  region    = "MAB",
  model_name = model_name,
  f1_thresh = f1conf_mab
)

img_all <- bind_rows(
  res_gb$img,
  res_mab$img
)

metrics_all <- dplyr::bind_rows(
  res_gb$metrics  |> dplyr::mutate(region = "GB"),
  res_mab$metrics |> dplyr::mutate(region = "MAB"),
  compute_combined_metrics(img_all, region_name = "All Survey Regions")
)

p_counts <- plot_image_level_fit(
  img_df = img_all, # |> dplyr::mutate(region = "All Surveys 2022"),
  metrics_df = metrics_all[1:2,],
  model_name = model_name,
  include_f1 = FALSE,
  plot_title = paste(model_name, " True vs. Predicted Count", sep = " -")
)
p_counts

p_counts_pooled <- plot_image_level_fit(
  img_df = img_all |> dplyr::mutate(region = "All Survey Regions"),
  metrics_df = metrics_all[3,],
  model_name = model_name,
  include_f1 = FALSE,
  plot_title = paste(model_name, " True vs. Predicted Count", sep = " -")
)
p_counts_pooled = p_counts_pooled + labs(title = "True vs. Predicted Count",
subtitle = "YOLOv12")

p_zoomed <- plot_image_level_fit_zoom(
  # img_df = img_all |> dplyr::mutate(region = "All Survey Regions"),
  img_df = img_all,
  metrics_df = metrics_all[1:2,],
  model_name = model_name,
  zoom_q=0.95
)
p_zoomed

save_cal(p_cal, id = paste("bydepth_",model_name,sep=""), 
         outdir = "../figures/diag1/", format = "png", width = 7)

save_cal(p_counts, id = paste("fit_",model_name,sep=""), 
         outdir = "../figures/", format = "png")

save_cal(p_counts_pooled, id = paste("fitpooled_",model_name,sep=""), 
         outdir = "../figures/diag1/", width = 5, height = 7, format = "png")

save_cal(p_zoomed, id = paste("fit_zoomedALL_",model_name,sep=""), 
         outdir = "../figures/diag1/", width = 5, height = 7, format = "png")
######
vis.gam(gams$GB,
        view = c("conf", "bottom_depth"),
        plot.type = "contour")
vis.gam(gams$GB,
        view = c("conf", "bottom_depth"),
        plot.type = "persp")
