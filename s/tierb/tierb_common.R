# ======================================================================
# tierb_common.R - shared set-up for the Tier B (independent 2026) scripts:
# configuration, 2026 metadata (roll offset removed), detections, the
# Calibration-fitting / Holdout split, and metric helpers.
# Sourced by tierb_calibration.R and tierb_effort.R (run from s/).
# ======================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(janitor)
  library(mgcv); library(ggplot2); library(purrr)
})
DET_ENGINE <- "gam"
source("./gamfunc.R")   # sources cvfunc.R

# ----------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------
DET_CSV   <- "../data/raw_model/2224scallop_yolo12n_slurm_598580/tierb2026/autotest_tierb2026_yolo12n.csv"
SPLIT_CSV <- "../data/raw/dataset_split_tierb2026.csv"
META_RDA  <- "../data/raw/allimg_2026.RData"
OUT_DIR   <- "../data/processed/tierb2026"
FIG_DIR   <- "../figures/tierb2026"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

N_FOLDS  <- 5
N_BOOT   <- 1000
SEED     <- 2026
ROLL_OFFSET_2026_DEG <- 46.12          # see claude/2026-fov-roll-offset-fix.md
MAIN_EMPTY_SHARE <- 0.124              # empties share of the 2022-2026 annotated set
N_MATCH_DRAWS    <- 200
MIN_POS_IMAGES_FIT <- 100              # min scallop images in Calibration-fitting to fit a region/platform

# Stage-1 forms: YOLOv12_formula_cv_2226.csv (chosen = TRUE)
S1_FORMS <- list(
  GB  = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) +
          te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + te(conf, bottom_depth, k = 5),
  MAB = truedetect ~ s(conf, k = 5)
)
S1_IDS <- c(GB = "M15_All_te", MAB = "M01_ConfOnly")

# Stage-2 YOLO-only candidates (same definitions as gammod.R)
S2_CANDIDATES <- list(
  M01_YoloRaw   = n_annotations ~ s(pred_yolo),
  M04_YoloLog   = n_annotations ~ s(log_yolo_pred),
  C01_Overlap   = n_annotations ~ s(log_yolo_pred) + s(overlap_frac, k = 5),
  C02_Spacing   = n_annotations ~ s(log_yolo_pred) + s(nn_ratio, k = 5),
  C03_LowConf   = n_annotations ~ s(log_yolo_pred) + s(lowconf_frac, k = 5),
  C04_BoxSize   = n_annotations ~ s(log_yolo_pred) + s(med_box_mm, k = 5),
  C05_teOverlap = n_annotations ~ te(log_yolo_pred, overlap_frac, k = c(5, 4)),
  C06_teBoxSize = n_annotations ~ te(log_yolo_pred, med_box_mm, k = c(5, 4)),
  C07_teLowConf = n_annotations ~ te(log_yolo_pred, lowconf_frac, k = c(5, 4)),
  C08_Structure = n_annotations ~ s(log_yolo_pred) + s(med_box_mm, k = 5) + s(lowconf_frac, k = 5) +
                    s(overlap_frac, k = 5) + s(nn_ratio, k = 5)
)
# Best YOLO-only stage-2 form in the main CV (synergy_cv_selection_2226.csv)
S2_TRANSFER <- c(GB = "C07_teLowConf", MAB = "C08_Structure")

if (!exists("LOG_NAME")) LOG_NAME <- "tierb_run_log.txt"
log_file <- file.path(OUT_DIR, LOG_NAME)
say <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE) }
cat("", file = log_file)
say("Tier B run (", LOG_NAME, ") ", format(Sys.time()))

# ----------------------------------------------------------------------
# 1. 2026 vehicle metadata (primary HabCam), roll offset removed
# ----------------------------------------------------------------------
split <- read.csv(SPLIT_CSV, stringsAsFactors = FALSE,
                  colClasses = c(name_cruise_token = "character", name_date_token = "character",
                                 name_time_token = "character")) %>%
  mutate(image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"),
         is_empty = as.logical(is_empty)) %>%
  distinct(image_id, .keep_all = TRUE)
say("Split images: ", nrow(split), " | ", paste(names(table(split$platform)), table(split$platform), collapse = ", "))

# Footprint (m^2) of a tilted camera over a flat seabed; same as FOV() in
# pred/allmeta1.R and s/datfunc.R (copied here so this script is standalone).
FOV <- function(altitude, roll, pitch) {
  DTOR <- pi / 180; focalLength <- 16 * 0.00133; PIXEL_SIZE <- 0.00000586
  sP <- sin(pitch * DTOR); cP <- cos(pitch * DTOR)
  sR <- sin(roll  * DTOR); cR <- cos(roll  * DTOR)
  m <- matrix(c(cP, 0, -sP,  sP * sR, cR, cP * sR,  sP * cR, -sR, cP * cR), nrow = 3, byrow = TRUE)
  ulX <- -1936 / 2 * PIXEL_SIZE; ulY <- 1216 / 2 * PIXEL_SIZE
  px <- c(ulX, ulX, -ulX, -ulX); py <- c(ulY, -ulY, -ulY, ulY)   # ul, ll, lr, ur
  X <- m[1, 1] * px + m[1, 2] * py + m[1, 3] * (-focalLength)
  Y <- m[2, 1] * px + m[2, 2] * py + m[2, 3] * (-focalLength)
  Z <- m[3, 1] * px + m[3, 2] * py + m[3, 3] * (-focalLength)
  X <- X * (altitude / Z); Y <- Y * (altitude / Z)
  j <- c(4, 1, 2, 3)
  -sum((X[j] + X) * (Y[j] - Y)) / 2
}
stopifnot(abs(FOV(2, 0, 0) - (1936 * 5.86e-6 * 2 / (16 * 0.00133)) * (1216 * 5.86e-6 * 2 / (16 * 0.00133))) < 1e-9)

load(META_RDA)   # allimg
allimg <- allimg %>%
  mutate(imagename = sub("\\.tif$", ".png", imagename),
         image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$")) %>%
  filter(image_id %in% split$image_id) %>%
  distinct(image_id, .keep_all = TRUE)
say("Split images with vehicle metadata: ", nrow(allimg))

meta26 <- allimg %>%
  rename(any_of(c(latitude = "Latitude", longitude = "Longitude",
                  vehicle_pitch_angle = "Pitch", vehicle_roll_angle = "Roll",
                  altitude = "Altitude", ctd_vehicle_depth_meter = "V_Depth",
                  ctd_temperature_celsius = "t", fluorometer_backscatter_ntu = "backscatter",
                  bottom_depth = "Bottom_Depth"))) %>%
  mutate(field_of_view_sq_meter = mapply(FOV, altitude, vehicle_roll_angle - ROLL_OFFSET_2026_DEG,
                                         vehicle_pitch_angle),
         millimeter_per_pixel = sqrt(field_of_view_sq_meter / (1936 * 1216)) * 1000) %>%
  select(image_id, latitude, longitude, altitude, bottom_depth, field_of_view_sq_meter,
         millimeter_per_pixel, ctd_temperature_celsius, fluorometer_backscatter_ntu,
         vehicle_roll_angle, vehicle_pitch_angle)
say(sprintf("FOV (m^2) median %.2f [IQR %.2f-%.2f]; mm/px median %.3f",
            median(meta26$field_of_view_sq_meter, na.rm = TRUE),
            quantile(meta26$field_of_view_sq_meter, .25, na.rm = TRUE),
            quantile(meta26$field_of_view_sq_meter, .75, na.rm = TRUE),
            median(meta26$millimeter_per_pixel, na.rm = TRUE)))

meta <- split %>%
  left_join(meta26, by = "image_id") %>%
  mutate(has_meta = !is.na(latitude) & !is.na(longitude),
         region = case_when(!has_meta ~ NA_character_, longitude >= -71 ~ "GB", TRUE ~ "MAB"))

# Sentinel values (<= -90 = sensor missing): NA, then cruise / region median (covariates only)
SENTINEL_VARS <- c("ctd_temperature_celsius", "fluorometer_backscatter_ntu", "bottom_depth",
                   "altitude", "field_of_view_sq_meter")
meta <- meta %>%
  mutate(cov_imputed = has_meta & if_any(all_of(SENTINEL_VARS), ~ is.na(.x) | .x <= -90),
         across(all_of(SENTINEL_VARS), ~ if_else(!is.na(.x) & .x <= -90, NA_real_, .x))) %>%
  group_by(name_cruise_token, region) %>%
  mutate(across(all_of(SENTINEL_VARS), ~ if_else(has_meta, coalesce(.x, median(.x, na.rm = TRUE)), .x))) %>%
  ungroup()

cov_tab <- meta %>% count(platform, has_meta, region, partition)
say("Metadata coverage (platform x region x partition):")
print(cov_tab); write.csv(cov_tab, file.path(OUT_DIR, "metadata_coverage.csv"), row.names = FALSE)
say("Images with an imputed covariate: ", sum(meta$cov_imputed, na.rm = TRUE))
n_prim_nometa <- sum(meta$platform == "primary" & !meta$has_meta)
if (n_prim_nometa > 0) say("!! ", n_prim_nometa, " primary images have no vehicle metadata; they are left out of the primary analysis.")

# ----------------------------------------------------------------------
# 2. Detections
# ----------------------------------------------------------------------
det_all <- read.csv(DET_CSV) %>%
  clean_names() %>%
  mutate(image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"),
         truedetect = tolower(as.character(truedetect)) == "true",
         spname = tolower(spname)) %>%
  filter(spname == "scallop", image_id %in% meta$image_id) %>%
  select(-imagename)
say("Detections (conf >= 0.1): ", nrow(det_all), " | true positives: ", sum(det_all$truedetect))

# Check: matched TPs never exceed annotations
chk <- det_all %>% group_by(image_id) %>% summarise(tp = sum(truedetect), .groups = "drop") %>%
  left_join(meta %>% select(image_id, n_annotations), by = "image_id")
if (any(chk$tp > chk$n_annotations)) say("!! ", sum(chk$tp > chk$n_annotations), " images have more TPs than annotations (check manifest vs eval labels).")

build_tables <- function(meta_df) {
  ids <- meta_df$image_id
  d <- det_all %>% filter(image_id %in% ids) %>% left_join(meta_df, by = "image_id")
  img <- meta_df %>%
    left_join(d %>% count(image_id, name = "n_auto"), by = "image_id") %>%
    mutate(n_auto = replace_na(n_auto, 0L))
  list(det = d, img = img)
}

# F1-optimal confidence threshold on Calibration-fitting (TP/FP from matching;
# FN = annotations not matched)
f1_threshold <- function(det, img) {
  n_gt <- sum(img$n_annotations)
  grid <- seq(0.10, 0.95, by = 0.01)
  tab <- map_dfr(grid, function(t) {
    k <- det$conf >= t; tp <- sum(det$truedetect & k); fp <- sum(!det$truedetect & k)
    p <- tp / max(tp + fp, 1); r <- tp / max(n_gt, 1)
    tibble(conf = t, precision = p, recall = r, f1 = ifelse(p + r > 0, 2 * p * r / (p + r), 0))
  })
  tab %>% filter(f1 == max(f1)) %>% slice_max(conf, n = 1)
}

# ----------------------------------------------------------------------
# 3. Metrics helpers
# ----------------------------------------------------------------------
metric_row <- function(y, yhat) {
  tibble(rmse = sqrt(mean((y - yhat)^2)), mae = mean(abs(y - yhat)),
         r2_1to1 = r2_1to1(y, yhat),
         pct_bias_total = 100 * (sum(yhat) - sum(y)) / sum(y),
         calibration_slope = unname(coef(lm(y ~ yhat))[2]),
         total_true = sum(y), total_pred = sum(yhat), n_images = length(y))
}
method_metrics <- function(df, methods, by = "region") {
  df %>% group_by(across(all_of(by))) %>%
    group_modify(function(d, k) map_dfr(names(methods), function(m)
      metric_row(d$n_annotations, d[[methods[[m]]]]) %>% mutate(method = m, .before = 1))) %>%
    ungroup()
}
matched_empty_metrics <- function(df, methods, share = MAIN_EMPTY_SHARE, draws = N_MATCH_DRAWS, by = "region") {
  set.seed(SEED)
  map_dfr(split(df, df[[by]]), function(d) {
    pos <- which(d$n_annotations > 0); emp <- which(d$n_annotations == 0)
    n_emp <- min(length(emp), round(share / (1 - share) * length(pos)))
    map_dfr(seq_len(draws), function(b) {
      dd <- d[c(pos, sample(emp, n_emp)), ]
      map_dfr(names(methods), function(m) metric_row(dd$n_annotations, dd[[methods[[m]]]]) %>%
                mutate(method = m, draw = b))
    }) %>%
      group_by(method) %>%
      summarise(across(c(rmse, mae, r2_1to1, pct_bias_total, calibration_slope),
                       list(med = median, lo = ~ quantile(.x, .025), hi = ~ quantile(.x, .975))),
                n_images = first(n_images), .groups = "drop") %>%
      mutate(!!by := d[[by]][1], empty_share = share, .before = 1)
  })
}
density_table <- function(df, methods, by = "region") {
  df %>% filter(!is.na(field_of_view_sq_meter)) %>%
    group_by(across(all_of(by))) %>%
    group_modify(function(d, k) {
      a <- sum(d$field_of_view_sq_meter)
      bind_rows(tibble(method = "Manual annotation", density_per_m2 = sum(d$n_annotations) / a),
                map_dfr(names(methods), function(m) tibble(method = m, density_per_m2 = sum(d[[methods[[m]]]]) / a))) %>%
        mutate(area_m2 = a, pct_diff_vs_manual = 100 * (density_per_m2 / density_per_m2[1] - 1))
    }) %>% ungroup()
}
haversine_m <- function(lat1, lon1, lat2, lon2) {
  r <- 6371000; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 + cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}

