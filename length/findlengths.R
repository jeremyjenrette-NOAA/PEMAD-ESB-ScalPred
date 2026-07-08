at_yolo = read.csv("../data/raw/model_output/2224scallop_yolo12n_261070/eval/autotest2224_yolo12n.csv") 
# mt_yolo = read.csv("../data/raw/2224scallop_yolo12n_261070/eval/mantest2224_yolo12n.csv")
yr2224 = rbind(
  read.csv("../../PEMAD-ESB-ScalTrain/data/raw/annotations_2022.csv"),
  read.csv("../../PEMAD-ESB-ScalTrain/data/raw/annotations_2023.csv"),
  read.csv("../../PEMAD-ESB-ScalTrain/data/raw/annotations_2024.csv")
)

meta = read.csv("../data/processed/metadata2224.csv")

at_yolo <- at_yolo %>%
  mutate(truedetect = if_else(truedetect == "True", TRUE, FALSE)) %>%
  filter(truedetect == TRUE)

###############################################################################

library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(png)
library(jpeg)
library(ggplot2)
library(grid)
library(magick)
library(clue)

source("./lengthfunc.R")

# ----------------------------
# build sourced GT table
# ----------------------------

yr_gt_src <- yr2224 %>%
  filter(geom_type == "line") %>%
  mutate(
    year = get_year(imagename),
    annotated_side = get_annotated_side(year),
    image_path = map2_chr(imagename, year, find_image_path)
  ) %>%
  filter(!is.na(image_path)) %>%
  mutate(
    parsed = map(geometry_text, parse_geometry_text),
    x1 = map_dbl(parsed, 1),
    y1 = map_dbl(parsed, 2),
    x2 = map_dbl(parsed, 3),
    y2 = map_dbl(parsed, 4)
  ) %>%
  select(-parsed)

# attach true image dimensions
dims_tbl <- yr_gt_src %>%
  distinct(image_path) %>%
  mutate(dims = map(image_path, read_img_dims_safe)) %>%
  unnest(dims)

yr_gt_src <- yr_gt_src %>%
  rowwise() %>%
  mutate(
    img_date = extract_img_date(imagename)
  ) %>%
  ungroup()

yr_gt_src <- yr_gt_src %>%
  left_join(dims_tbl, by = "image_path") %>%
  filter(img_readable, !is.na(img_width), !is.na(img_height)) %>%
  mutate(
    half_width = img_width / 2,
    
    # first remap from stereo -> split image coordinates
    x1_split_raw = case_when(
      annotated_side == "left"  ~ x1,
      annotated_side == "right" ~ x1 - half_width,
      TRUE ~ NA_real_
    ),
    x2_split_raw = case_when(
      annotated_side == "left"  ~ x2,
      annotated_side == "right" ~ x2 - half_width,
      TRUE ~ NA_real_
    ),
    
    y1_split = y1,
    y2_split = y2,
    
    # conditional calibration fix
    shift_px_applied = get_shift_px(year, img_date),
    
    x1_split = x1_split_raw + shift_px_applied,
    x2_split = x2_split_raw + shift_px_applied,
    
    split_side = annotated_side,
    split_width = half_width,
    split_height = img_height
  ) %>%
  mutate(
    # optional clipping after shift
    x1_split = pmin(pmax(x1_split, 0), split_width),
    x2_split = pmin(pmax(x2_split, 0), split_width)
  ) %>%
  mutate(
    gt_length_px = sqrt((x2_split - x1_split)^2 + (y2_split - y1_split)^2),
    gt_box = pmap(
      list(x1_split, y1_split, x2_split, y2_split),
      ~ lineseg2bb(..1, ..2, ..3, ..4)
    ),
    gt_TLx = map_dbl(gt_box, 1),
    gt_TLy = map_dbl(gt_box, 2),
    gt_BRx = map_dbl(gt_box, 3),
    gt_BRy = map_dbl(gt_box, 4)
  ) %>%
  select(-gt_box) %>%
  group_by(imagename) %>%
  mutate(man = row_number()) %>%
  ungroup() %>%
  filter(deprecated == "f")

dims_tbl %>%
  count(img_readable)

dims_tbl %>%
  filter(!img_readable) %>%
  head(20)

saveRDS(dims_tbl, "../data/raw/image_dimensions.rds")

# ----------------------------
# plot images with line annotations
# ----------------------------

plot_stereo_annotation("202403.20240514.225143959.8916.png")
plot_split_annotation("202403.20240514.225143959.8916.png")

plot_stereo_annotation("202404.20240627.040203881.4944.png")
plot_split_annotation("202404.20240627.040203881.4944.png")

# ----------------------------
# match annotations
# ----------------------------

at_tp <- at_yolo %>%
  filter(truedetect == TRUE) %>%
  mutate(
    pred_width_px  = BRx - TLx,
    pred_height_px = BRy - TLy,
    pred_length_px = (pred_width_px + pred_height_px) / 2
  ) %>%
  select(-c(man))

# restrict GT to the same valid images if needed
valid_images <- intersect(unique(at_tp$Imagename), unique(yr_gt_src$imagename))

at_tp2 <- at_tp %>%
  filter(Imagename %in% valid_images)

yr_gt2 <- yr_gt_src %>%
  filter(imagename %in% valid_images)

matched_pairs <- map_dfr(valid_images, function(img) {
  tp_sub <- at_tp2 %>% filter(Imagename == img)
  gt_sub <- yr_gt2 %>% filter(imagename == img)
  
  match_tp_to_gt_one_image(tp_sub, gt_sub, iou_tol = 1e-2)
})

meta2 <- meta %>%
  rename(Imagename = IMAGENAME)

matched_len_df <- at_tp2 %>%
  left_join(
    matched_pairs %>% select(Imagename, Detectid, man, iu_recorded, iu_recomputed, iu_abs_diff),
    by = c("Imagename", "Detectid")
  ) %>%
  left_join(
    yr_gt2,
    by = c("Imagename" = "imagename", "man" = "man")
  ) %>%
  left_join(meta2, by = "Imagename") %>%
  mutate(
    gt_length_mm   = gt_length_px * MILLIMETER_PER_PIXEL,
    pred_length_mm = pred_length_px * MILLIMETER_PER_PIXEL
  )


c(
  n_tp = nrow(at_tp2),
  n_matched = nrow(matched_pairs),
  prop_matched = nrow(matched_pairs) / nrow(at_tp2)
)

summary(matched_pairs$iu_abs_diff)

tp_counts <- at_tp2 %>%
  count(Imagename, name = "n_tp")

gt_counts <- yr_gt2 %>%
  count(imagename, name = "n_gt") %>%
  rename(Imagename = imagename)

tp_counts %>%
  left_join(gt_counts, by = "Imagename") %>%
  mutate(tp_gt_diff = n_tp - n_gt) %>%
  arrange(desc(tp_gt_diff)) %>%
  head(20)

ggplot(matched_pairs, aes(x = iu_recorded, y = iu_recomputed)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    title = "Recorded vs recomputed IoU",
    x = "Recorded iu",
    y = "Recomputed IoU"
  )

# ----------------------------
# Measuring error
# ----------------------------

matched_len_df <- matched_len_df %>%
  mutate(
    width_px  = BRx - TLx,
    height_px = BRy - TLy,
    
    # current method
    len_mean_px = (width_px + height_px) / 2,
    
    # alternative 1: max dimension (robust to orientation)
    len_max_px = pmax(width_px, height_px),
    
    # alternative 2: min dimension (probably bad, but test)
    len_min_px = pmin(width_px, height_px),
    
    # alternative 3: diagonal (overestimate)
    len_diag_px = sqrt(width_px^2 + height_px^2),
    
    # convert to mm
    len_mean_mm = len_mean_px * MILLIMETER_PER_PIXEL,
    len_max_mm  = len_max_px  * MILLIMETER_PER_PIXEL,
    len_min_mm  = len_min_px  * MILLIMETER_PER_PIXEL,
    len_diag_mm = len_diag_px * MILLIMETER_PER_PIXEL
  ) %>%
  filter(!is.na(iu_recorded))

results <- bind_rows(
  evaluate_length_method(matched_len_df, "len_mean_mm"),
  evaluate_length_method(matched_len_df, "len_max_mm"),
  evaluate_length_method(matched_len_df, "len_min_mm"),
  evaluate_length_method(matched_len_df, "len_diag_mm")
)

results

p_lengthfit = ggplot(matched_len_df, aes(x = gt_length_mm)) +
  geom_point(aes(y = len_mean_mm), alpha = 0.2, color = "black") +
  # geom_point(aes(y = len_max_mm), alpha = 0.2, color = "green") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    x = "True length (mm)",
    y = "Predicted length (mm)",
    title = "Comparison of length estimation methods"
  ) +
  theme_minimal()


p_lengtherror = matched_len_df %>%
  mutate(
    err_mean = len_mean_mm - gt_length_mm
  ) %>%
  pivot_longer(cols = c(err_mean),
               names_to = "method",
               values_to = "error") %>%
  ggplot(aes(x = gt_length_mm, y = error, color = method)) +
  geom_point(alpha = 0.2, col = "black") +
  geom_smooth(se = TRUE, col = "red") +
  # geom_line(aes(y=0), col = "red", lty = 2, lwd = 1) +
  labs(
    x = "True length (mm)",
    y = "Residual error (mm)",
    title = "Bias vs scallop size"
  ) +
  theme_minimal()

p_lengtherr_hist = matched_len_df %>%
  mutate(
    err_mean = len_mean_mm - gt_length_mm
  ) %>%
  pivot_longer(cols = c(err_mean),
               names_to = "method",
               values_to = "error") %>%
  ggplot(aes(x = error, fill = method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
  labs(
    x = "Error (mm)",
    y = "Count",
    title = "Distribution of length errors"
  ) +
  theme(
    legend.position = "none"
  )

img_name = "202203.20220517.084648378.99571.png"
p1 = plot_detection_length_debug(
  imagename = img_name,
  matched_df = matched_len_df
)
ggsave(p1, file = paste0("../data/images/",img_name), width = 7, height = 7)


ggsave(p_lengthfit, file = "../figures/2224_lengthfit.png", width = 6, height = 6)
ggsave(p_lengtherror, file = "../figures/2224_lengtherror.png", width = 6, height = 6)
# ggsave(p_lengtherr_hist, file = "../figures/2224_lengtherr_hist.png", width = 6, height = 6)
# ----------------------------
# Predicting error
# ----------------------------

cal_model <- lm(gt_length_mm ~ len_mean_mm, data = matched_len_df)

summary(cal_model)

matched_len_df <- matched_len_df %>%
  mutate(
    len_calibrated_mm = predict(cal_model)
  )

evaluate_length_method(matched_len_df, "len_calibrated_mm")

ggplot(matched_len_df, aes(x = gt_length_mm, y = len_calibrated_mm)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    title = "Calibrated length estimation",
    x = "True length (mm)",
    y = "Calibrated prediction (mm)"
  ) +
  theme_minimal()
