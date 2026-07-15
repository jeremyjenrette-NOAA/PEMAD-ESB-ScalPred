library(dplyr)
library(tidyr)
library(janitor)
library(stringr)

# ======================================================================
# 1. Master Metadata Prep
# ======================================================================
meta <- read.csv("../data/raw/dataset_split_crab.csv") %>%
  clean_names() %>%
  rename(n_annotations = total_annotations) %>% 
  mutate(
    image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"), 
    region = if_else(longitude >= -71, "GB", "MAB") %>% factor(levels = c("MAB", "GB")),
    is_test = as.logical(is_test)
  ) %>%
  distinct(image_id, .keep_all = TRUE)

# ======================================================================
# 2. Unified Processing Function
# ======================================================================
process_model <- function(det_csv_path, meta_df, dets_df) {
  
  # Ensure ground-truth names match expected lowercase formats
  dets_clean <- dets_df %>% clean_names()
  
  raw_at <- read.csv(det_csv_path)
  
  # Force coordinate names to lowercase to bypass janitor acronym rules
  names(raw_at)[tolower(names(raw_at)) == "tlx"] <- "tlx"
  names(raw_at)[tolower(names(raw_at)) == "tly"] <- "tly"
  names(raw_at)[tolower(names(raw_at)) == "brx"] <- "brx"
  names(raw_at)[tolower(names(raw_at)) == "bry"] <- "bry"
  
  at <- raw_at %>%
    clean_names() %>%
    mutate(
      image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"),
      truedetect = if_else(tolower(truedetect) == "true", TRUE, FALSE),
      spname = tolower(spname)
    ) %>%
    filter(spname == "crab") 
  
  # --- Match True Positives to 'dets' via IoU to recover original labels ---
  tp_dets <- at %>% filter(truedetect == TRUE)
  
  if (nrow(tp_dets) > 0) {
    tp_matched <- tp_dets %>%
      inner_join(
        dets_clean, 
        by = "imagename", 
        suffix = c("_det", "_gt"), 
        relationship = "many-to-many"
      ) %>%
      mutate(
        xa = pmax(tlx_det, tlx_gt),
        ya = pmax(tly_det, tly_gt),
        xb = pmin(brx_det, brx_gt),
        yb = pmin(bry_det, bry_gt),
        
        inter_area = pmax(0, xb - xa) * pmax(0, yb - ya),
        det_area = (brx_det - tlx_det) * (bry_det - tly_det),
        gt_area = (brx_gt - tlx_gt) * (bry_gt - tly_gt),
        
        iou = inter_area / (det_area + gt_area - inter_area)
      ) %>%
      group_by(detectid) %>%
      slice_max(iou, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(detectid, original_label = label)
    
    at <- at %>% left_join(tp_matched, by = "detectid")
  } else {
    at <- at %>% mutate(original_label = NA_character_)
  }
  
  # --- NEW: Compute ground truth species counts per image ---
  gt_counts <- dets_clean %>%
    mutate(
      image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"),
      species_cat = case_when(
        label == "rock_crab" ~ "n_rock_crab",
        label == "jonah_crab" ~ "n_jonah_crab",
        label %in% c("cancer_sp.", "cancer_sp") ~ "n_cancer_sp",
        TRUE ~ "n_other"
      )
    ) %>%
    count(image_id, species_cat) %>%
    pivot_wider(names_from = species_cat, values_from = n, values_fill = 0L)
  
  # Safely ensure all 4 target columns exist even if a category is missing entirely
  target_cols <- c("n_rock_crab", "n_jonah_crab", "n_cancer_sp", "n_other")
  for (col in target_cols) {
    if (!col %in% names(gt_counts)) gt_counts[[col]] <- 0L
  }
  
  # --- Isolate only the images used during evaluation ---
  test_meta_df <- meta_df %>% filter(is_test == TRUE) %>%
    select(-imagename)
  
  auto_counts <- at %>% count(image_id, name = "n_auto")
  
  # Combine counts into the image-level summary
  img_df <- test_meta_df %>%
    left_join(auto_counts, by = "image_id") %>%
    left_join(gt_counts, by = "image_id") %>%
    mutate(
      n_auto = replace_na(n_auto, 0L),
      # Efficiently replace NAs with 0 across all species count columns
      across(all_of(target_cols), ~ replace_na(.x, 0L)),
      man_density = n_annotations / field_of_view_sq_meter,
      auto_density = n_auto / field_of_view_sq_meter
    )
  
  det_df <- at %>%
    left_join(test_meta_df, by = "image_id")
  
  return(list(img = img_df, det = det_df))
}

# ======================================================================
# 3. Process Models and Bundle
# ======================================================================
print("Processing YOLO...")
yolo_data <- process_model("../data/raw/crab_eval_yolov12/autotest.csv", meta, dets)

print("Processing Cascade R-CNN...")
cas_data <- process_model("../data/raw/crab_eval_cascade/autotest2426_viame_cascade.csv", meta, dets)

model_results <- list(
  `YOLOv12` = yolo_data,
  `Cascade R-CNN` = cas_data
)

saveRDS(model_results, file = "../data/processed/crab_eval_2426.rds")
print("Success! Unified data saved with image-level species profiles.")