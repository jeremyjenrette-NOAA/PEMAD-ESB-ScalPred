library(dplyr)
library(tidyr)
library(janitor)
library(stringr)

# ======================================================================
# 1. Master Metadata Prep
# ======================================================================
meta <- read.csv("../data/raw/dataset_split_star.csv") %>%
  clean_names() %>%
  rename(n_annotations = total_annotations) %>% 
  mutate(
    image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"), 
    region = if_else(longitude >= -71, "GB", "MAB") %>% factor(levels = c("MAB", "GB")),
    
    # --- Robust Boolean Coercion for Python/Pandas Exported Strings ---
    is_empty = as.logical(as.character(is_empty)),
    is_train = as.logical(as.character(is_train)),
    is_test  = as.logical(as.character(is_test))
  ) %>%
  distinct(image_id, .keep_all = TRUE)

# ======================================================================
# 2. Unified Processing Function
# ======================================================================
process_model <- function(det_csv_path, meta_df, dets_df) {
  
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
      # --- Robust Boolean Coercion for YOLO/Cascade Predictions ---
      truedetect = as.logical(as.character(truedetect)),
      spname = tolower(spname)
    ) %>%
    # Filter for any of our three target crab species (using valid R syntax)
    filter(spname %in% c("jonah_crab", "rock_crab", "cancer_sp"))
  
  # --- Isolate only the images used during evaluation ---
  test_meta_df <- meta_df %>% 
    filter(is_test == TRUE) %>%
    select(-imagename)
  
  # Compute overall automated detection counts per image
  auto_total <- at %>% count(image_id, name = "n_auto")
  
  # Compute species-specific automated detection counts per image
  auto_species <- at %>%
    mutate(species_cat = paste0("n_auto_", spname)) %>%
    count(image_id, species_cat) %>%
    pivot_wider(names_from = species_cat, values_from = n, values_fill = 0L)
  
  # Structural safeguard: Ensure all target auto columns exist even if unpredicted
  auto_cols <- c("n_auto_jonah_crab", "n_auto_rock_crab", "n_auto_cancer_sp")
  for (col in auto_cols) {
    if (!col %in% names(auto_species)) auto_species[[col]] <- 0L
  }
  
  # Combine counts into the comprehensive image-level summary
  img_df <- test_meta_df %>%
    left_join(auto_total, by = "image_id") %>%
    left_join(auto_species, by = "image_id") %>%
    mutate(
      n_auto = replace_na(n_auto, 0L),
      across(all_of(auto_cols), ~ replace_na(.x, 0L)),
      man_density = n_annotations / field_of_view_sq_meter,
      auto_density = n_auto / field_of_view_sq_meter
    )
  
  # Join metadata directly to the detection-level frame
  det_df <- at %>%
    left_join(test_meta_df, by = "image_id")
  
  return(list(img = img_df, det = det_df))
}

# ======================================================================
# 3. Process Models and Bundle
# ======================================================================
print("Processing YOLO...")
yolo_data <- process_model("../data/raw/crab_eval_yolov12_multi/autotest2426_yolo12n.csv", meta, dets)

print("Processing Cascade R-CNN...")
cas_data <- process_model("../data/raw/crab_eval_cascade_multi/autotest2426_viame_cascade.csv", meta, dets)

<<<<<<< HEAD
# Combine into cohesive nested analytical lists
=======
>>>>>>> 099e1fe8751b6dd1e92581670bbe7f7234962132
model_results <- list(
  `YOLOv12` = yolo_data,
  `Cascade R-CNN` = cas_data
)

saveRDS(model_results, file = "../data/processed/crab_evalmulti_2426.rds")
print("Success! Multi-class data structures successfully unified with standardized logicals.")
