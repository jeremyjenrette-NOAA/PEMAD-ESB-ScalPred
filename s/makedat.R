library(dplyr)
library(tidyr)
library(janitor)
library(stringr)

# ======================================================================
# 1. Processing Function for Seal Detections
# ======================================================================
process_model <- function(det_csv_path) {
  
  raw_at <- read.csv(det_csv_path, stringsAsFactors = FALSE)
  
  # Force coordinate names to lowercase to bypass janitor acronym rules
  names(raw_at)[tolower(names(raw_at)) == "tlx"] <- "tlx"
  names(raw_at)[tolower(names(raw_at)) == "tly"] <- "tly"
  names(raw_at)[tolower(names(raw_at)) == "brx"] <- "brx"
  names(raw_at)[tolower(names(raw_at)) == "bry"] <- "bry"
  
  at <- raw_at %>%
    clean_names() %>%
    mutate(
      # Standardize boolean columns for R (converts "False"/"True" strings to FALSE/TRUE)
      truedetect = as.logical(toupper(as.character(truedetect))),
      
      # Normalize tiled image names (removes "_tile_X_Y" pattern from YOLO files)
      imagename_clean = str_remove(imagename, "_tile_\\d+_\\d+"),
      image_id = str_remove(imagename_clean, "\\.[A-Za-z0-9]+$"),
      
      # Standardize species labels
      spname = tolower(trimws(spname))
    ) %>%
    # Filter strictly for target seal classes
    filter(spname %in% c("adult", "pup"))
  
  # Compute overall automated detection counts per full image
  auto_total <- at %>% 
    count(image_id, name = "n_auto")
  
  # Compute species-specific automated detection counts per image
  auto_species <- at %>%
    mutate(species_cat = paste0("n_auto_", spname)) %>%
    count(image_id, species_cat) %>%
    pivot_wider(names_from = species_cat, values_from = n, values_fill = 0L)
  
  # Structural safeguard: Ensure target species columns exist even if zero detections
  auto_cols <- c("n_auto_adult", "n_auto_pup")
  for (col in auto_cols) {
    if (!col %in% names(auto_species)) auto_species[[col]] <- 0L
  }
  
  # Create image-level summary dataframe
  img_df <- auto_total %>%
    left_join(auto_species, by = "image_id")
  
  return(list(img = img_df, det = at))
}

# ======================================================================
# 2. Process Models and Save Output
# ======================================================================
print("Processing YOLOv12...")
yolo_data <- process_model("../data/raw/seal_eval_yolo/autotest26_yolo12n.csv")

print("Processing Cascade R-CNN...")
cas_data <- process_model("../data/raw/seal_eval_cascade/autotest26_viame_cascade.csv")

# Combine into cohesive nested analytical list structure
model_results <- list(
  `YOLOv12` = yolo_data,
  `Cascade R-CNN` = cas_data
)

# Save processed data structures
saveRDS(model_results, file = "../data/processed/seal_eval_2026.rds")
print("Success! Seal detection datasets successfully normalized and standardized.")
