library(dplyr)
library(janitor)
library(stringr)

# ======================================================================
# 1. Master Metadata Prep
# ======================================================================
meta <- read.csv("../data/raw/dataset_split_2226.csv") %>%
  clean_names() %>%
  mutate(
    image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"), 
    region = if_else(longitude >= -71, "GB", "MAB") %>% factor(levels = c("MAB", "GB")),
    # Robustly handle the boolean just in case pandas exported it as a string
    is_test = as.logical(is_test)
  ) %>%
  distinct(image_id, .keep_all = TRUE)

# ======================================================================
# 2. Unified Processing Function
# ======================================================================
process_model <- function(det_csv_path, meta_df) {
  
  at <- read.csv(det_csv_path) %>%
    clean_names() %>%
    mutate(
      image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"),
      truedetect = if_else(tolower(truedetect) == "true", TRUE, FALSE),
      spname = tolower(spname)
    ) %>%
    filter(spname == "scallop")
  
  # --- THE FIX: Isolate only the images used during evaluation ---
  test_meta_df <- meta_df %>% filter(is_test == TRUE) %>%
    select(-imagename)
  
  auto_counts <- at %>% count(image_id, name = "n_auto")
  
  img_df <- test_meta_df %>%
    left_join(auto_counts, by = "image_id") %>%
    mutate(
      n_auto = replace_na(n_auto, 0L),
      man_density = n_annotations / field_of_view_sq_meter,
      auto_density = n_auto / field_of_view_sq_meter
    )
  
  det_df <- at %>%
    left_join(
      test_meta_df,
      by = "image_id"
    )
  
  return(list(img = img_df, det = det_df))
}

# ======================================================================
# 3. Process Models and Bundle
# ======================================================================
print("Processing YOLO...")
yolo_data <- process_model("../data/raw/2226scallop_yolo12n_441288/eval/autotest2226_yolo12n.csv", meta)

print("Processing Cascade R-CNN...")
cas_data <- process_model("../data/raw/2226scallop_viame_441305/autotest2226_viame_cascade.csv", meta)

# Rename them right here:
model_results <- list(
  `YOLOv12` = yolo_data,
  `Cascade R-CNN` = cas_data
)

saveRDS(model_results, file = "../data/processed/eval_2226.rds")
print("Success! Unified data saved.")