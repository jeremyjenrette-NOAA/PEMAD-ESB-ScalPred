library(dplyr)
library(janitor)
library(stringr)
library(tidyr)
source("./cvfunc.R")

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

# Images with confirmed annotation errors (e.g. labelled empty but full of
# scallops). They are REMOVED, not relabelled: their detections carry
# truedetect labels from matching against the wrong ground truth, so fixing
# n_annotations alone would still feed wrong labels to the detection GAM.
# Add rows to the CSV (status = "confirmed") after visual review; rows with any
# other status are reported but kept. Re-annotating an image on the cluster and
# re-running the detector evaluation is the alternative to exclusion.
excl <- read.csv("../data/raw/annotation_exclusions_2226.csv", stringsAsFactors = FALSE)
excl_ids <- str_remove(trimws(excl$imagename[tolower(trimws(excl$status)) == "confirmed"]),
                       "\\.[A-Za-z0-9]+$")
cat("Annotation exclusions (confirmed):", length(excl_ids), "\n")
print(meta %>% filter(image_id %in% excl_ids) %>% count(dataset, region))
meta <- meta %>% filter(!image_id %in% excl_ids)

# Sentinel covariate values (-99.99 = sensor missing; MAB AUV cruise 202308 has
# them for CTD temperature and fluorometer backscatter). Left in place they
# enter s(ctd_temperature_celsius) etc. as real values near -100. They are set
# to NA and imputed with the cruise median (region median if the whole cruise
# is missing). Imputation uses covariates only, never labels. cov_imputed flags
# affected images.
SENTINEL_VARS <- c("ctd_temperature_celsius", "fluorometer_backscatter_ntu", "bottom_depth",
                   "altitude", "field_of_view_sq_meter")
meta <- meta %>%
  mutate(cov_imputed = if_any(any_of(SENTINEL_VARS), ~ !is.na(.x) & .x <= -90),
         across(any_of(SENTINEL_VARS), ~ if_else(.x <= -90, NA_real_, .x))) %>%
  group_by(name_cruise_token) %>%
  mutate(across(any_of(SENTINEL_VARS), ~ coalesce(.x, median(.x, na.rm = TRUE)))) %>%
  group_by(region) %>%
  mutate(across(any_of(SENTINEL_VARS), ~ coalesce(.x, median(.x, na.rm = TRUE)))) %>%
  ungroup()
cat("Images with imputed sentinel covariates:", sum(meta$cov_imputed), "\n")

# Cross-fitting folds for the detection-GAM stage (see cvfunc.R).
# Whole transect blocks, assigned within region, high-density blocks balanced.
# Assigned ONCE here so every detector shares the same folds.
N_FOLDS <- 5
meta <- assign_folds(meta, K = N_FOLDS, seed = 2226)

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
    filter(spname == "scallop", !image_id %in% excl_ids)
  
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
yolo_data <- process_model("../data/raw_model/2226data_yolo12n_slurm_598585/autotest2226_yolo12n.csv", meta)

print("Processing Cascade R-CNN...")
cas_data <- process_model("../data/raw_model/2226data_cascade_slurm_598695/autotest2226_viame_cascade.csv", meta)

model_results <- list(
  `YOLOv12` = yolo_data,
  `Cascade R-CNN` = cas_data
)

saveRDS(model_results, file = "../data/processed/eval_2226.rds")
print("Success! Unified data saved.")
