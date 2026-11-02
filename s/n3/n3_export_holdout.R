# ======================================================================
# n3_export_holdout.R - one Holdout test table with every workflow's prediction
# for error.R and fig_holdout.R (replaces synergy_holdout_eval_2226.rds).
#
# Predictions come from n3_ablation.R (primary rule: lowest CV RMSE);
# image metadata (covariates, position, cruise, crowding) from the table gammod.R saved.
# Run from s/:  Rscript n3/n3_export_holdout.R
# Writes ../data/processed/holdout_eval_2226_v2.rds
# ======================================================================
suppressPackageStartupMessages({ library(dplyr) })
n3  <- readRDS("../data/processed/n3/n3_results.rds")
old <- readRDS("../data/processed/synergy_holdout_eval_2226.rds") %>% filter(dataset == "test_GAM_test")

PRED <- c(pred_2s_yolo = "Two-stage, one detector: YOLOv12",
          pred_2s_yc   = "Two-stage, two detectors: YOLOv12 + Cascade R-CNN",
          pred_yolo    = "pred_Y",  pred_cascade = "pred_C",
          yolo_f1      = "f1_Y",    cascade_f1   = "f1_C")
h <- n3$hold %>% filter(dataset == "test_GAM_test")
stopifnot(all(PRED %in% names(h)))
new <- h %>% select(image_id, n_annotations_n3 = n_annotations, all_of(PRED))

# crowding features from YOLO boxes, if gammod's table lacks them
CROWD <- c("n_det_all", "nn_ratio", "overlap_frac", "lowconf_frac", "med_box_mm")
miss  <- setdiff(CROWD, names(old))
if (length(miss)) new <- new %>% left_join(h %>% select(image_id, any_of(paste0(miss, "_Y"))) %>%
                                             rename_with(~ sub("_Y$", "", .x), -image_id), by = "image_id")

ev <- old %>% select(-any_of(c(names(PRED), "pred_synergy"))) %>%
  inner_join(new, by = "image_id")
stopifnot(nrow(ev) == nrow(old), nrow(ev) == nrow(h),
          isTRUE(all.equal(ev$n_annotations, ev$n_annotations_n3)))
ev <- ev %>% select(-n_annotations_n3)
saveRDS(ev, "../data/processed/holdout_eval_2226_v2.rds")
cat("Holdout images:", nrow(ev), "| scallops:", sum(ev$n_annotations), "\n")
print(ev %>% group_by(region) %>%
        summarise(n_img = n(), manual = sum(n_annotations),
                  across(all_of(names(PRED)), ~ round(100 * (sum(.x) - manual) / manual, 1),
                         .names = "bias_{.col}")), width = Inf)
