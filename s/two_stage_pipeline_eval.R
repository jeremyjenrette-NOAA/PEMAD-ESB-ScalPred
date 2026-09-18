#!/usr/bin/env Rscript
# =============================================================================
# two_stage_pipeline_eval.R
#
# Performance exploration & status dashboard for the two-stage marine species
# pipeline (Stage 1: YOLO broad detector  ->  Stage 2: hierarchical taxonomic
# classifier). Built around the star24 (starfish) autotest run, but written
# generically enough to point at another taxon (e.g. crabdata2426) as long as
# the same file/column conventions are used.
#
# WHAT THIS SCRIPT DOES
#   1. Stage 1 (detector)   : precision/recall/AP for the broad "is there an
#                              animal here" detector, IoU behavior, false
#                              alarms on background-only images.
#   2. Stage 2 (classifier) : confusion matrices (species + genus), per-class
#                              precision/recall/F1, genus-only "restraint"
#                              behavior, confidence calibration.
#   3. End-to-end cascade   : per-species / per-genus AP using the classifier's
#                              own probability columns as the ranking score,
#                              rolled up into a single mAP headline number.
#   4. Workflow status      : known pipeline fixes, dataset composition, and
#                              a compact scorecard of all headline numbers.
#
# INPUT FILES (edit DATA_DIR below, or pass as first CLI arg)
#   autotest_two_stage_cascade.csv  - per-detection Stage1+Stage2 output (superset
#                                      of autotest.csv; used for everything)
#   fn.csv                          - ground-truth boxes the detector missed entirely
#   mantest.csv                     - full manual/GT annotation table (one row per
#                                      annotated animal), used for per-species GT totals
#   val_images24_yolov12.csv        - the complete list of evaluated images (includes
#                                      images with zero detections)
#   dataset_split_star.csv          - per-image train/test stratification + covariates
#
# OUTPUT
#   <OUTPUT_DIR>/*.png   - one PNG per analysis, plus a combined dashboard
#   <OUTPUT_DIR>/*.csv   - the underlying summary tables behind each plot
#   Console              - a printed scorecard at the end of the run
#
# USAGE
#   Rscript two_stage_pipeline_eval.R [data_dir] [output_dir]
#   (defaults: data_dir = ".", output_dir = "./pipeline_eval_output")
# =============================================================================

suppressWarnings(suppressMessages({
  
  # ---- package bootstrap -----------------------------------------------
  required_pkgs <- c("tidyverse", "scales", "patchwork")
  missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(required_pkgs, library, character.only = TRUE))
}))

# =============================================================================
# 0. CONFIG
# =============================================================================

args       <- commandArgs(trailingOnly = TRUE)
DATA_DIR   <- if (length(args) >= 1) args[[1]] else "../data/raw/star_twostage_yolo/"
OUTPUT_DIR <- if (length(args) >= 2) args[[2]] else file.path(DATA_DIR, "pipeline_eval_output")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

FILES <- list(
  cascade   = file.path(DATA_DIR, "autotest_two_stage_cascade.csv"),
  fn        = file.path(DATA_DIR, "fn.csv"),
  mantest   = file.path(DATA_DIR, "mantest.csv"),
  val_imgs  = file.path(DATA_DIR, "val_images24_yolov12.csv"),
  split     = file.path(DATA_DIR, "dataset_split_star.csv")
)

missing_files <- names(FILES)[!file.exists(unlist(FILES))]
if (length(missing_files) > 0) {
  stop(
    "Missing expected input file(s): ", paste(unlist(FILES[missing_files]), collapse = ", "),
    "\nSet DATA_DIR (arg 1) to the folder containing the pipeline CSVs."
  )
}

# The 6 taxonomic leaf classes for this run. Two are genus-only (no further
# species resolution exists for them yet) -- see taxonomy fix in project notes.
GENUS_ONLY_SPECIES <- c("henricia", "sclerasterias")
SPECIES_CLASSES    <- c("asterias_forbesi", "asterias_vulgaris", "astropecten_americanus",
                        "henricia", "leptasterias_tenera", "sclerasterias")
GENUS_CLASSES      <- c("Asterias", "Astropecten", "Henricia", "Leptasterias", "Sclerasterias")

# consistent plotting theme --------------------------------------------------
theme_pipeline <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.15)),
      plot.subtitle = element_text(color = "grey35"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}
theme_set(theme_pipeline())

PAL_OK   <- "#2a9d8f"
PAL_BAD  <- "#e76f51"
PAL_MID  <- "#e9c46a"
PAL_LINE <- "#264653"

save_plot <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 200, bg = "white")
  invisible(plot)
}

save_table <- function(df, name) {
  write_csv(df, file.path(OUTPUT_DIR, paste0(name, ".csv")))
  invisible(df)
}

# helper: turn a species class name into its genus, matching the capitalization
# convention used by the classifier's own `pred_genus` column ("Astropecten", etc.)
species_to_genus <- function(sp) {
  case_when(
    is.na(sp) ~ NA_character_,
    sp %in% GENUS_ONLY_SPECIES ~ str_to_title(sp),
    TRUE ~ str_to_title(str_split_fixed(sp, "_", 2)[, 1])
  )
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================

message("Loading data from: ", normalizePath(DATA_DIR))

# Booleans in these CSVs are written as "True"/"False". Coerce explicitly
# rather than trusting readr's type guessing, so results don't silently
# depend on how a given readr version parses mixed-case logicals.
as_logical_flexible <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  toupper(as.character(x)) %in% c("TRUE", "T", "1")
}

det <- read_csv(FILES$cascade, show_col_types = FALSE) %>%
  mutate(
    truedetect = as_logical_flexible(truedetect),
    species_is_resolved = as_logical_flexible(species_is_resolved),
    gt_genus = species_to_genus(gt_label)
  )

fn_tbl  <- read_csv(FILES$fn, show_col_types = FALSE)
man_tbl <- read_csv(FILES$mantest, show_col_types = FALSE)
val_imgs <- read_csv(FILES$val_imgs, show_col_types = FALSE)
split_tbl <- read_csv(FILES$split, show_col_types = FALSE) %>%
  mutate(
    is_train = as_logical_flexible(is_train),
    is_test = as_logical_flexible(is_test),
    is_empty = as_logical_flexible(is_empty)
  )

# ground-truth totals per species, taken directly from the manual annotation
# table (the authoritative source -- independent of anything the model did)
n_gt_by_species <- man_tbl %>%
  count(gt_label, name = "n_gt") %>%
  filter(gt_label %in% SPECIES_CLASSES)

n_gt_by_genus <- n_gt_by_species %>%
  mutate(genus = species_to_genus(gt_label)) %>%
  group_by(genus) %>%
  summarise(n_gt = sum(n_gt), .groups = "drop")

N_GT_TOTAL <- nrow(man_tbl)

# sanity cross-check: TP (truedetect & has gt_label) + FN (fn.csv rows) should
# reconcile against the manual annotation totals above. Warn (don't fail) if not.
reconcile <- det %>%
  filter(truedetect == TRUE) %>%
  count(gt_label, name = "n_tp") %>%
  full_join(fn_tbl %>% count(Spname, name = "n_fn"), by = c("gt_label" = "Spname")) %>%
  mutate(across(c(n_tp, n_fn), ~replace_na(., 0)), n_reconciled = n_tp + n_fn) %>%
  full_join(n_gt_by_species, by = "gt_label") %>%
  mutate(mismatch = n_reconciled != n_gt)

if (any(reconcile$mismatch, na.rm = TRUE)) {
  warning("GT totals do not fully reconcile between mantest.csv, autotest cascade, and fn.csv:\n",
          paste(capture.output(print(reconcile %>% filter(mismatch))), collapse = "\n"))
} else {
  message("GT reconciliation check passed: TP + FN == manual annotation totals for every class.")
}

# =============================================================================
# 2. SHARED METRIC HELPERS
# =============================================================================

#' Precision/recall curve + Average Precision (PASCAL VOC / COCO-style
#' monotonic-envelope integration) for a single class.
#'
#' @param score numeric ranking score (higher = more confident positive)
#' @param is_tp logical, TRUE where that scored item is a true positive
#' @param n_gt  total number of real positives for this class (TP + FN)
compute_pr_ap <- function(score, is_tp, n_gt) {
  ok <- !is.na(score)
  score <- score[ok]; is_tp <- is_tp[ok]
  
  if (n_gt == 0 || length(score) == 0) {
    return(list(curve = tibble(score = numeric(0), recall = numeric(0),
                               precision = numeric(0), precision_envelope = numeric(0)),
                ap = NA_real_))
  }
  
  ord <- order(score, decreasing = TRUE)
  is_tp_sorted <- is_tp[ord]
  score_sorted <- score[ord]
  
  cum_tp <- cumsum(is_tp_sorted)
  cum_fp <- cumsum(!is_tp_sorted)
  recall <- cum_tp / n_gt
  precision <- cum_tp / (cum_tp + cum_fp)
  
  # monotonic non-increasing envelope, integrated left-to-right over recall
  precision_env <- rev(cummax(rev(precision)))
  
  recall_pts    <- c(0, recall)
  precision_pts <- c(precision_env[1], precision_env)
  d_recall <- diff(recall_pts)
  ap <- sum(d_recall * precision_pts[-1])
  
  list(
    curve = tibble(score = score_sorted, recall = recall, precision = precision,
                   precision_envelope = precision_env),
    ap = ap
  )
}

#' Precision/recall/F1 at a fixed set of score thresholds (for picking an
#' operating point), independent of the full AP integration above.
threshold_sweep <- function(score, is_tp, n_gt, thresholds = seq(0.05, 0.95, by = 0.05)) {
  map_dfr(thresholds, function(t) {
    keep <- !is.na(score) & score >= t
    tp <- sum(is_tp[keep], na.rm = TRUE)
    fp <- sum(!is_tp[keep], na.rm = TRUE)
    fn <- n_gt - tp
    precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    recall    <- tp / n_gt
    f1 <- if (is.na(precision) || (precision + recall) == 0) NA_real_ else
      2 * precision * recall / (precision + recall)
    tibble(threshold = t, tp = tp, fp = fp, fn = fn,
           precision = precision, recall = recall, f1 = f1)
  })
}

#' Build a full (zero-filled) confusion matrix between two class vectors.
build_confusion <- function(truth, pred, classes) {
  tibble(truth = factor(truth, levels = classes), pred = factor(pred, levels = classes)) %>%
    count(truth, pred, .drop = FALSE) %>%
    group_by(truth) %>%
    mutate(row_total = sum(n), pct_of_truth = if_else(row_total > 0, n / row_total, 0)) %>%
    ungroup()
}

#' Per-class precision/recall/F1/support derived from a confusion matrix
#' produced by build_confusion().
per_class_metrics <- function(confusion) {
  totals_pred <- confusion %>%
    group_by(pred) %>%
    summarise(n_pred = sum(n), .groups = "drop") %>%
    mutate(pred = as.character(pred))
  confusion %>%
    filter(truth == pred) %>%
    transmute(class = as.character(truth), tp = n, support = row_total) %>%
    left_join(totals_pred, by = c("class" = "pred")) %>%
    mutate(
      precision = if_else(n_pred > 0, tp / n_pred, NA_real_),
      recall    = if_else(support > 0, tp / support, NA_real_),
      f1        = if_else(!is.na(precision) & !is.na(recall) & (precision + recall) > 0,
                          2 * precision * recall / (precision + recall), NA_real_)
    ) %>%
    select(class, support, precision, recall, f1)
}

# =============================================================================
# 3. STAGE 1 -- BROAD DETECTOR PERFORMANCE
# =============================================================================
# `det` has one row per candidate box the YOLO detector produced. truedetect
# indicates a correct localization (IoU-matched to a real GT box, regardless
# of species). n_gt for this stage is every real animal in the test set,
# whether the detector found it or not (TP + the entries in fn.csv).

message("\n--- Stage 1: detector ---")

stage1_pr  <- compute_pr_ap(det$Conf, det$truedetect, N_GT_TOTAL)
stage1_ap  <- stage1_pr$ap
stage1_sweep <- threshold_sweep(det$Conf, det$truedetect, N_GT_TOTAL)
save_table(stage1_sweep, "stage1_threshold_sweep")

best_f1_row <- stage1_sweep %>% filter(!is.na(f1)) %>% slice_max(f1, n = 1, with_ties = FALSE)

p_stage1_pr <- ggplot(stage1_pr$curve, aes(recall, precision_envelope)) +
  geom_line(color = PAL_LINE, linewidth = 1) +
  geom_ribbon(aes(ymin = 0, ymax = precision_envelope), fill = PAL_LINE, alpha = 0.08) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Stage 1 detector: precision-recall",
    subtitle = sprintf("AP = %.3f  |  %d ground-truth animals total", stage1_ap, N_GT_TOTAL),
    x = "Recall", y = "Precision"
  )
save_plot(p_stage1_pr, "01_stage1_pr_curve")

p_stage1_sweep <- stage1_sweep %>%
  pivot_longer(c(precision, recall, f1), names_to = "metric", values_to = "value") %>%
  ggplot(aes(threshold, value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = best_f1_row$threshold, linetype = "dashed", color = "grey40") +
  annotate("text", x = best_f1_row$threshold, y = 0.05,
           label = sprintf("best-F1 @ conf=%.2f", best_f1_row$threshold),
           hjust = -0.05, size = 3.2, color = "grey30") +
  scale_color_manual(values = c(precision = PAL_OK, recall = PAL_BAD, f1 = PAL_LINE)) +
  labs(title = "Stage 1 detector: metric vs. confidence threshold",
       x = "Confidence threshold", y = "Value", color = NULL)
save_plot(p_stage1_sweep, "02_stage1_threshold_sweep")

p_stage1_iou <- det %>%
  filter(truedetect == TRUE) %>%
  ggplot(aes(iu)) +
  geom_histogram(binwidth = 0.02, fill = PAL_LINE, boundary = 0) +
  geom_vline(xintercept = median(det$iu[det$truedetect], na.rm = TRUE),
             linetype = "dashed", color = PAL_BAD) +
  labs(title = "IoU distribution of matched (true-positive) detections",
       subtitle = sprintf("median IoU = %.2f", median(det$iu[det$truedetect], na.rm = TRUE)),
       x = "IoU with matched ground-truth box", y = "Count")
save_plot(p_stage1_iou, "03_stage1_iou_distribution")

# --- image-level behavior: coverage & false alarms on background images ----
images_with_gt   <- unique(man_tbl$Imagename)
images_with_det  <- unique(det$Imagename)
images_zero_det  <- setdiff(val_imgs$imagename, images_with_det)
images_fully_missed <- intersect(images_zero_det, images_with_gt)

empty_images <- split_tbl %>% filter(is_test == TRUE, total_annotations == 0) %>% pull(imagename)
fp_on_empty  <- det %>% filter(Imagename %in% empty_images, truedetect == FALSE)

image_level_summary <- tibble(
  metric = c("Test images (total)", "Test images with >=1 real animal (mantest)",
             "Test images the detector never fired on (0 proposals)",
             "...of which actually contained an animal (fully missed)",
             "Officially empty test images (no GT at all)",
             "...of which received >=1 false-positive detection",
             "False-positive detections landing on empty images"),
  value = c(nrow(val_imgs), length(images_with_gt), length(images_zero_det),
            length(images_fully_missed), length(empty_images),
            n_distinct(fp_on_empty$Imagename), nrow(fp_on_empty))
)
save_table(image_level_summary, "stage1_image_level_summary")

p_stage1_empty <- ggplot(
  tibble(
    group = c("Empty images\n(no FP)", "Empty images\n(>=1 FP)"),
    n = c(length(empty_images) - n_distinct(fp_on_empty$Imagename), n_distinct(fp_on_empty$Imagename))
  ),
  aes(group, n, fill = group)
) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, fontface = "bold") +
  scale_fill_manual(values = c(PAL_OK, PAL_BAD)) +
  labs(title = "False alarms on background-only images",
       subtitle = sprintf("%d of %d empty test images produced at least one false detection (%d FP boxes total)",
                          n_distinct(fp_on_empty$Imagename), length(empty_images), nrow(fp_on_empty)),
       x = NULL, y = "Number of images")
save_plot(p_stage1_empty, "04_stage1_false_alarms_on_empty_images")

# =============================================================================
# 4. STAGE 2 -- TAXONOMIC CLASSIFIER PERFORMANCE
# =============================================================================
# Classifier accuracy can only be judged on detections that were correctly
# localized AND have a known ground-truth label -- i.e. truedetect == TRUE.
# False positives (background boxes) have no real species to be right or
# wrong about, so they're excluded from this section (they're covered above).

message("--- Stage 2: taxonomic classifier ---")

tp_det <- det %>% filter(truedetect == TRUE, !is.na(gt_label))

species_confusion <- build_confusion(tp_det$gt_label, tp_det$stage2_species, SPECIES_CLASSES)
save_table(species_confusion, "stage2_species_confusion")

genus_confusion <- build_confusion(tp_det$gt_genus, tp_det$pred_genus, GENUS_CLASSES)
save_table(genus_confusion, "stage2_genus_confusion")

confusion_heatmap <- function(confusion, title, subtitle = NULL) {
  ggplot(confusion, aes(pred, truth, fill = pct_of_truth)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), size = 3.4,
              color = ifelse(confusion$pct_of_truth > 0.55, "white", "grey15")) +
    scale_fill_viridis_c(limits = c(0, 1), labels = percent, name = "% of true class") +
    labs(title = title, subtitle = subtitle, x = "Predicted", y = "Ground truth") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
}

p_species_confusion <- confusion_heatmap(
  species_confusion, "Stage 2 species-level confusion matrix",
  sprintf("%d true-positive detections with known species", nrow(tp_det))
)
save_plot(p_species_confusion, "05_stage2_confusion_species", width = 8.5, height = 7)

p_genus_confusion <- confusion_heatmap(
  genus_confusion, "Stage 2 genus-level confusion matrix",
  "Species collapsed to genus before comparing"
)
save_plot(p_genus_confusion, "06_stage2_confusion_genus", width = 7.5, height = 6.5)

species_metrics <- per_class_metrics(species_confusion) %>%
  mutate(is_genus_only = class %in% GENUS_ONLY_SPECIES)
save_table(species_metrics, "stage2_species_per_class_metrics")

stage2_accuracy <- mean(tp_det$gt_label == tp_det$stage2_species, na.rm = TRUE)
stage2_macro_f1 <- mean(species_metrics$f1, na.rm = TRUE)
stage2_balanced_acc <- mean(species_metrics$recall, na.rm = TRUE)

p_species_metrics <- species_metrics %>%
  pivot_longer(c(precision, recall, f1), names_to = "metric", values_to = "value") %>%
  mutate(class = fct_reorder(class, support)) %>%
  ggplot(aes(class, value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value))),
            position = position_dodge(width = 0.75), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c(precision = PAL_OK, recall = PAL_BAD, f1 = PAL_LINE)) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(title = "Stage 2 per-species precision / recall / F1",
       subtitle = sprintf("Overall accuracy = %.1f%%   |   macro F1 = %.2f   |   balanced accuracy = %.1f%%",
                          100 * stage2_accuracy, stage2_macro_f1, 100 * stage2_balanced_acc),
       x = NULL, y = NULL, fill = NULL)
save_plot(p_species_metrics, "07_stage2_per_class_metrics", width = 8.5, height = 5.5)

# --- genus-only "restraint" behavior ----------------------------------------
# henricia / sclerasterias have no resolved species below the genus. A model
# that correctly refuses to over-resolve these into a specific (wrong) species
# will predict the genus-only label itself (species_is_resolved == FALSE).
# Forcing a real species name here is a specific, avoidable error mode.
genus_only_tp <- tp_det %>% filter(gt_label %in% GENUS_ONLY_SPECIES)

restraint_summary <- genus_only_tp %>%
  mutate(behavior = if_else(species_is_resolved, "Wrongly over-resolved\nto a false species",
                            "Correctly kept as\ngenus-only")) %>%
  count(behavior)
save_table(restraint_summary, "stage2_genus_only_restraint")

restraint_rate <- if (nrow(genus_only_tp) == 0) NA_real_ else
  sum(!genus_only_tp$species_is_resolved) / nrow(genus_only_tp)

p_restraint <- ggplot(restraint_summary, aes(behavior, n, fill = behavior)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, fontface = "bold") +
  scale_fill_manual(values = setNames(c(PAL_OK, PAL_BAD),
                                      c("Correctly kept as\ngenus-only", "Wrongly over-resolved\nto a false species"))) +
  labs(title = "Genus-only classes: does the classifier over-resolve?",
       subtitle = sprintf("%d true positives where ground truth was genus-only (henricia/sclerasterias)\nrestraint rate = %.1f%%",
                          nrow(genus_only_tp), 100 * restraint_rate),
       x = NULL, y = "True positives")
save_plot(p_restraint, "08_stage2_genus_only_restraint", width = 6.5, height = 5.5)

# --- confidence calibration -------------------------------------------------
calibration <- tp_det %>%
  filter(!is.na(stage2_conf)) %>%
  mutate(correct = gt_label == stage2_species,
         conf_bin = cut(stage2_conf, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  group_by(conf_bin) %>%
  summarise(mean_conf = mean(stage2_conf), empirical_accuracy = mean(correct), n = n(), .groups = "drop") %>%
  filter(n > 0)
save_table(calibration, "stage2_calibration")

p_calibration <- ggplot(calibration, aes(mean_conf, empirical_accuracy)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(color = PAL_LINE) +
  geom_point(aes(size = n), color = PAL_LINE) +
  scale_size_continuous(name = "n detections", range = c(2, 8)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Stage 2 confidence calibration",
       subtitle = "Dashed line = perfect calibration. Points below it mean the classifier is overconfident.",
       x = "Mean predicted confidence (stage2_conf)", y = "Empirical accuracy")
save_plot(p_calibration, "09_stage2_calibration")

# =============================================================================
# 5. END-TO-END CASCADE PERFORMANCE (per-species / per-genus AP -> mAP)
# =============================================================================
# This is the metric that actually reflects what the pipeline delivers end to
# end: for a detection to count as a true positive of species X, it has to be
# BOTH correctly localized (truedetect) AND correctly classified (gt_label ==
# X). Following COCO/PASCAL-style per-category AP, the candidate pool for
# class X is every detection the pipeline's own argmax decision assigned to
# X (stage2_species == X), ranked by a combined "is this a real, correctly
# classified X" score (detector confidence x classifier confidence).
#
# NOTE: we deliberately do NOT rank by the raw, un-conditioned conf_<species>
# column across all 4740 detections. Because decoding is genus-conditioned
# (see project notes), a species' raw softmax value is only meaningful within
# its own genus -- using it globally would understate rare species nested
# inside a genus dominated by a more common one (e.g. asterias_forbesi living
# inside the Asterias genus alongside the much more common asterias_vulgaris).
# Scoring the model's actual argmax output is the fairer, deployment-faithful
# measure of "what does the cascade actually deliver."

message("--- End-to-end cascade: per-species / per-genus AP ---")

det <- det %>% mutate(combined_species_score = Conf * stage2_conf,
                      combined_genus_score = Conf * genus_conf)

species_ap_list <- map(SPECIES_CLASSES, function(sp) {
  candidates <- det %>% filter(stage2_species == sp)
  is_tp <- candidates$truedetect == TRUE & candidates$gt_label == sp
  n_gt  <- n_gt_by_species$n_gt[n_gt_by_species$gt_label == sp]
  n_gt  <- if (length(n_gt) == 0) 0 else n_gt
  compute_pr_ap(candidates$combined_species_score, is_tp, n_gt)
})
names(species_ap_list) <- SPECIES_CLASSES

species_ap_curves <- map_dfr(species_ap_list, "curve", .id = "class")
species_ap_table <- tibble(
  class = SPECIES_CLASSES,
  ap = map_dbl(species_ap_list, "ap"),
  n_gt = n_gt_by_species$n_gt[match(SPECIES_CLASSES, n_gt_by_species$gt_label)],
  n_candidates = map_dbl(species_ap_list, ~nrow(.x$curve))
)
species_mAP <- mean(species_ap_table$ap, na.rm = TRUE)
save_table(species_ap_table, "e2e_species_average_precision")

genus_ap_list <- map(GENUS_CLASSES, function(g) {
  candidates <- det %>% filter(pred_genus == g)
  is_tp <- candidates$truedetect == TRUE & candidates$gt_genus == g
  n_gt  <- n_gt_by_genus$n_gt[n_gt_by_genus$genus == g]
  n_gt  <- if (length(n_gt) == 0) 0 else n_gt
  compute_pr_ap(candidates$combined_genus_score, is_tp, n_gt)
})
names(genus_ap_list) <- GENUS_CLASSES
genus_ap_table <- tibble(
  class = GENUS_CLASSES,
  ap = map_dbl(genus_ap_list, "ap"),
  n_gt = n_gt_by_genus$n_gt[match(GENUS_CLASSES, n_gt_by_genus$genus)]
)
genus_mAP <- mean(genus_ap_table$ap, na.rm = TRUE)
save_table(genus_ap_table, "e2e_genus_average_precision")

p_species_pr_curves <- species_ap_curves %>%
  left_join(species_ap_table, by = "class") %>%
  mutate(facet_label = sprintf("%s\n(AP=%.2f, n=%d)", class, ap, n_gt)) %>%
  ggplot(aes(recall, precision_envelope)) +
  geom_line(color = PAL_LINE, linewidth = 1) +
  facet_wrap(~facet_label, ncol = 3) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "End-to-end per-species precision-recall (Stage 1 + Stage 2 combined)",
       subtitle = sprintf("mAP across %d species = %.3f", length(SPECIES_CLASSES), species_mAP),
       x = "Recall", y = "Precision")
save_plot(p_species_pr_curves, "10_e2e_pr_curves_species", width = 9, height = 6.5)

p_species_ap_bar <- species_ap_table %>%
  mutate(class = fct_reorder(class, ap)) %>%
  ggplot(aes(class, ap, fill = class %in% GENUS_ONLY_SPECIES)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f (n=%d)", ap, n_gt)), hjust = -0.05, size = 3.2) +
  scale_fill_manual(values = c(`TRUE` = PAL_MID, `FALSE` = PAL_LINE)) +
  geom_hline(yintercept = species_mAP, linetype = "dashed", color = PAL_BAD) +
  annotate("text", x = 0.7, y = species_mAP, label = sprintf("mAP = %.3f", species_mAP),
           vjust = -0.6, color = PAL_BAD, size = 3.3, hjust = 0) +
  coord_flip(ylim = c(0, 1.1)) +
  labs(title = "End-to-end Average Precision by species",
       subtitle = "Gold bars = genus-only classes (no species-level ground truth exists to fully resolve)",
       x = NULL, y = "Average Precision")
save_plot(p_species_ap_bar, "11_e2e_ap_by_species", width = 8, height = 5)

# =============================================================================
# 6. WORKFLOW STATUS
# =============================================================================
# A short, human-readable status board of where the pipeline stands. The bug
# entries below reflect the current known state of the codebase; edit this
# table as new issues are found/fixed so the dashboard stays current.

message("--- Workflow status ---")

pipeline_status <- tribble(
  ~component,                 ~issue,                                                      ~status,     ~note,
  "build_multiclass_yolo.py", "Class IDs assigned by alphabetically sorting CSV labels, disagreeing with taxonomy JSON ordering -> silent crop/gt_label cross-contamination", "Fixed", "Taxonomy JSON is now sole authority; data.yaml validated against it before exit",
  "extract_crops.py",         "Hardcoded to star24/ paths, unusable on other datasets without source edits",              "Fixed", "argparse added",
  "star_taxonomy.json / eval_two_stage_predictions.py", "Genus-only classes (henricia, sclerasterias) had no species output neuron -> guaranteed-wrong species classification for true genus-only positives", "Fixed", "Real species_id assigned with is_resolved_species flag; genus-conditioned decoding implemented",
  "train_classifier.py (weights)", "Output head shape changed by the genus-only taxonomy fix",                              "Retraining pending", "No code change needed -- picks up new species slots once taxonomy JSON is corrected"
)
save_table(pipeline_status, "workflow_pipeline_status")

status_colors <- c("Fixed" = PAL_OK, "Retraining pending" = PAL_MID, "Open" = PAL_BAD)

p_status <- pipeline_status %>%
  mutate(component = fct_rev(fct_inorder(component))) %>%
  ggplot(aes(x = 1, y = component, fill = status)) +
  geom_tile(width = 0.15, height = 0.7) +
  geom_text(aes(label = str_wrap(paste0(component, ": ", issue), 70)),
            x = 1.12, hjust = 0, size = 3.1, lineheight = 0.9) +
  scale_fill_manual(values = status_colors, name = NULL) +
  scale_x_continuous(limits = c(0.9, 3.2)) +
  labs(title = "Pipeline fix status", subtitle = NULL, x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = rel(1.15), hjust = 0))
save_plot(p_status, "12_workflow_pipeline_status", width = 10, height = 4.5)

# --- dataset composition / stratification (descriptive only) ---------------
# NOTE: dataset_split_star.csv also carries the covariates for downstream GAM
# density modeling. That modeling is explicitly out of scope here -- this
# section only checks that the train/test split looks sane and balanced.

split_counts <- split_tbl %>% count(dataset_level1, name = "n_images")
save_table(split_counts, "dataset_image_counts_by_split")

species_by_split <- split_tbl %>%
  select(dataset_level1, starts_with("n_")) %>%
  pivot_longer(-dataset_level1, names_to = "species", values_to = "n") %>%
  mutate(species = str_remove(species, "^n_")) %>%
  group_by(dataset_level1, species) %>%
  summarise(n = sum(n), .groups = "drop")
save_table(species_by_split, "dataset_species_counts_by_split")

p_dataset_composition <- species_by_split %>%
  ggplot(aes(fct_reorder(species, n, .fun = sum), n, fill = dataset_level1)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c(train = PAL_LINE, test = PAL_MID), name = "Split") +
  coord_flip() +
  labs(title = "Dataset composition: GT annotations by species and split",
       subtitle = sprintf("%d train images, %d test images",
                          split_counts$n_images[split_counts$dataset_level1 == "train"],
                          split_counts$n_images[split_counts$dataset_level1 == "test"]),
       x = NULL, y = "Number of annotated animals")
save_plot(p_dataset_composition, "13_dataset_composition", width = 8, height = 5)

p_dataset_spatial <- split_tbl %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  ggplot(aes(longitude, latitude, color = dataset_level1)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c(train = PAL_LINE, test = PAL_MID), name = "Split") +
  labs(title = "Spatial stratification of the train/test split",
       subtitle = "Context only -- density-based GAM stratification is out of scope for this script",
       x = "Longitude", y = "Latitude")
save_plot(p_dataset_spatial, "14_dataset_spatial_stratification", width = 7.5, height = 6)

# =============================================================================
# 7. SCORECARD
# =============================================================================

scorecard <- tribble(
  ~metric, ~value,
  "Stage 1 detector AP",                                          sprintf("%.3f", stage1_ap),
  "Stage 1 best-F1 operating point",                              sprintf("conf >= %.2f  (P=%.2f, R=%.2f, F1=%.2f)",
                                                                          best_f1_row$threshold, best_f1_row$precision,
                                                                          best_f1_row$recall, best_f1_row$f1),
  "Stage 1 false positives on empty images",                      sprintf("%d boxes across %d/%d empty images",
                                                                          nrow(fp_on_empty), n_distinct(fp_on_empty$Imagename),
                                                                          length(empty_images)),
  "Stage 1 fully-missed images (0 proposals, contained an animal)", as.character(length(images_fully_missed)),
  "Stage 2 species accuracy (on true positives)",                 sprintf("%.1f%%", 100 * stage2_accuracy),
  "Stage 2 macro F1 (on true positives)",                         sprintf("%.3f", stage2_macro_f1),
  "Stage 2 balanced accuracy (on true positives)",                sprintf("%.1f%%", 100 * stage2_balanced_acc),
  "Stage 2 genus-only restraint rate",                            sprintf("%.1f%%", 100 * restraint_rate),
  "End-to-end species mAP",                                       sprintf("%.3f", species_mAP),
  "End-to-end genus mAP",                                         sprintf("%.3f", genus_mAP),
  "Total ground-truth animals in test set",                       as.character(N_GT_TOTAL),
  "Total test images",                                            as.character(nrow(val_imgs))
)
save_table(scorecard, "00_scorecard")

message("\n==================== PIPELINE SCORECARD ====================")
print(scorecard, n = Inf)
message("==============================================================")

# --- combined dashboard (headline plots only) -------------------------------
dashboard <- (p_stage1_pr + p_species_confusion) /
  (p_species_ap_bar + p_status) +
  plot_annotation(title = "Two-stage pipeline -- status dashboard",
                  theme = theme(plot.title = element_text(face = "bold", size = 16)))
save_plot(dashboard, "00_dashboard_overview", width = 15, height = 11)

message("\nAll plots and tables written to: ", normalizePath(OUTPUT_DIR))
