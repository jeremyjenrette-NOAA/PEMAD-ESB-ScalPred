# ======================================================================
# n3_ablation.R - N3 (R2.5, N13): what does each calibration stage, and each
# added detector, contribute?
#
# Tiers, all scored on the same Holdout test images with block-bootstrap CIs:
#   a  F1 threshold count                         (per detector)
#   b  Stage 1 only: sum of detection-GAM p       (per detector)
#   c  Stage 2 on raw counts, no stage 1          (per detector)
#   d  Two-stage, one detector                    (per detector)
#   e  Two-stage, two detectors                   (each pair)
#   f  Two-stage, three detectors (+ RT-DETR-L)
# Within c-f the image-level form is chosen by block-fold CV RMSE on
# Calibration-fitting (lowest, as in gammod.R); fold SEs are reported so ties show.
#
# Concurvity-safe synergy forms (N13). The old forms used log_pred_diff =
# log1p(|C - Y|) next to log_yolo / log_cascade, which is a function of them
# (and drops the sign). Multi-detector forms here use a LEVEL plus SIGNED
# CONTRASTS of the log1p stage-1 sums:
#   two:   L = (l1 + l2)/2,            D = l1 - l2
#   three: L = (lY + lC + lR)/3,       H1 = lY - lC,  H2 = lR - (lY + lC)/2
# mgcv::concurvity ("worst") is reported for every selected model and for the
# old forms, for the response letter.
#
# RT-DETR-L has no stage 1 yet: its detection-GAM formula is chosen here by the
# same block CV + one-SE rule as gamstrat.R (cached), no random effects.
#
# Run from s/:  Rscript n3/n3_ablation.R     (~45-90 min the first time; RT-DETR CV cached)
# Needs: ../data/processed/eval_2226.rds, YOLOv12_deteval_2226.rds, CascadeR-CNN_deteval_2226.rds,
#        ../data/raw_model/2226data_rtdetr_l_gcp_vit_20260924_190833/autotest2226_rtdetr_l.csv
# Outputs: ../data/processed/n3/*.csv|rds, ../figures/ms_figures/2226_n3_ablation.pdf|png
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(stringr); library(janitor)
  library(mgcv); library(ggplot2); library(purrr) })
DET_ENGINE <- "gam"
source("./gamfunc.R"); source("./procfunc.R")
# Image-level selection rule: N3_SELECT=lowest (default; lowest CV RMSE) or
# N3_SELECT=one_se (simplest listed form within one SE of the lowest, as in stage 1)
SELECT_RULE <- Sys.getenv("N3_SELECT", "lowest")
SFX <- if (SELECT_RULE == "one_se") "_onese" else ""
OUT <- paste0("../data/processed/n3", SFX); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG <- "../figures/ms_figures"
N_FOLDS <- 5; N_BOOT <- 1000; SEED <- 2226
pr <- function(x) print(tibble::as_tibble(x), n = Inf, width = Inf)
logf <- file.path(OUT, "n3_log.txt"); cat("", file = logf)
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

# ----------------------------------------------------------------------
# 1. Stage 1 for RT-DETR-L (same data, folds and procedure as gamstrat.R)
# ----------------------------------------------------------------------
models <- readRDS("../data/processed/eval_2226.rds")
names(models)[names(models) == "YOLO"] <- "YOLOv12"
names(models)[names(models) == "Cascade"] <- "Cascade R-CNN"
img_meta <- models$YOLOv12$img %>% select(-any_of(c("n_auto", "man_density", "auto_density")))

rt_csv <- "../data/raw_model/2226data_rtdetr_l_gcp_vit_20260924_190833/autotest2226_rtdetr_l.csv"
at <- read.csv(rt_csv) %>% clean_names() %>%
  mutate(image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"),
         truedetect = tolower(as.character(truedetect)) == "true", spname = tolower(spname)) %>%
  filter(spname == "scallop", image_id %in% img_meta$image_id)
rt <- list(
  det = at %>% left_join(img_meta %>% select(-any_of("imagename")), by = "image_id"),
  img = img_meta %>% left_join(at %>% count(image_id, name = "n_auto"), by = "image_id") %>%
    mutate(n_auto = replace_na(n_auto, 0L))
)
say("RT-DETR-L detections on Calibration-set images: ", nrow(rt$det), " (TP ", sum(rt$det$truedetect), ")")

# gamstrat.R candidate set (verbatim)
candidate_forms <- list(
  M01_ConfOnly = truedetect ~ s(conf, k = 5),
  M02_Spatial  = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7),
  M03_Depth    = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5),
  M04_Optical  = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(field_of_view_sq_meter, k = 5),
  M05_Boxsize  = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(boxsize, k = 5),
  M06_Environ  = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(fluorometer_backscatter_ntu, k = 5),
  M07_teDepth   = truedetect ~ te(conf, bottom_depth, k = 5) + s(latitude, longitude, k = 7),
  M08_teAlt     = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5),
  M09_teOptics  = truedetect ~ te(conf, altitude, k = 5) + te(conf, field_of_view_sq_meter, k = 5) + s(latitude, longitude, k = 7),
  M10_teBoxsize = truedetect ~ te(conf, boxsize, k = 5) + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7),
  M11_PrevBest  = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7),
  M12_PB_Depth  = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(bottom_depth, k = 5),
  M13_PB_Backsc = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(fluorometer_backscatter_ntu, k = 5),
  M14_PB_Temp   = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(ctd_temperature_celsius, k = 5),
  M15_All_te    = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + te(conf, bottom_depth, k = 5)
)

rt_eval_file <- "../data/processed/RTDETR-L_deteval_2226.rds"
if (file.exists(rt_eval_file)) {
  say("Loading cached RT-DETR-L stage 1: ", rt_eval_file)
  rt_eval <- readRDS(rt_eval_file)
} else {
  rt$det <- add_re_factors(rt$det)
  cache <- "../data/processed/cache_formula_cv_RTDETR-L_2226.rds"
  fcv <- if (file.exists(cache)) readRDS(cache) else {
    say("RT-DETR-L stage-1 formula CV (15 candidates x 5 folds x 2 regions)...")
    x <- cv_select_detection_formula(rt$det, candidate_forms, K = N_FOLDS, engine = DET_ENGINE); saveRDS(x, cache); x }
  fsel <- choose_detection_formula(fcv)
  write.csv(fsel$one_se, file.path(OUT, "RTDETR-L_formula_cv_2226.csv"), row.names = FALSE)
  ids <- fsel$chosen
  say("RT-DETR-L stage 1 selected: GB = ", ids$GB, " | MAB = ", ids$MAB)
  forms <- list(GB = candidate_forms[[ids$GB]], MAB = candidate_forms[[ids$MAB]])
  tmp <- list(`RT-DETR-L` = list(img = rt$img %>% filter(dataset == "test_GAM_train"),
                                 det = rt$det %>% filter(dataset == "test_GAM_train")))
  best <- evaluate_pr_models(tmp, stratify_region = TRUE)$pr_all %>% group_by(model) %>%
    filter(f1 == max(f1, na.rm = TRUE)) %>% slice_max(conf, n = 1) %>% ungroup()
  f1_list <- list(GB = best$conf[grepl(" GB$", best$model)], MAB = best$conf[grepl(" MAB$", best$model)])
  say(sprintf("RT-DETR-L F1 thresholds: GB %.3f | MAB %.3f", f1_list$GB, f1_list$MAB))
  train_det <- det_fit_rows(rt$det %>% filter(dataset == "test_GAM_train"))
  gams <- lapply(c(GB = "GB", MAB = "MAB"), function(r) {
    f <- try_fit_det_gam(forms[[r]], train_det %>% filter(region == r), engine = DET_ENGINE, label = paste("RT-DETR final", r))
    if (!fit_ok(f)) stop("RT-DETR-L final stage-1 fit unusable in ", r); f })
  cf <- crossfit_detection_gams(rt$det, fixed_formulas = forms, K = N_FOLDS, nested_selection = FALSE, engine = DET_ENGINE)
  rt_eval <- test_calibration_gams(gams, cf$det, rt$img, f1_list, test_pred_source = "full")
  rt_eval$formulas <- ids
  saveRDS(rt_eval, rt_eval_file)
}
EV <- list(Y = readRDS("../data/processed/YOLOv12_deteval_2226.rds"),
           C = readRDS("../data/processed/CascadeR-CNN_deteval_2226.rds"),
           R = rt_eval)
DET_NAME <- c(Y = "YOLOv12", C = "Cascade R-CNN", R = "RT-DETR-L")
print(bind_rows(lapply(names(EV), function(k) EV[[k]]$metrics %>% mutate(detector = DET_NAME[[k]], .before = 1))) %>%
        select(detector, region, rmse_calibrated, pct_bias_total_calibrated, rmse_f1, pct_bias_total_f1))

# ----------------------------------------------------------------------
# 2. Image table: stage-1 sums (OOF on Calibration-fitting), raw and F1 counts,
#    crowding features per detector, level/contrast terms
# ----------------------------------------------------------------------
img <- EV$Y$img_eval %>% select(image_id, region, dataset, fold, n_annotations, transect_block_id)
for (k in names(EV)) {
  ie <- EV[[k]]$img_eval %>% select(image_id, pred = predicted_number, raw = raw_detection_number, f1 = predicted_f1_number)
  names(ie)[-1] <- paste0(names(ie)[-1], "_", k)
  img <- img %>% left_join(ie, by = "image_id")
  ft <- image_detection_features(EV[[k]]$det_eval)
  names(ft)[-1] <- paste0(names(ft)[-1], "_", k)
  img <- img %>% left_join(ft, by = "image_id")
}
img <- img %>% group_by(region) %>%
  mutate(across(starts_with(c("pred_", "raw_", "f1_", "n_det_all_", "overlap_frac_", "lowconf_frac_")), ~ coalesce(.x, 0)),
         across(starts_with("nn_ratio_"), ~ coalesce(.x, 5)),
         across(starts_with("med_box_mm_"), ~ coalesce(.x, median(.x, na.rm = TRUE)))) %>% ungroup() %>%
  mutate(l_Y = log1p(pred_Y), l_C = log1p(pred_C), l_R = log1p(pred_R),
         lraw_Y = log1p(raw_Y), lraw_C = log1p(raw_C), lraw_R = log1p(raw_R),
         L_YC = (l_Y + l_C) / 2, D_YC = l_Y - l_C,
         L_YR = (l_Y + l_R) / 2, D_YR = l_Y - l_R,
         L_CR = (l_C + l_R) / 2, D_CR = l_C - l_R,
         L3 = (l_Y + l_C + l_R) / 3, H1 = l_Y - l_C, H2 = l_R - (l_Y + l_C) / 2,
         # old (submitted) synergy covariates, for the concurvity comparison only
         log_yolo_pred = l_Y, log_cascade_pred = l_C, pred_diff = pred_C - pred_Y,
         log_pred_diff = log1p(abs(pred_C - pred_Y)))
stopifnot(!anyNA(img$l_R))

# ----------------------------------------------------------------------
# 3. Candidate image-level forms per tier and detector set
# ----------------------------------------------------------------------
fx <- function(...) as.formula(paste("n_annotations ~", paste(..., sep = " + ")))
single <- function(k) list(
  s_level      = fx(sprintf("s(l_%s)", k)),
  te_lowconf   = fx(sprintf("te(l_%s, lowconf_frac_%s, k = c(5, 4))", k, k)),
  te_box       = fx(sprintf("te(l_%s, med_box_mm_%s, k = c(5, 4))", k, k)),
  crowd        = fx(sprintf("s(l_%s)", k), sprintf("s(overlap_frac_%s, k = 5)", k), sprintf("s(nn_ratio_%s, k = 5)", k)),
  structure    = fx(sprintf("s(l_%s)", k), sprintf("s(med_box_mm_%s, k = 5)", k), sprintf("s(lowconf_frac_%s, k = 5)", k),
                    sprintf("s(overlap_frac_%s, k = 5)", k), sprintf("s(nn_ratio_%s, k = 5)", k)))
rawc <- function(k) list(
  s_raw      = fx(sprintf("s(lraw_%s)", k)),
  te_raw_low = fx(sprintf("te(lraw_%s, lowconf_frac_%s, k = c(5, 4))", k, k)))
pair <- function(a, b) { L <- paste0("L_", a, b); D <- paste0("D_", a, b); list(
  level          = fx(sprintf("s(%s)", L)),
  level_contrast = fx(sprintf("s(%s)", L), sprintf("s(%s, k = 5)", D)),
  te_level_contr = fx(sprintf("te(%s, %s, k = c(5, 4))", L, D)),
  te_lc_crowd    = fx(sprintf("te(%s, %s, k = c(5, 4))", L, D), "s(overlap_frac_Y, k = 5)", "s(nn_ratio_Y, k = 5)"),
  lc_lowconf     = fx(sprintf("s(%s)", L), sprintf("s(%s, k = 5)", D), "s(lowconf_frac_Y, k = 5)")) }
triple <- list(
  level          = fx("s(L3)"),
  level_contrast = fx("s(L3)", "s(H1, k = 5)", "s(H2, k = 5)"),
  te_L_H1        = fx("te(L3, H1, k = c(5, 4))", "s(H2, k = 5)"),
  lc_crowd       = fx("s(L3)", "s(H1, k = 5)", "s(H2, k = 5)", "s(overlap_frac_Y, k = 5)", "s(nn_ratio_Y, k = 5)"),
  lc_lowconf     = fx("s(L3)", "s(H1, k = 5)", "s(H2, k = 5)", "s(lowconf_frac_Y, k = 5)"))
SETS <- c(
  setNames(lapply(names(EV), rawc),   paste0("c|Stage 2 on raw counts|", DET_NAME)),
  setNames(lapply(names(EV), single), paste0("d|Two-stage, one detector|", DET_NAME)),
  list(`e|Two-stage, two detectors|YOLOv12 + Cascade R-CNN` = pair("Y", "C"),
       `e|Two-stage, two detectors|YOLOv12 + RT-DETR-L`     = pair("Y", "R"),
       `e|Two-stage, two detectors|Cascade R-CNN + RT-DETR-L` = pair("C", "R"),
       `f|Two-stage, three detectors|YOLOv12 + Cascade R-CNN + RT-DETR-L` = triple))

# ----------------------------------------------------------------------
# 4. Block-fold CV on Calibration-fitting (RMSE with fold SE), per region
# ----------------------------------------------------------------------
fit_tw <- function(f, d) tryCatch(gam(f, family = tw(link = "log"), data = d, method = "REML"), error = function(e) NULL)
cv_one <- function(f, d) {
  pred <- rep(NA_real_, nrow(d)); fr <- c()
  for (k in sort(unique(d$fold))) {
    g <- fit_tw(f, d[d$fold != k, ]); if (is.null(g)) return(NULL)
    ho <- d$fold == k; pred[ho] <- predict(g, d[ho, ], type = "response")
    fr <- c(fr, sqrt(mean((d$n_annotations[ho] - pred[ho])^2)))
  }
  tibble(cv_rmse = sqrt(mean((d$n_annotations - pred)^2)), cv_rmse_se = sd(fr) / sqrt(length(fr)),
         cv_pct_bias = 100 * (sum(pred) - sum(d$n_annotations)) / sum(d$n_annotations))
}
cal <- img %>% filter(dataset == "test_GAM_train", !is.na(fold))
cv_all <- imap_dfr(SETS, function(cands, key) {
  parts <- str_split_fixed(key, "\\|", 3)
  map_dfr(c("GB", "MAB"), function(r) {
    d <- cal %>% filter(region == r)
    imap_dfr(cands, function(f, id) {
      res <- cv_one(f, d)
      if (is.null(res)) res <- tibble(cv_rmse = NA_real_, cv_rmse_se = NA_real_, cv_pct_bias = NA_real_)
      res %>% mutate(tier = parts[1], tier_label = parts[2], set = parts[3], region = r, form = id,
                     formula = paste(deparse(f), collapse = " "), .before = 1)
    })
  })
})
cv_all <- cv_all %>% group_by(tier, set, region) %>%
  mutate(best = !is.na(cv_rmse) & cv_rmse == min(cv_rmse, na.rm = TRUE),
         within_1se = !is.na(cv_rmse) & cv_rmse <= min(cv_rmse, na.rm = TRUE) + cv_rmse_se[which(best)][1],
         pick = if (SELECT_RULE == "one_se") within_1se & cumsum(within_1se) == 1 else best) %>% ungroup()
write.csv(cv_all, file.path(OUT, "n3_cv_all.csv"), row.names = FALSE)
say("\nSelected image-level form per tier and detector set (rule: ", SELECT_RULE, "):")
sel <- cv_all %>% filter(pick)
pr(sel %>% select(tier, set, region, form, cv_rmse, cv_rmse_se, cv_pct_bias))

# ----------------------------------------------------------------------
# 5. Fit selected forms on Calibration-fitting, predict the Holdout test,
#    record concurvity
# ----------------------------------------------------------------------
hold <- img %>% filter(dataset == "test_GAM_test")
conc_worst <- function(g) { cc <- tryCatch(concurvity(g, full = TRUE), error = function(e) NULL)
  if (is.null(cc)) return(NA_real_); max(cc["worst", colnames(cc) != "para"]) }
conc <- list()
for (i in seq_len(nrow(sel))) {
  s1 <- sel[i, ]; meth <- paste0(s1$tier_label, ": ", s1$set)
  g <- fit_tw(SETS[[paste(s1$tier, s1$tier_label, s1$set, sep = "|")]][[s1$form]], cal %>% filter(region == s1$region))
  if (is.null(g)) { say("!! final fit failed: ", meth, " ", s1$region); next }
  if (!meth %in% names(hold)) hold[[meth]] <- NA_real_
  te <- hold$region == s1$region
  hold[[meth]][te] <- predict(g, hold[te, ], type = "response")
  conc[[length(conc) + 1]] <- tibble(method = meth, region = s1$region, form = s1$form, worst_concurvity = conc_worst(g),
                                     edf = sum(g$edf))
}
# Old (submitted) two-detector forms, for the concurvity comparison
OLD <- list(M10_YoloLogDiff = n_annotations ~ s(log_yolo_pred, log_pred_diff),
            M11_Additive    = n_annotations ~ s(log_yolo_pred) + s(log_cascade_pred) + s(log_pred_diff),
            M12_Tensor      = n_annotations ~ te(log_yolo_pred, log_cascade_pred),
            M13_TensorDiff  = n_annotations ~ ti(log_yolo_pred, log_cascade_pred) + log_pred_diff,
            M14_Champion    = n_annotations ~ s(log_yolo_pred, log_pred_diff) + s(log_cascade_pred, log_pred_diff),
            C10_TensorCrowd = n_annotations ~ te(log_yolo_pred, log_cascade_pred) + s(overlap_frac_Y, k = 5) + s(nn_ratio_Y, k = 5))
old_conc <- map_dfr(c("GB", "MAB"), function(r) imap_dfr(OLD, function(f, id) {
  g <- fit_tw(f, cal %>% filter(region == r))
  tibble(method = paste0("Submitted form ", id), region = r, form = id,
         worst_concurvity = if (is.null(g)) NA_real_ else conc_worst(g), edf = if (is.null(g)) NA_real_ else sum(g$edf)) }))
conc <- bind_rows(conc, old_conc)
write.csv(conc, file.path(OUT, "n3_concurvity.csv"), row.names = FALSE)
say("\nWorst-case concurvity (1 = a term is an exact function of the others; > 0.8 is a warning):")
pr(conc %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

# ----------------------------------------------------------------------
# 6. Holdout metrics with block-bootstrap CIs
# ----------------------------------------------------------------------
for (k in names(EV)) {
  hold[[paste0("F1 threshold: ", DET_NAME[[k]])]] <- hold[[paste0("f1_", k)]]
  hold[[paste0("Stage 1 only: ", DET_NAME[[k]])]] <- hold[[paste0("pred_", k)]]
}
TIER_ORDER <- c("F1 threshold", "Stage 1 only", "Stage 2 on raw counts", "Two-stage, one detector",
                "Two-stage, two detectors", "Two-stage, three detectors")
meth_cols <- names(hold)[str_detect(names(hold), paste0("^(", paste(TIER_ORDER, collapse = "|"), "): "))]
meth_cols <- meth_cols[order(match(str_extract(meth_cols, "^[^:]+"), TIER_ORDER), meth_cols)]
REF <- "Two-stage, two detectors: YOLOv12 + Cascade R-CNN"
boot <- block_bootstrap_metrics(hold, setNames(meth_cols, meth_cols), block = "transect_block_id",
                                ref = REF, B = N_BOOT, seed = SEED) %>%
  mutate(tier = str_extract(method, "^[^:]+"), set = str_remove(method, "^[^:]+: "))
write.csv(boot, file.path(OUT, "n3_holdout_metrics.csv"), row.names = FALSE)
saveRDS(list(img = img, hold = hold, cv = cv_all, conc = conc), file.path(OUT, "n3_results.rds"))
say("\nHoldout test (95% block-bootstrap intervals; d_rmse relative to ", REF, "):")
pr(boot %>% select(region, method, rmse, rmse_lo, rmse_hi, r2_1to1, pct_bias_total, bias_lo, bias_hi,
                   d_rmse_vs_ref, d_rmse_lo, d_rmse_hi) %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

# ----------------------------------------------------------------------
# 7. Figure: RMSE and total bias by tier (point + 95% interval)
# ----------------------------------------------------------------------
pdat <- boot %>%
  select(region, tier, set, rmse, rmse_lo, rmse_hi, pct_bias_total, bias_lo, bias_hi) %>%
  pivot_longer(c(rmse, pct_bias_total), names_to = "metric", values_to = "value") %>%
  mutate(lo = if_else(metric == "rmse", rmse_lo, bias_lo), hi = if_else(metric == "rmse", rmse_hi, bias_hi),
         metric = recode(metric, rmse = "Holdout RMSE (scallops per image)", pct_bias_total = "Holdout total bias (%)"),
         tier = factor(tier, levels = TIER_ORDER),
         set = factor(set, levels = rev(unique(set[order(tier)]))))
p <- ggplot(pdat, aes(value, set, colour = tier)) +
  geom_vline(data = tibble(metric = "Holdout total bias (%)", x = 0), aes(xintercept = x), linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.5, orientation = "y") +
  geom_point(size = 1.8) +
  facet_grid(tier ~ region + metric, scales = "free", space = "free_y", switch = "y",
             labeller = labeller(tier = label_wrap_gen(14))) +
  scale_colour_manual(values = c("#9e9e9e", "#6b7585", "#e87ba4", "#1baf7a", "#2a78d6", "#6250d6"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
                                  panel.grid.minor = element_blank())
ggsave(file.path(FIG, paste0("2226_n3_ablation", SFX, ".pdf")), p, width = 10, height = 6.5)
ggsave(file.path(FIG, paste0("2226_n3_ablation", SFX, ".png")), p, width = 10, height = 6.5, dpi = 300, bg = "white")
say("\nDone. Outputs in ", OUT)
