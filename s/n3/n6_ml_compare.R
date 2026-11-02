# ======================================================================
# n6_ml_compare.R - N6 (R2.4): would random forests or boosted trees do
# better than the GAMs in either calibration stage?
#
# YOLOv12 only (the recommended workflow). Same data, same block folds, same
# features as the GAMs, so only the model class changes:
#   stage 1 (one row per detection; is it a scallop?):
#     conf, altitude, field_of_view_sq_meter, boxsize, bottom_depth, latitude, longitude
#   stage 2 (one row per image; manual count):
#     stage-1 sum (pred), share below the F1 threshold, median box size,
#     overlap share, nearest-neighbour spacing, number of detections
# Workflows scored on the Holdout test (95% transect-block bootstrap):
#   GAM + GAM (reference, from n3_ablation.R) | GAM + RF | GAM + BRT | RF + RF | BRT + BRT
# Training-image stage-1 sums are out-of-fold (5 block folds), as for the GAMs.
# Hyperparameters are tuned by the same block-fold CV on Calibration-fitting
# (small grids; BRT depth x rounds on a grid). Also reports stage-1 log loss and
# stage-2 CV RMSE next to the GAM values.
#
# Run from s/:  Rscript n3/n6_ml_compare.R      (about 20-40 min on a laptop)
# Needs: ../data/processed/YOLOv12_deteval_2226.rds (gamstrat.R)
#        ../data/processed/n3/n3_results.rds       (n3_ablation.R)
#        R packages ranger and xgboost (install.packages(c("ranger", "xgboost")))
# Writes ../data/processed/n6_ml/*.csv, ../figures/ms_figures/2226_n6_ml_compare.pdf|png
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(ggplot2)
  library(ranger); library(xgboost) })
source("./gamfunc.R")   # image_detection_features(), add_crowding_features(), block_bootstrap_metrics()
OUT <- "../data/processed/n6_ml"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG <- "../figures/ms_figures"
SEED <- 2226; N_BOOT <- 1000; NTHREAD <- max(1, parallel::detectCores() - 1)
set.seed(SEED)
pr <- function(x) print(tibble::as_tibble(x), n = Inf, width = Inf)
# BRT rounds are tuned on a grid (no early stopping), so the code works with
# both the old (1.7) and new (>= 2) xgboost R interfaces
NROUNDS <- c(150, 400, 1000)
xgb_fit <- function(params, x, y, nrounds) xgb.train(params = params, data = xgb.DMatrix(as.matrix(x), label = y),
                                                     nrounds = nrounds, verbose = 0)
xgb_pred <- function(m, x) predict(m, xgb.DMatrix(as.matrix(x)))
logf <- file.path(OUT, "n6_log.txt"); cat("", file = logf)
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

# ----------------------------------------------------------------------
# 1. Data
# ----------------------------------------------------------------------
ev  <- readRDS("../data/processed/YOLOv12_deteval_2226.rds")
det <- ev$det_eval %>% mutate(region = as.character(region), y = as.integer(as.logical(truedetect)))
F1 <- c("conf", "altitude", "field_of_view_sq_meter", "boxsize", "bottom_depth", "latitude", "longitude")
stopifnot(all(c(F1, "fold", "dataset", "image_id", "pred_p") %in% names(det)))
# sentinel / missing covariates -> regional median (ranger needs complete rows)
det <- det %>% group_by(region) %>%
  mutate(across(all_of(F1), ~ { v <- if_else(.x <= -90, NA_real_, as.numeric(.x)); coalesce(v, median(v, na.rm = TRUE)) })) %>%
  ungroup()
n3  <- readRDS("../data/processed/n3/n3_results.rds")
img <- n3$img %>% select(image_id, region, dataset, fold, n_annotations, transect_block_id,
                         pred_Y, n_det_all_Y, lowconf_frac_Y, med_box_mm_Y, overlap_frac_Y, nn_ratio_Y) %>%
  mutate(region = as.character(region))
ref <- n3$hold %>% select(image_id, gam_gam = `Two-stage, one detector: YOLOv12`)
F2 <- c("pred", "lowconf_frac_Y", "med_box_mm_Y", "overlap_frac_Y", "nn_ratio_Y", "n_det_all_Y")
regions <- c("GB", "MAB")
say("Detections: ", nrow(det), " | Calibration-fitting images: ", sum(img$dataset == "test_GAM_train"),
    " | Holdout images: ", sum(img$dataset == "test_GAM_test"))

logloss <- function(y, p) { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); -mean(y * log(p) + (1 - y) * log(1 - p)) }
folds_of <- function(d) sort(unique(na.omit(d$fold)))

# ----------------------------------------------------------------------
# 2. Stage 1 learners: out-of-fold on Calibration-fitting, full fit for Holdout
# ----------------------------------------------------------------------
rf1_fit <- function(d, mns) ranger(x = d[, F1], y = factor(d$y, levels = c(0, 1)), probability = TRUE,
                                   num.trees = 500, min.node.size = mns, num.threads = NTHREAD, seed = SEED)
rf1_pred <- function(m, d) predict(m, d[, F1], num.threads = NTHREAD)$predictions[, "1"]
brt1_params <- function(depth) list(objective = "binary:logistic", eta = 0.05,
                                    max_depth = depth, subsample = 0.8, colsample_bytree = 0.8,
                                    min_child_weight = 5, nthread = NTHREAD)

stage1 <- function(r) {
  tr <- det %>% filter(region == r, dataset == "test_GAM_train", !is.na(fold))
  te <- det %>% filter(region == r, dataset == "test_GAM_test")
  ks <- folds_of(tr)
  # --- RF: tune min.node.size by block CV
  rf_grid <- c(10, 50, 200)
  rf_cv <- map_dfr(rf_grid, function(mns) {
    oof <- rep(NA_real_, nrow(tr))
    for (k in ks) { m <- rf1_fit(tr[tr$fold != k, ], mns); oof[tr$fold == k] <- rf1_pred(m, tr[tr$fold == k, ]) }
    tibble(mns = mns, logloss = logloss(tr$y, oof), oof = list(oof))
  })
  rb <- rf_cv[which.min(rf_cv$logloss), ]
  rf_full <- rf1_fit(tr, rb$mns)
  # --- BRT: tune depth x rounds by block CV
  brt_cv <- pmap_dfr(expand_grid(depth = c(2, 4, 6), nrounds = NROUNDS), function(depth, nrounds) {
    oof <- rep(NA_real_, nrow(tr))
    for (k in ks) { a <- tr$fold != k
      m <- xgb_fit(brt1_params(depth), tr[a, F1], tr$y[a], nrounds); oof[!a] <- xgb_pred(m, tr[!a, F1]) }
    tibble(depth = depth, nrounds = nrounds, logloss = logloss(tr$y, oof), oof = list(oof))
  })
  bb <- brt_cv[which.min(brt_cv$logloss), ]
  brt_full <- xgb_fit(brt1_params(bb$depth), tr[, F1], tr$y, bb$nrounds)
  say(sprintf("  %s stage 1: RF min.node.size %d (log loss %.4f) | BRT depth %d, %d rounds (log loss %.4f) | GAM OOF log loss %.4f",
              r, rb$mns, rb$logloss, bb$depth, bb$nrounds, bb$logloss, logloss(tr$y, tr$pred_p)))
  out <- bind_rows(
    tr %>% transmute(det_row = row_id, image_id, p_rf = rb$oof[[1]], p_brt = bb$oof[[1]]),
    te %>% transmute(det_row = row_id, image_id, p_rf = rf1_pred(rf_full, te), p_brt = xgb_pred(brt_full, te[, F1])))
  ll <- tibble(region = r,
               set = c("Calibration-fitting (out-of-fold)", "Holdout test"),
               GAM = c(logloss(tr$y, tr$pred_p), logloss(te$y, te$pred_p)),
               RF  = c(rb$logloss, logloss(te$y, rf1_pred(rf_full, te))),
               BRT = c(bb$logloss, logloss(te$y, xgb_pred(brt_full, te[, F1]))))
  list(p = out, ll = ll, tune = tibble(region = r, rf_min_node = rb$mns, brt_depth = bb$depth, brt_rounds = bb$nrounds))
}
det$row_id <- seq_len(nrow(det))
say("Stage 1 ...")
s1 <- lapply(setNames(regions, regions), stage1)
p1 <- bind_rows(lapply(s1, `[[`, "p"))
ll1 <- bind_rows(lapply(s1, `[[`, "ll"))
write.csv(ll1, file.path(OUT, "n6_stage1_logloss.csv"), row.names = FALSE)
say("\nStage 1 log loss (lower is better):"); pr(ll1 %>% mutate(across(c(GAM, RF, BRT), ~ round(.x, 4))))

# Stage-1 sums per image for each learner (OOF on Calibration-fitting)
sums <- p1 %>% group_by(image_id) %>% summarise(pred_rf = sum(p_rf), pred_brt = sum(p_brt), .groups = "drop")
img <- img %>% left_join(sums, by = "image_id") %>%
  mutate(pred_rf = coalesce(pred_rf, 0), pred_brt = coalesce(pred_brt, 0))

# ----------------------------------------------------------------------
# 3. Stage 2 learners on the image table
# ----------------------------------------------------------------------
X2 <- function(d, s1col) { x <- d[, F2[-1]]; x$pred <- d[[s1col]]; as.data.frame(x[, F2]) }
rf2_fit  <- function(x, y, mns, mtry) ranger(x = x, y = y, num.trees = 500, min.node.size = mns, mtry = mtry,
                                             num.threads = NTHREAD, seed = SEED)
brt2_params <- function(depth) list(objective = "reg:tweedie", tweedie_variance_power = 1.5,
                                    eta = 0.05, max_depth = depth, subsample = 0.8, colsample_bytree = 0.8,
                                    min_child_weight = 5, nthread = NTHREAD)
rmse <- function(y, p) sqrt(mean((y - p)^2))

stage2 <- function(r, learner, s1col) {
  tr <- img %>% filter(region == r, dataset == "test_GAM_train", !is.na(fold))
  te <- img %>% filter(region == r, dataset == "test_GAM_test")
  xtr <- X2(tr, s1col); xte <- X2(te, s1col); ks <- folds_of(tr)
  if (learner == "RF") {
    grid <- expand_grid(mns = c(5, 20, 50), mtry = c(2, 4))
    cv <- pmap_dfr(grid, function(mns, mtry) {
      oof <- rep(NA_real_, nrow(tr))
      for (k in ks) { m <- rf2_fit(xtr[tr$fold != k, ], tr$n_annotations[tr$fold != k], mns, mtry)
        oof[tr$fold == k] <- predict(m, xtr[tr$fold == k, ], num.threads = NTHREAD)$predictions }
      tibble(mns = mns, mtry = mtry, cv_rmse = rmse(tr$n_annotations, oof))
    })
    b <- cv[which.min(cv$cv_rmse), ]
    m <- rf2_fit(xtr, tr$n_annotations, b$mns, b$mtry)
    pred <- predict(m, xte, num.threads = NTHREAD)$predictions
    tune <- sprintf("min.node.size %d, mtry %d", b$mns, b$mtry)
  } else {
    cv <- pmap_dfr(expand_grid(depth = c(2, 3, 4), nrounds = NROUNDS), function(depth, nrounds) {
      oof <- rep(NA_real_, nrow(tr))
      for (k in ks) { a <- tr$fold != k
        m <- xgb_fit(brt2_params(depth), xtr[a, ], tr$n_annotations[a], nrounds); oof[!a] <- xgb_pred(m, xtr[!a, ]) }
      tibble(depth = depth, nrounds = nrounds, cv_rmse = rmse(tr$n_annotations, oof))
    })
    b <- cv[which.min(cv$cv_rmse), ]
    m <- xgb_fit(brt2_params(b$depth), xtr, tr$n_annotations, b$nrounds)
    pred <- xgb_pred(m, xte)
    tune <- sprintf("depth %d, %d rounds", b$depth, b$nrounds)
  }
  list(pred = tibble(image_id = te$image_id, pred = pred),
       info = tibble(region = r, learner = learner, stage1 = s1col, cv_rmse = b$cv_rmse, tuning = tune))
}

WF <- list("GAM + RF"  = c("RF",  "pred_Y"),  "GAM + BRT" = c("BRT", "pred_Y"),
           "RF + RF"   = c("RF",  "pred_rf"), "BRT + BRT" = c("BRT", "pred_brt"))
say("\nStage 2 ...")
hold <- img %>% filter(dataset == "test_GAM_test") %>% left_join(ref, by = "image_id")
info2 <- list()
for (w in names(WF)) {
  col <- paste0("wf_", gsub("[^A-Za-z]", "", w)); hold[[col]] <- NA_real_
  for (r in regions) {
    res <- stage2(r, WF[[w]][1], WF[[w]][2])
    i <- match(res$pred$image_id, hold$image_id); hold[[col]][i] <- res$pred$pred
    info2[[length(info2) + 1]] <- res$info %>% mutate(workflow = w, .before = 1)
    say(sprintf("  %s %s: CV RMSE %.3f (%s)", w, r, res$info$cv_rmse, res$info$tuning))
  }
}
info2 <- bind_rows(info2)
gam_cv <- n3$cv %>% filter(tier == "d", set == "YOLOv12", best) %>% select(region, gam_cv_rmse = cv_rmse, gam_form = form)
info2 <- info2 %>% left_join(gam_cv, by = "region")
write.csv(info2, file.path(OUT, "n6_stage2_cv.csv"), row.names = FALSE)
say("\nStage 2 CV RMSE on Calibration-fitting (GAM for comparison):"); pr(info2)

# ----------------------------------------------------------------------
# 4. Holdout comparison
# ----------------------------------------------------------------------
PRED <- c("GAM + GAM" = "gam_gam", "GAM + RF" = "wf_GAMRF", "GAM + BRT" = "wf_GAMBRT",
          "RF + RF" = "wf_RFRF", "BRT + BRT" = "wf_BRTBRT")
stopifnot(!anyNA(hold[, unname(PRED)]))
boot <- block_bootstrap_metrics(hold, PRED, block = "transect_block_id", ref = "GAM + GAM", B = N_BOOT, seed = SEED)
write.csv(boot, file.path(OUT, "n6_holdout_metrics.csv"), row.names = FALSE)
write.csv(hold %>% select(image_id, region, transect_block_id, n_annotations, all_of(unname(PRED))),
          file.path(OUT, "n6_holdout_predictions.csv"), row.names = FALSE)
say("\nHoldout test (95% block-bootstrap; d_rmse relative to GAM + GAM):")
pr(boot %>% select(region, method, rmse, rmse_lo, rmse_hi, pct_bias_total, bias_lo, bias_hi,
                   d_rmse_vs_ref, d_rmse_lo, d_rmse_hi) %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

pd <- boot %>% select(region, method, rmse, rmse_lo, rmse_hi, pct_bias_total, bias_lo, bias_hi) %>%
  pivot_longer(c(rmse, pct_bias_total), names_to = "metric", values_to = "value") %>%
  mutate(lo = if_else(metric == "rmse", rmse_lo, bias_lo), hi = if_else(metric == "rmse", rmse_hi, bias_hi),
         metric = recode(metric, rmse = "Holdout RMSE (scallops per image)", pct_bias_total = "Holdout total bias (%)"),
         region = recode(region, GB = "Georges Bank", MAB = "Mid-Atlantic Bight"),
         method = factor(method, levels = rev(names(PRED))))
p <- ggplot(pd, aes(value, method)) +
  geom_vline(data = tibble(metric = "Holdout total bias (%)", x = 0), aes(xintercept = x), linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", colour = "#2a78d6") +
  geom_point(colour = "#2a78d6", size = 1.8) +
  facet_grid(. ~ region + metric, scales = "free_x") +
  labs(x = NULL, y = "Stage 1 + stage 2 model") + theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(FIG, "2226_n6_ml_compare.pdf"), p, width = 10, height = 3)
ggsave(file.path(FIG, "2226_n6_ml_compare.png"), p, width = 10, height = 3, dpi = 300, bg = "white")
say("\nDone. Outputs in ", OUT)
