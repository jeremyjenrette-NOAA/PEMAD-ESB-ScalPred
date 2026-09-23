library(dplyr)
library(tidyr)
library(mgcv)

# ======================================================================
# Cross-fitting utilities for the two-stage (detection GAM -> synergy GAM)
# calibration pipeline.
#
# Why: the synergy GAM is trained on image-level sums of detection-GAM
# probabilities. If those sums come from a detection GAM that was fit on the
# same images, the synergy GAM learns from in-sample (over-fit) inputs but is
# applied to out-of-sample inputs at test time. Cross-fitting gives every
# test_GAM_train image an out-of-fold (OOF) prediction from a detection GAM
# that never saw that image, so stage-2 training inputs match test inputs.
#
# Design:
#   * Folds are assigned to whole transect blocks (transect_block_id), never
#     to individual detections or images, so neighbouring frames stay together.
#   * Folds are assigned separately within each region (GB / MAB), matching the
#     region-specific GAMs.
#   * High-density blocks are dealt out first and balanced across folds.
#   * The same fold key is used for every detector (assigned once in makedat.R).
#   * Formula selection is nested: within each fold the candidate formulas are
#     ranked by AIC on the fitting folds only, so the held-out fold never
#     influences which formula predicts it.
# ======================================================================

# Columns that must be non-missing for a detection to enter a detection-GAM fit.
# Kept identical to the filter used in gamfunc.R.
DET_REQUIRED_VARS <- c("truedetect", "conf", "bottom_depth", "latitude", "longitude",
                       "altitude", "field_of_view_sq_meter", "fluorometer_backscatter_ntu",
                       "ctd_temperature_celsius", "boxsize", "iu")

det_fit_rows <- function(det_df) {
  det_df %>% drop_na(all_of(DET_REQUIRED_VARS))
}

# ----------------------------------------------------------------------
# 1. Fold assignment (block-grouped, region-specific, density-balanced)
# ----------------------------------------------------------------------
assign_folds <- function(meta_df, K = 5, seed = 2226, fold_set = "test_GAM_train",
                         block_col = "transect_block_id") {
  
  # Highest density stratum = the density_bin with the largest lower bound
  bins     <- unique(na.omit(meta_df$density_bin))
  lower    <- suppressWarnings(as.numeric(sub("^[\\(\\[]\\s*([-0-9.eE]+)\\s*,.*$", "\\1", bins)))
  high_bin <- bins[which.max(lower)]
  cat("assign_folds: balancing high-density bin", high_bin, "across", K, "folds\n")
  
  set.seed(seed)
  blocks <- meta_df %>%
    filter(dataset == fold_set) %>%
    group_by(region, block = .data[[block_col]]) %>%
    summarise(n_img  = n(),
              n_high = sum(density_bin == high_bin, na.rm = TRUE),
              n_ann  = sum(n_annotations, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(u = runif(n())) %>%
    arrange(region, desc(n_high), desc(n_ann), u)
  
  # Greedy balancing within region: blocks holding high-density images go to the
  # fold with the fewest high-density images so far; all other blocks go to the
  # fold with the fewest images. Ties are broken at random.
  pick <- function(score) {
    cand <- which(score == min(score))
    if (length(cand) == 1) cand else sample(cand, 1)
  }
  blocks$fold <- NA_integer_
  for (r in unique(blocks$region)) {
    idx <- which(blocks$region == r)
    hi <- numeric(K); ni <- numeric(K)
    for (i in idx) {
      k <- if (blocks$n_high[i] > 0) pick(hi * 1e6 + ni) else pick(ni * 1e6 + hi)
      blocks$fold[i] <- k
      hi[k] <- hi[k] + blocks$n_high[i]
      ni[k] <- ni[k] + blocks$n_img[i]
    }
  }
  
  out <- meta_df %>%
    left_join(blocks %>% select(region, block, fold),
              by = setNames(c("region", "block"), c("region", block_col))) %>%
    # Only fold_set images carry a fold; test_GAM_test images from shared
    # blocks must stay NA so they are never used inside cross-fitting.
    mutate(fold = if_else(dataset == fold_set, fold, NA_integer_))
  
  cat("Fold balance (images per fold, by region x density bin):\n")
  print(out %>% filter(!is.na(fold)) %>% count(region, density_bin, fold) %>%
          pivot_wider(names_from = fold, values_from = n, names_prefix = "fold_"))
  out
}

# ----------------------------------------------------------------------
# 2. Fit candidate formulas, return the AIC table and the best fit
# ----------------------------------------------------------------------
fit_best_detection_gam <- function(train_data, candidate_formulas, verbose = FALSE) {
  best_fit <- NULL; best_aic <- Inf
  tab <- vector("list", length(candidate_formulas))
  for (j in seq_along(candidate_formulas)) {
    id  <- names(candidate_formulas)[j]
    if (verbose) cat("    Fitting:", id, "\n")
    fit <- gam(candidate_formulas[[j]], family = binomial(), data = train_data, method = "REML")
    a   <- AIC(fit)
    tab[[j]] <- tibble(Model_ID = id, AIC = a, Deviance_Explained = summary(fit)$dev.expl)
    if (a < best_aic) { best_aic <- a; best_fit <- fit; best_id <- id }
    rm(fit)
  }
  list(fit = best_fit, id = best_id, table = bind_rows(tab) %>% arrange(AIC))
}

# Convenience: pick the formula object with the lowest AIC from a
# compare_detection_gams() table.
pick_best_formula <- function(selection_table, candidate_formulas) {
  best_id <- selection_table$Model_ID[which.min(selection_table$AIC)]
  cat("  [Selected]", best_id, "\n")
  candidate_formulas[[best_id]]
}

# ----------------------------------------------------------------------
# 3. Cross-fit the detection GAMs -> out-of-fold pred_p for test_GAM_train
# ----------------------------------------------------------------------
#   det_df             detection table with `dataset`, `region`, `fold` columns
#   candidate_formulas named list of formulas (used when nested_selection = TRUE)
#   fixed_formulas     list(GB = , MAB = ) used when nested_selection = FALSE
#
# Returns list(det = det_df + pred_p_oof, selection = per-fold chosen formula)
crossfit_detection_gams <- function(det_df, candidate_formulas = NULL,
                                    fixed_formulas = NULL, K = 5,
                                    nested_selection = TRUE,
                                    fold_set = "test_GAM_train") {
  if (!"fold" %in% names(det_df)) {
    stop("No `fold` column in det_df. Re-run makedat.R (assign_folds) to rebuild eval_2226.rds.")
  }
  if (nested_selection && is.null(candidate_formulas)) stop("candidate_formulas required for nested selection.")
  if (!nested_selection && is.null(fixed_formulas)) stop("fixed_formulas required when nested_selection = FALSE.")
  
  det_df$pred_p_oof <- NA_real_
  sel <- list()
  t0 <- Sys.time()
  
  for (r in c("GB", "MAB")) {
    for (k in seq_len(K)) {
      in_set  <- det_df$dataset == fold_set & det_df$region == r & !is.na(det_df$fold)
      fit_idx <- which(in_set & det_df$fold != k)
      ho_idx  <- which(in_set & det_df$fold == k)
      if (length(ho_idx) == 0) next
      
      train_k <- det_fit_rows(det_df[fit_idx, ])
      
      if (nested_selection) {
        best <- fit_best_detection_gam(train_k, candidate_formulas)
        fit <- best$fit; chosen <- best$id; aic <- min(best$table$AIC)
      } else {
        fit <- gam(fixed_formulas[[r]], family = binomial(), data = train_k, method = "REML")
        chosen <- "fixed"; aic <- AIC(fit)
      }
      
      det_df$pred_p_oof[ho_idx] <- predict(fit, newdata = det_df[ho_idx, ],
                                           type = "response", na.action = na.pass)
      
      sel[[length(sel) + 1]] <- tibble(region = r, fold = k, selected = chosen, AIC = aic,
                                       n_fit_det = nrow(train_k), n_holdout_det = length(ho_idx))
      cat(sprintf("  %s fold %d/%d: %s (fit n = %d, holdout n = %d) [%.1f min elapsed]\n",
                  r, k, K, chosen, nrow(train_k), length(ho_idx),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      rm(fit)
    }
  }
  list(det = det_df, selection = bind_rows(sel))
}
