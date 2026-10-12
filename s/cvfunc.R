library(dplyr)
library(tidyr)
library(mgcv)

# ======================================================================
# Cross-fitting + random-effect utilities for the two-stage
# (detection GAM -> synergy GAM) calibration pipeline.
#
# Cross-fitting: every test_GAM_train image gets an out-of-fold (OOF)
# detection-GAM prediction from a model that never saw it, so the synergy
# GAM is trained on inputs with the same properties as its test inputs.
#   * Folds = whole transect blocks, assigned within region, high-density
#     blocks balanced (assign_folds(), run once in makedat.R).
#   * Formula selection is nested inside each fold (AIC on fitting folds).
#
# Random effects (detection level): detections are nested in images, images
# in transect blocks, blocks in cruises. Grouping terms are added as mgcv
# random-effect smooths:
#   s(image_f, bs = "re")   image-level conditions not captured by covariates
#   s(block_f, bs = "re")   persistent along-transect conditions
#   s(cruise_f, bs = "re")  survey/annotation/camera differences
#   s(conf, cruise_f, bs = "fs")  confidence calibration curve varying by cruise
# Prediction policy (predict_det_p()):
#   * Groups never seen in fitting (always true for held-out images/blocks)
#     are INTEGRATED OUT: p = E[plogis(eta + b)], b ~ N(0, sum of their
#     variances), by Gauss-Hermite quadrature. Predicting at b = 0 instead
#     would give the median, not the mean, probability and bias the summed
#     counts (Jensen's inequality).
#   * Groups listed in RE_CONDITION_ON (default: cruise) are conditioned on
#     when the level was seen in fitting, and integrated out otherwise.
# Specs are compared on OOF scores (compare_re_specs()), not AIC.
# ======================================================================

# ----------------------------------------------------------------------
# 0. Configuration defaults (override in the calling script if needed)
# ----------------------------------------------------------------------
DET_REQUIRED_VARS <- c("truedetect", "conf", "bottom_depth", "latitude", "longitude",
                       "altitude", "field_of_view_sq_meter", "fluorometer_backscatter_ntu",
                       "ctd_temperature_celsius", "boxsize", "iu")

# Listed simplest-first. fs terms use m = 1 (Pedersen et al. 2019, "GS" HGAM):
# per-cruise deviations from the global confidence curve are penalised toward
# flat, which avoids the "repeated 1-d smooths" concurvity warning.
RE_SPECS <- list(
  none               = character(0),
  cruise             = c("s(cruise_f, bs = 're')"),
  cruise_fs          = c("s(conf, cruise_f, bs = 'fs', m = 1, k = 5)"),
  image              = c("s(image_f, bs = 're')"),
  image_block        = c("s(image_f, bs = 're')", "s(block_f, bs = 're')"),
  image_block_cruise = c("s(image_f, bs = 're')", "s(block_f, bs = 're')", "s(cruise_f, bs = 're')"),
  image_cruise_fs    = c("s(image_f, bs = 're')", "s(conf, cruise_f, bs = 'fs', m = 1, k = 5)")
)
RE_CONDITION_ON <- c("cruise_f")   # grouping factors known at prediction time

# "gam" (REML, Laplace) is the default: in the 2026-10-05 run, bam(discrete = TRUE)
# returned non-converged fits (negative deviance explained) for the tensor-product
# formulas, and in the simulation self-check it under-estimated the image SD
# (0.885 vs 1.0; gam gave 1.007), which under-corrects the marginal sums.
# bam remains available for image/block REs where gam is too slow.
if (!exists("DET_ENGINE"))   DET_ENGINE   <- "gam"
if (!exists("DET_DISCRETE")) DET_DISCRETE <- FALSE   # only used when engine = "bam"
if (!exists("DET_THREADS"))  DET_THREADS  <- 1       # macOS R has no OpenMP

det_fit_rows <- function(det_df) {
  det_df %>% drop_na(all_of(DET_REQUIRED_VARS))
}

# Grouping factors used by the random-effect terms. Tolerates missing source
# columns (e.g. inference data without blocks): the factor is then NA, which
# predict_det_p() treats as an unseen level and integrates out.
add_re_factors <- function(df) {
  mk <- function(col) if (col %in% names(df)) factor(as.character(df[[col]])) else factor(rep(NA_character_, nrow(df)))
  df$image_f  <- mk("image_id")
  df$block_f  <- mk("transect_block_id")
  df$cruise_f <- mk("name_cruise_token")
  df
}

with_re <- function(formula, re_terms) {
  if (length(re_terms) == 0) return(formula)
  update(formula, as.formula(paste(". ~ . +", paste(re_terms, collapse = " + "))))
}

# re_terms may be a character vector (both regions) or list(GB = , MAB = )
re_for_region <- function(re_terms, region) {
  if (is.list(re_terms)) re_terms[[region]] else re_terms
}

fit_det_gam <- function(formula, data, engine = DET_ENGINE, nthreads = DET_THREADS,
                        discrete = DET_DISCRETE) {
  data <- droplevels(data)
  if (engine == "bam") {
    withCallingHandlers(
      bam(formula, family = binomial(), data = data, method = "fREML",
          discrete = discrete, nthreads = nthreads),
      warning = function(w) if (grepl("openMP", conditionMessage(w))) invokeRestart("muffleWarning"))
  } else {
    gam(formula, family = binomial(), data = data, method = "REML")
  }
}

# Deviance explained without summary() (summary.gam's RE tests can fail with a
# Cholesky error on large RE fits, which crashed the 2026-10-05 run).
dev_expl <- function(fit) 1 - fit$deviance / fit$null.deviance

# A fit is usable only if it converged, beats the intercept-only model, and
# raised no non-convergence warning. The warning check matters: check_engine.R
# showed bam(discrete = TRUE) returning converged = TRUE (and sometimes a
# positive deviance explained) for tensor-product fits whose iterations had in
# fact diverged ("algorithm did not converge", "fitted probabilities 0 or 1").
fit_ok <- function(fit) {
  !is.null(fit) && isTRUE(fit[["converged"]]) && !isTRUE(fit[["fit_warned"]]) &&
    is.finite(dev_expl(fit)) && dev_expl(fit) >= 0
}

# Warnings that make a fit unusable. For bam these include the separation
# message seen in check_engine.R. For gam only a PIRLS divergence counts:
# mgcv's outer-iteration notes ("step failure", "iteration limit") are common
# and usually benign (e.g. a smooth penalised to a straight line), so they are
# recorded in fit$fit_notes but do not reject the fit.
FIT_REJECT_PATTERN <- list(bam = "did not converge|numerically 0 or 1",
                           gam = "algorithm did not converge")

try_fit_det_gam <- function(formula, data, engine = DET_ENGINE, label = "", ...) {
  msgs <- character(0)
  fit <- tryCatch(
    withCallingHandlers(fit_det_gam(formula, data, engine = engine, ...),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
      }),
    error = function(e) { cat("  !! Fit failed", label, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)
  msgs <- unique(msgs)
  fit[["fit_notes"]]  <- paste(msgs, collapse = " | ")
  fit[["fit_warned"]] <- any(grepl(FIT_REJECT_PATTERN[[engine]], msgs))
  if (!fit_ok(fit)) {
    # cat() rather than warning(): shows up live in the log, even if R halts later
    cat(sprintf("  !! Unusable fit %s (converged = %s, dev. expl. = %.3f)%s\n", label,
                isTRUE(fit[["converged"]]), dev_expl(fit),
                if (length(msgs)) paste0(" warnings: ", fit[["fit_notes"]]) else ""))
  }
  fit
}

# ----------------------------------------------------------------------
# 1. Fold assignment (block-grouped, region-specific, density-balanced)
# ----------------------------------------------------------------------
assign_folds <- function(meta_df, K = 5, seed = 2226, fold_set = "test_GAM_train",
                         block_col = "transect_block_id") {
  
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
    mutate(fold = if_else(dataset == fold_set, fold, NA_integer_))
  
  cat("Fold balance (images per fold, by region x density bin):\n")
  print(out %>% filter(!is.na(fold)) %>% count(region, density_bin, fold) %>%
          pivot_wider(names_from = fold, values_from = n, names_prefix = "fold_"))
  out
}

# ----------------------------------------------------------------------
# 2. Random-effect bookkeeping and marginal prediction
# ----------------------------------------------------------------------
# Standard deviation of a bs = "re" term: sigma^2 = scale * S.scale / sp
# (same quantity mgcv::gam.vcomp reports; binomial scale = 1).
re_sd <- function(fit, sm) {
  sp <- fit$sp[sm$first.sp]
  s_scale <- if (is.null(sm$S.scale)) 1 else sm$S.scale
  as.numeric(sqrt(fit$sig2 * s_scale / sp))
}

re_term_info <- function(fit) {
  out <- list()
  for (sm in fit$smooth) {
    is_re <- inherits(sm, "random.effect")
    is_fs <- inherits(sm, "fs.interaction")
    if (!is_re && !is_fs) next
    fac <- if (is_re) sm$term[length(sm$term)] else sm$fterm
    out[[sm$label]] <- list(label = sm$label, factor = fac,
                            type = if (is_re) "re" else "fs",
                            sd = if (is_re) re_sd(fit, sm) else NA_real_)
  }
  out
}

# Gauss-Hermite rule (Golub-Welsch); E[f(Z)] ~ sum(w * f(sqrt(2) x)) / sqrt(pi)
gh_rule <- function(n = 20) {
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i / 2)
  J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}
GH_RULE <- gh_rule(20)

marginal_plogis <- function(eta, sd) {
  out <- plogis(eta)
  idx <- which(!is.na(eta) & sd > 0)
  if (length(idx) > 0) {
    z <- outer(eta[idx], rep(1, length(GH_RULE$x))) + outer(sd[idx], sqrt(2) * GH_RULE$x)
    out[idx] <- as.vector(plogis(z) %*% GH_RULE$w) / sqrt(pi)
  }
  out
}

# Link-scale prediction; falls back to non-discrete prediction for bam fits
# if discrete prediction rejects the request.
pred_link <- function(fit, nd, excl) {
  args <- list(fit, newdata = nd, type = "link",
               exclude = if (length(excl)) excl else NULL, na.action = na.pass)
  out <- tryCatch(do.call(predict, args), error = function(e) {
    if (inherits(fit, "bam")) do.call(predict, c(args, list(discrete = FALSE))) else stop(e)
  })
  as.numeric(out)
}

# Population-level p(true detection) for new data under the prediction policy
# described at the top of this file. Works for models with or without REs.
predict_det_p <- function(fit, newdata, condition_on = RE_CONDITION_ON, marginal = TRUE) {
  info <- re_term_info(fit)
  if (length(info) == 0) {
    return(as.numeric(predict(fit, newdata = newdata, type = "response", na.action = na.pass)))
  }
  if (!all(c("image_f", "block_f", "cruise_f") %in% names(newdata))) newdata <- add_re_factors(newdata)
  
  nd <- newdata
  n <- nrow(nd)
  excl <- character(0)          # terms integrated out for every row
  extra_var <- numeric(n)
  unknown <- list()             # conditioned terms with unseen levels, by row
  
  for (ti in info) {
    v   <- ti$factor
    lev <- levels(fit$model[[v]])
    x   <- as.character(nd[[v]])
    known <- !is.na(x) & x %in% lev
    x[!known] <- lev[1]                      # placeholder; term excluded below
    nd[[v]] <- factor(x, levels = lev)
    if (v %in% condition_on) {
      if (any(!known)) unknown[[ti$label]] <- !known
    } else {
      excl <- c(excl, ti$label)
      if (ti$type == "re") extra_var <- extra_var + ti$sd^2
    }
  }
  
  eta <- pred_link(fit, nd, excl)
  
  if (length(unknown) > 0) {
    unk <- do.call(cbind, unknown)
    pattern <- apply(unk, 1, paste, collapse = "")
    for (pt in unique(pattern[rowSums(unk) > 0])) {
      rows <- which(pattern == pt)
      labs <- colnames(unk)[unk[rows[1], ]]
      eta[rows] <- pred_link(fit, nd[rows, , drop = FALSE], c(excl, labs))
      for (l in labs) {
        if (info[[l]]$type == "re") {
          extra_var[rows] <- extra_var[rows] + info[[l]]$sd^2
        } else {
          warning("Unseen level for fs term ", l, ": term dropped without variance correction for ", length(rows), " rows.")
        }
      }
    }
  }
  
  if (marginal) marginal_plogis(eta, sqrt(extra_var)) else plogis(eta)
}

re_sd_table <- function(fit) {
  info <- re_term_info(fit)
  sds <- vapply(info, function(ti) ti$sd, numeric(1))
  sds <- sds[!is.na(sds)]
  if (length(sds) == 0) return(tibble())
  as_tibble(as.list(setNames(sds, paste0("sd_", gsub("^s\\(|\\)$", "", names(sds))))))
}

# ----------------------------------------------------------------------
# 3. Fit candidate formulas (+ RE terms), return the AIC table and best fit
# ----------------------------------------------------------------------
fit_best_detection_gam <- function(train_data, candidate_formulas, re_terms = character(0),
                                   engine = DET_ENGINE, verbose = FALSE) {
  best_fit <- NULL; best_aic <- Inf; best_id <- NA_character_
  tab <- vector("list", length(candidate_formulas))
  for (j in seq_along(candidate_formulas)) {
    id <- names(candidate_formulas)[j]
    if (verbose) cat("    Fitting:", id, "\n")
    fit <- try_fit_det_gam(with_re(candidate_formulas[[j]], re_terms), train_data, engine = engine, label = id)
    ok <- fit_ok(fit)
    a  <- if (is.null(fit)) NA_real_ else AIC(fit)
    tab[[j]] <- tibble(Model_ID = id, AIC = a,
                       Deviance_Explained = if (is.null(fit)) NA_real_ else dev_expl(fit),
                       edf = if (is.null(fit)) NA_real_ else sum(fit$edf), usable = ok)
    if (ok && a < best_aic) { best_aic <- a; best_fit <- fit; best_id <- id }
    rm(fit)
  }
  if (is.null(best_fit)) stop("No usable candidate fit.")
  list(fit = best_fit, id = best_id, table = bind_rows(tab) %>% arrange(AIC))
}

pick_best_formula <- function(selection_table, candidate_formulas) {
  best_id <- selection_table$Model_ID[which.min(selection_table$AIC)]
  cat("  [Selected]", best_id, "\n")
  candidate_formulas[[best_id]]
}

# ----------------------------------------------------------------------
# 4. Cross-fit the detection GAMs
# ----------------------------------------------------------------------
# Returns list(det = det_df + pred_p_oof [+ pred_p_foldavg], selection = per-fold table)
#   pred_p_oof     : test_GAM_train rows, from the model that did not see them
#   pred_p_foldavg : test_GAM_test rows, mean of the K fold models' predictions
crossfit_detection_gams <- function(det_df, candidate_formulas = NULL,
                                    fixed_formulas = NULL, K = 5,
                                    nested_selection = TRUE,
                                    re_terms = character(0),
                                    engine = DET_ENGINE,
                                    predict_test = TRUE,
                                    fold_set = "test_GAM_train",
                                    test_set = "test_GAM_test") {
  if (!"fold" %in% names(det_df)) {
    stop("No `fold` column in det_df. Re-run makedat.R (assign_folds) to rebuild eval_2226.rds.")
  }
  if (nested_selection && is.null(candidate_formulas)) stop("candidate_formulas required for nested selection.")
  if (!nested_selection && is.null(fixed_formulas)) stop("fixed_formulas required when nested_selection = FALSE.")
  
  det_df <- add_re_factors(det_df)
  det_df$pred_p_oof <- NA_real_
  if (predict_test) det_df$pred_p_foldavg <- NA_real_
  sel <- list()
  t0 <- Sys.time()
  
  for (r in c("GB", "MAB")) {
    re_r <- re_for_region(re_terms, r)
    te_idx <- which(det_df$dataset == test_set & det_df$region == r)
    te_acc <- numeric(length(te_idx)); te_n <- 0
    
    for (k in seq_len(K)) {
      in_set  <- det_df$dataset == fold_set & det_df$region == r & !is.na(det_df$fold)
      fit_idx <- which(in_set & det_df$fold != k)
      ho_idx  <- which(in_set & det_df$fold == k)
      if (length(ho_idx) == 0) next
      
      train_k <- det_fit_rows(det_df[fit_idx, ])
      
      if (nested_selection) {
        best <- fit_best_detection_gam(train_k, candidate_formulas, re_terms = re_r, engine = engine)
        fit <- best$fit; chosen <- best$id; aic <- min(best$table$AIC)
      } else {
        fit <- try_fit_det_gam(with_re(fixed_formulas[[r]], re_r), train_k, engine = engine,
                               label = sprintf("%s fold %d", r, k))
        if (is.null(fit)) stop(sprintf("Fit failed in %s fold %d", r, k))
        chosen <- "fixed"; aic <- AIC(fit)
      }
      
      det_df$pred_p_oof[ho_idx] <- predict_det_p(fit, det_df[ho_idx, ])
      if (predict_test && length(te_idx) > 0) {
        te_acc <- te_acc + predict_det_p(fit, det_df[te_idx, ]); te_n <- te_n + 1
      }
      
      # NB: bind only non-empty SD tables. bind_cols() recycles to the smallest
      # size, so binding a 0-row tibble (models without bs = "re" terms) silently
      # dropped the whole row - which emptied the selection table for "none" and
      # "cruise_fs" fits in the 2026-10-05 evening run.
      sd_tab <- re_sd_table(fit)
      sel_row <- tibble(region = r, fold = k, selected = chosen, AIC = aic,
               edf = sum(fit$edf), usable = fit_ok(fit),
               converged = isTRUE(fit[["converged"]]), dev_expl = dev_expl(fit),
               notes = if (is.null(fit[["fit_notes"]])) "" else fit[["fit_notes"]],
               n_fit_det = nrow(train_k), n_holdout_det = length(ho_idx))
      sel[[length(sel) + 1]] <- if (nrow(sd_tab) == 1) bind_cols(sel_row, sd_tab) else sel_row
      cat(sprintf("  %s fold %d/%d: %s (fit n = %d, holdout n = %d) [%.1f min elapsed]\n",
                  r, k, K, chosen, nrow(train_k), length(ho_idx),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      rm(fit)
    }
    if (predict_test && te_n > 0) det_df$pred_p_foldavg[te_idx] <- te_acc / te_n
  }
  list(det = det_df, selection = bind_rows(sel))
}

# ----------------------------------------------------------------------
# 5. OOF scoring and random-effect structure comparison
# ----------------------------------------------------------------------
r2_1to1 <- function(y, yhat) {
  ok <- !is.na(y) & !is.na(yhat)
  1 - sum((y[ok] - yhat[ok])^2) / sum((y[ok] - mean(y[ok]))^2)
}

# Stage-1 target is the TRUE-POSITIVE count per image (the detection GAM
# estimates which detections are real; missed scallops are stage 2's job).
score_oof_detection <- function(det_df, p_col = "pred_p_oof", set = "test_GAM_train") {
  d <- det_df %>%
    filter(dataset == set, !is.na(.data[[p_col]])) %>%
    mutate(p = pmin(pmax(.data[[p_col]], 1e-6), 1 - 1e-6), y = as.numeric(truedetect))
  det_s <- d %>% group_by(region) %>%
    summarise(log_loss = -mean(y * log(p) + (1 - y) * log(1 - p)),
              brier = mean((y - p)^2), n_det = n(), .groups = "drop")
  img_s <- d %>% group_by(region, image_id) %>%
    summarise(tp = sum(y), sp = sum(p), .groups = "drop") %>%
    group_by(region) %>%
    summarise(rmse_tp = sqrt(mean((sp - tp)^2)), r2_1to1_tp = r2_1to1(tp, sp),
              total_tp = sum(tp), total_pred = sum(sp),
              pct_bias_total = 100 * (sum(sp) - sum(tp)) / sum(tp), .groups = "drop")
  left_join(det_s, img_s, by = "region")
}

# Per-block OOF log loss (sum over detections) for paired, cluster-robust
# comparisons between candidates.
oof_block_loss <- function(det_df, p_col = "pred_p_oof", set = "test_GAM_train") {
  det_df %>%
    filter(dataset == set, !is.na(.data[[p_col]])) %>%
    mutate(p = pmin(pmax(.data[[p_col]], 1e-6), 1 - 1e-6), y = as.numeric(truedetect),
           ll = -(y * log(p) + (1 - y) * log(1 - p))) %>%
    group_by(region = as.character(region), block = as.character(transect_block_id)) %>%
    summarise(loss = sum(ll), n = n(), .groups = "drop")
}

# One-standard-error rule. For each region: find the candidate with the lowest
# OOF log loss; for every other candidate compute the paired per-detection loss
# difference and its block-clustered SE; keep candidates within 1 SE of the
# best and pick the least complex (lowest mean edf). Small differences that are
# noise therefore resolve to the simpler model instead of the noisiest winner.
#   scores: tibble(candidate, region, log_loss, complexity)
#   block_losses: named list candidate -> oof_block_loss() output
choose_one_se <- function(scores, block_losses) {
  out <- list(); tab <- list()
  for (r in unique(as.character(scores$region))) {
    sr_all <- scores %>% filter(as.character(region) == r)
    sr <- sr_all %>% filter(coalesce(all_usable, FALSE))
    if (nrow(sr) == 0) {
      cat("  !! No candidate had usable fits in every fold for", r,
          "- falling back to all candidates. Check the `usable`/`notes` columns.\n")
      sr <- sr_all
    }
    best <- sr$candidate[which.min(sr$log_loss)]
    bb <- block_losses[[best]] %>% filter(region == r)
    rows <- lapply(sr$candidate, function(cand) {
      cb <- block_losses[[cand]] %>% filter(region == r)
      j  <- inner_join(cb, bb, by = c("region", "block"), suffix = c("", "_best"))
      N  <- sum(j$n); B <- nrow(j)
      d  <- j$loss - j$loss_best
      tibble(candidate = cand, region = r,
             d_logloss = sum(d) / N,
             se = if (B > 1) sqrt(B * var(d)) / N else NA_real_)
    })
    t <- bind_rows(rows) %>% left_join(sr, by = c("candidate", "region")) %>%
      mutate(within_1se = d_logloss <= coalesce(se, 0))
    pick <- t %>% filter(within_1se) %>% arrange(complexity, log_loss) %>% slice(1)
    out[[r]] <- pick$candidate
    tab[[r]] <- t %>% mutate(best = candidate == best, chosen = candidate == pick$candidate)
  }
  list(chosen = out, table = bind_rows(tab))
}

# Cross-validated selection of the FIXED-EFFECT formula (no REs by default):
# every candidate is cross-fitted on the block folds and scored on OOF log
# loss; choice by the one-SE rule with mean edf as complexity. Because the
# folds are whole transect blocks, this respects the clustering that makes
# AIC over-favour complex formulas when detections are treated as independent.
cv_select_detection_formula <- function(det_df, candidate_formulas, K = 5,
                                        re_terms = character(0), engine = DET_ENGINE) {
  scores <- list(); bl <- list(); diag <- list()
  for (id in names(candidate_formulas)) {
    cat("\n-- Formula:", id, "--\n")
    cf <- tryCatch(
      crossfit_detection_gams(det_df, fixed_formulas = list(GB = candidate_formulas[[id]],
                                                            MAB = candidate_formulas[[id]]),
                              K = K, nested_selection = FALSE, re_terms = re_terms,
                              engine = engine, predict_test = FALSE),
      error = function(e) { warning("Formula ", id, " failed: ", conditionMessage(e)); NULL })
    if (is.null(cf)) next
    diag[[id]] <- cf$selection %>% mutate(candidate = id, .before = 1)
    ux <- cf$selection %>% group_by(region) %>%
      summarise(all_usable = all(usable), complexity = mean(edf), .groups = "drop")
    sc <- score_oof_detection(cf$det) %>% mutate(region = as.character(region)) %>%
      left_join(ux, by = "region") %>% mutate(candidate = id, .before = 1)
    scores[[id]] <- sc
    bl[[id]] <- oof_block_loss(cf$det)
  }
  scores <- bind_rows(scores)
  list(scores = scores, block_losses = bl, fold_diag = bind_rows(diag))
}

# Applies the one-SE rule to cv_select_detection_formula() output (kept
# separate so cached CV results can be re-scored without refitting).
choose_detection_formula <- function(cv) {
  pick <- choose_one_se(cv$scores, cv$block_losses)
  c(cv, list(one_se = pick$table, chosen = pick$chosen))
}

# Cross-fit a fixed formula under each RE spec and score OOF; choice by the
# one-SE rule (complexity = mean edf, so "none" wins ties).
compare_re_specs <- function(det_df, base_formulas, re_specs = RE_SPECS, K = 5,
                             engine = DET_ENGINE) {
  scores <- list(); sds <- list(); bl <- list()
  for (nm in names(re_specs)) {
    cat("\n-- RE spec:", nm, "--\n")
    cf <- tryCatch(
      crossfit_detection_gams(det_df, fixed_formulas = base_formulas, K = K,
                              nested_selection = FALSE, re_terms = re_specs[[nm]],
                              engine = engine, predict_test = FALSE),
      error = function(e) { warning("RE spec ", nm, " failed: ", conditionMessage(e)); NULL })
    if (is.null(cf)) next
    cx <- cf$selection %>% group_by(region) %>%
      summarise(all_usable = all(usable), complexity = mean(edf), .groups = "drop")
    scores[[nm]] <- score_oof_detection(cf$det) %>% mutate(region = as.character(region)) %>%
      left_join(cx, by = "region") %>% mutate(candidate = nm, .before = 1)
    sds[[nm]] <- cf$selection %>% mutate(re_spec = nm, .before = 1)
    bl[[nm]]  <- oof_block_loss(cf$det)
  }
  scores <- bind_rows(scores)
  pick <- choose_one_se(scores, bl)
  list(scores = scores, one_se = pick$table, sds = bind_rows(sds), chosen = pick$chosen)
}

# ----------------------------------------------------------------------
# 6. Variance-component report (inference, not prediction)
# ----------------------------------------------------------------------
# Fits the chosen fixed formula + image/block/cruise REs once per region on
# all of test_GAM_train and reports SDs with 95% CIs (mgcv::gam.vcomp). These
# document how strongly detections cluster within images/blocks/cruises, i.e.
# why a no-RE model's AIC/SEs are over-confident, even when conditioning on
# image/block does not improve predictions for NEW images.
# Uses bam(discrete = TRUE) by default: gam with ~3,000 image levels is very
# slow. bam's PQL-type estimates run low for binary data (simulation: image SD
# 0.885 vs true 1.0), so these SDs are conservative (understate clustering).
# bam(discrete) diverges on tensor-product formulas (check_engine.R), so pass an
# ADDITIVE fixed part here (gamstrat.R uses VARCOMP_FORMULA, default M05_Boxsize).
report_variance_components <- function(train_det, formulas, re_terms = RE_SPECS$image_block_cruise,
                                       engine = "bam", discrete = TRUE) {
  purrr::map_dfr(c("GB", "MAB"), function(r) {
    d <- train_det %>% filter(region == r)
    fit <- try_fit_det_gam(with_re(formulas[[r]], re_terms), d, engine = engine,
                           label = paste("varcomp", r), discrete = discrete)
    if (is.null(fit)) return(tibble(region = r))
    vc <- tryCatch({ utils::capture.output(v <- gam.vcomp(fit)); v }, error = function(e) NULL)
    if (is.list(vc) && !is.matrix(vc) && !is.null(vc$vc)) vc <- vc$vc
    if (is.null(vc) || !is.matrix(vc)) {
      sds <- vapply(re_term_info(fit), function(ti) ti$sd, numeric(1))
      return(tibble(region = r, term = names(sds), std_dev = sds, lower = NA_real_, upper = NA_real_))
    }
    keep <- rownames(vc) %in% names(re_term_info(fit))
    tibble(region = r, term = rownames(vc)[keep], std_dev = vc[keep, 1],
           lower = vc[keep, 2], upper = vc[keep, 3], engine = engine,
           n_det = nrow(d), dev_expl = dev_expl(fit))
  })
}
