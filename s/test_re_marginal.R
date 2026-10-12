# ======================================================================
# Self-check for the random-effect machinery in cvfunc.R.
# Run from s/:  Rscript test_re_marginal.R   (about 1-2 minutes)
# Simulates detections nested in images/blocks/cruises with KNOWN variance
# components, then checks that:
#   1. re_sd() recovers the simulated SDs (and agrees with mgcv::gam.vcomp)
#   2. Gauss-Hermite marginalisation matches Monte Carlo
#   3. On NEW images/blocks, summed marginal p is unbiased for the true
#      positive count, while plain plogis(eta at RE = 0) is biased
#   4. predict_det_p() handles unseen levels (incl. an unseen cruise) and fs terms
# ======================================================================
library(dplyr); library(tidyr); library(mgcv)
source("./cvfunc.R")
set.seed(1)

sim <- function(n_cruise = 5, n_block = 150, img_per_block = 6, det_per_img = 7,
                sd_img = 1.0, sd_block = 0.5, sd_cruise = 0.4, cruise_eff = NULL) {
  if (is.null(cruise_eff)) cruise_eff <- rnorm(n_cruise, 0, sd_cruise)
  blocks <- tibble(transect_block_id = paste0("b", seq_len(n_block)),
                   name_cruise_token = sample(paste0("c", seq_len(n_cruise)), n_block, TRUE),
                   b_block = rnorm(n_block, 0, sd_block))
  imgs <- blocks[rep(seq_len(n_block), each = img_per_block), ] %>%
    mutate(image_id = paste0(transect_block_id, "_i", sequence(rep(img_per_block, n_block))),
           b_img = rnorm(n(), 0, sd_img))
  dets <- imgs[rep(seq_len(nrow(imgs)), each = det_per_img), ] %>%
    mutate(conf = runif(n(), 0.05, 1),
           b_cruise = cruise_eff[as.integer(sub("c", "", name_cruise_token))],
           eta = -1 + 4 * (conf - 0.5) + b_img + b_block + b_cruise,
           p_true = plogis(eta),
           truedetect = rbinom(n(), 1, p_true) == 1)
  list(dets = dets, cruise_eff = cruise_eff)
}

tr <- sim()
te <- sim(cruise_eff = tr$cruise_eff)          # new images & blocks, same cruises
te$dets$image_id <- paste0("new_", te$dets$image_id)
te$dets$transect_block_id <- paste0("new_", te$dets$transect_block_id)

train <- add_re_factors(tr$dets); test <- add_re_factors(te$dets)
f <- truedetect ~ s(conf, k = 5)
re <- RE_SPECS$image_block_cruise

for (eng in c("bam", "gam")) {
  cat("\n=========== engine:", eng, "===========\n")
  fit <- fit_det_gam(with_re(f, re), train, engine = eng)
  
  # 1. Variance components
  info <- re_term_info(fit)
  est <- sapply(info, `[[`, "sd")
  cat("Estimated SDs:", round(est, 3), " (true: image 1.0, block 0.5, cruise ~0.4)\n")
  vc <- tryCatch(gam.vcomp(fit), error = function(e) NULL)
  if (!is.null(vc)) { cat("gam.vcomp for comparison:\n"); print(vc) }
  
  # 2. GH vs Monte Carlo
  eta0 <- c(-2, 0, 1.5); s0 <- c(1.2, 1.2, 0.3)
  mc <- sapply(seq_along(eta0), function(i) mean(plogis(eta0[i] + rnorm(2e5, 0, s0[i]))))
  cat("GH vs MC:", round(marginal_plogis(eta0, s0), 4), "|", round(mc, 4), "\n")
  
  # 3. Unbiased totals on new images/blocks
  p_marg  <- predict_det_p(fit, test)
  p_naive <- predict_det_p(fit, test, marginal = FALSE)
  cat(sprintf("Total TRUE positives: %d | expected (sum p_true): %.0f | marginal: %.0f | naive RE=0: %.0f\n",
              sum(test$truedetect), sum(test$p_true), sum(p_marg), sum(p_naive)))
  
  # 4. Unseen cruise -> integrated out, still close
  test2 <- test; test2$cruise_f <- factor("brand_new_cruise")
  cat(sprintf("Unseen-cruise total (marginal): %.0f\n", sum(predict_det_p(fit, test2))))
}

# fs term (calibration curve varying by cruise)
fit_fs <- fit_det_gam(with_re(f, RE_SPECS$image_cruise_fs), train, engine = "bam")
cat(sprintf("\nfs spec total on new images: %.0f (true TP %d)\n",
            sum(predict_det_p(fit_fs, test)), sum(test$truedetect)))
cat("\nChecks: SDs near truth; GH == MC to ~3 d.p.; marginal total near TP total;\n",
    "naive total visibly lower/higher than truth (that gap is the Jensen bias).\n")
