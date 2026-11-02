# ======================================================================
# plot_final_gams.R - what the recommended two-stage YOLOv12 GAMs actually do
#   Stage 1: detection confidence -> probability the detection is a scallop
#            (other covariates held at their medians), with binned observed TP rates
#   Stage 2: summed probabilities -> calibrated count per image
#            (GB: at low / median / high share of low-confidence detections;
#             MAB: other terms at their medians), with Calibration-fitting images
# Also writes the term tables (edf, p) for a slide.
# Run from s/:  Rscript n3/plot_final_gams.R
# Needs: ../data/processed/YOLOv12_detgams_2226.rds (gamstrat.R)
#        ../data/processed/n3/n3_results.rds        (n3_ablation.R)
# Writes ../figures/ms_figures/2226_final_gams_yolo.pdf|png, ../data/processed/n3/final_gam_terms.csv
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(mgcv); library(patchwork) })
FIG <- "../figures/ms_figures"; OUT <- "../data/processed/n3"
REG <- c(GB = "Georges Bank", MAB = "Mid-Atlantic Bight")

# ---------------- Stage 1 ----------------
s1 <- readRDS("../data/processed/YOLOv12_detgams_2226.rds")
s1_curve <- bind_rows(lapply(names(REG), function(r) {
  g <- s1[[r]]; mf <- g$model
  base <- mf %>% summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))
  nd <- as.data.frame(base)[rep(1, 200), , drop = FALSE] %>% mutate(conf = seq(0.1, 1, length.out = 200))
  pr <- predict(g, nd, type = "link", se.fit = TRUE)
  tibble(region = REG[[r]], conf = nd$conf, p = plogis(pr$fit),
         lo = plogis(pr$fit - 1.96 * pr$se.fit), hi = plogis(pr$fit + 1.96 * pr$se.fit),
         form = paste(deparse(formula(g), width.cutoff = 500), collapse = ""))
}))
s1_obs <- bind_rows(lapply(names(REG), function(r) {
  mf <- s1[[r]]$model
  y <- mf[[1]]; y <- if (is.logical(y)) as.numeric(y) else as.numeric(as.character(y))
  tibble(conf = mf$conf, y = y) %>% mutate(bin = cut(conf, seq(0.1, 1, by = 0.05), include.lowest = TRUE)) %>%
    group_by(bin) %>% summarise(conf = mean(conf), rate = mean(y), n = n(), .groups = "drop") %>%
    mutate(region = REG[[r]])
}))
p1 <- ggplot(s1_curve, aes(conf, p)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#7570B3", alpha = 0.25) +
  geom_line(colour = "#7570B3", linewidth = 0.9) +
  geom_point(data = s1_obs, aes(conf, rate, size = n), colour = "grey30", alpha = 0.7) +
  facet_wrap(~ region) + scale_size_area(max_size = 3, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) +
  labs(x = "YOLOv12 confidence", y = "P(detection is a scallop)",
       title = "Stage 1: confidence to detection probability",
       subtitle = "Line: GAM, other covariates at their medians (95% CI). Points: observed share of true positives per 0.05 confidence bin.") +
  theme(panel.grid.minor = element_blank())

# ---------------- Stage 2 ----------------
n3 <- readRDS(file.path(OUT, "n3_results.rds"))
cal <- n3$img %>% filter(dataset == "test_GAM_train")
F2 <- list(GB  = n_annotations ~ te(l_Y, lowconf_frac_Y, k = c(5, 4)),
           MAB = n_annotations ~ s(l_Y) + s(med_box_mm_Y, k = 5) + s(lowconf_frac_Y, k = 5) +
                                 s(overlap_frac_Y, k = 5) + s(nn_ratio_Y, k = 5))
s2 <- lapply(names(F2), function(r) gam(F2[[r]], family = tw(link = "log"), method = "REML",
                                        data = cal %>% filter(region == r))); names(s2) <- names(F2)
# Lines: GAM prediction across the OBSERVED range of Sigma p, with 95% CI.
#   GB: one line per tercile of the low-confidence share (images with Sigma p >= 1),
#       each drawn only over the Sigma p range of the images in that tercile (no extrapolation).
#   MAB: box size, low-confidence share, overlap and spacing at their medians
#       (images with detections); these four terms are near-linear and non-significant.
pred_ci <- function(g, nd) {
  pr <- predict(g, nd, type = "link", se.fit = TRUE)
  nd %>% mutate(fit = exp(pr$fit), lo = exp(pr$fit - 1.96 * pr$se.fit), hi = exp(pr$fit + 1.96 * pr$se.fit))
}
s2_curve_list <- lapply(setNames(names(REG), names(REG)), function(r) {
  d <- cal %>% filter(region == r); g <- s2[[r]]
  if (r == "GB") {
    dd <- d %>% filter(pred_Y >= 1) %>% mutate(ter = dplyr::ntile(lowconf_frac_Y, 3))
    bind_rows(lapply(1:3, function(k) {
      dk <- dd %>% filter(ter == k); v <- median(dk$lowconf_frac_Y)
      sp <- seq(0, max(dk$pred_Y), length.out = 200)
      pred_ci(g, tibble(pred_Y = sp, l_Y = log1p(sp), lowconf_frac_Y = v)) %>%
        mutate(level = sprintf("%s low-confidence share (%.2f)", c("Low", "Middle", "High")[k], v))
    }))
  } else {
    base <- d %>% filter(pred_Y > 0) %>% summarise(across(c(med_box_mm_Y, lowconf_frac_Y, overlap_frac_Y, nn_ratio_Y), median))
    sp <- seq(0, max(d$pred_Y), length.out = 200)
    pred_ci(g, as.data.frame(base)[rep(1, length(sp)), , drop = FALSE] %>% mutate(pred_Y = sp, l_Y = log1p(sp))) %>%
      mutate(level = "Other terms at their medians")
  }
})
s2_curve <- bind_rows(lapply(names(REG), function(r) {
  x <- s2_curve_list[[r]]; x$region <- REG[[r]]; x }))
pts <- cal %>% mutate(region = unname(REG[as.character(region)]))
lv <- unique(s2_curve$level)
cols <- setNames(c("#FB6A4A", "#7570B3", "#3F007D", "#2171B5")[seq_along(lv)], lv)
p2 <- ggplot(s2_curve, aes(pred_Y, fit)) +
  geom_point(data = pts, aes(pred_Y, n_annotations), colour = "grey60", size = 0.6, alpha = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = level), alpha = 0.15, colour = NA) +
  geom_line(aes(colour = level), linewidth = 0.9) +
  facet_wrap(~ region, scales = "free") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 25, 50, 100, 250)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 25, 50, 100, 250)) +
  scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, guide = "none") +
  theme_bw(base_size = 11) +
  labs(x = expression("Stage 1 expected count per image (" * Sigma * p[det] * ")"),
       y = "Final count (line) / manual count (points)",
       title = "Stage 2: stage 1 expected count to final calibrated count",
       subtitle = "Tweedie GAM, log link, 95% CI; log(1 + x) axes. Dashed 1:1 = no change from stage 1; above it the model adds missed scallops, below it removes false ones.") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

fig <- p1 / p2 + plot_layout(heights = c(1, 1.15))
ggsave(file.path(FIG, "2226_final_gams_yolo.pdf"), fig, width = 10, height = 9)
ggsave(file.path(FIG, "2226_final_gams_yolo.png"), fig, width = 10, height = 9, dpi = 300, bg = "white")

# ---------------- Term tables ----------------
terms <- bind_rows(
  bind_rows(lapply(names(REG), function(r) { st <- summary(s1[[r]])$s.table
    tibble(stage = "1 (binomial, logit)", region = r, term = rownames(st), edf = st[, "edf"], p = st[, ncol(st)],
           dev_expl = summary(s1[[r]])$dev.expl) })),
  bind_rows(lapply(names(REG), function(r) { st <- summary(s2[[r]])$s.table
    tibble(stage = "2 (Tweedie, log)", region = r, term = rownames(st), edf = st[, "edf"], p = st[, ncol(st)],
           dev_expl = summary(s2[[r]])$dev.expl) })))
write.csv(terms, file.path(OUT, "final_gam_terms.csv"), row.names = FALSE)
print(terms %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = Inf, width = Inf)
cat("\nStage 1 formulas:\n"); for (r in names(REG)) cat(" ", r, ":", deparse(formula(s1[[r]]), width.cutoff = 500), "\n")
