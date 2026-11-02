# ======================================================================
# n7_density.R - N7 (R2.6): scallop density (per m^2) for every workflow
#
# Every comparison in the paper is against absolute manual counts per image.
# Here counts become densities with each image's field of view (m^2):
#   area-weighted density  D = sum(count) / sum(FOV)    (the survey ratio estimator)
#   per-image density      d_i = count_i / FOV_i        (RMSE in scallops per m^2)
# by region, and by region x survey year for the three workflows a
# stock assessment would use. 95% intervals: transect-block bootstrap.
# For the area-weighted estimator the % bias equals the count bias by
# construction (same FOV in numerator and denominator); the table adds the
# absolute densities and their uncertainty, which is what R2.6 asks for.
#
# Run from s/:  Rscript n7_density.R
# Needs: ../data/processed/holdout_eval_2226_v2.rds (n3/n3_export_holdout.R)
# Writes ../data/processed/n7_density/{density_region.csv, density_region_year.csv}
#        tables/tableS6_density.tex
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr) })
OUT <- "../data/processed/n7_density"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
B <- 1000; set.seed(2226)
pr <- function(x) print(tibble::as_tibble(x), n = Inf, width = Inf)

M <- c("Manual count" = "n_annotations",
       "YOLOv12 F1 threshold" = "yolo_f1", "Cascade R-CNN F1 threshold" = "cascade_f1",
       "YOLOv12 stage 1 only" = "pred_yolo",
       "Two-stage YOLOv12" = "pred_2s_yolo", "Two-stage YOLOv12 + Cascade R-CNN" = "pred_2s_yc")
ev <- readRDS("../data/processed/holdout_eval_2226_v2.rds") %>%
  mutate(region = as.character(region),
         fov = if_else(field_of_view_sq_meter <= 0 | field_of_view_sq_meter <= -90, NA_real_, field_of_view_sq_meter),
         year = substr(as.character(name_cruise_token), 1, 4))
n_drop <- sum(is.na(ev$fov)); if (n_drop) cat("Images without a usable field of view (left out):", n_drop, "\n")
ev <- ev %>% filter(!is.na(fov))

dens_stats <- function(d) {
  a <- sum(d$fov)
  map_dfr(names(M), function(m) {
    y <- d[[M[[m]]]]
    tibble(method = m, D = sum(y) / a, mean_d = mean(y / d$fov),
           rmse_d = sqrt(mean((y / d$fov - d$n_annotations / d$fov)^2)))
  })
}
boot_dens <- function(d) {
  pt <- dens_stats(d)
  idx <- split(seq_len(nrow(d)), d$transect_block_id); nb <- length(idx)
  bs <- map_dfr(seq_len(B), function(b) dens_stats(d[unlist(idx[sample.int(nb, nb, TRUE)], use.names = FALSE), ]) %>%
                  mutate(b = b))
  # bias of each method's D relative to the manual D in the same resample
  bs <- bs %>% group_by(b) %>% mutate(bias = 100 * (D / D[method == "Manual count"] - 1)) %>% ungroup()
  ci <- bs %>% group_by(method) %>%
    summarise(D_lo = quantile(D, .025), D_hi = quantile(D, .975),
              bias_lo = quantile(bias, .025), bias_hi = quantile(bias, .975), .groups = "drop")
  pt %>% mutate(bias = 100 * (D / D[method == "Manual count"] - 1)) %>%
    left_join(ci, by = "method") %>%
    mutate(n_images = nrow(d), area_m2 = sum(d$fov), scallops = sum(d$n_annotations))
}

by_region <- ev %>% group_by(region) %>% group_modify(~ boot_dens(.x)) %>% ungroup() %>%
  mutate(method = factor(method, levels = names(M))) %>% arrange(region, method)
write.csv(by_region, file.path(OUT, "density_region.csv"), row.names = FALSE)
cat("\nDensity by region (scallops per m^2; bias vs manual, %):\n")
pr(by_region %>% mutate(across(c(D, D_lo, D_hi, mean_d, rmse_d), ~ round(.x, 4)), across(starts_with("bias"), ~ round(.x, 1))))

M_YR <- c("Manual count", "YOLOv12 F1 threshold", "Two-stage YOLOv12")
# No Mid-Atlantic Bight images were annotated in 2023: those cruise tokens are mislabelled
# (2022 or 2024), so they stay in the regional rows but get no survey-year row.
by_year <- ev %>% filter(!(region == "MAB" & year == "2023")) %>%
  group_by(region, year) %>% filter(sum(n_annotations) >= 50) %>%
  group_modify(~ boot_dens(.x) %>% filter(method %in% M_YR)) %>% ungroup() %>%
  mutate(method = factor(method, levels = M_YR)) %>% arrange(region, year, method)
write.csv(by_year, file.path(OUT, "density_region_year.csv"), row.names = FALSE)
cat("\nDensity by region and survey year (years with >= 50 Holdout scallops):\n")
pr(by_year %>% select(region, year, method, n_images, scallops, D, D_lo, D_hi, bias, bias_lo, bias_hi) %>%
     mutate(across(c(D, D_lo, D_hi), ~ round(.x, 4)), across(starts_with("bias"), ~ round(.x, 1))))

# ---------------- LaTeX table (supplement) ----------------
f3 <- function(x) formatC(x, format = "f", digits = 3)
f1s <- function(x) gsub("-", "$-$", formatC(x, format = "f", digits = 1))
row_tex <- function(r) {
  man <- r$method == "Manual count"
  sprintf("%s & %s & %s [%s, %s] & %s & %s \\\\", r$lab, r$method, f3(r$D), f3(r$D_lo), f3(r$D_hi),
          if (man) "--" else sprintf("%s [%s, %s]", f1s(r$bias), f1s(r$bias_lo), f1s(r$bias_hi)),
          if (man) "--" else f3(r$rmse_d))
}
tab_rows <- function(df, lab_fun) {
  out <- character(0); prev <- ""
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]; key <- lab_fun(r); r$lab <- if (key == prev) "" else key; prev <- key
    r$method <- as.character(r$method)
    if (i > 1 && r$lab != "") out <- c(out, "\\midrule")
    out <- c(out, row_tex(r))
  }
  out
}
reg_lab <- c(GB = "Georges Bank", MAB = "Mid-Atlantic Bight")
L <- c("\\begin{table}[htbp]", "\\centering",
       "\\caption{Scallop density on the Holdout test estimated by each workflow, by region and by region and survey year.}",
       "\\label{tab:density}", "\\resizebox{\\textwidth}{!}{%", "\\begin{tabular}{llccc}", "\\toprule",
       "\\textbf{Stratum} & \\textbf{Workflow} & \\textbf{Density (m$^{-2}$)} & \\textbf{Bias (\\%)} & \\textbf{RMSE (m$^{-2}$)} \\\\", "\\midrule",
       tab_rows(by_region, function(r) reg_lab[[r$region]]), "\\midrule",
       tab_rows(by_year, function(r) paste(reg_lab[[r$region]], r$year)),
       "\\bottomrule", "\\end{tabular}%", "}", "\\vspace{2mm}", "\\begin{flushleft}", "\\footnotesize",
       paste0("Density: total scallops divided by total imaged area (field of view summed over images). ",
              "Bias: difference from the manual density; it equals the total-count bias because both use the same area. ",
              "RMSE: per-image density error. Brackets: 95\\% transect-block bootstrap intervals (", B, " resamples). ",
              "Survey-year rows: years with at least 50 Holdout test scallops in the region; 97 Mid-Atlantic Bight images with a mislabelled cruise year are included in the regional rows only."),
       "\\end{flushleft}", "\\end{table}")
writeLines(L, "tables/tableS6_density.tex")
cat("\nWrote", OUT, "and tables/tableS6_density.tex\n")
