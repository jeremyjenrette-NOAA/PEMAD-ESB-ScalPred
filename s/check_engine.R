# ======================================================================
# Quick check (about 5-10 min): why did the tensor-product formulas fail
# under bam in the 2026-10-05 run? Refits a few Georges Bank formulas on all
# of test_GAM_train with three engines and reports convergence and deviance
# explained. A usable fit has converged = TRUE and dev_expl > 0.
# Run from s/:  Rscript check_engine.R
# ======================================================================
library(dplyr); library(tidyr); library(mgcv)
source("./cvfunc.R")
models <- readRDS("../data/processed/eval_2226.rds")
forms <- list(
  M05_Boxsize = truedetect ~ s(conf, k = 5) + s(latitude, longitude, k = 7) + s(bottom_depth, k = 5) + s(altitude, k = 5) + s(boxsize, k = 5),
  M12_PB_Depth = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(bottom_depth, k = 5),
  M15_All_te   = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) + te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + te(conf, bottom_depth, k = 5)
)
engines <- list(gam = list(engine = "gam", discrete = FALSE),
                bam = list(engine = "bam", discrete = FALSE),
                bam_discrete = list(engine = "bam", discrete = TRUE))
res <- list()
for (m in c("YOLOv12", "Cascade R-CNN")) {
  d <- det_fit_rows(models[[m]]$det %>% filter(dataset == "test_GAM_train", region == "GB"))
  for (f in names(forms)) for (e in names(engines)) {
    t0 <- Sys.time()
    fit <- tryCatch(fit_det_gam(forms[[f]], d, engine = engines[[e]]$engine, discrete = engines[[e]]$discrete),
                    error = function(err) NULL)
    res[[length(res) + 1]] <- tibble(model = m, formula = f, engine = e,
      converged = if (is.null(fit)) NA else isTRUE(fit$converged),
      dev_expl = if (is.null(fit)) NA_real_ else round(dev_expl(fit), 4),
      AIC = if (is.null(fit)) NA_real_ else round(AIC(fit)),
      secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    print(tail(bind_rows(res), 1))
  }
}
print(bind_rows(res), n = Inf)
