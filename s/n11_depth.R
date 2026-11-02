# ======================================================================
# n11_depth.R - R1.14: why does error rise between 70 and 100 m depth?
#
# Holdout test residuals (prediction - manual count per image) by bottom-
# depth bin, overall and within density classes. If the 70-100 m excess
# is a density effect (dense Georges Bank aggregations sit at 70-80 m),
# it should vanish within density classes.
#
# Run from s/:  Rscript n11_depth.R
# Needs: ../data/processed/holdout_eval_2226_v2.rds (n3/n3_export_holdout.R)
# Outputs: ../data/processed/n11_depth/*.csv, ../figures/ms_figures/2226_depth_density_error.pdf|png
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(purrr) })
OUT <- "../data/processed/n11_depth"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG <- "../figures/ms_figures"
METHODS <- c("Two-stage YOLOv12" = "pred_2s_yolo", "YOLOv12 F1" = "yolo_f1", "Cascade F1" = "cascade_f1",
             "YOLOv12 stage 1" = "pred_yolo", "Cascade stage 1" = "pred_cascade")
DEPTH_BREAKS <- c(0, 50, 60, 70, 80, 90, 100, Inf)
DEPTH_LABELS <- c("<50", "50-60", "60-70", "70-80", "80-90", "90-100", ">100")
pr <- function(x) print(tibble::as_tibble(x), n = Inf, width = Inf)

ev <- readRDS("../data/processed/holdout_eval_2226_v2.rds") %>%
  filter(dataset == "test_GAM_test") %>%
  mutate(region = as.character(region),
         bottom_depth = if_else(bottom_depth <= -90, NA_real_, bottom_depth),
         depth_bin = cut(bottom_depth, DEPTH_BREAKS, labels = DEPTH_LABELS, right = FALSE),
         density = cut(n_annotations, c(-Inf, 0, 2, 10, Inf), labels = c("0", "1-2", "3-10", ">10"))) %>%
  filter(!is.na(depth_bin))
stopifnot(all(METHODS %in% names(ev)))

# Block-clustered SE of a mean (transect blocks), as in error.R
cl_se <- function(e, cl) {
  k <- length(unique(cl)); if (k < 2) return(NA_real_)
  m <- mean(e); s <- tapply(e - m, cl, sum)
  sqrt(k / (k - 1) * sum(s^2)) / length(e)
}
long <- ev %>%
  pivot_longer(all_of(unname(METHODS)), names_to = "col", values_to = "pred") %>%
  mutate(method = factor(names(METHODS)[match(col, METHODS)], levels = names(METHODS)),
         err = pred - n_annotations)

summ <- function(d, ...) d %>% group_by(region, method, ...) %>%
  summarise(n_images = n(), scallops = sum(n_annotations), mean_count = mean(n_annotations),
            mean_err = mean(err), se = cl_se(err, transect_block_id), mae = mean(abs(err)),
            pct_bias = 100 * sum(err) / pmax(sum(n_annotations), 1), .groups = "drop")

by_depth   <- summ(long, depth_bin)
by_depth_d <- summ(long, depth_bin, density)
composition <- ev %>% group_by(region, depth_bin) %>%
  summarise(n_images = n(), scallops = sum(n_annotations), mean_count = mean(n_annotations),
            share_gt10 = mean(n_annotations > 10), .groups = "drop") %>%
  group_by(region) %>% mutate(share_of_region_scallops = scallops / sum(scallops)) %>% ungroup()

# Does depth still matter once density is accounted for? Within-density
# regression of error on depth bin (Two-stage YOLOv12, GB), block-clustered.
anova_like <- long %>% filter(method == "Two-stage YOLOv12") %>%
  group_by(region) %>%
  group_modify(function(d, k) {
    m0 <- lm(err ~ density, data = d); m1 <- lm(err ~ density + depth_bin, data = d)
    a <- anova(m0, m1)
    tibble(F_depth_given_density = a$F[2], p_depth_given_density = a$`Pr(>F)`[2],
           r2_density = summary(m0)$r.squared, r2_density_depth = summary(m1)$r.squared)
  }) %>% ungroup()

write.csv(composition, file.path(OUT, "depth_composition.csv"), row.names = FALSE)
write.csv(by_depth, file.path(OUT, "error_by_depth.csv"), row.names = FALSE)
write.csv(by_depth_d, file.path(OUT, "error_by_depth_density.csv"), row.names = FALSE)
write.csv(anova_like, file.path(OUT, "depth_given_density_test.csv"), row.names = FALSE)

cat("\nHoldout composition by depth:\n"); pr(composition %>% mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("\nMean error by depth (all images):\n")
pr(by_depth %>% filter(method %in% c("Two-stage YOLOv12", "YOLOv12 F1", "Cascade F1")) %>%
     mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("\nMean error by depth within density class (Two-stage YOLOv12):\n")
pr(by_depth_d %>% filter(method == "Two-stage YOLOv12") %>% mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("\nDoes depth add to density? (Two-stage YOLOv12 error ~ density + depth bin):\n"); pr(anova_like)

pd <- by_depth_d %>% filter(method %in% c("Two-stage YOLOv12", "YOLOv12 F1", "Cascade F1"), n_images >= 5)
p <- ggplot(pd, aes(depth_bin, mean_err, colour = method, group = method)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = mean_err - 1.96 * se, ymax = mean_err + 1.96 * se),
                  position = position_dodge(0.5), size = 0.25) +
  facet_grid(region ~ density, scales = "free_y", labeller = labeller(density = function(x) paste(x, "scallops"))) +
  scale_colour_manual(values = c("Two-stage YOLOv12" = "#2a78d6", "YOLOv12 F1" = "#1baf7a", "Cascade F1" = "#eb6834")) +
  labs(x = "Bottom depth (m)", y = "Mean error per image (predicted - manual)", colour = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(FIG, "2226_depth_density_error.pdf"), p, width = 8, height = 4.8)
ggsave(file.path(FIG, "2226_depth_density_error.png"), p, width = 8, height = 4.8, dpi = 300, bg = "white")
cat("\nDone. Tables in", OUT, "\n")
