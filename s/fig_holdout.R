# ======================================================================
# fig_holdout.R - Fig. 4 (predicted vs manual counts, the YOLOv12 calibration
# steps) and Fig. S3 (Holdout test totals by workflow), from the
# n3 Holdout table. Replaces sections 4-6 of gammod.R.
# Run from s/ after n3/n3_export_holdout.R:  Rscript fig_holdout.R
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(grid) })
ev <- readRDS("../data/processed/holdout_eval_2226_v2.rds")
FIG <- "../figures/ms_figures"

r2_11 <- function(y, x) 1 - sum((y - x)^2) / sum((y - mean(y))^2)
rmse  <- function(y, x) sqrt(mean((y - x)^2))

# ---------------- Fig. 4: rows = calibration steps (YOLOv12) ----------------
ROWS <- list(list(col = "yolo_f1",      title = expression("F"[1] ~ "threshold"),        colour = "#FB6A4A"),
             list(col = "pred_yolo",    title = expression("Stage 1 only:" ~ Sigma ~ p[det]), colour = "#CB181D"),
             list(col = "pred_2s_yolo", title = "Two-stage calibration",               colour = "#7570B3"))
zoom <- ev %>% group_by(region) %>% summarise(lim = max(5, quantile(n_annotations, 0.95)))
zl <- setNames(zoom$lim, zoom$region)

panel <- function(r, reg, zoomed, show_x, show_y, show_title) {
  d <- ev %>% filter(region == reg)
  lab <- sprintf("R²(1:1) = %.3f\nRMSE = %.2f", r2_11(d$n_annotations, d[[r$col]]), rmse(d$n_annotations, d[[r$col]]))
  p <- ggplot(d, aes(.data[[r$col]], n_annotations)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(colour = r$colour, size = 1.1, alpha = 0.5) +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.1, vjust = 1.3, size = 3.2) +
    theme_bw(base_size = 11) +
    labs(x = "Predicted count", y = "Manual count",
         title = if (show_title) (if (zoomed) "95% of images" else "All images") else NULL) +
    theme(plot.title = element_text(size = 10, hjust = 0.5), panel.grid.minor = element_blank())
  if (zoomed) p <- p + coord_cartesian(xlim = c(0, zl[[reg]]), ylim = c(0, zl[[reg]]))
  if (!show_x) p <- p + theme(axis.title.x = element_blank())
  if (!show_y) p <- p + theme(axis.title.y = element_blank())
  p
}
rows <- lapply(seq_along(ROWS), function(i) {
  r <- ROWS[[i]]; last <- i == length(ROWS); first <- i == 1
  row_lab <- wrap_elements(textGrob(r$title, rot = 90, gp = gpar(fontsize = 11)))
  row_lab | panel(r, "GB", FALSE, last, TRUE, first) | panel(r, "GB", TRUE, last, FALSE, first) |
    panel(r, "MAB", FALSE, last, FALSE, first) | panel(r, "MAB", TRUE, last, FALSE, first)
})
rows <- lapply(rows, function(x) x + plot_layout(widths = c(0.12, 1, 1, 1, 1)))
hdr <- wrap_elements(textGrob("")) |
  wrap_elements(textGrob("Georges Bank", gp = gpar(fontsize = 14, fontface = "bold"))) |
  wrap_elements(textGrob("Mid-Atlantic Bight", gp = gpar(fontsize = 14, fontface = "bold")))
hdr <- hdr + plot_layout(widths = c(0.12, 2, 2))
fig4 <- wrap_plots(c(list(hdr), rows), ncol = 1, heights = c(0.08, 1, 1, 1))
ggsave(file.path(FIG, "2226_holdout_1to1_yolo_steps.png"), fig4, width = 13, height = 10.5, dpi = 300, bg = "white")
ggsave(file.path(FIG, "2226_holdout_1to1_yolo_steps.pdf"), fig4, width = 13, height = 10.5,
       device = if (capabilities("cairo")) cairo_pdf else "pdf")

# ---------------- Fig. S3: totals by workflow and region ----------------
M <- c("YOLOv12 F1" = "yolo_f1", "YOLOv12 stage 1" = "pred_yolo", "Two-stage YOLOv12" = "pred_2s_yolo",
       "Cascade F1" = "cascade_f1", "Cascade stage 1" = "pred_cascade", "Two-stage YOLOv12 + Cascade" = "pred_2s_yc")
COL <- c("Manual count" = "grey40", "YOLOv12 F1" = "#FB6A4A", "YOLOv12 stage 1" = "#CB181D",
         "Two-stage YOLOv12" = "#7570B3", "Cascade F1" = "#6BAED6", "Cascade stage 1" = "#2171B5",
         "Two-stage YOLOv12 + Cascade" = "#3F007D")
tot <- ev %>% group_by(region) %>%
  summarise(`Manual count` = sum(n_annotations), across(all_of(unname(M)), sum), .groups = "drop") %>%
  rename(any_of(M)) %>%
  pivot_longer(-region, names_to = "Method", values_to = "Total") %>%
  group_by(region) %>% mutate(true = Total[Method == "Manual count"], err = 100 * (Total - true) / true) %>% ungroup() %>%
  mutate(Method = factor(Method, levels = names(COL)),
         region = factor(recode(region, GB = "Georges Bank", MAB = "Mid-Atlantic Bight"),
                         levels = c("Georges Bank", "Mid-Atlantic Bight")),
         lab = if_else(Method == "Manual count", format(round(Total), big.mark = ","),
                       sprintf("%s\n(%+.1f%%)", format(round(Total), big.mark = ","), err)))
figS3 <- ggplot(tot, aes(Method, Total, fill = Method)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_hline(aes(yintercept = true), linetype = "dashed") +
  geom_text(data = tot %>% group_by(region) %>% mutate(y_lab = 1.04 * max(Total)) %>% ungroup(),
            aes(y = y_lab, label = lab), vjust = 0, size = 3, lineheight = 0.9) +
  facet_wrap(~ region, scales = "free_y") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22)), labels = scales::comma) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Holdout test total (scallops)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
ggsave(file.path(FIG, "2226_final_abundance_summary.png"), figS3, width = 11, height = 5.5, dpi = 300, bg = "white")
OVL <- path.expand("~/Library/CloudStorage/Dropbox/Apps/Overleaf/NEFSC-scallop-detection/ms_figures")
if (dir.exists(OVL)) file.copy(file.path(FIG, c("2226_holdout_1to1_yolo_steps.pdf", "2226_final_abundance_summary.png")),
                               OVL, overwrite = TRUE)
cat("Wrote Fig. 4 (2226_holdout_1to1_yolo_steps) and Fig. S3 (2226_final_abundance_summary)\n")
