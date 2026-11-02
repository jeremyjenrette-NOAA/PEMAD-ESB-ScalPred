# ======================================================================
# fig2_design.R - Fig. 2: study regions, survey years, data splits and
# the density balance the stratification achieved.
#   (A) bathymetric overview with annotated image locations, coloured by survey year
#   (B) annotated images per region and year, stacked by data split
#       (scallop totals printed above each bar)
#   (C) share of images in each density class (0, 1-2, 3-10, >10 scallops)
#       for each split, by region - the classes used in Fig. S9
# Run from s/:  Rscript fig2_design.R
# Needs: ../data/raw/dataset_split_2226.csv; packages marmap, maps, patchwork
#        (marmap caches the NOAA bathymetry download in s/ after the first run)
# Writes ../figures/ms_figures/2226_study_design.pdf|png (copied to Overleaf if present)
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(maps); library(marmap); library(scales) })

SPLITS <- c(train = "Detector-training", test_GAM_train = "Calibration-fitting", test_GAM_test = "Holdout test")
SPLIT_COL <- c("Detector-training" = "#4D4D4D", "Calibration-fitting" = "#2C7FB8", "Holdout test" = "#D95F0E")
YEAR_COL <- c("2022" = "#E69F00", "2023" = "#56B4E9", "2024" = "#009E73", "2026" = "#CC79A7")
DENS <- c("0", "1-2", "3-10", ">10")
DENS_COL <- setNames(c("#F1EEF6", "#BDC9E1", "#74A9CF", "#0570B0"), DENS)

d <- read.csv("../data/raw/dataset_split_2226.csv") %>%
  filter(!is.na(longitude), !is.na(latitude), dataset %in% names(SPLITS)) %>%
  mutate(split = factor(SPLITS[dataset], levels = SPLITS),
         region = factor(ifelse(longitude >= -71, "Georges Bank", "Mid-Atlantic Bight"),
                         levels = c("Georges Bank", "Mid-Atlantic Bight")),
         year = factor(year),
         dens = factor(cut(n_annotations, c(-Inf, 0, 2, 10, Inf), labels = DENS), levels = DENS))
cat("Images:", nrow(d), "| scallops:", sum(d$n_annotations), "\n")
print(d %>% count(region, year, split) %>% pivot_wider(names_from = split, values_from = n, values_fill = 0))

th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold"))

# ---------------- (A) overview map ----------------
world <- map_data("world"); states <- map_data("state")
bathy <- fortify.bathy(getNOAA.bathy(lon1 = -78, lon2 = -64, lat1 = 36, lat2 = 44, resolution = 1, keep = TRUE))
box <- d %>% group_by(region) %>%
  summarise(xmin = min(longitude) - 0.2, xmax = max(longitude) + 0.2,
            ymin = min(latitude) - 0.2, ymax = max(latitude) + 0.2, .groups = "drop")
pA <- ggplot() +
  geom_raster(data = bathy %>% filter(z <= 0), aes(x, y, fill = z)) +
  scale_fill_gradientn(colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#EEF6FA"),
                       values = rescale(c(-4000, -2000, -500, -100, 0)), guide = "none") +
  geom_contour(data = bathy, aes(x, y, z = z), breaks = -100, colour = "#08519C", linewidth = 0.3, alpha = 0.7) +
  geom_polygon(data = world, aes(long, lat, group = group), fill = "grey88", colour = "grey55", linewidth = 0.2) +
  geom_polygon(data = states, aes(long, lat, group = group), fill = NA, colour = "grey60", linewidth = 0.15) +
  geom_point(data = d %>% arrange(year), aes(longitude, latitude, colour = year), size = 0.35, alpha = 0.6) +
  geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), colour = "black",
            fill = NA, linewidth = 0.35, linetype = 2) +
  geom_label(data = box, aes(x = xmin, y = ifelse(region == "Georges Bank", ymax, ymin), label = region,
                             vjust = ifelse(region == "Georges Bank", -0.25, 1.25)),
             hjust = 0, size = 2.8, fill = alpha("white", 0.8)) +
  scale_colour_manual(values = YEAR_COL, name = "Survey year") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_quickmap(xlim = c(-77, -65), ylim = c(36.5, 43), expand = FALSE) +
  labs(x = NULL, y = NULL, tag = "A") + th +
  theme(legend.position = "inside", legend.position.inside = c(0.985, 0.03), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(0.35, "cm"), axis.text = element_text(size = 7))

# ---------------- (B) images by region-year and split ----------------
cnt <- d %>% count(region, year, split, name = "images")
tot <- d %>% group_by(region, year) %>%
  summarise(images = n(), scallops = sum(n_annotations), .groups = "drop")
pB <- ggplot(cnt, aes(year, images, fill = split)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(data = tot, aes(year, images, label = paste0(comma(scallops), "\nscallops")),
            inherit.aes = FALSE, vjust = -0.15, size = 2.3, lineheight = 0.85, colour = "grey25") +
  facet_grid(. ~ region, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = SPLIT_COL, name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Survey year", y = "Annotated images", tag = "B") + th +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

# ---------------- (C) density-class balance across splits ----------------
comp <- d %>% count(region, split, dens) %>% group_by(region, split) %>%
  mutate(share = n / sum(n)) %>% ungroup()
pC <- ggplot(comp, aes(split, share, fill = dens)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(share >= 0.05, percent(share, accuracy = 1), "")),
            position = position_stack(vjust = 0.5), size = 2.2,
            colour = ifelse(comp$dens %in% c("3-10", ">10"), "white", "grey20")) +
  facet_grid(. ~ region) +
  scale_fill_manual(values = DENS_COL, name = "Scallops per image") +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  scale_x_discrete(labels = function(x) sub("-", "-\n", sub(" ", "\n", x))) +
  labs(x = NULL, y = "Share of images", tag = "C") + th +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

fig <- pA / (pB | pC) + plot_layout(heights = c(1.25, 1))
ggsave("../figures/ms_figures/2226_study_design.pdf", fig, width = 8.2, height = 8.2)
ggsave("../figures/ms_figures/2226_study_design.png", fig, width = 8.2, height = 8.2, dpi = 300, bg = "white")
ovl <- path.expand("~/Library/CloudStorage/Dropbox/Apps/Overleaf/NEFSC-scallop-detection/ms_figures")
if (dir.exists(ovl)) file.copy("../figures/ms_figures/2226_study_design.pdf", ovl, overwrite = TRUE)
cat("Wrote 2226_study_design.pdf/png\n")
