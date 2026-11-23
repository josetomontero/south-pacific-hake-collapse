# Shared figure styling and export helpers (ggplot2).
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# Segment colours between biomass breaks. Okabe-Ito hues, with the yellow and sky blue replaced
# by darker steps (#B8860B, #0072B2) so every colour reaches 3:1 contrast on white and adjacent
# pairs stay separable under deutan/protan/tritan vision (checked with the dataviz palette
# validator: all checks pass). Black marks the collapsed period (data ink, not a hue).
SEG_COLS <- c("1968-1974" = "#D55E00", "1975-1983" = "#009E73", "1984-1992" = "#B8860B",
              "1993-2001" = "#0072B2", "2002-2020" = "#000000")
seg_label <- function(year) cut(year, c(-Inf, 1974, 1983, 1992, 2001, Inf), labels = names(SEG_COLS))

PRED_COL <- "#C0392B"   # predicted trajectories
OBS_COL  <- "black"

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
          axis.title = element_text(size = base_size), axis.text = element_text(size = base_size - 1, colour = "black"),
          plot.tag = element_text(size = base_size + 1, face = "bold"),
          legend.text = element_text(size = base_size - 1), legend.title = element_blank(),
          legend.key.size = unit(0.4, "cm"), legend.background = element_blank(),
          plot.margin = margin(4, 6, 2, 2))
}

# Axis style of Figure 3 (R/08_Fig3.R): thin 0.3 lines and 0.12 cm ticks. Added to Figure 2 and to the
# phase portrait so that the three figures share the same axes.
theme_axes_fig3 <- theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
                         axis.ticks.length = unit(0.12, "cm"))

# Journal figure sizes: max 18 cm wide x 22 cm tall. Figures are written as PNG at 600 dpi only.
save_fig <- function(plot, name, width_cm, height_cm, dir = file.path(RES, "figs")) {
  ggsave(file.path(dir, paste0(name, ".png")), plot, width = width_cm, height = height_cm, units = "cm", dpi = 600, bg = "white")
  cat("wrote", paste0(name, ".png"), "in", sub(paste0(PROJ, "/"), "", dir), "\n")
}

# Figure numbers are not printed inside the files (they go in the file names and captions)
figure_label <- function(plot, n) plot
