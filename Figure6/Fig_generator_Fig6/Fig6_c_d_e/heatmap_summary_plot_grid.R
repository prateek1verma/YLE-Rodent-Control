############################################################
# plot_separate_heatmaps_from_summary_csv.R
#
# Reads a summary CSV with columns including:
#   rel, x_shred, rel_pct, x_shred_pct, rel_lab, x_shred_lab,
#   A_prob_elim, tte_to_plot, D_eq_final_norm_mean
#
# Saves three separate heatmaps:
#   1) Probability of elimination
#   2) Median/selected time to elimination
#   3) Mean total population after 50 years, normalised to K
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
  library(magick)
  library(showtext)
})

############################################################
# User settings
############################################################

font_add("Arial", "/Users/prateekverma/Documents/Fonts/arial.ttf")   # mac
# or path to Arial.ttf on your system

showtext_auto()

csv_file <- "heatmap_summary_plot_grid.csv"
out_dir  <- "heatmap_outputs"

# Axis variables
x_var <- "rel_lab"
y_var <- "x_shred_lab"

# Optional filters; use NULL to include all values
rel_range_pct     <- c(0.2, 4.0)   # e.g. c(0.2, 2.0), or NULL
x_shred_range     <- c(0, 1.0)     # e.g. c(0, 1), or NULL

# Tile size controls; tiles do NOT need to be square
tile_width_cm  <- 1.7
tile_height_cm <- 1.0

# Figure padding around heatmap body
extra_width_cm  <- 4.5
extra_height_cm <- 4.5

# Text sizes
base_size       <- 18
axis_text_size  <- 18
axis_title_size <- 22
title_size      <- 22
legend_text_size <- 16

# Axis tick density
x_tick_by <- 0.4    # show x labels every 0.2 release percentage units
y_tick_by <- 0.1    # show y labels every 0.1 x-shredder units

# Tile border/gap
border_colour <- "white"
border_width  <- 0.25

# Colour palettes
colour_prob <- "#006d2c"
colour_tte  <- "#D55E00"
colour_pop  <- "#0072B2"

# TTE display
mask_tte_when_poe_below <- TRUE
poe_threshold_for_tte   <- 0.5
highlight_tte_range     <- c(0, 5)  # NULL to disable

# Output options
save_png <- FALSE
save_pdf <- TRUE
dpi      <- 600

############################################################
# Helper functions
############################################################

lighten_hex <- function(col, amount = 0.90) {
  rgb(t(col2rgb(col) * (1 - amount) + 255 * amount) / 255)
}

num_from_label <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

get_limits_breaks <- function(x, fallback = c(0, 1), n = 5) {
  x <- x[is.finite(x)]
  lims <- if (length(x) > 0) range(x, na.rm = TRUE) else fallback

  if (!is.finite(diff(lims)) || diff(lims) == 0) {
    eps <- max(abs(lims[1]) * 1e-6, 1e-6)
    lims <- c(lims[1] - eps, lims[1] + eps)
  }

  brks <- scales::breaks_pretty(n = n)(lims)
  brks <- brks[brks >= lims[1] & brks <= lims[2]]

  list(lims = lims, breaks = brks)
}

make_axis <- function(levels, tick_by = NULL) {
  vals <- num_from_label(levels)

  if (is.null(tick_by)) {
    keep <- rep(TRUE, length(levels))
  } else {
    keep <- abs((vals / tick_by) - round(vals / tick_by)) < 1e-8
  }

  list(
    breaks = seq_along(levels)[keep],
    labels = levels[keep]
  )
}

build_tte_perimeter <- function(df, lo, hi) {
  inside <- df %>%
    filter(is.finite(tte_to_plot),
           tte_to_plot >= lo,
           tte_to_plot <= hi) %>%
    select(x_idx_plot, y_idx_plot)

  if (nrow(inside) == 0) return(NULL)

  left <- inside %>%
    anti_join(inside %>% transmute(x_idx_plot = x_idx_plot + 1, y_idx_plot),
              by = c("x_idx_plot", "y_idx_plot")) %>%
    transmute(x0 = x_idx_plot - 0.5, y0 = y_idx_plot - 0.5,
              x1 = x_idx_plot - 0.5, y1 = y_idx_plot + 0.5)

  right <- inside %>%
    anti_join(inside %>% transmute(x_idx_plot = x_idx_plot - 1, y_idx_plot),
              by = c("x_idx_plot", "y_idx_plot")) %>%
    transmute(x0 = x_idx_plot + 0.5, y0 = y_idx_plot - 0.5,
              x1 = x_idx_plot + 0.5, y1 = y_idx_plot + 0.5)

  bottom <- inside %>%
    anti_join(inside %>% transmute(x_idx_plot, y_idx_plot = y_idx_plot + 1),
              by = c("x_idx_plot", "y_idx_plot")) %>%
    transmute(x0 = x_idx_plot - 0.5, y0 = y_idx_plot - 0.5,
              x1 = x_idx_plot + 0.5, y1 = y_idx_plot - 0.5)

  top <- inside %>%
    anti_join(inside %>% transmute(x_idx_plot, y_idx_plot = y_idx_plot - 1),
              by = c("x_idx_plot", "y_idx_plot")) %>%
    transmute(x0 = x_idx_plot - 0.5, y0 = y_idx_plot + 0.5,
              x1 = x_idx_plot + 0.5, y1 = y_idx_plot + 0.5)

  bind_rows(left, right, bottom, top) %>% distinct()
}

base_theme <- function() {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      text = element_text(family = "Arial"),
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.0),
      legend.position = "right",
      legend.title = element_text(size = legend_text_size),
      legend.text = element_text(size = legend_text_size),
      plot.margin = margin(0, 0, 0, 0, unit = "mm")
    )
}

plot_heatmap <- function(df, fill_col, title, legend_title, high_colour,
                         limits = NULL, breaks = NULL,
                         add_tte_perimeter = FALSE) {

  if (is.null(limits) || is.null(breaks)) {
    scale_info <- get_limits_breaks(df[[fill_col]])
    limits <- scale_info$lims
    breaks <- scale_info$breaks
  }

  p <- ggplot(df, aes(x = x_idx_plot, y = y_idx_plot, fill = .data[[fill_col]])) +
    geom_tile(width = 1, height = 1,
              colour = border_colour, linewidth = border_width) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = x_axis$breaks,
      labels = x_axis$labels
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = y_axis$breaks,
      labels = y_axis$labels
    ) +
    scale_fill_gradient(
      low = lighten_hex(high_colour, 0.90),
      high = high_colour,
      name = legend_title,
      limits = limits,
      breaks = breaks,
      labels = label_number(accuracy = 0.1),
      oob = squish,
      na.value = "grey90"
    ) +
    coord_fixed(ratio = tile_width_cm / tile_height_cm,
                expand = FALSE, clip = "on") +
    labs(
      title = title,
      x = "% of released transgenic males\n at K per month",
      y = "X-shredder efficiency"
    ) +
    base_theme() +
    # guides(
    #   fill = guide_colorbar(
    #     direction = "horizontal",
    #     title.position = "top",
    #     barheight = unit(0.45, "cm"),
    #     barwidth = unit(9.5, "cm"),
    #     ticks.colour = "black"
    #   )
        guides(
          fill = guide_colorbar(
            direction = "vertical",
            title.position = "top",
            barheight = unit(11, "cm"),
            barwidth  = unit(0.5, "cm")
          )
    )
  
  

  if (add_tte_perimeter && !is.null(highlight_tte_range)) {
    seg <- build_tte_perimeter(df, highlight_tte_range[1], highlight_tte_range[2])

    if (!is.null(seg) && nrow(seg) > 0) {
      p <- p +
        geom_segment(
          data = seg,
          aes(x = x0, y = y0, xend = x1, yend = y1),
          inherit.aes = FALSE,
          linewidth = 1.1,
          lineend = "square",
          colour = high_colour,
          linetype = "dotted"
        )
    }
  }

  p
}

save_heatmap <- function(plot, filename_stub, nx, ny) {
  width_cm  <- nx * tile_width_cm  + extra_width_cm
  height_cm <- ny * tile_height_cm + extra_height_cm

  if (save_png) {
    ggsave(
      filename = file.path(out_dir, paste0(filename_stub, ".png")),
      plot = plot,
      width = width_cm,
      height = height_cm,
      units = "cm",
      dpi = dpi,
      bg = "white"
    )
  }

  if (save_pdf) {
    ggsave(
      filename = file.path(out_dir, paste0(filename_stub, ".pdf")),
      plot = plot,
      width = width_cm,
      height = height_cm,
      units = "cm",
      bg = "white",
      device = cairo_pdf
    )
  }
}

############################################################
# Read and prepare data
############################################################

if (!file.exists(csv_file)) stop("CSV file not found: ", csv_file)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grid <- read_csv(csv_file, show_col_types = FALSE)

required_cols <- c(
  "rel_pct", "x_shred", "rel_lab", "x_shred_lab",
  "A_prob_elim", "tte_to_plot", "D_eq_final_norm_mean"
)

missing_cols <- setdiff(required_cols, names(grid))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Optional filtering
grid <- grid %>%
  filter(
    if (is.null(rel_range_pct)) TRUE else rel_pct >= rel_range_pct[1] & rel_pct <= rel_range_pct[2],
    if (is.null(x_shred_range)) TRUE else x_shred >= x_shred_range[1] & x_shred <= x_shred_range[2]
  )

if (nrow(grid) == 0) stop("No rows left after filtering.")

# Rebuild contiguous plotting indices after filtering
x_levels <- grid %>%
  distinct(rel_lab) %>%
  mutate(x_num = num_from_label(rel_lab)) %>%
  arrange(x_num) %>%
  pull(rel_lab)

y_levels <- grid %>%
  distinct(x_shred_lab) %>%
  mutate(y_num = num_from_label(x_shred_lab)) %>%
  arrange(y_num) %>%
  pull(x_shred_lab)

grid <- grid %>%
  mutate(
    x_idx_plot = as.integer(factor(rel_lab, levels = x_levels)),
    y_idx_plot = as.integer(factor(x_shred_lab, levels = y_levels)),
    tte_to_plot = ifelse(
      mask_tte_when_poe_below & A_prob_elim < poe_threshold_for_tte,
      NA_real_,
      tte_to_plot
    )
  )

# x_axis <- make_axis(x_levels, tick_by = x_tick_by)
manual_x_labels <- c(0.2, 1.0, 2.0, 3.0, 4.0) #
idx <- match(manual_x_labels, as.numeric(x_levels))
x_axis <- list( breaks = idx[!is.na(idx)], labels = manual_x_labels[!is.na(idx)] )

y_axis <- make_axis(y_levels, tick_by = y_tick_by)


nx <- length(x_levels)
ny <- length(y_levels)

############################################################
# Create heatmaps
############################################################

prob_scale <- get_limits_breaks(grid$A_prob_elim, fallback = c(0, 1))
tte_scale  <- get_limits_breaks(grid$tte_to_plot, fallback = c(0, 1))
pop_scale  <- get_limits_breaks(grid$D_eq_final_norm_mean, fallback = c(0, 1))

p_prob <- plot_heatmap(
  df = grid,
  fill_col = "A_prob_elim",
  title = "Probability of elimination (PoE)\nof female population",
  legend_title = " ",
  high_colour = colour_prob,
  limits = prob_scale$lims,
  breaks = prob_scale$breaks
)

p_tte <- plot_heatmap(
  df = grid,
  fill_col = "tte_to_plot",
  title = "Median time to elimination of\nfemale population (when PoE>0.5)",
  legend_title = " ",
  high_colour = colour_tte,
  limits = tte_scale$lims,
  breaks = tte_scale$breaks,
  add_tte_perimeter = TRUE
) +
  ggtext::geom_richtext(
    x = 0.5 + 0.80*nx, y = 0.5 + 0.9*10,
    label = "Median TTE<br>&le; 5 yrs", fill = NA, label.color = NA, size = 0
  )

p_pop <- plot_heatmap(
  df = grid,
  fill_col = "D_eq_final_norm_mean",
  title = "Mean of total population after\n50 yrs normalised to K",
  legend_title = " ",
  high_colour = colour_pop,
  limits = pop_scale$lims,
  breaks = pop_scale$breaks
)

############################################################
# Save heatmaps
############################################################

save_heatmap(p_prob, "heatmap_probability_elimination", nx, ny)
save_heatmap(p_tte,  "heatmap_median_time_to_elimination", nx, ny)
save_heatmap(p_pop,  "heatmap_mean_total_population_norm", nx, ny)

# img1 <- image_read("heatmap_outputs/heatmap_probability_elimination.png") |> image_trim(fuzz = 10)
# image_write(img1, "heatmap_outputs/heatmap_probability_elimination.png")
# 
# img2 <- image_read("heatmap_outputs/heatmap_median_time_to_elimination.png") |> image_trim(fuzz = 10)
# image_write(img2, "heatmap_outputs/heatmap_median_time_to_elimination.png")
# 
# img3 <- image_read("heatmap_outputs/heatmap_mean_total_population_norm.png") |> image_trim(fuzz = 10)
# image_write(img3, "heatmap_outputs/heatmap_mean_total_population_norm.png")

# img1 <- image_read("heatmap_outputs/heatmap_probability_elimination.pdf") |> image_trim(fuzz = 10)
# image_write(img1, "heatmap_outputs/heatmap_probability_elimination.pdf")
# 
# img2 <- image_read("heatmap_outputs/heatmap_median_time_to_elimination.pdf") |> image_trim(fuzz = 10)
# image_write(img2, "heatmap_outputs/heatmap_median_time_to_elimination.pdf")
# 
# img3 <- image_read("heatmap_outputs/heatmap_mean_total_population_norm.pdf") |> image_trim(fuzz = 10)
# image_write(img3, "heatmap_outputs/heatmap_mean_total_population_norm.pdf")



cat("\nSaved heatmaps to: ", normalizePath(out_dir), "\n", sep = "")
