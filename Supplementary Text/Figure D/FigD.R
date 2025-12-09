############################################################
# heatmaps_from_gridcsv_multi.R — fast plotting from saved grids
# Stacks (A,B,C) for multiple CSVs; horizontal legends under each panel
############################################################

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
  library(magick)
  library(knitr)
  library(RColorBrewer)
  library(patchwork)   # <— require patchwork for layout/labels
})

# ------------------------------- #
# User settings                   #
# ------------------------------- #

# (1) EITHER: provide a vector of root dirs (each containing heatmap_summary_plot_grid.csv)
roots <- c(
  "mgdriveYLE_sweep_5yr_release_target_fertility"
)

# (2) OR: point directly to CSV paths (leave roots empty if you use csv_files)
csv_files <- NULL
# csv_files <- c("run1/heatmap_summary_plot_grid.csv", "run2/heatmap_summary_plot_grid.csv")

years  <- 10
tte_highlight_range_years <- c(0, 5)
fpt <- 20
label_fsize <- 32

gap_col <- "white"
gap_mm  <- 0.6

tile_colorA <- "#006d2c"
tile_colorB <- "#D55E00"
tile_colorC <- "#0072B2"

rescale_to_filtered <- TRUE

# windowing (per CSV, but we’ll compute shared B-scale from the union after filtering)
rel_range <- c(1, 10)
fy_range  <- c(0, 90)

# ------------------------------- #
# Helpers (unchanged from yours)  #
# ------------------------------- #

lighten_hex <- function(col, amt = 0.88) rgb(t(col2rgb(col) * (1 - amt) + 255 * amt) / 255)

compute_tile_size_cm <- function(nx, ny, min_cm = 0.9, max_cm = 2.8) {
  base <- 2.0; adj <- base * (12 / max(12, max(nx, ny)))
  pmin(max_cm, pmax(min_cm, adj))
}

get_limits_breaks <- function(vals, n = 5, fallback = c(0, 1), dec = 1) {
  v <- vals[is.finite(vals)]
  lims <- if (length(v)) range(v, na.rm = TRUE) else fallback
  if (!is.finite(diff(lims)) || diff(lims) == 0) {
    eps <- max(1e-6, abs(lims[1]) * 1e-6, 1e-6)
    lims <- c(lims[1] - eps, lims[2] + eps)
  }
  raw_breaks <- scales::breaks_pretty(n)(lims)
  brks <- sort(unique(c(lims[1], raw_breaks, lims[2])))
  brks <- brks[brks >= lims[1] & brks <= lims[2]]
  lower <- lims[1]; upper <- lims[2]
  above_lower <- suppressWarnings(min(brks[brks > lower]))
  if (is.finite(above_lower) && round(above_lower, dec) == round(lower, dec)) {
    brks <- brks[brks != above_lower]
  }
  below_upper <- suppressWarnings(max(brks[brks < upper]))
  if (is.finite(below_upper) && round(below_upper, dec) == round(upper, dec)) {
    brks <- brks[brks != below_upper]
  }
  list(lims = lims, breaks = brks)
}

# ---- your robust range filter (unchanged) ----
filter_grid_by_ranges <- function(grid,
                                  rel_range = NULL,
                                  fy_range  = NULL,
                                  tol = 0) {
  g <- grid; stopifnot(is.data.frame(g))
  pick_axis_col <- function(df, pref = c("rel_pct","rel","rel_lab"), axis = "rel") {
    for (nm in pref) if (nm %in% names(df)) return(nm)
    stop(sprintf("No usable %s axis column found (looked for: %s).",
                 axis, paste(pref, collapse = ", ")))
  }
  rel_col <- pick_axis_col(g, c("rel_pct","rel","rel_lab"), "rel")
  fy_col  <- pick_axis_col(g, c("fy_cost_pct","fy","fy_cost_lab"), "fy")
  make_normalizer <- function(vals) {
    v <- suppressWarnings(as.numeric(vals)); v <- v[is.finite(v)]
    if (!length(v)) return(function(x) x)
    mx <- max(v, na.rm = TRUE)
    if (mx > 1.5) function(x) x else function(x) x / 100
  }
  rel_norm <- if (grepl("_pct$", rel_col)) make_normalizer(g[[rel_col]]) else function(x) x
  fy_norm  <- if (grepl("_pct$", fy_col )) make_normalizer(g[[fy_col ]]) else function(x) x
  if (!is.null(rel_range)) {
    rr <- sort(rel_range); rr_use <- if (grepl("_pct$", rel_col)) rel_norm(rr) else rr
    if (tol > 0) {
      uniq_vals <- sort(unique(suppressWarnings(as.numeric(g[[rel_col]]))))
      snap <- function(b) { d <- abs(uniq_vals - b); if (any(d <= tol)) uniq_vals[which.min(d)] else b }
      rr_use <- c(snap(rr_use[1]), snap(rr_use[2]))
    }
    g <- dplyr::filter(g, .data[[rel_col]] >= rr_use[1], .data[[rel_col]] <= rr_use[2])
  }
  if (!is.null(fy_range)) {
    fr <- sort(fy_range); fr_use <- if (grepl("_pct$", fy_col)) fy_norm(fr) else fr
    if (tol > 0) {
      uniq_vals <- sort(unique(suppressWarnings(as.numeric(g[[fy_col]]))))
      snap <- function(b) { d <- abs(uniq_vals - b); if (any(d <= tol)) uniq_vals[which.min(d)] else b }
      fr_use <- c(snap(fr_use[1]), snap(fr_use[2]))
    }
    g <- dplyr::filter(g, .data[[fy_col]] >= fr_use[1], .data[[fy_col]] <= fr_use[2])
  }
  if (nrow(g) == 0) stop("Filter removed all rows.", call. = FALSE)
  g
}

auto_crop_png <- function(infile, outfile = infile, fuzz = 12) {
  img <- magick::image_read(infile) |> magick::image_trim(fuzz = fuzz)
  magick::image_write(img, path = outfile); invisible(outfile)
}

# ------------------------------- #
# Load all CSVs (NEW)             #
# ------------------------------- #

if (is.null(csv_files) || length(csv_files) == 0) {
  csv_files <- file.path(roots, "heatmap_summary_plot_grid.csv")
}
stopifnot(length(csv_files) >= 1, all(file.exists(csv_files)))

# Read, filter per-file, and also collect for global B-scale
all_items <- purrr::map(csv_files, function(f) {
  g <- readr::read_csv(f, show_col_types = FALSE)
  need_cols <- c("rel_lab","fy_cost_lab","x_idx","y_idx",
                 "A_prob_elim","tte_to_plot","D_eq_final_norm_mean","tte_stat_choice")
  stopifnot(all(need_cols %in% names(g)), nrow(g) > 0)
  g2 <- filter_grid_by_ranges(g, rel_range, fy_range)
  # rebuild contiguous indices per file
  rel_levels     <- g2 |> distinct(x_idx, rel_lab)     |> arrange(x_idx)     |> pull(rel_lab)
  fy_cost_levels <- g2 |> distinct(y_idx, fy_cost_lab) |> arrange(y_idx)     |> pull(fy_cost_lab)
  g2 <- g2 |>
    mutate(
      x_idx_f = as.integer(factor(rel_lab,     levels = rel_levels)),
      y_idx_f = as.integer(factor(fy_cost_lab, levels = fy_cost_levels))
    )
  list(file = f, grid = g2, rel_levels = rel_levels, fy_levels = fy_cost_levels)
})

# ---- compute global limits (A, B, C). Per request, force a single COMMON scale for B. ----
src_union <- dplyr::bind_rows(lapply(all_items, `[[`, "grid"))
sb_elim_global <- get_limits_breaks(src_union$A_prob_elim, n = 5, fallback = c(0,1))
sb_tte_global  <- get_limits_breaks(src_union$tte_to_plot, n = 5,
                                    fallback = range(src_union$tte_to_plot, na.rm = TRUE), dec = 0.1)
sb_tot_global  <- get_limits_breaks(src_union$D_eq_final_norm_mean, n = 5, fallback = c(0,1))

# ------------------------------- #
# Common theme & sizing           #
# ------------------------------- #

# size from the FIRST item (each row uses its own axes but similar sizing)
nx <- length(all_items[[1]]$rel_levels); ny <- length(all_items[[1]]$fy_levels)
tile_size_cm <- compute_tile_size_cm(nx, ny)
legend_h_cm  <- 0.5    # horizontal bar height (cm)
legend_w_cm  <- max(4, nx * 1)  # horizontal bar width (cm)
panel_pad_cm <- 2.0
title_axis_pad <- 4.0

w_tiles_cm   <- nx * tile_size_cm; h_tiles_cm <- ny * tile_size_cm
w_single     <- w_tiles_cm + panel_pad_cm
h_single     <- h_tiles_cm + title_axis_pad + 1.2  # + some room for bottom legend
fixed_ratio  <- nx / ny

axis_text_x <- element_text(size = fpt)
axis_text_y <- element_text(size = fpt)

base_theme <- theme_minimal(base_size = fpt) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = fpt),
        axis.text.x = axis_text_x,
        axis.text.y = axis_text_y,
        plot.title  = element_text(face = "plain", size = 20, hjust = 0),
        plot.margin = margin(0, 0, 0, 0, unit = "mm"),
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.5, hjust = 0.5, margin = margin(t = 6), size = 18)
  )

# ------------------------------- #
# Optional perimeter for TTE band #
# ------------------------------- #
build_perimeter <- function(grid_plot, lo, hi) {
  grid_df <- grid_plot %>%
    mutate(mask = ifelse(is.finite(tte_to_plot) & tte_to_plot >= lo & tte_to_plot <= hi, 1, 0)) %>%
    select(x_idx_f, y_idx_f, mask)
  inside_df <- grid_df %>% filter(mask == 1) %>% select(x_idx_f, y_idx_f)
  if (nrow(inside_df) == 0) return(NULL)
  left_edges <- inside_df %>%
    anti_join(inside_df %>% transmute(x_idx_f = x_idx_f + 1, y_idx_f), by = c("x_idx_f","y_idx_f")) %>%
    transmute(x0 = x_idx_f - 0.5, y0 = y_idx_f - 0.5, x1 = x_idx_f - 0.5, y1 = y_idx_f + 0.5)
  right_edges <- inside_df %>%
    anti_join(inside_df %>% transmute(x_idx_f = x_idx_f - 1, y_idx_f), by = c("x_idx_f","y_idx_f")) %>%
    transmute(x0 = x_idx_f + 0.5, y0 = y_idx_f - 0.5, x1 = x_idx_f + 0.5, y1 = y_idx_f + 0.5)
  bottom_edges <- inside_df %>%
    anti_join(inside_df %>% transmute(x_idx_f, y_idx_f = y_idx_f + 1), by = c("x_idx_f","y_idx_f")) %>%
    transmute(x0 = x_idx_f - 0.5, y0 = y_idx_f - 0.5, x1 = x_idx_f + 0.5, y1 = y_idx_f - 0.5)
  top_edges <- inside_df %>%
    anti_join(inside_df %>% transmute(x_idx_f, y_idx_f = y_idx_f - 1), by = c("x_idx_f","y_idx_f")) %>%
    transmute(x0 = x_idx_f - 0.5, y0 = y_idx_f + 0.5, x1 = x_idx_f + 0.5, y1 = y_idx_f + 0.5)
  bind_rows(left_edges, right_edges, bottom_edges, top_edges) %>% distinct()
}

# ------------------------------- #
# Panel builder for one CSV (NEW) #
# ------------------------------- #

make_panels_for_grid <- function(item) {
  grid_plot <- item$grid
  rel_levels <- item$rel_levels
  fy_levels  <- item$fy_levels
  
  # x tick every other label (same as your code)
  rel_levels_lab <- as.character(rel_levels)
  even_indices <- rel_levels %% 2 == 0
  new_breaks <- seq_along(rel_levels)[even_indices]
  new_labels <- rel_levels_lab[even_indices]
  
  perimeter_segments <- NULL
  if (!is.null(tte_highlight_range_years)) {
    perimeter_segments <- build_perimeter(grid_plot, min(tte_highlight_range_years), max(tte_highlight_range_years))
  }
  
  # Panel A
  pA <- ggplot(grid_plot, aes(x = x_idx, y = y_idx, fill = A_prob_elim)) +
    geom_tile(width = 1, height = 1, color = gap_col, linewidth = gap_mm) +
    scale_x_continuous(expand = c(0, 0), breaks = new_breaks, labels = new_labels) +
    scale_y_continuous(expand = c(0, 0), breaks = seq_along(fy_levels), labels = fy_levels) +
    scale_fill_gradient(
      low = lighten_hex(tile_colorA, 0.90), high = tile_colorA,
      name = "",
      limits = if (rescale_to_filtered) get_limits_breaks(grid_plot$A_prob_elim)$lims else sb_elim_global$lims,
      breaks = if (rescale_to_filtered) get_limits_breaks(grid_plot$A_prob_elim)$breaks else sb_elim_global$breaks,
      labels = label_number(accuracy = 0.1),
      oob = squish, na.value = "grey90"
    ) +
    labs(
      title = "Probability of elimination (PoE)\nof female population",
      x = "% of released YLE males\nat K (per month)",
      y = "% of life-span reduction in YLE\nmales relative to WT"
    ) +
    coord_fixed(ratio = fixed_ratio, expand = FALSE, clip = "on") +
    base_theme +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top",
                                 barheight = unit(legend_h_cm, "cm"),
                                 barwidth  = unit(legend_w_cm, "cm"),
                                 ticks.colour = "black"))
  
  # Panel B — COMMON scale across ALL CSVs
  pB <- ggplot(grid_plot, aes(x = x_idx, y = y_idx, fill = tte_to_plot)) +
    geom_tile(width = 1, height = 1, color = gap_col, linewidth = gap_mm) +
    scale_x_continuous(expand = c(0, 0), breaks = new_breaks, labels = new_labels) +
    scale_y_continuous(expand = c(0, 0), breaks = seq_along(fy_levels), labels = fy_levels) +
    scale_fill_gradient(
      low = lighten_hex(tile_colorB, 0.90), high = tile_colorB,
      name = "(in yrs)",
      limits = sb_tte_global$lims,
      breaks = sb_tte_global$breaks,
      labels = label_number(accuracy = 0.1),
      oob    = squish,
      na.value = "grey90"
    ) +
    labs(
      title = "Median time to elimination of\nfemale population (when PoE>0.5)",
      x = "% of released YLE males\nat K (per month)",
      y = NULL
    ) +
    coord_fixed(ratio = fixed_ratio, expand = FALSE, clip = "on") +
    base_theme +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top",
                                 barheight = unit(legend_h_cm, "cm"),
                                 barwidth  = unit(legend_w_cm, "cm"),
                                 ticks.colour = "black")) +
    ggtext::geom_richtext(
      x = 0.5 + 0.80*nx, y = 0.5 + 0.20*ny,
      label = "Median TTE<br>&le; 5 yrs", fill = NA, label.color = NA, size = 7
    )
  
  if (!is.null(tte_highlight_range_years) && !is.null(perimeter_segments) && nrow(perimeter_segments) > 0) {
    pB <- pB + geom_segment(data = perimeter_segments,
                            aes(x = x0, y = y0, xend = x1, yend = y1),
                            inherit.aes = FALSE,
                            linewidth   = 1.5,
                            lineend     = "square",
                            color       = "#D55E00",
                            linetype    = "dotted")
  }
  
  # Panel C
  pC <- ggplot(grid_plot, aes(x = x_idx, y = y_idx, fill = D_eq_final_norm_mean)) +
    geom_tile(width = 1, height = 1, color = gap_col, linewidth = gap_mm) +
    scale_x_continuous(expand = c(0, 0), breaks = new_breaks, labels = new_labels) +
    scale_y_continuous(expand = c(0, 0), breaks = seq_along(fy_levels), labels = fy_levels) +
    scale_fill_gradient(
      low = lighten_hex(tile_colorC, 0.90), high = tile_colorC,
      name = "",
      limits = if (rescale_to_filtered) get_limits_breaks(grid_plot$D_eq_final_norm_mean)$lims else sb_tot_global$lims,
      breaks = if (rescale_to_filtered) get_limits_breaks(grid_plot$D_eq_final_norm_mean)$breaks else sb_tot_global$breaks,
      labels = label_number(accuracy = 0.1),
      oob    = squish,
      na.value = "grey90"
    ) +
    labs(
      title = "Mean of total population after\n50 yrs normalised to K",
      x = "% of released YLE males\nat K (per month)",
      y = NULL
    ) +
    coord_fixed(ratio = fixed_ratio, expand = FALSE, clip = "on") +
    base_theme +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top",
                                 barheight = unit(legend_h_cm, "cm"),
                                 barwidth  = unit(legend_w_cm, "cm"),
                                 ticks.colour = "black"))
  
  list(pA = pA, pB = pB, pC = pC)
}

# ------------------------------- #
# Build all rows (NEW)            #
# ------------------------------- #

# ------------------------------- #
# Build all rows (already made)   #
# ------------------------------- #

rows_list <- lapply(all_items, make_panels_for_grid)  # as before

# ------------------------------- #
# Column-wise assembly (NEW)      #
# - collect legends per column -> 3 bars total
# - only top row shows each column title
# - only bottom row shows x-axis label
# ------------------------------- #

n_rows <- length(rows_list)

# grab column lists
colA <- lapply(rows_list, function(x) x$pA)
colB <- lapply(rows_list, function(x) x$pB)
colC <- lapply(rows_list, function(x) x$pC)

# helper to blank titles except top, and x label except bottom
apply_columnwise_tidy <- function(col_list, title_for_top, xlab_for_bottom) {
  out <- col_list
  for (i in seq_along(out)) {
    # remove titles except at top
    if (i != 1) {
      out[[i]] <- out[[i]] + labs(title = NULL)
    } else {
      out[[i]] <- out[[i]] + labs(title = title_for_top)
    }
    # remove x-axis title except at bottom
    if (i != length(out)) {
      out[[i]] <- out[[i]] + labs(x = NULL) +
        theme(axis.title.x = element_blank())
    } else {
      out[[i]] <- out[[i]] + labs(x = xlab_for_bottom)
    }
  }
  out
}

# keep your original titles/x labels text, but only where needed
colA <- apply_columnwise_tidy(
  colA,
  title_for_top   = "Probability of elimination (PoE)\nof female population",
  xlab_for_bottom = "% of released YLE males\nat K (per month)"
)
colB <- apply_columnwise_tidy(
  colB,
  title_for_top   = expression(atop("Median time to elimination (TTE) of",
                                    "female population (when PoE" >= "0.5)")),
  xlab_for_bottom = "% of released YLE males\nat K (per month)"
)
colC <- apply_columnwise_tidy(
  colC,
  title_for_top   = "Mean of total population after\n50 yrs normalised to K",
  xlab_for_bottom = "% of released YLE males\nat K (per month)"
)

# collect each column into a stacked patchwork, with ONE legend at the bottom
colA_combined <- wrap_plots(colA, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")
colB_combined <- wrap_plots(colB, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")
colC_combined <- wrap_plots(colC, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

# --- squeeze horizontal gaps between columns hard ---

# 1) strip inner plot margins everywhere (helps a lot)
strip_margins <- function(p) p & theme(
  plot.margin       = margin(0, 0, 0, 0, unit = "mm"),
  legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
  legend.margin     = margin(0, 0, 0, 0, unit = "mm")
)

colA_combined <- strip_margins(colA_combined)
colB_combined <- strip_margins(colB_combined)
colC_combined <- strip_margins(colC_combined)

# how much horizontal gap you want (relative, smaller = tighter)
# --- kill vertical whitespace between a↕e↕i (and b↕f↕j, c↕g↕k) ---

squash_vertical_gaps <- function(plots, gap_mm = 0.5) {
  n <- length(plots)
  out <- vector("list", n)
  for (i in seq_along(plots)) {
    top_mm    <- if (i == 1) 0          else 0          # no top padding
    bottom_mm <- if (i == n) 0          else gap_mm     # tiny gap before next plot
    out[[i]] <- plots[[i]] +
      theme(
        plot.margin = margin(top_mm, 0.5, bottom_mm, 0.5, unit = "mm"),
        legend.box.margin = margin(0, 0, 0, 0, "mm"),
        legend.margin     = margin(0, 0, 0, 0, "mm")
      )
  }
  # stack with no extra spacing added by patchwork
  wrap_plots(out, ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
}

colA_combined <- squash_vertical_gaps(colA, gap_mm = 0.5)
colB_combined <- squash_vertical_gaps(colB, gap_mm = 0.5)
colC_combined <- squash_vertical_gaps(colC, gap_mm = 0.5)

# final grid across columns (unchanged)
combined <- (colA_combined | colB_combined | colC_combined) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

# FINAL grid: 3 columns (A | B | C), multiple rows (one per CSV)
# combined <- (colA_combined | colB_combined | colC_combined)

# row-wise sequential tags for 3 columns per row:
# keep your make_row_tags()
make_row_tags <- function(row_index) {
  start <- (row_index - 1) * 3 + 1
  letters[start + 0:2]
}

# build row-wise, then vectorize the matrix column-wise
row_tags <- lapply(seq_len(n_rows), make_row_tags)  # list of length n_rows, each c(a,b,c) etc.
tag_mat  <- do.call(rbind, row_tags)                # n_rows x 3: 
# row1: a b c
# row2: d e f
# row3: g h i
all_tags <- as.vector(tag_mat)                      # -> a, d, g, ..., b, e, h, ..., c, f, i


combined <- combined + plot_annotation(tag_levels = list(all_tags)) &
  theme(plot.tag = element_text(size = label_fsize, face = "bold"))

# ------------------------------- #
# Save (same sizing logic)        #
# ------------------------------- #

# use the same width/height computation you had before
w_total_cm <- (w_tiles_cm + panel_pad_cm) * 1.9 
h_total_cm <- n_rows * (h_tiles_cm + title_axis_pad + 1.2)

out_root <- if (length(roots) >= 1) roots[[1]] else dirname(csv_files[[1]])
png_path <- file.path(out_root, "heatmaps_ABC_multi_from_gridcsv.png")
pdf_path <- file.path(out_root, "heatmaps_ABC_multi_from_gridcsv.pdf")

ggsave(png_path, plot = combined, width = w_total_cm, height = h_total_cm,
       units = "cm", dpi = 300, bg = "white")

auto_crop_png(png_path, fuzz = 12)
ggsave(pdf_path, plot = combined, device = "pdf", width = w_total_cm, height = h_total_cm,
       units = "cm", bg = "white")
knitr::plot_crop(pdf_path, quiet = TRUE)

cat("\nSaved combined figure with 3 shared colorbars (one per column):\n- ",
    basename(png_path), "\n- ", basename(pdf_path), "\n", sep = "")
