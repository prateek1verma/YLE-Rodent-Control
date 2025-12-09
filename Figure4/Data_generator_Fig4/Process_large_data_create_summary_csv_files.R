############################################################
# heatmaps_ABC_combined.R  —  compact, flexible, robust
# - Reads param-set/run folders of MGDrivE CSVs
# - Saves per-run and per-parameter summaries (A–F)
# - Plots three heatmaps: A (Pr[elim]), B (TTE), C (Eq. pop norm)
############################################################

rm(list = ls()); gc(); save.image()

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, ggplot2, readr, stringr, tidyr, purrr
  library(scales)
  library(grid)        # unit()
  library(magick)      # optional: auto crop PNG
  library(knitr)       # pdf crop
})

# ------------------------------- #
# User settings (edit as needed)  #
# ------------------------------- #

# Root folder containing parameter-set subfolders (each with 001, 002, ... runs)
root <- "mgdrivefRIDL_sweep_10yrs_heatmap"
stopifnot(dir.exists(root))

# Biology / time scaling
Neq               <- 10000        # carrying-capacity scale for normalization
fem_elim_threshold<- 1            # <= this considered eliminated
days_per_year     <- 365
release_start_day <- round(8*days_per_year) # day of first release; TTE measured from here
years             <- 10           # used in x-axis label only

# Equilibrium window: average over last 'equil_window_days' days
equil_window_days <- 365

# Panel-B behavior (user-controlled)
tte_stat             <- "median"   # "median" or "mean"
plot_tte_prob_min    <- 0.50       # require Pr[elim] >= this to show TTE in panel B
mask_below_threshold <- TRUE       # if TRUE, tiles below threshold show NA (grey)

# Optional highlight band on TTE heatmap (NULL to disable)
tte_highlight_range_years <- c(0, 5)

# Visual helpers
lighten_hex <- function(col, amt = 0.88) rgb(t(col2rgb(col) * (1 - amt) + 255 * amt) / 255)
outline_colorB   <- "#D55E00"
outline_size     <- 1.5
outline_linetype <- "dotted"

# Plot-window filters (in percent). Use NULL to include all.
rel_filter_pct     <- NULL   # e.g., c(1, 10); or NULL
fy_cost_filter_pct <- NULL       # e.g., c(0, 90); or NULL

# Colorbar scaling: TRUE uses observed range after filtering; FALSE uses fixed fallback
rescale_to_filtered <- TRUE


# ------------------------------- #
# Helpers                         #
# ------------------------------- #

read_ts_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"Time" %in% names(df)) stop(sprintf("'Time' column not found in %s", basename(path)))
  df <- df %>% mutate(Time = suppressWarnings(as.numeric(Time))) %>% arrange(Time)
  if (anyNA(df$Time)) stop(sprintf("Non-numeric 'Time' values in %s", basename(path)))
  df
}


# --- FLEXIBLE TAG PARSERS ---

# returns numeric value for a given key if present; NA_real_ otherwise
# Parse number from a parameter tag "rel0p05" -> 0.05
parse_num_from_tag <- function(tag, key) {
  # matches: rel0p02, rel1p00, fy0p40, pq0p00, etc. (order-independent)
  m <- stringr::str_match(tag, paste0("(^|_)", key, "([-0-9p]+)"))[,3]
  if (is.na(m)) return(NA_real_)
  as.numeric(sub("p", ".", m, fixed = TRUE))
}

# parse all known keys; supply defaults as needed
parse_param_tags <- function(tag,
                             defaults = list(rel = NA_real_, fy = NA_real_, pq = NA_real_,
                                             fl = NA_real_, fs = NA_real_,
                                             mu = NA_real_, j = NA_real_, c = NA_real_)) {
  vals <- list(
    rel = parse_num_from_tag(tag, "rel"),
    fy  = parse_num_from_tag(tag, "fy"),
    pq  = parse_num_from_tag(tag, "pq"),
    fl  = parse_num_from_tag(tag, "fl"),
    fs  = parse_num_from_tag(tag, "fs"),
    mu  = parse_num_from_tag(tag, "mu"),
    j   = parse_num_from_tag(tag, "j"),
    c   = parse_num_from_tag(tag, "c")
  )
  for (k in names(vals)) {
    if (is.na(vals[[k]]) && !is.null(defaults[[k]])) vals[[k]] <- defaults[[k]]
  }
  tibble::tibble(
    rel = vals$rel, fy = vals$fy, pq = vals$pq, fl = vals$fl, fs = vals$fs,
    mu = vals$mu, j = vals$j, c = vals$c
  )
}


# Per-run metric extraction (returns single-row tibble)
# - final_female, final_male, final_total, final_total_norm
# - elim (logical), tte_years (from first release)
# - eq_mean_total_norm: mean(total_adults) over last equil_window_days, normalized by Neq
run_metrics <- function(run_dir, fem_thresh, release_start_day, days_per_year, Neq, equil_window_days) {
  f_path <- list.files(run_dir, pattern = "^F_Aggregate_Run\\d+_Patch\\d+\\.csv$", full.names = TRUE)
  m_path <- list.files(run_dir, pattern = "^M_Run\\d+_Patch\\d+\\.csv$",         full.names = TRUE)
  if (length(f_path) == 0 || length(m_path) == 0) {
    warning(sprintf("Missing F or M csv in: %s", run_dir)); return(tibble())
  }
  fdf <- read_ts_csv(f_path[1]); mdf <- read_ts_csv(m_path[1])
  
  f_cols <- grep("^f", names(fdf), value = TRUE)
  m_cols <- grep("^m", names(mdf), value = TRUE)
  
  female_total <- if ("female_total" %in% names(fdf)) fdf$female_total else {
    if (length(f_cols)) rowSums(fdf[, f_cols, drop = FALSE]) else rep(NA_real_, nrow(fdf))
  }
  male_total <- if (length(m_cols)) rowSums(mdf[, m_cols, drop = FALSE]) else rep(NA_real_, nrow(mdf))
  
  TS <- inner_join(
    fdf %>% select(Time, female_total = !!sym(if ("female_total" %in% names(fdf)) "female_total" else f_cols[1])) %>%
      mutate(female_total = as.numeric(female_total)),
    mdf %>% select(Time, male_total = !!sym(m_cols[1])) %>%
      mutate(male_total = as.numeric(rowSums(mdf[, m_cols, drop = FALSE]))),
    by = "Time"
  ) %>% arrange(Time) %>%
    mutate(total_adults = male_total + female_total)
  
  final_female <- suppressWarnings(as.numeric(tail(TS$female_total, 1)))
  final_male   <- suppressWarnings(as.numeric(tail(TS$male_total,   1)))
  final_total  <- suppressWarnings(as.numeric(tail(TS$total_adults, 1)))
  
  elim <- is.finite(final_female) && final_female <= fem_thresh
  
  tte_years <- NA_real_
  if (elim) {
    idx0 <- which(TS$Time >= release_start_day)
    if (length(idx0)) {
      idx0 <- idx0[1]
      hit <- which(TS$female_total[idx0:nrow(TS)] <= fem_thresh)
      if (length(hit)) {
        hit_idx <- idx0 + hit[1] - 1
        tte_years <- (TS$Time[hit_idx] - release_start_day) / days_per_year
        if (!is.finite(tte_years)) tte_years <- NA_real_
      }
    }
  }
  
  tmax <- max(TS$Time, na.rm = TRUE)
  win_start <- max(0, tmax - equil_window_days)
  eq_mean_total <- TS %>%
    filter(Time >= win_start) %>%
    summarise(mean_tot = mean(total_adults, na.rm = TRUE)) %>%
    pull(mean_tot)
  
  tibble(
    run_dir, run_id = basename(run_dir),
    final_female, final_male, final_total,
    final_total_norm = final_total / Neq,
    elim = as.logical(elim),
    tte_years = tte_years,
    eq_mean_total_norm = as.numeric(eq_mean_total) / Neq
  )
}

# Legend limit helper (robust)
get_limits_breaks <- function(vals, n = 5, fallback = c(0, 1)) {
  v <- vals[is.finite(vals)]
  lims <- if (length(v)) range(v, na.rm = TRUE) else fallback
  if (!is.finite(diff(lims)) || diff(lims) == 0) {
    eps <- max(1e-6, abs(lims[1]) * 1e-6, 1e-6)
    lims <- c(lims[1] - eps, lims[2] + eps)
  }
  breaks <- scales::breaks_pretty(n)(lims)
  breaks <- sort(unique(c(lims[1], breaks, lims[2])))
  list(lims = lims, breaks = breaks[breaks >= lims[1] & breaks <= lims[2]])
}

# Auto crop utilities
auto_crop_png <- function(infile, outfile = infile, fuzz = 10) {
  stopifnot(requireNamespace("magick", quietly = TRUE))
  img <- magick::image_read(infile) |> magick::image_trim(fuzz = fuzz)
  magick::image_write(img, path = outfile); invisible(outfile)
}


# ------------------------------- #
# Scan parameter-set folders      #
# ------------------------------- #

param_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
stopifnot(length(param_dirs) > 0)

# 1) Build per-run table (FULL DATA you can analyze variations from)
per_run_tbl <- purrr::map_dfr(param_dirs, function(pdir) {
  tag <- basename(pdir)
  
  # NEW: parse all available tags, with defaults if you set them
  pars <- parse_param_tags(tag)
  rel <- pars$rel; fy <- pars$fy; pq <- pars$pq; fl <- pars$fl; fs <- pars$fs
  mu  <- pars$mu;  j  <- pars$j;  cval <- pars$c
  
  run_dirs <- list.dirs(pdir, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("/\\d{3}$", run_dirs)]
  if (!length(run_dirs)) return(tibble())
  
  purrr::map_dfr(run_dirs, function(rd) {
    rm <- run_metrics(rd, fem_elim_threshold, release_start_day, days_per_year, Neq, equil_window_days)
    if (!nrow(rm)) return(tibble())
    rm %>%
      mutate(
        tag = tag,
        rel = rel, fy = fy, pq = pq, fl = fl, fs = fs,
        mu = mu, j = j, c = cval
      )
  })
})


# Save per-run CSV (FULL data)
per_run_csv <- file.path(root, "heatmap_runs_per_run_full.csv")
if (nrow(per_run_tbl)) readr::write_csv(per_run_tbl, per_run_csv, na = "")

# 2) Summarize to per-parameter rows (A–F)
summ_tbl <- per_run_tbl %>%
  mutate(
    rel_pct     = 100 * rel,
    fy_cost_pct = pmax(0, pmin(100, 100 * (1 - fy)))
  ) %>%
  group_by(tag, rel, fy, pq, fl, fs, mu, j, c) %>%
  summarise(
    n_runs = n(),
    A_prob_elim          = mean(elim, na.rm = TRUE),                         # (A)
    B_mean_tte_cond      = if (any(elim)) mean(tte_years[elim], na.rm = TRUE) else NA_real_,  # (B)
    C_median_tte_cond    = if (any(elim)) median(tte_years[elim], na.rm = TRUE) else NA_real_,# (C)
    D_eq_final_norm_mean = mean(eq_mean_total_norm, na.rm = TRUE),           # (D) equilibrium mean, normalized
    E_sd_tte_cond        = if (sum(elim, na.rm = TRUE) > 1) sd(tte_years[elim], na.rm = TRUE) else NA_real_, # (E)
    F_iqr_tte_cond       = if (any(elim)) IQR(tte_years[elim], na.rm = TRUE) else NA_real_,    # (F)
    .groups = "drop"
  ) %>%
  mutate(
    rel_pct     = 100 * rel,
    fy_cost_pct = pmax(0, pmin(100, 100 * (1 - fy))),
    rel_lab     = sprintf("%d", as.integer(round(rel_pct))),
    fy_cost_lab = sprintf("%d", as.integer(round(fy_cost_pct)))
  )

# Column chosen for panel B (respect tte_stat and prob threshold)
summ_tbl <- summ_tbl %>%
  mutate(
    tte_stat_choice = match.arg(tte_stat, c("median", "mean")),
    tte_to_plot_raw = ifelse(tte_stat_choice == "median", C_median_tte_cond, B_mean_tte_cond),
    tte_to_plot = if (mask_below_threshold) ifelse(A_prob_elim >= plot_tte_prob_min, tte_to_plot_raw, NA_real_) else tte_to_plot_raw
  )

# ------------------------------- #
# Filtering for plotting          #
# ------------------------------- #

filtered_tbl <- summ_tbl %>%
  filter(is.finite(rel), is.finite(fy), n_runs > 0)

if (!is.null(rel_filter_pct)) {
  stopifnot(length(rel_filter_pct) == 2)
  filtered_tbl <- filtered_tbl %>%
    filter(rel_pct >= min(rel_filter_pct), rel_pct <= max(rel_filter_pct))
}
if (!is.null(fy_cost_filter_pct)) {
  stopifnot(length(fy_cost_filter_pct) == 2)
  filtered_tbl <- filtered_tbl %>%
    filter(fy_cost_pct >= min(fy_cost_filter_pct), fy_cost_pct <= max(fy_cost_filter_pct))
}
stopifnot(nrow(filtered_tbl) > 0)

# Ordered factor labels and integer grid coordinates
rel_levels     <- filtered_tbl %>% distinct(rel_pct, rel_lab) %>% arrange(rel_pct) %>% pull(rel_lab)
fy_cost_levels <- filtered_tbl %>% distinct(fy_cost_pct, fy_cost_lab) %>% arrange(fy_cost_pct) %>% pull(fy_cost_lab)

filtered_tbl <- filtered_tbl %>%
  mutate(
    x_idx = as.integer(factor(rel_lab,     levels = rel_levels)),
    y_idx = as.integer(factor(fy_cost_lab, levels = fy_cost_levels))
  )

# ------------------------------- #
# Optional perimeter for TTE band #
# ------------------------------- #

perimeter_segments <- NULL
if (!is.null(tte_highlight_range_years)) {
  lo <- min(tte_highlight_range_years); hi <- max(tte_highlight_range_years)
  
  grid_df <- tidyr::expand_grid(rel_lab = rel_levels, fy_cost_lab = fy_cost_levels) %>%
    left_join(filtered_tbl, by = c("rel_lab", "fy_cost_lab")) %>%
    mutate(
      x_idx = as.integer(factor(rel_lab, levels = rel_levels)),
      y_idx = as.integer(factor(fy_cost_lab, levels = fy_cost_levels)),
      mask  = ifelse(is.finite(tte_to_plot) & tte_to_plot >= lo & tte_to_plot <= hi, 1, 0)
    )
  
  inside_df <- grid_df %>% filter(mask == 1) %>% select(x_idx, y_idx)
  
  left_edges <- inside_df %>%
    anti_join(inside_df %>% transmute(x_idx = x_idx + 1, y_idx), by = c("x_idx", "y_idx")) %>%
    transmute(x0 = x_idx - 0.5, y0 = y_idx - 0.5, x1 = x_idx - 0.5, y1 = y_idx + 0.5)
  right_edges <- inside_df %>%
    anti_join(inside_df %>% transmute(x_idx = x_idx - 1, y_idx), by = c("x_idx", "y_idx")) %>%
    transmute(x0 = x_idx + 0.5, y0 = y_idx - 0.5, x1 = x_idx + 0.5, y1 = y_idx + 0.5)
  bottom_edges <- inside_df %>%
    anti_join(inside_df %>% transmute(x_idx, y_idx = y_idx + 1), by = c("x_idx", "y_idx")) %>%
    transmute(x0 = x_idx - 0.5, y0 = y_idx - 0.5, x1 = x_idx + 0.5, y1 = y_idx - 0.5)
  top_edges <- inside_df %>%
    anti_join(inside_df %>% transmute(x_idx, y_idx = y_idx - 1), by = c("x_idx", "y_idx")) %>%
    transmute(x0 = x_idx - 0.5, y0 = y_idx + 0.5, x1 = x_idx + 0.5, y1 = y_idx + 0.5)
  
  perimeter_segments <- bind_rows(left_edges, right_edges, bottom_edges, top_edges) %>% distinct()
}

# ------------------------------- #
# Sizing & themes                 #
# ------------------------------- #

nx <- length(rel_levels); ny <- length(fy_cost_levels)
compute_tile_size_cm <- function(nx, ny, min_cm = 0.9, max_cm = 2.8) {
  base <- 2.0; adj <- base * (12 / max(12, max(nx, ny)))
  pmin(max_cm, pmax(min_cm, adj))
}
tile_size_cm <- compute_tile_size_cm(nx, ny)
legend_h_cm  <- max(6, min(14, ny * tile_size_cm * 0.8))
legend_w_cm  <- 1.0
panel_pad_cm <- 2.5
title_axis_pad <- 4.5

w_tiles_cm <- nx * tile_size_cm; h_tiles_cm <- ny * tile_size_cm
w_single <- w_tiles_cm + panel_pad_cm + legend_w_cm
h_single <- h_tiles_cm + title_axis_pad
w_total_cm <- w_single * 3 + 2.0
h_total_cm <- h_single
fixed_ratio <- nx / ny

axis_text_x <- element_text(size = if (nx > 18) 9 else if (nx > 12) 10 else 14, angle = 0)
axis_text_y <- element_text(size = if (ny > 18) 9 else 14, angle = 0)

base_theme <- theme_minimal(base_size = 17) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.x = axis_text_x,
        axis.text.y = axis_text_y,
        plot.title = element_text(face = "bold"),
        legend.position = "right",
        plot.margin = margin(6, 6, 6, 6))

# Color scale bounds
sb_elim <- get_limits_breaks(filtered_tbl$A_prob_elim, n = 5, fallback = c(0, 1))
sb_tte  <- get_limits_breaks(filtered_tbl$tte_to_plot, n = 5,
                             fallback = range(filtered_tbl$tte_to_plot, na.rm = TRUE))
sb_tot  <- get_limits_breaks(filtered_tbl$D_eq_final_norm_mean, n = 5, fallback = c(0, 1))

# ------------------------------- #
# Panels                          #
# ------------------------------- #

# Panel A: Probability of female elimination
pA <- ggplot(filtered_tbl, aes(x = x_idx, y = y_idx, fill = A_prob_elim)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), breaks = seq_along(rel_levels), labels = rel_levels) +
  scale_y_continuous(expand = c(0, 0), breaks = seq_along(fy_cost_levels), labels = fy_cost_levels) +
  scale_fill_gradient(
    low = lighten_hex("#006d2c", 0.90), high = "#006d2c",
    name = "Probability of\nfemale elimination",
    limits = if (rescale_to_filtered) sb_elim$lims else c(0,1),
    breaks = sb_elim$breaks,
    labels = label_percent(accuracy = 1),
    oob = squish, na.value = "grey90"
  ) +
  labs(
    title = "Female elimination probability",
    x = sprintf("Release proportion of y-males\n(%% of adult males at K, per month for %d years)", years),
    y = "Fitness cost in lifespan reduction of y-males (%)"
  ) +
  coord_fixed(ratio = fixed_ratio, expand = FALSE, clip = "on") + base_theme +
  guides(fill = guide_colorbar(direction = "vertical",
                               barheight = unit(legend_h_cm, "cm"),
                               barwidth  = unit(legend_w_cm, "cm"),
                               ticks.colour = "black"))

# Panel B: TTE (median or mean), masked by Pr[elim] threshold if chosen
tte_label <- if (tte_stat == "median") "Median elimination\ntime (yrs)" else "Mean elimination\ntime (yrs)"
pB <- ggplot(filtered_tbl, aes(x = x_idx, y = y_idx, fill = tte_to_plot)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), breaks = seq_along(rel_levels), labels = rel_levels) +
  scale_y_continuous(expand = c(0, 0), breaks = seq_along(fy_cost_levels), labels = fy_cost_levels) +
  scale_fill_gradient(
    low = lighten_hex("#D55E00", 0.90), high = "#D55E00",
    name   = tte_label,
    limits = sb_tte$lims,
    breaks = sb_tte$breaks,
    labels = label_number(accuracy = 0.1),
    oob    = squish,
    na.value = "grey90"
  ) +
  labs(
    title = sprintf("Time to female elimination (%s%s)",
                    tte_stat,
                    if (mask_below_threshold) sprintf(", when Pr >= %.0f%%", 100*plot_tte_prob_min) else ""),
    x = sprintf("Release proportion of y-males\n(%% of adult males at K, per month for %d years)", years),
    y = "Fitness cost in lifespan reduction of y-males (%)"
  ) +
  coord_fixed(ratio = fixed_ratio, expand = FALSE, clip = "on") + base_theme +
  guides(fill = guide_colorbar(direction = "vertical",
                               barheight = unit(legend_h_cm, "cm"),
                               barwidth  = unit(legend_w_cm, "cm"),
                               ticks.colour = "black"))

# Optional perimeter outline for B
if (!is.null(tte_highlight_range_years) && !is.null(perimeter_segments) && nrow(perimeter_segments) > 0) {
  pB <- pB + geom_segment(data = perimeter_segments,
                          aes(x = x0, y = y0, xend = x1, yend = y1),
                          inherit.aes = FALSE,
                          linewidth   = outline_size,
                          lineend     = "square",
                          color       = outline_colorB,
                          linetype    = outline_linetype)
}

# Panel C: Equilibrium normalized population (mean across runs)
pC <- ggplot(filtered_tbl, aes(x = x_idx, y = y_idx, fill = D_eq_final_norm_mean)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), breaks = seq_along(rel_levels), labels = rel_levels) +
  scale_y_continuous(expand = c(0, 0), breaks = seq_along(fy_cost_levels), labels = fy_cost_levels) +
  scale_fill_gradient(
    low = lighten_hex("#0072B2", 0.90), high = "#0072B2",
    name   = "Normalized population\n(mean across runs)",
    limits = if (rescale_to_filtered) sb_tot$lims else c(0, 1),
    breaks = sb_tot$breaks,
    labels = label_number(accuracy = 0.1),
    oob    = squish,
    na.value = "grey90"
  ) +
  labs(
    title = "Total population after 50 yrs",
    x = sprintf("Release proportion of y-males\n(%% of adult males at K, per month for %d years)", years),
    y = "Fitness cost in lifespan reduction of y-males (%)"
  ) +
  coord_fixed(ratio = fixed_ratio, expand = FALSE, clip = "on") + base_theme +
  guides(fill = guide_colorbar(direction = "vertical",
                               barheight = unit(legend_h_cm, "cm"),
                               barwidth  = unit(legend_w_cm, "cm"),
                               ticks.colour = "black"))

# ------------------------------- #
# Arrange & Save                  #
# ------------------------------- #

.use_patchwork <- requireNamespace("patchwork", quietly = TRUE)
.use_cowplot   <- if (.use_patchwork) FALSE else requireNamespace("cowplot", quietly = TRUE)

if (.use_patchwork) {
  combined <- (pA | pB | pC) + patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 32, face = "bold"))
  print(combined)
} else if (.use_cowplot) {
  combined <- cowplot::plot_grid(pA, pB, pC, labels = c("a", "b", "c"), label_size = 32, nrow = 1, align = "hv")
  print(combined)
} else {
  warning("Install 'patchwork' or 'cowplot' to display the three heatmaps side-by-side with labels.")
  combined <- NULL
  print(pA); print(pB); print(pC)
}

# Save figures
png_path <- file.path(root, "heatmaps_ABC_combined.png")
pdf_path <- file.path(root, "heatmaps_ABC_combined.pdf")

if (!is.null(combined)) {
  ggsave(png_path, plot = combined, width = w_total_cm, height = h_total_cm, units = "cm", dpi = 600, bg = "white")
  auto_crop_png(png_path, fuzz = 12)
  
  ggsave(pdf_path, plot = combined, device = "pdf", width = w_total_cm, height = h_total_cm, units = "cm", bg = "white")
  knitr::plot_crop(pdf_path, quiet = TRUE)
}

# ------------------------------- #
# Save CSVs                       #
# ------------------------------- #

# FULL per-run table (variations readily analyzable)
readr::write_csv(per_run_tbl, file.path(root, "heatmap_runs_per_run_full.csv"), na = "")

# Parameter-level summary (A–F)
summary_param <- summ_tbl %>%
  transmute(
    tag, rel, fy, pq, fl, fs, mu, j, c,
    rel_pct = round(rel_pct, 2),
    fy_cost_pct = round(fy_cost_pct, 2),
    n_runs,
    A_prob_elim = A_prob_elim,             # (A)
    B_mean_time_years = B_mean_tte_cond,   # (B)
    C_median_time_years = C_median_tte_cond,                   # (C)
    D_final_norm_pop_equil = D_eq_final_norm_mean,             # (D)
    E_sd_mean_time_years = E_sd_tte_cond,                      # (E)
    F_iqr_time_years = F_iqr_tte_cond                          # (F)
  )
readr::write_csv(summary_param, file.path(root, "heatmap_summary_ABCDEF_per_param.csv"), na = "")

# Filtered grid used for plotting (with tile indices & chosen TTE stat)
summary_filtered <- filtered_tbl %>%
  select(tag, rel, fy, pq, fl, fs, mu, j, c,
         rel_pct, fy_cost_pct, rel_lab, fy_cost_lab, x_idx, y_idx,
         n_runs,
         A_prob_elim,
         B_mean_tte_cond, C_median_tte_cond, tte_stat_choice, tte_to_plot,
         D_eq_final_norm_mean)
readr::write_csv(summary_filtered, file.path(root, "heatmap_summary_plot_grid.csv"), na = "")

# Optional: long format for quick faceting or external plotting
summary_long <- summary_filtered %>%
  pivot_longer(cols = c(A_prob_elim, tte_to_plot, D_eq_final_norm_mean),
               names_to = "metric", values_to = "value")
readr::write_csv(summary_long, file.path(root, "heatmap_summary_plot_grid_long.csv"), na = "")

cat("\nSaved:\n",
    "- per-run: heatmap_runs_per_run_full.csv\n",
    "- per-parameter: heatmap_summary_ABCDEF_per_param.csv\n",
    "- plot-grid (filtered): heatmap_summary_plot_grid.csv (+ _long)\n", sep = "")
