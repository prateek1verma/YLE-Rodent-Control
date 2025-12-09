rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ------------------------------- #
# User settings                   #
# ------------------------------- #

root            <- "mgdrive_YLE_10yrs"   # parent folder with param-tag subfolders
years_window    <- 10                       # compute min/max over first N years
days_per_year   <- 365
fy_keep         <- 1                      # keep only fy ≈ 0.5
runs_expect     <- sprintf("%03d", 1:10)     # expected run folders
burn_in_years  <- 8
analysis_years <- 50
days_per_year  <- 365   # keep as before

# ------------------------------- #
# Helpers (your originals)        #
# ------------------------------- #

read_ts_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"Time" %in% names(df)) stop(sprintf("'Time' column not found in %s", basename(path)))
  df <- df %>% mutate(Time = suppressWarnings(as.numeric(Time))) %>% arrange(Time)
  if (anyNA(df$Time)) stop(sprintf("Non-numeric 'Time' values in %s", basename(path)))
  df
}

parse_num_from_tag <- function(tag, key) {
  m <- stringr::str_match(tag, paste0("(^|_)", key, "([-0-9p]+)"))[,3]
  if (is.na(m)) return(NA_real_)
  as.numeric(sub("p", ".", m, fixed = TRUE))
}

parse_param_tags <- function(tag,
                             defaults = list(rel = NA_real_, fy = NA_real_, pq = NA_real_,
                                             fl = NA_real_, mu = NA_real_, j = NA_real_, c = NA_real_)) {
  vals <- list(
    rel = parse_num_from_tag(tag, "rel"),
    fy  = parse_num_from_tag(tag, "fy"),
    pq  = parse_num_from_tag(tag, "pq"),
    fl  = parse_num_from_tag(tag, "fl"),
    mu  = parse_num_from_tag(tag, "mu"),
    j   = parse_num_from_tag(tag, "j"),
    c   = parse_num_from_tag(tag, "c")
  )
  for (k in names(vals)) if (is.na(vals[[k]]) && !is.null(defaults[[k]])) vals[[k]] <- defaults[[k]]
  tibble::tibble(rel = vals$rel, fy = vals$fy, pq = vals$pq, fl = vals$fl,
                 mu = vals$mu, j = vals$j, c = vals$c)
}

# ------------------------------- #
# File readers for totals         #
# ------------------------------- #

# --- READER: females + males + total ----------------------------------------
read_fem_male_totals <- function(run_dir) {
  f_path <- list.files(run_dir, pattern = "^F_.*Aggregate.*\\.csv$", full.names = TRUE)
  m_path <- list.files(run_dir, pattern = "^M_.*\\.csv$",               full.names = TRUE)
  
  if (length(f_path) == 0 || length(m_path) == 0) {
    stop(sprintf("Missing F or M csv in: %s (need both to compute female minima)", run_dir))
  }
  fdf <- read_ts_csv(f_path[1]); mdf <- read_ts_csv(m_path[1])
  
  f_cols <- grep("^f", names(fdf), value = TRUE)
  m_cols <- grep("^m", names(mdf), value = TRUE)
  
  female_total <- if ("female_total" %in% names(fdf)) fdf$female_total else {
    if (length(f_cols)) rowSums(fdf[, f_cols, drop = FALSE]) else rep(NA_real_, nrow(fdf))
  }
  male_total <- if (length(m_cols)) rowSums(mdf[, m_cols, drop = FALSE]) else rep(NA_real_, nrow(mdf))
  
  df <- dplyr::inner_join(
    fdf %>% dplyr::select(Time) %>% dplyr::mutate(female_total = as.numeric(female_total)),
    mdf %>% dplyr::select(Time) %>% dplyr::mutate(male_total   = as.numeric(male_total)),
    by = "Time"
  ) %>% dplyr::arrange(Time) %>%
    dplyr::mutate(total_adults = female_total + male_total)
  
  df %>% dplyr::select(Time, female_total, total_adults)
}

# --- Per-run extrema over first N years -------------------------------------
run_extrema_10y <- function(run_dir) {
  df <- read_fem_male_totals(run_dir)
  cutoff <- years_window * days_per_year
  df10 <- dplyr::filter(df, Time <= cutoff)
  tibble::tibble(
    run_dir        = run_dir,
    run_id         = basename(run_dir),
    max_total_10y  = max(df10$total_adults, na.rm = TRUE),   # series 1 (unchanged)
    min_female_10y = min(df10$female_total, na.rm = TRUE)    # series 2 (NEW)
  )
}

# --- Per-run extrema over [burn-in, burn-in+analysis] -----------------------
run_extrema_window <- function(run_dir) {
  df <- read_fem_male_totals(run_dir)   # from the previous version I gave you
  t0 <- burn_in_years  * days_per_year
  t1 <- (burn_in_years + analysis_years) * days_per_year
  dfw <- dplyr::filter(df, Time >= t0, Time <= t1)
  tibble::tibble(
    run_dir        = run_dir,
    run_id         = basename(run_dir),
    max_total_10y  = max(dfw$total_adults,  na.rm = TRUE),
    min_female_10y = min(dfw$female_total, na.rm = TRUE)
  )
}


# --- Build per-run table (uses your existing param scan above) ---------------
# (Re-run the per_run_tbl construction from your script after redefining run_extrema_10y)
stopifnot(dir.exists(root))
param_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
stopifnot(length(param_dirs) > 0)

per_run_tbl <- purrr::map_dfr(param_dirs, function(pdir) {
  tag <- basename(pdir)
  pars <- parse_param_tags(tag)
  if (!is.finite(pars$fy) || abs(pars$fy - fy_keep) > 1e-8) return(tibble::tibble())
  
  run_dirs <- list.dirs(pdir, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[basename(run_dirs) %in% runs_expect]
  if (!length(run_dirs)) return(tibble::tibble())
  
  purrr::map_dfr(run_dirs, function(rd) {
    # ext <- tryCatch(run_extrema_10y(rd), error = function(e) tibble::tibble())
    ext <- tryCatch(run_extrema_window(rd), error = function(e) tibble::tibble())
    if (!nrow(ext)) return(tibble::tibble())
    ext %>% dplyr::mutate(tag = tag, rel = pars$rel, fy = pars$fy)
  })
})

stopifnot(nrow(per_run_tbl) > 0)

# --- Summaries (rel as %, two stats: max total, min female) ------------------
per_run_long <- per_run_tbl %>%
  tidyr::pivot_longer(c(max_total_10y, min_female_10y),
                      names_to = "stat", values_to = "value") %>%
  dplyr::mutate(
    stat = dplyr::recode(stat,
                         max_total_10y  = "max_total",
                         min_female_10y = "min_female"),
    rel_pct = 100*rel,
    stat = factor(stat, levels = c("max_total", "min_female"))
  )

summ_tbl <- per_run_long %>%
  dplyr::group_by(rel_pct, stat) %>%
  dplyr::summarise(
    n    = dplyr::n(),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    se   = sd / sqrt(n),
    tval = qt(0.975, df = pmax(n - 1, 1)),
    lo   = mean - tval * se,
    hi   = mean + tval * se,
    .groups = "drop"
  )

# --- Plots: stacked vertically ----------------------------------------------
# --- Dual-axis plot: max total (left) & min females (right) ------------------

# split summaries
summ_max  <- dplyr::filter(summ_tbl, stat == "max_total")
summ_minF <- dplyr::filter(summ_tbl, stat == "min_female")

# linear mapping: y_right -> y_left  so curves share one panel
rng1   <- range(c(summ_max$lo,  summ_max$hi),  na.rm = TRUE)  # left axis span
rng2   <- range(c(summ_minF$lo, summ_minF$hi), na.rm = TRUE)  # right axis span
scale  <- diff(rng1) / diff(rng2)
offset <- rng1[1] - rng2[1] * scale

summ_minF_tr <- summ_minF %>%
  dplyr::mutate(
    mean_tr = mean * scale + offset,
    lo_tr   = lo   * scale + offset,
    hi_tr   = hi   * scale + offset
  )

xr <- range(summ_tbl$rel_pct, na.rm = TRUE)

# --- Font size handles -------------------------------------------------------
axis_title_size  <- 24  # x/y titles
axis_text_size   <- 24  # tick labels
legend_text_size <- 24
legend_title_size<- 24
title_size       <- 18  # overall plot title (optional)

# --- Dual-axis plot with shaded 95% CIs (CORRECTED START AT 10K) ---
summ_max  <- dplyr::filter(summ_tbl, stat == "max_total")
summ_minF <- dplyr::filter(summ_tbl, stat == "min_female")

xr <- range(summ_tbl$rel_pct, na.rm = TRUE)

# helper for "k" labels: 0, 25k, 50k, ...
k_lab <- function(x) paste0(round(x/1000, 1), "k")

# ==============================================================================
# NEW SCALING LOGIC: Map Right Axis [0, 5050] -> Left Axis [10000, Max_Red]
# ==============================================================================

# 1. Define where you want the Left Axis to start
y_left_min <- 10000

# 2. Define the max values for scaling
max_y_left  <- max(summ_max$mean, na.rm = TRUE)
# Use 5050 as the "ceiling" for the blue line (or use max(summ_minF$mean))
max_y_right <- 5050 

# 3. Calculate Slope (m)
# We want: 0 (Right) -> 10000 (Left)
#          max_y_right (Right) -> max_y_left (Left)
slope <- (max_y_left - y_left_min) / max_y_right

# 4. Apply Transformation: y_new = y_old * slope + offset
summ_minF_tr <- summ_minF %>%
  dplyr::mutate(
    mean_tr = mean * slope + y_left_min,
    lo_tr   = lo   * slope + y_left_min,
    hi_tr   = hi   * slope + y_left_min
  )

# 2. Plotting
p_dual <- ggplot() +
  # --- Reference Line ---
  geom_hline(yintercept = 10000, linetype = "dashed", color = "grey30") +
  
  # =======================================================
# LEFT SERIES: Max Total (Red)
# =======================================================
# 1. Error Bars
geom_errorbar(data = summ_max, 
              aes(x = rel_pct, ymin = lo, ymax = hi, color = "Max total"), 
              width = 0.5, linewidth = 0.6) + # width controls horizontal cap size
  # 2. Line
  geom_line(data = summ_max, 
            aes(x = rel_pct, y = mean, color = "Max total"), 
            linewidth = 1) +
  # 3. Points (plotted last to sit on top)
  geom_point(data = summ_max, 
             aes(x = rel_pct, y = mean, color = "Max total"), 
             size = 3) +
  
  # =======================================================
# RIGHT SERIES: Min Females (Blue) - Transformed
# =======================================================
# 1. Error Bars (using lo_tr and hi_tr)
geom_errorbar(data = summ_minF_tr, 
              aes(x = rel_pct, ymin = lo_tr, ymax = hi_tr, color = "Min females"), 
              width = 0.5, linewidth = 0.6) +
  # 2. Line
  geom_line(data = summ_minF_tr, 
            aes(x = rel_pct, y = mean_tr, color = "Min females"), 
            linewidth = 1) +
  # 3. Points
  geom_point(data = summ_minF_tr, 
             aes(x = rel_pct, y = mean_tr, color = "Min females"), 
             size = 3, shape = 17) +
  
  # =======================================================
# AXES & SCALES
# =======================================================
scale_x_continuous(
  name   = "Release size (% of adult sterile males released)",
  breaks = scales::breaks_pretty(n = 6),
  labels = function(x) sprintf("%.0f", x),
  limits = xr
) +
  scale_y_continuous(
    name   = "Peak value of total population",
    limits = c(y_left_min, NA),
    labels = k_lab,
    sec.axis = sec_axis(
      ~ (. - y_left_min) / slope,
      name   = "Minimum achieved total females",
      labels = k_lab
    )
  ) +
  scale_color_manual(
    NULL,
    values = c("Max total" = "#e41a1c", "Min females" = "#377eb8")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.margin           = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm"),
    legend.position       = c(0.80, 0.50),
    legend.background     = element_rect(fill = alpha("white", 0.2), colour = NA),
    legend.text           = element_text(size = legend_text_size),
    legend.title          = element_text(size = legend_title_size),
    axis.title.x          = element_text(size = axis_title_size, color = "black"),
    axis.title.y.left     = element_text(size = axis_title_size, color = "black", margin = margin(r = 5)),
    axis.title.y.right    = element_text(size = axis_title_size, color = "black"),
    axis.text.x           = element_text(size = axis_text_size),
    axis.text.y.left      = element_text(size = axis_text_size),
    axis.text.y.right     = element_text(size = axis_text_size),
    panel.border          = element_rect(colour = "black", fill = NA, size = 0.8),
    axis.line             = element_line(colour = "black"),
    axis.ticks            = element_line(colour = "black"),
    axis.ticks.length     = unit(3, "mm")
  )

print(p_dual)

## Optional save
out_dir <- root
ggsave(file.path(out_dir, "FigG.png"), p_dual, width = 10, height = 6, dpi = 300)



plot_tbl <- summ_tbl %>%
  dplyr::mutate(
    # y used on left axis
    y_left  = dplyr::if_else(stat == "max_total", mean, mean * scale),
    lo_left = dplyr::if_else(stat == "max_total", lo,   lo   * scale),
    hi_left = dplyr::if_else(stat == "max_total", hi,   hi   * scale),
    # original scale for the right axis (only meaningful for min_female)
    y_right = dplyr::if_else(stat == "min_female", mean, NA_real_)
  )

readr::write_csv(
  plot_tbl,
  file.path(out_dir, "YLE_extrema_dual_axis_data_all_scenarios.csv")
)
