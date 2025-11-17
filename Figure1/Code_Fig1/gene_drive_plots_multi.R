# ---- gene_drive_plots_multi.R ----
# Usage:
#   source("gene_drive_plots_multi.R")
#   make_run_patch_figure_single(run_id = "001", patch_id = "001", tstart = 3000)
# or run via Rscript after adding a call at the bottom.

suppressPackageStartupMessages({
  library(tidyverse)     # ggplot2, dplyr, readr, tidyr
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(dplyr)  # after any Bioconductor packages that define select()
})


# -----------------------------#
# Discover runs in a directory
# -----------------------------#
list_run_ids <- function(base_dir) {
  base_dir <- path.expand(base_dir)
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  # Accept "RunNNN" or "NNN"
  run_ids <- unique(c(
    sub("^Run", "", dirs[grepl("^Run\\d+$", dirs)]),
    dirs[grepl("^\\d+$", dirs)]
  ))
  if (length(run_ids) == 0) {
    stop(sprintf("No run folders found in %s (expecting 'NNN' or 'RunNNN').", base_dir))
  }
  sort(run_ids)
}

# -----------------------------------------------------#
# Build long series for *multiple* runs in one tibble
# -----------------------------------------------------#
build_series_long_multi <- function(base_dir, run_ids, patch_id, tstart) {
  # Use your existing per-run builder under the hood
  purrr::map_dfr(run_ids, function(rid) {
    built <- build_series_long(base_dir, rid, patch_id, tstart)
    built$series_long %>%
      mutate(
        run_id = rid,
        Years  = Time_rel / 365
      )
  })
}

#-----------------------------#
# I/O + path helpers
#-----------------------------#

# Read robustly (allow mixed types, ensure Time numeric), sorted by Time
read_ts_csv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"Time" %in% names(df)) {
    stop(sprintf("'Time' column not found in %s", path))
  }
  df <- df %>%
    mutate(Time = suppressWarnings(as.numeric(Time))) %>%
    arrange(Time)
  if (anyNA(df$Time)) {
    stop(sprintf("Non-numeric values found in 'Time' column of %s", path))
  }
  df
}

# Resolve whether the run folder is named "001" or "Run001"
resolve_run_dir <- function(base_dir, run_id) {
  base_dir <- path.expand(base_dir)
  candidates <- file.path(base_dir, c(run_id, sprintf("Run%s", run_id)))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop(
      sprintf(
        "Cannot find run folder. Tried:\n  %s\n  %s",
        candidates[1], candidates[2]
      )
    )
  }
  hit[1]
}

#-----------------------------#
# Summaries
#-----------------------------#

# Summarize male counts and (optionally) frequencies by Y vs y
summarize_males <- function(df) {
  maleY_cols <- grep("^mY", names(df), value = TRUE)
  maley_cols <- grep("^my", names(df), value = TRUE)
  
  if (length(maleY_cols) == 0 && length(maley_cols) == 0) {
    stop("No male genotype columns starting with 'mY' or 'my' found.")
  }
  
  df %>%
    mutate(
      male_Y     = if (length(maleY_cols)) rowSums(across(all_of(maleY_cols))) else 0,
      male_y     = if (length(maley_cols)) rowSums(across(all_of(maley_cols))) else 0,
      male_total = male_Y + male_y,
      freq_Y     = ifelse(male_total > 0, male_Y / male_total, NA_real_),
      freq_y     = ifelse(male_total > 0, male_y / male_total, NA_real_)
    ) %>%
    dplyr::select(Time, male_Y, male_y, male_total, freq_Y, freq_y)
}

# Summarize female totals from any fX* columns
summarize_females_from_df <- function(df) {
  female_cols <- grep("^fX", names(df), value = TRUE)
  if (length(female_cols) == 0) {
    stop("No female genotype columns starting with 'fX' found.")
  }
  df %>%
    transmute(
      Time,
      female_total = rowSums(across(all_of(female_cols)))
    )
}

# Try reading F_Aggregate_ file; otherwise fall back to female columns in the male CSV
get_female_totals <- function(base_dir, run_id, patch_id, male_df) {
  run_dir <- resolve_run_dir(base_dir, run_id)
  fpath   <- file.path(run_dir, sprintf("F_Aggregate_Run%s_Patch%s.csv", run_id, patch_id))
  
  if (file.exists(fpath)) {
    summarize_females_from_df(read_ts_csv(fpath))
  } else {
    summarize_females_from_df(male_df)
  }
}

# Join summarized males and females; add overall
join_male_female <- function(male_sum, female_sum) {
  inner_join(male_sum, female_sum, by = "Time") %>%
    mutate(overall_total = male_total + female_total)
}

#-----------------------------#
# Time filtering / reindexing
#-----------------------------#

# Filter to Time >= tstart and reindex to Time_rel = Time - tstart (x-axis starts at 0)
filter_and_reindex_time <- function(df, tstart) {
  out <- df %>%
    filter(Time >= tstart) %>%
    mutate(Time_rel = Time - tstart)
  
  if (nrow(out) == 0) {
    stop(sprintf("No rows at or after tstart = %s. Check your 'Time' range.", tstart))
  }
  out
}

# Build long-format series for plotting (male Y, male y, female total), with Time_rel
build_series_long <- function(base_dir, run_id, patch_id, tstart) {
  run_dir   <- resolve_run_dir(base_dir, run_id)
  male_path <- file.path(run_dir, sprintf("M_Run%s_Patch%s.csv", run_id, patch_id))
  
  male_raw   <- read_ts_csv(male_path)
  male_sum   <- summarize_males(male_raw)
  female_sum <- get_female_totals(base_dir, run_id, patch_id, male_raw)
  
  joined <- join_male_female(male_sum, female_sum) %>%
    filter_and_reindex_time(tstart = tstart)
  
  series_long <- joined %>%
    transmute(
      Time_rel,
      `male Y`       = male_Y,
      `male y`       = male_y,
      `female total` = female_total
    ) %>%
    pivot_longer(-Time_rel, names_to = "group", values_to = "count")
  
  list(
    run_dir     = run_dir,
    male_path   = male_path,
    data_joined = joined,
    series_long = series_long
  )
}


# -----------------------------------------------------#
# Plot: many thin trajectories + thick mean trajectory
# -----------------------------------------------------#
make_patch_figure_multi <- function(
    base_dir = "~/Library/CloudStorage/Box-Box/Research Projects/Rodent GD/MGDrivE-master/ExampleCode/mgdriveYLE_2patch",
    patch_id = "001",
    run_ids  = NULL,        # if NULL, auto-discover all runs in base_dir
    outfile  = NULL,
    width    = 9,
    height   = 6,
    dpi      = 600,
    y_log    = FALSE,
    tstart   = 3000,
    tEnd     = NULL,        # upper x-limit in YEARS (from tstart)
    Neq      = 5000,        # kept for backward compatibility; not used for y scaling
    alpha_many = 0.40,      # alpha for many-run lines
    lw_many    = 0.3,       # linewidth for many-run lines
    lw_mean    = 0.8        # linewidth for mean line
) {
  if (is.null(run_ids)) run_ids <- list_run_ids(base_dir)
  
  # Load & stack all runs
  all_series <- build_series_long_multi(base_dir, run_ids, patch_id, tstart)
  
  # Compute mean trajectory at each timepoint per group
  mean_series <- all_series %>%
    dplyr::group_by(group, Years) %>%
    dplyr::summarize(mean_count = mean(count, na.rm = TRUE), .groups = "drop")
  
  # Publication-friendly theme (BW, larger fonts)
  base_theme <- theme_bw(base_size = 18) +
    theme(
      legend.position      = c(0.5, 0.98),
      legend.justification = c(1, 1),
      legend.direction     = "vertical",
      legend.background    = element_rect(fill = alpha("white", 0.7), colour = NA),
      legend.title         = element_blank(),
      legend.text          = element_text(size = 16),
      plot.title           = element_text(face = "plain", hjust = 0, size = 18),
      axis.title           = element_text(face = "plain", size = 18),
      axis.text            = element_text(size = 16),
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_line(linewidth = 0.25)
    )
  
  # Palette
  okabe_ito <- c(
    "female total" = "#7F7F7F",  # grey #E69F00
    "male Y"       = "#1F77B4",  # blue 
    "male y"       = "#D55E00"   # orange-red #0072B2
  )
  
  # Build plot
  p <- ggplot() +
    # Many runs: thin translucent lines
    geom_line(
      data = all_series,
      aes(x = Years, y = count, color = group, group = interaction(run_id, group)),
      linewidth = lw_many, alpha = alpha_many
    ) +
    # Mean trajectory: thicker
    geom_line(
      data = mean_series,
      aes(x = Years, y = mean_count, color = group),
      linewidth = lw_mean, alpha = 0.9, linetype = "solid"
    ) +
    labs(
      title = sprintf("Patch %s — %d runs (tstart = %s)", patch_id, length(run_ids), tstart),
      x = "Time (years)",
      y = if (y_log) "Abundance (log10)" else "Abundance"
    ) +
    scale_color_manual(values = okabe_ito) +
    base_theme
  
  # X scale (optional right limit)
  x_breaks <- if (is.null(tEnd)) pretty(all_series$Years) else pretty(c(0, tEnd))
  p <- p + scale_x_continuous(
    breaks = x_breaks,
    limits = if (is.null(tEnd)) NULL else c(0, tEnd),
    expand = expansion(mult = c(0.01, 0.02))
  )
  
  # Y scale: auto ticks & auto upper limit (no manual limits)
  if (y_log) {
    p <- p + scale_y_log10(
      labels = scales::label_number(big.mark = ""),
      breaks = scales::breaks_log(n = 6, base = 10)
    )
    # optional: clamp lower bound for readability if needed
    # p <- p + coord_cartesian(ylim = c(1, NA))
  } else {
    p <- p + scale_y_continuous(
      labels = scales::label_number(big.mark = ""),
      breaks = scales::breaks_pretty(n = 6)
    ) +
      coord_cartesian(ylim = c(0, NA))  # start at 0, auto upper
  }
  
  # Output path
  if (is.null(outfile)) {
    outfile <- file.path(
      path.expand(base_dir),
      sprintf("Figure_Patch%s_MULTI_t%s.png", patch_id, tstart)
    )
  }
  
  ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  message(sprintf("Saved figure: %s", outfile))
  
  invisible(list(
    plot = p,
    outfile = outfile,
    runs = run_ids,
    series_all = all_series,
    series_mean = mean_series
  ))
}


# basedir <- "~/Library/CloudStorage/Box-Box/Research Projects/Rodent GD/MGDrivE-master/ExampleCode/mgdriveYLE_sweep_10yr_release1/rel0p08_fy0p80_pq0p00_fl1_mu1p00_j0p00_c1p00"
# basedir <- "~/Library/CloudStorage/Box-Box/Research Projects/Rodent GD/MGDrivE-master/ExampleCode/mgdriveYLE_sweep_10yr_release_recessive/rel0p10_fy0p80_pq0p00_fl1_mu1p00_j0p00_c1p00"


# Auto-detect all runs in base_dir
# make_patch_figure_multi(
#   base_dir = basedir,
#   patch_id = "001",
#   tstart   = round(8*365),
#   tEnd     = 50,         # show first 6 years after tstart
#   y_log    = FALSE,
#   width    = 10,
#   height   = 6,
# )

# # Or specify a subset of runs explicitly
# make_patch_figure_multi(
#   base_dir = ".../mgdriveYLE_2patch",
#   patch_id = "001",
#   run_ids  = c("001","002","003","004","005"),
#   tstart   = 3000,
#   tEnd     = NULL
# )
