# ---- gene_drive_plots.R ----
# Usage:
#   source("gene_drive_plots.R")
#   make_run_patch_figure_single(run_id = "001", patch_id = "001", tstart = 3000)
# or run via Rscript after adding a call at the bottom.

suppressPackageStartupMessages({
  library(tidyverse)     # ggplot2, dplyr, readr, tidyr
  library(patchwork)
  library(RColorBrewer)
  library(scales)
})

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
    select(Time, male_Y, male_y, male_total, freq_Y, freq_y)
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

#-----------------------------#
# Plot (single figure)
#-----------------------------#

make_run_patch_figure_single <- function(
    base_dir = "~/Library/CloudStorage/Box-Box/Research Projects/Rodent GD/MGDrivE-master/ExampleCode/mgdriveYLE_2patch",
    run_id   = "001",
    patch_id = "001",
    outfile  = NULL,
    width    = 9,
    height   = 6,
    dpi      = 300,
    y_log    = FALSE,      # set TRUE to use log10 y-axis
    tstart   = 3000,       # plot starting generation
    Neq      = 5000,       # suggested y-limit helper (linear scale)
    alpha_ln = 0.8,         # line transparency (0..1)
    tEnd     = NULL        #  upper limit for x-axis (in YEARS)
) {
  built <- build_series_long(base_dir, run_id, patch_id, tstart)
  run_dir     <- built$run_dir
  series_long <- built$series_long
  
  # Add years column when building long data
  series_long <- built$series_long %>%
    mutate(Years = Time_rel / 365)
  
  # Publication-friendly theme (BW, larger fonts)
  base_theme <- theme_bw(base_size = 18) +
    theme(
      legend.position    = c(0.98, 0.98),
      legend.justification= c(1, 1),
      legend.direction   = "vertical",
      legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
      legend.title       = element_blank(),
      legend.text        = element_text(size = 16),
      plot.title         = element_text(face = "plain", hjust = 0, size = 18),
      axis.title         = element_text(face = "plain", size = 18),
      axis.text          = element_text(size = 16),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(linewidth = 0.25)
      # panel.border       = element_rect(linewidth = 0.8, color = "black")
    )
  
  # Color palette: female = grey, Y males = blue, y males = orange/red
  okabe_ito <- c(
    "female total" = "#7F7F7F",  # medium grey
    "male Y"       = "#1F77B4",  # clear blue
    "male y"       = "#D55E00"   # reddish orange
  )
  
  # Plot with ColorBrewer colors + transparency; solid lines (journal-friendly)
  p <- ggplot(series_long, aes(Years, count, color = group)) +
    geom_line(linewidth = 1.1, alpha = alpha_ln) +
    labs(
      x = "Time (years)",
      y = if (y_log) "Abundance (log10)" else "Abundance"
    ) +
    # scale_color_brewer(palette = "Set1") +
    scale_color_manual(values = okabe_ito) +   # <- custom colorblind-safe palette
    # scale_color_viridis_d(option = "plasma") +
    base_theme
  
  # # X scale: break in whole years 
  # p <- p + scale_x_continuous(
  #   breaks = pretty(series_long$Years),  # nice whole-number breaks
  #   expand = expansion(mult = c(0.01, 0.02))
  # )  
  
  # ---------- X scale ----------
  # If tEnd is given, clamp x-axis to [0, tEnd] years; otherwise use data range.
  x_breaks <- if (is.null(tEnd)) pretty(series_long$Years) else pretty(c(0, tEnd))
  p <- p + scale_x_continuous(
    breaks = x_breaks,
    limits = if (is.null(tEnd)) NULL else c(0, tEnd),
    expand = expansion(mult = c(0.01, 0.02))
  )
  
  # Y scale: **no commas** in labels
  # Y scale: no commas in labels, fixed breaks every 1000
  if (y_log) {
    p <- p + scale_y_log10(labels = label_number(big.mark = ""))
  } else {
    p <- p + scale_y_continuous(
      labels = label_number(big.mark = ""),
      breaks = seq(0, Neq, by = round(Neq/5)),         # major ticks
      minor_breaks = seq(0, Neq, by = round(Neq/5)),    # minor ticks in between
      limits = c(0, Neq+round(0.5*Neq/5))                       # clip at Neq
    )
  }
  
  # Output path
  if (is.null(outfile)) {
    outfile <- file.path(run_dir, sprintf("Figure_Run%s_Patch%s_single_t%s.png", run_id, patch_id, tstart))
  }
  
  ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  message(sprintf("Saved figure: %s", outfile))
  
  invisible(list(plot = p, outfile = outfile, series = series_long))
}

#-----------------------------#
# Example (comment/uncomment)
#-----------------------------#
# make_run_patch_figure_single(run_id = "001", patch_id = "001", tstart = 3000, y_log = FALSE, alpha_ln = 0.8)
