# ---- gene_drive_plots_multi.R (updated) ----
# Usage:
#   source("gene_drive_plots_multi.R")
#   make_patch_figure_multi(
#     base_dir = "/path/to/rel.../your_scenario",
#     patch_id = "001",
#     tstart   = round(8*365),
#     center   = "median"  # or "mean"
#   )

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# -----------------------------#
# Discover runs in a directory
# -----------------------------#
list_run_ids <- function(base_dir) {
  base_dir <- path.expand(base_dir)
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  run_ids <- unique(c(
    sub("^Run", "", dirs[grepl("^Run\\d+$", dirs)]),
    dirs[grepl("^\\d+$", dirs)]
  ))
  if (!length(run_ids)) stop(sprintf("No run folders found in %s", base_dir))
  sort(run_ids)
}

#-----------------------------#
# I/O + path helpers
#-----------------------------#
read_ts_csv <- function(path) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"Time" %in% names(df)) stop(sprintf("'Time' column not found in %s", path))
  df %>% mutate(Time = suppressWarnings(as.numeric(Time))) %>% arrange(Time) -> df
  if (anyNA(df$Time)) stop(sprintf("Non-numeric values in 'Time' of %s", path))
  df
}

resolve_run_dir <- function(base_dir, run_id) {
  base_dir <- path.expand(base_dir)
  candidates <- file.path(base_dir, c(run_id, sprintf("Run%s", run_id)))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop(sprintf("Cannot find run folder for %s", run_id))
  hit[1]
}

#-----------------------------#
# Category summaries
#-----------------------------#
summarize_males_split <- function(df) {
  male_Y_cols <- grep("^mY", names(df), value = TRUE)  # WT males
  male_y_cols <- grep("^my", names(df), value = TRUE)  # YLE males
  if (!length(male_Y_cols) && !length(male_y_cols)) {
    stop("No male genotype columns starting with 'mY' or 'my' found.")
  }
  tibble(
    Time    = df$Time,
    `WT male` = if (length(male_Y_cols)) rowSums(df[male_Y_cols]) else 0,
    `YLE male`  = if (length(male_y_cols)) rowSums(df[male_y_cols]) else 0
  )
}

summarize_females_split <- function(df_like_male_csv) {
  # Try to locate a separate F_Aggregate* first; if not, we’ll rely on fX* columns
  female_cols_all <- grep("^fX", names(df_like_male_csv), value = TRUE)
  
  # If columns indicate genotype detail (e.g., fXWW, fXHW, fXHH), split by WT vs TG.
  # WT female: fXWW
  # TG female: any H-bearing (fXHW, fXWH, fXHH). We include WH just in case column order varies.
  has_genotype_detail <- any(grepl("^fX[WH]{2}$", female_cols_all))
  
  if (!length(female_cols_all)) {
    stop("No female genotype columns starting with 'fX' found.")
  }
  
  if (has_genotype_detail) {
    f_wt <- grep("^fXWW$", names(df_like_male_csv), value = TRUE)
    f_tg <- grep("^fX(HW|WH|HH)$", names(df_like_male_csv), value = TRUE)
    tibble(
      Time       = df_like_male_csv$Time,
      `WT female` = if (length(f_wt)) rowSums(df_like_male_csv[f_wt]) else 0,
      `TG female` = if (length(f_tg)) rowSums(df_like_male_csv[f_tg]) else 0
    )
  } else {
    # Fallback: if only totals exist, treat all as WT for compatibility (or all TG = 0)
    warning("Female genotype detail not found (fXWW / fXHW / fXHH). Assuming all fX* are WT females.")
    tibble(
      Time        = df_like_male_csv$Time,
      `WT female` = rowSums(df_like_male_csv[female_cols_all]),
      `TG female` = 0
    )
  }
}

get_female_totals_split <- function(base_dir, run_id, patch_id, male_df) {
  run_dir <- resolve_run_dir(base_dir, run_id)
  fpath   <- file.path(run_dir, sprintf("F_Aggregate_Run%s_Patch%s.csv", run_id, patch_id))
  if (file.exists(fpath)) {
    # When an F_Aggregate is present, it might be totals only; try to split if columns allow
    fdf <- read_ts_csv(fpath)
    summarize_females_split(fdf)
  } else {
    # Fall back to the male file, which often also includes female genotype columns
    summarize_females_split(male_df)
  }
}

#-----------------------------#
# Time filtering / long build
#-----------------------------#
filter_and_reindex_time <- function(df, tstart) {
  out <- df %>% filter(Time >= tstart) %>% mutate(Time_rel = Time - tstart)
  if (!nrow(out)) stop(sprintf("No rows at or after tstart = %s", tstart))
  out
}

build_series_long <- function(base_dir, run_id, patch_id, tstart) {
  run_dir   <- resolve_run_dir(base_dir, run_id)
  male_path <- file.path(run_dir, sprintf("M_Run%s_Patch%s.csv", run_id, patch_id))
  
  male_raw  <- read_ts_csv(male_path)
  males     <- summarize_males_split(male_raw)
  females   <- get_female_totals_split(base_dir, run_id, patch_id, male_raw)
  
  joined <- males %>%
    left_join(females, by = "Time") %>%
    filter_and_reindex_time(tstart = tstart)
  
  # Long format with four groups
  series_long <- joined %>%
    mutate(Years = Time_rel / 365) %>%
    dplyr::select(Years, `WT female`, `TG female`, `WT male`, `YLE male`) %>%
    pivot_longer(-Years, names_to = "group", values_to = "count")
  
  list(
    run_dir     = run_dir,
    male_path   = male_path,
    data_joined = joined,
    series_long = series_long
  )
}

build_series_long_multi <- function(base_dir, run_ids, patch_id, tstart) {
  purrr::map_dfr(run_ids, function(rid) {
    built <- build_series_long(base_dir, rid, patch_id, tstart)
    built$series_long %>% mutate(run_id = rid)
  })
}

#-----------------------------#
# CI summaries
#-----------------------------#
summarize_center_ci <- function(series_long, center = c("mean","median"), ci = 0.95) {
  center <- match.arg(center)
  a <- (1 - ci) / 2
  b <- 1 - a
  series_long %>%
    group_by(group, Years) %>%
    summarise(
      center = if (center == "mean") mean(count, na.rm = TRUE) else median(count, na.rm = TRUE),
      lower  = quantile(count, probs = a, na.rm = TRUE, type = 7),
      upper  = quantile(count, probs = b, na.rm = TRUE, type = 7),
      .groups = "drop"
    )
}

# ---- add/replace: palette and a small validator ----
palette4 <- c(
  "WT female" = "#E69F00",
  "TG female" = "#F4C97A",
  "WT male"   = "#009E73",
  "YLE male"  = "#0072B2"
)


# palette4 <- c(
#   "WT female" = "#CC79A7",
#   "TG female" = "#F4C97A",
#   "WT male"   = "#009E73",
#   "YLE male"  = "#56B4E9"
# )

.validate_groups <- function(requested, valid_names) {
  if (is.null(requested)) return(valid_names)
  missing <- setdiff(requested, valid_names)
  if (length(missing)) {
    stop(sprintf("Unknown group(s): %s\nValid groups are: %s",
                 paste(missing, collapse = ", "),
                 paste(valid_names, collapse = ", ")))
  }
  unique(requested)
}


# -----------------------------#
# Plot (center line + 95% CI)
# -----------------------------#
# ---- replace your make_patch_figure_multi with this version ----
make_patch_figure_multi <- function(
    base_dir,
    patch_id = "001",
    run_ids  = NULL,
    outfile  = NULL,
    width    = 10,
    height   = 6,
    dpi      = 600,
    y_log    = FALSE,
    tstart   = round(8*365),
    tEnd     = NULL,            # x-limit in YEARS (from tstart)
    alpha_ci = 0.25,            # transparency for CI ribbons
    lw_mean  = 1.8,             # center line width
    center   = c("mean","median"),
    show_runs   = FALSE,        # overlay thin per-run lines
    alpha_runs  = 0.15,
    lw_runs     = 0.25,
    groups_to_plot = NULL       # <-- NEW: e.g., c("WT female","YLE male")
) {
  if (is.null(run_ids)) run_ids <- list_run_ids(base_dir)
  
  all_series <- build_series_long_multi(base_dir, run_ids, patch_id, tstart)
  sum_ci     <- summarize_center_ci(all_series, center = center, ci = 0.95)
  
  # Determine which groups are available and subset
  valid_groups <- sort(unique(sum_ci$group))
  groups_use   <- .validate_groups(groups_to_plot, valid_groups)
  
  all_series <- dplyr::filter(all_series, group %in% groups_use)
  sum_ci     <- dplyr::filter(sum_ci,     group %in% groups_use)
  
  # Keep color/fill only for selected groups (order respected in legend)
  pal_use <- palette4[groups_use]
  
  base_theme <- theme_bw(base_size = 24) +
    theme(
      legend.position      = c(0.70, 0.98),
      legend.justification = c(0, 1),
      legend.direction     = "vertical",
      legend.background    = element_rect(fill = alpha("white", 0.7), colour = NA),
      legend.title         = element_blank(),
      legend.text          = element_text(size = 24),
      plot.title           = element_text(hjust = 0, size = 20),
      axis.title           = element_text(size = 28),
      axis.text            = element_text(size = 24),
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_line(linewidth = 0.25)
    )
  
  p <- ggplot()
  
  if (show_runs) {
    p <- p +
      geom_line(
        data = all_series,
        aes(x = Years, y = count, color = group, group = interaction(run_id, group)),
        linewidth = lw_runs, alpha = alpha_runs
      )
  }
  
  p <- p +
    geom_ribbon(
      data = sum_ci,
      aes(x = Years, ymin = lower, ymax = upper, fill = group),
      alpha = alpha_ci, colour = NA
    ) +
    geom_line(
      data = sum_ci,
      aes(x = Years, y = center, color = group),
      linewidth = lw_mean
    ) +
    labs(
      # title = sprintf("Patch %s — %d runs (tstart = %s days)", patch_id, length(run_ids), tstart),
      x = "Time (years)",
      y = if (y_log) "Abundance (log10)" else "Abundance"
    ) +
    scale_color_manual(values = pal_use, limits = groups_use, breaks = groups_use) +
    scale_fill_manual(values = pal_use,  limits = groups_use, breaks = groups_use) +
    base_theme
  
  # X scale / limits
  x_breaks <- if (is.null(tEnd)) pretty(all_series$Years) else pretty(c(0, tEnd))
  p <- p + scale_x_continuous(
    breaks = x_breaks,
    limits = if (is.null(tEnd)) NULL else c(0, tEnd),
    expand = expansion(mult = c(0.01, 0.02))
  )
  
  # Y scale
  if (y_log) {
    p <- p + scale_y_log10(
      labels = scales::label_number(big.mark = ""),
      breaks = scales::breaks_log(n = 6, base = 10)
    )
  } else {
    p <- p + scale_y_continuous(
      labels = scales::label_number(big.mark = ""),
      breaks = scales::breaks_pretty(n = 6)
    ) +
      coord_cartesian(ylim = c(0, NA))
  }
  
  if (is.null(outfile)) {
    outfile <- file.path(
      path.expand(base_dir),
      sprintf(
        "Figure_Patch%s_MULTI_t%s_%s_CI_%s.png",
        patch_id, tstart, match.arg(center),
        paste(gsub(" ", "", groups_use), collapse = "-")
      )
    )
  }
  
  ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  message(sprintf("Saved figure: %s", outfile))
  
  invisible(list(
    plot = p,
    outfile = outfile,
    runs = run_ids,
    series_all = all_series,
    summary_ci = sum_ci,
    groups = groups_use
  ))
}


# 
# outFolder <- "mgdriveYLE_nonideal_viablity_10runs_5yr_400rel"
# base_dir <- file.path(getwd(), outFolder)
# # source("gene_drive_plots_multi_new.R")
# make_patch_figure_multi(
#   base_dir = base_dir,
#   patch_id = "001",
#   tstart   = round(8*365),
#   tEnd     = 20,
#   center   = "mean",  # or "mean"
#   groups_to_plot = c("WT female","YLE male")
# )
