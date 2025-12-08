# ---- gene_drive_plots_multi_more.R (updated) ----
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
  
  # Add Years and Total population
  joined <- joined %>%
    mutate(
      Years             = Time_rel / 365,
      `Total population` = `WT female` + `TG female` + `WT male` + `YLE male`
    )
  
  # Long format with 5 groups (4 original + Total population)
  series_long <- joined %>%
    dplyr::select(
      Years,
      `WT female`, `TG female`,
      `WT male`,   `YLE male`,
      `Total population`
    ) %>%
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
  "YLE male"  = "#0072B2",
  "Total population" = "#CC79A7"   # or another contrasting colour
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

# ---------- color helpers ----------
lighten_hex <- function(col, amt = 0.88) {
  rgb(t(col2rgb(col) * (1 - amt) + 255 * amt) / 255)
}

# Build a per-scenario palette by lightening the base palette
make_scenario_palette <- function(base_palette, amt) {
  vapply(base_palette, lighten_hex, character(1), amt = amt)
}

# Validate scenario args
.validate_scenarios <- function(base_dirs, labels, lighten_amts) {
  stopifnot(length(base_dirs) == length(labels))
  if (is.null(lighten_amts)) {
    # first stays base (0), others progressively lighter
    lighten_amts <- c(0, rep(0.55, length(base_dirs) - 1))
  }
  if (length(lighten_amts) == 1L) lighten_amts <- rep(lighten_amts, length(base_dirs))
  stopifnot(length(lighten_amts) == length(base_dirs))
  lighten_amts
}

# ---------- multi-scenario summary builder ----------
# Returns:
#   list(data_all, ci_all)
# where:
#   - data_all: long trajectories for Total population across scenarios
#               (Years, scenario, run_id, count)
#   - ci_all:   center/CI across runs (Years, scenario, center, lower, upper)
build_multi_scenario_ci <- function(
    base_dirs,
    labels,
    patch_id = "001",
    run_ids = NULL,
    tstart = round(8*365),
    center = c("mean","median")
) {
  center <- match.arg(center)
  
  all_series <- list()
  all_ci     <- list()
  
  for (i in seq_along(base_dirs)) {
    bd  <- base_dirs[i]
    lab <- labels[i]
    
    rids <- if (is.null(run_ids)) list_run_ids(bd) else run_ids
    
    ser_i <- build_series_long_multi(bd, rids, patch_id, tstart) |>
      dplyr::filter(group == "Total population") |>
      dplyr::mutate(scenario = lab)
    
    ci_i <- ser_i |>
      dplyr::group_by(scenario, Years) |>
      dplyr::summarise(
        center = if (center == "mean") mean(count, na.rm = TRUE) else median(count, na.rm = TRUE),
        lower  = stats::quantile(count, 0.025, na.rm = TRUE, type = 7),
        upper  = stats::quantile(count, 0.975, na.rm = TRUE, type = 7),
        .groups = "drop"
      )
    
    all_series[[i]] <- ser_i
    all_ci[[i]]     <- ci_i
  }
  
  data_all <- dplyr::bind_rows(all_series) |>
    dplyr::mutate(scenario = factor(scenario, levels = labels))
  ci_all   <- dplyr::bind_rows(all_ci) |>
    dplyr::mutate(scenario = factor(scenario, levels = labels))
  
  list(
    data_all = data_all,
    ci_all   = ci_all
  )
}


# -----------------------------#
# Plot (center line + 95% CI)
# -----------------------------#
# ---- replace your make_patch_figure_multi with this version ----
# -----------------------------#
# Plot Total population (center line + 95% CI)
# -----------------------------#
make_patch_figure_multi_compare <- function(
    base_dirs,
    labels,
    patch_id = "001",
    run_ids = NULL,
    tstart = round(8*365),
    center = c("mean","median"),
    scenario_cols = NULL,        # user-defined colour vector, one per scenario
    outfile  = NULL,
    width    = 10,
    height   = 6,
    dpi      = 600,
    y_log    = FALSE,
    tEnd     = NULL,
    alpha_ci = 0.15,
    lw_mean  = 1.8,
    show_runs  = FALSE,
    alpha_runs = 0.10,
    lw_runs    = 0.25
) {
  center <- match.arg(center)
  
  built <- build_multi_scenario_ci(
    base_dirs = base_dirs,
    labels    = labels,
    patch_id  = patch_id,
    run_ids   = run_ids,
    tstart    = tstart,
    center    = center
  )
  
  ser   <- built$data_all
  sumci <- built$ci_all
  
  # ----- colour handling: one colour per scenario -----
  if (is.null(scenario_cols)) {
    # default colour-blind friendly vector for three scenarios
    scenario_cols <- c("#E69F00", "#009E73", "#0072B2")
    # scenario_cols <- c("#edf8b1","#7fcdbb","#2c7fb8")
  }
  if (length(scenario_cols) != length(labels)) {
    stop("Length of 'scenario_cols' must equal length of 'labels'.")
  }
  if (is.null(names(scenario_cols))) {
    names(scenario_cols) <- labels
  }
  # ensure order follows labels
  scenario_cols <- scenario_cols[labels]
  
  base_theme <- theme_bw(base_size = 24) +
    theme(
      legend.position      = c(0.35, 1.00),
      legend.justification = c(0, 1),
      legend.direction     = "vertical",
      legend.background    = element_rect(fill = alpha("white", 0.2), colour = NA),
      legend.title         = element_blank(),
      legend.text          = element_text(size = 24),
      axis.title           = element_text(size = 28),
      axis.text            = element_text(size = 24),
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_line(linewidth = 0.25)
    )
  
  p <- ggplot()
  
  # optional: thin per-run lines (still Total population only)
  if (show_runs) {
    p <- p +
      geom_line(
        data = ser,
        aes(x = Years, y = count, color = scenario,
            group = interaction(run_id, scenario)),
        linewidth = lw_runs, alpha = alpha_runs
      )
  }
  
  p <- p +
    geom_ribbon(
      data = sumci,
      aes(x = Years, ymin = lower, ymax = upper, fill = scenario),
      alpha = alpha_ci, colour = NA
    ) +
    geom_line(
      data = sumci,
      aes(x = Years, y = center, color = scenario),
      linewidth = lw_mean
    ) +
    labs(
      x = "Time (years)",
      y = if (y_log) "Total Population (log10)" else "Total Population"
    ) +
    scale_color_manual(values = scenario_cols) +
    scale_fill_manual(values  = scenario_cols) +
    base_theme
  
  # X scale / limits
  x_breaks <- if (is.null(tEnd)) pretty(ser$Years) else pretty(c(0, tEnd))
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
    ) + coord_cartesian(ylim = c(0, NA))
  }
  
  if (is.null(outfile)) {
    outfile <- file.path(
      normalizePath(getwd()),
      sprintf(
        "FIG_compare_TotalPop_Patch%s_t%s_%s_%s.png",
        patch_id, tstart, center,
        paste(gsub("[^A-Za-z0-9]+", "", labels), collapse = "-")
      )
    )
  }
  ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", outfile)
  
  invisible(list(plot = p, outfile = outfile, ci = sumci, scenario_cols = scenario_cols))
}



make_totals_figure_multi_compare <- function(
    base_dirs,
    labels,
    patch_id = "001",
    run_ids = NULL,
    tstart = round(8*365),
    center = c("mean","median"),
    totals_to_plot = c("both"),   # any of: "males","females","both"
    lighten_amts = NULL,          # per-scenario lightening of gray
    outfile  = NULL,
    width    = 9,
    height   = 5.5,
    dpi      = 600,
    y_log    = FALSE,
    tEnd     = NULL,
    alpha_ci = 0.22,
    lw_mean  = 1.05
) {
  center <- match.arg(center)
  totals_to_plot <- intersect(tolower(totals_to_plot), c("males","females","both"))
  if (!length(totals_to_plot)) stop("totals_to_plot must include 'males', 'females', and/or 'both'.")
  
  lighten_amts <- .validate_scenarios(base_dirs, labels, lighten_amts)
  
  # base greyscale for scenarios; then lighten per scenario
  base_gray <- "#000000"
  scen_cols <- vapply(lighten_amts, function(a) lighten_hex(base_gray, amt = a), character(1))
  names(scen_cols) <- labels
  
  all <- list(); all_ci <- list()
  for (i in seq_along(base_dirs)) {
    bd <- base_dirs[i]; lab <- labels[i]
    rids <- if (is.null(run_ids)) list_run_ids(bd) else run_ids
    
    ser_i <- build_series_long_multi(bd, rids, patch_id, tstart) |>
      dplyr::mutate(scenario = lab)
    
    # get totals from long groups
    totals_i <- ser_i |>
      dplyr::mutate(sex =
                      dplyr::case_when(
                        group %in% c("WT male","YLE male") ~ "males",
                        group %in% c("WT female","TG female") ~ "females",
                        TRUE ~ NA_character_
                      )
      ) |>
      dplyr::filter(!is.na(sex)) |>
      dplyr::group_by(scenario, Years, sex, run_id) |>
      dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
    
    # optionally add "both" = males+females per run
    if ("both" %in% totals_to_plot) {
      both_i <- totals_i |>
        dplyr::group_by(scenario, Years, run_id) |>
        dplyr::summarise(count = sum(count), .groups = "drop") |>
        dplyr::mutate(sex = "both")
      totals_i <- dplyr::bind_rows(totals_i, both_i)
    }
    
    totals_i <- dplyr::filter(totals_i, sex %in% totals_to_plot)
    all[[i]] <- totals_i
    
    ci_i <- totals_i |>
      dplyr::group_by(scenario, Years, sex) |>
      dplyr::summarise(
        center = if (center == "mean") mean(count, na.rm = TRUE) else median(count, na.rm = TRUE),
        lower  = stats::quantile(count, 0.025, na.rm = TRUE, type = 7),
        upper  = stats::quantile(count, 0.975, na.rm = TRUE, type = 7),
        .groups = "drop"
      )
    all_ci[[i]] <- ci_i
  }
  
  totals_all <- dplyr::bind_rows(all)
  ci_all     <- dplyr::bind_rows(all_ci)
  
  linemap <- c(males = "solid", females = "twodash", both = "longdash")[totals_to_plot]
  
  base_theme <- theme_bw(base_size = 24) +
    theme(
      legend.position      = c(0.35, 0.98),
      legend.title         = element_blank(),
      legend.text          = element_text(size = 24),
      axis.title           = element_text(size = 28),
      axis.text            = element_text(size = 24),
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_line(linewidth = 0.25)
    )
  
  p <- ggplot() +
    geom_ribbon(
      data = ci_all,
      aes(x = Years, ymin = lower, ymax = upper, fill = scenario, linetype = sex),
      alpha = alpha_ci, colour = NA, show.legend = FALSE
    ) +
    geom_line(
      data = ci_all,
      aes(x = Years, y = center, color = scenario, linetype = sex),
      linewidth = lw_mean
    ) +
    scale_color_manual(values = scen_cols) +
    scale_fill_manual(values = scen_cols) +
    scale_linetype_manual(values = linemap) +
    labs(x = "Time (years)", y = if (y_log) "Abundance (log10)" else "Abundance") +
    base_theme
  
  x_breaks <- if (is.null(tEnd)) pretty(totals_all$Years) else pretty(c(0, tEnd))
  p <- p + scale_x_continuous(
    breaks = x_breaks,
    limits = if (is.null(tEnd)) NULL else c(0, tEnd),
    expand = expansion(mult = c(0.01, 0.02))
  )
  
  if (y_log) {
    p <- p + scale_y_log10(
      labels = scales::label_number(big.mark = ""),
      breaks = scales::breaks_log(n = 6, base = 10)
    )
  } else {
    p <- p + scale_y_continuous(
      labels = scales::label_number(big.mark = ""),
      breaks = scales::breaks_pretty(n = 6)
    ) + coord_cartesian(ylim = c(0, NA))
  }
  
  if (is.null(outfile)) {
    outfile <- file.path(
      normalizePath(getwd()),
      sprintf("FIG_totals_compare_Patch%s_t%s_%s.png",
              patch_id, tstart, paste(gsub("[^A-Za-z0-9]+","",labels), collapse = "-"))
    )
  }
  ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved totals: ", outfile)
  
  invisible(list(plot = p, outfile = outfile, ci = ci_all))
}

base_dirs <- c(
  "rel0p01_fy1p00_pq0p00_fl0p00_fs1p00_mu0p97_j0p25_c0p93",
  "rel0p02_fy1p00_pq0p00_fl0p00_fs1p00_mu0p97_j0p25_c0p93",
  "rel0p08_fy1p00_pq0p00_fl0p00_fs1p00_mu0p97_j0p25_c0p93"
)

labels <- c("50 YLE/month", "100 YLE/month", "400 YLE/month")

# Okabe–Ito colour-blind-friendly, high contrast
# scenario_cols <- c(
#   "#E69F00",  # orange
#   "#009E73",  # bluish green
#   "#0072B2"   # blue
# )


scenario_cols <- c("#006A8EFF", "#B1283AFF", "#A8A6A7FF") # c("#BB9753FF","#061922FF","#008348FF") # c("#252525","#CC79A7","#009E73") # Green: "#009E73" # "#525252" "#252525"

# c("#66c2a5","#fc8d62","#8da0cb")
names(scenario_cols) <- labels

make_patch_figure_multi_compare(
  base_dirs     = base_dirs,
  labels        = labels,
  patch_id      = "001",
  tstart        = round(8*365),
  center        = "mean",
  tEnd          = 20,
  scenario_cols = scenario_cols
)


# make_patch_figure_multi_compare(
#   base_dirs = base_dirs,
#   labels    = labels,
#   patch_id  = "001",
#   tstart    = round(8*365),
#   center    = "mean",
#   tEnd     = 20,
#   groups_to_plot = c("WT female"), #c("WT female","WT male","YLE male","TG female"),
#   lighten_amts   = col1 # base, lighter, lightest
# )

# palette4 <- c(
#   "WT female" = "#E69F00",
#   "TG female" = "#F4C97A",
#   "WT male"   = "#009E73",
#   "YLE male"  = "#0072B2"
# )

# make_totals_figure_multi_compare(
#   base_dirs = base_dirs,
#   labels    = labels,
#   patch_id  = "001",
#   tstart    = round(8*365),
#   center    = "mean",
#   totals_to_plot = c("both"),
#   tEnd     = 20,
#   lighten_amts   = col1 # base, lighter, lightest
# )

