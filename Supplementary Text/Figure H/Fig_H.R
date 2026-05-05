suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

list_run_ids <- function(base_dir) {
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  run_ids <- unique(c(
    sub("^Run", "", dirs[grepl("^Run\\d+$", dirs)]),
    dirs[grepl("^\\d+$", dirs)]
  ))
  if (!length(run_ids)) stop(sprintf("No run folders found in %s", base_dir))
  sort(run_ids)
}

resolve_run_dir <- function(base_dir, run_id) {
  candidates <- file.path(base_dir, c(run_id, sprintf("Run%s", run_id)))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop(sprintf("Cannot find run folder for run %s", run_id))
  hit[1]
}

read_ts_csv <- function(path) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(Time = as.numeric(Time)) %>%
    arrange(Time)
}

extract_rel_label <- function(folder_name) {
  rel <- stringr::str_extract(folder_name, "rel[0-9]+p[0-9]+")
  rel <- stringr::str_replace(rel, "rel", "")
  rel <- stringr::str_replace(rel, "p", ".")
  paste0(as.numeric(rel) * 100, "% release")
}

get_YLE_freq_one_run <- function(base_dir, run_id, patch_id = "001", tstart = round(8 * 365)) {
  run_dir <- resolve_run_dir(base_dir, run_id)
  path <- file.path(run_dir, sprintf("M_Run%s_Patch%s.csv", run_id, patch_id))
  
  df <- read_ts_csv(path)
  
  yle_cols <- grep("^my", names(df), value = TRUE)
  wt_male_cols <- grep("^mY", names(df), value = TRUE)
  
  if (!length(yle_cols)) stop("No YLE male genotype columns starting with 'my' found.")
  if (!length(wt_male_cols)) stop("No WT male genotype columns starting with 'mY' found.")
  
  df %>%
    transmute(
      Time = Time,
      Years = (Time - tstart) / 365,
      YLE_male = rowSums(across(all_of(yle_cols))),
      WT_male  = rowSums(across(all_of(wt_male_cols))),
      total_male = YLE_male + WT_male,
      YLE_freq = if_else(total_male > 0, YLE_male / total_male, NA_real_)
    ) %>%
    filter(Time >= tstart)
}

build_YLE_freq_multi_folder <- function(
    folders,
    parent_dir = getwd(),
    patch_id = "001",
    run_ids = NULL,
    tstart = round(8 * 365)
) {
  purrr::map_dfr(folders, function(folder) {
    base_dir <- file.path(parent_dir, folder)
    ids <- if (is.null(run_ids)) list_run_ids(base_dir) else run_ids
    
    purrr::map_dfr(ids, function(rid) {
      get_YLE_freq_one_run(base_dir, rid, patch_id, tstart) %>%
        mutate(
          run_id = rid,
          scenario = folder,
          rel_label = extract_rel_label(folder)
        )
    })
  })
}

summarise_YLE_freq <- function(df, center = c("mean", "median"), ci = 0.95) {
  center <- match.arg(center)
  a <- (1 - ci) / 2
  b <- 1 - a
  
  df %>%
    group_by(scenario, rel_label, Years) %>%
    summarise(
      center = if (center == "mean") mean(YLE_freq, na.rm = TRUE) else median(YLE_freq, na.rm = TRUE),
      lower  = quantile(YLE_freq, a, na.rm = TRUE),
      upper  = quantile(YLE_freq, b, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_YLE_freq_by_release <- function(
    folders,
    parent_dir = getwd(),
    patch_id = "001",
    run_ids = NULL,
    tstart = round(8 * 365),
    tEnd = NULL,
    center = "mean",
    outfile = NULL,
    width = 10,
    height = 6,
    dpi = 600
) {
  all_runs <- build_YLE_freq_multi_folder(
    folders = folders,
    parent_dir = parent_dir,
    patch_id = patch_id,
    run_ids = run_ids,
    tstart = tstart
  )
  
  summary_df <- summarise_YLE_freq(all_runs, center = center)
  
  # palette_yle <- c(
  #   "2% release" = "#e41a1c",
  #   "4% release" = "#377eb8",
  #   "6% release" = "#4daf4a",
  #   "8% release" = "#984ea3"
  # )
  # 
  # palette_yle <- c(
  #   "2% release" = "#f1eef6",
  #   "4% release" = "#bdc9e1",
  #   "6% release" = "#74a9cf",
  #   "8% release" = "#0570b0"
  # )
  
  palette_yle <- c(
    "2% release" = "#D55E00",  # vermillion
    "4% release" = "#0072B2",  # blue
    "6% release" = "#009E73",  # bluish green
    "8% release" = "#CC79A7"   # reddish purple
  )
  
  #f1eef6
  #bdc9e1
  #74a9cf
  #0570b0
  
  p <- ggplot(summary_df, aes(x = Years, y = center, color = rel_label, fill = rel_label)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.30, color = NA) +
    geom_line(linewidth = 1.8) +
    scale_y_continuous(
      limits = c(0, 1.02),
      breaks = seq(0, 1, 0.1),
      labels = label_number(accuracy = 0.1)
    ) +
    scale_x_continuous(
      limits = if (is.null(tEnd)) NULL else c(0, tEnd),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_color_manual(values = palette_yle) +
    scale_fill_manual(values = palette_yle) +
    labs(
      x = "Time since release started (years)",
      y = "YLE male frequency",
      color = "Release",
      fill = "Release"
    ) +
    theme_bw(base_size = 22) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
      axis.title.y = element_text(margin = margin(r = 12)),
      axis.title.x = element_text(margin = margin(t = 10)),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.key.height = unit(1.2, "lines")
      )
  
  if (is.null(outfile)) {
    outfile <- file.path(parent_dir, "FigH.png")
  }
  
  ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  
  invisible(list(
    plot = p,
    all_runs = all_runs,
    summary = summary_df,
    outfile = outfile
  ))
}

folders <- c(
  "rel0p02_fy1p00_pq0p00_fl0p00_fs1p00_mu0p97_j0p25_c0p93",
  "rel0p04_fy1p00_pq0p00_fl0p00_fs1p00_mu0p97_j0p25_c0p93",
  "rel0p06_fy1p00_pq0p00_fl0p00_fs1p00_mu0p97_j0p25_c0p93",
  "rel0p08_fy1p00_pq0p00_fl0p00_fs1p00_mu0p97_j0p25_c0p93"
)

out <- plot_YLE_freq_by_release(
  folders = folders,
  parent_dir = getwd(),
  patch_id = "001",
  tstart = round(8 * 365),
  tEnd = 20,
  center = "mean"
)

out$plot