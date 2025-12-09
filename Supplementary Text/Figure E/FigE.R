rm(list = ls())
gc()

library(tidyverse)
library(grid)   # for unit()

# -----------------------------#
# Helper: parse K from tag     #
# e.g. "K1k" -> 1000, "K200k"->200000
# -----------------------------#
parse_K_from_tag <- function(tag) {
  tag <- as.character(tag)
  m <- regexec("^K?([0-9]+)([kKmM]?)$", tag)
  r <- regmatches(tag, m)[[1]]
  if (length(r) == 0) {
    stop(sprintf("Cannot parse K from scenario_tag '%s'", tag))
  }
  num <- as.numeric(r[2])
  mag <- r[3]
  if (is.na(num)) stop(sprintf("Numeric part of '%s' is NA", tag))
  if (mag %in% c("k", "K")) {
    num * 1e3
  } else if (mag %in% c("m", "M")) {
    num * 1e6
  } else {
    num
  }
}

# -----------------------------#
# Helper: robust CSV reader    #
# -----------------------------#
read_ts_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"Time" %in% names(df)) {
    stop(sprintf("'Time' column not found in %s", path))
  }
  df
}

# -----------------------------#
# Summarise one scenario       #
# Returns: Time, mean_norm, q2p5_norm,
#          q97p5_norm, label, scenario_tag, K
# -----------------------------#
library(stringr)

# helper: summarise one F or M file into total pop per Time & Run
summarise_sex_file <- function(path, sex_label) {
  dat <- read_ts_csv(path)
  
  if (!"Time" %in% names(dat)) {
    stop(sprintf("'Time' column not found in %s", path))
  }
  
  # extract run id from filename: ...Run001_...
  run_match <- str_match(basename(path), "[Rr]un(\\d+)")
  run_id <- suppressWarnings(as.integer(run_match[, 2]))
  if (is.na(run_id)) {
    run_id <- NA_integer_
  }
  
  num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, "Time")
  if (!length(num_cols)) {
    stop(sprintf("No numeric population columns found in %s", path))
  }
  
  pop <- rowSums(dat[, num_cols, drop = FALSE], na.rm = TRUE)
  
  tibble(
    Time      = dat$Time,
    Run       = run_id,
    Sex       = sex_label,
    Pop_sex   = pop
  )
}

summarise_scenario <- function(base_dir,
                               scenario_tag,
                               label,
                               color = NULL) {
  if (!dir.exists(base_dir)) {
    stop(sprintf("Base directory not found: %s (wd: %s)",
                 base_dir, getwd()))
  }
  
  scen_dir <- file.path(base_dir, scenario_tag)
  if (!dir.exists(scen_dir)) {
    stop(sprintf("Scenario directory not found: %s", scen_dir))
  }
  
  # F and M files: e.g. F_Aggregate_Run001_Patch001.csv, M_Run001_Patch001.csv
  files_F <- list.files(
    scen_dir,
    pattern    = "^F_.*Run\\d+_Patch\\d+\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  )
  files_M <- list.files(
    scen_dir,
    pattern    = "^M_.*Run\\d+_Patch\\d+\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  )
  
  if (!length(files_F) && !length(files_M)) {
    ex <- utils::head(list.files(scen_dir, recursive = TRUE), 30)
    stop(sprintf(
      "No F/M run files found under %s.\nExample files:\n%s",
      scen_dir,
      paste(ex, collapse = "\n")
    ))
  }
  
  df_F <- if (length(files_F)) {
    purrr::map_dfr(files_F, summarise_sex_file, sex_label = "F")
  } else {
    tibble(Time = numeric(0), Run = integer(0), Sex = character(0), Pop_sex = numeric(0))
  }
  
  df_M <- if (length(files_M)) {
    purrr::map_dfr(files_M, summarise_sex_file, sex_label = "M")
  } else {
    tibble(Time = numeric(0), Run = integer(0), Sex = character(0), Pop_sex = numeric(0))
  }
  
  all_runs_sex <- bind_rows(df_F, df_M)
  
  if (!nrow(all_runs_sex)) {
    stop(sprintf("After reading F/M files, no rows remained for %s", scen_dir))
  }
  
  # total pop per Time & Run across sexes (and implicitly across patches, because
  # multiple Patch files per run are all included in all_runs_sex)
  all_runs <- all_runs_sex %>%
    group_by(Time, Run) %>%
    summarise(
      Total_pop = sum(Pop_sex, na.rm = TRUE),
      .groups   = "drop"
    )
  
  K_val <- parse_K_from_tag(scenario_tag)
  
  all_runs %>%
    group_by(Time) %>%
    summarise(
      mean_pop  = mean(Total_pop, na.rm = TRUE),
      q2p5_pop  = quantile(Total_pop, probs = 0.025, na.rm = TRUE),
      q97p5_pop = quantile(Total_pop, probs = 0.975, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(
      mean_norm   = mean_pop  / K_val,
      q2p5_norm   = q2p5_pop  / K_val,
      q97p5_norm  = q97p5_pop / K_val,
      label       = label,
      scenario_tag = scenario_tag,
      K           = K_val
    )
}


# -----------------------------#
# Main plotting function       #
# -----------------------------#
plot_total_pop_byK <- function(base_dir,
                               scenario_tags,
                               labels = NULL,
                               colors = NULL,
                               time_in_years = TRUE,
                               tstart_days = NULL,   # used both for cropping and as time origin
                               tend_days   = NULL,
                               base_size   = 14) {
  
  if (is.null(labels)) labels <- scenario_tags
  if (length(labels) != length(scenario_tags)) {
    stop("Length of 'labels' must match length of 'scenario_tags'.")
  }
  if (!is.null(colors) && length(colors) != length(scenario_tags)) {
    stop("Length of 'colors' must match length of 'scenario_tags'.")
  }
  
  all_sum <- purrr::map_dfr(
    seq_along(scenario_tags),
    ~ summarise_scenario(
      base_dir     = base_dir,
      scenario_tag = scenario_tags[.x],
      label        = labels[.x],
      color        = if (!is.null(colors)) colors[.x] else NULL
    )
  )
  
  # crop by tstart/tend in days (optional)
  if (!is.null(tstart_days)) {
    all_sum <- all_sum %>% dplyr::filter(Time >= tstart_days)
  }
  if (!is.null(tend_days)) {
    all_sum <- all_sum %>% dplyr::filter(Time <= tend_days)
  }
  
  # ----- time origin at tstart_days (if provided) -----
  origin_days <- if (!is.null(tstart_days)) tstart_days else 0
  
  all_sum <- all_sum %>%
    mutate(
      scenario  = factor(label, levels = labels),
      Time_plot = if (time_in_years) {
        (Time - origin_days) / 365        # t = 0 at tstart_days, in years
      } else {
        Time - origin_days                # t = 0 at tstart_days, in days
      }
    )
  
  col_vec <- if (!is.null(colors)) stats::setNames(colors, labels) else NULL
  
  p <- ggplot(all_sum,
              aes(x = Time_plot,
                  y = mean_norm,
                  colour = scenario,
                  fill   = scenario)) +
    geom_ribbon(aes(ymin = q2p5_norm, ymax = q97p5_norm),
                alpha = 0.25, linewidth = 0) +
    geom_line(linewidth = 0.9) +
    labs(
      x = if (time_in_years) "Time since release (years)" else "Time since release (days)",
      y = "Total population / K",
      colour = NULL,
      fill   = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      panel.border      = element_rect(colour = "black", linewidth = 0.8),
      axis.line         = element_line(colour = "black", linewidth = 0.6),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(size = base_size + 3),
      axis.text         = element_text(size = base_size + 1),
      axis.ticks        = element_line(linewidth = 0.6),
      axis.ticks.length = grid::unit(0.3, "cm"),
      
      # legend inside plot
      legend.position      = c(0.70, 0.95),
      legend.justification = c("right", "top"),
      legend.background    = element_rect(fill = "white", colour = "black", linewidth = 0.3),
      legend.title         = element_text(size = base_size + 1),
      legend.text          = element_text(size = base_size),
      legend.key.size      = grid::unit(0.6, "cm"),
      
      plot.margin = grid::unit(c(0.8, 0.8, 0.5, 0.7), "cm")
    )
  
  if (!is.null(col_vec)) {
    p <- p +
      scale_colour_manual(values = col_vec) +
      scale_fill_manual(values   = col_vec)
  }
  
  p
}

# ------------------------------#
# Helper to save with DPI/size  #
# ------------------------------#
save_total_pop_byK_figure <- function(filename,
                                      base_dir,
                                      scenario_tags,
                                      labels,
                                      colors,
                                      time_in_years = TRUE,
                                      tstart_days   = NULL,
                                      tend_days     = NULL,
                                      width_in  = 7,
                                      height_in = 5,
                                      dpi       = 600) {
  
  p <- plot_total_pop_byK(
    base_dir      = base_dir,
    scenario_tags = scenario_tags,
    labels        = labels,
    colors        = colors,
    time_in_years = time_in_years,
    tstart_days   = tstart_days,
    tend_days     = tend_days,
    base_size     = 14
  )
  
  ggsave(filename,
         plot   = p,
         width  = width_in,
         height = height_in,
         units  = "in",
         dpi    = dpi)
}

# ----------------- example usage -----------------
base_dir <- "mgdriveYLE_sweep_5yr_release_target_fertility_byK"
scenario_tags <- c("K1k", "K5k", "K10k", "K50k", "K100k", "K200k")
labels  <- c("K = 1k", "K = 5k", "K = 10k", "K = 50k", "K = 100k", "K = 200k")
colors  <- c("#e41a1c", "#377eb8", "#4daf4a",
             "#984ea3", "#ff7f00", "#a65628")

p <- plot_total_pop_byK(base_dir, scenario_tags, labels, colors,
                        time_in_years = TRUE,
                        tstart_days   = 8*365,
                        tend_days     = 20*365,
                        base_size = 20)
print(p)

save_total_pop_byK_figure(
  filename      = "Fig_E.png",
  base_dir      = base_dir,
  scenario_tags = scenario_tags,
  labels        = labels,
  colors        = colors,
  time_in_years = TRUE,
  tstart_days   = 8*365,
  tend_days     = 28*365,
  width_in      = 6,
  height_in     = 4,
  dpi           = 300
)
