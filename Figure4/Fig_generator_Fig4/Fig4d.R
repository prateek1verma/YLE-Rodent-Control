suppressPackageStartupMessages({
  library(tidyverse)
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
# Safe total-pop calculator
#-----------------------------#
# Rule:
#  - If F_Aggregate exists: total = sum(male-only cols from M_ starting with 'm') + sum(all numeric in F_Aggregate except Time)
#  - Else: total = sum(all numeric columns in M_ except Time)  (many exports already contain both sexes)
sum_numeric_except <- function(df, exclude = "Time") {
  num <- names(df)[vapply(df, is.numeric, logical(1))]
  num <- setdiff(num, exclude)
  if (!length(num)) return(rep(0, nrow(df)))
  rowSums(df[, num, drop = FALSE], na.rm = TRUE)
}

build_total_series_for_run <- function(base_dir, run_id, patch_id, tstart) {
  run_dir   <- resolve_run_dir(base_dir, run_id)
  m_path    <- file.path(run_dir, sprintf("M_Run%s_Patch%s.csv", run_id, patch_id))
  f_path    <- file.path(run_dir, sprintf("F_Aggregate_Run%s_Patch%s.csv", run_id, patch_id))
  m_raw     <- read_ts_csv(m_path)
  
  if (file.exists(f_path)) {
    f_raw <- read_ts_csv(f_path)
    # males: prefer explicit male-only columns if present (start with 'm' or 'M')
    m_cols <- grep("^[mM]", names(m_raw), value = TRUE)
    male_total <- if (length(m_cols)) {
      rowSums(m_raw[, m_cols, drop = FALSE], na.rm = TRUE)
    } else {
      # fallback: assume M_ file already sums (avoid double counting by using only M_ numeric)
      sum_numeric_except(m_raw, exclude = "Time")
    }
    female_total <- sum_numeric_except(f_raw, exclude = "Time")
    total <- male_total + female_total
    Time  <- m_raw$Time  # both should share the same time grid
  } else {
    # No female file: take all numeric from M_ (except Time)
    total <- sum_numeric_except(m_raw, exclude = "Time")
    Time  <- m_raw$Time
  }
  
  tibble(Time = Time, total = total) %>%
    filter(Time >= tstart) %>%
    transmute(Years = (Time - tstart)/365, total = total)
}

#-----------------------------#
# CI summary for totals
#-----------------------------#
summarize_totals_ci <- function(df, center = c("mean","median"), ci = 0.95) {
  center <- match.arg(center)
  a <- (1 - ci)/2; b <- 1 - a
  df %>%
    group_by(Years) %>%
    summarise(
      center = if (center == "mean") mean(total, na.rm = TRUE) else median(total, na.rm = TRUE),
      lower  = quantile(total, probs = a, na.rm = TRUE, type = 7),
      upper  = quantile(total, probs = b, na.rm = TRUE, type = 7),
      .groups = "drop"
    )
}

# ---------- color helper ----------
lighten_hex <- function(col, amt = 0.88) {
  rgb(t(col2rgb(col) * (1 - amt) + 255 * amt) / 255)
}
.validate_scenarios <- function(base_dirs, labels, lighten_amts) {
  stopifnot(length(base_dirs) == length(labels))
  if (is.null(lighten_amts)) lighten_amts <- c(0, rep(0.55, length(base_dirs) - 1))
  if (length(lighten_amts) == 1L) lighten_amts <- rep(lighten_amts, length(base_dirs))
  stopifnot(length(lighten_amts) == length(base_dirs))
  lighten_amts
}

# Map user colors to labels; recycle or error if needed
.map_colors_to_labels <- function(labels, colors = NULL, lighten_amts = NULL, base_col = "#000000") {
  stopifnot(length(labels) >= 1)
  
  # If user provided colors, validate and recycle
  if (!is.null(colors)) {
    # validate colors (throws if invalid)
    invisible(lapply(colors, function(clr) grDevices::col2rgb(clr)))
    if (length(colors) < length(labels)) {
      colors <- rep(colors, length.out = length(labels))
    } else if (length(colors) > length(labels)) {
      warning("More colors than labels; extra colors will be ignored.")
      colors <- colors[seq_along(labels)]
    }
    names(colors) <- labels
    return(colors)
  }
  
  # Fallback: build from lighten_amts (as before)
  if (is.null(lighten_amts)) lighten_amts <- c(0, rep(0.55, length(labels) - 1))
  if (length(lighten_amts) == 1L) lighten_amts <- rep(lighten_amts, length(labels))
  stopifnot(length(lighten_amts) == length(labels))
  
  lighten_hex <- function(col, amt = 0.88) {
    rgb(t(col2rgb(col) * (1 - amt) + 255 * amt) / 255)
  }
  cols <- vapply(lighten_amts, function(a) lighten_hex(base_col, amt = a), character(1))
  names(cols) <- labels
  cols
}


#-----------------------------#
# Build per-scenario totals CI
#-----------------------------#
build_totals_multi_scenario_ci <- function(
    base_dirs, labels, patch_id = "001", run_ids = NULL, tstart = round(8*365),
    center = c("mean","median"), lighten_amts = NULL
) {
  center <- match.arg(center)
  lighten_amts <- .validate_scenarios(base_dirs, labels, lighten_amts)
  
  all_ci <- list(); pal <- c()
  base_col <- "#000000"
  scen_cols <- vapply(lighten_amts, function(a) lighten_hex(base_col, amt = a), character(1))
  names(scen_cols) <- labels
  
  for (i in seq_along(base_dirs)) {
    bd <- base_dirs[i]; lab <- labels[i]
    rids <- if (is.null(run_ids)) list_run_ids(bd) else run_ids
    totals_i <- purrr::map_dfr(rids, ~build_total_series_for_run(bd, .x, patch_id, tstart) %>%
                                 mutate(run_id = .x, scenario = lab))
    ci_i <- totals_i %>%
      summarize_totals_ci(center = center, ci = 0.95) %>%
      mutate(scenario = lab)
    all_ci[[i]] <- ci_i
  }
  list(ci = bind_rows(all_ci), cols = scen_cols)
}

#-----------------------------#
# Plot: mean total ±95% CI
#-----------------------------#
make_totals_mean_CI_compare <- function(
    base_dirs, labels,
    patch_id = "001",
    run_ids = NULL,
    tstart = round(8*365),
    center = c("mean","median"),
    lighten_amts = NULL,
    colors = NULL,                # <-- NEW: user palette (vector of hex or names)
    outfile = NULL,
    width = 9, height = 6, dpi = 600,
    y_log = FALSE, tEnd = NULL,
    legpos_tag = c(0.5,0.5),
    alpha_ci = 0.22, lw_mean = 1.4,
    csv_outfile = NULL          # <--- NEW
) {
  built <- build_totals_multi_scenario_ci(
    base_dirs, labels, patch_id, run_ids, tstart, center, lighten_amts
  )
  ci_all <- built$ci
  
  # ---- write CI summary to CSV ----
  if (is.null(csv_outfile)) {
    csv_outfile <- file.path(
      normalizePath(getwd()),
      sprintf("TOTALS_CI_Patch%s_t%s_%s.csv",
              patch_id, tstart,
              paste(gsub("[^A-Za-z0-9]+","", labels), collapse = "-"))
    )
  }
  readr::write_csv(ci_all, csv_outfile)
  message("Saved CI summary: ", csv_outfile)
  
  
  # ---- FIX 1: enforce scenario factor levels ----
  ci_all <- ci_all %>%
    dplyr::mutate(scenario = factor(scenario, levels = labels))
  
  # Build scenario colors: prefer user-supplied 'colors'
  # ---- FIX 2: named palette aligned with labels ----
  if (!is.null(colors)) {
    # if user forgot to name the vector, name by labels
    if (is.null(names(colors))) names(colors) <- labels
    scen_cols <- colors[labels]         # ensure same order as labels
  } else {
    scen_cols <- .map_colors_to_labels(labels, colors = NULL, lighten_amts = lighten_amts)
  }
  
  base_theme <- theme_bw(base_size = 24) +
    theme(
      legend.position      = legpos_tag,
      legend.justification = c(0,1),
      legend.title         = element_blank(),
      legend.text          = element_text(size = 24),
      axis.title           = element_text(size = 26),
      axis.text            = element_text(size = 22),
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_line(linewidth = 0.25),
      legend.background    = element_rect(fill   = scales::alpha("white", 0.5))  # 0–1, higher = more opaque
    )
  
  p <- ggplot() +
    geom_ribbon(
      data = ci_all,
      aes(x = Years, ymin = lower, ymax = upper, fill = scenario),
      alpha = alpha_ci, colour = NA
    ) +
    geom_line(
      data = ci_all,
      aes(x = Years, y = center, color = scenario),
      linewidth = lw_mean
    ) +
    scale_color_manual(values = scen_cols, breaks = labels) +
    scale_fill_manual(values = scen_cols, breaks = labels) +
    labs(x = "Time (years)", y = if (y_log) "Total (log10)" else "Total Population") +
    base_theme
  
  x_breaks <- if (is.null(tEnd)) pretty(ci_all$Years) else pretty(c(0, tEnd))
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
      breaks = scales::breaks_pretty(n = 6),
      labels = function(x) paste0(x/1000, "k") 
    ) + coord_cartesian(ylim = c(0, NA))
    
    # p <- p + scale_y_continuous(
    #   labels = scales::label_number(big.mark = ""),
    #   breaks = scales::breaks_pretty(n = 6)
    # ) + coord_cartesian(ylim = c(0, NA))
  }
  
  if (is.null(outfile)) {
    outfile <- file.path(
      normalizePath(getwd()),
      sprintf("FIG_total_only_Patch%s_t%s_%s.png",
              patch_id, tstart, paste(gsub("[^A-Za-z0-9]+","",labels), collapse = "-"))
    )
  }
  ggsave(outfile, p, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved totals: ", outfile)
  
  invisible(list(plot = p, outfile = outfile, ci = ci_all, colors = scen_cols))
}


base_dirs <- (c("mgdriveSIT_rel_100per_5000",
                "mgdrivefsRIDL_rel_10per_500", 
                "mgdriveYLE_rel_5per_250"))

labels_tag <- c("5000 SMR/month","500 fsRRDL/month", "250 YLE/month")


col_tags <- c("#CC79A7","#E69F00","#0072B2")
#          c("#E69F00", "#009E73", "#0072B2")
  # Color blind freindly: c("#000000", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")
  # Color Brewer source: c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#e41a1c")

make_totals_mean_CI_compare(
  base_dirs = base_dirs,
  labels    = labels_tag,
  patch_id  = "001",
  tstart    = round(8*365),
  center    = "mean",
  tEnd      = 20,
  legpos_tag = c(0.5,0.98),
  y_log = FALSE,
  colors    = col_tags,  # <-- user palette
  alpha_ci = 0.20, lw_mean = 1.4,
  csv_outfile = "totals_mean_CI_all_scenarios.csv"          # <--- NEW
)


## UNCOMMENT THIS FUNCTION to complie and write all the data of the plot in a csv file
# #!/usr/bin/env Rscript
# 
# suppressPackageStartupMessages({
#   library(tidyverse)
# })
# 
# # ---------------------------------------------------------
# # Function: analyse CI CSV from make_totals_mean_CI_compare
# # ---------------------------------------------------------
# # csv_file    : path to CI csv (with columns Years, center, lower, upper, scenario)
# # target_time : numeric, time (in years) at which you want total pop + CI
# # scenarios   : NULL = all scenarios; character vector = subset of scenarios
# # returns     : tibble with rows for target/max/min per scenario
# analyse_totals_ci <- function(csv_file,
#                               target_time,
#                               scenarios = NULL) {
#   if (!file.exists(csv_file)) {
#     stop("CSV not found: ", csv_file)
#   }
#   
#   df <- readr::read_csv(csv_file, show_col_types = FALSE)
#   
#   req_cols <- c("Years", "center", "lower", "upper", "scenario")
#   missing <- setdiff(req_cols, names(df))
#   if (length(missing)) {
#     stop("CSV is missing required columns: ", paste(missing, collapse = ", "))
#   }
#   
#   if (!is.null(scenarios)) {
#     df <- df %>% filter(scenario %in% scenarios)
#     if (nrow(df) == 0) {
#       stop("No rows left after filtering by scenarios: ",
#            paste(scenarios, collapse = ", "))
#     }
#   }
#   
#   # Split per scenario
#   split(df, df$scenario) %>%
#     imap_dfr(function(d, scen) {
#       d <- d %>% arrange(Years)
#       
#       # --- value at target_time (nearest time point) ---
#       idx_target <- which.min(abs(d$Years - target_time))
#       row_target <- d[idx_target, ]
#       
#       # --- max and min over time (by center) ---
#       idx_max <- which.max(d$center)
#       idx_min <- which.min(d$center)
#       
#       row_max <- d[idx_max, ]
#       row_min <- d[idx_min, ]
#       
#       tibble(
#         scenario = scen,
#         stat     = c("target_time", "max", "min"),
#         Years    = c(row_target$Years, row_max$Years, row_min$Years),
#         center   = c(row_target$center, row_max$center, row_min$center),
#         lower    = c(row_target$lower,  row_max$lower,  row_min$lower),
#         upper    = c(row_target$upper,  row_max$upper,  row_min$upper)
#       )
#     })
# }
# 
# # ---------------------------------------------------------
# # Example usage (edit these paths/values for your case)
# # ---------------------------------------------------------
# if (sys.nframe() == 0) {
#   # INPUTS: change as needed
#   csv_path    <- "totals_mean_CI_all_scenarios.csv"
#   target_time <- 10   # e.g. 10 years
#   # scenarios <- c("5000 SMR/month", "250 YLE/month")  # optional subset
#   scenarios <- NULL   # NULL = all scenarios in the file
#   
#   res <- analyse_totals_ci(
#     csv_file    = csv_path,
#     target_time = target_time,
#     scenarios   = scenarios
#   )
#   
#   cat("\nSummary of total population and CI:\n")
#   print(res, n = Inf)
#   
#   # Optionally, write results to CSV
#   out_csv <- sub("\\.csv$", paste0("_summary_t", target_time, ".csv"), csv_path)
#   readr::write_csv(res, out_csv)
#   cat("\nSaved summary to:", out_csv, "\n")
# }
