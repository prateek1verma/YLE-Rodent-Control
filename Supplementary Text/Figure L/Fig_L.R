rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
  library(patchwork)
})

# ------------------------------- #
# User settings
# ------------------------------- #

scenario_info <- tibble::tribble(
  ~strategy,       ~root,                      ~fy_keep,     ~x_label,
  "Boosted YLE (XS=0.5)",  "mgdrive_boosted_YLE_rel_xshred_0p50",   0.8,          "Release size (% of adult transgenic males released)",
  "Boosted YLE (XS=0.8)",  "mgdrive_boosted_YLE_rel_xshred_0p80",   0.8,          "Release size (% of adult transgenic males released)",
  "YLE",          "mgdrive_YLE_rel_high_reso_100rep",           0.8,          "Release size (% of adult transgenic males released)"
)

runs_expect    <- sprintf("%03d", 1:100)  # adjust if needed
days_per_year  <- 365
burn_in_years  <- 8
analysis_years <- 50
K              <- 10000

out_file_png <- "FigL_v1.png"
out_file_pdf <- "combined_peak_total_min_female_all_strategies.pdf"
out_file_csv <- "combined_peak_total_min_female_all_strategies.csv"

# ------------------------------- #
# Helpers
# ------------------------------- #

read_ts_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"Time" %in% names(df)) stop(sprintf("'Time' column not found in %s", basename(path)))
  df %>%
    mutate(Time = suppressWarnings(as.numeric(Time))) %>%
    arrange(Time)
}

parse_num_from_tag <- function(tag, key) {
  m <- stringr::str_match(tag, paste0("(^|_)", key, "([-0-9p]+)"))[, 3]
  if (is.na(m)) return(NA_real_)
  as.numeric(sub("p", ".", m, fixed = TRUE))
}

parse_param_tags <- function(tag) {
  tibble(
    rel = parse_num_from_tag(tag, "rel"),
    fy  = parse_num_from_tag(tag, "fy")
  )
}

read_fem_male_totals <- function(run_dir) {
  f_path <- list.files(run_dir, pattern = "^F_.*Aggregate.*\\.csv$", full.names = TRUE)
  m_path <- list.files(run_dir, pattern = "^M_.*\\.csv$", full.names = TRUE)
  
  if (length(f_path) == 0 || length(m_path) == 0) {
    stop(sprintf("Missing F or M csv in: %s", run_dir))
  }
  
  fdf <- read_ts_csv(f_path[1])
  mdf <- read_ts_csv(m_path[1])
  
  f_cols <- grep("^f", names(fdf), value = TRUE)
  m_cols <- grep("^m", names(mdf), value = TRUE)
  
  female_total <- if ("female_total" %in% names(fdf)) {
    fdf$female_total
  } else {
    rowSums(fdf[, f_cols, drop = FALSE], na.rm = TRUE)
  }
  
  male_total <- rowSums(mdf[, m_cols, drop = FALSE], na.rm = TRUE)
  
  inner_join(
    tibble(Time = fdf$Time, female_total = as.numeric(female_total)),
    tibble(Time = mdf$Time, male_total = as.numeric(male_total)),
    by = "Time"
  ) %>%
    arrange(Time) %>%
    mutate(total_adults = female_total + male_total)
}

run_extrema_window <- function(run_dir) {
  df <- read_fem_male_totals(run_dir)
  
  t0 <- burn_in_years * days_per_year
  t1 <- (burn_in_years + analysis_years) * days_per_year
  
  dfw <- df %>% filter(Time >= t0, Time <= t1)
  
  tibble(
    run_dir        = run_dir,
    run_id         = basename(run_dir),
    max_total      = max(dfw$total_adults, na.rm = TRUE),
    min_female     = min(dfw$female_total, na.rm = TRUE)
  )
}

read_strategy <- function(strategy, root, fy_keep = NA_real_) {
  stopifnot(dir.exists(root))
  
  param_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  
  map_dfr(param_dirs, function(pdir) {
    tag  <- basename(pdir)
    pars <- parse_param_tags(tag)
    
    if (!is.na(fy_keep)) {
      if (!is.finite(pars$fy) || abs(pars$fy - fy_keep) > 1e-8) return(tibble())
    }
    
    run_dirs <- list.dirs(pdir, recursive = FALSE, full.names = TRUE)
    run_dirs <- run_dirs[basename(run_dirs) %in% runs_expect]
    
    if (!length(run_dirs)) return(tibble())
    
    map_dfr(run_dirs, function(rd) {
      out <- tryCatch(run_extrema_window(rd), error = function(e) tibble())
      if (!nrow(out)) return(tibble())
      
      out %>%
        mutate(
          strategy = strategy,
          tag      = tag,
          rel      = pars$rel,
          fy       = pars$fy,
          rel_pct  = 100 * pars$rel
        )
    })
  })
}

k_lab <- function(x) {
  paste0(round(x / 1000, 1), "k")
}

# ------------------------------- #
# Read and summarise all scenarios
# ------------------------------- #

per_run_tbl <- pmap_dfr(
  scenario_info,
  function(strategy, root, fy_keep, x_label) {
    read_strategy(strategy, root, fy_keep)
  }
)

stopifnot(nrow(per_run_tbl) > 0)

summ_tbl <- per_run_tbl %>%
  pivot_longer(
    cols = c(max_total, min_female),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      max_total  = "Maximum total population",
      min_female = "Minimum female population"
    ),
    metric = factor(metric, levels = c("Maximum total population", "Minimum female population")),
    strategy = factor(strategy, levels = scenario_info$strategy)
  ) %>%
  group_by(strategy, rel_pct, metric) %>%
  summarise(
    n    = n(),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    se   = sd / sqrt(n),
    tval = qt(0.975, df = pmax(n - 1, 1)),
    lo   = mean - tval * se,
    hi   = mean + tval * se,
    .groups = "drop"
  )

readr::write_csv(summ_tbl, out_file_csv)

# ------------------------------- #
# Publication figure
# ------------------------------- #

# "Max total" = "#e41a1c", "Min females" = "#377eb8"

# cols <- c(
#   "Peak total population" = "#D55E00",
#   "Minimum female population" = "#0072B2"
# )

cols <- c(
  "Maximum total population" = "#e41a1c",
  "Minimum female population" = "#377eb8"
)


p <- ggplot(
  summ_tbl,
  aes(x = rel_pct, y = mean, color = metric, fill = metric, shape = metric)
) +
  geom_hline(yintercept = K, linetype = "dashed", color = "grey35", linewidth = 0.6) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
  # geom_errorbar(
  #   aes(ymin = lo, ymax = hi),
  #   width = 0.4,
  #   linewidth = 0.6,
  #   alpha = 0.8
  # ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8, stroke = 0.7) +
  facet_wrap(~ strategy, nrow = 1, scales = "free_x") +
  scale_color_manual(values = cols, name = NULL) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(
    breaks = breaks_pretty(n = 5),
    labels = label_number(accuracy = 1)
  ) +
  scale_y_continuous(
    labels = k_lab,
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    x = "Release size (% of adult males at K)",
    y = "Population size"
    # caption = "Lines show means across stochastic runs; ribbons show 95% confidence intervals. Dashed line denotes K = 10,000."
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 17),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.position = "top",
    legend.text = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.ticks.length = unit(2.5, "mm")
    # plot.caption = element_text(size = 11, hjust = 0),
    # plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )
 
print(p)

ggsave(out_file_png, p, width = 10, height = 4.5, dpi = 600, bg = "white")
ggsave(out_file_pdf, p, width = 10, height = 4.5, bg = "white")


# Assumes you already created `summ_tbl` from the combined script

K <- 10000

k_lab <- function(x) paste0(round(x / 1000, 1), "k")

strategy_cols <- c(
  "YLE"    = "#d95f02",
  "Boosted YLE (XS=0.8)" = "#7570b3",
  "Boosted YLE (XS=0.5)" = "#1b9e77"
)

#1b9e77
#d95f02
#7570b3

#1b9e77
#d95f02
#7570b3

#e41a1c
#377eb8
#4daf4a

make_metric_plot <- function(metric_name, ylab, y_limits = NULL, x_limits = NULL) {
  
  plot_df <- summ_tbl %>%
    filter(metric == metric_name) %>%
    mutate(strategy = factor(strategy, levels = c("YLE", "Boosted YLE (XS=0.8)","Boosted YLE (XS=0.5)")))
  
  ggplot(plot_df, aes(x = rel_pct, y = mean,
                      color = strategy, fill = strategy,
                      shape = strategy)) +
    geom_hline(yintercept = K, linetype = "dashed",
               color = "grey35", linewidth = 0.7) +
    geom_hline(yintercept = K/2, linetype = "dotted",
               color = "grey35", linewidth = 0.7) +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                alpha = 0.25, color = NA) +
    geom_line(linewidth = 1.0) +
    # geom_point(size = 3.2, stroke = 0.8) +
    geom_point(size = 2.5, stroke = 0.7) +
    scale_color_manual(values = strategy_cols) +
    scale_fill_manual(values = strategy_cols) +
    scale_shape_manual(values = c("YLE" = 16, "Boosted YLE (XS=0.8)" = 17,"Boosted YLE (XS=0.5)" = 17)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
    scale_y_continuous(labels = k_lab) +
    coord_cartesian(ylim = y_limits, xlim = x_limits) +
    labs(
      x = "Release size (% of adult males at K)",
      y = ylab
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "top",
      legend.title = element_blank()
    )
}

p1 <- make_metric_plot(
  "Maximum total population",
  "Maximum total population",
  y_limits = c(9900, 11000),
  x_limits = c(0, 4.5)
)

p2 <- make_metric_plot(
  "Minimum female population",
  "Minimum female population",
  y_limits = c(-100, 5200),
  x_limits = c(0, 4.5)
)

p2 <- p2 +
 # theme(legend.position = "none") +
  theme(
    legend.position = c(0.90, 0.90),   # x, y inside panel b legend.position = c(0.43, 0.90),
    # legend.justification = c("left", "bottom"),
    legend.justification = c("right", "top"),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.75),
      colour = NA
    ),
    legend.key.size = unit(0.55, "cm"),
    legend.text = element_text(size = 15)
  ) + 
  annotate("text", x = 3, y = 5000,
           label = "Initial females",
           hjust = 0, vjust = -0.5,
           size = 5, color = "grey30")

p1 <- p1 +
  theme(legend.position = "none") +
  # theme(
  #   legend.position = c(0.40, 0.95),   # x, y inside panel b legend.position = c(0.43, 0.90),
  #   # legend.justification = c("left", "bottom"),
  #   legend.justification = c("center", "top"),
  #   legend.direction = "horizontal",
  #   legend.title = element_blank(),
  #   legend.background = element_rect(
  #     fill = scales::alpha("white", 0.75),
  #     colour = NA
  #   ),
  #   legend.key.size = unit(0.55, "cm"),
  #   legend.text = element_text(size = 15)
  # ) + 
  annotate("text", x = 3, y = 10000,
         label = "Initial total (K)",
         hjust = 0, vjust = -0.5,
         size = 5, color = "grey30")

p_combined <- (p1 | p2) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 30, face = "bold")
  ) 


print(p_combined)

ggsave(
  "Fig_combined_two_panel.pdf",
  p_combined,
  width = 12,
  height = 5,
  bg = "white"
)

