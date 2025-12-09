# ==========================================
# Sensitivity analysis via BRT (per paper)
# ==========================================
rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(gbm)        # engine used by dismo::gbm.step
  library(dismo)      # gbm.step (paper's method)
  library(ggplot2)
})

# --------- paths / basic settings ----------
root <- "mgdriveYLE_sweep_Latin4"           # << change only if your root differs
index_csv <- file.path(root, "lhs_param_sets.csv")
stopifnot(dir.exists(root), file.exists(index_csv))

days_per_year     <- 365
release_start_day <- 8 * 365            # << match your simulations
fem_elim_threshold <- 1                 # final female total <= this => eradicated (successful)
Neq <- 10000

# --------- helpers ----------
read_ts_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  stopifnot("Time" %in% names(df))
  df <- df %>% mutate(Time = suppressWarnings(as.numeric(Time))) %>% arrange(Time)
  if (anyNA(df$Time)) stop(sprintf("Non-numeric 'Time' in %s", basename(path)))
  df
}

# Extract summary per run: event (elim) and event time (years) if successful.
# Also record follow-up length in days (from release start) for completeness.
summarize_run_female <- function(run_dir, fem_thresh, release_start_day, days_per_year) {
  f_path <- list.files(run_dir, pattern = "^F_Aggregate_Run\\d+_Patch\\d+\\.csv$", full.names = TRUE)
  if (length(f_path) == 0) {
    warning(sprintf("Missing F_Aggregate csv in: %s", run_dir))
    return(list(final_f = NA_real_, elim = NA, tte_years = NA_real_, t_end_days = NA_real_))
  }
  fdf <- read_ts_csv(f_path[1])
  
  female_cols <- grep("^f", names(fdf), value = TRUE)
  female_total <- if ("female_total" %in% names(fdf)) fdf$female_total else {
    if (length(female_cols)) rowSums(fdf[, female_cols, drop = FALSE]) else rep(NA_real_, nrow(fdf))
  }
  
  final_f <- suppressWarnings(as.numeric(female_total[length(female_total)]))
  elim    <- is.finite(final_f) && final_f <= fem_thresh
  
  # starting index at release
  idx0 <- which(fdf$Time >= release_start_day)
  if (!length(idx0)) {
    return(list(final_f = final_f, elim = elim, tte_years = NA_real_, t_end_days = NA_real_))
  }
  idx0 <- idx0[1]
  
  # event time if eradication occurs after release start
  hit <- which(female_total[idx0:length(female_total)] <= fem_thresh)
  if (length(hit)) {
    hit_idx <- idx0 + hit[1] - 1
    tte_years <- (fdf$Time[hit_idx] - release_start_day) / days_per_year
    t_end_days <- fdf$Time[hit_idx] - release_start_day
  } else {
    tte_years <- NA_real_
    t_end_days <- tail(fdf$Time, 1) - release_start_day
  }
  
  list(final_f = final_f, elim = elim, tte_years = tte_years, t_end_days = t_end_days)
}

# --------- assemble replicate-level dataset ----------
# Expect index to have columns: rel, fy, pq, fl, mu, j, c, tag, set_id
index <- readr::read_csv(index_csv, show_col_types = FALSE)

param_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
param_dirs <- param_dirs[basename(param_dirs) %in% index$tag]

dataset <- purrr::map_dfr(param_dirs, function(pdir){
  tag <- basename(pdir)
  run_dirs <- list.dirs(pdir, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("/\\d{3}$", run_dirs)]
  if (!length(run_dirs)) return(tibble())
  res <- purrr::map(run_dirs, ~ summarize_run_female(.x, fem_elim_threshold, release_start_day, days_per_year))
  finals_f <- vapply(res, `[[`, numeric(1), "final_f")
  elims    <- vapply(res, `[[`, logical(1), "elim")
  ttes     <- vapply(res, `[[`, numeric(1), "tte_years")
  tend     <- vapply(res, `[[`, numeric(1), "t_end_days")
  tibble(tag = tag, elim = as.integer(elims), t_elim_years = ttes, t_end_days = tend, final_f = finals_f)
}) %>%
  left_join(index, by = "tag") %>%
  filter(!is.na(elim))  # keep only completed replicates

# Save replicate-level BRT input (reproducibility)
readr::write_csv(dataset, file.path(root, "brt_input_replicate_level.csv"))

# ===============================
# Choose predictors + optional filters
# ===============================
selected_params <- c("rel","fy","p","q","fl","fs","j","mu","c")   # adjust if you want all 8
range_filter <- list(
  rel = c(0.01, 0.20),
  fy  = c(0.20, 1.00),
  p  = c(0.00, 0.20),
  q  = c(0.00, 0.20),
  fl = c(0.00, 0.00),
  fs = c(0.80, 1.00),
  mu  = c(0.80, 1.00),
  j   = c(0.00, 1.00),
  c   = c(0.80, 1.00)
)

# bounds <- tribble(
#   ~param, ~low,  ~high, ~scale, ~type,
#   "rel",   0.10,  0.10,  "lin",  "cont",   # release proportion of adult males
#   "fy",    0.20,  1.00,  "lin",  "cont",   # relative fitness of Y-males (0..1)
#   "p",     0.00,  0.50,  "lin",  "cont",   # Proportion of dominant NHEJ, repair rate p=q
#   "q",     0.00,  0.50,  "lin",  "cont",   # Proportion of dominant NHEJ, repair rate p=q
#   "fl",    0.50,  1.00,  "lin",  "cont",   # female lethality flag (0/1); map U(0,1) to {0,1}
#   "mu",    0.80,  1.00,   "lin",  "cont",  # Probability of mutation (in wild-type homozygotes) 
#   "j",     0.00,  1.00,  "lin",  "cont",   # Probability of NHEJ given joining (in heterozygotes)
#   "c",     0.80,  1.00,  "lin",  "cont"    # Probability of cleavage (in heterozygotes)
# )

# Apply row filters (only columns present)
if (length(range_filter)) {
  for (nm in intersect(names(range_filter), names(dataset))) {
    rng <- range_filter[[nm]]
    if (!is.null(rng) && length(rng) == 2 && all(is.finite(rng))) {
      dataset <- dataset %>% dplyr::filter(.data[[nm]] >= rng[1], .data[[nm]] <= rng[2])
    }
  }
}

keep_cols <- intersect(selected_params, c("rel","fy","p","q","fs","mu","j","c"))
stopifnot(length(keep_cols) >= 1)

# ----------------------------------------------------
# Utilities for gbm.step: zero-variance drop + col ids
# ----------------------------------------------------
# drop_nzv_and_build_x <- function(df, y_col) {
#   # make everything numeric *except* y, then drop NA rows
#   df2 <- df %>%
#     dplyr::mutate(!!y_col := as.numeric(.data[[y_col]])) %>%
#     dplyr::mutate(dplyr::across(-all_of(y_col), \(x) as.numeric(x))) %>%
#     tidyr::drop_na()
#   
#   # drop zero-variance predictors
#   nzv <- vapply(df2[names(df2) != y_col], function(x) sd(x) > 0, logical(1))
#   keep <- c(y_col, names(nzv)[nzv])
#   
#   # RETURN A BASE data.frame (not a tibble!)
#   df3 <- as.data.frame(df2[, keep, drop = FALSE])
#   
#   list(
#     data = df3,
#     gbm.y = which(names(df3) == y_col),
#     gbm.x = which(names(df3) != y_col)
#   )
# }

drop_nzv_and_build_x <- function(df, y_col) {
  # make everything numeric *except* y, then drop NA rows
  df2 <- df %>%
    dplyr::mutate(!!y_col := as.numeric(.data[[y_col]])) %>%
    dplyr::mutate(dplyr::across(-all_of(y_col), \(x) as.numeric(x))) %>%
    tidyr::drop_na()
  
  if(nrow(df2) == 0) stop("Filtering resulted in 0 rows of data. Check your range_filter bounds.")
  
  # drop zero-variance predictors (handle NA sd)
  nzv <- vapply(df2[names(df2) != y_col], function(x) {
    s <- sd(x)
    !is.na(s) && s > 0
  }, logical(1))
  
  keep <- c(y_col, names(nzv)[nzv])
  
  # RETURN A BASE data.frame
  df3 <- as.data.frame(df2[, keep, drop = FALSE])
  
  list(
    data = df3,
    gbm.y = which(names(df3) == y_col),
    gbm.x = which(names(df3) != y_col)
  )
}


# ===============================
# BINOMIAL BRT — probability
# (gbm.step, per paper)
# ===============================
mod_prob <- dataset %>%
  dplyr::select(elim, all_of(keep_cols)) %>%
  dplyr::mutate(elim = as.integer(elim))  # 0/1 numeric

prob_prep <- drop_nzv_and_build_x(mod_prob, "elim")

set.seed(1)
gbm_prob <- dismo::gbm.step(
  data            = prob_prep$data,               # base data.frame
  gbm.x           = prob_prep$gbm.x,
  gbm.y           = prob_prep$gbm.y,
  family          = "bernoulli",
  tree.complexity = 3,
  learning.rate   = 0.01,
  bag.fraction    = 0.75,
  n.folds         = 5,
  verbose         = FALSE
)

imp_prob <- gbm::summary.gbm(gbm_prob, plotit = FALSE) |>
  tibble::as_tibble() |>
  dplyr::transmute(param = var, rel_inf = rel.inf, outcome = "Probability of eradication")


# Save
readr::write_csv(imp_prob, file.path(root, "BRT_importance_probability_SELECTED.csv"))

# ==========================================
# POISSON BRT — time to eradication (successes only)
# (gbm.step, per paper)
# ==========================================
succ_reps <- dataset %>%
  dplyr::filter(elim == 1, is.finite(t_elim_years), t_elim_years >= 0) %>%
  dplyr::mutate(t_days = pmax(1L, as.integer(round(t_elim_years * days_per_year)))) %>%
  dplyr::select(t_days, all_of(keep_cols)) %>%
  dplyr::rename(y_time = t_days)

if (nrow(succ_reps) >= 10) {
  time_prep <- drop_nzv_and_build_x(succ_reps, "y_time")
  
  set.seed(2)
  gbm_time <- dismo::gbm.step(
    data            = time_prep$data,           # base data.frame
    gbm.x           = time_prep$gbm.x,
    gbm.y           = time_prep$gbm.y,
    family          = "poisson",
    tree.complexity = 3,
    learning.rate   = 0.01,
    bag.fraction    = 0.75,
    n.folds         = 5,
    verbose         = FALSE
  )
  
  imp_time <- gbm::summary.gbm(gbm_time, plotit = FALSE) |>
    tibble::as_tibble() |>
    dplyr::transmute(param = var, rel_inf = rel.inf, outcome = "Time to eradication")
  readr::write_csv(imp_time, file.path(root, "BRT_importance_time_SELECTED.csv"))
} else {
  message("Not enough successes for the Poisson time-to-eradication model; skipping.")
}


# ===============================
# Define the Expression Map
# ===============================
# This maps your short codes (param) to the fancy plot expressions
expression_map <- c(
  "rel" = "Release Size",
  "fy"  = expression(paste("YLE Fitness (", italic(f)[y], ")")),
  "p"   = expression(paste("Het. Resistance (", italic(p), ")")),
  "q"   = expression(paste("Hom. Resistance (", italic(q), ")")),
  "fs"  = expression(paste("Fem. Sterility (", italic(f)[S], ")")), # Assuming fs is sterility
  "mu"  = expression(paste("Hom. Cleavage (", italic(mu), ")")),
  "j"   = expression(paste("NHEJ Fraction (", italic(j), ")")),
  "c"   = expression(paste("Het. Cleavage (", italic(c), ")"))
)

# ===============================
# Prepare Data
# ===============================
imp_all <- dplyr::bind_rows(imp_prob, imp_time) %>%
  # Create the numeric sign value for the mirror effect
  dplyr::mutate(sign_val = ifelse(outcome == "Probability of eradication", -rel_inf, rel_inf)) %>%
  # Filter to keep only the parameters present in 'keep_cols'
  # (Assuming 'keep_cols' contains the short names like "fy", "c", etc.)
  dplyr::filter(param %in% keep_cols)

# ===============================
# Order Factors (Crucial for sorting)
# ===============================
# Calculate total influence to sort the bars
order_levels <- imp_all %>%
  dplyr::group_by(param) %>% 
  dplyr::summarise(tot = sum(rel_inf), .groups = "drop") %>%
  dplyr::arrange(tot) %>% 
  dplyr::pull(param)

# Set the factor levels on the 'param' column directly
imp_all$param <- factor(imp_all$param, levels = order_levels)

# =======================================================
# Paper-Quality Mirrored Bar Plot (Tornado Plot)
# =======================================================

max_abs <- 100  # full 0–100% range on each side

p <- ggplot(imp_all, aes(x = sign_val, y = param, fill = outcome)) +
  # 1. Bars: Add a thin white border for better definition
  geom_col(width = 0.7, color = "white", size = 0.2) +
  
  # 2. Reference Line: Make it subtle
  geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
  
  # 3. Y-Axis: Use your expression map for proper math notation
  scale_y_discrete(labels = expression_map) +
  
  # 4. X-Axis: Symmetrical expansion to give breathing room
  scale_x_continuous(
    labels = function(x) abs(x), 
    expand = expansion(mult = c(0.1, 0.1)) # 10% padding on sides
  ) +
  
  # 5. Colors: High contrast Red/Blue (ColorBrewer Set1/RdBu inspired)
  scale_fill_manual(values = c("Probability of eradication" = "#B2182B",
                               "Time to eradication"        = "#2166AC")) +
  
  # 6. Labels: Professional casing
  labs(
    x = "Relative influence (%)", 
    y = NULL,
    # Titles are often removed in papers in favor of the figure caption. 
    # If you keep it, keep it simple.
    title = NULL 
  ) +
  
  # 7. Theme: Start with classic and override for 'Paper Look'
  theme_classic(base_size = 18) + # Base size 18 is great for papers
  theme(
    # --- Legend ---
    # Position: c(0,0) is bottom-left, c(1,1) is top-right.
    # Since 'rel' (top bar) is huge, put legend at the BOTTOM (0.15)
    legend.position = c(0.8, 0.15), 
    legend.justification = c(0.5, 0.5),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey80", size = 0.3),
    legend.margin = margin(6, 6, 6, 6),
    legend.text = element_text(size = 14),
    
    # --- Axis ---
    axis.text.y = element_text(color = "black", margin = margin(r = 5)), # Math expressions
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    
    # --- Gridlines ---
    # Optional: Add faint vertical grids to help read percentages
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    
    # --- Plot Margins ---
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# Save with dimensions suitable for a full-width or half-width figure
ggsave(file.path(root, "BRT_relative_influence_mirrored_PUBLICATION.png"), p,
       width = 9, height = 6, dpi = 600) # 600 dpi for high-quality print
# ===============================
# Plot: mirrored relative influence (optional)
# ===============================
# label_map <- c(
#   rel = "Release Size",
#   fy  = "f_y (y-males fitness)",
#   p  = "p (NHEJ dom. 1)",
#   q  = "q (NHEJ dom. 2)",
#   fs  = "fL (female sterility)",
#   mu  = "mu (mutation in W/W)",
#   j   = "j (NHEJ, het.)",
#   c   = "c (cleavage, het.)"
# )
# 
# imp_all <- dplyr::bind_rows(imp_prob, imp_time) %>%
#   dplyr::mutate(label = dplyr::recode(param, !!!label_map),
#                 sign_val = ifelse(outcome == "Probability of eradication", -rel_inf, rel_inf)) %>%
#   # keep only labels for predictors we actually used
#   dplyr::filter(label %in% dplyr::recode(keep_cols, !!!label_map))
# 
# order_levels <- imp_all %>%
#   dplyr::group_by(label) %>% dplyr::summarise(tot = sum(rel_inf), .groups = "drop") %>%
#   dplyr::arrange(tot) %>% dplyr::pull(label)
# imp_all$label <- factor(imp_all$label, levels = order_levels)
# 
# p <- ggplot(imp_all, aes(x = sign_val, y = label, fill = outcome)) +
#   geom_col(width = 0.7) +
#   geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
#   scale_x_continuous(labels = function(x) abs(x), expand = expansion(mult = c(0.05, 0.05))) +
#   scale_fill_manual(values = c("Probability of eradication" = "#B2182B",
#                                "Time to eradication"       = "#2166AC")) +
#   labs(x = "Relative influence (%)", y = NULL,
#        title = "Relative influence (Boosted Regression Tree)") +
#   theme_classic(base_size = 12) +
#   theme(legend.title = element_blank(),
#         plot.title = element_text(face = "bold"))
# ggsave(file.path(root, "BRT_relative_influence_mirrored_SELECTED.png"), p,
#        width = 8, height = 5, dpi = 300)
# print(p)

# ===============================
# Optional: save aggregated view (not required for 1 run/param set)
# ===============================
agg <- dataset %>%
  group_by(tag, across(all_of(keep_cols))) %>%
  summarise(
    trials       = dplyr::n(),  # <--- Change n() to dplyr::n()
    successes    = sum(elim == 1, na.rm = TRUE),
    prop_elim    = if (trials > 0) successes / trials else NA_real_,
    mean_t_years = if (successes > 0) mean(t_elim_years[elim == 1 & is.finite(t_elim_years)], na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  filter(is.finite(prop_elim))
