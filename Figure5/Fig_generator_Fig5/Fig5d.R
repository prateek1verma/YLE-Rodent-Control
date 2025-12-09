## ============================
##  Generalised non-WT barplot
## ============================
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggnewscale) # Required for dual legends
})

#-----------------------------#
# Helpers
#-----------------------------#

# Discover run directories that look like "001", "002", ... or "Run001", "Run002", ...
list_run_dirs <- function(base_dir) {
  base_dir <- path.expand(base_dir)
  if (!dir.exists(base_dir)) stop("Base directory not found: ", base_dir)
  dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  runs <- dirs[basename(dirs) %>% str_detect("^(\\d{3}|Run\\d{3})$")]
  if (length(runs) == 0) stop("No run folders found under: ", base_dir)
  runs
}

# Read a time series CSV and ensure numeric Time
read_ts_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"Time" %in% names(df)) stop("'Time' column not found in: ", path)
  df <- df %>%
    mutate(Time = suppressWarnings(as.numeric(Time))) %>%
    arrange(Time)
  if (anyNA(df$Time)) stop("Non-numeric values in 'Time' column of: ", path)
  df
}

# Extract run_id and patch_id from filename "M_Run001_Patch003.csv" or "F_Run001_Patch003.csv"
parse_ids_from_filename <- function(file) {
  bn <- basename(file)
  # Matches M_Run001_Patch001.csv, M_Aggregate_Run001_Patch001.csv,
  # and same for F_...
  m <- str_match(bn, "^[MF](?:_Aggregate)?_Run(\\d{3})_Patch(\\d{3})\\.csv$")
  if (is.na(m[1, 1])) return(NULL)
  list(run_id = m[1, 2], patch_id = m[1, 3])
}


# Compute non-WT abundance statistic (timepoint or peak) given male & female files
# Compute abundance statistic (timepoint or peak) for a set of genotypes
# geno_vec: character vector of column names to use (across M/F files).
#           If NULL, falls back to "all non-WT" = all cols except mYWW, fXWW.
stat_from_files <- function(male_file,
                            female_file,
                            mode = c("timepoint", "peak"),
                            t_years = NULL,
                            days_per_year = 365,
                            geno_vec = NULL) {
  mode <- match.arg(mode)
  
  dm <- read_ts_csv(male_file)
  df <- read_ts_csv(female_file)
  
  if (!is.null(geno_vec)) {
    geno_m <- intersect(geno_vec, names(dm))
    geno_f <- intersect(geno_vec, names(df))
  } else {
    # original non-WT behaviour
    geno_m <- setdiff(names(dm), c("Time", "mYWW", "fXWW"))
    geno_f <- setdiff(names(df), c("Time", "mYWW", "fXWW"))
  }
  
  if (length(geno_m) == 0 && length(geno_f) == 0) {
    stop("No selected genotype columns found in:\n  ", male_file, "\n  ", female_file)
  }
  
  if (mode == "timepoint") {
    if (is.null(t_years)) stop("t_years must be supplied for mode = 'timepoint'.")
    t_target <- t_years * days_per_year
    
    im <- which.min(abs(dm$Time - t_target))
    i_f <- which.min(abs(df$Time - t_target))
    
    vm <- if (length(geno_m)) rowSums(dm[im, geno_m, drop = FALSE]) else 0
    vf <- if (length(geno_f)) rowSums(df[i_f, geno_f, drop = FALSE]) else 0
    
    return(as.numeric(vm + vf))
  }
  
  if (mode == "peak") {
    vm_t <- if (length(geno_m)) rowSums(dm[, geno_m, drop = FALSE]) else 0
    vf_t <- if (length(geno_f)) rowSums(df[, geno_f, drop = FALSE]) else 0
    return(max(vm_t + vf_t))
  }
  
  stop("Unknown mode: ", mode)
}

#-----------------------------#
# Main: multi-design barplot
#-----------------------------#

# base_dirs  : character vector of design folders
# labels_tag : labels for legend (same length as base_dirs)
# col_tags   : colors for designs (same length as base_dirs)
# mode       : "timepoint" (needs t_years) or "peak"
plot_nonwt_multidesign_by_patch <- function(
    base_dirs,
    labels_tag,
    col_tags,
    mode = c("timepoint", "peak"),
    t_years = NULL,
    days_per_year = 365,
    add_errorbars = TRUE,
    alpha_fill = 1.0,
    outfile = NULL,
    genotype_map = NULL
) {
  
  mode <- match.arg(mode)
  
  # 1. Data Loading & Processing
  # -------------------------------------------------------
  # (Identical to previous logic, just collapsing for brevity)
  all_design_data <- vector("list", length(base_dirs))
  
  for (k in seq_along(base_dirs)) {
    bdir <- base_dirs[k]
    lab  <- labels_tag[k]
    
    geno_vec <- if (!is.null(genotype_map)) genotype_map[[lab]] else NULL
    
    run_dirs <- list_run_dirs(bdir)
    male_files <- map(run_dirs, ~ list.files(.x, pattern = "^M(?:_Aggregate)?_Run\\d{3}_Patch\\d{3}\\.csv$", full.names = TRUE)) %>% unlist()
    
    if (length(male_files) == 0) stop("No male files found in: ", bdir)
    
    df_design <- map_dfr(male_files, function(mf) {
      ids <- parse_ids_from_filename(mf)
      if (is.null(ids)) return(NULL)
      
      bn   <- basename(mf)
      dirn <- dirname(mf)
      core <- sub("^M(?:_Aggregate)?_", "", bn)
      ff1 <- file.path(dirn, paste0("F_Aggregate_", core))
      ff2 <- file.path(dirn, paste0("F_",           core))
      ff  <- if (file.exists(ff1)) ff1 else if (file.exists(ff2)) ff2 else stop("Female file missing")
      
      stat_val <- stat_from_files(mf, ff, mode, t_years, days_per_year, geno_vec)
      
      tibble(design = lab, run_id = ids$run_id, patch_id = ids$patch_id, value = stat_val)
    })
    all_design_data[[k]] <- df_design
  }
  
  df_all <- bind_rows(all_design_data)
  
  # Aggregate stats
  df_summary <- df_all %>%
    group_by(design, patch_id) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val   = sd(value,  na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(
      patch_num = as.integer(patch_id),
      patch_lab = paste0("Patch ", patch_num)
    )
  
  # 2. Manual Positioning Calculation (With Gap)
  # -------------------------------------------------------
  design_order <- labels_tag
  # Configuration
  bar_width <- 0.12        # Width of individual bars
  group_gap <- 0.05        # Space between the Blue group and Red group
  
  # We assume the input 'labels_tag' is ordered: 
  # First 3 = Blue Group (Mig=2), Last 3 = Red Group (Mig=1)
  
  # Calculate offsets relative to the Patch Center (0)
  # The "Red" group starts at +gap/2
  # The "Blue" group ends at -gap/2
  
  # Red Offsets (Right side)
  r1 <- (group_gap / 2) + (bar_width / 2)
  r2 <- (group_gap / 2) + (bar_width * 1.5)
  r3 <- (group_gap / 2) + (bar_width * 2.5)
  offsets_red <- c(r1, r2, r3)
  
  # Blue Offsets (Left side) - Mirror of Red
  b3 <- -((group_gap / 2) + (bar_width / 2))
  b2 <- -((group_gap / 2) + (bar_width * 1.5))
  b1 <- -((group_gap / 2) + (bar_width * 2.5))
  offsets_blue <- c(b1, b2, b3)
  
  # Combine and name them to match 'design' levels
  all_offsets <- c(offsets_blue, offsets_red)
  names(all_offsets) <- labels_tag 
  
  # Apply offsets
  df_summary <- df_summary %>%
    mutate(
      x_pos = patch_num + all_offsets[as.character(design)],
      group_id = ifelse(design %in% labels_tag[1:3], "Blue", "Red")
    )
  
  # 3. Aesthetics & Legends
  # -------------------------------------------------------
  
  # Common Legend Labels (Reused for both blocks)
  # This creates the text: "YLE fitness, f_y = X"
  common_labels <- c(
    expression(paste("YLE fitness, ", f[y] == 1)),
    expression(paste("YLE fitness, ", f[y] == 0.8)),
    expression(paste("YLE fitness, ", f[y] == 0.5))
  )
  
  # Subset Data
  df_blue <- df_summary %>% filter(group_id == "Blue")
  df_red  <- df_summary %>% filter(group_id == "Red")
  
  # Subset Colors
  cols_blue <- col_tags[1:3]
  cols_red  <- col_tags[4:6]
  
  names(cols_blue) <- design_order[1:3]
  names(cols_red)  <- design_order[4:6]
  
  y_lab <- if (mode == "timepoint") paste0("Mean invasive abundance\n(t = ", t_years, " years)") else "Mean peak \nYLE males abundance"
  
  base_theme <- theme_classic(base_size = 20) +
    theme(
      # Stack legends vertically in top-right
      legend.position       = c(0.98, 0.98),
      legend.justification  = c("right", "top"),
      legend.box            = "vertical", # Stack the two legends
      legend.background     = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.margin         = margin(6, 10, 6, 6),
      
      # Text Styling
      legend.title          = element_text(size = 19, face = "bold", hjust = 0), # Header size
      legend.text           = element_text(size = 17, hjust = 0),
      legend.key.size       = unit(0.5, "cm"),
      legend.spacing.y      = unit(0.25, "cm"), # Gap between legend items
      
      axis.title            = element_text(size = 24, face = "bold", color = "black"),
      axis.text             = element_text(size = 24, color = "black"),
      axis.line             = element_line(linewidth = 1.0),
      axis.ticks            = element_line(linewidth = 1.0),
      panel.grid.major.y    = element_line(color = "grey85", linewidth = 0.5),
      panel.grid.major.x    = element_blank()
    )
  
  p <- ggplot() +
    # --- BLUE BARS (Mig = 2) ---
    geom_col(
      data = df_blue,
      aes(x = x_pos, y = mean_val, fill = design),
      width = bar_width, # Use the fixed width defined above
      color = "black", linewidth = 0.4
    ) +
    # Error bars for Blue
    geom_errorbar(
      data = df_blue,
      aes(x = x_pos, ymin = pmax(0, mean_val - sd_val), ymax = mean_val + sd_val),
      width = bar_width * 0.4, linewidth = 0.6, color = "grey20"
    ) +
    # Scale for Blue
    scale_fill_manual(
      name = "Migration rate = 2 mice/day", # Legend Header 1
      values = cols_blue,
      labels = common_labels # Using the math expressions
    ) +
    
    # --- RESET SCALE ---
    new_scale_fill() + 
    
    # --- RED BARS (Mig = 1) ---
    geom_col(
      data = df_red,
      aes(x = x_pos, y = mean_val, fill = design),
      width = bar_width, # Use the fixed width defined above
      color = "black", linewidth = 0.4
    ) +
    # Error bars for Red
    geom_errorbar(
      data = df_red,
      aes(x = x_pos, ymin = pmax(0, mean_val - sd_val), ymax = mean_val + sd_val),
      width = bar_width * 0.4, linewidth = 0.6, color = "grey20"
    ) +
    # Scale for Red
    scale_fill_manual(
      name = "Migration rate = 1 mouse/day", # Legend Header 2
      values = cols_red,
      labels = common_labels # Same labels, different colors
    ) +
    
    # --- Layout ---
    scale_x_continuous(
      breaks = sort(unique(df_summary$patch_num)),
      labels = function(x) paste("Patch", x)
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 6),
      labels = label_number(big.mark = ","),
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(x = NULL, y = y_lab) +
    base_theme
  
  if (is.null(outfile)) {
    suffix <- if (mode == "timepoint") paste0("_t", t_years, "yr") else "_peak"
    outfile <- file.path(path.expand(getwd()), paste0("Figure_SplitLegend", suffix, ".png"))
  }
  
  ggsave(outfile, p, width = 10, height = 6, dpi = 300, bg = "white")
  message("Saved plot: ", outfile)
  
  invisible(list(plot = p, outfile = outfile))
}

# =======================================================
# EXECUTION
# =======================================================

base_dirs <- c(
  "mgdriveYLE_rel_5per_250_5patches_mig_2daily_fy1",
  "mgdriveYLE_rel_5per_250_5patches_mig_2daily",
  "mgdriveYLE_rel_5per_250_5patches_mig_2daily_fy0p5",
  "mgdriveYLE_rel_5per_250_5patches_mig_1daily_fy1",
  "mgdriveYLE_rel_5per_250_5patches_mig_1daily",
  "mgdriveYLE_rel_5per_250_5patches_mig_1daily_fy0p5"
)

labels_tag <- c(
  "2 mouse/day, f_y = 1", "2 mouse/day, f_y = 0.8", "2 mouse/day, f_y = 0.5",
  "1 mouse/day, f_y = 1", "1 mouse/day, f_y = 0.8", "1 mouse/day, f_y = 0.5"
)

# col_tags <- c(
#   "#2B6595", "#377EB8", "#A1B9D9",  # Blue Scale
#   "#B80F15", "#E41A1C", "#FF9294"   # Red Scale
# )

# Paste this into your code
# col_tags <- c(
#   # Migration rate = 2 (Steel Blue)
#   "#2171B5",  # High
#   "#6BAED6",  # Mid
#   "#BDD7E7",  # Low
# 
#   # Migration rate = 1 (Crimson Red)
#   "#CB181D",  # High
#   "#FB6A4A",  # Mid
#   "#FCAE91"   # Low
# )

# Paste this into your code
col_tags <- c(
  # Migration rate = 2 (Purple Gradient)
  "#542788",  # High Fitness (Deep Purple)
  "#756BB1",  # Mid Fitness (Medium Purple)
  "#BCBDDC",  # Low Fitness (Lavender)
  
  # Migration rate = 1 (Emerald/Green Gradient)
  "#006D2C",  # High Fitness (Deep Emerald)
  "#31A354",  # Mid Fitness (Grass Green)
  "#A1D99B"   # Low Fitness (Pale Green)
)

# # Paste this into your code
# col_tags <- c(
#   # Migration rate = 2 (Blue Gradient)
#   "#08519C",  # High Fitness (Dark Blue)
#   "#3182BD",  # Mid Fitness (Medium Blue)
#   "#9ECAE1",  # Low Fitness (Light Blue)
#   
#   # Migration rate = 1 (Orange Gradient)
#   "#A63603",  # High Fitness (Dark Orange/Rust)
#   "#E6550D",  # Mid Fitness (Vibrant Orange)
#   "#FDAE6B"   # Low Fitness (Pale Orange)
# )

genotype_map <- list(
  "2 mouse/day, f_y = 1"   = c("myWW","myHW","myRW","myHH","myHR","myRR"),
  "2 mouse/day, f_y = 0.8" = c("myWW","myHW","myRW","myHH","myHR","myRR"),
  "2 mouse/day, f_y = 0.5" = c("myWW","myHW","myRW","myHH","myHR","myRR"),
  "1 mouse/day, f_y = 1"   = c("myWW","myHW","myRW","myHH","myHR","myRR"),
  "1 mouse/day, f_y = 0.8" = c("myWW","myHW","myRW","myHH","myHR","myRR"),
  "1 mouse/day, f_y = 0.5" = c("myWW","myHW","myRW","myHH","myHR","myRR")
)

res <- plot_nonwt_multidesign_by_patch(
  base_dirs = base_dirs,
  labels_tag = labels_tag,
  col_tags = col_tags,
  mode = "peak",
  add_errorbars = TRUE,
  genotype_map = genotype_map
)