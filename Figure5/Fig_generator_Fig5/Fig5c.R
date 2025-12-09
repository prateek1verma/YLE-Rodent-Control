## ============================
##  Generalised non-WT barplot
## ============================
suppressPackageStartupMessages({
  library(tidyverse)   # ggplot2, dplyr, tidyr, readr, stringr, purrr
  library(scales)
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
    alpha_fill = 0.9,
    outfile = NULL,
    genotype_map = NULL
) {
  
  mode <- match.arg(mode)
  
  if (length(base_dirs) != length(labels_tag) ||
      length(base_dirs) != length(col_tags)) {
    stop("base_dirs, labels_tag, and col_tags must have the same length.")
  }
  
  all_design_data <- vector("list", length(base_dirs))
  
  for (k in seq_along(base_dirs)) {
    bdir <- base_dirs[k]
    lab  <- labels_tag[k]
    
    # pick genotypes for this design if supplied
    geno_vec <- NULL
    if (!is.null(genotype_map)) {
      if (!lab %in% names(genotype_map)) {
        stop("No genotype entry in 'genotype_map' for design label: ", lab)
      }
      geno_vec <- genotype_map[[lab]]
    }
    
    run_dirs <- list_run_dirs(bdir)
    
    # All male files in this design
    male_files <- map(run_dirs, ~ list.files(
      .x,
      pattern = "^M(?:_Aggregate)?_Run\\d{3}_Patch\\d{3}\\.csv$",
      full.names = TRUE
    )) %>% unlist()
    
    if (length(male_files) == 0) {
      stop("No male files (M*_Run###_Patch###.csv) found in: ", bdir)
    }
    
    df_design <- map_dfr(male_files, function(mf) {
      ids <- parse_ids_from_filename(mf)
      if (is.null(ids)) return(NULL)
      
      bn   <- basename(mf)
      dirn <- dirname(mf)
      core <- sub("^M(?:_Aggregate)?_", "", bn)
      
      ff1 <- file.path(dirn, paste0("F_Aggregate_", core))
      ff2 <- file.path(dirn, paste0("F_",           core))
      
      if (file.exists(ff1)) {
        ff <- ff1
      } else if (file.exists(ff2)) {
        ff <- ff2
      } else {
        stop("Matching female file not found for: ", mf)
      }
      
      stat_val <- stat_from_files(
        male_file    = mf,
        female_file  = ff,
        mode         = mode,
        t_years      = t_years,
        days_per_year = days_per_year,
        geno_vec     = geno_vec
      )
      
      tibble(
        design   = lab,
        run_id   = ids$run_id,
        patch_id = ids$patch_id,
        value    = stat_val
      )
    })
    all_design_data[[k]] <- df_design
  }
  
  df_all <- bind_rows(all_design_data)
  if (nrow(df_all) == 0) stop("No valid data found.")
  
  # Aggregate across runs
  df_summary <- df_all %>%
    dplyr::group_by(design, patch_id) %>%
    dplyr::summarise(
      n_runs   = dplyr::n(),
      mean_val = mean(value, na.rm = TRUE),
      sd_val   = sd(value,  na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(
      patch_num = as.integer(patch_id),
      patch_lab = paste0("Patch ", patch_num)
    ) %>%
    arrange(patch_num, design)
  
  # -------------------------------------------------------
  # AESTHETICS UPDATE: Theme & Labels
  # -------------------------------------------------------
  
  # Use classic theme as base (cleaner than bw), customize from there
  base_theme <- theme_classic(base_size = 20) +
    theme(
      # Legend at bottom
      legend.position       = "bottom",
      legend.direction      = "horizontal",
      legend.title          = element_blank(),
      legend.text           = element_text(size = 18),
      legend.key.size       = unit(0.7, "cm"),
      
      # Axis typography
      axis.title            = element_text(size = 24, face = "bold", color = "black"),
      axis.text             = element_text(size = 24, color = "black"),
      
      # Grid lines: Keep horizontal for reading values, remove vertical
      panel.grid.major.y    = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.major.x    = element_blank(),
      panel.grid.minor      = element_blank(),
      
      # Box and background
      axis.line             = element_line(linewidth = 0.8, color = "black"),
      plot.margin           = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  # Factor order
  patch_levels  <- df_summary %>% distinct(patch_lab, patch_num) %>%
    arrange(patch_num) %>% pull(patch_lab)
  design_levels <- labels_tag
  
  df_summary <- df_summary %>%
    mutate(
      patch_lab = factor(patch_lab, levels = patch_levels),
      design    = factor(design, levels = design_levels)
    )
  
  pd <- position_dodge(width = 0.8)
  
  # Fixed Label Logic: Removed the duplicate overwrite
  y_lab <- if (mode == "timepoint") {
    paste0("Mean transgene carrier abundance\n(t = ", t_years, " years)")
  } else {
    "Mean peak transgene \n carrier abundance"
  }
  
  p <- ggplot(df_summary,
              aes(x = patch_lab,
                  y = mean_val,
                  fill = design)) +
    # Added black border (color="black") for better definition
    geom_col(position = pd, alpha = alpha_fill, width = 0.7, color = "black", linewidth = 0.3) +
    
    {if (add_errorbars)
      geom_errorbar(
        aes(ymin = pmax(0, mean_val - sd_val),
            ymax = mean_val + sd_val),
        position = pd,
        width = 0.2,
        linewidth = 0.6,
        color = "grey20" # Slightly softer black for error bars
      )
    } +
    labs(
      x = NULL, # "Patch Number" is usually redundant if labels say "Patch 1"
      y = y_lab
    ) +
    scale_fill_manual(values = setNames(col_tags, labels_tag)) +
    
    # Fix floating bars: expand=c(0, 0.1) starts exactly at 0 and gives 10% headroom
    scale_y_continuous(
      labels = function(x) ifelse(x == 0, "0", paste0(x / 1000, "k")),
      expand = expansion(mult = c(0, 0.1)) 
    ) + 
    base_theme + 
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) # Arrange legend in 2 rows
  
  
  # Output logic
  if (is.null(outfile)) {
    suffix <- if (mode == "timepoint") paste0("_t", t_years, "yr") else "_peak"
    outfile <- file.path(path.expand(getwd()),
                         paste0("Figure5c_intact_homing_By_Patch_multidesign", suffix, ".png"))
  }
  
  ggsave(outfile, p, width = 10, height = 7, dpi = 300, bg = "white")
  message("Saved plot: ", outfile)
  
  invisible(list(summary = df_summary, runs = df_all, plot = p, outfile = outfile))
}


base_dirs <- c("mgdrive_Homing_5patch_250_mig_1daily",
               "mgdrive_Homing_non_ideal_5patch_250_mig_1daily",
               "mgdrive_Xshredder_ideal_5patch_250_mig_1daily",
               "mgdrive_Xshredder_nonideal_5patch_250_mig_1daily", 
               "mgdriveYLE_rel_5per_250_5patches_mig_1daily")

labels_tag <- c("Ideal HGD","Non-ideal HGD", "Ideal X-Shredder",
                "Non-ideal X-Shredder", "YLE (250/month)")

col_tags <- c("#377eb8", "#4daf4a", "#984ea3", 
              "#ff7f00","#e41a1c")

## 1) At 20 years (non-WT abundance)
# res_time <- plot_nonwt_multidesign_by_patch(
#   base_dirs   = base_dirs,
#   labels_tag  = labels_tag,
#   col_tags    = col_tags,
#   mode        = "timepoint",
#   t_years     = 50,
#   add_errorbars = TRUE
# )

genotype_map <- list(
  "Ideal HGD"           = c("mYHW","mYHH","mYHR","mYHN",
                            "fXHW","fXHH","fXHR","fXHN"),
  "Non-ideal HGD"       = c("mYHW","mYHH","mYHR","mYHN",
                            "fXHW","fXHH","fXHR","fXHN"),
  "Ideal X-Shredder"    = c("mXA","mRA"),
  "Non-ideal X-Shredder"= c("mXA","mRA"),
  "YLE (250/month)"       = c("myWW","myHW","myRW","myHH","myHR","myRR")
)


res_peak <- plot_nonwt_multidesign_by_patch(
  base_dirs    = base_dirs,
  labels_tag   = labels_tag,
  col_tags     = col_tags,
  mode         = "peak",
  add_errorbars = TRUE,
  genotype_map = genotype_map
)

# res_time <- plot_nonwt_multidesign_by_patch(
#   base_dirs    = base_dirs,
#   labels_tag   = labels_tag,
#   col_tags     = col_tags,
#   mode         = "timepoint",
#   t_years      = 50,
#   add_errorbars = TRUE,
#   genotype_map = genotype_map
# )
