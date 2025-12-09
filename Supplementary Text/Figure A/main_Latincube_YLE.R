# ============================
# LHS-driven MGDrivE YLE runs
# ============================
rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(MGDrivEmouse)
  library(tidyverse)
  library(lhs)
  library(progressr)
  library(future.apply)
})

# ---------- top-level config ----------
outFolder <- "mgdriveYLE_sweep_Latin4"
if (dir.exists(outFolder)) unlink(outFolder, recursive = TRUE, force = TRUE)
dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)
cat("New folder:", normalizePath(outFolder), "\n")

# Simulation horizon & replicates
tMax   <- 58 * 365       # days
nRep   <- 1             # replicates per parameter set
tstart <- 8 * 365        # release start day
trelease <- 12*10      # 10 years

# Demography (edit as needed)
bioParameters <- list(
  betaK = 6, litters = 7.5, tGest = 19, tNursing = 23, tAdo = 37,
  muAd = (1/690), theta = 22.4
)

sitesNumber <- 1
Neq  <- 10000
cap  <- Neq
moveMat <- matrix(1, 1, 1)
batchMigration <- basicBatchMigration(batchProbs = 0, sexProbs = c(.5, .5), numPatches = sitesNumber)

# Toxicology (same as before; keep or remove)
exposure  <- 1.26
expTime   <- 7
toxCurve  <- (1/(1+exp(-3.2*(log(exposure)-log(1.774)))))/expTime
toxTime   <- c(10^6)

# Inheritance cube source
source("generate_YLE_inheritance_cube.R")

# ---------- PARAMETERS & LHS ----------
# Define bounds (EDIT these to your study ranges)
bounds <- tribble(
  ~param, ~low,  ~high, ~scale, ~type,
  "rel",   0.01,  0.20,  "lin",  "cont",   # release proportion of adult males
  "fy",    0.20,  1.00,  "lin",  "cont",   # relative fitness of Y-males
  "p",     0.00,  0.20,  "lin",  "cont",   # NHEJ parameter p
  "q",     0.00,  0.20,  "lin",  "cont",   # NHEJ parameter q
  "fl",    0.00,  0.00,  "lin",  "cont",   # female lethality (fixed at 0)
  "fs",    0.80,  1.00,  "lin",  "cont",   # female sterility (0.8–1.0)
  "mu",    0.80,  1.00,  "lin",  "cont",   # mutation probability
  "j",     0.00,  1.00,  "lin",  "cont",   # NHEJ prob. given joining
  "c",     0.80,  1.00,  "lin",  "cont"    # cleavage probability
)


n_sets <- 10   # LHS parameter sets (EDIT)
set.seed(2025)

lhs_raw <- randomLHS(n_sets, nrow(bounds))

scale_col <- function(z, low, high, scale){
  if (scale == "log") return(exp(log(low) + z*(log(high)-log(low))))
  low + z*(high - low)
}

param_sets <- map2_dfc(
  as.data.frame(lhs_raw), seq_len(nrow(bounds)),
  ~ scale_col(.x, bounds$low[.y], bounds$high[.y], bounds$scale[.y])
) %>%
  set_names(bounds$param) %>%
  mutate(
    set_id = row_number()
  ) %>%
  relocate(set_id)

# quick sanity
# summary(param_sets$fl)

# Optional: freeze any parameter (e.g., keep pq=0 for all)
# param_sets$pq <- 0
# param_sets$c <- 1
# param_sets$mu <- 1

# ---------- tagging helpers ----------
num_tag <- function(x) gsub("\\.", "p", formatC(x, format = "f", digits = 4))

build_tag <- function(rel, fy, p, q, fl, fs, mu, j, c) {
  paste0(
    "rel", num_tag(rel),
    "_fy", num_tag(fy),
    "_p",  num_tag(p),
    "_q",  num_tag(q),
    "_fl", num_tag(fl),
    "_fs", num_tag(fs),
    "_mu", num_tag(mu),
    "_j",  num_tag(j),
    "_c",  num_tag(c)
  )
}


# ---------- single-parameter-set runner ----------
run_param_set <- function(row){
  rel_prop <- row$rel
  fy       <- row$fy
  p_val    <- row$p
  q_val    <- row$q
  fl_val   <- row$fl   # fixed at 0 from bounds
  fs_val   <- row$fs   # LHS: 0.8–1
  mu_val   <- row$mu
  j_val    <- row$j
  c_val    <- row$c
  
  # --- inheritance cube ---
  alleleFitness_vec <- c(W = 1, H = 1, R = 1, Y = 1, X = 1, y = fy)
  cube <<- generate_YLE_inheritance_cube(
    c = c_val, j = j_val, mu = mu_val,
    p = p_val, q = q_val,                  # FIX: use p,q not undefined pq
    female_lethality  = fl_val,            # = 0
    female_sterility  = fs_val,            # ∈ [0.8, 1]
    Haploinsufficient = TRUE,
    muAd              = 1/690,
    alleleFitness     = alleleFitness_vec
  )
  
  # --- releases (adult males) ---
  releasesParameters_SP <- list(
    releasesStart     = tstart,
    releasesNumber    = trelease,
    releasesInterval  = 30,
    releaseProportion = round(cap * rel_prop * 0.5)
  )
  ReleasesVector_SP <- generateReleaseVector(
    driveCube = cube,
    releasesParameters = releasesParameters_SP
  )
  
  patchReleases <- replicate(
    n = sitesNumber,
    expr = list(
      maleReleases         = NULL,
      femaleReleases       = NULL,
      eggReleases          = NULL,
      matedFemaleReleases  = NULL
    ),
    simplify = FALSE
  )
  patchReleases[[1]]$maleReleases <- ReleasesVector_SP
  
  # --- MGDrivE parameters ---
  netPar <- parameterizeMGDrivE(
    runID = 1, simTime = tMax, nPatch = sitesNumber,
    beta = bioParameters$betaK, litters = bioParameters$litters,
    tGest = bioParameters$tGest, tNursing = bioParameters$tNursing,
    tAdo = bioParameters$tAdo,
    k = cap, theta = bioParameters$theta, inheritanceCube = cube,
    AdPopRatio_M = matrix(c(1, 0), 1, 2, dimnames = list(NULL, c("mYWW", "fXWW"))),
    AdPopRatio_F = matrix(c(0, 1), 1, 2, dimnames = list(NULL, c("mYWW", "fXWW")))
  )
  
  # --- folders ---
  tag <- build_tag(rel_prop, fy, p_val, q_val, fl_val, fs_val, mu_val, j_val, c_val)
  param_dir <- file.path(outFolder, tag)
  if (dir.exists(param_dir)) unlink(param_dir, recursive = TRUE, force = TRUE)
  dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)
  
  run_dirs <- file.path(param_dir, formatC(1:nRep, width = 3, format = "d", flag = "0"))
  invisible(lapply(run_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  
  # --- stochastic engine & runs ---
  setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)
  MGDrivESim <- Network$new(
    params          = netPar,
    driveCube       = cube,
    patchReleases   = patchReleases,
    migrationMale   = moveMat,
    migrationFemale = moveMat,
    migrationBatch  = batchMigration,
    directory       = run_dirs,
    verbose         = TRUE,
    toxMort         = toxCurve,
    toxInt          = toxTime
  )
  MGDrivESim$multRun(verbose = TRUE)
  
  # --- post-process each replicate ---
  for (i in seq_len(nRep)) {
    rd <- run_dirs[i]
    splitOutput(readDir = rd, remFile = TRUE, verbose = FALSE)
    aggregateFemales(
      readDir   = rd,
      genotypes = cube$genotypesID,
      remFile   = TRUE,
      verbose   = FALSE
    )
  }
  
  tag
}


# ---------- parallel execution ----------
workers <- max(1, round(parallel::detectCores() * 0.25))

handlers(global = TRUE)
handlers("txtprogressbar")   # or "progress" or "pbcol" for fancier output
plan(multisession, workers = workers)
cat("Using workers:", workers, "\n")

# Keep an index of all sets (for traceability)
index <- param_sets %>%
  mutate(
    tag = pmap_chr(
      select(., rel, fy, p, q, fl, fs, mu, j, c),
      ~ build_tag(..1, ..2, ..3, ..4, ..5, ..6, ..7, ..8, ..9)
    )
  )
readr::write_csv(index, file.path(outFolder, "lhs_param_sets.csv"))


with_progress({
  p <- progressor(along = 1:nrow(param_sets))
  results <- future_lapply(split(param_sets, param_sets$set_id), function(row){
    out <- run_param_set(row)
    p()    # tick progress
    out
  })
})

# Run all sets
# tags <- future_lapply(split(param_sets, param_sets$set_id), \(row) run_param_set(row))
# cat("Completed", length(tags), "parameter sets.\n")
