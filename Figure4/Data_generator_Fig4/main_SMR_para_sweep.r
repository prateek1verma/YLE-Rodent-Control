####################
# Load libraries
####################
rm(list = ls()); gc(); save.image()   # overwrite .RData with empty workspace
suppressPackageStartupMessages({
  library(MGDrivEmouse)
  library(BBmisc)
  library(parallel)
  library(pbapply)   
  library(pbmcapply)
})

####################
# Top-level output
####################
outFolder <- "mgdriveSIT_sweep_10yr"
if (dir.exists(outFolder)) {
  unlink(outFolder, recursive = TRUE, force = TRUE)
  cat("Folder deleted:", outFolder, "\n")
}
dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)
cat("New folder created:", outFolder, "\n")

############################
# Simulation horizon & reps
############################
tMax  <- 58 * 365       # days
nRep  <- 10            # Monte Carlo iterations per parameter set
tstart <- 8 * 365       # release start (days since 0)

####################
# Biology & demography
####################
bioParameters <- list(
  betaK = 6, litters = 7.5, tGest = 19, tNursing = 23, tAdo = 37,
  muAd = (1/690), theta = 22.4
)

sitesNumber <- 1
Neq <- 10000
cap <- Neq
moveMat <- matrix(1, 1, 1)  # no dispersal
batchMigration <- basicBatchMigration(batchProbs = 0, sexProbs = c(.5, .5), numPatches = sitesNumber)

####################
# Parameter grids
####################
# Release proportion of ADULT MALES: 1%..10%
rel_grid   <- seq(0, 15, length.out = 31) # 10^(seq(log10(0.1), log10(20), length.out = 10))
  # seq(10, 2000, by = 0.01)

# Fitness of y-males (fy): set what you want here
fy_grid <- seq(1, 1, by = 0.3) # c(1.00, 0.80, 0.60, 0.40, 0.2, 0.1)   # example; 1=no cost, <1 reduces fitness

# Repair rates (p = q), 0..0.5
# pq_grid    <- seq(0.00, 0.50, by = 0.10)
pq_grid <- 0

# Female lethality flag
fl_grid    <- c(0)

# Constants for cube but INCLUDED in folder names (change if desired)
mu_val <- 1.0
j_val  <- 0.0
c_val  <- 1.0

####################
# Source inheritance cube
####################
# source("generate_YLE_inheritance_cube_v2.R")
source("generate_SIT_inheritance_cube.R")

####################
# Toxicology (kept from your script)
####################
exposure  <- 1.26
expTime   <- 7
toxCurve  <- (1/(1+exp(-3.2*(log(exposure)-log(1.774)))))/expTime
toxTime   <- c(10^6)  # effectively absent

####################
# Helpers
####################
# Pretty/portable number tags (0.1 -> p10)
# -- helper to make "0.10" -> "p10" etc.
num_tag <- function(x) gsub("\\.", "p", formatC(x, format = "f", digits = 2))

# build folder tag
build_tag <- function(rel, fy, pq, fl, mu, j, c) {
  paste0(
    "rel", num_tag(rel),
    "_fy", num_tag(fy),
    "_pq", num_tag(pq),
    "_fl", fl,
    "_mu", num_tag(mu),
    "_j",  num_tag(j),
    "_c",  num_tag(c)
  )
}

run_param_set <- function(rel_prop, fy, pq, fl) {
  # ===== inheritance cube for this param-set =====
  fy_val <- calcOmega(mu = bioParameters$muAd, lifespanReduction = fy)
  cube <<- generate_SIT_inheritance_cube(sterility = 1)
  cube$omega["my"] <<- fy_val
  
  
  # ===== releases (adult males; rel_prop of total adult pop × 0.5) =====
  releasesParameters_SP <- list(
    releasesStart    = tstart,
    releasesNumber   = round(12*10),
    releasesInterval = 30,
    releaseProportion= round(cap * rel_prop * 0.5)
  )
  ReleasesVector_SP <- generateReleaseVector(driveCube = cube, releasesParameters = releasesParameters_SP)
  patchReleases <- replicate(
    n = sitesNumber,
    expr = list(maleReleases = NULL, femaleReleases = NULL, eggReleases = NULL, matedFemaleReleases = NULL),
    simplify = FALSE
  )
  patchReleases[[1]]$maleReleases <- ReleasesVector_SP
  
  # ===== parameterized network (unchanged) =====
  netPar <- parameterizeMGDrivE(
    runID = 1, simTime = tMax, nPatch = sitesNumber,
    beta = bioParameters$betaK, litters = bioParameters$litters,
    tGest = bioParameters$tGest, tNursing = bioParameters$tNursing, tAdo = bioParameters$tAdo,
    k = cap, theta = bioParameters$theta, inheritanceCube = cube,
    AdPopRatio_M = matrix(c(1, 0), 1, 2, dimnames = list(NULL, c("mY", "fX"))),
    AdPopRatio_F = matrix(c(0, 1), 1, 2, dimnames = list(NULL, c("mY", "fX")))
  )
  
  # ===== folders: mgdriveYLE_sweep/<tag>/<rep> =====
  tag <- build_tag(rel_prop, fy, pq, fl, mu_val, j_val, c_val)
  param_dir <- file.path(outFolder, tag)
  
  # clean & recreate the param-set folder so only fresh results live there
  if (dir.exists(param_dir)) unlink(param_dir, recursive = TRUE, force = TRUE)
  dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)
  
  # one subfolder per replicate: 001, 002, ...
  run_dirs <- file.path(param_dir, formatC(1:nRep, width = 3, format = "d", flag = "0"))
  invisible(lapply(run_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  
  # ===== stochastic engine + run =====
  setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)
  
  MGDrivESim <- Network$new(
    params = netPar,
    driveCube = cube,
    patchReleases = patchReleases,
    migrationMale = moveMat,
    migrationFemale = moveMat,
    migrationBatch = batchMigration,
    directory = run_dirs,         # <- each replicate writes into its own subfolder
    verbose = TRUE,
    toxMort = toxCurve,
    toxInt  = toxTime
  )
  MGDrivESim$multRun(verbose = TRUE)
  
  # ===== post-process: split & aggregate per replicate =====
  for (i in seq_len(nRep)) {
    rd <- run_dirs[i]
    # This generates "M_Run%03d_Patch001.csv"
    splitOutput(readDir = rd, remFile = TRUE, verbose = FALSE)
    # This generates "F_Aggregate_Run%03d_Patch001.csv"
    aggregateFemales(readDir = rd, genotypes = cube$genotypesID, remFile = TRUE, verbose = TRUE)
    
    # quick sanity: ensure exactly the two CSVs exist for Patch001
    m_ok <- length(list.files(rd, pattern = "^M_Run\\d+_Patch001\\.csv$", full.names = TRUE)) >= 1
    f_ok <- length(list.files(rd, pattern = "^F_Aggregate_Run\\d+_Patch001\\.csv$", full.names = TRUE)) >= 1
    if (!m_ok || !f_ok) {
      warning(sprintf("Missing expected CSV(s) in %s", rd))
    }
  }
  
  cat("Finished:", tag, "\n")
  invisible(tag)
}


####################
# Build the parameter grid
####################
param_grid <- expand.grid(
  rel = rel_grid,
  fy = fy_grid,
  pq = pq_grid,
  fl = fl_grid,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

cat("Total parameter sets:", nrow(param_grid), "\n")

####################
# Run (parallel if possible)
####################
is_windows <- identical(.Platform$OS.type, "windows")
cores <- if (is_windows) 1 else max(1L, parallel::detectCores() - 4L)

cat("Using cores:", cores, "\n")

if (cores > 1) {
  invisible(pbmclapply(
    X = seq_len(nrow(param_grid)),
    FUN = function(ii) {
      with(param_grid[ii, ], run_param_set(rel, fy, pq, fl))
    },
    mc.cores = cores
  ))
} else {
  invisible(pblapply(
    X = seq_len(nrow(param_grid)),
    FUN = function(ii) {
      with(param_grid[ii, ], run_param_set(rel, fy, pq, fl))
    }
  ))
}


cat("All simulations completed. Output root:", normalizePath(outFolder), "\n")
