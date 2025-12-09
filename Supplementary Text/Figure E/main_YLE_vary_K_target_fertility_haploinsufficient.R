####################
# Load libraries
####################
rm(list = ls()); gc(); save.image()
suppressPackageStartupMessages({
  library(MGDrivE)
  library(BBmisc)
  library(parallel)
  library(pbapply)
  library(pbmcapply)
})

####################
# Top-level output
####################
outFolder <- "mgdriveYLE_sweep_5yr_release_target_viablity_byK"
if (dir.exists(outFolder)) {
  unlink(outFolder, recursive = TRUE, force = TRUE)
  cat("Folder deleted:", outFolder, "\n")
}
dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)
cat("New folder created:", outFolder, "\n")

############################
# Simulation horizon & reps
############################
tMax   <- 58 * 365   # days
nRep   <- 100         # Monte Carlo iterations per K
tstart <- 8 * 365    # release start (days since 0)

####################
# Biology & demography (fixed)
####################
bioParameters <- list(
  betaK   = 6.0,
  litters = 7.5,
  tGest   = 19,
  tNursing= 23,
  tAdo    = 37,
  muAd    = 1/690,
  theta   = 22.4
)

sitesNumber <- 1
moveMat     <- matrix(1, 1, 1)  # no dispersal
batchMigration <- basicBatchMigration(
  batchProbs = 0,
  sexProbs  = c(.5, .5),
  numPatches= sitesNumber
)

####################
# K values to explore
####################
K_grid <- c(1e3, 5e3, 1e4, 5e4, 1e5, 2e5)

####################
# Fixed genetic / release parameters
####################
rel_fixed <- 0.08   # 8% of adults per month (example; change as needed)
pq <- 0
fl    <- 0
fs    <- 1
fy <- 0.8
mu_val <- 0.97
j_val  <- 0.25
c_val  <- 0.93

####################
# Source inheritance cube
####################
source("generate_YLE_inheritance_cube.R")

####################
# Toxicology
####################
exposure  <- 1.26
expTime   <- 7
toxCurve  <- (1/(1+exp(-3.2*(log(exposure)-log(1.774)))))/expTime
toxTime   <- c(10^6)  # effectively absent

####################
# Helpers
####################
K_tag <- function(K) {
  if (K %% 1000 == 0) {
    paste0("K", K/1000, "k")  # 1000 -> K1k, 50000 -> K50k
  } else {
    paste0("K", K)
  }
}

run_for_K <- function(K_val) {
  cat("Starting K =", K_val, "\n")
  
  # ===== inheritance cube for fixed parameter set =====
  source("generate_YLE_inheritance_cube.R")
  alleleFitness_vec <- c(W = 1, H = 1, R = 1, Y = 1, X = 1, y = fy)
  cube <<- generate_YLE_inheritance_cube(
    c = c_val, j = j_val, mu = mu_val, p = pq, q = pq,
    female_lethality = fl,
    female_sterility = fs,
    Haploinsufficient = TRUE,
    muAd = 1/690,
    alleleFitness = alleleFitness_vec
  )
  
  # ===== releases (adult males; rel_fixed of total adult pop × 0.5) =====
  releasesParameters_SP <- list(
    releasesStart     = tstart,
    releasesNumber    = round(12*10),       # 10 years monthly
    releasesInterval  = 30,
    releaseProportion = round(K_val * rel_fixed * 0.5)
  )
  
  ReleasesVector_SP <- generateReleaseVector(
    driveCube          = cube,
    releasesParameters = releasesParameters_SP
  )
  
  patchReleases <- replicate(
    n      = sitesNumber,
    expr   = list(
      maleReleases         = NULL,
      femaleReleases       = NULL,
      eggReleases          = NULL,
      matedFemaleReleases  = NULL
    ),
    simplify = FALSE
  )
  patchReleases[[1]]$maleReleases <- ReleasesVector_SP
  
  # ===== parameterized network =====
  netPar <- parameterizeMGDrivE(
    runID       = 1,
    simTime     = tMax,
    nPatch      = sitesNumber,
    beta        = bioParameters$betaK,
    litters     = bioParameters$litters,
    tGest       = bioParameters$tGest,
    tNursing    = bioParameters$tNursing,
    tAdo        = bioParameters$tAdo,
    k           = K_val,
    theta       = bioParameters$theta,
    inheritanceCube = cube,
    AdPopRatio_M = matrix(c(1, 0), 1, 2,
                          dimnames = list(NULL, c("mYWW", "fXWW"))),
    AdPopRatio_F = matrix(c(0, 1), 1, 2,
                          dimnames = list(NULL, c("mYWW", "fXWW")))
  )
  
  # ===== folders: outFolder/Kxxk/rep =====
  tag      <- K_tag(K_val)
  param_dir <- file.path(outFolder, tag)
  
  if (dir.exists(param_dir))
    unlink(param_dir, recursive = TRUE, force = TRUE)
  dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)
  
  run_dirs <- file.path(
    param_dir,
    formatC(1:nRep, width = 3, format = "d", flag = "0")
  )
  invisible(lapply(run_dirs, dir.create, recursive = TRUE,
                   showWarnings = FALSE))
  
  # ===== stochastic engine + run =====
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
  
  # ===== post-process: split & aggregate per replicate =====
  for (i in seq_len(nRep)) {
    rd <- run_dirs[i]
    
    splitOutput(
      readDir = rd,
      remFile = TRUE,
      verbose = FALSE
    )
    
    aggregateFemales(
      readDir    = rd,
      genotypes  = cube$genotypesID,
      remFile    = TRUE,
      verbose    = TRUE
    )
    
    m_ok <- length(list.files(
      rd,
      pattern = "^M_Run\\d+_Patch001\\.csv$",
      full.names = TRUE
    )) >= 1
    f_ok <- length(list.files(
      rd,
      pattern = "^F_Aggregate_Run\\d+_Patch001\\.csv$",
      full.names = TRUE
    )) >= 1
    
    if (!m_ok || !f_ok) {
      warning(sprintf("Missing expected CSV(s) in %s", rd))
    }
  }
  
  cat("Finished K =", K_val, " -> ", tag, "\n")
  invisible(tag)
}

####################
# Run (parallel if possible) over K only
####################
is_windows <- identical(.Platform$OS.type, "windows")
cores <- if (is_windows) 1 else max(1L, 6) # parallel::detectCores()/2

cat("Using cores:", cores, "\n")
cat("K values:", paste(K_grid, collapse = ", "), "\n")

if (cores > 1) {
  invisible(pbmclapply(
    X = K_grid,
    FUN = run_for_K,
    mc.cores = cores
  ))
} else {
  invisible(pblapply(
    X = K_grid,
    FUN = run_for_K
  ))
}

cat("All simulations completed. Output root:",
    normalizePath(outFolder), "\n")
