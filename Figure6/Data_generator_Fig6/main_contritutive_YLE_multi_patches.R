####################
# Load libraries
####################
rm(list = ls())
gc()

library(MGDrivEmouse2)
# library(MGDrivE) 
library(BBmisc)
library(parallel)

####################
# User settings
####################
outFolder <- "mgdriveboostedYLE_resistance_rel_1per_5patches_mig_2daily_fy0p8_xshred_0p80"

tMax        <- 58 * 365
nRep        <- 6
sitesNumber <- 5
Neq         <- 10000

daily_migration <- 2
mig_rate <- daily_migration / Neq

rel          <- 0.01
fy           <- 0.8
x_shred_val  <- 0.8
Tperiod      <- 10
tstart       <- 8 * 365

cap <- rep(Neq, sitesNumber)

folderNames <- file.path(
  outFolder,
  formatC(seq_len(nRep), width = 3, format = "d", flag = "0")
)

####################
# Output folders
####################
if (dir.exists(outFolder)) {
  unlink(outFolder, recursive = TRUE, force = TRUE)
}

dir.create(outFolder, showWarnings = FALSE)
invisible(lapply(folderNames, dir.create, recursive = TRUE, showWarnings = FALSE))

####################
# Biological parameters
####################
bioParameters <- list(
  betaK    = 6,
  litters  = 7.5,
  tGest    = 19,
  tNursing = 23,
  tAdo     = 37,
  muAd     = 1 / 690,
  theta    = 22.4
)

####################
# Movement matrix
####################
make_tri_diag <- function(n, migration_rate) {
  if (n == 1) {
    return(matrix(1, nrow = 1, ncol = 1))
  }
  
  mat <- matrix(0, nrow = n, ncol = n)
  
  for (i in seq_len(n - 1)) {
    mat[i, i + 1] <- migration_rate
    mat[i + 1, i] <- migration_rate
  }
  
  diag(mat) <- 1 - rowSums(mat)
  mat
}

moveMat <- make_tri_diag(sitesNumber, mig_rate)

batchMigration <- basicBatchMigration(
  batchProbs = 0,
  sexProbs   = c(0.5, 0.5),
  numPatches = sitesNumber
)

####################
# Inheritance cube
####################
source("generate_boosted_YLE_constitutive_inheritance_cube.R")

cube <- generate_boosted_YLE_constitutive_inheritance_cube(
  c = 0.93,
  j = 0.25,
  mu = 0.97,
  p = 0,
  q = 0,
  female_lethality = 0,
  female_sterility = 1,
  Haploinsufficient = TRUE,
  muAd = bioParameters$muAd,
  alleleFitness = c(W = 1, H = 1, R = 1, Y = 1, X = 1, y = fy, A = 0.8, B = 1),
  x_shred = x_shred_val
)

genotypes <- cube$genotypesID

stopifnot(all(c("mYWWBB", "fXWWBB") %in% genotypes))

####################
# Releases
####################
releaseParameters <- list(
  releasesStart     = tstart,
  releasesNumber    = round(12 * Tperiod),
  releasesInterval  = 30,
  releaseProportion = round(cap[1] * rel * 0.5)
)

releaseVector <- generateReleaseVector(
  driveCube = cube,
  releasesParameters = releaseParameters
)

patchReleases <- replicate(
  n = sitesNumber,
  expr = list(
    maleReleases        = NULL,
    femaleReleases      = NULL,
    eggReleases         = NULL,
    matedFemaleReleases = NULL
  ),
  simplify = FALSE
)

release_patches <- 1L

for (p in release_patches) {
  patchReleases[[p]]$maleReleases <- releaseVector
}

####################
# Initial adult population
####################
AdPopRatio_M <- cbind(
  mYWWBB = rep(1, sitesNumber),
  fXWWBB = rep(0, sitesNumber)
)

AdPopRatio_F <- cbind(
  mYWWBB = rep(0, sitesNumber),
  fXWWBB = rep(1, sitesNumber)
)

####################
# Toxicology
####################
exposure <- 1.26
expTime  <- 7

toxCurve <- (1 / (1 + exp(-3.2 * (log(exposure) - log(1.774))))) / expTime
toxTime  <- c(1e6)

####################
# One replicate: simulate + post-process
####################
run_one_rep <- function(i) {
  
  setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)
  
  netPar_i <- parameterizeMGDrivE(
    runID           = i,
    simTime         = tMax,
    nPatch          = sitesNumber,
    beta            = bioParameters$betaK,
    litters         = bioParameters$litters,
    tGest           = bioParameters$tGest,
    tNursing        = bioParameters$tNursing,
    tAdo            = bioParameters$tAdo,
    k               = cap,
    theta           = bioParameters$theta,
    inheritanceCube = cube,
    AdPopRatio_M    = AdPopRatio_M,
    AdPopRatio_F    = AdPopRatio_F
  )
  
  sim_i <- Network$new(
    params          = netPar_i,
    driveCube       = cube,
    patchReleases   = patchReleases,
    migrationMale   = moveMat,
    migrationFemale = moveMat,
    migrationBatch  = batchMigration,
    directory       = folderNames[i],
    verbose         = FALSE,
    toxMort         = toxCurve,
    toxInt          = toxTime
  )
  
  sim_i$oneRun(verbose = FALSE)
  
  splitOutput(
    readDir = folderNames[i],
    remFile = TRUE,
    verbose = FALSE
  )
  
  aggregateFemales(
    readDir   = folderNames[i],
    genotypes = genotypes,
    remFile   = TRUE,
    verbose   = FALSE
  )
  
  return(folderNames[i])
}

####################
# Run replicates in parallel
####################
nCores <- max(1, detectCores() - 4)

results <- mclapply(
  X = seq_len(nRep),
  FUN = run_one_rep,
  mc.cores = min(nCores, nRep)
)

print(unlist(results))