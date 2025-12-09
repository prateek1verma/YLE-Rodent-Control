####################
# Load libraries
####################
rm(list = ls())
gc()
save.image()  # overwrites .RData with empty workspace

library(MGDrivE)
library(BBmisc)

####################
# User-controlled settings
####################
# "mgdrive_Homing_5patch_250_mig_5daily" result by mistake in the folder here
outFolder   <- "mgdrive_Homing_non_ideal_5patch_250_mig_5daily"
tMax        <- 58 * 365              # simulation time (days)
nRep        <- 100                   # number of Monte Carlo iterations
sitesNumber <- 5                     # <<< choose any number of patches
Neq         <- 10000                 # equilibrium per patch (if uniform)
rel <- 0.05

# carrying capacities per patch (length must equal sitesNumber)
cap <- rep(Neq, sitesNumber)
tstart <- 8 * 365                    # start of releases (days)

####################
# Output folder
####################
if (dir.exists(outFolder)) {
  unlink(outFolder, recursive = TRUE, force = TRUE)
  cat("Folder deleted:", outFolder, "\n")
} else {
  cat("Folder does not exist:", outFolder, "\n")
}

dir.create(outFolder)
cat("New folder created:", outFolder, "\n")

folderNames <- file.path(
  outFolder,
  formatC(x = 1:nRep, width = 3, format = "d", flag = "0")
)

####################
# Biological parameters
####################
bioParameters <- list(
  betaK   = 6,
  litters = 7.5,
  tGest   = 19,
  tNursing= 23,
  tAdo    = 37,
  muAd    = 1 / 690,
  theta   = 22.4
)

####################
# Helper: tri-diagonal movement matrix for line of patches
####################
triDiag <- function(upper, lower) {
  n <- length(upper) + 1
  retMat <- matrix(0, nrow = n, ncol = n)
  idx    <- 1:length(upper)
  retMat[cbind(idx + 1, idx)] <- lower
  retMat[cbind(idx, idx + 1)] <- upper
  diag(retMat) <- 1 - rowSums(retMat)
  retMat
}

####################
# Movement and batch migration
####################
moveMat <- triDiag(
  upper = rep.int(0.0001, times = sitesNumber - 1),
  lower = rep.int(0.0001, times = sitesNumber - 1)
)

batchMigration <- basicBatchMigration(
  batchProbs = 0,
  sexProbs   = c(.5, .5),
  numPatches = sitesNumber
)

####################
# Inheritance cube (Homing)
####################
source("generate_Homing_inheritance_cube.R")
# Non - ideal cube specification, Now call the function
alleleFitness_vec <- c(W = 1, H = 0.8, R = 1, N = 0.5)
cube <- generate_Homing_inheritance_cube(
  pC = 0.97,   # 95% cutting efficiency
  pN = 0.25,   # 25% NHEJ when cut occurs
  pL = 5/6,     # 500/6 % of NHEJ events cause loss-of-function
  alleleFitness = alleleFitness_vec  # Drive and R have costs
)

omega_calculated <- sapply(cube$omega, function(l) {
  calcOmega(mu = bioParameters$muAd, lifespanReduction = l)
})
cube$omega <- omega_calculated

genotypes <- cube$genotypesID

####################
# Releases
####################
# Empty release list, one element per patch
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

# Release parameters (here: same schedule, amount based on patch 1 capacity)
releasesParameters_SP <- list(
  releasesStart    = tstart,
  releasesNumber   = round(1),
  releasesInterval = 30,
  releaseProportion = round(cap[1] * rel * 0.5)  # 5% of K, half male
)

# Toxicology
exposure  <- 1.26
expTime   <- 7
toxCurve  <- (1 / (1 + exp(-3.2 * (log(exposure) - log(1.774))))) / expTime
toxTime   <- c(1e6)

# Release vector for a single patch
ReleasesVector_SP <- generateReleaseVector(
  driveCube          = cube,
  releasesParameters = releasesParameters_SP
)

# Choose which patches receive releases (here: only patch 1)
release_patches <- 1L
for (p in release_patches) {
  patchReleases[[p]]$maleReleases <- ReleasesVector_SP
  # patchReleases[[p]]$femaleReleases <- ReleasesVector_SP  # if needed
}

####################
# Initial adult population ratios for ANY number of patches
####################
# Each row = patch; columns = initial adult genotypes
# Here: all patches start with mYWW males and fXWW females
AdPopRatio_M <- cbind(
  mYWW = rep(1, sitesNumber),
  fXWW = rep(0, sitesNumber)
)

AdPopRatio_F <- cbind(
  mYWW = rep(0, sitesNumber),
  fXWW = rep(1, sitesNumber)
)

####################
# Combine parameters and run
####################
netPar <- parameterizeMGDrivE(
  runID          = 1,
  simTime        = tMax,
  nPatch         = sitesNumber,
  beta           = bioParameters$betaK,
  litters        = bioParameters$litters,
  tGest          = bioParameters$tGest,
  tNursing       = bioParameters$tNursing,
  tAdo           = bioParameters$tAdo,
  k              = cap,
  theta          = bioParameters$theta,
  inheritanceCube= cube,
  AdPopRatio_M   = AdPopRatio_M,
  AdPopRatio_F   = AdPopRatio_F
)

setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)

MGDrivESim <- Network$new(
  params          = netPar,
  driveCube       = cube,
  patchReleases   = patchReleases,
  migrationMale   = moveMat,
  migrationFemale = moveMat,
  migrationBatch  = batchMigration,
  directory       = folderNames,
  verbose         = TRUE,
  toxMort         = toxCurve,
  toxInt          = toxTime
)

MGDrivESim$multRun(verbose = TRUE)

####################
# Post-analysis
####################
for (i in 1:nRep) {
  splitOutput(
    readDir = folderNames[i],
    remFile = TRUE,
    verbose = FALSE
  )
  aggregateFemales(
    readDir   = folderNames[i],
    genotypes = cube$genotypesID,
    remFile   = TRUE,
    verbose   = TRUE
  )
}

