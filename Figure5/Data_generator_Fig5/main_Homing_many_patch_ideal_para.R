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
outFolder   <- "mgdrive_Homing_5patch_250_mig_1daily"
tMax        <- 58 * 365              # simulation time (days)
nRep        <- 100                   # number of Monte Carlo iterations
sitesNumber <- 5                     # <<< choose any number of patches
Neq         <- 10000                 # equilibrium per patch (if uniform)

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
# Inheritance cube (YLE)
####################
source("generate_Homing_inheritance_cube.R")

dayOmega <- calcOmega(mu = bioParameters$muAd, lifespanReduction = 1)
alleleFitness_vec <- c(W = 1, H = 1, R = 1, N = 1)

cube <- generate_Homing_inheritance_cube(
  pC = 1,           # Probability of successful cut
  pN = 0,           # Probability of NHEJ (given cut occurs)
  pL = 0,            # Probability of loss-of-function after NHEJ
  alleleFitness = alleleFitness_vec # Optional per-allele fitness costs
)
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
  releaseProportion = round(cap[1] * 0.05 * 0.5)  # 5% of K, half male
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

plotMGDrivESingle(
  readDir  = folderNames[1],
  totalPop = TRUE,
  lwd      = 1,
  alpha    = .8
)

plotMGDrivEMult(
  readDir = outFolder,
  lwd     = 0.35,
  alpha   = 0.70
)
