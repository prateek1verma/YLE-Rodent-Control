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
outFolder   <- "mgdrive_Xshredder_nonideal_5patch_250_mig_1daily"
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
# Source the function (this loads the function definition)
# source("generate_XShredder_inheritance_cube.R")
source("cube_XShredderY.R")

# Non Ideal cube
cube <- cube_XShredderY(cX = 0.97, crX = 0.25)
# Cas-9 carriers: mXA + mXB + mRA + mRB
fertility_cost <- 0.2
cube$s[c("mXA","mXB","mRA","mRB")] <- 1-fertility_cost

# ## Ideal cube
# cube <- cube_XShredderY(cX = 1, crX = 0)
# # Cas-9 carriers: mXA + mXB + mRA + mRB
# fertility_cost <- 0
# cube$s[c("mXA","mXB","mRA","mRB")] <- 1-fertility_cost
# genotypes <- cube$genotypesID

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
  mXY = rep(1, sitesNumber),
  fXX = rep(0, sitesNumber)
)

AdPopRatio_F <- cbind(
  mXY = rep(0, sitesNumber),
  fXX = rep(1, sitesNumber)
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
