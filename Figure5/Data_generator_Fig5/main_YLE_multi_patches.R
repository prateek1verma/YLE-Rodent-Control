####################
# Load libraries
####################
rm(list = ls())
gc()
save.image(file = paste0("workspace_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData"))  # overwrites .RData with empty workspace

library(MGDrivE)
library(BBmisc)

####################
# User-controlled settings
####################
outFolder   <- "mgdriveYLE_rel_5per_250_5patches_mig_5daily_fy0p5"
tMax        <- 58 * 365              # simulation time (days)
nRep        <- 100                    # number of Monte Carlo iterations
sitesNumber <- 5                     # <<< choose any number of patches
Neq         <- 10000                 # equilibrium per patch (if uniform)
daily_migration <- 5
mig_rate <- daily_migration/Neq
rel <- 0.05
fy <- 0.5
Tperiod <- 10

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
  upper = rep.int(mig_rate, times = sitesNumber - 1),
  lower = rep.int(mig_rate, times = sitesNumber - 1)
)

batchMigration <- basicBatchMigration(
  batchProbs = 0,
  sexProbs   = c(.5, .5),
  numPatches = sitesNumber
)

####################
# Inheritance cube (YLE)
####################
source("generate_YLE_inheritance_cube.R")
pq <- 0
fl    <- 0
fs    <- 1
mu_val <- 0.97
j_val  <- 0.25
c_val  <- 0.93
alleleFitness_vec <- c(W = 1, H = 1, R = 1, Y = 1, X = 1, y = fy)
cube <- generate_YLE_inheritance_cube(
  c = c_val, j = j_val, mu = mu_val, p = pq, q = pq,
  female_lethality = fl,
  female_sterility = fs,
  Haploinsufficient = TRUE,
  muAd = 1/690,
  alleleFitness = alleleFitness_vec
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
  releasesNumber   = round(12 * Tperiod),
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
