####################
# Load libraries
####################
rm(list = ls())
gc()
save.image()  # overwrites .RData with empty workspace

library(MGDrivEmouse) # I have renamed the "MouseGD - MGDrivE" R package/library by Brown et al to "MGDrivEmouse"
# library(MGDrivE) # Use this if you have installed MouseGd from Github and have not not renamed the R package.
library(BBmisc)  # Library for plotting

####################
# Output Folder
####################
# Assign the folder name to a variable
outFolder <- "mgdriveYLE_viablity_100runs_5yr_50rel"

# Check if the folder exists before deleting
if (dir.exists(outFolder)) {
  unlink(outFolder, recursive = TRUE, force = TRUE)
  cat("Folder deleted:", outFolder, "\n")
} else {
  cat("Folder does not exist:", outFolder, "\n")
}

# Create a new, empty folder
dir.create(outFolder)

# Print a message to confirm the action
cat("New folder created:", outFolder, "\n")
####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 58*365

# number of Monte Carlo iterations
nRep <- 2 # change to 100 runs to replicate results of figure 1

# each Monte Carlo iteration gets its own folder
folderNames <- file.path(outFolder,
                         formatC(x = 1:nRep, width = 3, format = "d", flag = "0"))

# biological parameters
bioParameters <- list(betaK= 6, litters = 7.5, tGest=19, tNursing=23, tAdo=37, muAd = (1/690), theta = 22.4)

sitesNumber <- 1 # number of patches
Neq <- 10000
cap <- Neq # adult carrying capacity for small patch and big patch respectively
tstart <- 8*365

# a 1-node network where mosquitoes do not leave
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)
sitesNumber <- nrow(moveMat)

# batch migration is disabled by setting the probability to 0
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5),
                                      numPatches=sitesNumber)

####################
# Basic Inheritance pattern
####################
# Source the function (this loads the function definition)
source("generate_YLE_inheritance_cube.R")

# Now call the function
dayOmega <- calcOmega(mu = bioParameters$muAd, lifespanReduction = 1) # lifespan reduction of 0%
alleleFitness_vec <- c(W = 1, H = 1, R = 1, Y=1, X=1, y=dayOmega)
cube <- generate_YLE_inheritance_cube(c = 1, j = 0, mu = 1, p = 0, q = 0., female_lethality = 1, alleleFitness = alleleFitness_vec)
# cube <- generate_YLE_inheritance_cube(c = 0.93, j = 0.043, mu = 1, p = 0.1, q = 0.1, alleleFitness = alleleFitness_vec)
# c = 0.93, j = 0.043, mu = 1, p = 0, q = 0,
genotypes <- cube$genotypesID

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
releasesParameters_SP <- list(releasesStart=tstart,
                              releasesNumber=round(12*5),
                              releasesInterval=30,
                              releaseProportion=round(cap*0.01*0.5))

# specify toxicological parameters
exposure  = 1.26 # estimated average dose of toxicant in mice
expTime = 7 # length of toxicity experiment used to generate exposure-response curve

# daily average toxicant-induced mortality in mice
# calculated from exposure-response equation and the
# length of the corresponding toxicity test:

toxCurve = (1/(1+exp(-3.2*(log(exposure)-log(1.774)))))/expTime 

# the time interval at which the toxicant is
# present in the system, assuming no chemical degradation,
# assuming chemical CAN be removed:

# toxTime = c((releasesParameters_SP$releasesStart
#              #+ releasesParameters_BP$releasesNumber*releasesParameters_BP$releasesInterval
# ),(releasesParameters_SP$releasesStart +
#      #releasesParameters_BP$releasesNumber*releasesParameters_BP$releasesInterval
#      + 56)) 
toxTime = c(10^6)

# generate small patch release vector
ReleasesVector_SP <- generateReleaseVector(driveCube=cube,
                                           releasesParameters=releasesParameters_SP)

# put releases into the proper place in the release list
#  This specifies the releases for the small patch (SP) and big patch (BP) respectively 
patchReleases[[1]]$maleReleases <- ReleasesVector_SP
# patchReleases[[1]]$femaleReleases <- ReleasesVector_SP



####################
# Combine parameters and run!
####################
# setup parameters for the network. This builds a list of parameters required for
# every population in the network.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                              beta=bioParameters$betaK, litters = bioParameters$litters,
                              tGest=bioParameters$tGest, 
                              tNursing=bioParameters$tNursing, tAdo=bioParameters$tAdo,
                              k = cap, theta = bioParameters$theta, inheritanceCube = cube,
                              AdPopRatio_M = matrix(c(1,0),1,2, dimnames = list(NULL,c("mYWW","fXWW"))),
                              AdPopRatio_F = matrix(c(0,1),1,2, dimnames = list(NULL,c("mYWW","fXWW"))))


# set MGDrivE to run stochastic
setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=folderNames,
                          verbose = TRUE,
                          toxMort = toxCurve,
                          toxInt = toxTime)
# run simulation
MGDrivESim$multRun(verbose = TRUE)


####################
# Post Analysis
####################
# First, split output by patch
# Second, aggregate females by their mate choice
for(i in 1:nRep){
  splitOutput(readDir = folderNames[i], remFile = TRUE, verbose = FALSE)
  aggregateFemales(readDir = folderNames[i], genotypes = cube$genotypesID,
                   remFile = TRUE, verbose = TRUE)
}


# plot output of first run to see effect
plotMGDrivESingle(readDir=folderNames[1],totalPop = TRUE,lwd=1,alpha=.8)

# # plot all repetitions together
# png(paste(outFolder, ".png", sep = ""), width = 1000, height = 750)
# plotMGDrivEMult(readDir=outFolder,lwd=0.35,alpha=0.75)
# dev.off()


# plot all nRep repetitions together
# plotMGDrivEMult(readDir = outFolder, lwd = 0.35, alpha = 0.70)

# --- USAGE EXAMPLE ---

# # outFolder <- "mgdriveYLE_viablity_100runs" #"mgdriveYLE_viablity_20runs_2yr"
# base_dir <- file.path(getwd(), outFolder)
# source("gene_drive_plots_multi.R")
# # Auto-detect all runs in base_dir
# make_patch_figure_multi(
#   base_dir = base_dir,
#   patch_id = "001",
#   tstart   = round(tstart),
#   tEnd     = 10,         # show first 6 years after tstart
#   y_log    = FALSE,
#   Neq      = 0.5*Neq
# )
# 
# 
# source("gene_drive_plots_multi_tot.R")
# make_patch_figure_multi_tot(
#   base_dir = base_dir,
#   patch_id = "001",
#   outfile  = NULL,
#   tstart   = round(tstart),
#   tEnd     = 10,         # show first 6 years after tstart
#   y_log    = FALSE
# )
# 
# 
# source("gene_drive_plots_multi_new.R")
# make_patch_figure_multi(
#   base_dir = base_dir,
#   patch_id = "001",
#   tstart   = round(tstart),
#   tEnd     = 20,
#   center   = "mean",  # or "mean"
#   groups_to_plot = c("WT female","YLE male")
# )
