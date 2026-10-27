###############################################################################
# Early-Stage Female-Specific RIDL (fsRIDL) Inheritance Cube Generator
# Model: Sex-explicit 6-genotype framework (mY... / fX...)
# Lethality Stage: Embryonic / Egg Stage (via Viability Tensor 'tau')
###############################################################################

generate_early_fsRIDL_cube <- function(female_lethality = 1.0, alleleFitness = NULL) {
  # Genotype state space (3 male genotypes, 3 female genotypes)
  male_genotypes   <- c("mYWW", "mYHW", "mYHH")
  female_genotypes <- c("fXWW", "fXHW", "fXHH")
  genotypes        <- c(male_genotypes, female_genotypes)
  nG               <- length(genotypes)
  
  # --- Internal Helper Functions ---
  split_geno <- function(g) {
    core <- substring(g, 3, 4) # Extract "WW", "HW", or "HH"
    c(substr(core, 1, 1), substr(core, 2, 2))
  }
  
  gametes <- function(g) {
    a <- split_geno(g)
    if (a[1] == a[2]) setNames(1, a[1]) else setNames(c(0.5, 0.5), sort(a))
  }
  
  # --- 1. Inheritance Tensor (ih): 6 x 6 x 6 ---
  cube_ih <- array(0, dim = c(nG, nG, nG), dimnames = list(genotypes, genotypes, genotypes))
  
  for (f in female_genotypes) {
    gf <- gametes(f)
    for (m in male_genotypes) {
      gm <- gametes(m)
      for (aF in names(gf)) {
        for (aM in names(gm)) {
          son <- paste0("mY", paste0(sort(c(aF, aM)), collapse = ""))
          dau <- paste0("fX", paste0(sort(c(aF, aM)), collapse = ""))
          
          # Mendelian gametic combination with 1:1 primary sex ratio
          cube_ih[f, m, son] <- cube_ih[f, m, son] + 0.5 * gf[aF] * gm[aM]
          cube_ih[f, m, dau] <- cube_ih[f, m, dau] + 0.5 * gf[aF] * gm[aM]
        }
      }
    }
  }
  
  # --- 2. Viability Mask (tau): Egg-Stage Embryonic Lethality ---
  cube_tau <- array(1, dim = c(nG, nG, nG), dimnames = list(genotypes, genotypes, genotypes))
  
  # Set egg survival for transgenic females (fXHW and fXHH)
  surv_female_H <- 1 - female_lethality
  cube_tau[, , "fXHW"] <- surv_female_H
  cube_tau[, , "fXHH"] <- surv_female_H
  
  # --- 3. Fecundity / Egg Laying Mask (phi) ---
  phi <- setNames(c(rep(0, 3), rep(1, 3)), genotypes)
  
  # --- 4. Adult Fitness / Longevity (omega) ---
  omega <- setNames(rep(1.0, nG), genotypes)
  if (!is.null(alleleFitness)) {
    omega_W <- if (!is.na(alleleFitness["W"])) alleleFitness["W"] else 1.0
    omega_H <- if (!is.na(alleleFitness["H"])) alleleFitness["H"] else 1.0
    for (g in genotypes) {
      alleles <- split_geno(g)
      omega[g] <- if ("H" %in% alleles) omega_H else omega_W
    }
  }
  
  # --- 5. Mating Competitiveness (eta) ---
  eta <- matrix(1.0, nrow = nG, ncol = nG, dimnames = list(genotypes, genotypes))
  
  # --- 6. Post-Hatch Survival (xiF, xiM) ---
  xiF <- setNames(rep(1.0, nG), genotypes)
  xiM <- setNames(rep(1.0, nG), genotypes)
  
  return(list(
    ih          = cube_ih,
    tau         = cube_tau,
    genotypesID = genotypes,
    genotypesN  = nG,
    wildType    = "fXWW",   # Scalar string required by MGDrivE internals
    releaseType = "mYHH",   # Default release genotype (Homozygous Males)
    eta         = eta,
    phi         = phi,
    omega       = omega,
    xiF         = xiF,
    xiM         = xiM,
    s           = setNames(c(1, 1, 1, 1, surv_female_H, surv_female_H), genotypes)
  ))
}