generate_SIT_inheritance_cube <- function(sterility = 1) {
  
  # Genotype sets (diploid locus; sex denoted by prefix mY / fX)
  male_genotypes   <- c("mY", "my")
  female_genotypes <- c("fX")
  genotypes <- c(male_genotypes, female_genotypes)
  
  # --- inheritance cube (ih): parents→offspring genotype probabilities ---
  cube_ih <- array(0, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                   dimnames = list(genotypes, genotypes, genotypes))

  cube_ih["fX", "mY", "mY"] <- 0.5
  cube_ih["fX", "mY", "fX"] <-  0.5
  cube_ih["fX", "my", "mY"] <- 0.5 
  cube_ih["fX", "my", "fX"] <-  0.5

  
  # --- viability cube (tau): here genotype-only viability (same for all parent pairs) ---
  cube_tau <- array(0, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                    dimnames = list(genotypes, genotypes, genotypes))
  
  # Female-specific lethality for any H-bearing genotype
  # --- viability cube (tau): pgSIT is not a zygotic-lethal system in the field ---
  # Females carrying both components are culled in the factory (not released),
  # so in-field zygotes don't have genotype-based lethal effects here.
  cube_tau <- array(1, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                    dimnames = list(genotypes, genotypes, genotypes))
  
  # --- fertility mask (phi): 0 for males, 1 for females (as in your original) ---
  phi <- c(rep(0, length(male_genotypes)), rep(1, length(female_genotypes)))
  names(phi) <- genotypes
  
  # --- genotype-level fitness (omega): multiplicative from alleleFitness (optional) ---
  omega <- setNames(rep(1, length(genotypes)), genotypes)
  
  # --- additional vectors to mirror your structure ---
  xiF <- setNames(rep(1, length(genotypes)), genotypes)
  xiM <- setNames(rep(1, length(genotypes)), genotypes)
  eta <- matrix(1, nrow = length(genotypes), ncol = length(genotypes),
                dimnames = list(genotypes, genotypes))
  
  # Per-genotype survival shortcut (s): mirror tau’s genotype effect for convenience
  s <- setNames(rep(1, length(genotypes)), genotypes)
  s["my"] <- 1 - sterility   # e.g., sterility=1 -> s["my"]=0 (no siring)
  
  list(
    ih = cube_ih,
    tau = cube_tau,
    genotypesID = genotypes,
    genotypesN = length(genotypes),
    wildType = c("mY", "fX"),
    eta = eta,
    phi = phi,
    omega = omega,
    xiF = xiF,
    xiM = xiM,
    s = s,
    releaseType = "my"   # standard SIT release: heterozygous males
  )
}

# example
# cube_SIT <- generate_SIT_inheritance_cube(sterility = 1) 
