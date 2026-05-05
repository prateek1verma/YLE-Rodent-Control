generate_fRIDL_inheritance_cube <- function(female_lethality = 1, alleleFitness = NULL) {
  # female_lethality: 1 = all H-bearing females die (field); 0 = full rescue (lab)
  if (is.null(alleleFitness)) {
    alleleFitness <- c(W = 1, H = 1, Y=1, X=1)  # per-allele multiplicative fitness (can customize)
  }
  
  # Genotype sets (diploid locus; sex denoted by prefix mY / fX)
  alleleSet <- c("W", "H")
  male_genotypes   <- c("mYWW", "mYHW", "mYHH")
  female_genotypes <- c("fXWW", "fXHW", "fXHH")
  genotypes <- c(male_genotypes, female_genotypes)
  
  # --- helpers ---
  split_geno_to_alleles <- function(geno) {
    core <- substr(geno, 3, 5)         # e.g., "WW", "HW", "HH"
    a1 <- substr(core, 1, 1); a2 <- substr(core, 2, 2)
    if (!(a1 %in% alleleSet && a2 %in% alleleSet)) stop("Invalid genotype: ", geno)
    c(a1, a2)
  }
  get_sex_prefix <- function(geno) {
    if (startsWith(geno, "mY")) return("mY")
    if (startsWith(geno, "fX")) return("fX")
    stop("Genotype must start with 'mY' or 'fX'")
  }
  
  # Mendelian gamete production (no drive/editing)
  gametes_from_geno <- function(geno) {
    alleles <- split_geno_to_alleles(geno)
    if (alleles[1] == alleles[2]) {
      setNames(1, alleles[1])
    } else {
      setNames(c(0.5, 0.5), sort(alleles))
    }
  }
  
  # --- inheritance cube (ih): parents→offspring genotype probabilities ---
  cube_ih <- array(0, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                   dimnames = list(genotypes, genotypes, genotypes))
  
  for (f in female_genotypes) {
    gf <- gametes_from_geno(f)
    for (m in male_genotypes) {
      gm <- gametes_from_geno(m)
      
      for (aF in names(gf)) for (aM in names(gm)) {
        # sons (mY)
        son <- paste0("mY", paste0(sort(c(aF, aM)), collapse = ""))
        # daughters (fX)
        dau <- paste0("fX", paste0(sort(c(aF, aM)), collapse = ""))
        
        cube_ih[f, m, son] <- cube_ih[f, m, son] + 0.5 * gf[aF] * gm[aM]
        cube_ih[f, m, dau] <- cube_ih[f, m, dau] + 0.5 * gf[aF] * gm[aM]
      }
    }
  }
  
  # --- viability cube (tau): here genotype-only viability (same for all parent pairs) ---
  cube_tau <- array(1, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                    dimnames = list(genotypes, genotypes, genotypes))
  
  # Female-specific lethality for any H-bearing genotype
  lethal_genos_f <- c("fXHW", "fXHH")
  surv_female_H  <- 1 - female_lethality  # 0 in field; 1 in lab
  for (g in lethal_genos_f) {
    cube_tau[ , , g] <- surv_female_H
  }
  
  # --- fertility mask (phi): 0 for males, 1 for females (as in your original) ---
  phi <- c(rep(0, length(male_genotypes)), rep(1, length(female_genotypes)))
  names(phi) <- genotypes
  
  # # --- genotype-level fitness (omega): multiplicative from alleleFitness (optional) ---
  # calcGenotypeFitness <- function(geno, fitnessVals) {
  #   alleles <- split_geno_to_alleles(geno)
  #   prod(fitnessVals[alleles])
  # }
  # omega <- setNames(sapply(genotypes, calcGenotypeFitness, fitnessVals = alleleFitness), genotypes)
  # 
  # Non-multiplicative genotype fitness:
  # - WW -> omega_W  (taken from alleleFitness["W"], default 1)
  # - HW or HH -> omega_H (taken from alleleFitness["H"], default 1)
  calcGenotypeFitness <- function(geno, fitnessVals) {
    # defaults if not present
    omega_W <- if (!is.na(fitnessVals["W"])) fitnessVals["W"] else 1
    omega_H <- if (!is.na(fitnessVals["H"])) fitnessVals["H"] else omega_W
    
    # strip sex prefix ("mY"/"fX"): alleles are the remaining two chars
    alle <- substring(geno, 3)  # "WW", "HW", or "HH"
    
    if (alle == "WW") return(omega_W)
    if (grepl("H", alle)) return(omega_H)  # applies to HW or HH
    omega_W
  }
  
  # Named vector of genotype fitnesses (ω)
  omega <- setNames(
    vapply(genotypes, calcGenotypeFitness, numeric(1), fitnessVals = alleleFitness),
    genotypes
  )
  
  # --- additional vectors to mirror your structure ---
  xiF <- setNames(rep(1, length(genotypes)), genotypes)
  xiM <- setNames(rep(1, length(genotypes)), genotypes)
  eta <- matrix(1, nrow = length(genotypes), ncol = length(genotypes),
                dimnames = list(genotypes, genotypes))
  
  # Per-genotype survival shortcut (s): mirror tau’s genotype effect for convenience
  s <- setNames(rep(1, length(genotypes)), genotypes)
  s[lethal_genos_f] <- surv_female_H
  
  list(
    ih = cube_ih,
    tau = cube_tau,
    genotypesID = genotypes,
    genotypesN = length(genotypes),
    wildType = c("mYWW", "fXWW"),
    eta = eta,
    phi = phi,
    omega = omega,
    xiF = xiF,
    xiM = xiM,
    s = s,
    releaseType = "mYHH"   # standard f-RIDL release: heterozygous males
  )
}

# example
# cube_fRIDL_field <- generate_fRIDL_inheritance_cube(female_lethality = 1)   # field
cube_fRIDL_lab   <- generate_fRIDL_inheritance_cube(female_lethality = 1, alleleFitness = c(W = 1, H = 0.9, Y=1, X=1))   # tetracycline rescue
