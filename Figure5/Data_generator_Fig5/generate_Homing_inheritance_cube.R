generate_Homing_inheritance_cube <- function(
    pC = 1,           # Probability of successful cut
    pN = 0,           # Probability of NHEJ (given cut occurs)
    pL = 0,            # Probability of loss-of-function after NHEJ
    alleleFitness = NULL # Optional per-allele fitness costs
) {
  
  # =========================================================================
  # HOMING SUPPRESSION DRIVE SYSTEM
  # =========================================================================
  # Target: Haplosufficient female fertility gene (autosomal)
  # Alleles:
  #   W = Wild-type (functional)
  #   H = Drive allele (nonfunctional - contains CRISPR machinery)
  #   R = Resistant functional (NHEJ repair that maintains function)
  #   N = Resistant nonfunctional (NHEJ repair that disrupts function)
  #
  # Fertility effects:
  #   - Females need at least 1 functional allele (W or R) to be fertile
  #   - Females with HH, HN, or NN are infertile
  #   - Males are unaffected (gene not required for male fertility)
  # =========================================================================
  
  # Default fitness (no cost)
  if (is.null(alleleFitness)) {
    alleleFitness <- c(W = 1, H = 1, R = 1, N = 1)
  }
  
  # Define allele set
  alleleSet <- c("W", "H", "R", "N")
  
  # Generate all possible genotypes (autosomal, so same structure for both sexes)
  # Prefix: mY = male, fX = female
  # Create all combinations including order (WH and HW are same, but we'll handle this)
  all_combinations <- expand.grid(alleleSet, alleleSet, stringsAsFactors = FALSE)
  # Sort alleles alphabetically to ensure consistent naming (e.g., HW becomes HW)
  genotype_cores <- unique(apply(all_combinations, 1, function(x) paste0(sort(x), collapse = "")))
  
  male_genotypes   <- paste0("mY", genotype_cores)
  female_genotypes <- paste0("fX", genotype_cores)
  genotypes <- c(male_genotypes, female_genotypes)
  
  # =========================================================================
  # HELPER FUNCTIONS
  # =========================================================================
  
  # Extract alleles from genotype string
  split_geno_to_alleles <- function(geno) {
    core <- substr(geno, 3, nchar(geno))
    if (nchar(core) != 2) stop("Invalid genotype format: ", geno)
    a1 <- substr(core, 1, 1)
    a2 <- substr(core, 2, 2)
    if (!(a1 %in% alleleSet && a2 %in% alleleSet)) {
      stop("Invalid alleles in genotype: ", geno)
    }
    c(a1, a2)
  }
  
  # Get sex prefix
  get_sex_prefix <- function(geno) {
    if (startsWith(geno, "mY")) return("mY")
    if (startsWith(geno, "fX")) return("fX")
    stop("Genotype must start with 'mY' or 'fX'")
  }
  
  # =========================================================================
  # GAMETE PRODUCTION WITH HOMING DRIVE MECHANICS
  # =========================================================================
  
  gametes_from_geno <- function(geno) {
    alleles <- split_geno_to_alleles(geno)
    a1 <- alleles[1]
    a2 <- alleles[2]
    
    # Homozygotes: simple Mendelian segregation
    if (a1 == a2) {
      return(setNames(1, a1))
    }
    
    # Heterozygote with Drive allele: HOMING OCCURS
    # Drive (H) attempts to convert wild-type (W) allele
    if (setequal(alleles, c("W", "H"))) {
      # Homing drive mechanics:
      # - Successful homing: pC * (1 - pN) → produces H
      # - NHEJ with function retained: pC * pN * (1 - pL) → produces R
      # - NHEJ with function lost: pC * pN * pL → produces N
      # - No cut: (1 - pC) → produces W
      
      p_homing <- pC * (1 - pN)          # H copies itself
      p_resist_functional <- pC * pN * (1 - pL)  # R created
      p_resist_nonfunctional <- pC * pN * pL     # N created
      p_no_cut <- 1 - pC                  # W transmitted
      
      return(c(
        W = 0.5 * p_no_cut,              # Wild-type from non-drive chromosome
        H = 0.5 + 0.5 * p_homing,        # Drive + successful homing
        R = 0.5 * p_resist_functional,   # Functional resistant
        N = 0.5 * p_resist_nonfunctional # Nonfunctional resistant
      ))
    }
    
    # Heterozygote: Drive with resistant alleles
    # Drive can still attempt to cut resistant alleles, but typically
    # resistant alleles are immune. We'll assume no further homing.
    if (a1 == "H" || a2 == "H") {
      # Simple Mendelian segregation (resistant alleles block homing)
      return(setNames(c(0.5, 0.5), sort(alleles)))
    }
    
    # All other heterozygotes: standard Mendelian segregation
    return(setNames(c(0.5, 0.5), sort(alleles)))
  }
  
  # =========================================================================
  # BUILD INHERITANCE CUBE (ih)
  # =========================================================================
  
  cube_ih <- array(0, 
                   dim = c(length(genotypes), length(genotypes), length(genotypes)),
                   dimnames = list(genotypes, genotypes, genotypes))
  
  for (f in female_genotypes) {
    gf <- gametes_from_geno(f)
    
    for (m in male_genotypes) {
      gm <- gametes_from_geno(m)
      
      # Cross all possible gamete combinations
      for (aF in names(gf)) {
        for (aM in names(gm)) {
          # Create offspring genotype (sorted alphabetically)
          offspring_alleles <- paste0(sort(c(aF, aM)), collapse = "")
          
          # Sons and daughters (1:1 sex ratio)
          son <- paste0("mY", offspring_alleles)
          dau <- paste0("fX", offspring_alleles)
          
          prob <- gf[aF] * gm[aM]
          
          cube_ih[f, m, son] <- cube_ih[f, m, son] + 0.5 * prob
          cube_ih[f, m, dau] <- cube_ih[f, m, dau] + 0.5 * prob
        }
      }
    }
  }
  
  # =========================================================================
  # BUILD VIABILITY CUBE (tau)
  # =========================================================================
  
  cube_tau <- array(1, 
                    dim = c(length(genotypes), length(genotypes), length(genotypes)),
                    dimnames = list(genotypes, genotypes, genotypes))
  
  # Female fertility loss: need at least one functional allele (W or R)
  # Infertile genotypes: HH, HN, NN (and NH which is same as HN)
  infertile_female_genos <- c("fXHH", "fXHN", "fXNN")
  
  # Infertile females cannot produce offspring (viability = 0)
  for (g in infertile_female_genos) {
    cube_tau[g, , ] <- 0  # No offspring from infertile females
  }
  
  # =========================================================================
  # FERTILITY MASK (phi)
  # =========================================================================
  
  # phi: 0 for males, 1 or 0 for females
  phi <- setNames(rep(1, length(genotypes)), genotypes)
  phi[male_genotypes] <- rep(0, length(male_genotypes))
  
  # =========================================================================
  # GENOTYPE FITNESS (omega)
  # =========================================================================
  
  calcGenotypeFitness <- function(geno, fitnessVals) {
    alleles <- split_geno_to_alleles(geno)
    prod(fitnessVals[alleles])
  }
  
  omega <- setNames(
    sapply(genotypes, calcGenotypeFitness, fitnessVals = alleleFitness),
    genotypes
  )
  
  # =========================================================================
  # ADDITIONAL PARAMETERS
  # =========================================================================
  
  # Sex-specific viability modifiers (default: no effect)
  xiF <- setNames(rep(1, length(genotypes)), genotypes)
  xiM <- setNames(rep(1, length(genotypes)), genotypes)
  
  # Mating compatibility matrix (default: all compatible)
  eta <- matrix(1, 
                nrow = length(genotypes), 
                ncol = length(genotypes),
                dimnames = list(genotypes, genotypes))
  
  # Per-genotype survival (mirrors fertility for females)
  s <- setNames(rep(1, length(genotypes)), genotypes)
  
  # =========================================================================
  # RETURN INHERITANCE CUBE OBJECT
  # =========================================================================
  
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
    releaseType = "mYHH",  # Release homozygous drive males
    
    # Store parameters for reference
    parameters = list(
      pC = pC,
      pN = pN,
      pL = pL,
      alleleFitness = alleleFitness
    ),
    
    # Metadata
    description = "Homing suppression drive targeting haplosufficient female viability gene"
  )
}

# =========================================================================
# EXAMPLE USAGE
# =========================================================================

# High efficiency drive (field conditions)
# cube <- generate_Homing_inheritance_cube(
#   pC = 1,   # 95% cutting efficiency
#   pN = 0.0,   # 1% NHEJ when cut occurs
#   pL = 1/3,     # 100/3 % of NHEJ events cause loss-of-function
#   alleleFitness = c(W = 1, H = 1, R = 1, N = 0.6)  # Drive and R have costs
# )

# # Lower efficiency drive
# cube_homing_low <- generate_Homing_inheritance_cube(
#   pC = 0.80,   # 80% cutting efficiency
#   pN = 0.10,   # 10% NHEJ when cut occurs
#   pL = 0.85    # 85% of NHEJ events cause loss-of-function
# )
# 
# # With fitness costs on drive allele
# cube_homing_cost <- generate_Homing_inheritance_cube(
#   pC = 0.95,
#   pN = 0.05,
#   pL = 0.9,
#   alleleFitness = c(W = 1, H = 0.95, R = 0.98, N = 1)  # Drive and R have costs
# )

# # Print summary
# print(paste("Total genotypes:", cube_homing_field$genotypesN))
# print(paste("Release genotype:", cube_homing_field$releaseType))
# print(paste("Wild-type genotypes:", paste(cube_homing_field$wildType, collapse = ", ")))