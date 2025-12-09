generate_YLE_inheritance_cube <- function(c = 0.93, j = 0.25, mu = 0.97, p = 0, q = 0, female_lethality = 1, female_sterility = 1, Haploinsufficient = TRUE, muAd = 1/690, alleleFitness = NULL) {
  
  if (is.null(alleleFitness)) {
    alleleFitness <- c(W = 1, H = 1, R = 1, Y=1, X=1, y=0.8)
  }
  
  # Define genotypes
  male_genotypes <- c("mYWW", "mYHW", "mYRW", "mYHH", "mYHR", "mYRR",
                      "myWW", "myHW", "myRW", "myHH", "myHR", "myRR")
  female_genotypes <- c("fXWW", "fXHW", "fXRW", "fXHH", "fXHR", "fXRR")
  genotypes <- c(male_genotypes, female_genotypes)
  
  # Initialize inheritance cube
  cube_ih <- array(0, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                   dimnames = list(genotypes, genotypes, genotypes))
  
  alleleSet <- c("W", "H", "R")
  
  # Function to extract alleles from genotype string
  split_geno_to_alleles <- function(geno) {
    for (a1 in alleleSet) {
      if (startsWith(geno, a1)) {
        a2 <- substring(geno, nchar(a1) + 1)
        if (a2 %in% alleleSet) {
          return(c(a1, a2))
        }
      }
    }
    stop("Invalid genotype: ", geno)
  }
  
  
  # Function: extract Y type ("Y" or "y")
  get_Y_type <- function(geno) {
    if (startsWith(geno, "mY")) return("Y")
    if (startsWith(geno, "my")) return("y")
    if (startsWith(geno, "fX")) return("X")
    stop("Genotype not recognized: must start with 'mY' or 'my' or 'fX'")
  }
  
  # Inheritance logic
  # Function: Generate gametes from a genotype
  gametes_from_geno <- function(geno) {
    focus_geno <- substr(geno, 3, 5)
    Y_type <- get_Y_type(geno)
    alleles <- split_geno_to_alleles(focus_geno)
    
    # Case: y with WW (WT homozygote)
    if (setequal(alleles, c("W","W")) && Y_type=="y") {
      return(c(
        "W" = (1-mu),
        "H" = mu*(1-q),        # resistant “H”
        "R" = mu*q     # resistant “R”
      ))
    }
    
    # Case: y with WH (heterozygote)
    if (setequal(alleles,c("W","H")) && Y_type=="y") {
      return(c(
        "W" = 0.5*(1-c),                  # uncut W survives
        "H" = 0.5 + 0.5*c*(1-j) + 0.5*j*c*(1-p), # HDR to H or NHEJ to H
        "R" = 0.5*(j*c*p)             # NHEJ to R
      ))
    }
    
    # Case: y with WR
    if (setequal(alleles,c("W","R")) && Y_type=="y") {
      return(c(
        "W" = 0.5*(1-c), # uncut W survives
        "H" = 0.5*(j*c*(1-p)), # NHEJ to H
        "R" = 0.5 + 0.5*c*(1-j) + 0.5*c*j*p # HDR to R or NHEJ to R
      ))
    }
    
    # Default Mendelian segregation
    # Homo zygotes
    if (setequal(alleles, c("H", "H"))) {return(c("H" = 1))}
    if (setequal(alleles, c("R", "R"))) {return(c("R" = 1))}
    if (setequal(alleles, c("W", "W"))) {return(c("W" = 1))}
    
    # Hetero zygotes 
    return(setNames(rep(0.5, 2), alleles))
  }
  
  for (f in female_genotypes) {
    for (m in male_genotypes) {
      gf <- gametes_from_geno(f)
      gm <- gametes_from_geno(m)
      Y_type <- get_Y_type(m)
      
      for (a1 in names(gf)) {
        for (a2 in names(gm)) {
          offspring_son <- paste0(c("m",Y_type,sort(c(a1, a2))), collapse = "")
          offspring_daughter <- paste0(c("fX",sort(c(a1, a2))), collapse = "")
          if (offspring_son %in% genotypes) {
            cube_ih[f, m, offspring_son] <- cube_ih[f, m, offspring_son] + 0.5*gf[a1]*gm[a2]
          }
          if (offspring_daughter %in% genotypes) {
            cube_ih[f, m, offspring_daughter] <- cube_ih[f, m, offspring_daughter] + 0.5*gf[a1]*gm[a2]
          }          
        }
      }
    }
  }
  
  
  # Viability cube
  cube_tau <- array(1, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                    dimnames = list(genotypes, genotypes, genotypes))
  
  # Female-specific lethality for any L-bearing genotype
  if (Haploinsufficient == TRUE) {
    lethal_genos_f <- c("fXHW", "fXHH", "fXHR", "fXRR")
    sterile_genos_f <- c("fXHW", "fXHH", "fXHR", "fXRR")
  } else {
    lethal_genos_f <- c("fXHR", "fXHH", "fXRR")
    sterile_genos_f <- c("fXHR", "fXHH", "fXRR")
  }
  surv_female_L  <- 1 - female_lethality  # 0 in field; 1 in lab
  for (g in lethal_genos_f) {
    cube_tau[ , , g] <- surv_female_L
  }
  
  
  # Fertility mask (phi): 0 for males, 1 or 0 for females
  phi <- c(rep(0, length(male_genotypes)),rep(1, length(female_genotypes)))
  
  fs <- 1 - female_sterility  # FOR Haploinsufficient == TRUE/FALSE
  
  # Fertility (s): sex-specific fitness
  
  if (Haploinsufficient == TRUE) {
    # sterile: fXHW, fXHH, fXHR, fXRR  
    #      WW HW RW HH HR RR
    s <- c(1, 1, 1, 1, 1, 1,   # mY...
           1, 1, 1, 1, 1, 1,   # my...
           1, fs, 1, fs, fs, fs)  # fX...
  } else {
    # sterile: fXHR, fXHH, fXRR  
    #      WW HW RW HH HR RR 
    s <- c(1, 1, 1, 1, 1, 1,   # mY...
           1, 1, 1, 1, 1, 1,   # my...
           1, 1, 1, fs, fs, fs)  # fX...
  }
  
  names(s) <- genotypes
  
  # Fitness (omega): default 1
  omegaNew <- setNames(rep(1, length(genotypes)), genotypes)
  
  # Genotype fitness
  calcGenotypeFitness <- function(geno, fitnessVals) {
    focus_geno <- substr(geno, 3, 5)
    Sextype <- get_Y_type(geno)
    alleles <- split_geno_to_alleles(focus_geno)
    alleles_tot <- c(alleles,Sextype)
    Eff_fit <- prod(fitnessVals[alleles_tot])
    max(0,(1 - muAd / Eff_fit)/(1-muAd)) # analytical close form formula
  }
  
  omegaNew <- setNames(sapply(genotypes, calcGenotypeFitness, fitnessVals = alleleFitness), genotypes)
  
  
  list(
    ih = cube_ih,
    tau = cube_tau,
    genotypesID = genotypes,
    genotypesN = length(genotypes),
    wildType = c("mYWW", "fXWW"),
    eta = matrix(1, nrow = length(genotypes), ncol = length(genotypes),
                 dimnames = list(genotypes, genotypes)),
    phi = phi,
    omega = omegaNew,
    xiF = setNames(rep(1, length(genotypes)), genotypes),
    xiM = setNames(rep(1, length(genotypes)), genotypes),
    s = s,
    releaseType = "myHH"
  )
}

# cube_t1 <- generate_YLE_inheritance_cube(c = 0.93, j = 0.25, mu = 0.97, p = 0, q = 0, female_lethality = 1, female_sterility = 1, Haploinsufficient = TRUE, muAd = 1/690, alleleFitness = NULL)
# cube_t2 <- generate_YLE_inheritance_cube(c = 0.93, j = 0.25, mu = 0.97, p = 0, q = 0, female_lethality = 1, female_sterility = 1, Haploinsufficient = FALSE, muAd = 1/690, alleleFitness = NULL)

