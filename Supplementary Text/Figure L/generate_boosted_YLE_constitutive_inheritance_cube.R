generate_boosted_YLE_constitutive_inheritance_cube <- function(c = 0.93, j = 0.25, mu = 0.97, p = 0, q = 0, 
                                                    female_lethality = 0, female_sterility = 1, 
                                                    Haploinsufficient = TRUE, muAd = 1/690, 
                                                    alleleFitness = NULL, x_shred = 0.95) {
    
    # 1. Setup Fitness and Genotypes
    if (is.null(alleleFitness)) {
      alleleFitness <- c(W = 1, H = 1, R = 1, Y = 1, X = 1, y = 0.8, A = 1, B = 1)
    }
    
    # Define components for 3-locus system: Sex/YLE + Target + Autosomal Shredder
    sex_types <- c("mY", "my", "fX")
    target_genotypes <- c("WW", "HW", "RW", "HH", "HR", "RR")
    shredder_genotypes <- c("AA", "AB", "BB")
    
    genotypes <- c()
    for(s in sex_types) {
      for(t in target_genotypes) {
        for(sh in shredder_genotypes) {
          genotypes <- c(genotypes, paste0(s, t, sh))
        }
      }
    }
    
    male_genotypes <- genotypes[startsWith(genotypes, "m")]
    female_genotypes <- genotypes[startsWith(genotypes, "f")]
    
    # Initialize cube
    cube_ih <- array(0, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                     dimnames = list(genotypes, genotypes, genotypes))
    
    # 2. Helper Functions
    get_Y_type <- function(geno) substr(geno, 2, 2)
    
    split_geno <- function(geno) {
      list(
        target = c(substr(geno, 3, 3), substr(geno, 4, 4)),
        shredder = c(substr(geno, 5, 5), substr(geno, 6, 6))
      )
    }
    
    # Fixed gametes function: Returns a 3x2 matrix (Target x Shredder)
    gametes_from_geno <- function(geno) {
      Y_type <- get_Y_type(geno)
      parts <- split_geno(geno)
      
      # Target Locus (YLE Logic)
      t_alleles <- parts$target
      t_probs <- c(W=0, H=0, R=0)
      
      if (Y_type == "y") {
        # Drive/Editor active logic
        if (all(t_alleles == "W")) {
          t_probs <- c(W = (1-mu), H = mu*(1-q), R = mu*q)
        } else if (setequal(t_alleles, c("W", "H"))) {
          t_probs <- c(W = 0.5*(1-c), H = 0.5 + 0.5*c*(1-j) + 0.5*j*c*(1-p), R = 0.5*(j*c*p))
        } else if (setequal(t_alleles, c("W", "R"))) {
          t_probs <- c(W = 0.5*(1-c), H = 0.5*(j*c*(1-p)), R = 0.5 + 0.5*c*(1-j) + 0.5*c*j*p)
        } else {
          if(all(t_alleles=="H")) t_probs["H"]=1 else if(all(t_alleles=="R")) t_probs["R"]=1 else {t_probs["H"]=0.5; t_probs["R"]=0.5}
        }
      } else {
        # Standard Mendelian for mY and fX
        if(t_alleles[1] == t_alleles[2]) t_probs[t_alleles[1]] = 1
        else { t_probs[t_alleles[1]] = 0.5; t_probs[t_alleles[2]] = 0.5 }
      }
      
      # Shredder Locus (Mendelian A/B)
      sh_alleles <- parts$shredder
      sh_probs <- c(A=0, B=0)
      if(sh_alleles[1] == sh_alleles[2]) sh_probs[sh_alleles[1]] = 1
      else { sh_probs["A"] = 0.5; sh_probs["B"] = 0.5 }
      
      return(outer(t_probs, sh_probs)) # Returns matrix with proper names
    }
    
    # 3. Fill Inheritance Cube
    for (f in female_genotypes) {
      gf_mat <- gametes_from_geno(f)
      for (m in male_genotypes) {
        gm_mat <- gametes_from_geno(m)
        
        # Determine sex-ratio bias based on Autosomal Shredder (A allele)
        m_shredder <- split_geno(m)$shredder
        prob_Y <- if (any(m_shredder == "A")) 1/(2 - x_shred) else 0.5
        prob_X <- 1 - prob_Y
        
        Y_type <- get_Y_type(m) # 'Y' or 'y'
        
        # Iterate over non-zero probability gametes
        target_indices <- which(rowSums(gf_mat) > 0 | rowSums(gm_mat) > 0)
        shred_indices  <- which(colSums(gf_mat) > 0 | colSums(gm_mat) > 0)
        
        for (f_t in rownames(gf_mat)) {
          for (f_s in colnames(gf_mat)) {
            if (gf_mat[f_t, f_s] == 0) next
            
            for (m_t in rownames(gm_mat)) {
              for (m_s in colnames(gm_mat)) {
                if (gm_mat[m_t, m_s] == 0) next
                
                prob_joint <- gf_mat[f_t, f_s] * gm_mat[m_t, m_s]
                
                # Construct offspring strings
                off_t <- paste0(sort(c(f_t, m_t)), collapse="")
                off_s <- paste0(sort(c(f_s, m_s)), collapse="")
                
                son <- paste0("m", Y_type, off_t, off_s)
                dau <- paste0("fX", off_t, off_s)
                
                cube_ih[f, m, son] <- cube_ih[f, m, son] + prob_joint * prob_Y
                cube_ih[f, m, dau] <- cube_ih[f, m, dau] + prob_joint * prob_X
              }
            }
          }
        }
      }
    }
    
    # 4. Viability and Fitness
    cube_tau <- array(1, dim = c(length(genotypes), length(genotypes), length(genotypes)),
                      dimnames = list(genotypes, genotypes, genotypes))
    
    # Fertility Mask (phi) and Sterility Logic
    phi <- setNames(rep(0, length(genotypes)), genotypes)
    s <- setNames(rep(1, length(genotypes)), genotypes)
    fs <- 1 - female_sterility
    
    for (g in genotypes) {
      if (startsWith(g, "f")) {
        phi[g] <- 1
        target <- substr(g, 3, 4)
        is_sterile <- if (Haploinsufficient) target %in% c("HW", "HH", "HR", "RR") else target %in% c("HH", "HR", "RR")
        if (is_sterile) {
          cube_tau[,,g] <- 1 - female_lethality
          s[g] <- fs
        }
      }
    }
    
    # Calculate final fitness (Omega)
    calcGenotypeFitness <- function(geno, fitnessVals) {
      Sex <- substr(geno, 2, 2)
      Tar <- c(substr(geno, 3, 3), substr(geno, 4, 4))
      Shr <- c(substr(geno, 5, 5), substr(geno, 6, 6))
      eff_fit <- prod(fitnessVals[c(Sex, Tar, Shr)])
      return(max(0, (1 - muAd / eff_fit) / (1 - muAd)))
    }
  
    
    omegaNew <- sapply(genotypes, calcGenotypeFitness, fitnessVals = alleleFitness)
    
    return(list(
      ih = cube_ih, 
      tau = cube_tau, 
      genotypesID = genotypes,
      genotypesN = length(genotypes),
      omega = omegaNew, 
      eta = matrix(1, nrow = length(genotypes), ncol = length(genotypes),
                   dimnames = list(genotypes, genotypes)),
      phi = phi, 
      s = s, 
      xiF = setNames(rep(1, length(genotypes)), genotypes),
      xiM = setNames(rep(1, length(genotypes)), genotypes),
      wildType = c("mYWWBB", "fXWWBB"), 
      releaseType = "myHHAA"
    ))
  }
  