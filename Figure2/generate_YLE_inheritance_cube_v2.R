#' Generate inheritance and fitness cubes for a Y-linked genome editor (YLE)
#'
#' This function constructs the full set of inheritance, viability, fertility and fitness
#' components required by MouseGD/MGDrivE for a Y-linked genome editor targeting
#' a single autosomal locus with three alleles:
#'   - W: wild-type
#'   - H: "drive-favoured" edited allele
#'   - R: resistant allele
#'
#' Male genotypes carry either:
#'   - Y  : wild-type Y chromosome  (prefix "mY")
#'   - y  : YLE-bearing Y chromosome (prefix "my")
#'
#' Female genotypes carry two X chromosomes (prefix "fX").
#'
#' Editing occurs only in YLE males (my***), following these rules:
#'   - myWW: both W alleles can be edited with probability mu, producing H or R
#'           according to probabilities (1 - q) and q.
#'   - myWH: the W allele is cleaved with probability c. Cleaved alleles are repaired
#'           by:
#'             * HDR with prob (1 - j), copying the homologous H allele
#'             * NHEJ with prob j, producing H with prob (1 - p) or R with prob p
#'   - myWR: analogous to myWH but the homologous allele is R.
#'   - All other genotypes follow Mendelian segregation.
#'
#' Female-specific lethality and sterility can be applied to genotypes carrying H/R,
#' under either haploinsufficient or (partially) haplosufficient assumptions.
#'
#' @param c  Cleavage probability of W in heterozygous yW(H/R) males.
#' @param j  Fraction of cleaved alleles repaired by NHEJ (vs HDR).
#' @param mu Editing probability in myWW homozygotes (both W alleles).
#' @param p  Fraction of NHEJ outcomes that generate R rather than H.
#' @param q  Fraction of edits in myWW that generate R rather than H.
#' @param female_lethality  Degree of female lethality (0 = none, 1 = full) for
#'        specified H/R-bearing genotypes.
#' @param female_sterility  Degree of female sterility (0 = none, 1 = full) for
#'        specified H/R-bearing genotypes.
#' @param Haploinsufficient Logical; if TRUE, a single edited allele (H/R) in females
#'        is sufficient to trigger lethality/sterility (haploinsufficient target).
#'        If FALSE, only HH/HR/RR females are affected (haplosufficient target).
#' @param muAd Baseline adult daily mortality used in the genotype fitness formula.
#' @param alleleFitness Named vector of per-allele fitness multipliers with entries
#'        c(W, H, R, Y, X, y). If NULL, defaults to c(W=1, H=1, R=1, Y=1, X=1, y=0.8).
#'
#' @return A list with components:
#'   \item{ih}{3D inheritance cube (parent_female x parent_male x offspring) with
#'             genotype-level transmission probabilities.}
#'   \item{tau}{3D viability cube (same dimensions as ih) giving genotype survival
#'              after zygote formation (here: female-specific lethality).}
#'   \item{genotypesID}{Character vector of genotype labels (rows/cols of cubes).}
#'   \item{genotypesN}{Number of genotypes.}
#'   \item{wildType}{Vector giving the wild-type male and female genotypes.}
#'   \item{eta}{Mating structure matrix (here all 1, i.e. random mating).}
#'   \item{phi}{Sex-specific fertility mask (0 for males, 1 for females).}
#'   \item{omega}{Genotype fitness vector (adult survival scaling).}
#'   \item{xiF, xiM}{Female and male mating fitness (here all 1).}
#'   \item{s}{Genotype-specific fertility scaling (used to encode female sterility).}
#'   \item{releaseType}{Genotype released (here, "myHH").}
#'
generate_YLE_inheritance_cube <- function(
    c = 0.93,
    j = 0.25,
    mu = 0.97,
    p = 0,
    q = 0,
    female_lethality  = 1,
    female_sterility  = 1,
    Haploinsufficient = TRUE,
    muAd              = 1/690,
    alleleFitness     = NULL
) {
  # Default per-allele fitness multipliers if none supplied
  if (is.null(alleleFitness)) {
    alleleFitness <- c(W = 1, H = 1, R = 1, Y = 1, X = 1, y = 0.8)
  }
  
  # ------------------------------------------------------------------
  # Genotype definitions
  # ------------------------------------------------------------------
  # Male genotypes: mY.. = wild-type Y; my.. = YLE Y; two autosomal alleles
  male_genotypes <- c(
    "mYWW", "mYHW", "mYRW", "mYHH", "mYHR", "mYRR",
    "myWW", "myHW", "myRW", "myHH", "myHR", "myRR"
  )
  # Female genotypes: fX.., two autosomal alleles
  female_genotypes <- c("fXWW", "fXHW", "fXRW", "fXHH", "fXHR", "fXRR")
  genotypes <- c(male_genotypes, female_genotypes)
  
  # ------------------------------------------------------------------
  # Inheritance cube (ih): P(offspring genotype | female, male)
  # ------------------------------------------------------------------
  cube_ih <- array(
    0,
    dim      = c(length(genotypes), length(genotypes), length(genotypes)),
    dimnames = list(genotypes, genotypes, genotypes)
  )
  
  alleleSet <- c("W", "H", "R")
  
  # Helper: extract autosomal alleles (2-char vector) from genotype label
  split_geno_to_alleles <- function(geno) {
    # autosomal part is always positions 3–5, e.g. "mYWW" -> "WW"
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
  
  # Helper: identify sex/Y-type from genotype prefix
  # Returns "Y" (wild-type Y), "y" (YLE Y) or "X" (female)
  get_Y_type <- function(geno) {
    if (startsWith(geno, "mY")) return("Y")
    if (startsWith(geno, "my")) return("y")
    if (startsWith(geno, "fX")) return("X")
    stop("Genotype not recognized: must start with 'mY', 'my', or 'fX'")
  }
  
  # ------------------------------------------------------------------
  # Gamete formation with Y-linked editing
  # ------------------------------------------------------------------
  gametes_from_geno <- function(geno) {
    # Extract autosomal alleles and sex/Y-type
    focus_geno <- substr(geno, 3, 5)
    Y_type     <- get_Y_type(geno)
    alleles    <- split_geno_to_alleles(focus_geno)
    
    # --- Editing only for YLE males ("my..") -------------------------
    
    # Case 1: y with WW (homozygous wild-type target)
    #   Each W is edited with prob mu, giving H or R.
    if (setequal(alleles, c("W", "W")) && Y_type == "y") {
      return(c(
        "W" = (1 - mu),        # no edit
        "H" = mu * (1 - q),    # edited to H
        "R" = mu * q           # edited to R
      ))
    }
    
    # Case 2: y with WH (heterozygote, homologous H)
    #   The W allele is cut with prob c, then:
    #     - HDR (1 - j): becomes H
    #     - NHEJ  j    : H with (1 - p) or R with p
    if (setequal(alleles, c("W", "H")) && Y_type == "y") {
      return(c(
        "W" = 0.5 * (1 - c),                                   # uncut W
        "H" = 0.5 + 0.5 * c * (1 - j) + 0.5 * j * c * (1 - p), # H from HDR or NHEJ
        "R" = 0.5 * (j * c * p)                                # R from NHEJ
      ))
    }
    
    # Case 3: y with WR (heterozygote, homologous R)
    #   Analogous logic but homologous allele is R.
    if (setequal(alleles, c("W", "R")) && Y_type == "y") {
      return(c(
        "W" = 0.5 * (1 - c),                      # uncut W
        "H" = 0.5 * (j * c * (1 - p)),            # H from NHEJ
        "R" = 0.5 + 0.5 * c * (1 - j) +           # R from HDR
          0.5 * c * j * p                     # R from NHEJ
      ))
    }
    
    # ----------------------------------------------------------------
    # Otherwise: Mendelian segregation
    # ----------------------------------------------------------------
    # Homozygotes
    if (setequal(alleles, c("H", "H"))) return(c("H" = 1))
    if (setequal(alleles, c("R", "R"))) return(c("R" = 1))
    if (setequal(alleles, c("W", "W"))) return(c("W" = 1))
    
    # Heterozygotes: 50/50 segregation
    return(setNames(rep(0.5, 2), alleles))
  }
  
  # ------------------------------------------------------------------
  # Fill inheritance cube: for each female × male, combine gametes
  # ------------------------------------------------------------------
  for (f in female_genotypes) {
    for (m in male_genotypes) {
      # Gamete distributions from each parent
      gf <- gametes_from_geno(f)
      gm <- gametes_from_geno(m)
      Y_type <- get_Y_type(m)  # determines son Y-chromosome type
      
      for (a1 in names(gf)) {
        for (a2 in names(gm)) {
          # Sons: mY.. or my.. depending on sire Y_type
          offspring_son      <- paste0(c("m", Y_type, sort(c(a1, a2))), collapse = "")
          # Daughters: always fX..
          offspring_daughter <- paste0(c("fX", sort(c(a1, a2))), collapse = "")
          
          # Each mating produces 50% sons and 50% daughters
          if (offspring_son %in% genotypes) {
            cube_ih[f, m, offspring_son] <-
              cube_ih[f, m, offspring_son] + 0.5 * gf[a1] * gm[a2]
          }
          if (offspring_daughter %in% genotypes) {
            cube_ih[f, m, offspring_daughter] <-
              cube_ih[f, m, offspring_daughter] + 0.5 * gf[a1] * gm[a2]
          }
        }
      }
    }
  }
  
  # ------------------------------------------------------------------
  # Viability cube (tau): here used for female-specific lethality
  # ------------------------------------------------------------------
  cube_tau <- array(
    1,
    dim      = c(length(genotypes), length(genotypes), length(genotypes)),
    dimnames = list(genotypes, genotypes, genotypes)
  )
  
  # Choose which female genotypes are affected by lethality/sterility
  if (Haploinsufficient) {
    # Any female carrying H or R is affected
    lethal_genos_f  <- c("fXHW", "fXHH", "fXHR", "fXRR")
    sterile_genos_f <- c("fXHW", "fXHH", "fXHR", "fXRR")
  } else {
    # Only HH/HR/RR females affected (haplosufficient)
    lethal_genos_f  <- c("fXHR", "fXHH", "fXRR")
    sterile_genos_f <- c("fXHR", "fXHH", "fXRR")
  }
  
  # Survival multiplier for lethal genotypes
  surv_female_L <- 1 - female_lethality  # 0 in field; 1 in lab
  
  # Apply viability reduction to all matings producing these female genotypes
  for (g in lethal_genos_f) {
    cube_tau[ , , g] <- surv_female_L
  }
  
  # ------------------------------------------------------------------
  # Fertility mask (phi) and genotype-level fertility scaling (s)
  # ------------------------------------------------------------------
  # phi: 0 for males (do not directly produce offspring), 1 for females
  phi <- c(
    rep(0, length(male_genotypes)),
    rep(1, length(female_genotypes))
  )
  
  # Fertility reduction factor for affected females
  fs <- 1 - female_sterility
  
  # s: genotype-specific fertility scaling
  if (Haploinsufficient) {
    # sterile: fXHW, fXHH, fXHR, fXRR
    #      WW HW RW HH HR RR
    s <- c(
      1, 1, 1, 1, 1, 1,   # mY...
      1, 1, 1, 1, 1, 1,   # my...
      1, fs, 1, fs, fs, fs  # fX...
    )
  } else {
    # sterile: fXHR, fXHH, fXRR
    #      WW HW RW HH HR RR
    s <- c(
      1, 1, 1, 1, 1, 1,   # mY...
      1, 1, 1, 1, 1, 1,   # my...
      1, 1, 1, fs, fs, fs  # fX...
    )
  }
  names(s) <- genotypes
  
  # ------------------------------------------------------------------
  # Genotype fitness (omega): adult survival adjustment
  # ------------------------------------------------------------------
  # Helper: compute genotype-specific fitness from alleleFitness
  calcGenotypeFitness <- function(geno, fitnessVals) {
    focus_geno <- substr(geno, 3, 5)
    SexType    <- get_Y_type(geno)        # "Y", "y", or "X"
    alleles    <- split_geno_to_alleles(focus_geno)
    alleles_tot <- c(alleles, SexType)    # autosomal + sex/Y-type
    Eff_fit     <- prod(fitnessVals[alleles_tot])
    
    # Closed-form mapping from allele fitness to adult mortality scaling
    # Ensures omega = 1 for baseline (Eff_fit = 1) and declines with reduced fitness.
    max(0, (1 - muAd / Eff_fit) / (1 - muAd))
  }
  
  omegaNew <- setNames(
    sapply(genotypes, calcGenotypeFitness, fitnessVals = alleleFitness),
    genotypes
  )
  
  # ------------------------------------------------------------------
  # Assemble output list in MouseGD/MGDrivE format
  # ------------------------------------------------------------------
  list(
    ih          = cube_ih,
    tau         = cube_tau,
    genotypesID = genotypes,
    genotypesN  = length(genotypes),
    wildType    = c("mYWW", "fXWW"),
    eta         = matrix(1, nrow = length(genotypes), ncol = length(genotypes),
                         dimnames = list(genotypes, genotypes)),
    phi         = phi,
    omega       = omegaNew,
    xiF         = setNames(rep(1, length(genotypes)), genotypes),
    xiM         = setNames(rep(1, length(genotypes)), genotypes),
    s           = s,
    releaseType = "myHH"
  )
}

# Example usage:
# cube_t1 <- generate_YLE_inheritance_cube(Haploinsufficient = TRUE)
# cube_t2 <- generate_YLE_inheritance_cube(Haploinsufficient = FALSE)
