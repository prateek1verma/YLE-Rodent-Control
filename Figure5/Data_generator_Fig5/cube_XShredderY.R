# -------------------------------------------------------------------
# Y-linked X-Shredder inheritance cube (sex-tagged, CODE-2 style)
# -------------------------------------------------------------------
cube_XShredderY <- function(cX = 1, crX = 0, cB = 0,
                                   eta = NULL, phi = NULL, omega = NULL,
                                   xiF = NULL, xiM = NULL, s = NULL) {
  # sanity
  if (any(c(cX, crX, cB) > 1) || any(c(cX, crX, cB) < 0)) {
    stop("cX, crX, cB must be in [0,1].")
  }
  
  # sex-tagged genotype IDs (female first, then male)
  fem <- c("fXX","fXR","fRR")
  mal <- c("mXY","mXA","mXB","mRY","mRA","mRB")
  gtype <- c(fem, mal)
  n <- length(gtype)
  
  # 3D inheritance array: parents x parents x offspring (all genotypes listed)
  ih <- array(0, dim = c(n, n, n), dimnames = list(gtype, gtype, gtype))
  
  # ---- mapping of original CODE 1 rules to sex-tagged names ----
  # helper aliases for brevity
  fXX <- "fXX"; fXR <- "fXR"; fRR <- "fRR"
  mXY <- "mXY"; mXA <- "mXA"; mXB <- "mXB"; mRY <- "mRY"; mRA <- "mRA"; mRB <- "mRB"
  
  # Original: tMatrix['XX','XY',c('XX','XY')] <- c(1,1)/2
  ih[fXX, mXY, c(fXX, mXY)] <- c(1,1)/2
  
  # Original: tMatrix['XX','XA',c('XX','XR','XA','XB')] <- c(1-cX, cX*crX, 1-cB, cB)/(2 + cX*(crX-1))
  ih[fXX, mXA, c(fXX, fXR, mXA, mXB)] <- c(1 - cX, cX*crX, 1 - cB, cB) / (2 + cX*(crX - 1))
  
  # Original: tMatrix['XX','XB',c('XX','XB')] <- c(1,1)/2
  ih[fXX, mXB, c(fXX, mXB)] <- c(1,1)/2
  
  # Original: tMatrix['XX','RY',c('XR','XY')] <- c(1,1)/2
  ih[fXX, mRY, c(fXR, mXY)] <- c(1,1)/2
  
  # Original: tMatrix['XX','RA',c('XR','XA','XB')] <- c(1,1-cB,cB)/2
  ih[fXX, mRA, c(fXR, mXA, mXB)] <- c(1, 1 - cB, cB)/2
  
  # Original: tMatrix['XX','RB',c('XR','XB')] <- c(1,1)/2
  ih[fXX, mRB, c(fXR, mXB)] <- c(1,1)/2
  
  # Original: tMatrix['XR','XY',c('XX','XR','XY','RY')] <- rep(1,4)/4
  ih[fXR, mXY, c(fXX, fXR, mXY, mRY)] <- 1/4
  
  # Original:
  # tMatrix['XR','XA',c('XX','XR','XA','XB','RR','RA','RB')] <-
  #   c(1-cX, cX*crX + 1-cX, 1-cB, cB, cX*crX, 1-cB, cB)/(4 + 2*cX*(crX-1))
  ih[fXR, mXA, c(fXX, fXR, mXA, mXB, fRR, mRA, mRB)] <-
    c(1 - cX, cX*crX + 1 - cX, 1 - cB, cB, cX*crX, 1 - cB, cB) / (4 + 2*cX*(crX - 1))
  
  # Original: tMatrix['XR','XB',c('XX','XR','XB','RB')] <- rep(1,4)/4
  ih[fXR, mXB, c(fXX, fXR, mXB, mRB)] <- 1/4
  
  # Original: tMatrix['XR','RY',c('XR','RR','XY','RY')] <- rep(1,4)/4
  ih[fXR, mRY, c(fXR, fRR, mXY, mRY)] <- 1/4
  
  # Original:
  # tMatrix['XR','RA',c('XR','XA','XB','RR','RA','RB')] <- c(1,1-cB,cB, 1,1-cB,cB)/4
  ih[fXR, mRA, c(fXR, mXA, mXB, fRR, mRA, mRB)] <- c(1, 1 - cB, cB, 1, 1 - cB, cB)/4
  
  # Original: tMatrix['XR','RB',c('XR','RR','XB','RB')] <- rep(1,4)/4
  ih[fXR, mRB, c(fXR, fRR, mXB, mRB)] <- 1/4
  
  # Original: tMatrix['RR','XY',c('XR','RY')] <- c(1,1)/2
  ih[fRR, mXY, c(fXR, mRY)] <- c(1,1)/2
  
  # Original: tMatrix['RR','XA',c('XR','RR','RA','RB')] <- c(1-cX, cX*crX, 1-cB, cB)/(2 + cX*(crX-1))
  ih[fRR, mXA, c(fXR, fRR, mRA, mRB)] <- c(1 - cX, cX*crX, 1 - cB, cB) / (2 + cX*(crX - 1))
  
  # Original: tMatrix['RR','XB',c('XR','RB')] <- c(1,1)/2
  ih[fRR, mXB, c(fXR, mRB)] <- c(1,1)/2
  
  # Original: tMatrix['RR','RY',c('RR','RY')] <- c(1,1)/2
  ih[fRR, mRY, c(fRR, mRY)] <- c(1,1)/2
  
  # Original: tMatrix['RR','RA',c('RR','RA','RB')] <- c(1,1-cB,cB)/2
  ih[fRR, mRA, c(fRR, mRA, mRB)] <- c(1, 1 - cB, cB)/2
  
  # Original: tMatrix['RR','RB',c('RR','RB')] <- c(1,1)/2
  ih[fRR, mRB, c(fRR, mRB)] <- c(1,1)/2
  
  # numerical cleanliness
  ih[ih < .Machine$double.eps] <- 0
  
  # viability mask (all 1s by default, as in CODE 1)
  tau <- array(1, dim = c(n, n, n), dimnames = list(gtype, gtype, gtype))
  
  # genotype-specific modifiers (defaults mirror CODE 1)
  if (!is.null(phi)) {
    stop("This cube has a sex-specific phi; edit after construction only if you are sure.")
  }
  phi_vec <- setNames(c(1,1,1, 0,0,0,0,0,0), gtype)  # females=1, males=0
  
  eta_mat <- if (is.null(eta)) {
    matrix(1, n, n, dimnames = list(gtype, gtype))
  } else eta
  omega_vec <- if (is.null(omega)) setNames(rep(1, n), gtype) else omega
  xiF_vec   <- if (is.null(xiF))   setNames(rep(1, n), gtype) else xiF
  xiM_vec   <- if (is.null(xiM))   setNames(rep(1, n), gtype) else xiM
  s_vec     <- if (is.null(s))     setNames(rep(1, n), gtype) else s
  
  list(
    ih = ih,
    tau = tau,
    genotypesID = gtype,
    genotypesN = n,
    wildType = c("fXX","mXY"),
    eta = eta_mat,
    phi = phi_vec,
    omega = omega_vec,
    xiF = xiF_vec,
    xiM = xiM_vec,
    s = s_vec,
    releaseType = "mXA"   # YLE release genotype (sex-tagged)
  )
}

# ---- example ----
# cube <- cube_XShredderY_tagged(cX = 1, crX = 0, cB = 0)
