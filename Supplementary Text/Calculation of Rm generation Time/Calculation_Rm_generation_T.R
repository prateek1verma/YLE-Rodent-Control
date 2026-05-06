## Basic reproductive number and generation time calculations

## Baseline parameters
t_Gest <- 19
t_Nurs <- 23
t_Ado  <- 37
t_Ad   <- 690

mu_Ad  <- 1 / 690
L      <- 7.5
Q      <- 6
K      <- 10000
Theta  <- 22.4

## Life-cycle timing
A_mat <- t_Nurs + t_Ado
A0    <- A_mat + t_Gest

## Daily female offspring production
fec <- (L / 365) * (Q / 2)

## Adult survival
s_low <- (1 - mu_Ad) * (1 / 2)^(1 / Theta)
s_K   <- (1 - mu_Ad) * (1 / 4)^(1 / Theta)

## Basic reproductive number at low density
R_m <- (fec * s_low^(A0 - A_mat)) / (1 - s_low)

## Euler-Lotka function
euler_lotka <- function(r, s, fec, A0, A_mat) {
  fec * s^(A0 - A_mat) * exp(-r * A0) / (1 - s * exp(-r)) - 1
}

## Solve for r
safe_uniroot <- function(s, fec, A0, A_mat) {
  singular <- log(s)
  
  lower <- singular + 1e-6
  upper <- 0.1
  
  uniroot(
    euler_lotka,
    interval = c(lower, upper),
    s = s,
    fec = fec,
    A0 = A0,
    A_mat = A_mat
  )$root
}

r_low <- safe_uniroot(s_low, fec, A0, A_mat)
r_K   <- safe_uniroot(s_K, fec, A0, A_mat)

## Generation time
generation_time <- function(r, s, fec, A0, A_mat) {
  x <- s * exp(-r)
  fec * s^(A0 - A_mat) * exp(-r * A0) *
    (A0 / (1 - x) + x / (1 - x)^2)
}

T_low <- generation_time(r_low, s_low, fec, A0, A_mat)
T_K   <- generation_time(r_K, s_K, fec, A0, A_mat)

## Output
results <- data.frame(
  quantity = c(
    "A_mat",
    "A0",
    "f",
    "s_low",
    "s_K",
    "R_m",
    "r_low",
    "T_low_days",
    "T_low_years",
    "r_K",
    "T_K_days",
    "T_K_years"
  ),
  value = c(
    A_mat,
    A0,
    fec,
    s_low,
    s_K,
    R_m,
    r_low,
    T_low,
    T_low / 365,
    r_K,
    T_K,
    T_K / 365
  )
)

print(results, digits = 6)
