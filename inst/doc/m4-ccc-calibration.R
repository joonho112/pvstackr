## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----common-load--------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # external Rubin / BRR–Fay target (M2's object); Sigma_tgt = T_MI

c(M = tg$M, R = tg$R, fay_k = tg$fay_k)   # 2, 4, 0.5
is.null(get_draws(fit))                   # TRUE -- no retained draws to calibrate

## ----center-control-----------------------------------------------------------
fit$control$center   # "target" -- the reportable convention (c = beta_bar)

## ----ccc-invariant------------------------------------------------------------
set.seed(5104)                      # M4 is the only Method vignette with randomness
p <- 3L; S <- 4000L

## an arbitrary raw draw cloud with a non-trivial mean and covariance
Sigma_raw_true <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
raw <- sweep(matrix(rnorm(S * p), S, p) %*% chol(Sigma_raw_true),
             2L, c(10, -3, 0.5), `+`)

## target moments to impose (stand in for c = beta_bar and Sigma_tgt = T_MI)
c_tgt     <- c(458, 47, 2)
Sigma_tgt <- matrix(c(4.0, 0.5, 0.0,
                      0.5, 2.0, 0.3,
                      0.0, 0.3, 1.5), 3, 3, byrow = TRUE)

## EQ-CCC: beta_cal = c + L_tgt L_raw^{-1} (beta - colMeans(raw))
L_tgt <- t(chol(Sigma_tgt))         # lower-triangular Cholesky factor of the target
L_raw <- t(chol(cov(raw)))          # lower-triangular Cholesky factor of the (empirical) raw cov
A     <- L_tgt %*% solve(L_raw)     # the affine transform L_tgt L_raw^{-1}

cal <- sweep(raw, 2L, colMeans(raw), `-`) %*% t(A)   # recenter, then linear map
cal <- sweep(cal, 2L, c_tgt, `+`)                    # shift to the target center

all.equal(unname(colMeans(cal)), c_tgt)      # TRUE  (max|diff| ~ 1e-13)
all.equal(unname(cov(cal)),      Sigma_tgt)  # TRUE  (max|diff| ~ 1e-15)

## ----ccc-diagnostics----------------------------------------------------------
dg <- get_diagnostics(fit)

names(dg)                                                    # preflight, stack_fit, stack_fit_warnings, ccc
dg$ccc[c("center_status", "delta_c_rel", "delta_c_max")]

## ----session-info-------------------------------------------------------------
sessionInfo()

