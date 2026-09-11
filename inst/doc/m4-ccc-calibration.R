## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----common-load--------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # the target of the fit, made by pv_brr_target()

is.null(get_draws(fit))

## ----center-control-----------------------------------------------------------
fit$control$center   # "target": the calibrated draws are centered at beta_bar

## ----ccc-invariant------------------------------------------------------------
set.seed(5104)                      # fixed seed, so the simulated draws are the same in every run
p <- 3L; S <- 4000L

## simulated draws of p = 3 coefficients with a chosen mean and covariance
Sigma_raw_true <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
raw <- sweep(matrix(rnorm(S * p), S, p) %*% chol(Sigma_raw_true),
             2L, c(10, -3, 0.5), `+`)

## target mean and covariance for the check, in place of beta_bar and T_MI
c_tgt     <- c(458, 47, 2)
Sigma_tgt <- matrix(c(4.0, 0.5, 0.0,
                      0.5, 2.0, 0.3,
                      0.0, 0.3, 1.5), 3, 3, byrow = TRUE)

## the map: beta_cal = c_tgt + A (beta - colMeans(raw)), with A = L_tgt L_raw^{-1}
L_tgt <- t(chol(Sigma_tgt))         # lower Cholesky factor of the target covariance
L_raw <- t(chol(cov(raw)))          # lower Cholesky factor of the covariance of the draws
A     <- L_tgt %*% solve(L_raw)     # the calibration matrix A

cal <- sweep(raw, 2L, colMeans(raw), `-`) %*% t(A)   # center the draws, then multiply by A
cal <- sweep(cal, 2L, c_tgt, `+`)                    # add the target mean

all.equal(unname(colMeans(cal)), c_tgt)      # TRUE (largest difference about 1e-13)
all.equal(unname(cov(cal)),      Sigma_tgt)  # TRUE (largest difference about 1e-14)

## ----ccc-diagnostics----------------------------------------------------------
dg <- get_diagnostics(fit)

names(dg)
dg$ccc[c("center_status", "delta_c_rel", "delta_c_max")]

## ----session-info-------------------------------------------------------------
sessionInfo()

