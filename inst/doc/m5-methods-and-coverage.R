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

## ----k-hat-check--------------------------------------------------------------
## pvstackr's check of the Pareto k-hat values of a stack_psis fit, written in
## base R to illustrate it; pvstackr applies it inside pv_fit_stack_psis().
## A value passes when it is finite and strictly below the cut-off.
k_passes <- function(k, threshold = 0.7) is.finite(k) & k < threshold

k <- c(0.20, 0.50, 0.69, 0.70, 0.95, NA)   # made-up k-hat values
data.frame(
  k            = k,
  passes_0.7   = k_passes(k),                               # the default cut-off
  passes_0.667 = k_passes(k, threshold = 1 - 1 / log10(1000))
)

## ----coverage-flag------------------------------------------------------------
est <- get_estimates(fit)

unique(est$interval_role)
unique(est$coverage_claim_allowed)
est[, c("term", "interval_role", "df_method", "coverage_claim_allowed")]

## ----session-info-------------------------------------------------------------
sessionInfo()

