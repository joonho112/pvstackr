## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----common-load--------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # the BRR–Fay target of the fit, made by pv_brr_target()

## ----stacked-wls--------------------------------------------------------------
pisa_tiny <- read.csv(system.file("extdata", "pisa_tiny.csv", package = "pvstackr"))
X <- model.matrix(~ x + female, data = pisa_tiny)     # the design matrix
w <- pisa_tiny$W_FSTUWT / mean(pisa_tiny$W_FSTUWT)    # normalized weights
M <- 2
stacked_y <- c(pisa_tiny$PV1READ, pisa_tiny$PV2READ)  # one copy per plausible value
stacked_X <- rbind(X, X)
b_stacked <- lm.wfit(stacked_X, stacked_y, rep(w / M, M))$coefficients
max(abs(unname(b_stacked) - unname(tg$beta_bar)))     # rounding error only

b_equal <- lm.wfit(stacked_X, stacked_y, rep(1 / M, M * nrow(X)))$coefficients
max(abs(unname(b_equal) - unname(tg$beta_bar)))       # rows weighted 1/M only

## ----estimate-from-target-----------------------------------------------------
est <- get_estimates(fit)

all.equal(est$estimate, unname(tg$beta_bar))   # TRUE: the estimate is copied from beta_bar
data.frame(
  term             = est$term,
  stacked_estimate = est$estimate,
  rubin_beta_bar   = unname(tg$beta_bar)
)

## ----session-info-------------------------------------------------------------
sessionInfo()

