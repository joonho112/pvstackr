## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----common-load--------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # external Rubin / BRR–Fay target (M2's object)

c(M = tg$M, R = tg$R, fay_k = tg$fay_k)   # 2, 4, 0.5

## ----point-identity-----------------------------------------------------------
est <- get_estimates(fit)

all.equal(est$estimate, unname(tg$beta_bar))   # TRUE  -- EQ-THM22 realized on the fixture
data.frame(
  term             = est$term,
  stacked_estimate = est$estimate,
  rubin_beta_bar   = unname(tg$beta_bar)
)

## ----se-is-external-----------------------------------------------------------
all.equal(est$se, unname(sqrt(diag(tg$T_MI))))   # TRUE  -- the SE is M2's target, not the bridge's

## ----session-info-------------------------------------------------------------
sessionInfo()

