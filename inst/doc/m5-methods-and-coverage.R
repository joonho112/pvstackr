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

## ----psis-classifier----------------------------------------------------------
## PSIS Pareto-k threshold classification (EQ-PSIS), pure base-R
classify_k <- function(k) {
  cut(k, breaks = c(-Inf, 0.5, 0.7, Inf),
      labels = c("good", "borderline", "unreliable"), right = FALSE)
}

data.frame(
  k       = c(0.04, 0.20, 0.49, 0.50, 0.69, 0.70),
  verdict = classify_k(c(0.04, 0.20, 0.49, 0.50, 0.69, 0.70))
)

## ----coverage-flag------------------------------------------------------------
est <- get_estimates(fit)

unique(est$interval_role)            # "descriptive_classic_rubin"
unique(est$coverage_claim_allowed)   # FALSE -- even though method == "stack_direct"
est[, c("term", "interval_role", "df_method", "coverage_claim_allowed")]

## ----session-info-------------------------------------------------------------
sessionInfo()

