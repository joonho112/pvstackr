## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----load-fixture-------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit

## ----get-draws----------------------------------------------------------------
get_draws(fit)

## ----get-estimates------------------------------------------------------------
est <- get_estimates(fit)

names(est)

## ----estimate-core------------------------------------------------------------
est[, c("term", "estimate", "se", "df", "df_method")]

## ----estimate-provenance------------------------------------------------------
est[, c("term", "parameter_scope", "target_source", "target_hash")]

## ----interval-metadata--------------------------------------------------------
est[, c("term", "interval_role", "df_method",
        "df_complete", "coverage_claim_allowed")]

## ----fmi----------------------------------------------------------------------
tg <- get_target(fit)

data.frame(
  term = tg$fe_names,
  fmi  = round(as.numeric(tg$fmi), 3),
  riv  = round(as.numeric(tg$riv), 2),
  df   = round(as.numeric(tg$df),  2)
)

## ----target-fields------------------------------------------------------------
diag(tg$T_MI)        # total target variance per coefficient
round(tg$df, 3)      # degrees of freedom per coefficient
tg$target_hash       # content fingerprint (matches the estimate table)

## ----target-identity----------------------------------------------------------
all.equal(est$estimate, unname(tg$beta_bar))            # estimates == pooled beta
all.equal(est$se,       unname(sqrt(diag(tg$T_MI))))    # SEs == sqrt(diag T_MI)

## ----get-draws-again----------------------------------------------------------
get_draws(fit)

## ----get-diagnostics----------------------------------------------------------
dg <- get_diagnostics(fit)

names(dg)

## ----ccc-peek-----------------------------------------------------------------
dg$ccc[c("center_status", "delta_c_rel", "delta_c_max")]

## ----session-info-------------------------------------------------------------
sessionInfo()

