## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----design-------------------------------------------------------------------
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

design <- pv_design(
  pisa_tiny,
  formula     = OUTCOME ~ x + female,   # OUTCOME stands for each plausible value
  pv_suffix   = "READ",                 # finds PV1READ and PV2READ
  expected_M  = 2L,                     # stop unless 2 plausible values are found
  expected_R  = 4L,                     # stop unless 4 replicate weights are found
  id_cols     = "CNTSTUID"             # the column that identifies each row
)

design

## ----design-fields------------------------------------------------------------
design$pv_cols
design$weight_col
design$rep_weight_cols

## ----detect-columns-----------------------------------------------------------
detect_pisa_pv_columns(pisa_tiny, suffix = "READ")
detect_pisa_brr_replicate_weights(pisa_tiny)

## ----target-------------------------------------------------------------------
target <- pv_brr_target(
  pisa_tiny,
  formula         = OUTCOME ~ x + female,
  pv_cols         = design$pv_cols,
  weight_col      = design$weight_col,
  rep_weight_cols = design$rep_weight_cols,
  fay_k           = design$fay_k,
  id_cols         = design$id_cols
)

target

## ----target-fields------------------------------------------------------------
target$beta
target$T_MI
sqrt(diag(target$T_MI))   # these become the reported standard errors

## ----target-group-term, error = TRUE------------------------------------------
try({
pv_brr_target(
  pisa_tiny,
  formula         = OUTCOME ~ x + (1 | CNTSCHID),
  pv_cols         = design$pv_cols,
  weight_col      = design$weight_col,
  rep_weight_cols = design$rep_weight_cols,
  fay_k           = design$fay_k
)
})

## ----example-fit--------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit

fit

## ----fit-bundled-engine, eval = FALSE-----------------------------------------
# fit <- pv_fit(
#   data    = pisa_tiny,
#   formula = OUTCOME ~ x + female,
#   target  = target,              # from pv_brr_target() above
#   method  = "stack_direct",
#   control = pv_control(method = "stack_direct", backend = "brms")
# )

## ----fit-own-functions, eval = FALSE------------------------------------------
# fit <- pv_fit(
#   data              = pisa_tiny,
#   formula           = OUTCOME ~ x + female,
#   target            = target,
#   method            = "stack_direct",
#   control           = pv_control(method = "stack_direct"),
#   fit_function      = your_fit_function,      # fits the stacked model
#   draws_function    = your_draws_function,    # returns the draws of that fit
#   diagnose_function = your_diagnose_function, # returns its sampler diagnostics
#   cache_dir         = NULL                    # no cache file for your function
# )

## ----summary------------------------------------------------------------------
summary(fit)

## ----get-estimates------------------------------------------------------------
est <- get_estimates(fit)

est[, c("term", "estimate", "se", "df",
        "conf_low", "conf_high",
        "interval_role", "coverage_claim_allowed")]

## ----own-functions-sketch, eval = FALSE---------------------------------------
# # fit_function: fit the stacked model and return the fitted object. pvstackr
# # calls it with the named arguments formula, data, family, prior, chains, iter,
# # warmup, cores, seed, backend, file and file_refit, plus the elements of
# # additional_args. There is no weights argument: the weights are in the
# # formula term weights(.pvstackr_weight) and in the data column
# # .pvstackr_weight.
# my_fit_function <- function(formula, data, ...) {
#   # ... call your Bayesian engine here ...
# }
# 
# # draws_function: return the draws of that fit as a numeric matrix with one
# # row per draw. Name the fixed-effect columns as in target$fe_names
# # (b_Intercept, b_x, ...) or after the stacked formula (b_pvstackrMM001, ...),
# # or pass param_map.
# my_draws_function <- function(backend_fit, ...) {
#   # ... extract a draws matrix from backend_fit ...
# }
# 
# # diagnose_function: return the sampler diagnostics as a named list with six
# # values. If it is missing, fails or leaves out a value, the fit is blocked.
# my_diagnose_function <- function(backend_fit, ...) {
#   # ... return list(rhat_max = , ess_bulk_min = , ess_tail_min = ,
#   #                 divergences = , chains = , post_warmup_draws_per_chain = )
# }
# 
# fit <- pv_fit(
#   data              = your_real_pisa_data,
#   formula           = OUTCOME ~ x + female,
#   target            = your_target,
#   method            = "stack_direct",
#   control           = pv_control(method = "stack_direct"),
#   fit_function      = my_fit_function,
#   draws_function    = my_draws_function,
#   diagnose_function = my_diagnose_function,
#   cache_dir         = NULL
# )

