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
  formula     = OUTCOME ~ x + female,   # OUTCOME is a placeholder (see below)
  pv_suffix   = "READ",                 # matches PV1READ / PV2READ
  expected_M  = 2L,                     # assert exactly 2 plausible values
  expected_R  = 4L,                     # assert exactly 4 replicate weights
  id_cols     = "CNTSTUID"             # column(s) that identify a unique row
)

design

## ----design-fields------------------------------------------------------------
design$pv_cols
design$weight_col
design$rep_weight_cols

## ----detect-helpers-----------------------------------------------------------
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

## ----target-internals---------------------------------------------------------
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

## ----fit-cached---------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit

fit

## ----fit-live-shape, eval = FALSE---------------------------------------------
# fit <- pv_fit(
#   data           = pisa_tiny,
#   formula        = OUTCOME ~ x + female,
#   target         = target,                                  # from Section 3
#   method         = "stack_direct",
#   control        = pv_control(
#                      method  = "stack_direct",
#                      backend = "injected"                   # backend policy
#                    ),
#   fit_function   = your_fit_function,                       # the modelling backend
#   draws_function = your_draws_function                      # extract posterior draws
# )

## ----summary------------------------------------------------------------------
summary(fit)

## ----get-estimates------------------------------------------------------------
est <- get_estimates(fit)

est[, c("term", "estimate", "se", "df",
        "conf_low", "conf_high",
        "interval_role", "coverage_claim_allowed")]

## ----coef-figure, fig.width = 7, fig.height = 3, fig.cap = "Slope coefficients from the cached synthetic stack_direct fit, with 95% descriptive intervals. The intercept is omitted because its scale (~458 score points) would dominate the axis. The dashed line marks zero (no effect). These are illustrative synthetic values, not real PISA estimates.", fig.alt = "A horizontal dot-and-interval plot of two slope coefficients from the synthetic fixture fit. The coefficient on x is about 47 score points with a narrow interval well to the right of zero. The coefficient on female is about 2 score points with a wide interval spanning zero from roughly minus 42 to plus 46. A dashed vertical reference line is drawn at zero."----
slopes <- est[est$term != "b_Intercept", ]
slopes <- slopes[order(slopes$term), ]

y    <- seq_len(nrow(slopes))
xlim <- range(c(slopes$conf_low, slopes$conf_high, 0))

op <- par(mar = c(4.5, 7, 1, 1))
plot(
  slopes$estimate, y,
  xlim = xlim, ylim = c(0.5, nrow(slopes) + 0.5),
  yaxt = "n", ylab = "",
  xlab = "Coefficient (synthetic reading-score points)",
  pch = 19, cex = 1.4, col = "#1f6f9c"
)
abline(v = 0, lty = 2, col = "grey50")
segments(slopes$conf_low, y, slopes$conf_high, y, lwd = 2, col = "#1f6f9c")
points(slopes$estimate, y, pch = 19, cex = 1.4, col = "#1f6f9c")
axis(2, at = y, labels = slopes$term, las = 1)
par(op)

## ----honest-columns-----------------------------------------------------------
est[, c("term", "interval_role", "coverage_claim_allowed")]

## ----live-backend-sketch, eval = FALSE----------------------------------------
# # fit_function: fit the stacked model on the prepared (N * M)-row data and return
# # whatever the backend produces (e.g. a brms or cmdstanr fit object).
# my_fit_function <- function(formula, data, weights, ...) {
#   # ... call your Bayesian engine here ...
# }
# 
# # draws_function: return a posterior draws matrix whose columns are the
# # fixed-effect parameters CCC will calibrate (b_* naming, or supply param_map).
# my_draws_function <- function(backend_fit, ...) {
#   # ... extract a draws matrix from backend_fit ...
# }
# 
# fit <- pv_fit(
#   data           = your_real_pisa_data,
#   formula        = OUTCOME ~ x + female,
#   target         = your_target,
#   method         = "stack_direct",
#   control        = pv_control(method = "stack_direct", backend = "injected"),
#   fit_function   = my_fit_function,
#   draws_function = my_draws_function
# )

## ----session-info-------------------------------------------------------------
sessionInfo()

