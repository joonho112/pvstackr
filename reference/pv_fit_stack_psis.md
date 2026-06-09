# Fit the PSIS-Reweighted Stacked Method

`pv_fit_stack_psis()` implements the `stack_psis` API for this package
stage: one stacked draw source, supplied/precomputed or injected
per-plausible-value PSIS weights, Pareto-k diagnostics, model-based
Rubin pooling of PSIS-weighted fixed-effect summaries, and explicit
gating when Pareto-k diagnostics fail.

## Usage

``` r
pv_fit_stack_psis(
  data = NULL,
  formula = NULL,
  pv_cols = NULL,
  control = pv_control(method = "stack_psis"),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  stack_fit = NULL,
  stacked_draws = NULL,
  param_map = NULL,
  psis_weights = NULL,
  pareto_k = NULL,
  log_ratios = NULL,
  psis_function = NULL,
  fallback = c("block", "warn"),
  weight_col = NULL,
  rep_weight_cols = NULL,
  fay_k = 0.5,
  id_cols = NULL,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL,
  allow_m1 = FALSE,
  cache_dir = "cache",
  cache_stem = "pvstackr-stack-psis",
  additional_args = list()
)
```

## Arguments

- data:

  Optional analysis data frame for the injected stacked-fit route.

- formula:

  Optional formula with `OUTCOME` on the left-hand side for the injected
  stacked-fit route.

- pv_cols:

  Plausible-value columns.

- control:

  A
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  object with `method = "stack_psis"`.

- family, prior:

  Optional backend arguments passed to an injected stacked fit.

- fit_function, draws_function:

  Injected stacked backend fit and draw extractor.

- stack_fit:

  Optional existing `pvstackr_stack_fit`.

- stacked_draws:

  Optional stacked posterior draw matrix.

- param_map:

  Optional explicit draw-column map. Supply `fe_names` or `fe_idx` to
  identify fixed-effect columns, and optional `vc_names` or `vc_idx` for
  nuisance variance-component columns. Use `vc_names = character()` to
  drop all nuisance columns. Explicit maps are recommended when backend
  draw names do not follow the automatic `b_*` fixed-effect convention,
  including distributional names such as `b_sigma_*`.

- psis_weights:

  Optional normalized or unnormalized PSIS weight matrix, one column per
  plausible value.

- pareto_k:

  Pareto-k diagnostics aligned to plausible values.

- log_ratios:

  Optional log-ratio matrix used with `psis_function`, or
  self-normalized when `pareto_k` is supplied.

- psis_function:

  Optional function returning `weights` and `pareto_k`.

- fallback:

  Behavior when Pareto-k exceeds the threshold: `"block"` or `"warn"`.

- weight_col, rep_weight_cols, fay_k, id_cols:

  Optional design metadata.

- df_method, df_complete, allow_m1:

  Rubin pooling options.

- cache_dir, cache_stem:

  Cache metadata for injected stacked fits.

- additional_args:

  Additional named arguments passed to injected fits.

## Value

A `pvstackr_fit` object with `method = "stack_psis"`.

## Details

This function does not depend on `loo`, `brms`, or `cmdstanr` directly
and does not run a live PSIS routine from the `loo` package by default.
Users may supply PSIS weights plus Pareto-k values, or inject a
`psis_function(log_ratios)` that returns `weights` and `pareto_k`. In
v0.1, reportable output is fixed-effect-only; group terms such as
`(1 | school)` and `(1 || school)` are rejected.

`stack_psis` interval metadata is diagnostic/reference vocabulary, not a
design-coverage claim. Classic Rubin pooling is labeled
`interval_role = "psis_classic_rubin"`; Barnard-Rubin pooling is labeled
`interval_role = "psis_barnard_rubin"` and uses the supplied
`df_complete`. Both roles set `coverage_claim_allowed = FALSE` because
the pooled covariance is model-based weighted covariance of the stacked
draws. The pooling provenance is labeled `stack_psis_rubin_pooling`.
PSIS diagnostics record their weight source as `supplied_weights`,
`psis_function`, or `log_ratios_self_normalized`; the last path
self-normalizes the log ratios and does not by itself run Pareto
smoothing.

If any Pareto-k exceeds `control$psis_k_threshold`, the default
`fallback = "block"` returns a blocked `pvstackr_fit` with diagnostics
but no reportable estimates. `fallback = "warn"` retains estimates with
warning status for diagnostic comparison workflows.

## Reportable scope and coverage

In this package stage, reportable output is **fixed-effect-only**;
variance components are fit but not calibrated to the target. Coverage
claims are enabled **only** for `stack_direct` rows backed by the
external Rubin/BRR-Fay target
(`interval_role = "coverage_barnard_rubin"`,
`coverage_claim_allowed = TRUE`). `stack_psis` intervals are
**descriptive/reference** even with Barnard-Rubin degrees of freedom.
"One stacked fit" describes the computational topology, not a
benchmarked speed claim. A small Pareto-\\\hat k\\ does **not** certify
correct variance: reusing one draw cloud correlates the per-PV
imputations, so intervals can run narrow; \\\hat k\\ is
specification-dependent. Treat `stack_psis` as a cross-check, never the
deliverable.

## References

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5), 1413-1432.

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72), 1-58.

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.*
Wiley.

## See also

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md);
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).

Other pvstackr-fitting:
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)

## Examples

``` r
set.seed(1)
# Injected stack_psis demo: a stacked draw matrix plus supplied PSIS weights
# and Pareto-k (the function does not run `loo` by default). Fixed-effect
# columns are named like `b_*`; weights have one column per plausible value.
M <- 2L
stacked_draws <- matrix(
  rnorm(400 * 2), ncol = 2,
  dimnames = list(NULL, c("b_Intercept", "b_x"))
)
psis_weights <- matrix(1 / 400, nrow = 400, ncol = M)  # equal weights
pareto_k     <- rep(0.2, M)                             # all "good" (< 0.5)
fit_psis <- pv_fit_stack_psis(
  stacked_draws = stacked_draws,
  pv_cols       = paste0("PV", seq_len(M)),
  psis_weights  = psis_weights,
  pareto_k      = pareto_k,
  control       = pv_control(method = "stack_psis")
)
fit_psis                    # a stack_psis pvstackr_fit
#> pvstackr fit
#>   method: stack_psis
#>   status: ok
#>   fixed effects: 2
#>   target: none
#>   draws: not retained
#>   diagnostics: psis, pooling, weighted
#>   interval note: intervals are descriptive rather than coverage-claimable.
get_estimates(fit_psis)     # interval_role = "psis_*"; coverage_claim_allowed = FALSE
#>          term    estimate       se std.error         df df_method df_complete
#> 1 b_Intercept  0.03808867 4.335938  4.335938 4.5036e+15   classic          NA
#> 2         b_x -0.07074459 4.828052  4.828052 4.5036e+15   classic          NA
#>   conf_level  conf_low conf_high  conf.low conf.high      interval_role
#> 1       0.95 -8.460193  8.536371 -8.460193  8.536371 psis_classic_rubin
#> 2       0.95 -9.533552  9.392063 -9.533552  9.392063 psis_classic_rubin
#>   coverage_claim_allowed parameter_scope            target_source target_hash
#> 1                  FALSE    fixed_effect stack_psis_rubin_pooling    90cfbed1
#> 2                  FALSE    fixed_effect stack_psis_rubin_pooling    90cfbed1
#>             pooling_source pooling_hash psis_status pareto_k_max
#> 1 stack_psis_rubin_pooling     90cfbed1          ok          0.2
#> 2 stack_psis_rubin_pooling     90cfbed1          ok          0.2
#>   psis_k_threshold
#> 1              0.7
#> 2              0.7
```
