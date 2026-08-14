# Fit the PSIS-Reweighted Stacked Method

`pv_fit_stack_psis()` implements the `stack_psis` API for this package
stage: one stacked draw source, externally produced per-plausible-value
weights and Pareto-k diagnostics, model-based Rubin pooling of weighted
fixed-effect summaries, and explicit provenance and diagnostic gates.

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
  psis_producer = NULL,
  psis_producer_version = NULL,
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

  Optional normalized or unnormalized external weight matrix, one column
  per plausible value. The numeric matrix alone does not establish that
  PSIS smoothing occurred.

- pareto_k:

  Pareto-k diagnostics aligned to plausible values.

- log_ratios:

  Optional log-ratio matrix used with `psis_function`, or
  self-normalized when `pareto_k` is supplied.

- psis_function:

  Optional function returning `weights` and `pareto_k`.

- psis_producer, psis_producer_version:

  Optional bounded scalar strings identifying the external PSIS producer
  and its version. Both are required together for supplied or injected
  weights to be reportable. They record a caller declaration;
  package-verified `loo` execution is deferred.

- fallback:

  Requested behavior when Pareto-k fails: `"block"` or the legacy
  `"warn"`. Both choices now fail closed; `"warn"` cannot make failed
  PSIS output reportable.

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
Users may supply weights plus Pareto-k values, or inject a
`psis_function(log_ratios)` that returns `weights` and `pareto_k`.
Numeric weights or a function name cannot prove that Pareto smoothing
occurred. Reportable output therefore also requires both `psis_producer`
and `psis_producer_version`; these fields are caller-declared
provenance, not package verification. In v0.2, reportable output is
fixed-effect-only; group terms such as `(1 | school)` and
`(1 || school)` are rejected.

`stack_psis` interval metadata is diagnostic/reference vocabulary, not a
design-coverage claim. Classic Rubin pooling is labeled
`interval_role = "psis_classic_rubin"`; Barnard-Rubin pooling is labeled
`interval_role = "psis_barnard_rubin"` and uses the supplied
`df_complete`. Both roles set `coverage_claim_allowed = FALSE` because
the pooled covariance is model-based weighted covariance of the stacked
draws. The pooling provenance is labeled `stack_psis_rubin_pooling`.
Diagnostics separate the input route (`supplied_psis_weights`,
`injected_psis_function`, or `self_normalized_log_ratios`), the Pareto-k
source, and the weight method. The self-normalized path does not run
Pareto smoothing and is always blocked from estimates. Supplied or
injected weights without producer/version provenance are likewise
diagnostic-only.

Every path records per-PV weight-concentration diagnostics after column
normalization: `weight_ess_iid = 1 / sum(w^2)`, its draw-count fraction,
and the largest normalized weight. This is a Kish-style iid weight
diagnostic, not MCMC ESS, the `loo` relative-efficiency diagnostic, or
an autocorrelation-adjusted PSIS effective sample size. These quantities
explain weight concentration but do not relax or replace the Pareto-k
gate.

Reportability uses an immutable fail-closed gate: every plausible value
must have a finite Pareto-k strictly below `control$psis_k_threshold`,
and that threshold cannot exceed `0.7`. Failed, unevaluated, unsmoothed,
or provenance-incomplete input always returns a blocked `pvstackr_fit`
with no reportable estimates, raw weights, or draws; bounded ESS and
maximum-weight diagnostics remain available. The legacy
`fallback = "warn"` argument remains accepted for call compatibility and
emits a deprecation warning, but it cannot relax the gate and is
recorded only as the requested policy.

When `control$return_draws = TRUE`, retained normalized weights are used
to recompute the ESS, maximum-weight, and weighted-summary diagnostics
during deep validation
(`weight_diagnostic_authority = "retained_weights_recomputed"`). When
weights are intentionally redacted, these bounded diagnostics are
protected by the package-owned validation stamp and cross-field
feasibility checks (`"owned_stamp_bounded_projection"`), but the
original weight vector cannot be independently reconstructed from the
compact object. This is an explicit portability and threat-model
boundary, not a stronger statistical claim.

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
# Fail-closed stack_psis demo: equal placeholder weights and arbitrary
# Pareto-k values do not establish that PSIS ran. This example deliberately
# omits producer/version, so it remains diagnostic-only.
# Fixed-effect
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
fit_psis                    # blocked: provenance_incomplete
#> pvstackr fit
#>   method: stack_psis
#>   status: blocked
#>   fixed effects: 0
#>   target: none
#>   draws: not retained
#>   diagnostics: psis, redaction
#>   reason codes: psis_weight_provenance_incomplete
#>   warnings: 1
get_diagnostics(fit_psis)$psis[c("status", "weight_method")]
#> $status
#> [1] "provenance_incomplete"
#> 
#> $weight_method
#> [1] "unspecified_external"
#> 
```
