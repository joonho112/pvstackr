# Fit the Direct Stacked Plausible-Value Model

`pv_fit_direct()` runs the v0.1 `stack_direct` path: fixed-effect-only
compatibility preflight, one stacked plausible-value fit, CCC
calibration to an external BRR-Fay/Rubin target, and assembly of a
reportable `pvstackr_fit` object.

## Usage

``` r
pv_fit_direct(
  data,
  formula,
  target,
  control = pv_control(method = "stack_direct"),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  param_map = NULL,
  diagnose_function = NULL,
  log_lik_function = NULL,
  extract_log_lik = FALSE,
  cache_dir = "cache",
  cache_stem = "pvstackr-stack-direct",
  additional_args = list()
)
```

## Arguments

- data:

  Analysis data frame.

- formula:

  Two-sided formula with `OUTCOME` on the left-hand side.

- target:

  A `pvstackr_brr_target` object from
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).

- control:

  A
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  object with `method = "stack_direct"`.

- family:

  Optional backend family object passed to the injected fit function.

- prior:

  Optional backend prior object. Explicit priors are allowed but
  reported as warning-level diagnostics because the current identity
  result is scoped to MLE/flat-prior fixed-effect regimes.

- fit_function:

  Injected backend fitting function for this package stage.

- draws_function:

  Function extracting posterior draws from the injected backend fit.

- param_map:

  Optional explicit draw-column map passed to the stack-fit layer.
  Supply `fe_names` or `fe_idx` to identify fixed-effect columns, and
  optional `vc_names` or `vc_idx` for nuisance variance-component
  columns. Use `vc_names = character()` to drop all nuisance columns.
  Explicit maps are recommended when backend draw names do not follow
  the automatic `b_*` fixed-effect convention, including distributional
  names such as `b_sigma_*`.

- diagnose_function:

  Optional backend diagnostic extractor.

- log_lik_function:

  Optional log-likelihood extractor.

- extract_log_lik:

  Whether to extract log-likelihood draws.

- cache_dir, cache_stem:

  Cache location metadata passed to the injected fit function.

- additional_args:

  Additional named arguments passed to the injected fit function.

## Value

A `pvstackr_fit` object.

## Details

In v0.1, `stack_direct` is fixed-effect-only and requires an external
`pvstackr_brr_target` from
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
The formula RHS and derived fixed-effect names must match the target
exactly. Group terms such as `(1 | school)` are rejected until a
two-level target engine is implemented.

The returned `pvstackr_fit` has `status = "ok"`, `"warning"`, or
`"blocked"`. Yellow CCC center separation and explicit priors produce
warning-status fits; red CCC center separation blocks reportable
estimates. CCC calibration-matrix conditioning is also gated by
`kappa_A`: values at or above `1e6` produce warning-status fits, and
values at or above `1e8` block reportable estimates. Reportable
estimates and retained draws are fixed-effect-only.

With the default `control$center = "target"`, the reportable
fixed-effect estimate table is target-based: `estimate` is the external
Rubin/BRR-Fay coefficient after CCC centering, `se` is
`sqrt(diag(target$T_MI))`, and the degrees of freedom and interval
metadata are inherited from `target`. Therefore `df_method`,
`df_complete`, `interval_role`, and `coverage_claim_allowed` describe
the external target policy, not a residual df or interval policy
estimated by the stacked backend fit. When the target uses classic Rubin
df, the reported df are the target's classic Rubin imputation df; when
it uses Barnard-Rubin df, `df_complete` is the explicit complete-data df
supplied to
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
The stacked fit supplies the draw cloud used for calibration, retained
calibrated draws when requested, and center-separation agreement
diagnostics; it is not a replacement source for the headline
Rubin/BRR-Fay fixed-effect numbers.

Reportable `pv_fit_direct()` output requires
`control$center = "target"`. The lower-level CCC convention
`center = "posterior"` is diagnostic and exploratory: it leaves the
output fixed-effect center at the raw stacked posterior mean while
pairing that center with the external target covariance.
`pv_fit_direct()` does not expose that hybrid as a reportable
`stack_direct` estimate table.

See also
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).

## Reportable scope and coverage

In this package stage, reportable output is **fixed-effect-only**;
variance components are fit but not calibrated to the target. Coverage
claims are enabled **only** for `stack_direct` rows backed by the
external Rubin/BRR-Fay target with Barnard-Rubin df and an explicit
`df_complete` (`interval_role = "coverage_barnard_rubin"`,
`coverage_claim_allowed = TRUE`); with classic Rubin df the rows are
descriptive (`interval_role = "descriptive_classic_rubin"`,
`coverage_claim_allowed = FALSE`), as in the bundled cached example.
`per_pv` and `stack_psis` intervals are **descriptive/reference** even
with Barnard-Rubin degrees of freedom. "One stacked fit" describes the
computational topology, not a benchmarked speed claim.

## References

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.*
Wiley.

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948-955.

## See also

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md),
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md);
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).

Other pvstackr-fitting:
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)

## Examples

``` r
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

# Build the design and the external BRR-Fay target on the fixture.
design <- pv_design(
  pisa_tiny, formula = OUTCOME ~ x + female,
  pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
)
target <- pv_brr_target(
  pisa_tiny, formula = OUTCOME ~ x + female,
  pv_cols = design$pv_cols, weight_col = design$weight_col,
  rep_weight_cols = design$rep_weight_cols, fay_k = design$fay_k,
  id_cols = design$id_cols
)

# A live stack_direct fit needs injected/precomputed backend adapters
# (e.g. brms/cmdstanr); not run in checks.
if (FALSE) { # \dontrun{
fit <- pv_fit_direct(
  data = pisa_tiny, formula = OUTCOME ~ x + female,
  target = target,
  fit_function = my_backend_adapter, draws_function = my_draws_adapter
)
} # }

# Inspect the reportable fixed-effect surface via the bundled cached fit.
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
if (nzchar(path)) {
  fit <- readRDS(path)$fit     # a stack_direct pvstackr_fit
  head(get_estimates(fit))     # reportable fixed-effect table
}
#>          term   estimate        se std.error       df df_method df_complete
#> 1 b_Intercept 457.894088 1.2873118 1.2873118 1.021194   classic          NA
#> 2         b_x  46.883361 0.3717929 0.3717929 1.402308   classic          NA
#> 3    b_female   2.143702 3.5550309 3.5550309 1.013730   classic          NA
#>   conf_level  conf_low conf_high  conf.low conf.high             interval_role
#> 1       0.95 442.31804 473.47013 442.31804 473.47013 descriptive_classic_rubin
#> 2       0.95  44.41457  49.35215  44.41457  49.35215 descriptive_classic_rubin
#> 3       0.95 -41.60687  45.89428 -41.60687  45.89428 descriptive_classic_rubin
#>   coverage_claim_allowed parameter_scope          target_source target_hash
#> 1                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
#> 2                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
#> 3                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
```
