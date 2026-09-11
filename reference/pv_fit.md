# Fit a model to plausible-value data

`pv_fit()` fits a regression model to assessment data whose outcome is
given as plausible values, as in PISA, and returns estimates, standard
errors and intervals for the fixed effects. Each method has its own
function, whose help page lists the arguments that you pass through
`...`:
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
for `"stack_direct"`,
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)
for `"per_pv"` and
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)
for `"stack_psis"`.

## Usage

``` r
pv_fit(
  data,
  formula,
  target = NULL,
  method = "stack_direct",
  control = NULL,
  ...
)
```

## Arguments

- data:

  A data frame with the plausible-value columns, the weight columns and
  the variables in `formula`, or `NULL` when the method starts from
  precomputed draws. For `"stack_direct"` it must be the data used to
  build `target` (the same rows in the same order, with the same
  plausible values, weights and covariates); otherwise `pv_fit()` stops
  with an error.

- formula:

  A two-sided formula with the placeholder `OUTCOME` on the left-hand
  side, for example `OUTCOME ~ x + female`. `OUTCOME` stands for the
  plausible values, so the model is written once. Use `NULL` when the
  method starts from precomputed draws. `"stack_direct"` and
  `"stack_psis"` stop with an error on random-effect terms such as
  `(1 | school)`. For `"stack_direct"` it must be the formula used to
  build `target`.

- target:

  For `"stack_direct"`, the `pvstackr_brr_target` object returned by
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md);
  `pv_fit()` stops with an error if it is `NULL`. Ignored by `"per_pv"`
  and `"stack_psis"`. Default `NULL`.

- method:

  The fitting method: `"stack_direct"` (default), `"per_pv"` or
  `"stack_psis"`. See Details.

- control:

  A
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  object, or `NULL` (default) to use `pv_control(method = method)`. That
  default selects no fitting engine, so `"stack_direct"` stops with an
  error unless you set `pv_control(backend = "brms")` or pass your own
  `fit_function` (see Details). `pv_fit()` stops with an error if
  `control$method` differs from `method`.

- ...:

  Further arguments passed unchanged to the method's function (see
  Details).

## Value

A `pvstackr_fit` object. Read it with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
rather than with `$`; they stop with an error if the object was changed
after it was created. Its `status` is `"ok"`, `"warning"` (estimates are
returned and `warnings` says why) or `"blocked"` (the estimate table is
empty); `reason_codes` gives the reasons for a warning or a block. A
blocked fit does not stop with an error; in code, read the status as
`summary(fit)$status` and the reasons as `summary(fit)$reason_codes`
([`summary()`](https://rdrr.io/r/base/summary.html) checks the fit as
the reading functions do).
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
returns one row per fixed effect with columns such as `term`,
`estimate`, `se`, `df`, `conf_low`, `conf_high`, `interval_role` and
`coverage_claim_allowed`.

## Details

### Methods

- `"stack_direct"` (default) fits one model to the stacked data (one
  copy of the data per plausible value), then applies the Cholesky
  calibration correction (CCC): the fixed-effect draws are transformed
  so that their mean and covariance equal those of `target`. The target
  comes from
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
  which fits a survey-weighted regression to each plausible value,
  estimates the sampling covariance of the coefficients with the BRR-Fay
  replicate weights, and combines the results with Rubin's rules.

- `"per_pv"` fits one model per plausible value and combines the results
  with Rubin's rules.

- `"stack_psis"` reweights one stacked fit toward each plausible value
  with importance weights and Pareto k-hat values that you supply
  (pvstackr checks them but does not compute them), then combines the
  reweighted results with Rubin's rules.

Arguments in `...` are passed unchanged to the method's function.
Besides `data` and `formula`, which a method needs when it fits a model,
each method needs:

- `"stack_direct"`: `target` and a model-fitting engine. Either set
  `pv_control(backend = "brms")`, which fits the stacked model with brms
  functions bundled with pvstackr (this needs the brms and posterior
  packages; if cmdstanr is installed, CmdStan must be configured or the
  fit stops with an error; without cmdstanr, brms uses rstan), or pass
  your own `fit_function`, `draws_function` and `diagnose_function`
  (without `diagnose_function` the fit is blocked). The plausible-value
  and weight columns are taken from `target`. See
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
  also for `cache_dir`.

- `"per_pv"`: `pv_cols` with your own `fit_function` and
  `draws_function`, which are called once per plausible value and do not
  receive the survey weights automatically, or `per_pv_draws` computed
  elsewhere. There is no bundled engine for this method. See
  [`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md).

- `"stack_psis"`: one source of stacked draws (`stacked_draws`,
  `stack_fit`, or `fit_function` with `draws_function` and `pv_cols`);
  importance weights with Pareto k-hat values (`psis_weights` and
  `pareto_k`, or your `psis_function` applied to `log_ratios`); and
  `psis_producer` with `psis_producer_version`, which name the program
  that produced the weights. Without the last two, or when a k-hat value
  is not below `control$psis_k_threshold` (default 0.7), the fit is
  blocked. See
  [`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md).

When a method starts from draws computed elsewhere, pass `data = NULL`
and `formula = NULL` explicitly; `pv_fit()` has no defaults for them.

### What is reported

Only the fixed effects are reported: the intercept and the slope
coefficients. Other model parameters, such as the residual standard
deviation, are not reported, and `"stack_direct"` does not calibrate
them. `"stack_direct"` requires `control$center = "target"` (the
default): its reported estimates, standard errors and degrees of freedom
are those of the target, and the stacked fit supplies the calibrated
draws and the diagnostics that compare the fit with the target.

Each row of the estimate table has an interval (`conf_low`, `conf_high`)
at the level `conf_level` of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
(default 0.95): the estimate plus or minus a t quantile with `df`
degrees of freedom times `se`. `coverage_claim_allowed` says whether
pvstackr's reporting rule lets you read the interval as a confidence
interval with nominal coverage (`TRUE`) or labels it descriptive
(`FALSE`), and `interval_role` names the case. The label records how the
interval was built; it does not certify its coverage. pvstackr sets both
columns by this rule:

- `"stack_direct"` with a target built with Barnard-Rubin degrees of
  freedom (`df_method = "barnard_rubin"` and a `df_complete` value in
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)):
  `interval_role = "coverage_barnard_rubin"` and
  `coverage_claim_allowed = TRUE`. These are the only intervals that the
  rule lets you read as confidence intervals with nominal coverage.

- `"stack_direct"` with the classic Rubin degrees of freedom
  (`df_method = "classic"`, the default): `"descriptive_classic_rubin"`
  and `FALSE`. The bundled example fit, computed from synthetic data, is
  of this kind.

- `"per_pv"` (`"reference_classic_rubin"`, `"reference_barnard_rubin"`)
  and `"stack_psis"` (`"psis_classic_rubin"`, `"psis_barnard_rubin"`):
  always `FALSE`, because their variance within each plausible value
  comes from the model's posterior draws, not from the replicate
  weights.

These labels come from the target and the method, not from the status: a
fit with status `"warning"` keeps them, and a blocked fit has no rows.

## See also

[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
to build the target,
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
to compare fits, and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
for the parts of a fit.

Other pvstackr-fitting:
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)

## Examples

``` r
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

# Declare the columns and build the BRR-Fay target from the bundled
# synthetic data.
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

# Fit with the bundled brms engine. This samples with Stan, so it is not
# run here.
if (FALSE) { # \dontrun{
fit <- pv_fit(
  data = pisa_tiny, formula = OUTCOME ~ x + female,
  target = target, method = "stack_direct",
  control = pv_control(method = "stack_direct", backend = "brms")
)

# With your own engine, pass all three functions (see ?pv_fit_direct).
# Without diagnose_function the fit is blocked.
fit <- pv_fit(
  data = pisa_tiny, formula = OUTCOME ~ x + female,
  target = target, method = "stack_direct",
  fit_function = my_fit_function, draws_function = my_draws_function,
  diagnose_function = my_diagnose_function, cache_dir = NULL
)
} # }

# Read the example fit that ships with pvstackr. It was made from the same
# data with fitting functions that return draws around the target, not
# with a sampler.
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
if (nzchar(path)) {
  fit <- readRDS(path)$fit     # a stack_direct fit (class pvstackr_fit)
  print(fit)                   # method, status and interval note
  head(get_estimates(fit))     # the fixed-effect estimate table
}
#> pvstackr fit
#>   method: stack_direct
#>   status: ok
#>   fixed effects: 3
#>   target: external_brr_fay_rubin
#>   draws: not retained
#>   diagnostics: preflight, sampler, sampler_gate, stack_fit, stack_fit_warnings, ccc
#>   interval note: intervals are descriptive rather than coverage-claimable.
#>          term   estimate        se std.error       df df_method df_complete
#> 1 b_Intercept 457.894088 1.2873118 1.2873118 1.021194   classic          NA
#> 2         b_x  46.883361 0.3717929 0.3717929 1.402308   classic          NA
#> 3    b_female   2.143702 3.5550309 3.5550309 1.013730   classic          NA
#>   conf_level  conf_low conf_high  conf.low conf.high             interval_role
#> 1       0.95 442.31804 473.47013 442.31804 473.47013 descriptive_classic_rubin
#> 2       0.95  44.41457  49.35215  44.41457  49.35215 descriptive_classic_rubin
#> 3       0.95 -41.60687  45.89428 -41.60687  45.89428 descriptive_classic_rubin
#>   coverage_claim_allowed parameter_scope          target_source
#> 1                  FALSE    fixed_effect external_brr_fay_rubin
#> 2                  FALSE    fixed_effect external_brr_fay_rubin
#> 3                  FALSE    fixed_effect external_brr_fay_rubin
#>                                                               target_hash
#> 1 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
#> 2 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
#> 3 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
```
