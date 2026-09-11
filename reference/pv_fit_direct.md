# Fit a stacked model and calibrate it to a BRR-Fay target

`pv_fit_direct()` runs the `"stack_direct"` method of
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md).
It fits one model to the stacked data (one copy of `data` per plausible
value), checks the sampler diagnostics of that fit, and then applies the
Cholesky calibration correction (CCC): the fixed-effect draws are
transformed so that their mean and covariance equal the estimates and
the covariance (`T_MI`) of `target`. The reported estimates, standard
errors and degrees of freedom are those of `target`.

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

  The data frame used to build `target`: the same rows in the same
  order, with the same plausible values, weights and covariates.
  Otherwise `pv_fit_direct()` stops with an error.

- formula:

  The formula used to build `target`, with the placeholder `OUTCOME` on
  the left-hand side, such as `OUTCOME ~ x + female` (see
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)).
  Random-effect terms such as `(1 | school)` stop with an error.

- target:

  The `pvstackr_brr_target` object returned by
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  for `data` and `formula`.

- control:

  A
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  object with `method = "stack_direct"`. The default is
  `pv_control(method = "stack_direct")`, and `NULL` gives the same. Its
  `center` must be `"target"`.

- family:

  `NULL` (default), `"gaussian"`, or a family object or list with
  elements `family` and `link`, such as
  [`gaussian()`](https://rdrr.io/r/stats/family.html). Only the Gaussian
  family with the identity link is accepted; any other family or link
  stops with an error. The fitting function always receives
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html) made by
  pvstackr, not the object you pass.

- prior:

  `NULL` (default) or a prior table: a data frame with character columns
  `class` and `coef`, such as the result of the brms function
  `set_prior()`. Only priors of class `"b"` or `"sigma"` with an empty
  `coef` are accepted; a prior on the intercept, on one coefficient or
  of another class stops with an error. A `"b"` prior is given to the
  slopes only, as in the original formula, so the intercept keeps the
  default prior. Any explicit prior, even a flat one, gives the fit the
  status `"warning"` (`explicit_prior_warning`).

- fit_function:

  `NULL` (default) or your function that fits the stacked model. `NULL`
  selects the bundled brms engine when `control$backend` is `"brms"` and
  stops with an error otherwise. The function is called with the
  arguments `formula`, `data`, `family`, `prior`, `chains`, `iter`,
  `warmup`, `cores`, `seed`, `backend`, `file` and `file_refit` (the
  stacked formula and data, the settings from `control` and the cache
  file; see Details and `cache_dir`), plus the elements of
  `additional_args`. There is no `weights` argument: the weights are in
  `formula` and in the column `.pvstackr_weight` of `data`, and a
  function that ignores them fits an unweighted model. It may return any
  object that `draws_function` and `diagnose_function` accept.

- draws_function:

  Your function that takes the object returned by the fitting function
  and returns its draws: a numeric matrix or data frame with one row per
  draw (at least two), finite values and unique column names. It may be
  `NULL` (default) only with the bundled engine. The fixed-effect
  columns must be named as in `target$fe_names` (such as `b_Intercept`)
  or after the columns of the stacked formula (`b_pvstackrMM001`,
  `b_pvstackrMM002`, ...); otherwise give `param_map` and list all fixed
  effects in the order of `target$fe_names`.

- param_map:

  `NULL` (default) or a named list that says which columns of the draws
  are the fixed effects. With `NULL`, columns whose names start with
  `b_` are the fixed effects; `sigma` and `sd_*` columns are other model
  parameters, which are not reported, and all other columns, such as
  `lp__`, are dropped. Give a list with `fe_names` (names) or `fe_idx`
  (positions) when the fixed-effect columns have other names or when
  other columns start with `b_`, such as distributional parameters named
  `b_sigma_*`. Its optional elements `vc_names` and `vc_idx` select the
  other model parameters and change nothing that is reported.

- diagnose_function:

  Your function that takes the object returned by the fitting function
  and returns its sampler diagnostics as a named list with the elements
  `rhat_max` (the largest R-hat), `ess_bulk_min` and `ess_tail_min` (the
  smallest bulk and tail effective sample sizes, totals over all
  chains), `divergences` (the number of divergent transitions), `chains`
  and `post_warmup_draws_per_chain`; the list may also sit in an element
  `sampler`. For the per-chain check, pvstackr divides each effective
  sample size by `chains`; a returned `ess_bulk_per_chain_min` or
  `ess_tail_per_chain_min` must equal that ratio. If the function is
  missing or fails, or a value is missing or invalid, the fit is blocked
  (`sampler_diagnostics_incomplete`), and `diagnostic_reason_codes` in
  `get_diagnostics(fit)$sampler` says why. The bundled engine ignores
  this argument: it takes R-hat and the effective sample sizes over all
  variables of the fit from the posterior package and the number of
  divergent transitions from the brms sampler record.

- log_lik_function:

  `NULL` (default) or your function that takes the object returned by
  the fitting function and returns the pointwise log-likelihood: a
  finite numeric matrix with one row per draw and one column per row of
  the stacked data. It is called only when `extract_log_lik = TRUE`; the
  bundled engine does not supply one.

- extract_log_lik:

  Whether to call `log_lik_function`, which is then required: `TRUE` or
  `FALSE` (default). The matrix is kept in the fit's `stack_fit` only
  with `pv_control(keep_log_lik = TRUE)`, for your own use: no pvstackr
  function uses it. `keep_log_lik = TRUE` stops with an error unless
  `extract_log_lik = TRUE`.

- cache_dir, cache_stem:

  The folder and the file name (without extension) for a cached stacked
  fit; defaults `"cache"` and `"pvstackr-stack-direct"`. `cache_stem`
  must be a plain file name, without a path separator. For the bundled
  brms engine, pvstackr creates the folder, verifies that it is
  writable, and uses brms `file_refit = "on_change"`. Set
  `cache_dir = NULL` to disable file caching. Your own `fit_function`
  receives `file = file.path(path.expand(cache_dir), cache_stem)` and
  `file_refit = "on_change"`, or `file = NULL` and
  `file_refit = "never"` when `cache_dir = NULL`; pvstackr does not
  create the folder for it.

  brms refits a cached model only when its Stan code or data change, not
  when `chains`, `iter`, `warmup` or `seed` change. The sampler check
  then blocks a cached fit whose number of chains or of draws per chain
  differs from `control`, and a changed seed returns the cached draws
  without a warning. Use `cache_dir = NULL` or a new `cache_stem` when
  you vary these settings, as in a seed-stability check.

- additional_args:

  A named list of further arguments for the fitting function,
  [`list()`](https://rdrr.io/r/base/list.html) by default. They are
  passed to your `fit_function` or, with the bundled engine, to the brms
  function `brm()`, for example
  `list(control = list(adapt_delta = 0.95))`. They may not repeat the
  arguments that pvstackr sets (see `fit_function`).

## Value

A `pvstackr_fit` object with `method = "stack_direct"`. Read it with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md);
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
explains the status and the estimate table, and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
describes every part of the object. The estimate table has one row per
fixed effect, with the estimates, standard errors and degrees of freedom
of `target`.
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
returns the calibrated fixed-effect draws (one row per draw of the
stacked fit, one column per fixed effect) when `control$return_draws` is
`TRUE`, the default.
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
returns `preflight` (the check that `data`, `formula` and `target`
match), `sampler` and `sampler_gate` (the sampler diagnostics and their
check), `stack_fit`, `stack_fit_warnings` (notes from the stacked fit,
such as dropped draw columns) and `ccc` (the calibration diagnostics).
For a blocked fit it returns only `preflight`, `sampler`,
`sampler_gate`, `redaction` (what was removed) and, when the calibration
check blocked the fit, `ccc`, with its values grouped as `center`,
`conditioning`, `residual` and `prior` (for example
`get_diagnostics(fit)$ccc$center$delta_c_max`).

## Details

### The stacked fit

The stacked data hold one copy of `data` for each plausible value, with
that plausible value as the outcome. Each row is weighted by the final
survey weight named in `target`, divided by its mean and by the number
of plausible values, and the model is a Gaussian linear model with the
identity link. pvstackr adds the columns of the design matrix of
`formula`, the intercept included, to the stacked data as ordinary
covariates, so the fitting function receives the formula
`.pvstackr_y | weights(.pvstackr_weight) ~ 0 + pvstackrMM001 + ...`, in
which the weights enter through
[`weights()`](https://rdrr.io/r/stats/weights.html). With the bundled
brms engine and `prior = NULL`, brms therefore gives every fixed effect,
the intercept included, its default flat prior. The residual standard
deviation `sigma` is estimated but neither calibrated nor reported.

### What is reported

Only the fixed effects are reported. By pvstackr's reporting rule, the
intervals of a `stack_direct` fit can be read as confidence intervals
with nominal coverage only when `target` uses Barnard-Rubin degrees of
freedom;
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
gives the full rule.

`pv_fit_direct()` requires `control$center = "target"` (the default of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)),
so the estimate table holds the numbers of the target: `estimate` is the
target estimate, `se` is `sqrt(diag(target$T_MI))`, and `df`,
`df_method`, `df_complete`, `interval_role` and `coverage_claim_allowed`
are copied from `target`: the degrees of freedom follow the rule chosen
in
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
and are not residual degrees of freedom estimated from the stacked fit.
The stacked fit supplies the draws that are calibrated and the
diagnostics described below. The other value, `center = "posterior"`,
would keep the mean of the stacked draws and calibrate only their
covariance; it is accepted by
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
but `pv_fit_direct()` stops with an error when it is set.

### Fitting engine

With `pv_control(backend = "brms")` and `fit_function = NULL`, the brms
engine bundled with pvstackr fits the stacked model. It needs the brms
and posterior packages, and
[`pv_backend_brms_fit_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
describes how it chooses the Stan backend. A `draws_function` that you
pass replaces its draws step, but a `diagnose_function` is ignored: the
bundled engine computes the sampler diagnostics itself. Otherwise, as
with the default `backend = "none"`, pass all three functions. A
`fit_function` that you pass is used whatever the `backend` value, which
is only passed on to it. There is no argument for draws computed
elsewhere.

### Checks and status

Before the model is fitted, `pv_fit_direct()` stops with an error when
`data` or `formula` differ from those used to build `target`, when
`formula` has a random-effect term such as `(1 | school)`, when `family`
is not Gaussian with the identity link, or when `prior` is not accepted.
After the stacked fit, pvstackr's checks give the fit the status `"ok"`,
`"warning"` or `"blocked"`:

- Sampler diagnostics, checked first. The six values listed under
  `diagnose_function`, from the bundled engine or from your function,
  must all be available, and `chains` and `post_warmup_draws_per_chain`
  must equal `chains` and `iter - warmup` in `control`. The largest
  R-hat and the smallest bulk and tail effective sample sizes can give a
  warning or a block, and any divergent transition blocks the fit. A fit
  blocked here is not calibrated.

- Calibration. `delta_c_max`, the largest difference between a target
  estimate and the mean of the stacked draws in target standard errors,
  and `kappa_A`, the condition number of the calibration matrix, can
  each give a warning or a block.

- Priors. Any explicit `prior` gives a warning.

The thresholds are listed under "Status and checks" in
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).
The status is `"ok"` when no check gave a warning or a block; it does
not show that the sampler converged. A fit with status `"warning"` has
estimates, and `reason_codes` and `warnings` name the checks. A blocked
fit has an empty estimate table and no draws; it keeps the target and
the diagnostics that explain the block.

## References

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.*
Wiley.

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948-955.

## See also

[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
to build the target,
[`pv_backend_brms_fit_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
for the bundled brms engine, and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
for the parts of a fit and the thresholds of the checks.

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

# A live fit samples with Stan, so it is not run here. With
# backend = "brms" and no fit_function, the bundled brms engine is used.
if (FALSE) { # \dontrun{
fit <- pv_fit_direct(
  data = pisa_tiny, formula = OUTCOME ~ x + female, target = target,
  control = pv_control(method = "stack_direct", backend = "brms")
)

# With the default backend = "none", fit_function and draws_function are
# required, and without diagnose_function the fit is blocked.
fit <- pv_fit_direct(
  data = pisa_tiny, formula = OUTCOME ~ x + female, target = target,
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
  head(get_estimates(fit))     # the target's estimates, se and df
}
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
