# Reweight stacked draws toward each plausible value

`pv_fit_stack_psis()` runs the `"stack_psis"` method of
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md).
It takes one set of posterior draws from a model fitted to the stacked
data (one copy of the data per plausible value) and, for each plausible
value, importance weights for these draws and a Pareto k-hat value,
computed outside pvstackr, for example by Pareto smoothed importance
sampling (PSIS) with the loo package. It reweights the draws toward the
posterior of each plausible value and combines the weighted results with
Rubin's rules. pvstackr checks the weights and the k-hat values but does
not compute them.

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

  A data frame, or `NULL` (default). Required with `fit_function`, which
  receives the stacked data (see `fit_function`). Otherwise it is
  optional: given together with `formula`, it is used only to record the
  design (the fit's `design`), and the plausible-value names must then
  be columns of `data`.

- formula:

  A two-sided formula with the placeholder `OUTCOME` on the left-hand
  side, such as `OUTCOME ~ x + female` (see
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)),
  or `NULL` (default). Required with `fit_function`. A random-effect
  term such as `(1 | school)` or a
  [`weights()`](https://rdrr.io/r/stats/weights.html) term stops the
  function with an error, also when the formula only serves to record
  the design.

- pv_cols:

  A character vector with the names of the plausible values, at least
  two, or `NULL` (default). Required with `fit_function`: the
  plausible-value columns of `data`. Otherwise it is optional and names
  the columns of the weights; if the weights have column names, they
  must be the same. Without `pv_cols`, the column names of the weights
  are used, or `PV1`, `PV2`, ... when they have none.

- control:

  A
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  object with `method = "stack_psis"`. The default is
  `pv_control(method = "stack_psis")`, and `NULL` gives the same. Its
  `psis_k_threshold` sets the k-hat cut-off, and its sampler settings
  and `backend` are only passed on to `fit_function`.

- family, prior:

  Passed unchanged to `fit_function`; pvstackr does not check them.
  Default `NULL`.

- fit_function:

  Your function that fits the model to the stacked data, or `NULL`
  (default). pvstackr stacks `data` (one copy per plausible value, with
  the plausible value in the column `.pvstackr_y`) and calls the
  function once with the arguments `formula`, `data`, `family`, `prior`,
  `chains`, `iter`, `warmup`, `cores`, `seed`, `backend`, `file` and
  `file_refit` (the stacked formula and data, the settings from
  `control` and the cache file; see `cache_dir`), plus the elements of
  `additional_args`. The formula has the right-hand side of `formula`,
  for example `.pvstackr_y | weights(.pvstackr_weight) ~ x + female` for
  `OUTCOME ~ x + female`. The row weight `.pvstackr_weight` is
  `weight_col` divided by its mean and by the number of plausible
  values, or 1/M without `weight_col`; a function that ignores it fits
  an unweighted model. There is no bundled engine for this method, also
  with `pv_control(backend = "brms")`.

- draws_function:

  Your function that takes the object returned by `fit_function` and
  returns its draws: a numeric matrix or data frame with one row per
  draw (at least two), finite values and unique column names. Required
  with `fit_function`.

- stack_fit:

  A stacked-fit object of class `pvstackr_stack_fit` that holds its
  stacked draws, or `NULL` (default). The `stack_fit` component of a fit
  returned by pvstackr keeps no stacked draws and stops the function
  with an error, so in practice give `stacked_draws` or `fit_function`
  instead.

- stacked_draws:

  Posterior draws of a model fitted to the stacked data, computed
  elsewhere, or `NULL` (default): a numeric matrix or data frame with
  one row per draw (at least two), finite values and unique column
  names.

- param_map:

  `NULL` (default) or a named list that says which columns of the
  stacked draws are the fixed effects, by name (`fe_names`) or by
  position (`fe_idx`), as for
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md).
  With `NULL`, the columns whose names start with `b_` are the fixed
  effects; with `stack_fit`, the selection recorded in `stack_fit` is
  used. Only the fixed effects are used, and their column names become
  the terms of the estimate table, so they must start with `b_`. Give a
  list to leave out other columns that start with `b_`, such as
  distributional parameters named `b_sigma_*`. With `fit_function` or
  `stack_fit`, give the columns by name (`fe_names`), because positions
  can point to other columns there.

- psis_weights:

  A numeric matrix of importance weights, one row per stacked draw and
  one column per plausible value, or `NULL` (default). The values must
  be finite and not negative, and every column must have a positive sum;
  pvstackr divides each column by its sum. Give it with `pareto_k`.

- pareto_k:

  The Pareto k-hat values, one per plausible value: a numeric vector
  named by plausible value, or in the order of the weight columns.
  Required with `psis_weights`, and with `log_ratios` when there is no
  `psis_function`; leave it `NULL` (default) with `psis_function`, which
  returns it.

- log_ratios:

  A numeric matrix of log importance ratios for reweighting the stacked
  draws toward the posterior of each plausible value, with one row per
  stacked draw, one column per plausible value and finite values; or
  `NULL` (default). With `psis_function`, pvstackr passes it to that
  function. With `pareto_k` alone, pvstackr turns the ratios into
  normalized weights without smoothing and the fit is blocked
  (`psis_smoothing_not_applied`), so this gives only the diagnostics in
  `get_diagnostics(fit)$psis`.

- psis_function:

  Your function that takes `log_ratios` and returns a list with
  `weights` (a matrix of the same dimensions, with the columns in the
  same order) and `pareto_k`, or `NULL` (default). For example, a
  function that applies the function `psis()` of the loo package to each
  column.

- psis_producer, psis_producer_version:

  The name and the version of the program that produced the weights and
  the k-hat values, such as `"loo"` and
  `as.character(packageVersion("loo"))`: single strings of at most 128
  and 64 bytes, without leading or trailing spaces or control
  characters. Give both or neither; without them the fit is blocked
  (`psis_weight_provenance_incomplete`). pvstackr records them in the
  estimate table and in `get_diagnostics(fit)$psis` as your statement of
  how the weights were made, but it does not check them. They cannot be
  given with weights that pvstackr makes from `log_ratios` without
  `psis_function`.

- fallback:

  `"block"` (default) or `"warn"`. Both block a fit that fails the check
  (see Details). `"warn"` is kept so that older code still runs: it
  gives a deprecation warning and is only recorded, as
  `fallback_requested` in `get_diagnostics(fit)$psis`.

- weight_col, rep_weight_cols, fay_k, id_cols:

  The final weight column, the replicate-weight columns, the Fay
  coefficient and the row-identifier columns, as in
  [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md).
  `weight_col` sets the row weights for `fit_function` (see there). With
  `data` and `formula`, all four are checked and recorded in the fit's
  `design`; the replicate weights and `fay_k` are not used in the
  calculation. Defaults `NULL`, `NULL`, `0.5` and `NULL`.

- df_method:

  The rule for the degrees of freedom: `"classic"` (default) or
  `"barnard_rubin"`, which needs `df_complete`. Both rules are defined
  in
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
  With either rule the intervals are descriptive.

- df_complete:

  For `df_method = "barnard_rubin"`, the complete-data degrees of
  freedom, which you state: one positive number for all fixed effects,
  or one value per fixed effect, preferably named by the fixed-effect
  columns (such as `b_Intercept`); an unnamed vector is taken in the
  order of the fixed-effect columns. Leave it `NULL` (default) with
  `"classic"`: a value there changes nothing but the `df_complete`
  column of the estimate table.

- allow_m1:

  Has no effect: a `stack_psis` fit needs weights and k-hat values for
  at least two plausible values. Default `FALSE`.

- cache_dir, cache_stem:

  The folder and the file name that are passed to `fit_function`;
  defaults `"cache"` and `"pvstackr-stack-psis"`. It receives
  `file = file.path(path.expand(cache_dir), cache_stem)` and
  `file_refit = "on_change"`, or `file = NULL` and
  `file_refit = "never"` when `cache_dir = NULL`. pvstackr does not
  create the folder.

- additional_args:

  A named list of further arguments passed to `fit_function`,
  [`list()`](https://rdrr.io/r/base/list.html) by default. They may not
  repeat the arguments that pvstackr sets (see `fit_function`).

## Value

A `pvstackr_fit` object with `method = "stack_psis"`. Read it with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md);
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
and
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
return `NULL` for this method, and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
describes every part of the object. The estimate table has one row per
fixed effect with the combined estimate, its standard error and its
degrees of freedom; besides the columns of every method, it has
`pooling_source`, `pooling_hash` and the columns `psis_status`,
`pareto_k_max`, `psis_k_threshold`, `psis_source`, `pareto_k_source`,
`weight_method`, `psis_producer` and `psis_producer_version`.
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
returns `psis` (the k-hat values and their check, the sources of the
weights and the weight-concentration values), `pooling` (the
Rubin's-rules combination, such as `U_bar`, `B` and `T_MI`) and
`weighted` (the weighted estimate and covariance of each plausible value
and, with `control$return_draws = TRUE`, the fixed-effect draws and the
normalized weights). For a blocked fit the estimate table is empty, and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
returns only `psis` and `redaction`, which lists what was removed.

## Details

### Inputs

The stacked draws come from exactly one of `stacked_draws`,
`fit_function` with `draws_function`, or `stack_fit`. The weights come
from `psis_weights` with `pareto_k`, or from your `psis_function`, which
pvstackr calls on `log_ratios` and which returns both. The weights have
one row per stacked draw and one column per plausible value. Name the
program that produced them with `psis_producer` and
`psis_producer_version`: the numbers alone cannot show that the weights
were Pareto-smoothed.

PSIS stabilizes importance weights with a generalized Pareto
distribution fitted to the upper tail of the importance ratios; the
Pareto k-hat, the estimated shape parameter of that distribution, shows
how reliable the reweighted estimates are (Vehtari et al. 2024).
pvstackr itself does not run Pareto smoothing or any other part of PSIS:
it fits no Pareto distribution and estimates no k-hat.

### Calculation

For each plausible value \\m\\, pvstackr divides the weights by their
sum and computes the weighted mean \\\hat\beta_m\\ and the weighted
covariance matrix \\U_m\\ of the fixed-effect draws (see `param_map`),
with the covariance divided by \\1 - \sum w^2\\. Rubin's rules then
combine the \\M\\ results as in
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md):
the total covariance is \\T\_{\mathrm{MI}} = \bar U + (1 + 1/M) B\\, the
standard errors are the square roots of its diagonal, and the degrees of
freedom follow `df_method`.

### Checks and status

By pvstackr's rule, the fit is blocked (status `"blocked"`, with no
estimates, weights or draws) unless every plausible value has a finite
k-hat strictly below `control$psis_k_threshold` (default 0.7, the
largest value that
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
accepts) and the weights come with `psis_producer` and
`psis_producer_version`. The reason code says which part failed:
`psis_k_not_evaluated` (a k-hat is missing or not finite),
`psis_k_too_high`, `psis_weight_provenance_incomplete` (no program
named) or `psis_smoothing_not_applied` (weights that pvstackr made from
`log_ratios` without smoothing). Otherwise the status is `"ok"`; a
`stack_psis` fit never has status `"warning"`, and pvstackr does not
check the sampler diagnostics of the stacked fit.

Vehtari et al. (2024) recommend the cut-off \\\min(1 - 1/\log\_{10} S,
0.7)\\ for \\S\\ draws, for example 0.67 for 1000 draws; to use it, set
`psis_k_threshold` in
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).

### What is reported

As for every method, only the fixed effects are reported
([`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
gives the scope and the full rule for the intervals). The variance
within each plausible value comes from the weighted posterior draws, not
from replicate weights, so every row of the estimate table has
`interval_role = "psis_classic_rubin"` or `"psis_barnard_rubin"` and
`coverage_claim_allowed = FALSE`. A k-hat below the cut-off lets
pvstackr report the estimates, but the intervals remain descriptive.
pvstackr provides `stack_psis` as a check on a stacked fit, not as a
replacement for `stack_direct`.

`get_diagnostics(fit)$psis` also records, for a blocked fit too, how
concentrated the normalized weights are, and
`weight_diagnostic_authority` says whether the weights were kept
(`"retained_weights_recomputed"`, for a fit with estimates and
`control$return_draws = TRUE`) or removed
(`"owned_stamp_bounded_projection"`). These values do not change the
status;
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
explains them.

## References

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72), 1-58.

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.*
Wiley.

## See also

[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
to compare fits,
[`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
for `stack_psis` fits saved by earlier versions of pvstackr, and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
for the parts of a fit.

Other pvstackr-fitting:
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)

## Examples

``` r
set.seed(1)
# Simulated stacked draws with the fixed-effect columns b_Intercept and
# b_x, equal placeholder weights with one column per plausible value, and
# k-hat values below the 0.7 cut-off. These weights were not made by
# PSIS, so no program is named, and pvstackr blocks the fit.
M <- 2L
stacked_draws <- matrix(
  rnorm(400 * 2), ncol = 2,
  dimnames = list(NULL, c("b_Intercept", "b_x"))
)
psis_weights <- matrix(1 / 400, nrow = 400, ncol = M)  # equal weights
pareto_k     <- rep(0.2, M)                             # below 0.7
fit_psis <- pv_fit_stack_psis(
  stacked_draws = stacked_draws,
  pv_cols       = paste0("PV", seq_len(M)),
  psis_weights  = psis_weights,
  pareto_k      = pareto_k,
  control       = pv_control(method = "stack_psis")
)
# Status "blocked", reason code psis_weight_provenance_incomplete:
fit_psis
#> pvstackr fit
#>   method: stack_psis
#>   status: blocked
#>   fixed effects: 0
#>   target: none
#>   draws: not retained
#>   diagnostics: psis, redaction
#>   reason codes: psis_weight_provenance_incomplete
#>   warnings: 1
# The PSIS check has status "provenance_incomplete", and weight_method is
# "unspecified_external" because no program was named:
get_diagnostics(fit_psis)$psis[c("status", "weight_method")]
#> $status
#> [1] "provenance_incomplete"
#> 
#> $weight_method
#> [1] "unspecified_external"
#> 
```
