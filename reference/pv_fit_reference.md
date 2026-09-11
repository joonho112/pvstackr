# Combine one fit per plausible value with Rubin's rules

`pv_fit_reference()` runs the `"per_pv"` method of
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md).
It takes the posterior draws of one model fit per plausible value,
either computed elsewhere (`per_pv_draws`) or returned by your own
fitting functions (`fit_function` and `draws_function`, called once per
plausible value), and combines the fixed effects with Rubin's rules. The
variance within each plausible value is the covariance of that fit's
posterior draws, so the intervals are always descriptive.

## Usage

``` r
pv_fit_reference(
  data = NULL,
  formula = NULL,
  pv_cols = NULL,
  per_pv_draws = NULL,
  control = pv_control(method = "per_pv"),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  param_map = NULL,
  weight_col = NULL,
  rep_weight_cols = NULL,
  fay_k = 0.5,
  id_cols = NULL,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL,
  cache_dir = "cache",
  cache_stem = "pvstackr-per-pv",
  additional_args = list()
)
```

## Arguments

- data:

  A data frame with the plausible-value columns and the variables in
  `formula`, or `NULL` (default). Required with `fit_function`, which
  receives it unchanged. With `per_pv_draws` it is optional: given
  together with `formula`, it is used only to record the design (the
  fit's `design`), and the plausible-value names must then be columns of
  `data`.

- formula:

  A two-sided formula with the placeholder `OUTCOME` on the left-hand
  side, such as `OUTCOME ~ x + female` (see
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)),
  or `NULL` (default). Required with `fit_function`. A
  [`weights()`](https://rdrr.io/r/stats/weights.html) term stops the
  function with an error; all other terms, random-effect terms included,
  are passed to `fit_function` as written.

- pv_cols:

  A character vector with the names of the plausible-value columns, at
  least two, or `NULL` (default). Required with `fit_function`: the
  columns of `data` that are fitted one at a time. With `per_pv_draws`
  it is optional and names the elements of the list; if the list has
  names, they must be the same. Without `pv_cols`, the names of
  `per_pv_draws` are used, or `PV1`, `PV2`, ... for an unnamed list.

- per_pv_draws:

  A list of draws computed elsewhere, one element per plausible value
  (at least two), or `NULL` (default). Each element is a numeric matrix
  or data frame with one row per posterior draw (at least two), finite
  values and unique column names. The fixed-effect columns (see
  `param_map`) must have the same names, in the same order, in every
  element. Give either `per_pv_draws` or `fit_function`, not both.

- control:

  A
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  object with `method = "per_pv"`. The default is
  `pv_control(method = "per_pv")`, and `NULL` gives the same. Its
  `conf_level` sets the interval level; its sampler settings and
  `backend` are only passed on to `fit_function`.

- family, prior:

  Passed unchanged to `fit_function`; pvstackr does not check them.
  Default `NULL`.

- fit_function:

  Your function that fits the model to one plausible value, or `NULL`
  (default) when you give `per_pv_draws`. It is called once per
  plausible value with the arguments `formula` (with `OUTCOME` replaced
  by that plausible-value column), `data`, `family`, `prior`, `chains`,
  `iter`, `warmup`, `cores`, `seed`, `backend`, `file` and `file_refit`
  (the settings from `control` and the cache file; see `cache_dir`),
  plus the elements of `additional_args`. There is no `weights` argument
  (see Details). It may return any object that `draws_function` accepts.

- draws_function:

  Your function that takes the object returned by `fit_function` and
  returns its draws, in the form described for `per_pv_draws`. Required
  with `fit_function`.

- param_map:

  `NULL` (default) or a named list that says which columns of the draws
  are the fixed effects, by name (`fe_names`) or by position (`fe_idx`),
  as for
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md).
  With `NULL`, the columns whose names start with `b_` are the fixed
  effects. Only the fixed effects are used, and their column names
  become the terms of the estimate table, so they must start with `b_`.
  Give a list to leave out other columns that start with `b_`, such as
  distributional parameters named `b_sigma_*`.

- weight_col:

  The name of the final weight column, or `NULL` (default). When `data`
  is used, pvstackr checks that the weights are positive and records the
  column in the fit's `design`, but it does not use the weights: they
  are not passed to `fit_function`.

- rep_weight_cols, fay_k, id_cols:

  Replicate-weight columns, the Fay coefficient and row-identifier
  columns, as in
  [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md).
  With `data` and `formula`, they are checked and recorded in the fit's
  `design`; the calculation does not use them. Defaults `NULL`, `0.5`
  and `NULL`.

- df_method:

  The rule for the degrees of freedom: `"classic"` (default) or
  `"barnard_rubin"`, which needs `df_complete`. Both rules are defined
  in
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
  With either rule the intervals are descriptive.

- df_complete:

  For `df_method = "barnard_rubin"`, the complete-data degrees of
  freedom, which you state: the degrees of freedom that the analysis
  would have if the outcome were observed directly. Give one positive
  number for all fixed effects, or a vector with one value per fixed
  effect, named by the fixed-effect columns (such as `b_Intercept`); an
  unnamed vector of several values is an error. Default `NULL`; giving a
  value with `df_method = "classic"` is an error.

- cache_dir, cache_stem:

  The folder and the start of the file name that are passed to
  `fit_function`; defaults `"cache"` and `"pvstackr-per-pv"`. For each
  plausible value, `fit_function` receives
  `file = file.path(path.expand(cache_dir), paste0(cache_stem, "-", pv))`,
  where `pv` is the plausible-value column (such as
  `"cache/pvstackr-per-pv-PV1READ"`), and `file_refit = "on_change"`;
  with `cache_dir = NULL` it receives `file = NULL` and
  `file_refit = "never"`. They are meant for the arguments of the same
  names of the brms function `brm()`. pvstackr does not create the
  folder.

- additional_args:

  A named list of further arguments passed to `fit_function` for every
  plausible value, [`list()`](https://rdrr.io/r/base/list.html) by
  default. They may not repeat the arguments that pvstackr sets (see
  `fit_function`).

## Value

A `pvstackr_fit` object with `method = "per_pv"`. Read it with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md);
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
returns `NULL` for this method, and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
describes every part of the object. The estimate table has one row per
fixed effect with the combined estimate, its standard error and its
degrees of freedom; besides the columns of every method, it has
`pooling_source` and `pooling_hash`.
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
returns the combined result, a list of class `pvstackr_reference_pool`
with `beta` (the estimates), `U_bar`, `B`, `T_MI`, `se`, `df`, `fmi` and
related fields as in
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
and `per_pv` (the mean `beta` and the covariance `U` of the draws of
each plausible value).
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
returns `reference` (where the draws came from, the plausible values,
the numbers of draws and, with `control$return_draws = TRUE`, the
fixed-effect draws of each plausible value in `per_pv_draws`) and
`pooling` (the combined quantities).

## Details

### Calculation

For each plausible value \\m\\, pvstackr takes the fixed-effect columns
of the draws (see `param_map`) and computes their mean \\\hat\beta_m\\
and their covariance matrix \\U_m\\. Rubin's rules then combine the
\\M\\ results as in
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md):
the estimate \\\bar\beta\\ is the mean of the \\\hat\beta_m\\, and the
total covariance is \\T\_{\mathrm{MI}} = \bar U + (1 + 1/M) B\\, where
\\\bar U\\ is the mean of the \\U_m\\ and \\B\\ is the covariance matrix
of the \\\hat\beta_m\\. The standard errors are the square roots of the
diagonal of \\T\_{\mathrm{MI}}\\, and the degrees of freedom follow
`df_method`.

### What is reported

As for every method, only the fixed effects are reported
([`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
gives the scope and the full rule for the intervals). Each \\U_m\\ is
the model-based posterior covariance of the draws, not a design-based
variance from replicate weights, so every row of the estimate table has
`interval_role = "reference_classic_rubin"` or
`"reference_barnard_rubin"` and `coverage_claim_allowed = FALSE`. The
replicate weights play no part in this method.

### Fitting functions and status

pvstackr itself fits no model for this method, and there is no bundled
engine: `pv_control(backend = "brms")` does not supply one. With
`fit_function`, pvstackr replaces `OUTCOME` in `formula` by each
plausible-value column in turn (for example `PV1READ ~ x + female`),
calls `fit_function` once for each, and passes each result to
`draws_function`. It does not check the model further: random-effect
terms such as `(1 | school)`, `family` and `prior` reach `fit_function`
unchanged. The survey weights are not passed as an argument to
`fit_function`; a weighted fit must take them from `data` inside that
function or from `additional_args`.

pvstackr collects no sampler diagnostics (R-hat, effective sample sizes,
divergent transitions) from these fits and has no other check for this
method, so the status of a `per_pv` fit is always `"ok"`. Check the
convergence of each fit yourself.

## References

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.*
Wiley.

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948-955.

## See also

[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
to compare fits, and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
for the parts of a fit.

Other pvstackr-fitting:
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)

## Examples

``` r
# Draws computed elsewhere: one matrix per plausible value, one row per
# draw, with the fixed-effect columns b_Intercept and b_x. Real draws
# would come from a model fitted to each plausible value; these are
# simulated to show the input.
set.seed(1)
make_draws <- function(n, b0, bx) {
  cbind(
    b_Intercept = rnorm(n, b0, 0.5),
    b_x         = rnorm(n, bx, 0.2)
  )
}
per_pv_draws <- list(
  PV1READ = make_draws(200, 1.0, 0.30),
  PV2READ = make_draws(200, 1.1, 0.28)
)
fit_ref <- pv_fit_reference(
  per_pv_draws = per_pv_draws,
  pv_cols      = c("PV1READ", "PV2READ"),
  control      = pv_control(method = "per_pv")
)
fit_ref                     # method, status "ok" and interval note
#> pvstackr fit
#>   method: per_pv
#>   status: ok
#>   fixed effects: 2
#>   target: per_pv_rubin_draws
#>   draws: not retained
#>   diagnostics: reference, pooling
#>   interval note: intervals are descriptive rather than coverage-claimable.
get_estimates(fit_ref)      # interval_role "reference_classic_rubin"
#>          term  estimate        se std.error        df df_method df_complete
#> 1 b_Intercept 1.0484282 0.5038295 0.5038295 8103.9103   classic          NA
#> 2         b_x 0.2840975 0.2143992 0.2143992  704.1008   classic          NA
#>   conf_level    conf_low conf_high    conf.low conf.high
#> 1       0.95  0.06079298 2.0360634  0.06079298 2.0360634
#> 2       0.95 -0.13684085 0.7050359 -0.13684085 0.7050359
#>             interval_role coverage_claim_allowed parameter_scope
#> 1 reference_classic_rubin                  FALSE    fixed_effect
#> 2 reference_classic_rubin                  FALSE    fixed_effect
#>        target_source target_hash     pooling_source pooling_hash
#> 1 per_pv_rubin_draws    c226f7bb per_pv_rubin_draws     c226f7bb
#> 2 per_pv_rubin_draws    c226f7bb per_pv_rubin_draws     c226f7bb
```
