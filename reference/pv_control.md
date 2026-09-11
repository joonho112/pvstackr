# Settings for the fitting functions

`pv_control()` collects the settings used by
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
and by
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)
and
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md):
the method, the sampler settings and engine passed to the model-fitting
function, the interval level, the Pareto k-hat cut-off of
`"stack_psis"`, and what a fitted object keeps. It checks each value and
returns them in one object.

## Usage

``` r
pv_control(
  method = "stack_direct",
  chains = 4L,
  iter = 2000L,
  warmup = NULL,
  cores = 1L,
  seed = NULL,
  backend = "none",
  conf_level = 0.95,
  psis_k_threshold = 0.7,
  center = "target",
  allow_target_nearpd = FALSE,
  return_draws = TRUE,
  keep_data = FALSE,
  keep_backend_fit = FALSE,
  keep_log_lik = FALSE,
  verbose = FALSE
)

# S3 method for class 'pvstackr_control'
print(x, ...)
```

## Arguments

- method:

  The method these settings are for: `"stack_direct"` (default),
  `"per_pv"` or `"stack_psis"`. It must match the method you fit with:
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
  stops with an error if it differs from its `method` argument, and each
  method function accepts only its own method.

- chains:

  Number of Markov chains, a whole number of at least 1. Default `4L`.

- iter:

  Number of iterations per chain, warmup included, a whole number of at
  least 2. Default `2000L`.

- warmup:

  Number of warmup iterations per chain, a whole number from 0 to
  `iter - 1`. `NULL` (default) uses `floor(iter / 2)`, which is 1000 for
  the default `iter`.

- cores:

  Number of cores the model-fitting function may use, a whole number of
  at least 1. Default `1L`.

- seed:

  The random seed, a whole number of at least 0, or `NULL` (default) to
  leave the seed to the model-fitting function.

- backend:

  The engine for `"stack_direct"`: `"none"` (default), `"injected"`,
  `"brms"` or `"cmdstanr"`. Only `"brms"` selects one: when you give no
  `fit_function`, `"stack_direct"` then uses the bundled brms engine
  (see Details). With any other value, `"stack_direct"` needs your own
  `fit_function`, `draws_function` and `diagnose_function` (see
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)).
  `"per_pv"` and `"stack_psis"` have no bundled engine and always need
  your own functions or draws computed elsewhere, whatever the value.

- conf_level:

  The level of the intervals in the estimate table, a number strictly
  between 0 and 1. Default `0.95`. It applies to all three methods; the
  `conf_level` stored in a
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  object does not change the reported intervals.

- psis_k_threshold:

  The cut-off for the Pareto k-hat values (Pareto shape estimates) that
  come with the importance weights of `"stack_psis"`, a number greater
  than 0 and at most 0.7. Default `0.7`. pvstackr blocks a
  `"stack_psis"` fit unless the k-hat of every plausible value is finite
  and below this cut-off, so a smaller value makes the check stricter;
  values above 0.7 are not allowed.

- center:

  Where the calibrated fixed-effect draws are centered: `"target"`
  (default) or `"posterior"`. With `"target"`, the Cholesky calibration
  correction (CCC) gives the draws the mean and covariance of the
  target, whose estimates are the ones reported. `"stack_direct"`
  requires this value. `"posterior"` would keep the mean of the stacked
  draws and calibrate only their covariance; it is accepted here, but no
  fitting function uses it:
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  stops with an error for it. `"per_pv"` and `"stack_psis"` ignore
  `center`.

- allow_target_nearpd:

  Whether a target covariance matrix that is not positive definite may
  be repaired. pvstackr has no such repair and stops with an error on
  such a target, so the value must be `FALSE` (default); `TRUE` stops
  with an error.

- return_draws:

  Whether a fit keeps its fixed-effect draws (see Details). Default
  `TRUE`.

- keep_data:

  Whether a fit keeps the data frame. Default `FALSE`: the fit then
  keeps only a description of the data, such as column names and
  checksums.

- keep_backend_fit:

  Whether a fit keeps the object returned by the model-fitting function,
  such as a brms fit (for `"per_pv"`, one per plausible value). This
  object may contain the data, so the fitting functions stop with an
  error when `keep_backend_fit = TRUE` and `keep_data = FALSE`. Default
  `FALSE`.

- keep_log_lik:

  Whether a fit keeps the log-likelihood draws. Only `"stack_direct"`
  can extract them:
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  does so when given `extract_log_lik = TRUE` and a `log_lik_function`.
  Set `TRUE` only then: otherwise a `"stack_direct"` fit, or a
  `"stack_psis"` fit from a `fit_function`, stops with an error. Default
  `FALSE`.

- verbose:

  Stored in the object but not used by any pvstackr function. Default
  `FALSE`. For progress messages while the target is computed, use the
  `verbose` argument of
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).

- x:

  A `pvstackr_control` object from `pv_control()`.

- ...:

  Ignored.

## Value

A `pvstackr_control` object: a list of the 16 settings in the order of
the arguments, with `warmup` filled in and whole numbers stored as
integers. Pass it as `control` to
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
or to a method function. To change a setting, call `pv_control()` again
rather than editing the list: the fitting functions stop with an error
if a value was replaced by one of another type, as `ctrl$chains <- 2`
does (it stores a double).
[`print()`](https://rdrr.io/r/base/print.html) shows `method`,
`backend`, `iter`, `warmup` and `chains`, notes that target repair is
not supported, and returns the object invisibly.

## Details

### Which settings each method uses

`chains`, `iter`, `warmup`, `cores`, `seed` and `backend` are passed to
the function that fits the model: the bundled brms engine or your own
`fit_function`. They have no effect when a method starts from draws
computed elsewhere. For `"stack_direct"`, pvstackr also checks that the
sampler diagnostics report `chains` chains with `iter - warmup` draws
each after warmup, and blocks the fit (status `"blocked"`) if they do
not.

`center` is used only by `"stack_direct"` and `psis_k_threshold` only by
`"stack_psis"`. `conf_level` and the settings that decide what a fitted
object keeps apply to all three methods.

### The bundled brms engine

With `backend = "brms"` and no `fit_function`, `"stack_direct"` fits the
stacked model with the brms function `brm()`, through
[`pv_backend_brms_fit_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md).
This needs the brms and posterior packages. If cmdstanr is installed,
CmdStan must be configured, otherwise the fit stops with an error; if
cmdstanr is not installed, brms uses rstan. The other `backend` values
select no engine: they are passed to your `fit_function` as its
`backend` argument and recorded in the fit.

### What a fitted object keeps

With `return_draws = TRUE` (the default), a fit keeps its fixed-effect
draws: for `"stack_direct"` the calibrated draws, read with
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md);
for `"per_pv"` the draws of each plausible value; for `"stack_psis"` the
stacked draws with the normalized weights for each plausible value. The
last two are in
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md),
and
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
returns `NULL` for them. Draws of other parameters, such as `sigma`, are
not kept.

`keep_data`, `keep_backend_fit` and `keep_log_lik` (all `FALSE` by
default) keep the data, the object returned by the model-fitting
function and the log-likelihood draws; each makes a saved fit larger. A
blocked fit keeps none of these parts: its `control` records all four
settings as `FALSE`, whatever was requested.
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
describes the parts of a fitted object.

## See also

[`pv_backend_brms_fit_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
for the bundled brms engine;
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
for the parts of a fitted object.

Other pvstackr-fitting:
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)

## Examples

``` r
# Default settings, for method = "stack_direct".
ctrl <- pv_control()
ctrl
#> pvstackr control
#>   method: stack_direct
#>   backend: none
#>   iter/warmup/chains: 2000/1000/4
#>   target repair: unsupported (disabled)

# A control's method must match the method you fit with.
ctrl_psis <- pv_control(method = "stack_psis", psis_k_threshold = 0.7)
ctrl_psis$method
#> [1] "stack_psis"

# Settings that make a fit keep more. keep_backend_fit = TRUE needs
# keep_data = TRUE, and keep_log_lik = TRUE needs extract_log_lik = TRUE
# and a log_lik_function in pv_fit_direct().
ctrl_heavy <- pv_control(
  keep_data = TRUE, keep_backend_fit = TRUE, keep_log_lik = TRUE
)
c(ctrl_heavy$keep_data, ctrl_heavy$keep_backend_fit, ctrl_heavy$keep_log_lik)
#> [1] TRUE TRUE TRUE
```
