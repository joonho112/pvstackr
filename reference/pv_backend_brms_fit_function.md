# Bundled brms engine for the stacked fit

These three functions are the brms engine that
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
uses when `pv_control(backend = "brms")` is set and no `fit_function` is
given: `pv_backend_brms_fit_function()` fits the stacked model with
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html),
`pv_backend_brms_draws_function()` returns the draws of that fit, and
`pv_backend_brms_sampler_diagnostics()` returns the sampler diagnostics
that pvstackr checks. They are exported so that you can pass them to
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
yourself and replace one of them with your own function.

## Usage

``` r
pv_backend_brms_fit_function(
  formula,
  data,
  family,
  prior,
  chains,
  iter,
  warmup,
  cores,
  seed,
  backend,
  file,
  file_refit,
  ...
)

pv_backend_brms_draws_function(fit, ...)

pv_backend_brms_sampler_diagnostics(
  fit,
  draws_array_function = NULL,
  summary_function = NULL,
  nuts_function = NULL
)
```

## Arguments

- formula:

  The stacked formula built by
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
  such as
  `.pvstackr_y | weights(.pvstackr_weight) ~ 0 + pvstackrMM001 + ...`.

- data:

  The stacked data built by
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md):
  one copy of the data per plausible value, with the columns used in
  `formula`.

- family:

  A family object, or `NULL` for
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html).
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  always passes
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html).

- prior:

  `NULL` or a prior table, passed to `brm()` unchanged.
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  describes the priors it accepts.

- chains, iter, warmup, cores, seed:

  Sampler settings passed to `brm()`;
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  takes them from
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).

- backend:

  The Stan backend for `brm()`. `"cmdstanr"` and `"rstan"` are used as
  given; any other value is replaced by the backend that the bundled
  engine would choose (see "Stan backend").

- file, file_refit:

  The cache file (a path without extension) and the refit rule, passed
  to `brm()`.
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  sets them from `cache_dir` and `cache_stem`:
  `file_refit = "on_change"`, or `file = NULL` and
  `file_refit = "never"` when caching is off.

- ...:

  Further arguments. `pv_backend_brms_fit_function()` passes them to
  `brm()`
  ([`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  passes the elements of `additional_args` here);
  `pv_backend_brms_draws_function()` ignores them.

- fit:

  The object returned by the fitting function: a `brmsfit` for the
  functions documented here. With all three replacement functions below,
  `pv_backend_brms_sampler_diagnostics()` accepts any object that they
  accept.

- draws_array_function, summary_function, nuts_function:

  Functions that replace, in this order,
  [`posterior::as_draws_array()`](https://mc-stan.org/posterior/reference/draws_array.html),
  the
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
  call and
  [`brms::nuts_params()`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html),
  or `NULL` (default) to use those. With all three, the diagnostics can
  be computed without the posterior and brms packages, for example for a
  fit from another program. `draws_array_function(fit)` must return an
  array with the dimensions draws, chains and variables;
  `summary_function(draws)` a data frame with the columns `rhat`,
  `ess_bulk` and `ess_tail` and one row per variable; and
  `nuts_function(fit)` a data frame with the columns `Parameter` and
  `Value`, whose rows with `Parameter == "divergent__"` hold 1 for a
  divergent transition and 0 otherwise.

## Value

`pv_backend_brms_fit_function()` returns the `brmsfit` object from
`brm()`.

`pv_backend_brms_draws_function()` returns a base R numeric matrix with
one row per draw after warmup, over all chains, and one column for each
variable whose name starts with `b_` and for `sigma`. For the stacked
model the `b_` columns are the fixed effects of the stacked formula
(`b_pvstackrMM001`, `b_pvstackrMM002`, ...), which pvstackr renames to
the names used in the target, such as `b_Intercept`.

`pv_backend_brms_sampler_diagnostics()` returns the list that a fit made
with the bundled engine keeps as `get_diagnostics(fit)$sampler`:
`rhat_max`, `ess_bulk_min`, `ess_tail_min`, `ess_bulk_per_chain_min` and
`ess_tail_per_chain_min` (the two ESS values divided by `chains`),
`divergences`, `chains`, `post_warmup_draws_per_chain`,
`diagnostic_source` (`"bundled_brms_posterior_and_nuts"`),
`diagnostic_complete` (`TRUE` when every value was extracted) and
`diagnostic_reason_codes` (the reasons when it is `FALSE`).

## Details

### Replacing one of the functions

To change one step (fitting the stacked model, extracting its draws or
computing its sampler diagnostics), pass all three functions to
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
as `fit_function`, `draws_function` and `diagnose_function`: your own
function for that step and the functions documented here for the others.
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
describes what each function receives and returns, and what happens when
one is missing. The draws and diagnostics functions here expect a
`brmsfit`; for a fit from another program,
`pv_backend_brms_sampler_diagnostics()` can still compute the
diagnostics with its three replacement functions. With
`pv_control(backend = "brms")` and no `fit_function`, only the draws
function can be replaced: a `diagnose_function` that you pass is not
used, and the fit records that it was ignored.

When you pass the functions yourself, pvstackr checks no package and
creates no cache folder before it calls `fit_function`. A missing
posterior package is found only after sampling, when the draws are
extracted; create `cache_dir` first or set `cache_dir = NULL`. A fit
made this way has `get_diagnostics(fit)$sampler$diagnostic_source` equal
to `"injected_diagnose_function"`, also when your `diagnose_function` is
`pv_backend_brms_sampler_diagnostics()`.

### Stan backend

`pv_backend_brms_fit_function()` passes `backend` to `brm()`. For the
bundled engine, pvstackr chooses it before sampling: `"cmdstanr"` when
the cmdstanr package is installed and CmdStan is configured, and
`"rstan"` when cmdstanr is not installed. If cmdstanr is installed but
CmdStan is not configured, the fit stops with an error instead of using
rstan. When you pass the function yourself, it receives the `backend`
value of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md):
`"cmdstanr"` is passed to `brm()` unchanged, and any other value is
replaced by the choice described above. A fit that fails is not repeated
with the other backend.

### Sampler diagnostics

`pv_backend_brms_sampler_diagnostics()` computes R-hat and the bulk and
tail effective sample sizes (ESS) of every variable of the fit with
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html),
and reports the largest R-hat and the smallest ESS values, which are
totals over all chains. The number of divergent transitions comes from
[`brms::nuts_params()`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html),
and `chains` and `post_warmup_draws_per_chain` from the dimensions of
the draws. When the posterior or brms package is missing, or a value
cannot be extracted, it does not stop with an error: it returns the list
with missing values, `diagnostic_complete = FALSE` and reason codes such
as `posterior_namespace_unavailable`, and
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
then blocks the fit. The limits that pvstackr applies to these values
are listed under "Status and checks" in
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).

## See also

[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
for the arguments `fit_function`, `draws_function` and
`diagnose_function`, and
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
for `backend`.

## Examples

``` r
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)
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

# Pass the three functions yourself and replace any one of them with your
# own. A live fit samples with Stan, so it is not run here.
if (FALSE) { # \dontrun{
fit <- pv_fit_direct(
  data = pisa_tiny, formula = OUTCOME ~ x + female, target = target,
  fit_function = pv_backend_brms_fit_function,
  draws_function = pv_backend_brms_draws_function,
  diagnose_function = pv_backend_brms_sampler_diagnostics,
  cache_dir = NULL
)
} # }

# The sampler diagnostics record. Stand-in functions replace the posterior
# and brms functions, so this runs without those packages.
diagnostics <- pv_backend_brms_sampler_diagnostics(
  fit = NULL,
  draws_array_function = function(fit) array(0, dim = c(1000, 4, 3)),
  summary_function = function(draws) {
    data.frame(
      rhat = c(1.001, 1.004, 1.002),
      ess_bulk = c(3900, 4100, 3700),
      ess_tail = c(2900, 3100, 3000)
    )
  },
  nuts_function = function(fit) {
    data.frame(Parameter = "divergent__", Value = rep(0, 4000))
  }
)
str(diagnostics, digits.d = 4)
#> List of 11
#>  $ rhat_max                   : num 1.004
#>  $ ess_bulk_min               : num 3700
#>  $ ess_tail_min               : num 2900
#>  $ ess_bulk_per_chain_min     : num 925
#>  $ ess_tail_per_chain_min     : num 725
#>  $ divergences                : int 0
#>  $ chains                     : int 4
#>  $ post_warmup_draws_per_chain: int 1000
#>  $ diagnostic_source          : chr "bundled_brms_posterior_and_nuts"
#>  $ diagnostic_complete        : logi TRUE
#>  $ diagnostic_reason_codes    : chr(0) 
```
