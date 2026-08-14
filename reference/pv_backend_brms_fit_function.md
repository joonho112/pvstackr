# Bundled brms adapter for the stacked fit

The three functions that `pv_control(backend = "brms")` installs when no
adapter is injected: `pv_backend_brms_fit_function()` fits the prepared
stacked model, `pv_backend_brms_draws_function()` extracts its draws,
and `pv_backend_brms_sampler_diagnostics()` reports its sampler
diagnostics.

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

  Prepared stacked model formula, with the stacked outcome and weight
  columns already bound.

- data:

  Prepared stacked data frame.

- family:

  Response family, or `NULL` for
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html).

- prior:

  Prior specification passed through to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html), or
  `NULL`.

- chains, iter, warmup, cores, seed:

  Sampler settings from
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).

- backend:

  Resolved Stan backend, `"cmdstanr"` or `"rstan"`.

- file, file_refit:

  Cache location and refit policy for
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- ...:

  Further arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) by
  `pv_backend_brms_fit_function()`; accepted and ignored by
  `pv_backend_brms_draws_function()`.

- fit:

  The object returned by the fit function.

- draws_array_function, summary_function, nuts_function:

  Optional replacements for
  [`posterior::as_draws_array()`](https://mc-stan.org/posterior/reference/draws_array.html),
  the
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
  call that produces R-hat and ESS, and
  [`brms::nuts_params()`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html).
  Supplying these lets an engine report diagnostics without depending on
  those packages.

## Value

`pv_backend_brms_fit_function()` returns a `brmsfit`.
`pv_backend_brms_draws_function()` returns a base numeric matrix of the
`b_*` and `sigma` draws. `pv_backend_brms_sampler_diagnostics()` returns
a normalized sampler diagnostics record: maximum R-hat, minimum bulk and
tail ESS, divergence count, chains, post-warmup draws per chain, a
completeness flag, and reason codes.

## Details

They are exported so that an injected adapter and the bundled backend
run the same fitting, extraction, and diagnostic code rather than two
implementations that have to be kept in step. Given the same data and
seed, the two routes produce the same fit, the same draws matrix, and
the same sampler numbers.

Injecting the three is not the same as selecting the bundled backend,
and the recorded provenance says so. An injected fit records
`adapter_source = "injected"`, `resolved_backend = "injected"`,
`engine_id = "injected_fit_function"`, an empty `package_versions`, and
a `cmdstan_state` that was never evaluated; its sampler evidence is
labelled `injected_diagnose_function` rather than
`bundled_brms_posterior_and_nuts`. Two things the bundled route does for
you are the adapter's own responsibility when injected: checking that
brms, posterior, and a usable Stan backend are present before sampling
starts, and creating the cache directory. `cache_dir` is passed through
as given, so create it first or pass `cache_dir = NULL`.

A reportable injected fit needs all three. An adapter that supplies only
a fit and a draws function leaves the sampler evidence incomplete, and
the fit is blocked rather than reported.

Under `pv_control(backend = "brms")` a supplied `draws_function`
replaces the bundled one, but a supplied `diagnose_function` does not:
the bundled sampler record always wins, and the fit records that the
override was ignored. Supplying `fit_function` switches to the injected
route as a whole.

`pv_backend_brms_draws_function()` returns the fixed-effect (`b_*`) and
residual-scale (`sigma`) draws as a plain base matrix, the draw format
the calibration layer expects. A
[`posterior::draws_matrix`](https://mc-stan.org/posterior/reference/draws_matrix.html)
is converted rather than passed through, because the calibration
validator compares draw blocks by value and a classed matrix would fail
that comparison on class alone.

`pv_backend_brms_sampler_diagnostics()` reports maximum R-hat, minimum
bulk and tail ESS, divergence count, chains, and post-warmup draws per
chain. When a required namespace is missing or the extraction fails, it
returns an explicit incomplete record instead of silently omitting the
evidence, and the fit is blocked downstream.

## Backend resolution

cmdstanr is used only when both its namespace and a configured CmdStan
are available, and rstan is selected only when cmdstanr is absent. A fit
failure is not retried against the other backend. The bundled route
resolves this before calling the fit function; reached through an
injected adapter, the fit function applies the same policy itself, so
`pv_control(backend = )` may name any accepted value.

## See also

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
and
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
for how an adapter is selected.
