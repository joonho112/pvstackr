# Access pvstackr Draws

Return retained top-level reportable draws from a fit. Methods that do
not synthesize a single reportable draw matrix, or fits created with
`return_draws = FALSE`, return `NULL`.

## Usage

``` r
get_draws(x, ...)

# Default S3 method
get_draws(x, ...)

# S3 method for class 'pvstackr_fit'
get_draws(x, ...)

# S3 method for class 'pvstackr_legacy_psis_inspection'
get_draws(x, ...)
```

## Arguments

- x:

  A pvstackr fit object.

- ...:

  Reserved for future extensions.

## Value

The retained reportable draw matrix, or `NULL` when draws were not
retained or the method does not synthesize top-level reportable draws.
Per-PV reference draws and the PSIS proposal/weight pair, when retained,
remain available in diagnostics rather than through this top-level
reportable-draw accessor. An inspection-only legacy PSIS object fails
explicitly instead of returning historical draws.

## Details

This accessor returns only the synthesized top-level reportable draw
matrix. Per-PV reference draws (`per_pv`) and the PSIS fixed-effect
proposal/weight pair (`stack_psis`), when retained, are not surfaced
here; they live in the fit's diagnostics (see
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)),
not in this top-level reportable-draw accessor.

## See also

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md);
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md).

Other pvstackr-accessors:
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)

## Examples

``` r
path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
if (nzchar(path)) {
  fit <- readRDS(path)$fit
  get_draws(fit)
}
#> NULL
```
