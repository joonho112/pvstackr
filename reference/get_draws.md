# Get the calibrated draws of a fit

`get_draws()` returns the calibrated fixed-effect draws of a
`stack_direct` fit: the draws of the stacked fit after the Cholesky
calibration correction (CCC), which gives them the mean and covariance
of the target.

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

  A fit (class `pvstackr_fit`).

- ...:

  Ignored.

## Value

A numeric matrix with one row per draw of the stacked fit and one column
per fixed effect, named as in the estimate table (such as
`b_Intercept`). `NULL` for a fit made with `return_draws = FALSE` (such
as the bundled example fit), for a blocked fit, and for `per_pv` and
`stack_psis` fits.

`get_draws()` stops with an error if the fit was changed after it was
created, and for the inspection object that
[`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
makes from a `stack_psis` fit of an earlier pvstackr version, which has
no draws.

## Details

A fit keeps these draws only with `return_draws = TRUE`, the default of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).
`per_pv` and `stack_psis` fits keep their draws in
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
instead; "The fit object" in
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
says where.

## See also

[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
for how the draws are calibrated.

Other pvstackr-accessors:
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)

## Examples

``` r
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
if (nzchar(path)) {
  fit <- readRDS(path)$fit
  get_draws(fit)   # NULL: the example fit was saved without its draws
}
#> NULL
```
