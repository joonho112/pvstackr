# Get the target of a fit

`get_target()` returns the target that the estimates of a fit come from.
For a `stack_direct` fit this is the BRR-Fay target from
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
to which the stacked fit was calibrated; for a `per_pv` fit it is the
Rubin's-rules combination of the draws of the per-plausible-value fits.
A `stack_psis` fit has no target object.

## Usage

``` r
get_target(x, ...)

# Default S3 method
get_target(x, ...)

# S3 method for class 'pvstackr_fit'
get_target(x, ...)

# S3 method for class 'pvstackr_brr_target'
get_target(x, ...)

# S3 method for class 'pvstackr_legacy_psis_inspection'
get_target(x, ...)
```

## Arguments

- x:

  A fit (class `pvstackr_fit`), or a target from
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).

- ...:

  Ignored.

## Value

For a `stack_direct` fit, its `pvstackr_brr_target` object, also when
the fit is blocked.
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
lists its elements, such as `beta`, `T_MI`, `df` and `fmi` (the fraction
of missing information, which the estimate table does not have). For a
`per_pv` fit, a `pvstackr_reference_pool` object with the same kind of
elements. For a `stack_psis` fit, `NULL`: its estimate table has
`target_source = "stack_psis_rubin_pooling"`, but that is only a label
for the combined result, not a target object. Given a target,
`get_target()` checks it and returns it unchanged.

`get_target()` stops with an error if the fit or target was changed
after it was created. For the inspection object that
[`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
makes from a `stack_psis` fit of an earlier pvstackr version, it returns
`NULL`.

## See also

[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
for the elements of a target.

Other pvstackr-accessors:
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)

## Examples

``` r
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
if (nzchar(path)) {
  fit <- readRDS(path)$fit
  get_target(fit)
}
#> pvstackr BRR-Fay target
#>   fixed effects: 3
#>   plausible values: 2
#>   replicate weights: 4
#>   fay_k: 0.5
#>   df method: classic
#>   interval role: descriptive_classic_rubin
#>   source: external_brr_fay_rubin
```
