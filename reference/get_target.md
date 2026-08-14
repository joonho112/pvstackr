# Access pvstackr Targets

Return the formal target object carried by a fit. Estimate-row
provenance labels such as `target_source` are separate from this
accessor and do not imply that every method has a formal target object.

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

  A pvstackr fit or target object.

- ...:

  Reserved for future extensions.

## Value

The target object used by a fit, or `NULL` when the method has no target
component. Estimate-row `target_source` labels are provenance metadata
and may be present even when a method, such as `stack_psis`, does not
carry a formal target object. A legacy PSIS inspection object also
returns `NULL`.

## See also

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md);
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md).

Other pvstackr-accessors:
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)

## Examples

``` r
path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
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
