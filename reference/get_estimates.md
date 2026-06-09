# Access pvstackr Estimates

Return the reportable fixed-effect estimate table for a fit, or the
aligned estimate table for a method comparison. The returned columns
depend on the method contract but keep interval and provenance fields
intact.

## Usage

``` r
get_estimates(x, ...)

# Default S3 method
get_estimates(x, ...)

# S3 method for class 'pvstackr_fit'
get_estimates(x, ...)

# S3 method for class 'pvstackr_method_comparison'
get_estimates(x, ...)
```

## Arguments

- x:

  A pvstackr object.

- ...:

  Reserved for future extensions.

## Value

A data frame of reportable fixed-effect estimates. Fit estimate tables
include interval/provenance columns such as `df_method`, `df_complete`,
`interval_role`, `coverage_claim_allowed`, `target_source`,
`target_hash`, `pooling_source`, and `pooling_hash` when those columns
are part of the method contract.

## See also

[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md);
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md).

Other pvstackr-accessors:
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)

## Examples

``` r
path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
if (nzchar(path)) {
  fit <- readRDS(path)$fit
  head(get_estimates(fit))
}
#>          term   estimate        se std.error       df df_method df_complete
#> 1 b_Intercept 457.894088 1.2873118 1.2873118 1.021194   classic          NA
#> 2         b_x  46.883361 0.3717929 0.3717929 1.402308   classic          NA
#> 3    b_female   2.143702 3.5550309 3.5550309 1.013730   classic          NA
#>   conf_level  conf_low conf_high  conf.low conf.high             interval_role
#> 1       0.95 442.31804 473.47013 442.31804 473.47013 descriptive_classic_rubin
#> 2       0.95  44.41457  49.35215  44.41457  49.35215 descriptive_classic_rubin
#> 3       0.95 -41.60687  45.89428 -41.60687  45.89428 descriptive_classic_rubin
#>   coverage_claim_allowed parameter_scope          target_source target_hash
#> 1                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
#> 2                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
#> 3                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
```
