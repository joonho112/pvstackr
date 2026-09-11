# Get the estimate table of a fit

`get_estimates()` returns the table of fixed-effect estimates of a fit,
with their standard errors, degrees of freedom and intervals. For a
method comparison it returns the table that lines up the estimates of
the compared fits.

## Usage

``` r
get_estimates(x, ...)

# Default S3 method
get_estimates(x, ...)

# S3 method for class 'pvstackr_fit'
get_estimates(x, ...)

# S3 method for class 'pvstackr_method_comparison'
get_estimates(x, ...)

# S3 method for class 'pvstackr_legacy_psis_inspection'
get_estimates(x, ...)
```

## Arguments

- x:

  A fit (class `pvstackr_fit`) from
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
  or a method function, or a method comparison from
  [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md).

- ...:

  Ignored.

## Value

A data frame with one row per fixed effect. For a `stack_direct` fit it
has 17 columns: `term`, `estimate`, `se`, `std.error`, `df`,
`df_method`, `df_complete`, `conf_level`, `conf_low`, `conf_high`,
`conf.low`, `conf.high`, `interval_role`, `coverage_claim_allowed`,
`parameter_scope`, `target_source` and `target_hash`. `per_pv` and
`stack_psis` fits add further columns. The fraction of missing
information is not a column; for `stack_direct` and `per_pv` fits it is
`get_target(fit)$fmi`. A blocked fit gives an empty data frame. For a
method comparison, the data frame has one row per method and fixed
effect (see
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)).

`get_estimates()` stops with an error if the fit or comparison was
changed after it was created, and for the inspection object that
[`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
makes from a `stack_psis` fit of an earlier pvstackr version, which has
no estimates.

## Details

Only the fixed effects are reported. `coverage_claim_allowed` says
whether pvstackr's reporting rule lets you read the interval of a row as
a confidence interval with nominal coverage;
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
gives the rule, and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
describes every column.

## See also

[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
for the columns of the estimate table.

Other pvstackr-accessors:
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)

## Examples

``` r
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
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
#>   coverage_claim_allowed parameter_scope          target_source
#> 1                  FALSE    fixed_effect external_brr_fay_rubin
#> 2                  FALSE    fixed_effect external_brr_fay_rubin
#> 3                  FALSE    fixed_effect external_brr_fay_rubin
#>                                                               target_hash
#> 1 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
#> 2 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
#> 3 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
```
