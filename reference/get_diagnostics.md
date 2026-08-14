# Access pvstackr Diagnostics

Return the structured diagnostics list stored on a fit or method
comparison. Diagnostics keep method-specific details that are
intentionally not flattened into the reportable estimate table.

## Usage

``` r
get_diagnostics(x, ...)

# Default S3 method
get_diagnostics(x, ...)

# S3 method for class 'pvstackr_fit'
get_diagnostics(x, ...)

# S3 method for class 'pvstackr_method_comparison'
get_diagnostics(x, ...)

# S3 method for class 'pvstackr_legacy_psis_inspection'
get_diagnostics(x, ...)
```

## Arguments

- x:

  A pvstackr object.

- ...:

  Reserved for future extensions.

## Value

A structured diagnostics list. For a legacy PSIS inspection object this
is bounded Pareto-k evidence plus its redaction record; historical
pooling, weights, estimates, and draws are never returned.

## Details

Diagnostic keys are method-specific; inspect
`names(get_diagnostics(x))`. A current `stack_direct` fit carries
top-level `sampler` and `sampler_gate` records plus `preflight`,
`stack_fit`, `stack_fit_warnings`, and `ccc`; sampler-blocked fits
retain only slim preflight/sampler/gate evidence and the independently
valid external target, rebuilt from an exact recursive allowlist with no
formula object in preflight and a safe formula environment on the target
snapshot. A legacy cached stack-direct fit may predate the sampler keys.
A `per_pv` fit carries `reference` and `pooling`; a `stack_psis` fit
carries `psis` (status and Pareto-k), `pooling`, and `weighted`. A
method comparison instead carries comparison-level keys such as
`agreement`, `method_diagnostics`, `timing`, and `target_overlap` (the
shared-provenance summary).

## See also

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md);
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md).

Other pvstackr-accessors:
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)

## Examples

``` r
path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
if (nzchar(path)) {
  fit <- readRDS(path)$fit
  names(get_diagnostics(fit))
}
#> [1] "preflight"          "sampler"            "sampler_gate"      
#> [4] "stack_fit"          "stack_fit_warnings" "ccc"               
```
