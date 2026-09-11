# Get the diagnostics of a fit or comparison

`get_diagnostics()` returns the diagnostics list of a fit or a method
comparison: the results of pvstackr's checks and the details of a method
that are not in the estimate table.

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

  A fit (class `pvstackr_fit`) or a method comparison from
  [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md).

- ...:

  Ignored.

## Value

A named list whose elements depend on the object:

- A `stack_direct` fit: `preflight` (the check that `data`, `formula`
  and `target` match), `sampler` and `sampler_gate` (the sampler
  diagnostics and pvstackr's check of them), `stack_fit` and
  `stack_fit_warnings` (records and notes from the stacked fit) and
  `ccc` (the calibration diagnostics, such as `delta_c_max` and
  `kappa_A`). A blocked fit has only `preflight`, `sampler`,
  `sampler_gate`, `redaction` (what was removed) and, when the
  calibration check blocked it, `ccc`, with its values grouped as
  `center`, `conditioning`, `residual` and `prior`.

- A `per_pv` fit: `reference` (a record of the per-plausible-value fits
  and, when kept, their draws) and `pooling` (their Rubin's-rules
  combination).

- A `stack_psis` fit: `psis` (the Pareto k-hat values, the decision and
  the weight diagnostics), `pooling` and `weighted` (the weighted result
  of each plausible value and, when kept, the stacked draws and the
  normalized weights). A blocked fit has only `psis` and `redaction`.

- A method comparison: `reference_method`, `methods`, `statuses`,
  `blocked_methods`, `warning_methods`, `agreement`,
  `method_diagnostics`, `timing` and `target_overlap` (see
  [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)).

- The inspection object that
  [`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  makes from a `stack_psis` fit of an earlier pvstackr version: `psis`
  (the Pareto k-hat values and the decision) and `redaction` (what was
  removed).

`get_diagnostics()` stops with an error if the fit or comparison was
changed after it was created.

## See also

[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
and
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
for the elements of a fit's diagnostics; "Status and checks" in
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
for the thresholds of the checks and their reason codes.

Other pvstackr-accessors:
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)

## Examples

``` r
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
if (nzchar(path)) {
  fit <- readRDS(path)$fit
  names(get_diagnostics(fit))
}
#> [1] "preflight"          "sampler"            "sampler_gate"      
#> [4] "stack_fit"          "stack_fit_warnings" "ccc"               
```
