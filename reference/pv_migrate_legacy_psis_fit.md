# Convert a Legacy PSIS Fit to a Safe Inspection Object

`pv_migrate_legacy_psis_fit()` returns a current, fully validated
`stack_psis` fit unchanged. A non-current or invalid historical
`stack_psis` fit is never promoted to reportable status: it is projected
to a separately allocated inspection-only object containing only bounded
Pareto-k decision evidence and an explicit redaction record. The source
object is not modified.

## Usage

``` r
pv_migrate_legacy_psis_fit(fit)

# S3 method for class 'pvstackr_legacy_psis_inspection'
print(x, ...)

# S3 method for class 'pvstackr_legacy_psis_inspection'
summary(object, ...)

# S3 method for class 'summary.pvstackr_legacy_psis_inspection'
print(x, ...)
```

## Arguments

- fit:

  A current or historical `pvstackr_fit` with `method = "stack_psis"`.

- x, object:

  A `pvstackr_legacy_psis_inspection` object, or its summary for the
  summary print method.

- ...:

  Ignored.

## Value

The unchanged current fit when full current validation succeeds;
otherwise a `pvstackr_legacy_psis_inspection` object with no reportable
estimates, draws, pooling, weights, backend, or data payload.

## Details

Inspection objects allow
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
and compact print/summary methods.
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
and
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
fail explicitly, because an historical warning, blocked, unevaluated,
incomplete, or `k >= 0.7` result cannot be grandfathered into current
reportable output. When historical diagnostic and control thresholds
differ, the effective inspection gate is the minimum of both valid
declarations and the immutable `0.7` ceiling.
