# Inspect a stack_psis fit saved by an earlier version of pvstackr

`pv_migrate_legacy_psis_fit()` returns a `stack_psis` fit unchanged when
it passes the checks of the current version of pvstackr. Any other
`stack_psis` fit, such as one saved by an earlier version (for example
with status `"warning"`) or one that was changed after it was created,
becomes an object for inspection only (class
`pvstackr_legacy_psis_inspection`). It keeps the Pareto k-hat values and
their check, but no estimates, weights or draws. The fit that you pass
is not modified.

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

  A fit with `method = "stack_psis"`, made by the current or an earlier
  version of pvstackr, or an inspection object made by this function,
  which is checked and returned unchanged.

- x, object:

  An inspection object; for the
  [`print()`](https://rdrr.io/r/base/print.html) method of a summary,
  the list returned by
  [`summary()`](https://rdrr.io/r/base/summary.html).

- ...:

  Ignored.

## Value

`pv_migrate_legacy_psis_fit()` returns the fit itself when it passes the
current checks, also when it is blocked. Otherwise it returns a list of
class `pvstackr_legacy_psis_inspection` with the elements `method`
(`"stack_psis"`), `inspection_only` (`TRUE`), `reportable` (`FALSE`),
`reason_codes` (`"legacy_psis_inspection_only"`), `diagnostics` (`psis`
and `redaction`; see Details), `source` (the class, format version and
status of the old fit, `"unknown"` where they cannot be read),
`schema_version` (the format version of the inspection object),
`provenance` (a record of how it was made) and `warnings`.

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly; it shows the status of the old fit and `evidence_status`.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a list of
class `summary.pvstackr_legacy_psis_inspection` with `method`, `status`
(`"inspection_only"`), `source_status` (the status of the old fit),
`evidence_status`, `bad_pv_cols`, `reportable`, `diagnostic_keys`,
`schema_version` and `source_inspection` (the inspection object).

## Details

Earlier versions of pvstackr gave estimates for some `stack_psis` fits
that the current version blocks, so the inspection object does not carry
any numbers over:
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
and
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
stop with an error on it, and
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
returns `NULL`. To get estimates, fit the model again with
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md).
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md),
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) work on the
inspection object.

`get_diagnostics(x)$psis` holds the k-hat values found in the old fit
(`pareto_k`), the largest of them (`pareto_k_max`), the plausible values
whose k-hat is not finite or not below the cut-off (`bad_pv_cols`), and
the result of the check (`evidence_status`):

- `"not_evaluated"`: the old fit does not hold a finite k-hat value for
  each of at least two plausible values (`complete` is `FALSE`).

- `"failed"`: a k-hat value is not below the cut-off.

- `"legacy_unsafe"`: all k-hat values are below the cut-off, but the fit
  does not pass the current checks for another reason.

The cut-off, `effective_threshold`, is the smallest of 0.7
(`hard_threshold`) and the thresholds, where present, recorded in the
old fit's diagnostics (`diagnostic_threshold`) and settings
(`control_threshold`), so an old fit is never judged against a threshold
above 0.7. `fallback_requested` is the old `fallback` setting
(`"block"`, `"warn"` or `"unknown"`). `get_diagnostics(x)$redaction`
lists what was removed.

`pv_migrate_legacy_psis_fit()` stops with an error for a fit of another
method.
