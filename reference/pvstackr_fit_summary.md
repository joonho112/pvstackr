# Display Methods for pvstackr Fits

Compact console [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) for a
[pvstackr_fit](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
object. The [`print()`](https://rdrr.io/r/base/print.html) method shows
the method id, fit status, fixed-effect count, target source, draw
retention, diagnostics keys, and the interval note (set when reportable
intervals are descriptive rather than coverage-claimable). The
[`summary()`](https://rdrr.io/r/base/summary.html) method builds a
structured `summary.pvstackr_fit` object; its
[`print()`](https://rdrr.io/r/base/print.html) method adds the compact
estimate columns (`term`, `estimate`, `se`, `df`, `conf_low`,
`conf_high`, where present).

## Usage

``` r
# S3 method for class 'pvstackr_fit'
print(x, ...)

# S3 method for class 'pvstackr_fit'
summary(object, ...)

# S3 method for class 'summary.pvstackr_fit'
print(x, ...)
```

## Arguments

- ...:

  Ignored.

- object, x:

  A `pvstackr_fit` object (or its summary, for the
  `print.summary.pvstackr_fit` method).

## Value

The `print` methods return their input invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`summary.pvstackr_fit` list with fields:

- `method`, `status`:

  Method id and fit status.

- `n_terms`, `terms`:

  Number and names of reportable fixed-effect terms.

- `has_target`, `target_source`:

  Whether a formal target object is carried, and the estimate-row
  target-source provenance label.

- `has_draws`, `draw_dim`:

  Whether reportable draws were retained, and their `c(nrow, ncol)`
  dimension (`c(0L, 0L)` when not retained).

- `diagnostic_keys`:

  Names of the entries in the diagnostics list.

- `interval_note`:

  Set when reportable intervals are descriptive rather than
  coverage-claimable; `NA` otherwise.

- `reason_codes`, `warnings`:

  Status reason codes and any captured fit warnings.

- `estimates`, `diagnostics`:

  The reportable fixed-effect estimate table and the structured
  diagnostics list.

- `summary_schema_version`, `source_validation`,
  `source_reportability_fit`, `validation`:

  The current summary schema, deep-valid compact source fit and its
  stamp, and the owned-summary SHA-256 record used before printing.

## See also

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
