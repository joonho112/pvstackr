# Display Methods for pvstackr Method Comparisons

Compact console [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) for a
`pvstackr_method_comparison` object. The
[`print()`](https://rdrr.io/r/base/print.html) method shows the
reference method, the compared method labels, the fixed-effect count,
any blocked or warning methods, the provenance note (the shared-target
caveat), and the interval note. The
[`summary()`](https://rdrr.io/r/base/summary.html) method builds a
structured `summary.pvstackr_method_comparison` object; its
[`print()`](https://rdrr.io/r/base/print.html) method also renders the
agreement table.

## Usage

``` r
# S3 method for class 'pvstackr_method_comparison'
print(x, ...)

# S3 method for class 'pvstackr_method_comparison'
summary(object, ...)

# S3 method for class 'summary.pvstackr_method_comparison'
print(x, ...)
```

## Arguments

- ...:

  Ignored.

- object, x:

  A `pvstackr_method_comparison` object (or its summary, for the
  `print.summary.pvstackr_method_comparison` method).

## Value

The `print` methods return their input invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`summary.pvstackr_method_comparison` list with fields:

- `reference_method`:

  The method used as the comparison reference.

- `methods`, `method_labels`, `n_methods`:

  The compared method ids, their display labels, and the number of
  compared methods.

- `n_terms`:

  Number of distinct reportable fixed-effect terms.

- `blocked_methods`, `warning_methods`:

  Methods whose reportable output was blocked, and methods that fit with
  warnings.

- `interval_note`:

  Set when reportable intervals are descriptive rather than
  coverage-claimable; `NA` otherwise.

- `provenance_note`:

  Shared-provenance caveat (close agreement is not independent
  corroboration when compared methods share a target source, pooling
  source, or estimand construction); `NA` otherwise.

- `estimate_table`, `diagnostic_table`:

  The aligned cross-method estimate table and the per-method diagnostic
  table.

- `agreement`:

  The descriptive cross-method agreement table.

- `timing`:

  Per-method timing metadata.

- `summary_schema_version`, `source_validation`,
  `source_reportability_comparison`, `validation`:

  The current summary schema, deep-valid compact source comparison and
  its stamp, and the owned-summary SHA-256 record used before printing.

## See also

[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
