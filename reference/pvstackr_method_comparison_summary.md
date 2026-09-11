# Print and summarize a method comparison

[`print()`](https://rdrr.io/r/base/print.html) shows a short overview of
a comparison made by
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md).
[`summary()`](https://rdrr.io/r/base/summary.html) returns the same
overview, together with the tables of the comparison, as a list of class
`summary.pvstackr_method_comparison`; its
[`print()`](https://rdrr.io/r/base/print.html) method also shows the
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

  A comparison (class `pvstackr_method_comparison`) from
  [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md);
  for the [`print()`](https://rdrr.io/r/base/print.html) method of a
  summary, the list returned by
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns a
list of class `summary.pvstackr_method_comparison` with these elements:

- `reference_method`: the label of the reference fit.

- `methods`, `method_labels` and `n_methods`: the method of each fit
  (named by label), the labels and the number of fits.

- `n_terms`: the number of fixed effects in the comparison.

- `blocked_methods` and `warning_methods`: the labels of the fits with
  status `"blocked"` or `"warning"`.

- `interval_note` and `provenance_note`: the texts of the two note
  lines; `interval_note` is `NA` when all intervals are
  coverage-claimable.

- `estimate_table`, `diagnostic_table`, `agreement` and `timing`: the
  tables described in
  [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md).

- `schema_version`, `summary_schema_version`, `source_validation`,
  `source_reportability_comparison` and `validation`: the format
  versions of the comparison and of the summary; a copy of the
  comparison without its fits and the checksum of that copy; and a
  SHA-256 checksum of the summary, which
  [`print()`](https://rdrr.io/r/base/print.html) recomputes to detect
  changes.

## Details

The overview shows the label of the reference fit, the label and the
method of each fit (only their number in the summary), the number of
fixed effects, and the labels of blocked fits and of fits with status
`"warning"` when there are any. A line beginning "provenance note:"
follows for every comparison: it says that agreement is descriptive and
that a shared target, combined result or source label is not independent
confirmation
([`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
explains this). A line beginning "interval note:" is added when some or
all intervals are descriptive.

The three methods stop with an error if the comparison or the summary
was changed after it was created, and for a comparison saved by pvstackr
0.1.x, which has to be made again from current fits. A summary that
contains a `stack_psis` fit and was saved by an earlier version of
pvstackr also stops with an error when printed (see
[`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)).

## See also

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
to read a comparison.
