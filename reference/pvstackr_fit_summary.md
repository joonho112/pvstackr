# Print and summarize a fit

[`print()`](https://rdrr.io/r/base/print.html) shows a short overview of
a fit (class `pvstackr_fit`, see
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)).
[`summary()`](https://rdrr.io/r/base/summary.html) returns the same
overview, together with the estimate table and the diagnostics, as a
list of class `summary.pvstackr_fit`; its
[`print()`](https://rdrr.io/r/base/print.html) method also shows the
estimates.

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

  A fit (class `pvstackr_fit`); for the
  [`print()`](https://rdrr.io/r/base/print.html) method of a summary,
  the list returned by
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns a
list of class `summary.pvstackr_fit` with these elements:

- `method`, `status`, `reason_codes` and `warnings`: as in the fit.

- `n_terms` and `terms`: the number and the names of the fixed effects
  in the estimate table (`0` and
  [`character()`](https://rdrr.io/r/base/character.html) for a blocked
  fit).

- `has_target` and `target_source`: whether the fit has a target object
  ([`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
  does not return `NULL`), and its `target_source` label, which is
  `"none"` for a `stack_psis` fit.

- `has_draws` and `draw_dim`: whether the fit has calibrated draws
  ([`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
  does not return `NULL`), and their numbers of rows and columns
  (`c(0L, 0L)` when it has none).

- `diagnostic_keys`: the names of the list that
  [`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
  returns.

- `interval_note`: the text of the line beginning "interval note:", or
  `NA` when there is no such line.

- `estimates` and `diagnostics`: the estimate table and the diagnostics
  list.

- `schema_version`, `summary_schema_version`, `source_validation`,
  `source_reportability_fit` and `validation`: the format versions of
  the fit and of the summary; the checksum of the fit that the summary
  was made from; for a `stack_psis` fit, a copy of the fit without its
  stacked draws and weights, from which `estimates` and `diagnostics`
  are taken (`NULL` for the other methods); and a SHA-256 checksum of
  the summary, which [`print()`](https://rdrr.io/r/base/print.html)
  recomputes to detect changes.

## Details

The overview has one line each for the method, the status, the number of
fixed effects, the target (its `target_source` label, or "none" when
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
returns `NULL`, as for a `stack_psis` fit), the calibrated draws (their
numbers of rows and columns, or "not retained" when
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
returns `NULL`) and the names of the diagnostics. A line beginning
"interval note:" follows when some or all intervals are descriptive, and
the reason codes and the number of warnings are shown when there are
any. The [`print()`](https://rdrr.io/r/base/print.html) method of the
summary adds the columns `term`, `estimate`, `se`, `df`, `conf_low` and
`conf_high` of the estimate table.

The three methods stop with an error if the fit or the summary was
changed after it was created, and for a fit saved by pvstackr 0.1.x,
which has to be made again. A summary of a `stack_psis` fit that was
saved by an earlier version of pvstackr also stops with an error when
printed (see
[`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)).

## See also

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
to read a fit.
