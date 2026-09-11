# Find PISA-style plausible-value columns

`detect_pisa_pv_columns()` returns the names of the columns that hold
the plausible values of one subject, such as `PV1READ`, `PV2READ`, ...,
`PV10READ`. A column matches when its whole name is `prefix`, a number
and `suffix`; the match is case-sensitive.
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
and
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
call this function when you do not give `pv_cols`.

## Usage

``` r
detect_pisa_pv_columns(data, prefix = "PV", suffix = "", expected_M = NULL)
```

## Arguments

- data:

  A data frame. Only its column names are used.

- prefix:

  The text before the number. Default `"PV"`.

- suffix:

  The text after the number, such as `"MATH"` for `PV1MATH`, `PV2MATH`,
  and so on. The default `""` matches only bare names such as `PV1`.

- expected_M:

  The number of plausible-value columns you expect, a whole number. If a
  different number is found, the function stops with an error, which
  catches a wrong suffix or an incomplete set of columns. `NULL`
  (default) skips the check.

## Value

A character vector of column names in numeric order: `PV2READ` comes
before `PV10READ`, unlike in alphabetical order.

## Details

PISA files name the plausible values by subject, for example `PV1MATH`
or `PV1READ`, so give the subject as the suffix: `suffix = "MATH"` here,
or `pv_suffix = "MATH"` in
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
and
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
The default `suffix = ""` matches only bare names such as `PV1` and
`PV2`, never `PV1MATH`; if the data also contain subject-suffixed
plausible values, the function warns.

The numbers must run from 1 without gaps or repeats (`PV1` and `PV01`
both count as 1); otherwise the function stops with an error. It also
stops when no column matches, and the error message then lists any
subject suffixes found in the data.

## See also

[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
and
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
which call this function.

Other pvstackr-detection:
[`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md)

## Examples

``` r
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)
# Subject-suffixed plausible values: pass the subject suffix explicitly.
detect_pisa_pv_columns(pisa_tiny, suffix = "READ", expected_M = 2L)
#> [1] "PV1READ" "PV2READ"
```
