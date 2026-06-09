# Detect PISA-Style Plausible-Value Columns

Modern PISA plausible values are often subject-suffixed, for example
`PV1MATH`, `PV2MATH`, or `PV1READ`, `PV2READ`. Pass the subject/domain
suffix explicitly, such as `suffix = "MATH"` here or
`pv_suffix = "MATH"` in
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
and
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
The default `suffix = ""` is reserved for data whose plausible values
are intentionally named as bare `PV1`, `PV2`, and so on. Detection is
anchored, so `suffix = ""` does not match `PV1MATH`. Use `expected_M` to
guard against selecting the wrong subject or an incomplete
plausible-value set.

## Usage

``` r
detect_pisa_pv_columns(data, prefix = "PV", suffix = "", expected_M = NULL)
```

## Arguments

- data:

  Data frame containing plausible-value columns.

- prefix:

  Character scalar. Column prefix before the numeric plausible-value
  index. Default `"PV"`.

- suffix:

  Character scalar. Column suffix after the numeric plausible-value
  index, such as `"MATH"` for columns named `PV1MATH`, `PV2MATH`, and so
  on. Default `""` matches bare `PV1`, `PV2`, ... only; detection is
  anchored, so the default does not match subject-suffixed columns.

- expected_M:

  Optional integer scalar. Expected plausible-value count; if supplied,
  detection errors unless exactly `expected_M` columns are found. If
  `NULL` (default), the count is not checked.

## Value

Character vector of column names in natural numeric order.

## See also

[`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md);
build a design with
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
or a target with
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).

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
