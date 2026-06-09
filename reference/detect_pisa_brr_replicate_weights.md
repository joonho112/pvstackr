# Detect PISA-Style BRR Replicate-Weight Columns

PISA balanced-repeated-replication (BRR) replicate weights are named as
a prefix followed by a numeric replicate index, for example
`W_FSTURWT1`, `W_FSTURWT2`, and so on. This finds the `<prefix><r>`
columns and returns them in natural numeric order using the default
prefix `"W_FSTURWT"`. It is the replicate-weight companion to
[`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md).
Use `expected_R` to guard against an incomplete replicate set or the
wrong prefix.

## Usage

``` r
detect_pisa_brr_replicate_weights(
  data,
  prefix = "W_FSTURWT",
  expected_R = NULL
)
```

## Arguments

- data:

  Data frame containing replicate-weight columns.

- prefix:

  Character scalar. Column prefix before the numeric replicate index.
  Default `"W_FSTURWT"`.

- expected_R:

  Optional integer scalar. Expected replicate-weight count; if supplied,
  detection errors unless exactly `expected_R` columns are found. If
  `NULL` (default), the count is not checked.

## Value

Character vector of column names in natural numeric order.

## See also

[`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md);
build a design with
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
or a target with
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).

Other pvstackr-detection:
[`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)

## Examples

``` r
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)
detect_pisa_brr_replicate_weights(pisa_tiny, expected_R = 4L)
#> [1] "W_FSTURWT1" "W_FSTURWT2" "W_FSTURWT3" "W_FSTURWT4"
```
