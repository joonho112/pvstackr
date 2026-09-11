# Find PISA-style replicate-weight columns

`detect_pisa_brr_replicate_weights()` returns the names of the
replicate-weight columns, such as `W_FSTURWT1`, ..., `W_FSTURWT80` in
PISA 2022 files, where they are balanced repeated replication (BRR)
weights made with Fay's method. A column matches when its whole name is
`prefix` followed by a number; the match is case-sensitive.
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
and
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
call this function when you do not give `rep_weight_cols`.

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

  A data frame. Only its column names are used.

- prefix:

  The text before the number. Default `"W_FSTURWT"`, the PISA name.

- expected_R:

  The number of replicate-weight columns you expect, a whole number, for
  example `80` for PISA 2022. If a different number is found, the
  function stops with an error. `NULL` (default) skips the check.

## Value

A character vector of column names in numeric order (`W_FSTURWT1`,
`W_FSTURWT2`, ..., `W_FSTURWT80`).

## Details

The numbers must run from 1 without gaps or repeats; otherwise, or when
no column matches, the function stops with an error. Only the names are
checked:
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
and
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
check the weights themselves.

## See also

[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
and
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
which call this function.

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
