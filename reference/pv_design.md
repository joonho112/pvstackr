# Declare a PISA-Style Plausible-Value Design

`pv_design()` records plausible-value columns, final weights, BRR
replicate weights, Fay coefficient, row-support metadata, and stable
hashes without running a model or assembling a target. It is the
applied-user wrapper around the internal `pvstackr_design` object.

## Usage

``` r
pv_design(
  data,
  formula,
  pv_cols = NULL,
  weight_col = "W_FSTUWT",
  rep_weight_cols = NULL,
  fay_k = 0.5,
  pv_prefix = "PV",
  pv_suffix = "",
  rep_weight_prefix = "W_FSTURWT",
  expected_M = NULL,
  expected_R = NULL,
  id_cols = NULL,
  roles = list()
)

# S3 method for class 'pvstackr_design'
print(x, ...)
```

## Arguments

- data:

  Data frame containing plausible-value, weight, and replicate weight
  columns. Must have unique, non-empty column names.

- formula:

  Formula using `OUTCOME` as the plausible-value placeholder, for
  example `OUTCOME ~ x + female`. Do not embed
  [`weights()`](https://rdrr.io/r/stats/weights.html); pass weight
  columns explicitly.

- pv_cols:

  Optional character vector of plausible-value columns. If `NULL`
  (default), columns are detected with
  [`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)
  using `pv_prefix`/`pv_suffix`.

- weight_col:

  Character scalar naming the final survey-weight column. Defaults to
  PISA's `"W_FSTUWT"`.

- rep_weight_cols:

  Optional character vector of BRR replicate-weight columns. If `NULL`
  (default), columns are detected with
  [`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md)
  using `rep_weight_prefix`.

- fay_k:

  Numeric scalar Fay coefficient used by the BRR replicate design. Must
  satisfy `0 <= fay_k < 1`. Default `0.5`.

- pv_prefix, pv_suffix:

  Character scalars giving the prefix and suffix used when detecting
  plausible values. Modern PISA files commonly use subject-suffixed
  plausible values such as `PV1MATH`; set `pv_suffix = "MATH"` or the
  relevant subject suffix instead of relying on the bare `PV1`, `PV2`,
  ... default. Defaults `"PV"` and `""`.

- rep_weight_prefix:

  Character scalar prefix used when detecting replicate weights. Default
  `"W_FSTURWT"`.

- expected_M:

  Optional integer scalar. Expected plausible-value count; checked
  against detected or supplied `pv_cols`. If `NULL` (default), the count
  is not checked.

- expected_R:

  Optional integer scalar. Expected replicate-weight count; checked
  against detected or supplied `rep_weight_cols`. If `NULL` (default),
  the count is not checked.

- id_cols:

  Optional character vector of columns that jointly identify unique
  rows. If `NULL` (default), row support is tracked by row number.

- roles:

  Optional named list of extra design roles, merged over the defaults
  recorded by the wrapper. Default
  [`list()`](https://rdrr.io/r/base/list.html).

- x:

  A `pvstackr_design` object.

- ...:

  Ignored.

## Value

A `pvstackr_design` object: a list recording the declared design and its
provenance. User-facing fields include:

- `n`, `M`, `R`:

  Row count, number of plausible values, and number of BRR replicate
  weights.

- `formula`, `formula_string`:

  The design formula (with the `OUTCOME` placeholder) and its string
  form.

- `pv_cols`, `weight_col`, `rep_weight_cols`:

  Resolved plausible-value, final-weight, and replicate-weight column
  names.

- `fay_k`:

  Fay coefficient for the BRR replicate design.

- `id_cols`, `row_support`:

  Row-identifier columns and the row-support descriptor.

- `design_hash`, `row_support_hash`, `pv_value_hash`,
  `weight_design_hash`:

  Stable content hashes of the design and its components.

- `roles`, `provenance`:

  Design roles and wrapper provenance.

The compact [`print()`](https://rdrr.io/r/base/print.html) method shows
`n`, `formula_string`, `M`, `weight_col`, `R`, `fay_k`, and
`design_hash`.

## See also

Detect inputs with
[`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)
and
[`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md);
assemble the external target with
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md);
fit a method with
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md).

## Examples

``` r
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)
design <- pv_design(
  pisa_tiny, formula = OUTCOME ~ x + female,
  pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
)
design
#> pvstackr design
#>   rows: 12
#>   formula: OUTCOME ~ x + female
#>   plausible values: 2
#>   final weight: W_FSTUWT
#>   replicate weights: 4
#>   fay_k: 0.5
#>   design hash: fa78c04b
```
