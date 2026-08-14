# Revalidate a Legacy BRR-Fay Target Against Its Original Inputs

`pv_revalidate_brr_target()` upgrades a schema-0.1 BRR-Fay target only
after rebuilding its model, replicate-weight target, target content, and
binding manifest from the complete original data and formula. Legacy
objects without those inputs remain inspection-only and cannot be used
for reportable fitting.

## Usage

``` r
pv_revalidate_brr_target(
  target,
  data = NULL,
  formula = NULL,
  conf_level = NULL
)
```

## Arguments

- target:

  A `pvstackr_brr_target` object.

- data:

  The complete original data frame used to create a schema-0.1 target.

- formula:

  The original model formula.

- conf_level:

  Confidence level used by the legacy target. This must be supplied when
  it is not retained on the source object.

## Value

A separately allocated, validated schema-0.2 `pvstackr_brr_target`, or
the unchanged schema-0.2 input.

## Details

Already bound schema-0.2 targets are strictly validated and returned
unchanged.
