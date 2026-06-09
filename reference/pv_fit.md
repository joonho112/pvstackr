# Fit a pvstackr Method

`pv_fit()` is the generic public fitting entry point. In this package
stage, `method = "stack_direct"`, `method = "per_pv"`, and
`method = "stack_psis"` are implemented.

## Usage

``` r
pv_fit(
  data,
  formula,
  target = NULL,
  method = "stack_direct",
  control = NULL,
  ...
)
```

## Arguments

- data:

  Analysis data frame.

- formula:

  Two-sided formula with `OUTCOME` on the left-hand side.

- target:

  A method-specific target object. For `method = "stack_direct"`, this
  must be a `pvstackr_brr_target` from
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  and is required. Ignored for `"per_pv"` and `"stack_psis"`. Default
  `NULL`.

- method:

  Public method identifier. Character scalar; one of `"stack_direct"`,
  `"stack_psis"`, or `"per_pv"`. Default `"stack_direct"`.

- control:

  Optional
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  object. If `NULL`, one is built with `pv_control(method = method)`. If
  supplied, `control$method` must equal `method`. Default `NULL`.

- ...:

  Additional arguments forwarded to the dispatched fitter:
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  for `"stack_direct"`,
  [`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)
  for `"per_pv"`, or
  [`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)
  for `"stack_psis"`.

## Value

A `pvstackr_fit` object. Read it with the accessors
([`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md))
rather than by `$`-indexing; user-facing components include `estimates`,
`draws`, `target`, `diagnostics`, `status`, `reason_codes`, `warnings`,
and `method`.

## Details

### Method dispatch

`pv_fit()` validates `method` and `control`, then forwards to a
method-specific fitter:

- `"stack_direct"` (default) dispatches to
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  and **requires** a `pvstackr_brr_target` `target` from
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md);
  `pv_fit()` errors if `target` is `NULL`.

- `"per_pv"` dispatches to
  [`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)
  (the orthodox per-PV reference); `target` is ignored.

- `"stack_psis"` dispatches to
  [`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)
  (the PSIS-reweighted stacked path); `target` is ignored.

If `control` is `NULL` it is constructed with
`pv_control(method = method)`; otherwise `control$method` **must** equal
`method` or the call errors. Any arguments in `...` are forwarded
verbatim to the dispatched fitter, so consult that fitter's signature
for the available extras (for example, the injected backend adapters of
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
or the injected PSIS weights of
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)).

In this package stage, only `stack_direct` output is coverage-claimable:
its intervals are backed by the external Rubin/BRR-Fay target. `per_pv`
and `stack_psis` intervals are descriptive/reference, even with
Barnard-Rubin degrees of freedom.

## See also

[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md),
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md);
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)

Other pvstackr-fitting:
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)

## Examples

``` r
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

# Build the design and the external BRR-Fay target on the fixture.
design <- pv_design(
  pisa_tiny, formula = OUTCOME ~ x + female,
  pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
)
target <- pv_brr_target(
  pisa_tiny, formula = OUTCOME ~ x + female,
  pv_cols = design$pv_cols, weight_col = design$weight_col,
  rep_weight_cols = design$rep_weight_cols, fay_k = design$fay_k,
  id_cols = design$id_cols
)

# A live stack_direct fit needs an injected/precomputed backend:
if (FALSE) { # \dontrun{
fit <- pv_fit(
  data = pisa_tiny, formula = OUTCOME ~ x + female,
  target = target, method = "stack_direct",
  fit_function = my_backend_adapter, draws_function = my_draws_adapter
)
} # }

# Inspect the reportable object surface via the bundled cached fit instead.
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
if (nzchar(path)) {
  fit <- readRDS(path)$fit     # a stack_direct pvstackr_fit
  fit                          # compact console print
  head(get_estimates(fit))     # reportable fixed-effect table
}
#>          term   estimate        se std.error       df df_method df_complete
#> 1 b_Intercept 457.894088 1.2873118 1.2873118 1.021194   classic          NA
#> 2         b_x  46.883361 0.3717929 0.3717929 1.402308   classic          NA
#> 3    b_female   2.143702 3.5550309 3.5550309 1.013730   classic          NA
#>   conf_level  conf_low conf_high  conf.low conf.high             interval_role
#> 1       0.95 442.31804 473.47013 442.31804 473.47013 descriptive_classic_rubin
#> 2       0.95  44.41457  49.35215  44.41457  49.35215 descriptive_classic_rubin
#> 3       0.95 -41.60687  45.89428 -41.60687  45.89428 descriptive_classic_rubin
#>   coverage_claim_allowed parameter_scope          target_source target_hash
#> 1                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
#> 2                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
#> 3                  FALSE    fixed_effect external_brr_fay_rubin    4a4d40f8
```
