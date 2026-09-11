# Compare fits made with different methods

`pv_compare_methods()` lines up the fixed-effect estimates of two or
more fits, for example a `stack_direct` fit and a `per_pv` fit of the
same model, and computes how far each fit's estimates are from those of
a reference fit. It also records the status, the kind of interval and
the source labels of each fit, and the time that each fit took if you
supply it. It does not fit any model itself.

## Usage

``` r
pv_compare_methods(
  ...,
  fits = NULL,
  reference_method = NULL,
  timings = NULL,
  include_fits = FALSE
)
```

## Arguments

- ...:

  Two or more fits (class `pvstackr_fit`) from
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
  or the method functions, optionally named, or a single list of fits.

- fits:

  A list of fits, as an alternative to `...`; default `NULL`. Give the
  fits either in `...` or in `fits`, not both.

- reference_method:

  The label or the method of the reference fit, or `NULL` (default) to
  choose it by the rule in Details. It must match exactly one fit, and
  that fit must not be blocked.

- timings:

  Elapsed times in seconds that you measured for the fits, for example
  with [`system.time()`](https://rdrr.io/r/base/system.time.html): a
  numeric vector or list named by fit label or method, with finite
  values of 0 or more. pvstackr does not time the fits itself. With
  `NULL` (default), and for a fit whose label and method are not among
  the names, the time is `NA`.

- include_fits:

  `TRUE` to keep the compared fits in the element `fits` of the result;
  default `FALSE`.

## Value

A list of class `pvstackr_method_comparison`. Read it with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md);
its [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods are described
in
[pvstackr_method_comparison_summary](https://joonho112.github.io/pvstackr/reference/pvstackr_method_comparison_summary.md).

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
returns the estimate table, a data frame with one row per fit and fixed
effect:

- `method`, `method_label`, `term` and `status`: the method and the
  label of the fit, the fixed effect and the status of the fit.

- `estimate`, `se` (repeated as `std.error`), `df`, `df_method`,
  `df_complete`, `conf_level`, `conf_low`, `conf_high`, `interval_role`
  and `coverage_claim_allowed`: copied from the estimate table of the
  fit
  ([pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
  describes them).

- `psis_source`, `pareto_k_source`, `weight_method`, `psis_producer` and
  `psis_producer_version`: how the weights of a `stack_psis` fit were
  obtained; `NA` for the other methods.

- `reference_method`, `reference_estimate` and `reference_se`: the label
  of the reference fit and its estimate and standard error.

- `estimate_diff`, `se_ratio`, `abs_z_diff` and `agreement_band`: the
  comparison with the reference fit (see "Reading the agreement").

- `reason_codes`: the reason codes of the fit, separated by commas.

[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
returns a list:

- `reference_method`: the label of the reference fit.

- `methods` and `statuses`: the method and the status of each fit, named
  by label.

- `blocked_methods` and `warning_methods`: the labels of the fits with
  status `"blocked"` or `"warning"`.

- `agreement`: a data frame with one row per fit that gives
  `max_abs_z_diff` and `max_abs_log_se_ratio` (the largest `abs_z_diff`
  and the largest absolute log of `se_ratio` over its fixed effects),
  `n_terms` (the number of fixed effects in the comparison) and
  `n_available` (the number of them with an estimate from this fit).

- `method_diagnostics`: a data frame with one row per fit. Its columns
  are `method`, `method_label`, `status`, `reason_codes` and
  `warning_count` (the number of warnings); `n_terms` and `n_available`;
  `df_method`, `interval_role`, `coverage_claim_allowed` and
  `n_descriptive_intervals` (the number of rows with
  `coverage_claim_allowed = FALSE`); the source labels and checksums
  `target_source`, `target_hash`, `pooling_source` and `pooling_hash`
  (see "Source labels" in
  [pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md));
  `center` (`"target"` for a `stack_direct` fit, `NA` otherwise); the
  five flags of `target_overlap` for this fit; `psis_status`,
  `pareto_k_max`, the five weight-source columns of the estimate table,
  `weight_ess_iid_min`, `weight_ess_fraction_min` and
  `max_normalized_weight_max` (the smallest or largest value over the
  plausible values), which describe a `stack_psis` fit and are `NA` for
  the other methods; and `n_fits` and `elapsed_seconds`. Here
  `coverage_claim_allowed` is `TRUE` only when all available rows of the
  fit have `TRUE`, `FALSE` when any has `FALSE`, and `NA` for a fit
  without estimates. `target_source` is only a label: a `stack_psis` fit
  has no target object.

- `timing`: a data frame with one row per fit that gives `n_fits`, the
  number of models fitted (the number of plausible values for `per_pv`,
  1 for `stack_direct` and `stack_psis`, `NA` for a blocked
  `stack_direct` fit), and `elapsed_seconds`, the time from `timings`
  (`NA` without one).

- `target_overlap`: whether fits share a target or a combined result
  (see "Reading the agreement").

The result also holds this list as `diagnostics`; the tables as
`estimate_table` (and `table`), `diagnostic_table` (the same as
`method_diagnostics`), `agreement` and `timing`; `reference_method`,
`methods` and `method_labels`; and the compared fits as `fits` when
`include_fits = TRUE` (otherwise `NULL`). Its other elements are
described in "Technical details".

## Details

Each fit is known by a label: its argument name, such as `per_pv` in
`pv_compare_methods(stack_direct = fit_a, per_pv = fit_b)`, or its
method when it has no name. Repeated labels get a suffix, as in
`"per_pv_1"`. The estimate table has one row per fit and fixed effect,
for every fixed effect that appears in any of the fits; a fit without
that fixed effect has `NA` in the row.

Differences are taken from the reference fit, which `reference_method`
names. By default it is the first `per_pv` fit that is not blocked or,
if there is none, the first fit that is not blocked. A blocked fit
cannot be the reference, but it keeps its rows, with `NA` values and its
reason codes, so that it stays visible in the comparison.

`pv_compare_methods()` stops with an error if it gets fewer than two
fits, an object that is not a fit or a fit that was changed after it was
created, and when all fits are blocked.

## Reading the agreement

`estimate_diff` is the estimate minus the reference estimate, `se_ratio`
is the standard error divided by the reference standard error, and
`abs_z_diff` is the absolute difference divided by
`sqrt(se^2 + reference_se^2)`. `agreement_band` groups `abs_z_diff` by
pvstackr's cut-offs: `"close"` below 0.1, `"moderate"` from 0.1 to below
0.5, and `"different"` from 0.5 on. It is `"not_available"` for a
blocked fit and where either estimate is missing. The rows of the
reference fit itself have a difference of 0 and are `"close"`.

These values describe how far apart the fits are; they are not a
statistical test. The fits usually come from the same data and plausible
values, and a `stack_psis` fit may reweight the draws of the same
stacked model, so agreement between them is not independent confirmation
of a result. Where fits share their numbers, agreement follows by
construction: a `stack_direct` fit reports the estimates and standard
errors of its target, so two `stack_direct` fits calibrated to the same
target have the same estimates and standard errors, whatever their
draws. `get_diagnostics(x)$target_overlap` records such cases:

- `shared_target_hash` and `shared_pooling_hash`: `TRUE` when two or
  more fits have the same `target_hash` (or `pooling_hash`).

- `shared_external_target`: `TRUE` only when two or more fits have
  `target_source = "external_brr_fay_rubin"` and the same `target_hash`,
  that is, when they are `stack_direct` fits calibrated to the same
  target from
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
  A blocked `stack_direct` fit keeps its target, so it counts as well.

- `shares_reference_target` and `shares_reference_pooling`: `TRUE` when
  a fit other than the reference fit has the `target_hash` (or
  `pooling_hash`) of the reference fit.

- `shared_target_sources`: the `target_source` labels that occur in two
  or more fits.

- `independence_caveat_required`: `TRUE` when any of these flags is
  `TRUE` or `shared_target_sources` is not empty.

- `independence_caveat`: the text of the line beginning "provenance
  note:", which [`print()`](https://rdrr.io/r/base/print.html) shows for
  every comparison.

Each row keeps the interval columns of its fit, so an interval is read
by the rule of that fit: only `stack_direct` rows whose target uses
Barnard-Rubin degrees of freedom have `coverage_claim_allowed = TRUE`
([`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
gives the rule). [`print()`](https://rdrr.io/r/base/print.html) adds a
line beginning "interval note:" when some or all intervals are
descriptive.

## Technical details

The result keeps what the comparison uses from each fit
(`source_fit_projection`), with the fit's checksum
(`source_fit_validation`) and, for a `stack_psis` fit, a copy of the fit
without its stacked draws and weights (`source_fit_reportability`).
`validation` holds a SHA-256 checksum of the comparison, `created_at`
the time it was made, `schema_version` its format version and
`provenance` a record of how it was made; `warnings` is empty. The
functions that read a comparison stop with an error if it was changed
after it was created or was saved by pvstackr 0.1.x; make it again from
current fits.

## See also

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
to make the fits,
[pvstackr_method_comparison_summary](https://joonho112.github.io/pvstackr/reference/pvstackr_method_comparison_summary.md)
for [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`vignette("a4-comparing-methods", package = "pvstackr")`](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.md)
for a comparison of the three methods.

## Examples

``` r
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
if (nzchar(path)) {
  fit_direct <- readRDS(path)$fit          # the bundled stack_direct fit

  # A per_pv fit made from random draws with the same fixed-effect names
  # (b_Intercept, b_x, b_female). It is not fitted to any data, so the
  # differences below only show the layout of the output.
  set.seed(1)
  fe_names <- c("b_Intercept", "b_x", "b_female")
  draw_block <- function() {
    matrix(rnorm(200 * 3), ncol = 3, dimnames = list(NULL, fe_names))
  }
  fit_per_pv <- pv_fit_reference(
    per_pv_draws = list(PV1READ = draw_block(), PV2READ = draw_block()),
    control      = pv_control(method = "per_pv")
  )

  cmp <- pv_compare_methods(stack_direct = fit_direct, per_pv = fit_per_pv)
  print(cmp)                   # reference fit, labels and the two notes
  get_estimates(cmp)[, c("method_label", "term", "estimate", "se",
                         "abs_z_diff", "agreement_band")]
}
#> pvstackr method comparison
#>   reference: per_pv
#>   methods: stack_direct=stack_direct, per_pv=per_pv
#>   fixed effects: 3
#>   provenance note: Agreement bands are descriptive; shared target, pooling, or source metadata should not be read as independent corroboration.
#>   interval note: intervals are descriptive rather than coverage-claimable.
#>   method_label        term     estimate        se  abs_z_diff agreement_band
#> 1 stack_direct b_Intercept 457.89408827 1.2873118 278.8090409      different
#> 2 stack_direct         b_x  46.88336051 0.3717929  42.4029899      different
#> 3 stack_direct    b_female   2.14370179 3.5550309   0.5953395      different
#> 4       per_pv b_Intercept  -0.03206137 1.0200128   0.0000000          close
#> 5       per_pv         b_x   0.02385441 1.0406796   0.0000000          close
#> 6       per_pv    b_female  -0.06147476 1.0401232   0.0000000          close
```
