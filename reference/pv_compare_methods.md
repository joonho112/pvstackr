# Compare pvstackr Method Fits

`pv_compare_methods()` compares already-created `pvstackr_fit` objects
by aligning their fixed-effect estimate tables, computing differences
against a reference method, and recording agreement and timing
diagnostics.

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

  Two or more `pvstackr_fit` objects, or a single list of fits.

- fits:

  Optional explicit list of `pvstackr_fit` objects.

- reference_method:

  Optional method label or method id used as the comparison reference.
  Defaults to the first non-blocked `per_pv` fit when available,
  otherwise the first non-blocked fit.

- timings:

  Optional named numeric vector of elapsed seconds, named by method
  label or method id.

- include_fits:

  Whether to retain the original fit objects in the comparison object.

## Value

A `pvstackr_method_comparison` object with top-level `estimate_table`,
`diagnostic_table`, `agreement`, and `timing` fields. `estimate_table`
includes method-level interval and PSIS provenance metadata alongside
the aligned fixed-effect estimates. `diagnostic_table` includes
method-level target/pooling provenance, PSIS source and bounded
weight-concentration summaries, and shared-provenance flags.

## Details

The returned object keeps blocked methods visible. Blocked fits receive
rows for every aligned fixed-effect term, but their reportable numeric
comparison fields are `NA` and their reason codes remain available.

The comparison estimate table preserves interval metadata from each fit:
`df_method`, `df_complete`, `conf_level`, `interval_role`, and
`coverage_claim_allowed`. Default print methods keep their tables
compact and emit a one-line note when one or more compared intervals are
descriptive rather than coverage-claimable. In the method-level
diagnostic table, `coverage_claim_allowed` is `TRUE` only when all
available reportable rows for that method are coverage-claimable,
`FALSE` when any available row is descriptive, and `NA` when the method
has no available estimate rows.

Method diagnostics also carry target and pooling provenance where
available: `target_source`, `target_hash`, `pooling_source`,
`pooling_hash`, and shared hash flags. The nested
`diagnostics$target_overlap` summary records whether any compared
methods share external targets, target hashes, pooling hashes,
reference-target hashes, or target-source families. Agreement
diagnostics are descriptive; close agreement is not automatic
independent corroboration when compared methods share a target hash,
pooling hash, target source, or estimand construction. `target_source`
is provenance vocabulary, not necessarily a formal target object:
`stack_psis` rows use `stack_psis_rubin_pooling` even though
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
returns `NULL` for `stack_psis` fits. For `stack_psis`, the aligned
estimate rows additionally retain `psis_source`, `pareto_k_source`,
`weight_method`, `psis_producer`, and `psis_producer_version`. The
method diagnostic table adds the minimum per-PV Kish-style iid weight
ESS and ESS fraction, plus the maximum normalized weight. These fields
preserve the source fit's bounded provenance and concentration
diagnostics; they are not MCMC ESS or a substitute for Pareto-k gating.

Current comparison objects record deep-valid compact PSIS source fits,
canonical per-source reportability projections, each source authority's
validation stamp, and an owned-payload SHA-256 stamp. Derived tables and
diagnostics are recomputed from those projections during validation. A
pre-marker serialized comparison containing `stack_psis` is
inspection-only and must be rebuilt from current validated fits; its
saved numeric table is not grandfathered. Independently of the stamp,
warning-status PSIS rows are forbidden and blocked PSIS rows must carry
no reportable numeric or pooling metadata.

## Interpreting agreement

Agreement diagnostics are **descriptive**, not independent
corroboration. Close agreement does **not** confirm a result when the
compared methods share a target hash, pooling hash, target source, or
estimand construction: the methods are then reading the same
information, so concordance is expected and uninformative about
correctness. The nested `diagnostics$target_overlap` summary flags this,
and the print methods emit a one-line provenance note when overlap is
detected. Only `stack_direct` rows backed by the external Rubin/BRR-Fay
target are coverage-claimable; `per_pv` and `stack_psis` rows are
descriptive/reference even with Barnard-Rubin degrees of freedom.

## See also

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)

## Examples

``` r
path <- system.file(
  "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
)
if (nzchar(path)) {
  fit_direct <- readRDS(path)$fit          # cached stack_direct pvstackr_fit

  # Build a per_pv reference fit sharing the same fixed-effect names
  # (b_Intercept, b_x, b_female) from small injected posterior draw clouds.
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
  cmp                          # descriptive agreement + provenance notes
  head(get_estimates(cmp))     # aligned per-method estimate rows
}
#>         method method_label        term status     estimate        se std.error
#> 1 stack_direct stack_direct b_Intercept     ok 457.89408827 1.2873118 1.2873118
#> 2 stack_direct stack_direct         b_x     ok  46.88336051 0.3717929 0.3717929
#> 3 stack_direct stack_direct    b_female     ok   2.14370179 3.5550309 3.5550309
#> 4       per_pv       per_pv b_Intercept     ok  -0.03206137 1.0200128 1.0200128
#> 5       per_pv       per_pv         b_x     ok   0.02385441 1.0406796 1.0406796
#> 6       per_pv       per_pv    b_female     ok  -0.06147476 1.0401232 1.0401232
#>             df df_method df_complete conf_level             interval_role
#> 1 1.021194e+00   classic          NA       0.95 descriptive_classic_rubin
#> 2 1.402308e+00   classic          NA       0.95 descriptive_classic_rubin
#> 3 1.013730e+00   classic          NA       0.95 descriptive_classic_rubin
#> 4 5.759265e+03   classic          NA       0.95   reference_classic_rubin
#> 5 1.642548e+06   classic          NA       0.95   reference_classic_rubin
#> 6 8.726226e+05   classic          NA       0.95   reference_classic_rubin
#>   coverage_claim_allowed psis_source pareto_k_source weight_method
#> 1                  FALSE        <NA>            <NA>          <NA>
#> 2                  FALSE        <NA>            <NA>          <NA>
#> 3                  FALSE        <NA>            <NA>          <NA>
#> 4                  FALSE        <NA>            <NA>          <NA>
#> 5                  FALSE        <NA>            <NA>          <NA>
#> 6                  FALSE        <NA>            <NA>          <NA>
#>   psis_producer psis_producer_version   conf_low  conf_high reference_method
#> 1          <NA>                  <NA> 442.318045 473.470132           per_pv
#> 2          <NA>                  <NA>  44.414570  49.352151           per_pv
#> 3          <NA>                  <NA> -41.606872  45.894276           per_pv
#> 4          <NA>                  <NA>  -2.031670   1.967547           per_pv
#> 5          <NA>                  <NA>  -2.015842   2.063551           per_pv
#> 6          <NA>                  <NA>  -2.100082   1.977132           per_pv
#>   reference_estimate reference_se estimate_diff  se_ratio  abs_z_diff
#> 1        -0.03206137     1.020013    457.926150 1.2620545 278.8090409
#> 2         0.02385441     1.040680     46.859506 0.3572597  42.4029899
#> 3        -0.06147476     1.040123      2.205177 3.4178942   0.5953395
#> 4        -0.03206137     1.020013      0.000000 1.0000000   0.0000000
#> 5         0.02385441     1.040680      0.000000 1.0000000   0.0000000
#> 6        -0.06147476     1.040123      0.000000 1.0000000   0.0000000
#>   agreement_band reason_codes
#> 1      different             
#> 2      different             
#> 3      different             
#> 4          close             
#> 5          close             
#> 6          close             
```
