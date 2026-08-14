# Fit the Per-PV Bayesian/Backend Reference Method

`pv_fit_reference()` implements the `per_pv` reference path: one backend
fit or posterior-draw source per plausible value, fixed-effect summary
extraction, and Rubin pooling with model-based within-PV covariance from
the selected fixed-effect draws.

## Usage

``` r
pv_fit_reference(
  data = NULL,
  formula = NULL,
  pv_cols = NULL,
  per_pv_draws = NULL,
  control = pv_control(method = "per_pv"),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  param_map = NULL,
  weight_col = NULL,
  rep_weight_cols = NULL,
  fay_k = 0.5,
  id_cols = NULL,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL,
  cache_dir = "cache",
  cache_stem = "pvstackr-per-pv",
  additional_args = list()
)
```

## Arguments

- data:

  Optional analysis data frame. Required for the injected fitting route.

- formula:

  Optional two-sided formula with `OUTCOME` on the left-hand side.
  Required for the injected fitting route.

- pv_cols:

  Plausible-value columns. Required for the injected fitting route;
  optional alignment labels for `per_pv_draws`.

- per_pv_draws:

  Optional list of per-PV posterior draw matrices.

- control:

  A
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  object with `method = "per_pv"`.

- family:

  Optional backend family object passed to the injected fit function.

- prior:

  Optional backend prior object passed to the injected fit function.

- fit_function:

  Injected backend fitting function.

- draws_function:

  Function extracting posterior draws from each injected backend fit.

- param_map:

  Optional explicit draw-column map passed to the fixed-effect
  extraction layer. Use this when backend draw columns do not follow the
  package's automatic `b_*` fixed-effect naming convention. Supply
  `fe_names` or `fe_idx` to identify fixed-effect columns, and optional
  `vc_names` or `vc_idx` for nuisance variance-component columns. Use
  `vc_names = character()` to drop all nuisance columns, including
  distributional names such as `b_sigma_*`.

- weight_col:

  Optional main weight column recorded in the design object and
  validated when the injected fitting route is used. It is not passed as
  an automatic backend argument.

- rep_weight_cols:

  Optional replicate-weight columns recorded in the design object. They
  are not used by `per_pv` pooling and are not passed as automatic
  backend arguments.

- fay_k:

  Fay coefficient recorded in the design object when replicate weights
  are supplied.

- id_cols:

  Optional row identifier columns recorded in the design object.

- df_method:

  Rubin degrees-of-freedom method. The default `"classic"` is
  descriptive only; `"barnard_rubin"` requires `df_complete`.

- df_complete:

  Complete-data degrees of freedom used when
  `df_method = "barnard_rubin"`. Supply a positive scalar for all fixed
  effects or a named numeric vector aligned by fixed-effect name.
  Unnamed length-p vectors are rejected to avoid positional ambiguity.

- cache_dir, cache_stem:

  Cache location metadata passed to injected fits.

- additional_args:

  Additional named arguments passed to each injected fit.

## Value

A `pvstackr_fit` object with `method = "per_pv"`.

## Details

This reference method accepts either already-computed per-PV posterior
draws via `per_pv_draws`, or an injected `fit_function`/`draws_function`
pair. The injected route rewrites the `OUTCOME` placeholder to each
plausible-value column and calls the backend once per plausible value.
Current reportable output is fixed-effect-only, selected from draw
columns named like `b_*` or by an explicit `param_map`. The within-PV
variance component is the posterior draw covariance of the selected
fixed-effect columns in each plausible value. `weight_col`,
`rep_weight_cols`, `fay_k`, and `id_cols` are recorded as design
metadata for this path; BRR/Fay replicate weights are not used to form
that within-PV covariance or an external target. These weight arguments
are not passed automatically to an injected backend. If a backend fit
should use weights, handle that explicitly inside `fit_function`, in the
data supplied to the backend, or through backend-specific
`additional_args`.

`per_pv` interval metadata is source-specific. Classic Rubin pooling is
labeled `interval_role = "reference_classic_rubin"`; Barnard-Rubin
pooling is labeled `interval_role = "reference_barnard_rubin"` and uses
the explicit `df_complete` supplied by the user. Both roles set
`coverage_claim_allowed = FALSE` because the within-PV covariance is
model-based posterior-draw covariance, not an external design replicate
variance.

## Reportable scope and coverage

In this package stage, reportable output is **fixed-effect-only**;
variance components are fit but not calibrated to the target. Coverage
claims are enabled **only** for `stack_direct` rows backed by the
external Rubin/BRR-Fay target
(`interval_role = "coverage_barnard_rubin"`,
`coverage_claim_allowed = TRUE`). `per_pv` intervals are
**descriptive/reference** even with Barnard-Rubin degrees of freedom.
"One stacked fit" describes the computational topology, not a
benchmarked speed claim. The within-PV variance is the **model-based**
posterior covariance of the selected fixed-effect draws, not a BRR/Fay
replicate variance.

## References

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.*
Wiley.

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948-955.

## See also

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md),
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md);
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md),
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).

Other pvstackr-fitting:
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)

## Examples

``` r
# Injected per-PV draws: one posterior draw matrix per plausible value, with
# fixed-effect columns named like `b_*` (auto-selected). Real fits would come
# from a backend; here we use small synthetic draw clouds to show the shape.
set.seed(1)
make_draws <- function(n, b0, bx) {
  cbind(
    b_Intercept = rnorm(n, b0, 0.5),
    b_x         = rnorm(n, bx, 0.2)
  )
}
per_pv_draws <- list(
  PV1READ = make_draws(200, 1.0, 0.30),
  PV2READ = make_draws(200, 1.1, 0.28)
)
fit_ref <- pv_fit_reference(
  per_pv_draws = per_pv_draws,
  pv_cols      = c("PV1READ", "PV2READ"),
  control      = pv_control(method = "per_pv")
)
fit_ref                     # a per_pv pvstackr_fit
#> pvstackr fit
#>   method: per_pv
#>   status: ok
#>   fixed effects: 2
#>   target: per_pv_rubin_draws
#>   draws: not retained
#>   diagnostics: reference, pooling
#>   interval note: intervals are descriptive rather than coverage-claimable.
get_estimates(fit_ref)      # interval_role = "reference_*"; coverage_claim_allowed = FALSE
#>          term  estimate        se std.error        df df_method df_complete
#> 1 b_Intercept 1.0484282 0.5038295 0.5038295 8103.9103   classic          NA
#> 2         b_x 0.2840975 0.2143992 0.2143992  704.1008   classic          NA
#>   conf_level    conf_low conf_high    conf.low conf.high
#> 1       0.95  0.06079298 2.0360634  0.06079298 2.0360634
#> 2       0.95 -0.13684085 0.7050359 -0.13684085 0.7050359
#>             interval_role coverage_claim_allowed parameter_scope
#> 1 reference_classic_rubin                  FALSE    fixed_effect
#> 2 reference_classic_rubin                  FALSE    fixed_effect
#>        target_source target_hash     pooling_source pooling_hash
#> 1 per_pv_rubin_draws    c226f7bb per_pv_rubin_draws     c226f7bb
#> 2 per_pv_rubin_draws    c226f7bb per_pv_rubin_draws     c226f7bb
```
