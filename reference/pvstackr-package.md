# pvstackr: Bayesian Plausible-Value Analysis with One Calibrated Stacked Fit

pvstackr fits survey-weighted linear regression models to assessment
data whose outcome is given as plausible values, as in PISA, and reports
estimates, standard errors and intervals for the fixed effects. The
default method, `stack_direct`, fits one Bayesian model to the stacked
data (one copy of the data per plausible value) and then applies the
Cholesky calibration correction (CCC): the fixed-effect draws are
transformed to have the mean and covariance of a design-based target,
whose estimates, standard errors and degrees of freedom are the ones
reported. The other methods, `per_pv` (one fit per plausible value) and
`stack_psis` (one stacked fit reweighted toward each plausible value),
work from draws, fitting functions or importance weights that you
supply.

## Details

### Workflow

1.  Optionally,
    [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
    checks and records the plausible-value, final-weight and
    replicate-weight columns and the Fay coefficient. By default it
    finds PISA-style column names with
    [`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)
    and
    [`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md).
    This step can be skipped:
    [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
    takes the same column arguments and finds the PISA-style columns
    itself.

2.  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
    computes the target: a survey-weighted regression for each plausible
    value, the sampling covariance of the coefficients from the BRR-Fay
    replicate weights, and Rubin's rules across plausible values.

3.  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
    fits a method with settings from
    [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).
    Only `stack_direct` has a bundled engine, brms, selected with
    `pv_control(backend = "brms")`.

4.  Optionally,
    [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
    compares fits made with different methods.

5.  [`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
    [`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
    [`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
    and
    [`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
    read the results. A fit with `status` `"blocked"` has no estimates;
    `reason_codes` says why.

### Scope

pvstackr is built for the fixed effects of a survey-weighted linear
model.
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
`stack_direct` and `stack_psis` stop with an error on a random-effect
term such as `(1 | school)`. Only the fixed-effect coefficients
(intercept and slopes) are reported; other parameters, such as the
residual standard deviation `sigma`, are neither calibrated nor
reported.

By pvstackr's reporting rule, only a `stack_direct` fit whose target
uses Barnard-Rubin degrees of freedom (`df_method = "barnard_rubin"` in
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md))
has coverage-claimable intervals, which can be read as confidence
intervals with nominal coverage. All other intervals are descriptive,
including those of a `stack_direct` fit whose target uses the default
classic Rubin degrees of freedom.
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
gives the full rule.

### Articles

The Applied articles, starting with
[`vignette("a1-getting-started", package = "pvstackr")`](https://joonho112.github.io/pvstackr/articles/a1-getting-started.md),
show how to run an analysis and read the results. The Method articles,
starting with
[`vignette("m1-foundations-and-notation", package = "pvstackr")`](https://joonho112.github.io/pvstackr/articles/m1-foundations-and-notation.md),
give the formulas and assumptions behind the target, the stacked fit and
the calibration. The companion methods preprint (Lee, Williams and
Savitsky 2026,
[doi:10.5281/zenodo.22407935](https://doi.org/10.5281/zenodo.22407935) )
describes the stacked fit and its calibration.

## See also

Design:
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md),
[`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md),
[`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md).
Target:
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
Fit:
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md),
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).
Compare:
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md).
Read:
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
Returned objects and their fields:
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]

Authors:

- JoonHo Lee <jlee296@ua.edu>
  ([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]
