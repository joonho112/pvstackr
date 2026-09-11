#' pvstackr: Bayesian Plausible-Value Analysis with One Calibrated Stacked Fit
#'
#' pvstackr fits survey-weighted linear regression models to assessment data
#' whose outcome is given as plausible values, as in PISA, and reports
#' estimates, standard errors and intervals for the fixed effects. The default
#' method, `stack_direct`, fits one Bayesian model to the stacked data (one
#' copy of the data per plausible value) and then applies the Cholesky
#' calibration correction (CCC): the fixed-effect draws are transformed to
#' have the mean and covariance of a design-based target, whose estimates,
#' standard errors and degrees of freedom are the ones reported. The other
#' methods, `per_pv` (one fit per plausible value) and `stack_psis` (one
#' stacked fit reweighted toward each plausible value), work from draws,
#' fitting functions or importance weights that you supply.
#'
#' @details
#' ## Workflow
#'
#' 1. Optionally, [pv_design()] checks and records the plausible-value,
#'    final-weight and replicate-weight columns and the Fay coefficient. By
#'    default it finds PISA-style column names with [detect_pisa_pv_columns()]
#'    and [detect_pisa_brr_replicate_weights()]. This step can be skipped:
#'    [pv_brr_target()] takes the same column arguments and finds the
#'    PISA-style columns itself.
#' 2. [pv_brr_target()] computes the target: a survey-weighted regression for
#'    each plausible value, the sampling covariance of the coefficients from
#'    the BRR-Fay replicate weights, and Rubin's rules across plausible values.
#' 3. [pv_fit()] fits a method with settings from [pv_control()]. Only
#'    `stack_direct` has a bundled engine, brms, selected with
#'    `pv_control(backend = "brms")`.
#' 4. Optionally, [pv_compare_methods()] compares fits made with different
#'    methods.
#' 5. [get_estimates()], [get_target()], [get_draws()] and [get_diagnostics()]
#'    read the results. A fit with `status` `"blocked"` has no estimates;
#'    `reason_codes` says why.
#'
#' ## Scope
#'
#' pvstackr is built for the fixed effects of a survey-weighted linear model.
#' [pv_brr_target()], `stack_direct` and `stack_psis` stop with an error on a
#' random-effect term such as `(1 | school)`. Only the fixed-effect
#' coefficients (intercept and slopes) are reported; other parameters, such as
#' the residual standard deviation `sigma`, are neither calibrated nor
#' reported.
#'
#' By pvstackr's reporting rule, only a `stack_direct` fit whose target uses
#' Barnard-Rubin degrees of freedom (`df_method = "barnard_rubin"` in
#' [pv_brr_target()]) has coverage-claimable intervals, which can be read as
#' confidence intervals with nominal coverage. All other intervals are
#' descriptive, including those of a `stack_direct` fit whose target uses the
#' default classic Rubin degrees of freedom. [pv_fit()] gives the full rule.
#'
#' ## Articles
#'
#' The Applied articles, starting with
#' `vignette("a1-getting-started", package = "pvstackr")`, show how to run an
#' analysis and read the results. The Method articles, starting with
#' `vignette("m1-foundations-and-notation", package = "pvstackr")`, give the
#' formulas and assumptions behind the target, the stacked fit and the
#' calibration. The companion methods preprint (Lee, Williams and Savitsky
#' 2026, \doi{10.5281/zenodo.22407935}) describes the stacked fit and its
#' calibration.
#'
#' @seealso
#' Design: [pv_design()], [detect_pisa_pv_columns()],
#'   [detect_pisa_brr_replicate_weights()]. Target: [pv_brr_target()]. Fit:
#'   [pv_fit()], [pv_fit_direct()], [pv_fit_reference()], [pv_fit_stack_psis()],
#'   [pv_control()]. Compare: [pv_compare_methods()]. Read: [get_estimates()],
#'   [get_target()], [get_draws()], [get_diagnostics()]. Returned objects and
#'   their fields: [pvstackr_object_contracts].
#'
#' @keywords internal
"_PACKAGE"
