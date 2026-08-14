#' pvstackr: Stacked-Fit Calibration to Rubin/BRR-Fay Fixed-Effect Targets for Plausible Values
#'
#' `pvstackr` currently implements the `stack_direct` path: one stacked
#' plausible-value fit plus fixed-effect CCC calibration to a design-based
#' external Rubin/BRR-Fay target. It also implements the `per_pv` reference
#' path: one backend fit or posterior-draw source per plausible value plus
#' model-based Rubin pooling, and the `stack_psis` diagnostic/reference path:
#' one stacked draw source plus caller-declared external PSIS weights, model-based
#' Rubin pooling, and Pareto-k gating.
#'
#' @details
#' ## Typical workflow
#'
#' A reportable analysis moves through the package API in five stages: declare
#' the plausible-value design with [pv_design()]; assemble the external
#' Rubin/BRR-Fay fixed-effect target with [pv_brr_target()]; fit a method with
#' [pv_fit()] (default `method = "stack_direct"`); optionally line the methods
#' up side by side with [pv_compare_methods()]; and read the results with the
#' accessor family [get_estimates()], [get_target()], [get_draws()], and
#' [get_diagnostics()].
#'
#' ## Where to learn more
#'
#' The documentation is organized as two vignette tracks. The Applied track is
#' workflow-first and starts at
#' `vignette("a1-getting-started", package = "pvstackr")`. The Method track is
#' the rigorous treatment of the notation, target, stacked bridge, and
#' calibration, and starts at
#' `vignette("m1-foundations-and-notation", package = "pvstackr")`.
#'
#' ## Scope and honesty
#'
#' Reportable output is fixed-effect-only (variance components are fit but not
#' calibrated); among the three methods only `stack_direct`, whose intervals
#' rest on the external Rubin/BRR-Fay target, is coverage-claimable, while
#' `per_pv` and `stack_psis` are descriptive/reference paths.
#'
#' @seealso
#' Design: [pv_design()], [detect_pisa_pv_columns()],
#'   [detect_pisa_brr_replicate_weights()]. Target: [pv_brr_target()]. Fit:
#'   [pv_fit()], [pv_fit_direct()], [pv_fit_reference()], [pv_fit_stack_psis()],
#'   [pv_control()]. Compare: [pv_compare_methods()]. Read: [get_estimates()],
#'   [get_target()], [get_draws()], [get_diagnostics()]. Contracts:
#'   [pvstackr_object_contracts].
#'
#' @keywords internal
"_PACKAGE"
