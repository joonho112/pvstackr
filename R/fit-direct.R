pv_fit_direct_estimates <- function(ccc, target, conf_level) {
  fe_names <- target$fe_names
  estimate <- ccc$psi_hat[fe_names]
  se <- sqrt(diag(ccc$Sigma_target))[fe_names]
  df <- target$df[fe_names]
  if (!is.numeric(estimate) || any(!is.finite(estimate))) {
    pv_abort("Direct-fit calibrated estimates must be finite.")
  }
  if (!is.numeric(se) || any(!is.finite(se)) || any(se <= 0)) {
    pv_abort("Direct-fit target standard errors must be finite and positive.")
  }
  if (!is.numeric(df) || any(is.na(df)) || any(df <= 0)) {
    pv_abort("Direct-fit target degrees of freedom must be positive.")
  }
  df_complete <- target$df_complete
  if (is.null(df_complete)) {
    df_complete <- stats::setNames(rep(NA_real_, length(fe_names)), fe_names)
  } else {
    df_complete <- df_complete[fe_names]
  }
  alpha <- 1 - conf_level
  crit <- stats::qt(1 - alpha / 2, df = df)
  data.frame(
    term = fe_names,
    estimate = unname(estimate),
    se = unname(se),
    std.error = unname(se),
    df = unname(df),
    df_method = target$df_method,
    df_complete = unname(df_complete),
    conf_level = conf_level,
    conf_low = unname(estimate - crit * se),
    conf_high = unname(estimate + crit * se),
    conf.low = unname(estimate - crit * se),
    conf.high = unname(estimate + crit * se),
    interval_role = target$interval_role,
    coverage_claim_allowed = target$coverage_claim_allowed,
    parameter_scope = "fixed_effect",
    target_source = target$target_source,
    target_hash = target$target_hash,
    stringsAsFactors = FALSE
  )
}

pv_fit_direct_status <- function(stack_fit, ccc) {
  reason_codes <- character()
  warnings <- character()

  prior_policy <- stack_fit$meta$prior_policy %||% list()
  if (isTRUE(prior_policy$explicit_prior_warning)) {
    prior_reason_code <- prior_policy$reason_code
    if (!is.character(prior_reason_code) ||
        length(prior_reason_code) != 1L ||
        is.na(prior_reason_code) ||
        !nzchar(prior_reason_code)) {
      prior_reason_code <- "explicit_prior_warning"
    }
    reason_codes <- c(reason_codes, prior_reason_code)
    warnings <- c(warnings, prior_policy$warning)
  }

  center_status <- ccc$diagnostics$center_status %||% "ok"
  if (identical(center_status, "warning")) {
    reason_codes <- c(reason_codes, ccc$diagnostics$center_reason_code %||% "center_separation_warning")
    warnings <- c(warnings, ccc$warnings)
  }
  if (identical(center_status, "blocked")) {
    reason_codes <- c(reason_codes, ccc$diagnostics$center_reason_code %||% "center_separation_blocked")
    warnings <- c(warnings, ccc$warnings)
    status <- "blocked"
  } else {
    status <- if (length(reason_codes) > 0L) "warning" else "ok"
  }

  conditioning_status <- ccc$diagnostics$conditioning_status %||% "ok"
  if (identical(conditioning_status, "warning")) {
    reason_codes <- c(reason_codes, ccc$diagnostics$conditioning_reason_code %||% "ccc_conditioning_warning")
    warnings <- c(warnings, ccc$warnings)
  }
  if (identical(conditioning_status, "blocked")) {
    reason_codes <- c(reason_codes, ccc$diagnostics$conditioning_reason_code %||% "ccc_conditioning_blocked")
    warnings <- c(warnings, ccc$warnings)
    status <- "blocked"
  } else if (!identical(status, "blocked") && length(reason_codes) > 0L) {
    status <- "warning"
  }

  reason_codes <- unique(reason_codes[nzchar(reason_codes)])
  warnings <- unique(warnings[nzchar(warnings)])
  if (identical(status, "warning") && length(warnings) == 0L) {
    warnings <- "Direct fit produced warning-level diagnostics."
  }
  list(
    status = status,
    reason_codes = reason_codes,
    warnings = warnings
  )
}

pv_fit_not_implemented <- function(method) {
  pv_abort(sprintf(
    "`method = \"%s\"` is recognized but not implemented in this package stage. Use `method = \"stack_direct\"`, `method = \"per_pv\"`, `method = \"stack_psis\"`, `pv_fit_direct()`, `pv_fit_reference()`, or `pv_fit_stack_psis()`.",
    method
  ))
}

#' Fit a pvstackr Method
#'
#' `pv_fit()` is the generic public fitting entry point. In this package stage,
#' `method = "stack_direct"`, `method = "per_pv"`, and
#' `method = "stack_psis"` are implemented.
#'
#' @details
#' ## Method dispatch
#' `pv_fit()` validates `method` and `control`, then forwards to a
#' method-specific fitter:
#' - `"stack_direct"` (default) dispatches to [pv_fit_direct()] and **requires**
#'   a `pvstackr_brr_target` `target` from [pv_brr_target()]; `pv_fit()` errors
#'   if `target` is `NULL`.
#' - `"per_pv"` dispatches to [pv_fit_reference()] (the orthodox per-PV
#'   reference); `target` is ignored.
#' - `"stack_psis"` dispatches to [pv_fit_stack_psis()] (the PSIS-reweighted
#'   stacked path); `target` is ignored.
#'
#' If `control` is `NULL` it is constructed with `pv_control(method = method)`;
#' otherwise `control$method` **must** equal `method` or the call errors. Any
#' arguments in `...` are forwarded verbatim to the dispatched fitter, so consult
#' that fitter's signature for the available extras (for example, the injected
#' backend adapters of [pv_fit_direct()], or the injected PSIS weights of
#' [pv_fit_stack_psis()]).
#'
#' In this package stage, only `stack_direct` output is coverage-claimable: its
#' intervals are backed by the external Rubin/BRR-Fay target. `per_pv` and
#' `stack_psis` intervals are descriptive/reference, even with Barnard-Rubin
#' degrees of freedom.
#'
#' @param data Analysis data frame.
#' @param formula Two-sided formula with `OUTCOME` on the left-hand side.
#' @param target A method-specific target object. For `method = "stack_direct"`,
#'   this must be a `pvstackr_brr_target` from [pv_brr_target()] and is required.
#'   Ignored for `"per_pv"` and `"stack_psis"`. Default `NULL`.
#' @param method Public method identifier. Character scalar; one of
#'   `"stack_direct"`, `"stack_psis"`, or `"per_pv"`. Default `"stack_direct"`.
#' @param control Optional [pv_control()] object. If `NULL`, one is built with
#'   `pv_control(method = method)`. If supplied, `control$method` must equal
#'   `method`. Default `NULL`.
#' @param ... Additional arguments forwarded to the dispatched fitter:
#'   [pv_fit_direct()] for `"stack_direct"`, [pv_fit_reference()] for `"per_pv"`,
#'   or [pv_fit_stack_psis()] for `"stack_psis"`.
#'
#' @returns A `pvstackr_fit` object. Read it with the accessors
#'   ([get_estimates()], [get_target()], [get_draws()], [get_diagnostics()])
#'   rather than by `$`-indexing; user-facing components include `estimates`,
#'   `draws`, `target`, `diagnostics`, `status`, `reason_codes`, `warnings`, and
#'   `method`.
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#'
#' # Build the design and the external BRR-Fay target on the fixture.
#' design <- pv_design(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
#' )
#' target <- pv_brr_target(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_cols = design$pv_cols, weight_col = design$weight_col,
#'   rep_weight_cols = design$rep_weight_cols, fay_k = design$fay_k,
#'   id_cols = design$id_cols
#' )
#'
#' # A live stack_direct fit needs an injected/precomputed backend:
#' \dontrun{
#' fit <- pv_fit(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female,
#'   target = target, method = "stack_direct",
#'   fit_function = my_backend_adapter, draws_function = my_draws_adapter
#' )
#' }
#'
#' # Inspect the reportable object surface via the bundled cached fit instead.
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit     # a stack_direct pvstackr_fit
#'   fit                          # compact console print
#'   head(get_estimates(fit))     # reportable fixed-effect table
#' }
#' @family pvstackr-fitting
#' @seealso [pv_fit_direct()], [pv_fit_reference()], [pv_fit_stack_psis()],
#'   [pv_control()]; [pv_brr_target()], [pv_compare_methods()], [get_estimates()]
#' @export
pv_fit <- function(
  data,
  formula,
  target = NULL,
  method = "stack_direct",
  control = NULL,
  ...
) {
  method <- pv_validate_method(method)
  control <- if (is.null(control)) pv_control(method = method) else pv_validate_control(control)
  if (!identical(control$method, method)) {
    pv_abort("`control$method` must match `method`.")
  }
  if (identical(method, "per_pv")) {
    return(pv_fit_reference(
      data = data,
      formula = formula,
      control = control,
      ...
    ))
  }
  if (identical(method, "stack_psis")) {
    return(pv_fit_stack_psis(
      data = data,
      formula = formula,
      control = control,
      ...
    ))
  }
  if (!identical(method, "stack_direct")) {
    pv_fit_not_implemented(method)
  }
  if (is.null(target)) {
    pv_abort("`target` is required for `method = \"stack_direct\"`.")
  }
  pv_fit_direct(
    data = data,
    formula = formula,
    target = target,
    control = control,
    ...
  )
}

#' Fit the Direct Stacked Plausible-Value Model
#'
#' `pv_fit_direct()` runs the v0.1 `stack_direct` path: fixed-effect-only
#' compatibility preflight, one stacked plausible-value fit, CCC calibration to
#' an external BRR-Fay/Rubin target, and assembly of a reportable
#' `pvstackr_fit` object.
#'
#' @details
#' In v0.1, `stack_direct` is fixed-effect-only and requires an external
#' `pvstackr_brr_target` from [pv_brr_target()]. The formula RHS and derived
#' fixed-effect names must match the target exactly. Group terms such as
#' `(1 | school)` are rejected until a two-level target engine is implemented.
#'
#' The returned `pvstackr_fit` has `status = "ok"`, `"warning"`, or
#' `"blocked"`. Yellow CCC center separation and explicit priors produce
#' warning-status fits; red CCC center separation blocks reportable estimates.
#' CCC calibration-matrix conditioning is also gated by `kappa_A`: values at or
#' above `1e6` produce warning-status fits, and values at or above `1e8` block
#' reportable estimates. Reportable estimates and retained draws are
#' fixed-effect-only.
#'
#' With the default `control$center = "target"`, the reportable fixed-effect
#' estimate table is target-based: `estimate` is the external Rubin/BRR-Fay
#' coefficient after CCC centering, `se` is `sqrt(diag(target$T_MI))`, and the
#' degrees of freedom and interval metadata are inherited from `target`.
#' Therefore `df_method`, `df_complete`, `interval_role`, and
#' `coverage_claim_allowed` describe the external target policy, not a residual
#' df or interval policy estimated by the stacked backend fit. When the target
#' uses classic Rubin df, the reported df are the target's classic Rubin
#' imputation df; when it uses Barnard-Rubin df, `df_complete` is the explicit
#' complete-data df supplied to [pv_brr_target()]. The stacked fit supplies the
#' draw cloud used for calibration, retained calibrated draws when requested,
#' and center-separation agreement diagnostics; it is not a replacement source
#' for the headline Rubin/BRR-Fay fixed-effect numbers.
#'
#' Reportable `pv_fit_direct()` output requires `control$center = "target"`.
#' The lower-level CCC convention `center = "posterior"` is diagnostic and
#' exploratory: it leaves the output fixed-effect center at the raw stacked
#' posterior mean while pairing that center with the external target covariance.
#' `pv_fit_direct()` does not expose that hybrid as a reportable
#' `stack_direct` estimate table.
#'
#' See also [pvstackr_object_contracts].
#'
#' @param data Analysis data frame.
#' @param formula Two-sided formula with `OUTCOME` on the left-hand side.
#' @param target A `pvstackr_brr_target` object from [pv_brr_target()].
#' @param control A [pv_control()] object with `method = "stack_direct"`.
#' @param family Optional backend family object passed to the injected fit
#'   function.
#' @param prior Optional backend prior object. Explicit priors are allowed but
#'   reported as warning-level diagnostics because the current identity result
#'   is scoped to MLE/flat-prior fixed-effect regimes.
#' @param fit_function Injected backend fitting function for this package stage.
#' @param draws_function Function extracting posterior draws from the injected
#'   backend fit.
#' @param param_map Optional explicit draw-column map passed to the stack-fit
#'   layer. Supply `fe_names` or `fe_idx` to identify fixed-effect columns, and
#'   optional `vc_names` or `vc_idx` for nuisance variance-component columns.
#'   Use `vc_names = character()` to drop all nuisance columns. Explicit maps
#'   are recommended when backend draw names do not follow the automatic `b_*`
#'   fixed-effect convention, including distributional names such as
#'   `b_sigma_*`.
#' @param diagnose_function Optional backend diagnostic extractor.
#' @param log_lik_function Optional log-likelihood extractor.
#' @param extract_log_lik Whether to extract log-likelihood draws.
#' @param cache_dir,cache_stem Cache location metadata passed to the injected
#'   fit function.
#' @param additional_args Additional named arguments passed to the injected fit
#'   function.
#'
#' @returns A `pvstackr_fit` object.
#' @section Reportable scope and coverage:
#' In this package stage, reportable output is **fixed-effect-only**; variance
#' components are fit but not calibrated to the target. Coverage claims are
#' enabled **only** for `stack_direct` rows backed by the external Rubin/BRR-Fay
#' target with Barnard-Rubin df and an explicit `df_complete`
#' (`interval_role = "coverage_barnard_rubin"`, `coverage_claim_allowed = TRUE`);
#' with classic Rubin df the rows are descriptive
#' (`interval_role = "descriptive_classic_rubin"`, `coverage_claim_allowed = FALSE`),
#' as in the bundled cached example. `per_pv` and `stack_psis` intervals are
#' **descriptive/reference** even with Barnard-Rubin degrees of freedom. "One
#' stacked fit" describes the computational topology, not a benchmarked speed
#' claim.
#'
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#'
#' # Build the design and the external BRR-Fay target on the fixture.
#' design <- pv_design(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
#' )
#' target <- pv_brr_target(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_cols = design$pv_cols, weight_col = design$weight_col,
#'   rep_weight_cols = design$rep_weight_cols, fay_k = design$fay_k,
#'   id_cols = design$id_cols
#' )
#'
#' # A live stack_direct fit needs injected/precomputed backend adapters
#' # (e.g. brms/cmdstanr); not run in checks.
#' \dontrun{
#' fit <- pv_fit_direct(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female,
#'   target = target,
#'   fit_function = my_backend_adapter, draws_function = my_draws_adapter
#' )
#' }
#'
#' # Inspect the reportable fixed-effect surface via the bundled cached fit.
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit     # a stack_direct pvstackr_fit
#'   head(get_estimates(fit))     # reportable fixed-effect table
#' }
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
#'
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#' @family pvstackr-fitting
#' @seealso [pv_fit()], [pv_fit_reference()], [pv_fit_stack_psis()],
#'   [pv_control()]; [pv_brr_target()], [pv_compare_methods()], [get_estimates()],
#'   [pvstackr_object_contracts].
#' @export
pv_fit_direct <- function(
  data,
  formula,
  target,
  control = pv_control(method = "stack_direct"),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  param_map = NULL,
  diagnose_function = NULL,
  log_lik_function = NULL,
  extract_log_lik = FALSE,
  cache_dir = "cache",
  cache_stem = "pvstackr-stack-direct",
  additional_args = list()
) {
  control <- if (is.null(control)) pv_control(method = "stack_direct") else pv_validate_control(control)
  if (!identical(control$method, "stack_direct")) {
    pv_abort("`control$method` must be `stack_direct` for `pv_fit_direct()`.")
  }
  if (!identical(control$center, "target")) {
    pv_abort("Reportable `pv_fit_direct()` output requires `control$center = \"target\"`; `center = \"posterior\"` is reserved for CCC diagnostic/exploratory checks.")
  }
  preflight <- pv_stack_direct_preflight(
    data = data,
    formula = formula,
    target = target
  )
  design <- new_pvstackr_design(
    data = data,
    formula = formula,
    pv_cols = target$pv_cols,
    weight_col = target$weight_col,
    rep_weight_cols = target$rep_weight_cols,
    fay_k = target$fay_k,
    id_cols = target$id_cols,
    roles = list(
      outcome_placeholder = "OUTCOME",
      method = "stack_direct"
    ),
    provenance = list(
      source = "pv_fit_direct",
      target_hash = target$target_hash
    )
  )
  stack_fit <- pv_stack_fit(
    data = data,
    formula = formula,
    pv_cols = target$pv_cols,
    weight_col = target$weight_col,
    control = control,
    family = family,
    prior = prior,
    fit_function = fit_function,
    draws_function = draws_function,
    param_map = param_map,
    diagnose_function = diagnose_function,
    log_lik_function = log_lik_function,
    extract_log_lik = extract_log_lik,
    cache_dir = cache_dir,
    cache_stem = cache_stem,
    additional_args = additional_args
  )
  ccc <- ccc_calibrate(
    draws = stack_fit$stacked_draws,
    target = target,
    param_map = stack_fit$param_map,
    center = control$center
  )
  status_info <- pv_fit_direct_status(stack_fit, ccc)
  estimates <- if (identical(status_info$status, "blocked")) {
    data.frame()
  } else {
    pv_fit_direct_estimates(ccc, target, control$conf_level)
  }
  draws <- if (isTRUE(control$return_draws) && !identical(status_info$status, "blocked")) {
    ccc$draws_fe_cal
  } else {
    NULL
  }

  new_pvstackr_fit(
    method = "stack_direct",
    design = design,
    target = target,
    stack_fit = stack_fit,
    ccc = ccc,
    estimates = estimates,
    draws = draws,
    diagnostics = list(
      preflight = preflight,
      stack_fit = stack_fit$diagnostics,
      stack_fit_warnings = stack_fit$warnings,
      ccc = ccc$diagnostics
    ),
    status = status_info$status,
    control = control,
    reason_codes = status_info$reason_codes,
    provenance = list(
      wrapper_function = "pv_fit_direct",
      target_hash = target$target_hash,
      ccc_target_hash = ccc$target_hash,
      stack_fit_long_data_hash = stack_fit$weight_summary$long_data_hash
    ),
    warnings = status_info$warnings
  )
}
