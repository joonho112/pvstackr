#' Construct pvstackr Fitting Controls
#'
#' `pv_control()` validates package-level options used by the fitting,
#' calibration, diagnostics, and object-retention layers, returning a frozen
#' `pvstackr_control` object that every `pv_fit*()` entry point consumes. It
#' centralizes backend policy, interval level, the calibration centering
#' convention, and how heavy a fitted object is allowed to be.
#'
#' @details
#' ## Centering and reportable output
#' The `center` convention decides whether a fit is *reportable* or merely
#' *diagnostic*. Reportable `stack_direct` output requires `center = "target"`:
#' fixed-effect draws are CCC-calibrated so their mean and covariance match the
#' external Rubin/BRR-Fay target, and the estimate table inherits the target's
#' standard errors, degrees of freedom, and interval metadata. `center =
#' "posterior"` is diagnostic and exploratory only: it leaves fixed-effect draws
#' at the raw stacked posterior mean while still computing target-covariance
#' diagnostics against the external target. [pv_fit_direct()] rejects a control
#' with `center = "posterior"` rather than emit a non-reportable `stack_direct`
#' estimate table; treat `"posterior"` as a CCC check, never a deliverable.
#'
#' ## Object retention
#' The retention flags govern how much a fitted `pvstackr_fit` carries, and the
#' defaults keep fits light. `return_draws` (default `TRUE`) retains the
#' reportable fixed-effect draws used downstream by [get_draws()]; the remaining
#' flags default to `FALSE`. `keep_data` retains the user data frame,
#' `keep_backend_fit` retains the heavy backend fit object, and `keep_log_lik`
#' retains log-likelihood draws. Enable the `FALSE`-by-default flags only when
#' you need the extra payload (for example, re-extraction or model checking), as
#' each materially increases the size of the saved object.
#'
#' @param method Public method identifier. Character scalar; must be one of
#'   `"stack_direct"`, `"stack_psis"`, or `"per_pv"`. Default `"stack_direct"`.
#'   A control's `method` must equal the `method` passed to [pv_fit()] and to the
#'   dispatched fitter.
#' @param chains Number of MCMC chains requested by a live backend. Integer-valued
#'   scalar, `>= 1`.
#' @param iter Total iterations per chain. Integer-valued scalar, `>= 2`.
#' @param warmup Warmup iterations per chain. Integer-valued scalar, `>= 0` and
#'   strictly less than `iter`. If `NULL`, defaults to `floor(iter / 2)`.
#' @param cores Number of cores requested by a live backend. Integer-valued
#'   scalar, `>= 1`.
#' @param seed Optional random seed. Integer-valued scalar `>= 0`, or `NULL` for
#'   no fixed seed.
#' @param backend Backend policy. Character scalar; one of `"none"`,
#'   `"injected"`, `"brms"`, or `"cmdstanr"`. Default `"none"`. In this package
#'   stage the backend fit is injected or precomputed.
#' @param conf_level Confidence or credible-interval level for report tables.
#'   Numeric scalar in `(0, 1)`. Default `0.95`.
#' @param psis_k_threshold Pareto-k warning threshold for `stack_psis`. Numeric
#'   scalar in `(0, 1]`. Default `0.7`, matching the Vehtari et al. unreliable-k
#'   cutoff.
#' @param center Calibration centering convention, `"target"` or `"posterior"`.
#'   Reportable `stack_direct` output requires `"target"`. `"posterior"` is
#'   reserved for CCC diagnostic/exploratory checks: it leaves fixed-effect draws
#'   at the raw stacked posterior mean while still computing target-covariance
#'   diagnostics against the external target. Default `"target"`.
#' @param allow_target_nearpd Reserved for future target-covariance repair.
#'   Logical scalar; must be `FALSE`. Automatic target repair is not currently
#'   supported. Default `FALSE`.
#' @param return_draws Logical scalar. Whether fitted objects retain reportable
#'   fixed-effect draws (read via [get_draws()]). Default `TRUE`.
#' @param keep_data Logical scalar. Whether fitted objects may retain the user
#'   data frame. Default `FALSE` (fits stay light).
#' @param keep_backend_fit Logical scalar. Whether fitted objects may retain the
#'   heavy backend fit object. Default `FALSE` (fits stay light).
#' @param keep_log_lik Logical scalar. Whether fitted objects may retain
#'   log-likelihood draws. Default `FALSE` (fits stay light).
#' @param verbose Logical scalar. Whether functions emit progress messages.
#'   Default `FALSE`.
#'
#' @returns A validated `pvstackr_control` object: a named list of the resolved
#'   options above, with class `c("pvstackr_control", "list")`. Pass it to a
#'   `pv_fit*()` function via the `control` argument; print it for a compact
#'   summary.
#' @examples
#' # Default controls target reportable stack_direct output.
#' ctrl <- pv_control()
#' ctrl
#'
#' # A control's method must match the method you fit with.
#' ctrl_psis <- pv_control(method = "stack_psis", psis_k_threshold = 0.7)
#' ctrl_psis$method
#'
#' # Retain the heavy backend fit and log-likelihood draws when you need them.
#' ctrl_heavy <- pv_control(keep_backend_fit = TRUE, keep_log_lik = TRUE)
#' c(ctrl_heavy$keep_backend_fit, ctrl_heavy$keep_log_lik)
#' @family pvstackr-fitting
#' @seealso [pv_fit()], [pv_fit_direct()], [pv_fit_reference()],
#'   [pv_fit_stack_psis()]
#' @export
pv_control <- function(
  method = "stack_direct",
  chains = 4L,
  iter = 2000L,
  warmup = NULL,
  cores = 1L,
  seed = NULL,
  backend = "none",
  conf_level = 0.95,
  psis_k_threshold = 0.7,
  center = "target",
  allow_target_nearpd = FALSE,
  return_draws = TRUE,
  keep_data = FALSE,
  keep_backend_fit = FALSE,
  keep_log_lik = FALSE,
  verbose = FALSE
) {
  method <- pv_validate_method(method)
  chains <- pv_assert_scalar_number(chains, "chains", integer = TRUE, lower = 1)
  iter <- pv_assert_scalar_number(iter, "iter", integer = TRUE, lower = 2)

  if (is.null(warmup)) {
    warmup <- floor(iter / 2)
  }
  warmup <- pv_assert_scalar_number(warmup, "warmup", integer = TRUE, lower = 0)
  if (warmup >= iter) {
    pv_abort("`warmup` must be smaller than `iter`.")
  }

  cores <- pv_assert_scalar_number(cores, "cores", integer = TRUE, lower = 1)
  seed <- pv_assert_scalar_number(seed, "seed", integer = TRUE, lower = 0, allow_null = TRUE)
  backend <- pv_validate_backend(backend)
  conf_level <- pv_assert_probability(conf_level, "conf_level")
  psis_k_threshold <- pv_assert_scalar_number(
    psis_k_threshold,
    "psis_k_threshold",
    lower = 0,
    upper = 1,
    inclusive_lower = FALSE,
    inclusive_upper = TRUE
  )
  center <- pv_validate_center(center)
  allow_target_nearpd <- pv_validate_target_repair_control(allow_target_nearpd)
  return_draws <- pv_assert_scalar_logical(return_draws, "return_draws")
  keep_data <- pv_assert_scalar_logical(keep_data, "keep_data")
  keep_backend_fit <- pv_assert_scalar_logical(keep_backend_fit, "keep_backend_fit")
  keep_log_lik <- pv_assert_scalar_logical(keep_log_lik, "keep_log_lik")
  verbose <- pv_assert_scalar_logical(verbose, "verbose")

  control <- list(
    method = method,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    backend = backend,
    conf_level = conf_level,
    psis_k_threshold = psis_k_threshold,
    center = center,
    allow_target_nearpd = allow_target_nearpd,
    return_draws = return_draws,
    keep_data = keep_data,
    keep_backend_fit = keep_backend_fit,
    keep_log_lik = keep_log_lik,
    verbose = verbose
  )
  class(control) <- c("pvstackr_control", "list")
  pv_validate_control(control)
}

#' @rdname pv_control
#' @param x A `pvstackr_control` object.
#' @param ... Ignored.
#' @export
print.pvstackr_control <- function(x, ...) {
  x <- pv_validate_control(x)
  cat("pvstackr control\n")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  backend: ", x$backend, "\n", sep = "")
  cat("  iter/warmup/chains: ", x$iter, "/", x$warmup, "/", x$chains, "\n", sep = "")
  cat("  target repair: unsupported (disabled)\n", sep = "")
  invisible(x)
}
