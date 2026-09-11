#' Settings for the fitting functions
#'
#' `pv_control()` collects the settings used by [pv_fit()] and by
#' [pv_fit_direct()], [pv_fit_reference()] and [pv_fit_stack_psis()]: the
#' method, the sampler settings and engine passed to the model-fitting
#' function, the interval level, the Pareto k-hat cut-off of `"stack_psis"`,
#' and what a fitted object keeps. It checks each value and returns them in
#' one object.
#'
#' @details
#' ## Which settings each method uses
#'
#' `chains`, `iter`, `warmup`, `cores`, `seed` and `backend` are passed to the
#' function that fits the model: the bundled brms engine or your own
#' `fit_function`. They have no effect when a method starts from draws
#' computed elsewhere. For `"stack_direct"`, pvstackr also checks that the
#' sampler diagnostics report `chains` chains with `iter - warmup` draws each
#' after warmup, and blocks the fit (status `"blocked"`) if they do not.
#'
#' `center` is used only by `"stack_direct"` and `psis_k_threshold` only by
#' `"stack_psis"`. `conf_level` and the settings that decide what a fitted
#' object keeps apply to all three methods.
#'
#' ## The bundled brms engine
#'
#' With `backend = "brms"` and no `fit_function`, `"stack_direct"` fits the
#' stacked model with the brms function `brm()`, through
#' [pv_backend_brms_fit_function()]. This needs the brms and posterior
#' packages. If cmdstanr is installed, CmdStan must be configured, otherwise
#' the fit stops with an error; if cmdstanr is not installed, brms uses rstan.
#' The other `backend` values select no engine: they are passed to your
#' `fit_function` as its `backend` argument and recorded in the fit.
#'
#' ## What a fitted object keeps
#'
#' With `return_draws = TRUE` (the default), a fit keeps its fixed-effect
#' draws: for `"stack_direct"` the calibrated draws, read with [get_draws()];
#' for `"per_pv"` the draws of each plausible value; for `"stack_psis"` the
#' stacked draws with the normalized weights for each plausible value. The
#' last two are in [get_diagnostics()], and [get_draws()] returns `NULL` for
#' them. Draws of other parameters, such as `sigma`, are not kept.
#'
#' `keep_data`, `keep_backend_fit` and `keep_log_lik` (all `FALSE` by
#' default) keep the data, the object returned by the model-fitting function
#' and the log-likelihood draws; each makes a saved fit larger. A blocked fit
#' keeps none of these parts: its `control` records all four settings as
#' `FALSE`, whatever was requested. [pvstackr_object_contracts] describes the
#' parts of a fitted object.
#'
#' @param method The method these settings are for: `"stack_direct"`
#'   (default), `"per_pv"` or `"stack_psis"`. It must match the method you fit
#'   with: [pv_fit()] stops with an error if it differs from its `method`
#'   argument, and each method function accepts only its own method.
#' @param chains Number of Markov chains, a whole number of at least 1.
#'   Default `4L`.
#' @param iter Number of iterations per chain, warmup included, a whole number
#'   of at least 2. Default `2000L`.
#' @param warmup Number of warmup iterations per chain, a whole number from 0
#'   to `iter - 1`. `NULL` (default) uses `floor(iter / 2)`, which is 1000 for
#'   the default `iter`.
#' @param cores Number of cores the model-fitting function may use, a whole
#'   number of at least 1. Default `1L`.
#' @param seed The random seed, a whole number of at least 0, or `NULL`
#'   (default) to leave the seed to the model-fitting function.
#' @param backend The engine for `"stack_direct"`: `"none"` (default),
#'   `"injected"`, `"brms"` or `"cmdstanr"`. Only `"brms"` selects one: when
#'   you give no `fit_function`, `"stack_direct"` then uses the bundled brms
#'   engine (see Details). With any other value, `"stack_direct"` needs your
#'   own `fit_function`, `draws_function` and `diagnose_function` (see
#'   [pv_fit()]). `"per_pv"` and `"stack_psis"` have no bundled engine and
#'   always need your own functions or draws computed elsewhere, whatever the
#'   value.
#' @param conf_level The level of the intervals in the estimate table, a
#'   number strictly between 0 and 1. Default `0.95`. It applies to all three
#'   methods; the `conf_level` stored in a [pv_brr_target()] object does not
#'   change the reported intervals.
#' @param psis_k_threshold The cut-off for the Pareto k-hat values (Pareto
#'   shape estimates) that come with the importance weights of
#'   `"stack_psis"`, a number greater than 0 and at most 0.7. Default `0.7`.
#'   pvstackr blocks a `"stack_psis"` fit unless the k-hat of every plausible
#'   value is finite and below this cut-off, so a smaller value makes the check
#'   stricter; values above 0.7 are not allowed.
#' @param center Where the calibrated fixed-effect draws are centered:
#'   `"target"` (default) or `"posterior"`. With `"target"`, the Cholesky
#'   calibration correction (CCC) gives the draws the mean and covariance of
#'   the target, whose estimates are the ones reported. `"stack_direct"`
#'   requires this value. `"posterior"` would keep the mean of the stacked
#'   draws and calibrate only their covariance; it is accepted here, but no
#'   fitting function uses it: [pv_fit_direct()] stops with an error for it.
#'   `"per_pv"` and `"stack_psis"` ignore `center`.
#' @param allow_target_nearpd Whether a target covariance matrix that is not
#'   positive definite may be repaired. pvstackr has no such repair and stops
#'   with an error on such a target, so the value must be `FALSE` (default);
#'   `TRUE` stops with an error.
#' @param return_draws Whether a fit keeps its fixed-effect draws (see
#'   Details). Default `TRUE`.
#' @param keep_data Whether a fit keeps the data frame. Default `FALSE`: the
#'   fit then keeps only a description of the data, such as column names and
#'   checksums.
#' @param keep_backend_fit Whether a fit keeps the object returned by the
#'   model-fitting function, such as a brms fit (for `"per_pv"`, one per
#'   plausible value). This object may contain the data, so the fitting
#'   functions stop with an error when `keep_backend_fit = TRUE` and
#'   `keep_data = FALSE`. Default `FALSE`.
#' @param keep_log_lik Whether a fit keeps the log-likelihood draws. Only
#'   `"stack_direct"` can extract them: [pv_fit_direct()] does so when given
#'   `extract_log_lik = TRUE` and a `log_lik_function`. Set `TRUE` only then:
#'   otherwise a `"stack_direct"` fit, or a `"stack_psis"` fit from a
#'   `fit_function`, stops with an error. Default `FALSE`.
#' @param verbose Stored in the object but not used by any pvstackr function.
#'   Default `FALSE`. For progress messages while the target is computed, use
#'   the `verbose` argument of [pv_brr_target()].
#'
#' @returns A `pvstackr_control` object: a list of the 16 settings in the order
#'   of the arguments, with `warmup` filled in and whole numbers stored as
#'   integers. Pass it as `control` to [pv_fit()] or to a method function. To
#'   change a setting, call `pv_control()` again rather than editing the list:
#'   the fitting functions stop with an error if a value was replaced by one of
#'   another type, as `ctrl$chains <- 2` does (it stores a double). `print()`
#'   shows `method`, `backend`, `iter`, `warmup` and `chains`, notes that
#'   target repair is not supported, and returns the object invisibly.
#' @examples
#' # Default settings, for method = "stack_direct".
#' ctrl <- pv_control()
#' ctrl
#'
#' # A control's method must match the method you fit with.
#' ctrl_psis <- pv_control(method = "stack_psis", psis_k_threshold = 0.7)
#' ctrl_psis$method
#'
#' # Settings that make a fit keep more. keep_backend_fit = TRUE needs
#' # keep_data = TRUE, and keep_log_lik = TRUE needs extract_log_lik = TRUE
#' # and a log_lik_function in pv_fit_direct().
#' ctrl_heavy <- pv_control(
#'   keep_data = TRUE, keep_backend_fit = TRUE, keep_log_lik = TRUE
#' )
#' c(ctrl_heavy$keep_data, ctrl_heavy$keep_backend_fit, ctrl_heavy$keep_log_lik)
#' @family pvstackr-fitting
#' @seealso [pv_backend_brms_fit_function()] for the bundled brms engine;
#'   [pvstackr_object_contracts] for the parts of a fitted object.
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
  psis_k_threshold <- pv_validate_psis_k_threshold(psis_k_threshold)
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

pv_fit_blocked_control <- function(control) {
  out <- pv_validate_control(control)
  out$return_draws <- FALSE
  out$keep_data <- FALSE
  out$keep_backend_fit <- FALSE
  out$keep_log_lik <- FALSE
  pv_validate_control(out)
}

#' @rdname pv_control
#' @param x A `pvstackr_control` object from `pv_control()`.
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
