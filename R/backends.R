# Bundled backend adapters -----------------------------------------------------
#
# pv_fit(method = "stack_direct") accepts injected `fit_function` /
# `draws_function` pairs so that any Bayesian engine can drive the stacked fit.
# For the common case the package bundles a brms adapter, selected with
# `pv_control(backend = "brms")` and no `fit_function`. The bundled adapter is
# equivalent to the two-function adapter documented in the vignettes: it fits
# the prepared stacked formula with `brms::brm()` (Gaussian family unless one
# is supplied) and returns the fixed-effect and residual-scale draws as a plain
# base matrix.
#
# brms remains in Suggests: the adapter checks availability at call time and
# the package installs and runs without it.

pv_backend_package_version <- function(package) {
  tryCatch(
    as.character(utils::packageVersion(package)),
    error = function(e) NA_character_
  )
}

pv_backend_cmdstan_state <- function() {
  namespace_available <- requireNamespace("cmdstanr", quietly = TRUE)
  if (!namespace_available) {
    return(list(
      namespace_available = FALSE,
      package_version = NA_character_,
      cmdstan_configured = FALSE,
      cmdstan_version = NA_character_,
      cmdstan_path_basename = NA_character_,
      state_reason = "cmdstanr_namespace_unavailable",
      toolchain_checked = FALSE
    ))
  }

  path_error <- NULL
  path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) {
      path_error <<- conditionMessage(e)
      NA_character_
    }
  )
  version <- tryCatch(
    cmdstanr::cmdstan_version(error_on_NA = FALSE),
    error = function(e) NULL
  )
  version <- if (length(version) == 1L && !is.na(version)) {
    as.character(version)
  } else {
    NA_character_
  }
  configured <- length(path) == 1L && !is.na(path) && nzchar(path) &&
    dir.exists(path) && !is.na(version) && nzchar(version)

  list(
    namespace_available = TRUE,
    package_version = pv_backend_package_version("cmdstanr"),
    cmdstan_configured = configured,
    cmdstan_version = version,
    cmdstan_path_basename = if (configured) basename(path) else NA_character_,
    state_reason = if (configured) {
      "cmdstanr_namespace_and_cmdstan_configured"
    } else if (!is.null(path_error)) {
      "cmdstanr_namespace_without_configured_cmdstan"
    } else {
      "cmdstanr_namespace_with_invalid_cmdstan_state"
    },
    toolchain_checked = FALSE
  )
}

pv_backend_resolve_brms_engine <- function(cmdstan_state) {
  required <- c(
    "namespace_available", "package_version", "cmdstan_configured",
    "cmdstan_version", "cmdstan_path_basename", "state_reason",
    "toolchain_checked"
  )
  if (!is.list(cmdstan_state) || !all(required %in% names(cmdstan_state))) {
    pv_abort("Internal CmdStan state is incomplete.")
  }
  if (!is.logical(cmdstan_state$namespace_available) ||
      length(cmdstan_state$namespace_available) != 1L ||
      is.na(cmdstan_state$namespace_available) ||
      !is.logical(cmdstan_state$cmdstan_configured) ||
      length(cmdstan_state$cmdstan_configured) != 1L ||
      is.na(cmdstan_state$cmdstan_configured) ||
      !is.logical(cmdstan_state$toolchain_checked) ||
      length(cmdstan_state$toolchain_checked) != 1L ||
      is.na(cmdstan_state$toolchain_checked)) {
    pv_abort("Internal CmdStan state has invalid logical fields.")
  }
  if (isTRUE(cmdstan_state$cmdstan_configured) &&
      !isTRUE(cmdstan_state$namespace_available)) {
    pv_abort("Internal CmdStan state cannot be configured without the cmdstanr namespace.")
  }
  if (isTRUE(cmdstan_state$namespace_available) &&
      !isTRUE(cmdstan_state$cmdstan_configured)) {
    pv_abort(paste0(
      "The cmdstanr package is installed, but a working CmdStan installation ",
      "is not configured. Configure CmdStan with cmdstanr::set_cmdstan_path() ",
      "or remove cmdstanr to use brms's rstan backend. pvstackr does not ",
      "silently retry with rstan after selecting cmdstanr."
    ))
  }

  resolved_backend <- if (isTRUE(cmdstan_state$cmdstan_configured)) {
    "cmdstanr"
  } else {
    "rstan"
  }
  list(
    resolved_backend = resolved_backend,
    selection_policy = "cmdstanr_when_namespace_and_cmdstan_configured_else_rstan",
    selection_reason = if (identical(resolved_backend, "cmdstanr")) {
      "configured_cmdstan_selected"
    } else {
      "cmdstanr_namespace_absent_rstan_selected"
    }
  )
}

pv_backend_brms_engine_spec <- function(cmdstan_state = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    pv_abort("`pv_control(backend = \"brms\")` requires the brms package; install brms or inject a `fit_function`.")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    pv_abort("The bundled brms backend requires the posterior package before sampling can begin.")
  }
  if (is.null(cmdstan_state)) {
    cmdstan_state <- pv_backend_cmdstan_state()
  }
  resolution <- pv_backend_resolve_brms_engine(cmdstan_state)
  list(
    adapter_source = "bundled",
    adapter_id = "bundled_brms",
    requested_backend = "brms",
    resolved_backend = resolution$resolved_backend,
    engine_id = paste0("bundled_brms_", resolution$resolved_backend),
    selection_policy = resolution$selection_policy,
    selection_reason = resolution$selection_reason,
    package_versions = list(
      brms = pv_backend_package_version("brms"),
      posterior = pv_backend_package_version("posterior"),
      cmdstanr = cmdstan_state$package_version,
      cmdstan = cmdstan_state$cmdstan_version,
      # brms owns its rstan dependency; pvstackr does not directly query it.
      rstan = NA_character_
    ),
    cmdstan_state = cmdstan_state
  )
}

pv_backend_injected_engine_spec <- function(requested_backend) {
  requested_backend <- pv_assert_scalar_string(requested_backend, "requested_backend")
  list(
    adapter_source = "injected",
    adapter_id = "injected_fit_function",
    requested_backend = requested_backend,
    resolved_backend = "injected",
    engine_id = "injected_fit_function",
    selection_policy = "caller_supplied_fit_function",
    selection_reason = "fit_function_supplied",
    package_versions = list(),
    cmdstan_state = list(
      namespace_available = NA,
      cmdstan_configured = NA,
      toolchain_checked = FALSE,
      state_reason = "not_evaluated_for_injected_adapter"
    )
  )
}

#' Bundled brms adapter for the stacked fit
#'
#' The three functions that `pv_control(backend = "brms")` installs when no
#' adapter is injected: `pv_backend_brms_fit_function()` fits the prepared
#' stacked model, `pv_backend_brms_draws_function()` extracts its draws, and
#' `pv_backend_brms_sampler_diagnostics()` reports its sampler diagnostics.
#'
#' They are exported so that the bundled backend and the injected route are the
#' same code path rather than two implementations that have to be kept in step.
#' Passing all three to [pv_fit()] as `fit_function`, `draws_function`, and
#' `diagnose_function` reproduces `backend = "brms"` exactly, and each one can
#' be replaced individually to attach a different engine.
#'
#' A reportable fit needs all three. An adapter that supplies only a fit and a
#' draws function leaves the sampler evidence incomplete, and the fit is then
#' blocked rather than reported.
#'
#' @section Backend resolution:
#' Resolution happens before `pv_backend_brms_fit_function()` is called:
#' cmdstanr is used only when both its namespace and a configured CmdStan are
#' available, and rstan is selected only when cmdstanr is absent. A fit failure
#' is not retried against the other backend.
#'
#' @param formula Prepared stacked model formula, with the stacked outcome and
#'   weight columns already bound.
#' @param data Prepared stacked data frame.
#' @param family Response family, or `NULL` for `stats::gaussian()`.
#' @param prior Prior specification passed through to `brms::brm()`, or `NULL`.
#' @param chains,iter,warmup,cores,seed Sampler settings from [pv_control()].
#' @param backend Resolved Stan backend, `"cmdstanr"` or `"rstan"`.
#' @param file,file_refit Cache location and refit policy for `brms::brm()`.
#' @param ... Further arguments passed to `brms::brm()`.
#'
#' @returns A `brmsfit`.
#'
#' @seealso [pv_fit()] and [pv_control()] for how an adapter is selected.
#' @family backend adapters
#' @export
pv_backend_brms_fit_function <- function(formula, data, family, prior, chains,
                                         iter, warmup, cores, seed, backend,
                                         file, file_refit, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    pv_abort("`pv_control(backend = \"brms\")` requires the brms package; install brms or inject a `fit_function`.")
  }
  if (!is.character(backend) || length(backend) != 1L || is.na(backend) ||
      !backend %in% c("cmdstanr", "rstan")) {
    pv_abort("The resolved bundled brms backend must be `cmdstanr` or `rstan`.")
  }
  brms::brm(
    formula = formula,
    data = data,
    family = if (is.null(family)) stats::gaussian() else family,
    prior = prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    backend = backend,
    file = file,
    file_refit = file_refit,
    ...
  )
}

#' @rdname pv_backend_brms_fit_function
#'
#' @param fit The object returned by the fit function.
#'
#' @details
#' `pv_backend_brms_draws_function()` returns the fixed-effect (`b_*`) and
#' residual-scale (`sigma`) draws as a plain base matrix, the draw format the
#' calibration layer expects. A `posterior::draws_matrix` is converted rather
#' than passed through, because the calibration validator compares draw blocks
#' by value and a classed matrix would fail that comparison on class alone.
#'
#' @export
pv_backend_brms_draws_function <- function(fit, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    pv_abort("The bundled brms backend requires the posterior package to extract draws.")
  }
  dm <- posterior::as_draws_matrix(fit)
  dm <- dm[, grep("^(b_|sigma$)", colnames(dm)), drop = FALSE]
  matrix(as.numeric(dm), nrow = nrow(dm),
         dimnames = list(NULL, colnames(dm)))
}

#' @rdname pv_backend_brms_fit_function
#'
#' @param draws_array_function,summary_function,nuts_function Optional
#'   replacements for `posterior::as_draws_array()`, the
#'   `posterior::summarise_draws()` call that produces R-hat and ESS, and
#'   `brms::nuts_params()`. Supplying these lets an engine report diagnostics
#'   without depending on those packages.
#'
#' @details
#' `pv_backend_brms_sampler_diagnostics()` reports maximum R-hat, minimum bulk
#' and tail ESS, divergence count, chains, and post-warmup draws per chain. When
#' a required namespace is missing or the extraction fails, it returns an
#' explicit incomplete record instead of silently omitting the evidence, and the
#' fit is blocked downstream.
#'
#' @export
pv_backend_brms_sampler_diagnostics <- function(
  fit,
  draws_array_function = NULL,
  summary_function = NULL,
  nuts_function = NULL
) {
  source <- "bundled_brms_posterior_and_nuts"
  if (is.null(draws_array_function)) {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      return(pv_sampler_diagnostics_incomplete(
        source,
        "posterior_namespace_unavailable"
      ))
    }
    draws_array_function <- posterior::as_draws_array
  }
  if (is.null(summary_function)) {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      return(pv_sampler_diagnostics_incomplete(
        source,
        "posterior_namespace_unavailable"
      ))
    }
    summary_function <- function(draws) {
      posterior::summarise_draws(
        draws,
        rhat = posterior::rhat,
        ess_bulk = posterior::ess_bulk,
        ess_tail = posterior::ess_tail
      )
    }
  }
  if (is.null(nuts_function)) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      return(pv_sampler_diagnostics_incomplete(
        source,
        "brms_namespace_unavailable"
      ))
    }
    nuts_function <- brms::nuts_params
  }

  failed <- FALSE
  payload <- tryCatch({
    draws <- draws_array_function(fit)
    dims <- dim(draws)
    if (length(dims) != 3L || any(!is.finite(dims)) || any(dims <= 0L)) {
      stop("invalid draws array")
    }
    summary <- summary_function(draws)
    if (!is.data.frame(summary) ||
        !all(c("rhat", "ess_bulk", "ess_tail") %in% names(summary)) ||
        nrow(summary) == 0L) {
      stop("invalid posterior summary")
    }
    nuts <- nuts_function(fit)
    if (!is.data.frame(nuts) ||
        !all(c("Parameter", "Value") %in% names(nuts))) {
      stop("invalid NUTS parameters")
    }
    divergent <- nuts$Value[nuts$Parameter == "divergent__"]
    list(
      rhat_max = max(summary$rhat),
      ess_bulk_min = min(summary$ess_bulk),
      ess_tail_min = min(summary$ess_tail),
      divergences = if (length(divergent) == 0L) NA_real_ else sum(divergent),
      chains = dims[[2L]],
      post_warmup_draws_per_chain = dims[[1L]]
    )
  }, error = function(e) {
    failed <<- TRUE
    NULL
  })
  if (failed) {
    return(pv_sampler_diagnostics_incomplete(
      source,
      "diagnostic_extraction_failed"
    ))
  }
  pv_sampler_diagnostics_normalize(payload, diagnostic_source = source)
}
