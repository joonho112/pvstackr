# Bundled backend adapters -----------------------------------------------------
#
# pv_fit(method = "stack_direct") accepts injected `fit_function`,
# `draws_function` and `diagnose_function` arguments, so that any Bayesian
# engine can drive the stacked fit. For the common case the package bundles a
# brms adapter, selected with `pv_control(backend = "brms")` and no
# `fit_function`. Its three functions are exported, so an injected adapter can
# reuse them: it fits the prepared stacked formula with `brms::brm()` (Gaussian
# family unless one is supplied), returns the fixed-effect and residual-scale
# draws as a plain base matrix, and reports the sampler diagnostics.
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

#' Bundled brms engine for the stacked fit
#'
#' These three functions are the brms engine that [pv_fit_direct()] uses when
#' `pv_control(backend = "brms")` is set and no `fit_function` is given:
#' `pv_backend_brms_fit_function()` fits the stacked model with `brms::brm()`,
#' `pv_backend_brms_draws_function()` returns the draws of that fit, and
#' `pv_backend_brms_sampler_diagnostics()` returns the sampler diagnostics
#' that pvstackr checks. They are exported so that you can pass them to
#' [pv_fit_direct()] yourself and replace one of them with your own function.
#'
#' @details
#' ## Replacing one of the functions
#'
#' To change one step (fitting the stacked model, extracting its draws or
#' computing its sampler diagnostics), pass all three functions to
#' [pv_fit_direct()] as `fit_function`, `draws_function` and
#' `diagnose_function`: your own function for that step and the functions
#' documented here for the others. [pv_fit_direct()] describes what each
#' function receives and returns, and what happens when one is missing. The
#' draws and diagnostics functions here expect a `brmsfit`; for a fit from
#' another program, `pv_backend_brms_sampler_diagnostics()` can still
#' compute the diagnostics with its three replacement functions. With
#' `pv_control(backend = "brms")` and no `fit_function`, only the draws
#' function can be replaced: a `diagnose_function` that you pass is not
#' used, and the fit records that it was ignored.
#'
#' When you pass the functions yourself, pvstackr checks no package and
#' creates no cache folder before it calls `fit_function`. A missing
#' posterior package is found only after sampling, when the draws are
#' extracted; create `cache_dir` first or set `cache_dir = NULL`. A fit made
#' this way has `get_diagnostics(fit)$sampler$diagnostic_source` equal to
#' `"injected_diagnose_function"`, also when your `diagnose_function` is
#' `pv_backend_brms_sampler_diagnostics()`.
#'
#' ## Stan backend
#'
#' `pv_backend_brms_fit_function()` passes `backend` to `brm()`. For the
#' bundled engine, pvstackr chooses it before sampling: `"cmdstanr"` when the
#' cmdstanr package is installed and CmdStan is configured, and `"rstan"` when
#' cmdstanr is not installed. If cmdstanr is installed but CmdStan is not
#' configured, the fit stops with an error instead of using rstan. When you
#' pass the function yourself, it receives the `backend` value of
#' [pv_control()]: `"cmdstanr"` is passed to `brm()` unchanged, and any other
#' value is replaced by the choice described above. A fit that fails is not
#' repeated with the other backend.
#'
#' @param formula The stacked formula built by [pv_fit_direct()], such as
#'   `.pvstackr_y | weights(.pvstackr_weight) ~ 0 + pvstackrMM001 + ...`.
#' @param data The stacked data built by [pv_fit_direct()]: one copy of the
#'   data per plausible value, with the columns used in `formula`.
#' @param family A family object, or `NULL` for `stats::gaussian()`.
#'   [pv_fit_direct()] always passes `stats::gaussian()`.
#' @param prior `NULL` or a prior table, passed to `brm()` unchanged.
#'   [pv_fit_direct()] describes the priors it accepts.
#' @param chains,iter,warmup,cores,seed Sampler settings passed to `brm()`;
#'   [pv_fit_direct()] takes them from [pv_control()].
#' @param backend The Stan backend for `brm()`. `"cmdstanr"` and `"rstan"` are
#'   used as given; any other value is replaced by the backend that the
#'   bundled engine would choose (see "Stan backend").
#' @param file,file_refit The cache file (a path without extension) and the
#'   refit rule, passed to `brm()`. [pv_fit_direct()] sets them from
#'   `cache_dir` and `cache_stem`: `file_refit = "on_change"`, or
#'   `file = NULL` and `file_refit = "never"` when caching is off.
#' @param ... Further arguments. `pv_backend_brms_fit_function()` passes them
#'   to `brm()` ([pv_fit_direct()] passes the elements of `additional_args`
#'   here); `pv_backend_brms_draws_function()` ignores them.
#'
#' @returns
#' `pv_backend_brms_fit_function()` returns the `brmsfit` object from `brm()`.
#'
#' `pv_backend_brms_draws_function()` returns a base R numeric matrix with one
#' row per draw after warmup, over all chains, and one column for each
#' variable whose name starts with `b_` and for `sigma`. For the stacked model
#' the `b_` columns are the fixed effects of the stacked formula
#' (`b_pvstackrMM001`, `b_pvstackrMM002`, ...), which pvstackr renames to the
#' names used in the target, such as `b_Intercept`.
#'
#' `pv_backend_brms_sampler_diagnostics()` returns the list that a fit made
#' with the bundled engine keeps as `get_diagnostics(fit)$sampler`:
#' `rhat_max`, `ess_bulk_min`, `ess_tail_min`, `ess_bulk_per_chain_min` and
#' `ess_tail_per_chain_min` (the two ESS values divided by `chains`),
#' `divergences`, `chains`, `post_warmup_draws_per_chain`,
#' `diagnostic_source` (`"bundled_brms_posterior_and_nuts"`),
#' `diagnostic_complete` (`TRUE` when every value was extracted) and
#' `diagnostic_reason_codes` (the reasons when it is `FALSE`).
#'
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
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
#' # Pass the three functions yourself and replace any one of them with your
#' # own. A live fit samples with Stan, so it is not run here.
#' \dontrun{
#' fit <- pv_fit_direct(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female, target = target,
#'   fit_function = pv_backend_brms_fit_function,
#'   draws_function = pv_backend_brms_draws_function,
#'   diagnose_function = pv_backend_brms_sampler_diagnostics,
#'   cache_dir = NULL
#' )
#' }
#'
#' # The sampler diagnostics record. Stand-in functions replace the posterior
#' # and brms functions, so this runs without those packages.
#' diagnostics <- pv_backend_brms_sampler_diagnostics(
#'   fit = NULL,
#'   draws_array_function = function(fit) array(0, dim = c(1000, 4, 3)),
#'   summary_function = function(draws) {
#'     data.frame(
#'       rhat = c(1.001, 1.004, 1.002),
#'       ess_bulk = c(3900, 4100, 3700),
#'       ess_tail = c(2900, 3100, 3000)
#'     )
#'   },
#'   nuts_function = function(fit) {
#'     data.frame(Parameter = "divergent__", Value = rep(0, 4000))
#'   }
#' )
#' str(diagnostics, digits.d = 4)
#'
#' @seealso [pv_fit_direct()] for the arguments `fit_function`,
#'   `draws_function` and `diagnose_function`, and [pv_control()] for
#'   `backend`.
#' @family pvstackr-backends
#' @export
pv_backend_brms_fit_function <- function(formula, data, family, prior, chains,
                                         iter, warmup, cores, seed, backend,
                                         file, file_refit, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    pv_abort("`pv_control(backend = \"brms\")` requires the brms package; install brms or inject a `fit_function`.")
  }
  if (!is.character(backend) || length(backend) != 1L || is.na(backend) ||
      !backend %in% c("cmdstanr", "rstan")) {
    # Reached through an injected adapter, where `backend` is whatever
    # `pv_control()` was given rather than an already-resolved Stan backend.
    # Resolve it here with the same policy the bundled route uses, so both
    # routes pick the engine the same way.
    backend <- pv_backend_resolve_brms_engine(
      pv_backend_cmdstan_state()
    )$resolved_backend
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
#' @param fit The object returned by the fitting function: a `brmsfit` for
#'   the functions documented here. With all three replacement functions
#'   below, `pv_backend_brms_sampler_diagnostics()` accepts any object that
#'   they accept.
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
#' @param draws_array_function,summary_function,nuts_function Functions that
#'   replace, in this order, `posterior::as_draws_array()`, the
#'   `posterior::summarise_draws()` call and `brms::nuts_params()`, or `NULL`
#'   (default) to use those. With all three, the diagnostics can be computed
#'   without the posterior and brms packages, for example for a fit from
#'   another program. `draws_array_function(fit)` must return an array with
#'   the dimensions draws, chains and variables; `summary_function(draws)` a
#'   data frame with the columns `rhat`, `ess_bulk` and `ess_tail` and one
#'   row per variable; and `nuts_function(fit)` a data frame with the columns
#'   `Parameter` and `Value`, whose rows with `Parameter == "divergent__"`
#'   hold 1 for a divergent transition and 0 otherwise.
#'
#' @details
#' ## Sampler diagnostics
#'
#' `pv_backend_brms_sampler_diagnostics()` computes R-hat and the bulk and
#' tail effective sample sizes (ESS) of every variable of the fit with
#' `posterior::summarise_draws()`, and reports the largest R-hat and the
#' smallest ESS values, which are totals over all chains. The number of
#' divergent transitions comes from `brms::nuts_params()`, and `chains` and
#' `post_warmup_draws_per_chain` from the dimensions of the draws. When the
#' posterior or brms package is missing, or a value cannot be extracted, it
#' does not stop with an error: it returns the list with missing values,
#' `diagnostic_complete = FALSE` and reason codes such as
#' `posterior_namespace_unavailable`, and [pv_fit_direct()] then blocks the
#' fit. The limits that pvstackr applies to these values are listed under
#' "Status and checks" in [pvstackr_object_contracts].
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
