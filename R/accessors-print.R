# Accessors check the fit's stored checksum instead of repeating every
# method-specific check. Fits that keep a backend fit object cannot use this
# shortcut, so `validate_pvstackr_fit()` runs the full check for them.
pv_validate_fit_for_access <- function(x) {
  validate_pvstackr_fit(x, tier = "cheap")
}

#' Get the estimate table of a fit
#'
#' `get_estimates()` returns the table of fixed-effect estimates of a fit,
#' with their standard errors, degrees of freedom and intervals. For a method
#' comparison it returns the table that lines up the estimates of the
#' compared fits.
#'
#' @details
#' Only the fixed effects are reported. `coverage_claim_allowed` says whether
#' pvstackr's reporting rule lets you read the interval of a row as a
#' confidence interval with nominal coverage; [pv_fit()] gives the rule, and [pvstackr_object_contracts]
#' describes every column.
#'
#' @param x A fit (class `pvstackr_fit`) from [pv_fit()] or a method function,
#'   or a method comparison from [pv_compare_methods()].
#' @param ... Ignored.
#'
#' @returns A data frame with one row per fixed effect. For a `stack_direct`
#'   fit it has 17 columns: `term`, `estimate`, `se`, `std.error`, `df`,
#'   `df_method`, `df_complete`, `conf_level`, `conf_low`, `conf_high`,
#'   `conf.low`, `conf.high`, `interval_role`, `coverage_claim_allowed`,
#'   `parameter_scope`, `target_source` and `target_hash`. `per_pv` and
#'   `stack_psis` fits add further columns. The fraction of missing
#'   information is not a column; for `stack_direct` and `per_pv` fits it is
#'   `get_target(fit)$fmi`. A blocked fit gives an empty data frame. For a
#'   method comparison, the data frame has one row per method and fixed
#'   effect (see [pv_compare_methods()]).
#'
#'   `get_estimates()` stops with an error if the fit or comparison was
#'   changed after it was created, and for the inspection object that
#'   [pv_migrate_legacy_psis_fit()] makes from a `stack_psis` fit of an
#'   earlier pvstackr version, which has no estimates.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit
#'   head(get_estimates(fit))
#' }
#' @family pvstackr-accessors
#' @seealso [pvstackr_object_contracts] for the columns of the estimate table.
#' @export
get_estimates <- function(x, ...) {
  UseMethod("get_estimates")
}

#' @rdname get_estimates
#' @export
get_estimates.default <- function(x, ...) {
  pv_abort("No `get_estimates()` method is available for this object.")
}

#' @rdname get_estimates
#' @export
get_estimates.pvstackr_fit <- function(x, ...) {
  pv_validate_fit_for_access(x)
  x$estimates
}

#' @rdname get_estimates
#' @export
get_estimates.pvstackr_method_comparison <- function(x, ...) {
  validate_pvstackr_method_comparison(x)
  x$estimate_table
}

#' Get the target of a fit
#'
#' `get_target()` returns the target that the estimates of a fit come from.
#' For a `stack_direct` fit this is the BRR-Fay target from [pv_brr_target()],
#' to which the stacked fit was calibrated; for a `per_pv` fit it is the
#' Rubin's-rules combination of the draws of the per-plausible-value fits. A
#' `stack_psis` fit has no target object.
#'
#' @param x A fit (class `pvstackr_fit`), or a target from [pv_brr_target()].
#' @param ... Ignored.
#'
#' @returns For a `stack_direct` fit, its `pvstackr_brr_target` object, also
#'   when the fit is blocked. [pv_brr_target()] lists its elements, such as
#'   `beta`, `T_MI`, `df` and `fmi` (the fraction of missing information,
#'   which the estimate table does not have). For a `per_pv` fit, a
#'   `pvstackr_reference_pool` object with the same kind of elements. For a
#'   `stack_psis` fit, `NULL`: its estimate table has
#'   `target_source = "stack_psis_rubin_pooling"`, but that is only a label
#'   for the combined result, not a target object. Given a target,
#'   `get_target()` checks it and returns it unchanged.
#'
#'   `get_target()` stops with an error if the fit or target was changed
#'   after it was created. For the inspection object that
#'   [pv_migrate_legacy_psis_fit()] makes from a `stack_psis` fit of an
#'   earlier pvstackr version, it returns `NULL`.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit
#'   get_target(fit)
#' }
#' @family pvstackr-accessors
#' @seealso [pv_brr_target()] for the elements of a target.
#' @export
get_target <- function(x, ...) {
  UseMethod("get_target")
}

#' @rdname get_target
#' @export
get_target.default <- function(x, ...) {
  pv_abort("No `get_target()` method is available for this object.")
}

#' @rdname get_target
#' @export
get_target.pvstackr_fit <- function(x, ...) {
  pv_validate_fit_for_access(x)
  x$target
}

#' @rdname get_target
#' @export
get_target.pvstackr_brr_target <- function(x, ...) {
  validate_pvstackr_brr_target(x)
  x
}

#' Get the calibrated draws of a fit
#'
#' `get_draws()` returns the calibrated fixed-effect draws of a `stack_direct`
#' fit: the draws of the stacked fit after the Cholesky calibration
#' correction (CCC), which gives them the mean and covariance of the target.
#'
#' @details
#' A fit keeps these draws only with `return_draws = TRUE`, the default of
#' [pv_control()]. `per_pv` and `stack_psis` fits keep their draws in
#' [get_diagnostics()] instead; "The fit object" in
#' [pvstackr_object_contracts] says where.
#'
#' @param x A fit (class `pvstackr_fit`).
#' @param ... Ignored.
#'
#' @returns A numeric matrix with one row per draw of the stacked fit and one
#'   column per fixed effect, named as in the estimate table (such as
#'   `b_Intercept`). `NULL` for a fit made with `return_draws = FALSE` (such
#'   as the bundled example fit), for a blocked fit, and for `per_pv` and
#'   `stack_psis` fits.
#'
#'   `get_draws()` stops with an error if the fit was changed after it was
#'   created, and for the inspection object that
#'   [pv_migrate_legacy_psis_fit()] makes from a `stack_psis` fit of an
#'   earlier pvstackr version, which has no draws.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit
#'   get_draws(fit)   # NULL: the example fit was saved without its draws
#' }
#' @family pvstackr-accessors
#' @seealso [pv_fit_direct()] for how the draws are calibrated.
#' @export
get_draws <- function(x, ...) {
  UseMethod("get_draws")
}

#' @rdname get_draws
#' @export
get_draws.default <- function(x, ...) {
  pv_abort("No `get_draws()` method is available for this object.")
}

#' @rdname get_draws
#' @export
get_draws.pvstackr_fit <- function(x, ...) {
  pv_validate_fit_for_access(x)
  x$draws
}

#' Get the diagnostics of a fit or comparison
#'
#' `get_diagnostics()` returns the diagnostics list of a fit or a method
#' comparison: the results of pvstackr's checks and the details of a method
#' that are not in the estimate table.
#'
#' @param x A fit (class `pvstackr_fit`) or a method comparison from
#'   [pv_compare_methods()].
#' @param ... Ignored.
#'
#' @returns A named list whose elements depend on the object:
#'
#' - A `stack_direct` fit: `preflight` (the check that `data`, `formula` and
#'   `target` match), `sampler` and `sampler_gate` (the sampler diagnostics
#'   and pvstackr's check of them), `stack_fit` and `stack_fit_warnings`
#'   (records and notes from the stacked fit) and `ccc` (the calibration
#'   diagnostics, such as `delta_c_max` and `kappa_A`). A blocked fit has
#'   only `preflight`, `sampler`, `sampler_gate`, `redaction` (what was
#'   removed) and, when the calibration check blocked it, `ccc`, with its
#'   values grouped as `center`, `conditioning`, `residual` and `prior`.
#' - A `per_pv` fit: `reference` (a record of the per-plausible-value fits
#'   and, when kept, their draws) and `pooling` (their Rubin's-rules
#'   combination).
#' - A `stack_psis` fit: `psis` (the Pareto k-hat values, the decision and
#'   the weight diagnostics), `pooling` and `weighted` (the weighted result
#'   of each plausible value and, when kept, the stacked draws and the
#'   normalized weights). A blocked fit has only `psis` and `redaction`.
#' - A method comparison: `reference_method`, `methods`, `statuses`,
#'   `blocked_methods`, `warning_methods`, `agreement`,
#'   `method_diagnostics`, `timing` and `target_overlap` (see
#'   [pv_compare_methods()]).
#' - The inspection object that [pv_migrate_legacy_psis_fit()] makes from a
#'   `stack_psis` fit of an earlier pvstackr version: `psis` (the Pareto
#'   k-hat values and the decision) and `redaction` (what was removed).
#'
#' `get_diagnostics()` stops with an error if the fit or comparison was
#' changed after it was created.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit
#'   names(get_diagnostics(fit))
#' }
#' @family pvstackr-accessors
#' @seealso [pv_fit_direct()] and [pvstackr_object_contracts] for the
#'   elements of a fit's diagnostics; "Status and checks" in
#'   [pvstackr_object_contracts] for the thresholds of the checks and their
#'   reason codes.
#' @export
get_diagnostics <- function(x, ...) {
  UseMethod("get_diagnostics")
}

#' @rdname get_diagnostics
#' @export
get_diagnostics.default <- function(x, ...) {
  pv_abort("No `get_diagnostics()` method is available for this object.")
}

#' @rdname get_diagnostics
#' @export
get_diagnostics.pvstackr_fit <- function(x, ...) {
  pv_validate_fit_for_access(x)
  x$diagnostics
}

#' @rdname get_diagnostics
#' @export
get_diagnostics.pvstackr_method_comparison <- function(x, ...) {
  validate_pvstackr_method_comparison(x)
  x$diagnostics
}

pv_estimate_print_columns <- function(estimates) {
  columns <- intersect(c("term", "estimate", "se", "df", "conf_low", "conf_high"), names(estimates))
  estimates[, columns, drop = FALSE]
}

pv_fit_target_source <- function(x) {
  if (!is.null(x$target$target_source)) {
    return(x$target$target_source)
  }
  if (!is.null(x$diagnostics$pooling$target_source)) {
    return(x$diagnostics$pooling$target_source)
  }
  if (!is.null(x$diagnostics$psis)) {
    return("none")
  }
  NA_character_
}

pv_draw_dim_label <- function(draws) {
  if (is.null(draws)) {
    return("not retained")
  }
  paste0(nrow(draws), " x ", ncol(draws))
}

pv_interval_note <- function(estimates) {
  if (!is.data.frame(estimates) || !"coverage_claim_allowed" %in% names(estimates) ||
      nrow(estimates) == 0L) {
    return(NA_character_)
  }
  coverage <- estimates$coverage_claim_allowed
  coverage <- coverage[!is.na(coverage)]
  if (length(coverage) == 0L || all(coverage)) {
    return(NA_character_)
  }
  if (all(!coverage)) {
    return("intervals are descriptive rather than coverage-claimable.")
  }
  "some intervals are descriptive rather than coverage-claimable."
}

pv_print_interval_note <- function(note) {
  if (is.character(note) && length(note) == 1L && !is.na(note) && nzchar(note)) {
    cat("  interval note: ", note, "\n", sep = "")
  }
}

pv_comparison_provenance_note <- function(diagnostics) {
  note <- diagnostics$target_overlap$independence_caveat %||% NA_character_
  if (is.character(note) && length(note) == 1L && !is.na(note) && nzchar(note)) {
    return(note)
  }
  NA_character_
}

pv_print_provenance_note <- function(note) {
  if (is.character(note) && length(note) == 1L && !is.na(note) && nzchar(note)) {
    cat("  provenance note: ", note, "\n", sep = "")
  }
}

pv_fit_summary_reportability_fit <- function(x) {
  if (!identical(x$method, "stack_psis")) {
    return(NULL)
  }
  validate_pvstackr_fit(x, tier = "deep")
  if (identical(x$status, "blocked")) {
    return(x)
  }
  projected_control <- x$control
  projected_control$return_draws <- FALSE
  projected_control$keep_data <- FALSE
  projected_control$keep_backend_fit <- FALSE
  projected_control$keep_log_lik <- FALSE
  projected_control <- pv_validate_control(projected_control)
  projected_diagnostics <- x$diagnostics
  projected_diagnostics$weighted["proposal_draws"] <- list(NULL)
  projected_diagnostics$weighted["weights"] <- list(NULL)
  projected_diagnostics$psis$weight_diagnostic_authority <-
    "owned_stamp_bounded_projection"
  new_pvstackr_fit(
    method = "stack_psis",
    estimates = x$estimates,
    diagnostics = projected_diagnostics,
    status = x$status,
    control = projected_control,
    reason_codes = x$reason_codes,
    provenance = list(
      wrapper_function = "pv_fit_stack_psis",
      stacked_source = "stacked_draws",
      psis_source = x$diagnostics$psis$source,
      pooling_hash = x$diagnostics$pooling$pooling_hash
    ),
    warnings = x$warnings
  )
}

pv_fit_summary <- function(x) {
  pv_validate_fit_for_access(x)
  source_reportability_fit <- pv_fit_summary_reportability_fit(x)
  validation_source <- source_reportability_fit %||% x
  estimates <- validation_source$estimates
  diagnostics <- validation_source$diagnostics
  list(
    method = x$method,
    status = x$status,
    reason_codes = x$reason_codes,
    warnings = x$warnings,
    n_terms = nrow(estimates),
    terms = if ("term" %in% names(estimates)) as.character(estimates$term) else character(),
    has_target = !is.null(x$target),
    has_draws = !is.null(x$draws),
    draw_dim = if (is.null(x$draws)) c(0L, 0L) else dim(x$draws),
    target_source = pv_fit_target_source(x),
    diagnostic_keys = names(diagnostics),
    interval_note = pv_interval_note(estimates),
    estimates = estimates,
    diagnostics = diagnostics,
    schema_version = x$schema_version,
    summary_schema_version = "0.2.0",
    source_validation = list(
      schema_version = validation_source$validation$schema_version,
      stamp = validation_source$validation$stamp
    ),
    source_reportability_fit = source_reportability_fit,
    validation = NULL
  )
}

pv_summary_validation_sentinel <- function() {
  paste0("sha256:", strrep("0", 64L))
}

pv_summary_validation_schema <- function(kind = c("fit", "comparison")) {
  kind <- match.arg(kind)
  paste0("pvstackr_", kind, "_summary_validation_v1")
}

pv_summary_validation_record <- function(
  kind = c("fit", "comparison"),
  stamp = pv_summary_validation_sentinel()
) {
  kind <- match.arg(kind)
  list(
    schema_version = pv_summary_validation_schema(kind),
    policy_id = "source_stamp_plus_owned_summary_sha256_v1",
    canonicalizer_id = "r_xdr_v2_summary_owned_payload_v1",
    stamp = stamp
  )
}

pv_summary_validation_digest <- function(x, kind = c("fit", "comparison")) {
  kind <- match.arg(kind)
  projected <- x
  projected$validation$stamp <- pv_summary_validation_sentinel()
  if (identical(kind, "fit") &&
      is.list(projected$diagnostics) &&
      is.list(projected$diagnostics$reference) &&
      !is.null(projected$diagnostics$reference$backend_fits)) {
    backend_fits <- projected$diagnostics$reference$backend_fits
    projected$diagnostics$reference$backend_fits <- list(
      marker = "pvstackr::opaque_backend_fits",
      length = length(backend_fits),
      names = names(backend_fits)
    )
  }
  payload <- pv_validation_payload_bytes(projected, "Summary validation payload")
  bytes <- c(
    charToRaw(paste0("pvstackr-", kind, "-summary-validation-v1")),
    as.raw(0L),
    payload
  )
  paste0("sha256:", digest::digest(bytes, algo = "sha256", serialize = FALSE))
}

pv_summary_issue_validation_stamp <- function(x, kind = c("fit", "comparison")) {
  kind <- match.arg(kind)
  x$validation <- pv_summary_validation_record(kind)
  x$validation$stamp <- pv_summary_validation_digest(x, kind)
  x
}

pv_validate_summary_stamp <- function(x, kind = c("fit", "comparison")) {
  kind <- match.arg(kind)
  validation <- x$validation
  if (!is.list(validation) ||
      !identical(
        validation,
        pv_summary_validation_record(kind, validation$stamp %||% "")
      ) || !grepl("^sha256:[0-9a-f]{64}$", validation$stamp) ||
      !identical(validation$stamp, pv_summary_validation_digest(x, kind))) {
    pv_abort("Summary validation record or owned-payload stamp is invalid.")
  }
  invisible(validation)
}

pv_validate_summary_source <- function(source, expected_schema) {
  fields <- c("schema_version", "stamp")
  if (!is.list(source) || !identical(names(source), fields) ||
      !identical(attributes(source), list(names = fields)) ||
      !identical(source$schema_version, expected_schema) ||
      !is.character(source$stamp) || length(source$stamp) != 1L ||
      is.na(source$stamp) || !grepl("^sha256:[0-9a-f]{64}$", source$stamp)) {
    pv_abort("Summary source-validation record is noncanonical.")
  }
  invisible(source)
}

pv_validate_stack_psis_fit_summary <- function(x) {
  if (!identical(x$method, "stack_psis")) {
    return(invisible(x))
  }
  if (identical(x$status, "warning")) {
    pv_abort("Warning-status stack_psis summaries are legacy unsafe and inspection-only.")
  }
  psis <- x$diagnostics$psis
  if (!is.list(psis)) {
    pv_abort("stack_psis summary requires PSIS decision evidence.")
  }
  if (identical(x$status, "ok")) {
    threshold <- psis$threshold
    pareto_k <- psis$pareto_k
    if (!identical(psis$status, "ok") || !is.numeric(threshold) ||
        length(threshold) != 1L || !is.finite(threshold) || threshold <= 0 ||
        threshold > pv_legacy_psis_hard_threshold() ||
        !is.numeric(pareto_k) || length(pareto_k) < 2L ||
        any(!is.finite(pareto_k)) || any(pareto_k >= threshold) ||
        x$n_terms < 1L || !is.data.frame(x$estimates) ||
        nrow(x$estimates) != x$n_terms) {
      pv_abort("Reportable stack_psis summary fails the immutable Pareto-k firewall.")
    }
    return(invisible(x))
  }
  if (!identical(x$status, "blocked") ||
      !psis$status %in% c(
        "failed", "not_evaluated", "unsmoothed", "provenance_incomplete"
      ) ||
      x$n_terms != 0L || !is.data.frame(x$estimates) ||
      nrow(x$estimates) != 0L || !identical(x$has_draws, FALSE)) {
    pv_abort("Blocked stack_psis summary must contain no reportable numeric output.")
  }
  invisible(x)
}

pv_validate_fit_summary_for_print <- function(x) {
  current_fields <- c(
    "method", "status", "reason_codes", "warnings", "n_terms", "terms",
    "has_target", "has_draws", "draw_dim", "target_source",
    "diagnostic_keys", "interval_note", "estimates", "diagnostics",
    "schema_version", "summary_schema_version", "source_validation",
    "source_reportability_fit", "validation"
  )
  is_current <- all(c(
    "summary_schema_version", "source_validation", "validation"
  ) %in% names(x))
  if (!is_current) {
    legacy_fields <- c(
      "method", "status", "reason_codes", "warnings", "n_terms", "terms",
      "has_target", "has_draws", "draw_dim", "target_source",
      "diagnostic_keys", "interval_note", "estimates", "diagnostics",
      "schema_version"
    )
    attrs <- attributes(x)
    if (!is.list(x) || !identical(names(x), legacy_fields) ||
        !identical(names(attrs), c("names", "class")) ||
        !identical(attrs$names, legacy_fields) ||
        !identical(attrs$class, c("summary.pvstackr_fit", "list")) ||
        !is.character(x$method) || length(x$method) != 1L ||
        is.na(x$method) || !x$method %in% pv_allowed_methods()) {
      pv_abort("Legacy fit summary schema is unrecognized and cannot be printed.")
    }
    if (identical(x$method, "stack_psis")) {
      pv_abort(
        "Legacy stack_psis summaries are inspection-only and must be rebuilt from a current validated fit."
      )
    }
    return(invisible(x))
  }
  root_attrs <- attributes(x)
  if (!identical(names(x), current_fields) ||
      !identical(names(root_attrs), c("names", "class")) ||
      !identical(root_attrs$names, current_fields) ||
      !identical(root_attrs$class, c("summary.pvstackr_fit", "list")) ||
      !identical(x$summary_schema_version, "0.2.0")) {
    pv_abort("Current fit summary fields, order, class, and schema must be exact.")
  }
  pv_validate_summary_source(x$source_validation, pv_fit_validation_schema())
  if (identical(x$method, "stack_psis")) {
    source <- x$source_reportability_fit
    if (!inherits(source, "pvstackr_fit")) {
      pv_abort("Current stack_psis summary requires a validated source reportability fit.")
    }
    validate_pvstackr_fit(source, tier = "deep")
    if (!identical(
          x$source_validation,
          list(
            schema_version = source$validation$schema_version,
            stamp = source$validation$stamp
          )
        ) || !identical(x$method, source$method) ||
        !identical(x$status, source$status) ||
        !identical(x$reason_codes, source$reason_codes) ||
        !identical(x$warnings, source$warnings) ||
        !identical(x$estimates, source$estimates) ||
        !identical(x$diagnostics, source$diagnostics) ||
        !identical(x$n_terms, nrow(source$estimates)) ||
        !identical(
          x$terms,
          if ("term" %in% names(source$estimates)) {
            as.character(source$estimates$term)
          } else {
            character()
          }
        ) || !identical(x$has_target, !is.null(source$target)) ||
        !identical(x$has_draws, !is.null(source$draws)) ||
        !identical(
          x$draw_dim,
          if (is.null(source$draws)) c(0L, 0L) else dim(source$draws)
        ) || !identical(x$target_source, pv_fit_target_source(source)) ||
        !identical(x$diagnostic_keys, names(source$diagnostics)) ||
        !identical(x$interval_note, pv_interval_note(source$estimates)) ||
        !identical(x$schema_version, source$schema_version)) {
      pv_abort("Current stack_psis summary must exactly reproduce its validated source reportability fit.")
    }
  } else if (!is.null(x$source_reportability_fit)) {
    pv_abort("Only stack_psis summaries may carry a source reportability fit.")
  }
  pv_validate_summary_stamp(x, "fit")
  pv_validate_stack_psis_fit_summary(x)
  invisible(x)
}

#' Print and summarize a fit
#'
#' `print()` shows a short overview of a fit (class `pvstackr_fit`, see
#' [pv_fit()]). `summary()` returns the same overview, together with the
#' estimate table and the diagnostics, as a list of class
#' `summary.pvstackr_fit`; its `print()` method also shows the estimates.
#'
#' @details
#' The overview has one line each for the method, the status, the number of
#' fixed effects, the target (its `target_source` label, or "none" when
#' [get_target()] returns `NULL`, as for a `stack_psis` fit), the calibrated
#' draws (their numbers of rows and columns, or "not retained" when
#' [get_draws()] returns `NULL`) and the names of the diagnostics. A line
#' beginning "interval note:" follows when some or all intervals are
#' descriptive, and the reason codes and the number of warnings are shown
#' when there are any. The `print()` method of the summary adds the columns
#' `term`, `estimate`, `se`, `df`, `conf_low` and `conf_high` of the
#' estimate table.
#'
#' The three methods stop with an error if the fit or the summary was
#' changed after it was created, and for a fit saved by pvstackr 0.1.x,
#' which has to be made again. A summary of a `stack_psis` fit that was saved
#' by an earlier version of pvstackr also stops with an error when printed
#' (see [pv_migrate_legacy_psis_fit()]).
#'
#' @param object,x A fit (class `pvstackr_fit`); for the `print()` method of
#'   a summary, the list returned by `summary()`.
#' @param ... Ignored.
#'
#' @returns `print()` returns its input invisibly. `summary()` returns a list
#'   of class `summary.pvstackr_fit` with these elements:
#'
#'   - `method`, `status`, `reason_codes` and `warnings`: as in the fit.
#'   - `n_terms` and `terms`: the number and the names of the fixed effects
#'     in the estimate table (`0` and `character()` for a blocked fit).
#'   - `has_target` and `target_source`: whether the fit has a target object
#'     ([get_target()] does not return `NULL`), and its `target_source`
#'     label, which is `"none"` for a `stack_psis` fit.
#'   - `has_draws` and `draw_dim`: whether the fit has calibrated draws
#'     ([get_draws()] does not return `NULL`), and their numbers of rows and
#'     columns (`c(0L, 0L)` when it has none).
#'   - `diagnostic_keys`: the names of the list that [get_diagnostics()]
#'     returns.
#'   - `interval_note`: the text of the line beginning "interval note:", or
#'     `NA` when there is no such line.
#'   - `estimates` and `diagnostics`: the estimate table and the diagnostics
#'     list.
#'   - `schema_version`, `summary_schema_version`, `source_validation`,
#'     `source_reportability_fit` and `validation`: the format versions of
#'     the fit and of the summary; the checksum of the fit that the summary
#'     was made from; for a `stack_psis` fit, a copy of the fit without its
#'     stacked draws and weights, from which `estimates` and `diagnostics`
#'     are taken (`NULL` for the other methods); and a SHA-256 checksum of
#'     the summary, which `print()` recomputes to detect changes.
#' @seealso [get_estimates()] and [get_diagnostics()] to read a fit.
#' @name pvstackr_fit_summary
NULL

#' @rdname pvstackr_fit_summary
#' @export
print.pvstackr_fit <- function(x, ...) {
  pv_validate_fit_for_access(x)
  n_terms <- nrow(x$estimates)
  cat("pvstackr fit\n")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  status: ", x$status, "\n", sep = "")
  cat("  fixed effects: ", n_terms, "\n", sep = "")
  cat("  target: ", if (!is.null(x$target)) pv_fit_target_source(x) else "none", "\n", sep = "")
  cat("  draws: ", pv_draw_dim_label(x$draws), "\n", sep = "")
  if (length(names(x$diagnostics)) > 0L) {
    cat("  diagnostics: ", paste(names(x$diagnostics), collapse = ", "), "\n", sep = "")
  }
  pv_print_interval_note(pv_interval_note(x$estimates))
  if (length(x$reason_codes) > 0L) {
    cat("  reason codes: ", paste(x$reason_codes, collapse = ", "), "\n", sep = "")
  }
  if (length(x$warnings) > 0L) {
    cat("  warnings: ", length(x$warnings), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname pvstackr_fit_summary
#' @export
summary.pvstackr_fit <- function(object, ...) {
  summary <- pv_fit_summary(object)
  class(summary) <- c("summary.pvstackr_fit", "list")
  summary <- pv_summary_issue_validation_stamp(summary, "fit")
  pv_validate_fit_summary_for_print(summary)
  summary
}

#' @rdname pvstackr_fit_summary
#' @export
print.summary.pvstackr_fit <- function(x, ...) {
  pv_validate_fit_summary_for_print(x)
  cat("pvstackr fit summary\n")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  status: ", x$status, "\n", sep = "")
  cat("  fixed effects: ", x$n_terms, "\n", sep = "")
  cat("  target: ", if (isTRUE(x$has_target)) x$target_source else "none", "\n", sep = "")
  cat("  draws: ", if (isTRUE(x$has_draws)) paste0(x$draw_dim[1L], " x ", x$draw_dim[2L]) else "not retained", "\n", sep = "")
  if (length(x$diagnostic_keys) > 0L) {
    cat("  diagnostics: ", paste(x$diagnostic_keys, collapse = ", "), "\n", sep = "")
  }
  pv_print_interval_note(x$interval_note)
  if (length(x$reason_codes) > 0L) {
    cat("  reason codes: ", paste(x$reason_codes, collapse = ", "), "\n", sep = "")
  }
  if (length(x$warnings) > 0L) {
    cat("  warnings: ", length(x$warnings), "\n", sep = "")
  }
  if (x$n_terms > 0L) {
    print(pv_estimate_print_columns(x$estimates), row.names = FALSE)
  }
  invisible(x)
}

pv_comparison_summary <- function(x) {
  validate_pvstackr_method_comparison(x)
  estimates <- x$estimate_table
  diagnostics <- x$diagnostics
  source_reportability_comparison <- NULL
  if (!is.null(x$validation)) {
    source_reportability_comparison <- x
    source_reportability_comparison["fits"] <- list(NULL)
    source_reportability_comparison <-
      pv_comparison_issue_validation_stamp(source_reportability_comparison)
    validate_pvstackr_method_comparison(source_reportability_comparison)
  }
  out <- list(
    reference_method = x$reference_method,
    methods = x$methods,
    method_labels = x$method_labels,
    n_methods = length(x$method_labels),
    n_terms = length(unique(estimates$term)),
    blocked_methods = diagnostics$blocked_methods,
    warning_methods = diagnostics$warning_methods,
    interval_note = pv_interval_note(estimates),
    provenance_note = pv_comparison_provenance_note(diagnostics),
    estimate_table = estimates,
    diagnostic_table = x$diagnostic_table,
    agreement = x$agreement,
    timing = x$timing,
    schema_version = x$schema_version
  )
  if (is.null(x$validation)) {
    return(out)
  }
  c(out, list(
    summary_schema_version = "0.2.0",
    source_validation = list(
      schema_version = source_reportability_comparison$validation$schema_version,
      stamp = source_reportability_comparison$validation$stamp
    ),
    source_reportability_comparison = source_reportability_comparison,
    validation = NULL
  ))
}

pv_validate_comparison_summary_for_print <- function(x) {
  current_fields <- c(
    "reference_method", "methods", "method_labels", "n_methods", "n_terms",
    "blocked_methods", "warning_methods", "interval_note",
    "provenance_note", "estimate_table", "diagnostic_table", "agreement",
    "timing", "schema_version", "summary_schema_version",
    "source_validation", "source_reportability_comparison", "validation"
  )
  is_current <- all(c(
    "summary_schema_version", "source_validation", "validation"
  ) %in% names(x))
  if (!is_current) {
    legacy_fields <- c(
      "reference_method", "methods", "method_labels", "n_methods", "n_terms",
      "blocked_methods", "warning_methods", "interval_note",
      "provenance_note", "estimate_table", "diagnostic_table", "agreement",
      "timing", "schema_version"
    )
    attrs <- attributes(x)
    valid_methods <- is.character(x$methods) && length(x$methods) >= 2L &&
      !anyNA(x$methods) && all(x$methods %in% pv_allowed_methods())
    valid_tables <- is.data.frame(x$estimate_table) &&
      "method" %in% names(x$estimate_table) &&
      is.data.frame(x$diagnostic_table) &&
      "method" %in% names(x$diagnostic_table)
    if (!is.list(x) || !identical(names(x), legacy_fields) ||
        !identical(names(attrs), c("names", "class")) ||
        !identical(attrs$names, legacy_fields) ||
        !identical(
          attrs$class,
          c("summary.pvstackr_method_comparison", "list")
        ) || !valid_methods || !valid_tables) {
      pv_abort("Legacy method-comparison summary schema is unrecognized and cannot be printed.")
    }
    observed_methods <- c(
      unname(x$methods),
      as.character(x$estimate_table$method),
      as.character(x$diagnostic_table$method)
    )
    if (any(observed_methods == "stack_psis", na.rm = TRUE)) {
      pv_abort(
        paste(
          "Legacy method-comparison summaries containing stack_psis are",
          "inspection-only and must be rebuilt from current validated fits."
        )
      )
    }
    return(invisible(x))
  }
  root_attrs <- attributes(x)
  if (!identical(names(x), current_fields) ||
      !identical(names(root_attrs), c("names", "class")) ||
      !identical(root_attrs$names, current_fields) ||
      !identical(
        root_attrs$class,
        c("summary.pvstackr_method_comparison", "list")
      ) || !identical(x$summary_schema_version, "0.2.0")) {
    pv_abort("Current method-comparison summary fields, order, class, and schema must be exact.")
  }
  pv_validate_summary_source(
    x$source_validation,
    pv_comparison_validation_schema()
  )
  source <- x$source_reportability_comparison
  if (!inherits(source, "pvstackr_method_comparison")) {
    pv_abort("Current method-comparison summary requires a validated source reportability comparison.")
  }
  validate_pvstackr_method_comparison(source)
  expected_interval_note <- pv_interval_note(source$estimate_table)
  expected_provenance_note <- pv_comparison_provenance_note(source$diagnostics)
  if (!identical(
        x$source_validation,
        list(
          schema_version = source$validation$schema_version,
          stamp = source$validation$stamp
        )
      ) || !is.null(source$fits) ||
      !identical(x$reference_method, source$reference_method) ||
      !identical(x$methods, source$methods) ||
      !identical(x$method_labels, source$method_labels) ||
      !identical(x$n_methods, length(source$method_labels)) ||
      !identical(x$n_terms, length(unique(source$estimate_table$term))) ||
      !identical(x$blocked_methods, source$diagnostics$blocked_methods) ||
      !identical(x$warning_methods, source$diagnostics$warning_methods) ||
      !identical(x$interval_note, expected_interval_note) ||
      !identical(x$provenance_note, expected_provenance_note) ||
      !identical(x$estimate_table, source$estimate_table) ||
      !identical(x$diagnostic_table, source$diagnostic_table) ||
      !identical(x$agreement, source$agreement) ||
      !identical(x$timing, source$timing) ||
      !identical(x$schema_version, source$schema_version)) {
    pv_abort("Current method-comparison summary must exactly reproduce its validated source comparison.")
  }
  pv_validate_summary_stamp(x, "comparison")
  pv_validate_stack_psis_derived_tables(
    x$estimate_table,
    x$diagnostic_table,
    context = "method-comparison summary"
  )
  invisible(x)
}

#' Print and summarize a method comparison
#'
#' `print()` shows a short overview of a comparison made by
#' [pv_compare_methods()]. `summary()` returns the same overview, together
#' with the tables of the comparison, as a list of class
#' `summary.pvstackr_method_comparison`; its `print()` method also shows the
#' agreement table.
#'
#' @details
#' The overview shows the label of the reference fit, the label and the
#' method of each fit (only their number in the summary), the number of
#' fixed effects, and the labels of blocked fits and of fits with status
#' `"warning"` when there are any. A line beginning "provenance note:"
#' follows for every comparison: it says that agreement is descriptive and
#' that a shared target, combined result or source label is not independent
#' confirmation ([pv_compare_methods()] explains this). A line beginning
#' "interval note:" is added when some or all intervals are descriptive.
#'
#' The three methods stop with an error if the comparison or the summary was
#' changed after it was created, and for a comparison saved by pvstackr
#' 0.1.x, which has to be made again from current fits. A summary that
#' contains a `stack_psis` fit and was saved by an earlier version of
#' pvstackr also stops with an error when printed (see
#' [pv_migrate_legacy_psis_fit()]).
#'
#' @param object,x A comparison (class `pvstackr_method_comparison`) from
#'   [pv_compare_methods()]; for the `print()` method of a summary, the list
#'   returned by `summary()`.
#' @param ... Ignored.
#'
#' @returns `print()` returns its input invisibly. `summary()` returns a list
#'   of class `summary.pvstackr_method_comparison` with these elements:
#'
#'   - `reference_method`: the label of the reference fit.
#'   - `methods`, `method_labels` and `n_methods`: the method of each fit
#'     (named by label), the labels and the number of fits.
#'   - `n_terms`: the number of fixed effects in the comparison.
#'   - `blocked_methods` and `warning_methods`: the labels of the fits with
#'     status `"blocked"` or `"warning"`.
#'   - `interval_note` and `provenance_note`: the texts of the two note
#'     lines; `interval_note` is `NA` when all intervals are
#'     coverage-claimable.
#'   - `estimate_table`, `diagnostic_table`, `agreement` and `timing`: the
#'     tables described in [pv_compare_methods()].
#'   - `schema_version`, `summary_schema_version`, `source_validation`,
#'     `source_reportability_comparison` and `validation`: the format
#'     versions of the comparison and of the summary; a copy of the
#'     comparison without its fits and the checksum of that copy; and a
#'     SHA-256 checksum of the summary, which `print()` recomputes to detect
#'     changes.
#' @seealso [get_estimates()] and [get_diagnostics()] to read a comparison.
#' @name pvstackr_method_comparison_summary
NULL

#' @rdname pvstackr_method_comparison_summary
#' @export
print.pvstackr_method_comparison <- function(x, ...) {
  validate_pvstackr_method_comparison(x)
  cat("pvstackr method comparison\n")
  cat("  reference: ", x$reference_method, "\n", sep = "")
  cat("  methods: ", paste(paste0(x$method_labels, "=", unname(x$methods)), collapse = ", "), "\n", sep = "")
  cat("  fixed effects: ", length(unique(x$estimate_table$term)), "\n", sep = "")
  if (length(x$diagnostics$blocked_methods) > 0L) {
    cat("  blocked: ", paste(x$diagnostics$blocked_methods, collapse = ", "), "\n", sep = "")
  }
  if (length(x$diagnostics$warning_methods) > 0L) {
    cat("  warnings: ", paste(x$diagnostics$warning_methods, collapse = ", "), "\n", sep = "")
  }
  pv_print_provenance_note(pv_comparison_provenance_note(x$diagnostics))
  pv_print_interval_note(pv_interval_note(x$estimate_table))
  invisible(x)
}

#' @rdname pvstackr_method_comparison_summary
#' @export
summary.pvstackr_method_comparison <- function(object, ...) {
  summary <- pv_comparison_summary(object)
  class(summary) <- c("summary.pvstackr_method_comparison", "list")
  if ("validation" %in% names(summary)) {
    summary <- pv_summary_issue_validation_stamp(summary, "comparison")
  }
  pv_validate_comparison_summary_for_print(summary)
  summary
}

#' @rdname pvstackr_method_comparison_summary
#' @export
print.summary.pvstackr_method_comparison <- function(x, ...) {
  pv_validate_comparison_summary_for_print(x)
  cat("pvstackr method comparison summary\n")
  cat("  reference: ", x$reference_method, "\n", sep = "")
  cat("  methods: ", x$n_methods, "\n", sep = "")
  cat("  fixed effects: ", x$n_terms, "\n", sep = "")
  if (length(x$blocked_methods) > 0L) {
    cat("  blocked: ", paste(x$blocked_methods, collapse = ", "), "\n", sep = "")
  }
  if (length(x$warning_methods) > 0L) {
    cat("  warnings: ", paste(x$warning_methods, collapse = ", "), "\n", sep = "")
  }
  pv_print_provenance_note(x$provenance_note)
  pv_print_interval_note(x$interval_note)
  print(x$agreement, row.names = FALSE)
  invisible(x)
}
