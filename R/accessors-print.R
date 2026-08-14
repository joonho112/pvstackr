# Accessors verify the current owned payload stamp without repeating the full
# method-specific semantic recomputation. Opaque-backend fits are marked
# ineligible and `validate_pvstackr_fit()` automatically falls back to deep.
pv_validate_fit_for_access <- function(x) {
  validate_pvstackr_fit(x, tier = "cheap")
}

#' Access pvstackr Estimates
#'
#' Return the reportable fixed-effect estimate table for a fit, or the aligned
#' estimate table for a method comparison. The returned columns depend on the
#' method contract but keep interval and provenance fields intact.
#'
#' @param x A pvstackr object.
#' @param ... Reserved for future extensions.
#'
#' @returns A data frame of reportable fixed-effect estimates. Fit estimate
#'   tables include interval/provenance columns such as `df_method`,
#'   `df_complete`, `interval_role`, `coverage_claim_allowed`,
#'   `target_source`, `target_hash`, `pooling_source`, and `pooling_hash` when
#'   those columns are part of the method contract. An inspection-only legacy
#'   PSIS object fails explicitly instead of returning historical numeric
#'   output.
#'
#' @examples
#' path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit
#'   head(get_estimates(fit))
#' }
#' @family pvstackr-accessors
#' @seealso [get_target()], [get_draws()], [get_diagnostics()]; [pv_fit()],
#'   [pv_compare_methods()].
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

#' Access pvstackr Targets
#'
#' Return the formal target object carried by a fit. Estimate-row provenance
#' labels such as `target_source` are separate from this accessor and do not
#' imply that every method has a formal target object.
#'
#' @param x A pvstackr fit or target object.
#' @param ... Reserved for future extensions.
#'
#' @returns The target object used by a fit, or `NULL` when the method has no
#'   target component. Estimate-row `target_source` labels are provenance
#'   metadata and may be present even when a method, such as `stack_psis`, does
#'   not carry a formal target object. A legacy PSIS inspection object also
#'   returns `NULL`.
#'
#' @examples
#' path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit
#'   get_target(fit)
#' }
#' @family pvstackr-accessors
#' @seealso [get_estimates()], [get_draws()], [get_diagnostics()];
#'   [pv_brr_target()], [pv_fit()].
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

#' Access pvstackr Draws
#'
#' Return retained top-level reportable draws from a fit. Methods that do not
#' synthesize a single reportable draw matrix, or fits created with
#' `return_draws = FALSE`, return `NULL`.
#'
#' @details
#' This accessor returns only the synthesized top-level reportable draw matrix.
#' Per-PV reference draws (`per_pv`) and the PSIS fixed-effect proposal/weight
#' pair (`stack_psis`), when retained, are not surfaced here; they live in the
#' fit's diagnostics (see [get_diagnostics()]), not in this top-level
#' reportable-draw accessor.
#'
#' @param x A pvstackr fit object.
#' @param ... Reserved for future extensions.
#'
#' @returns The retained reportable draw matrix, or `NULL` when draws were not
#'   retained or the method does not synthesize top-level reportable draws.
#'   Per-PV reference draws and the PSIS proposal/weight pair, when retained,
#'   remain available in diagnostics rather than through this top-level
#'   reportable-draw accessor. An inspection-only legacy PSIS object fails
#'   explicitly instead of returning historical draws.
#'
#' @examples
#' path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit
#'   get_draws(fit)
#' }
#' @family pvstackr-accessors
#' @seealso [get_estimates()], [get_target()], [get_diagnostics()]; [pv_fit()].
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

#' Access pvstackr Diagnostics
#'
#' Return the structured diagnostics list stored on a fit or method comparison.
#' Diagnostics keep method-specific details that are intentionally not flattened
#' into the reportable estimate table.
#'
#' @details
#' Diagnostic keys are method-specific; inspect `names(get_diagnostics(x))`. A
#' current `stack_direct` fit carries top-level `sampler` and `sampler_gate`
#' records plus `preflight`, `stack_fit`, `stack_fit_warnings`, and `ccc`;
#' sampler-blocked fits retain only slim preflight/sampler/gate evidence and the
#' independently valid external target, rebuilt from an exact recursive
#' allowlist with no formula object in preflight and a safe formula environment
#' on the target snapshot. A legacy cached stack-direct fit may
#' predate the sampler keys. A `per_pv` fit carries `reference` and `pooling`;
#' a `stack_psis` fit
#' carries `psis` (status and Pareto-k), `pooling`, and `weighted`. A method
#' comparison instead carries comparison-level keys such as `agreement`,
#' `method_diagnostics`, `timing`, and `target_overlap` (the shared-provenance
#' summary).
#'
#' @param x A pvstackr object.
#' @param ... Reserved for future extensions.
#'
#' @returns A structured diagnostics list. For a legacy PSIS inspection object
#'   this is bounded Pareto-k evidence plus its redaction record; historical
#'   pooling, weights, estimates, and draws are never returned.
#'
#' @examples
#' path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit
#'   names(get_diagnostics(fit))
#' }
#' @family pvstackr-accessors
#' @seealso [get_estimates()], [get_target()], [get_draws()]; [pv_fit()],
#'   [pv_compare_methods()].
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
  payload <- serialize(
    projected,
    NULL,
    ascii = FALSE,
    xdr = TRUE,
    version = 2L
  )
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

#' Display Methods for pvstackr Fits
#'
#' Compact console `print()` and `summary()` for a [pvstackr_fit][pv_fit]
#' object. The `print()` method shows the method id, fit status, fixed-effect
#' count, target source, draw retention, diagnostics keys, and the interval
#' note (set when reportable intervals are descriptive rather than
#' coverage-claimable). The `summary()` method builds a structured
#' `summary.pvstackr_fit` object; its `print()` method adds the compact
#' estimate columns (`term`, `estimate`, `se`, `df`, `conf_low`, `conf_high`,
#' where present).
#'
#' @param object,x A `pvstackr_fit` object (or its summary, for the
#'   `print.summary.pvstackr_fit` method).
#' @param ... Ignored.
#'
#' @returns The `print` methods return their input invisibly. `summary()`
#'   returns a `summary.pvstackr_fit` list with fields:
#'   \describe{
#'     \item{`method`, `status`}{Method id and fit status.}
#'     \item{`n_terms`, `terms`}{Number and names of reportable fixed-effect
#'       terms.}
#'     \item{`has_target`, `target_source`}{Whether a formal target object is
#'       carried, and the estimate-row target-source provenance label.}
#'     \item{`has_draws`, `draw_dim`}{Whether reportable draws were retained,
#'       and their `c(nrow, ncol)` dimension (`c(0L, 0L)` when not retained).}
#'     \item{`diagnostic_keys`}{Names of the entries in the diagnostics list.}
#'     \item{`interval_note`}{Set when reportable intervals are descriptive
#'       rather than coverage-claimable; `NA` otherwise.}
#'     \item{`reason_codes`, `warnings`}{Status reason codes and any captured
#'       fit warnings.}
#'     \item{`estimates`, `diagnostics`}{The reportable fixed-effect estimate
#'       table and the structured diagnostics list.}
#'     \item{`summary_schema_version`, `source_validation`,
#'       `source_reportability_fit`, `validation`}{The current summary schema,
#'       deep-valid compact source fit and its stamp, and the owned-summary
#'       SHA-256 record used before printing.}
#'   }
#' @seealso [pv_fit()], [get_estimates()], [get_diagnostics()].
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

#' Display Methods for pvstackr Method Comparisons
#'
#' Compact console `print()` and `summary()` for a
#' `pvstackr_method_comparison` object. The `print()` method shows the
#' reference method, the compared method labels, the fixed-effect count, any
#' blocked or warning methods, the provenance note (the shared-target caveat),
#' and the interval note. The `summary()` method builds a structured
#' `summary.pvstackr_method_comparison` object; its `print()` method also
#' renders the agreement table.
#'
#' @param object,x A `pvstackr_method_comparison` object (or its summary, for
#'   the `print.summary.pvstackr_method_comparison` method).
#' @param ... Ignored.
#'
#' @returns The `print` methods return their input invisibly. `summary()`
#'   returns a `summary.pvstackr_method_comparison` list with fields:
#'   \describe{
#'     \item{`reference_method`}{The method used as the comparison reference.}
#'     \item{`methods`, `method_labels`, `n_methods`}{The compared method ids,
#'       their display labels, and the number of compared methods.}
#'     \item{`n_terms`}{Number of distinct reportable fixed-effect terms.}
#'     \item{`blocked_methods`, `warning_methods`}{Methods whose reportable
#'       output was blocked, and methods that fit with warnings.}
#'     \item{`interval_note`}{Set when reportable intervals are descriptive
#'       rather than coverage-claimable; `NA` otherwise.}
#'     \item{`provenance_note`}{Shared-provenance caveat (close agreement is not
#'       independent corroboration when compared methods share a target source,
#'       pooling source, or estimand construction); `NA` otherwise.}
#'     \item{`estimate_table`, `diagnostic_table`}{The aligned cross-method
#'       estimate table and the per-method diagnostic table.}
#'     \item{`agreement`}{The descriptive cross-method agreement table.}
#'     \item{`timing`}{Per-method timing metadata.}
#'     \item{`summary_schema_version`, `source_validation`,
#'       `source_reportability_comparison`, `validation`}{The current summary
#'       schema, deep-valid compact source comparison and its stamp, and the
#'       owned-summary SHA-256 record used before printing.}
#'   }
#' @seealso [pv_compare_methods()], [get_estimates()], [get_diagnostics()].
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
