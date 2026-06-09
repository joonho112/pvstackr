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
#'   those columns are part of the method contract.
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
  validate_pvstackr_fit(x)
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
#'   not carry a formal target object.
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
  validate_pvstackr_fit(x)
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
#' Per-PV reference draws (`per_pv`) and PSIS importance weights (`stack_psis`),
#' when retained, are not surfaced here; they live in the fit's diagnostics
#' (see [get_diagnostics()]), not in this top-level reportable-draw accessor.
#'
#' @param x A pvstackr fit object.
#' @param ... Reserved for future extensions.
#'
#' @returns The retained reportable draw matrix, or `NULL` when draws were not
#'   retained or the method does not synthesize top-level reportable draws.
#'   Per-PV reference draws and PSIS weights, when retained, remain available in
#'   diagnostics rather than through this top-level reportable-draw accessor.
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
  validate_pvstackr_fit(x)
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
#' `stack_direct` fit carries `preflight`, `stack_fit`, `stack_fit_warnings`, and
#' `ccc`; a `per_pv` fit carries `reference` and `pooling`; a `stack_psis` fit
#' carries `psis` (status and Pareto-k), `pooling`, and `weighted`. A method
#' comparison instead carries comparison-level keys such as `agreement`,
#' `method_diagnostics`, `timing`, and `target_overlap` (the shared-provenance
#' summary).
#'
#' @param x A pvstackr object.
#' @param ... Reserved for future extensions.
#'
#' @returns A structured diagnostics list.
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
  validate_pvstackr_fit(x)
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

pv_fit_summary <- function(x) {
  estimates <- get_estimates(x)
  diagnostics <- get_diagnostics(x)
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
    schema_version = x$schema_version
  )
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
#'   }
#' @seealso [pv_fit()], [get_estimates()], [get_diagnostics()].
#' @name pvstackr_fit_summary
NULL

#' @rdname pvstackr_fit_summary
#' @export
print.pvstackr_fit <- function(x, ...) {
  validate_pvstackr_fit(x)
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
  summary
}

#' @rdname pvstackr_fit_summary
#' @export
print.summary.pvstackr_fit <- function(x, ...) {
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
  estimates <- get_estimates(x)
  diagnostics <- get_diagnostics(x)
  list(
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
  summary
}

#' @rdname pvstackr_method_comparison_summary
#' @export
print.summary.pvstackr_method_comparison <- function(x, ...) {
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
