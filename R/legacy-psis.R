pv_legacy_psis_inspection_schema <- function() {
  "pvstackr_legacy_psis_inspection_v1"
}

pv_legacy_psis_hard_threshold <- function() {
  0.7
}

pv_legacy_psis_plain_named_list <- function(x) {
  is.list(x) && !is.data.frame(x) && !is.object(x) &&
    !is.null(names(x)) && !anyNA(names(x)) && all(nzchar(names(x))) &&
    !anyDuplicated(names(x)) &&
    identical(attributes(x), list(names = names(x)))
}

pv_legacy_psis_scalar_character <- function(x, allowed, fallback = "unknown") {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !is.null(attributes(x)) || !x %in% allowed) {
    return(fallback)
  }
  unname(x)
}

pv_legacy_psis_extract_pareto <- function(x) {
  attrs <- attributes(x)
  valid <- is.numeric(x) && length(x) >= 2L && length(x) <= 1000L &&
    identical(names(attrs), "names") &&
    is.character(attrs$names) && length(attrs$names) == length(x) &&
    !anyNA(attrs$names) && all(nzchar(attrs$names)) &&
    all(nchar(attrs$names, type = "bytes") <= 128L) &&
    !anyDuplicated(attrs$names)
  if (!valid) {
    return(stats::setNames(numeric(), character()))
  }
  stats::setNames(as.numeric(x), attrs$names)
}

pv_legacy_psis_extract_threshold <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 ||
      x > 1 || !is.null(attributes(x))) {
    return(NA_real_)
  }
  unname(as.numeric(x))
}

pv_legacy_psis_validate_optional_threshold <- function(x, label) {
  missing <- is.numeric(x) && length(x) == 1L && is.na(x) &&
    is.null(attributes(x))
  if (missing) {
    return(NA_real_)
  }
  out <- pv_legacy_psis_extract_threshold(x)
  if (is.na(out)) {
    pv_abort(sprintf("Legacy PSIS inspection `%s` is noncanonical.", label))
  }
  out
}

pv_legacy_psis_redaction <- function() {
  list(
    status = "withheld",
    policy = "legacy_psis_inspection_only",
    withheld = c(
      "estimates", "draws", "design", "target", "stack_fit", "ccc",
      "pooling", "weighted", "weights", "proposal_draws", "backend_fit",
      "prepared_data", "log_lik"
    )
  )
}

pv_legacy_psis_inspection_provenance <- function() {
  list(
    function_name = "pv_migrate_legacy_psis_fit",
    package = "pvstackr",
    policy = "inspection_only_no_reportable_payload"
  )
}

pv_legacy_psis_inspection_diagnostics <- function(source) {
  diagnostics <- source[["diagnostics"]]
  if (!pv_legacy_psis_plain_named_list(diagnostics)) {
    diagnostics <- list()
  }
  psis <- diagnostics[["psis"]]
  if (!pv_legacy_psis_plain_named_list(psis)) {
    psis <- list()
  }
  pareto_k <- pv_legacy_psis_extract_pareto(psis[["pareto_k"]])
  diagnostic_threshold <- pv_legacy_psis_extract_threshold(psis[["threshold"]])
  control_threshold <- NA_real_
  control <- source[["control"]]
  if (is.list(control)) {
    control <- unclass(control)
    control_threshold <- pv_legacy_psis_extract_threshold(
      control[["psis_k_threshold"]]
    )
  }
  threshold_candidates <- c(
    pv_legacy_psis_hard_threshold(),
    diagnostic_threshold,
    control_threshold
  )
  effective_threshold <- min(threshold_candidates[!is.na(threshold_candidates)])
  complete <- length(pareto_k) >= 2L && all(is.finite(pareto_k))
  bad_pv_cols <- if (length(pareto_k) == 0L) {
    character()
  } else {
    names(pareto_k)[
      !is.finite(pareto_k) | pareto_k >= effective_threshold
    ]
  }
  evidence_status <- if (!complete) {
    "not_evaluated"
  } else if (length(bad_pv_cols) > 0L) {
    "failed"
  } else {
    "legacy_unsafe"
  }
  fallback_requested <- pv_legacy_psis_scalar_character(
    psis[["fallback_requested"]],
    c("block", "warn")
  )
  if (identical(fallback_requested, "unknown")) {
    fallback_requested <- pv_legacy_psis_scalar_character(
      psis[["fallback"]],
      c("block", "warn")
    )
  }

  list(
    psis = list(
      evidence_status = evidence_status,
      pareto_k = pareto_k,
      diagnostic_threshold = diagnostic_threshold,
      control_threshold = control_threshold,
      hard_threshold = pv_legacy_psis_hard_threshold(),
      effective_threshold = effective_threshold,
      pareto_k_max = if (complete) unname(max(pareto_k)) else NA_real_,
      bad_pv_cols = bad_pv_cols,
      complete = complete,
      fallback_requested = fallback_requested
    ),
    redaction = pv_legacy_psis_redaction()
  )
}

pv_validate_legacy_psis_inspection <- function(x) {
  fields <- c(
    "method", "inspection_only", "reportable", "reason_codes",
    "diagnostics", "source", "schema_version", "provenance", "warnings"
  )
  root_attrs <- attributes(x)
  if (!is.list(x) || !identical(names(x), fields) ||
      !identical(names(root_attrs), c("names", "class")) ||
      !identical(root_attrs$names, fields) ||
      !identical(
        root_attrs$class,
        c("pvstackr_legacy_psis_inspection", "list")
      )) {
    pv_abort("Legacy PSIS inspection fields, order, class, and attributes must be exact.")
  }
  if (!identical(x$method, "stack_psis") ||
      !identical(x$inspection_only, TRUE) ||
      !identical(x$reportable, FALSE) ||
      !identical(x$reason_codes, "legacy_psis_inspection_only") ||
      !identical(x$schema_version, pv_legacy_psis_inspection_schema()) ||
      !identical(x$warnings, paste(
        "Legacy stack_psis output is inspection-only; reportable estimates",
        "and draws were not migrated."
      )) ||
      !identical(x$provenance, pv_legacy_psis_inspection_provenance())) {
    pv_abort("Legacy PSIS inspection root metadata must use the exact non-reportable schema.")
  }
  source_fields <- c("object_class", "schema_version", "status")
  source_scalars <- if (is.list(x$source)) {
    vapply(
      x$source,
      function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
          is.null(attributes(value))
      },
      logical(1)
    )
  } else {
    logical()
  }
  if (!is.list(x$source) || !identical(names(x$source), source_fields) ||
      !identical(attributes(x$source), list(names = source_fields)) ||
      !identical(x$source$object_class, "pvstackr_fit") ||
      length(source_scalars) != length(source_fields) || any(!source_scalars) ||
      !x$source$schema_version %in% c("0.1.0", "0.2.0", "unknown") ||
      !x$source$status %in% c("ok", "warning", "blocked", "unknown")) {
    pv_abort("Legacy PSIS inspection source metadata must be canonical scalars.")
  }
  if (!is.list(x$diagnostics) ||
      !identical(names(x$diagnostics), c("psis", "redaction")) ||
      !identical(
        attributes(x$diagnostics),
        list(names = c("psis", "redaction"))
      ) || !identical(x$diagnostics$redaction, pv_legacy_psis_redaction())) {
    pv_abort("Legacy PSIS inspection diagnostics must use the exact slim schema.")
  }
  psis <- x$diagnostics$psis
  psis_fields <- c(
    "evidence_status", "pareto_k", "diagnostic_threshold",
    "control_threshold", "hard_threshold", "effective_threshold",
    "pareto_k_max", "bad_pv_cols", "complete", "fallback_requested"
  )
  canonical_evidence_status <- is.list(psis) &&
    is.character(psis$evidence_status) &&
    length(psis$evidence_status) == 1L && !is.na(psis$evidence_status) &&
    is.null(attributes(psis$evidence_status))
  canonical_fallback <- is.list(psis) &&
    is.character(psis$fallback_requested) &&
    length(psis$fallback_requested) == 1L &&
    !is.na(psis$fallback_requested) &&
    is.null(attributes(psis$fallback_requested))
  if (!is.list(psis) || !identical(names(psis), psis_fields) ||
      !identical(attributes(psis), list(names = psis_fields)) ||
      !canonical_evidence_status || !canonical_fallback ||
      !psis$evidence_status %in% c("failed", "not_evaluated", "legacy_unsafe") ||
      !identical(psis$hard_threshold, pv_legacy_psis_hard_threshold()) ||
      !is.logical(psis$complete) || length(psis$complete) != 1L ||
      is.na(psis$complete) || !is.null(attributes(psis$complete)) ||
      !psis$fallback_requested %in% c("block", "warn", "unknown")) {
    pv_abort("Legacy PSIS inspection evidence fields are noncanonical.")
  }
  pareto_k <- pv_legacy_psis_extract_pareto(psis$pareto_k)
  if (!identical(psis$pareto_k, pareto_k)) {
    pv_abort("Legacy PSIS inspection Pareto-k evidence is noncanonical.")
  }
  diagnostic_threshold <- pv_legacy_psis_validate_optional_threshold(
    psis$diagnostic_threshold,
    "diagnostic_threshold"
  )
  control_threshold <- pv_legacy_psis_validate_optional_threshold(
    psis$control_threshold,
    "control_threshold"
  )
  threshold_candidates <- c(
    pv_legacy_psis_hard_threshold(),
    diagnostic_threshold,
    control_threshold
  )
  expected_effective_threshold <- min(
    threshold_candidates[!is.na(threshold_candidates)]
  )
  expected_complete <- length(pareto_k) >= 2L && all(is.finite(pareto_k))
  expected_bad <- if (length(pareto_k) == 0L) character() else names(pareto_k)[
    !is.finite(pareto_k) | pareto_k >= expected_effective_threshold
  ]
  expected_status <- if (!expected_complete) {
    "not_evaluated"
  } else if (length(expected_bad) > 0L) {
    "failed"
  } else {
    "legacy_unsafe"
  }
  expected_max <- if (expected_complete) unname(max(pareto_k)) else NA_real_
  if (!identical(psis$complete, expected_complete) ||
      !identical(psis$bad_pv_cols, expected_bad) ||
      !identical(psis$evidence_status, expected_status) ||
      !identical(psis$effective_threshold, expected_effective_threshold) ||
      !identical(psis$pareto_k_max, expected_max) ||
      any(!vapply(
        psis[c(
          "evidence_status", "diagnostic_threshold", "control_threshold",
          "hard_threshold", "effective_threshold", "pareto_k_max",
          "bad_pv_cols", "fallback_requested"
        )],
        function(value) is.null(attributes(value)),
        logical(1)
      ))) {
    pv_abort("Legacy PSIS inspection evidence must reproduce its Pareto-k decision.")
  }
  invisible(x)
}

#' Convert a Legacy PSIS Fit to a Safe Inspection Object
#'
#' `pv_migrate_legacy_psis_fit()` returns a current, fully validated
#' `stack_psis` fit unchanged. A non-current or invalid historical
#' `stack_psis` fit is never promoted to reportable status: it is projected to
#' a separately allocated inspection-only object containing only bounded
#' Pareto-k decision evidence and an explicit redaction record. The source
#' object is not modified.
#'
#' Inspection objects allow [get_diagnostics()] and compact print/summary
#' methods. [get_estimates()] and [get_draws()] fail explicitly, because an
#' historical warning, blocked, unevaluated, incomplete, or `k >= 0.7` result
#' cannot be grandfathered into current reportable output. When historical
#' diagnostic and control thresholds differ, the effective inspection gate is
#' the minimum of both valid declarations and the immutable `0.7` ceiling.
#'
#' @param fit A current or historical `pvstackr_fit` with
#'   `method = "stack_psis"`.
#' @param x,object A `pvstackr_legacy_psis_inspection` object, or its summary
#'   for the summary print method.
#' @param ... Ignored.
#'
#' @return The unchanged current fit when full current validation succeeds;
#'   otherwise a `pvstackr_legacy_psis_inspection` object with no reportable
#'   estimates, draws, pooling, weights, backend, or data payload.
#' @export
pv_migrate_legacy_psis_fit <- function(fit) {
  if (inherits(fit, "pvstackr_legacy_psis_inspection")) {
    pv_validate_legacy_psis_inspection(fit)
    return(fit)
  }
  if (!inherits(fit, "pvstackr_fit") || !is.list(fit)) {
    pv_abort("`fit` must be a current or historical pvstackr_fit object.")
  }
  source <- unclass(fit)
  if (!identical(source[["method"]], "stack_psis")) {
    pv_abort("`pv_migrate_legacy_psis_fit()` accepts only stack_psis fits.")
  }
  is_current <- tryCatch(
    {
      validate_pvstackr_fit(fit, tier = "deep")
      TRUE
    },
    error = function(error) FALSE
  )
  if (is_current) {
    return(fit)
  }

  out <- list(
    method = "stack_psis",
    inspection_only = TRUE,
    reportable = FALSE,
    reason_codes = "legacy_psis_inspection_only",
    diagnostics = pv_legacy_psis_inspection_diagnostics(source),
    source = list(
      object_class = "pvstackr_fit",
      schema_version = pv_legacy_psis_scalar_character(
        source[["schema_version"]],
        c("0.1.0", "0.2.0")
      ),
      status = pv_legacy_psis_scalar_character(
        source[["status"]],
        c("ok", "warning", "blocked")
      )
    ),
    schema_version = pv_legacy_psis_inspection_schema(),
    provenance = pv_legacy_psis_inspection_provenance(),
    warnings = paste(
      "Legacy stack_psis output is inspection-only; reportable estimates",
      "and draws were not migrated."
    )
  )
  class(out) <- c("pvstackr_legacy_psis_inspection", "list")
  pv_validate_legacy_psis_inspection(out)
  out
}

#' @rdname get_estimates
#' @export
get_estimates.pvstackr_legacy_psis_inspection <- function(x, ...) {
  pv_validate_legacy_psis_inspection(x)
  pv_abort("Inspection-only legacy PSIS objects cannot expose reportable estimates.")
}

#' @rdname get_draws
#' @export
get_draws.pvstackr_legacy_psis_inspection <- function(x, ...) {
  pv_validate_legacy_psis_inspection(x)
  pv_abort("Inspection-only legacy PSIS objects cannot expose reportable draws.")
}

#' @rdname get_target
#' @export
get_target.pvstackr_legacy_psis_inspection <- function(x, ...) {
  pv_validate_legacy_psis_inspection(x)
  NULL
}

#' @rdname get_diagnostics
#' @export
get_diagnostics.pvstackr_legacy_psis_inspection <- function(x, ...) {
  pv_validate_legacy_psis_inspection(x)
  x$diagnostics
}

#' @rdname pv_migrate_legacy_psis_fit
#' @export
print.pvstackr_legacy_psis_inspection <- function(x, ...) {
  pv_validate_legacy_psis_inspection(x)
  cat("pvstackr legacy PSIS inspection\n")
  cat("  status: inspection_only\n")
  cat("  source status: ", x$source$status, "\n", sep = "")
  cat(
    "  PSIS evidence: ", x$diagnostics$psis$evidence_status, "\n",
    sep = ""
  )
  cat("  reportable estimates/draws: withheld\n")
  invisible(x)
}

#' @rdname pv_migrate_legacy_psis_fit
#' @export
summary.pvstackr_legacy_psis_inspection <- function(object, ...) {
  pv_validate_legacy_psis_inspection(object)
  out <- list(
    method = object$method,
    status = "inspection_only",
    source_status = object$source$status,
    evidence_status = object$diagnostics$psis$evidence_status,
    bad_pv_cols = object$diagnostics$psis$bad_pv_cols,
    reportable = FALSE,
    diagnostic_keys = names(object$diagnostics),
    schema_version = object$schema_version,
    source_inspection = object
  )
  class(out) <- c("summary.pvstackr_legacy_psis_inspection", "list")
  out
}

#' @rdname pv_migrate_legacy_psis_fit
#' @export
print.summary.pvstackr_legacy_psis_inspection <- function(x, ...) {
  fields <- c(
    "method", "status", "source_status", "evidence_status",
    "bad_pv_cols", "reportable", "diagnostic_keys", "schema_version",
    "source_inspection"
  )
  attrs <- attributes(x)
  if (!is.list(x) || !identical(names(x), fields) ||
      !identical(names(attrs), c("names", "class")) ||
      !identical(attrs$names, fields) ||
      !identical(attrs$class, c(
        "summary.pvstackr_legacy_psis_inspection", "list"
      )) || !identical(x$method, "stack_psis") ||
      !identical(x$status, "inspection_only") ||
      !identical(x$reportable, FALSE) ||
      !is.character(x$source_status) || length(x$source_status) != 1L ||
      is.na(x$source_status) ||
      !x$source_status %in% c("ok", "warning", "blocked", "unknown") ||
      !is.character(x$evidence_status) || length(x$evidence_status) != 1L ||
      is.na(x$evidence_status) ||
      !x$evidence_status %in% c("failed", "not_evaluated", "legacy_unsafe") ||
      !is.character(x$bad_pv_cols) || anyNA(x$bad_pv_cols) ||
      !identical(x$diagnostic_keys, c("psis", "redaction")) ||
      !identical(x$schema_version, pv_legacy_psis_inspection_schema())) {
    pv_abort("Legacy PSIS inspection summary is noncanonical.")
  }
  pv_validate_legacy_psis_inspection(x$source_inspection)
  if (!identical(x$method, x$source_inspection$method) ||
      !identical(x$source_status, x$source_inspection$source$status) ||
      !identical(
        x$evidence_status,
        x$source_inspection$diagnostics$psis$evidence_status
      ) || !identical(
        x$bad_pv_cols,
        x$source_inspection$diagnostics$psis$bad_pv_cols
      ) || !identical(
        x$diagnostic_keys,
        names(x$source_inspection$diagnostics)
      ) || !identical(x$schema_version, x$source_inspection$schema_version)) {
    pv_abort("Legacy PSIS inspection summary must match its validated source inspection.")
  }
  cat("pvstackr legacy PSIS inspection summary\n")
  cat("  source status: ", x$source_status, "\n", sep = "")
  cat("  PSIS evidence: ", x$evidence_status, "\n", sep = "")
  cat("  reportable estimates/draws: withheld\n")
  invisible(x)
}
