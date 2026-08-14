# Normalized sampler diagnostics ---------------------------------------------

pv_sampler_diagnostic_required_fields <- function() {
  c(
    "rhat_max",
    "ess_bulk_min",
    "ess_tail_min",
    "ess_bulk_per_chain_min",
    "ess_tail_per_chain_min",
    "divergences",
    "chains",
    "post_warmup_draws_per_chain",
    "diagnostic_source",
    "diagnostic_complete"
  )
}

pv_sampler_diagnostic_allowed_sources <- function() {
  c(
    "bundled_brms_posterior_and_nuts",
    "injected_diagnose_function",
    "injected_diagnose_function_absent"
  )
}

pv_sampler_diagnostic_reason_allowed <- function(reason) {
  fixed <- c(
    "diagnostic_incomplete_unspecified",
    "diagnostic_payload_not_list",
    "diagnostic_extractor_not_supplied",
    "diagnostic_extraction_failed",
    "diagnostic_counts_from_control_not_observed",
    "diagnostic_per_chain_mismatch",
    "posterior_namespace_unavailable",
    "brms_namespace_unavailable"
  )
  fields <- paste(c(
    "rhat_max", "ess_bulk_min", "ess_tail_min",
    "ess_bulk_per_chain_min", "ess_tail_per_chain_min",
    "divergences", "chains", "post_warmup_draws_per_chain"
  ), collapse = "|")
  reason %in% fixed || grepl(
    paste0(
      "^diagnostic_(missing|nonfinite|noninteger|nonpositive|negative)_(",
      fields,
      ")$"
    ),
    reason
  )
}

pv_sampler_diagnostics_incomplete <- function(
  diagnostic_source,
  reason_codes,
  chains = NA_integer_,
  post_warmup_draws_per_chain = NA_integer_
) {
  diagnostic_source <- pv_assert_scalar_string(
    diagnostic_source,
    "diagnostic_source"
  )
  reason_codes <- unique(as.character(reason_codes))
  reason_codes <- reason_codes[!is.na(reason_codes) & nzchar(reason_codes)]
  if (length(reason_codes) == 0L) {
    reason_codes <- "diagnostic_incomplete_unspecified"
  }
  list(
    rhat_max = NA_real_,
    ess_bulk_min = NA_real_,
    ess_tail_min = NA_real_,
    ess_bulk_per_chain_min = NA_real_,
    ess_tail_per_chain_min = NA_real_,
    divergences = NA_integer_,
    chains = as.integer(chains),
    post_warmup_draws_per_chain = as.integer(post_warmup_draws_per_chain),
    diagnostic_source = diagnostic_source,
    diagnostic_complete = FALSE,
    diagnostic_reason_codes = reason_codes
  )
}

pv_sampler_scalar_number <- function(payload, field, integerish = FALSE,
                                     positive = FALSE, nonnegative = FALSE) {
  value <- payload[[field]]
  reason <- character()
  if (is.null(value)) {
    return(list(value = if (integerish) NA_integer_ else NA_real_,
                reason = paste0("diagnostic_missing_", field)))
  }
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    return(list(value = if (integerish) NA_integer_ else NA_real_,
                reason = paste0("diagnostic_nonfinite_", field)))
  }
  if (integerish && value != floor(value)) {
    reason <- paste0("diagnostic_noninteger_", field)
  }
  if (positive && value <= 0) {
    reason <- paste0("diagnostic_nonpositive_", field)
  }
  if (nonnegative && value < 0) {
    reason <- paste0("diagnostic_negative_", field)
  }
  list(
    value = if (length(reason) > 0L) {
      if (integerish) NA_integer_ else NA_real_
    } else if (integerish) {
      as.integer(value)
    } else {
      as.numeric(value)
    },
    reason = reason
  )
}

pv_sampler_diagnostics_normalize <- function(
  payload,
  diagnostic_source,
  default_chains = NA_integer_,
  default_post_warmup_draws_per_chain = NA_integer_
) {
  diagnostic_source <- pv_assert_scalar_string(
    diagnostic_source,
    "diagnostic_source"
  )
  if (!is.list(payload)) {
    return(pv_sampler_diagnostics_incomplete(
      diagnostic_source,
      "diagnostic_payload_not_list",
      chains = default_chains,
      post_warmup_draws_per_chain = default_post_warmup_draws_per_chain
    ))
  }

  default_reasons <- character()
  if (is.null(payload$chains) && length(default_chains) == 1L &&
      is.finite(default_chains)) {
    payload$chains <- default_chains
    default_reasons <- c(
      default_reasons,
      "diagnostic_counts_from_control_not_observed"
    )
  }
  if (is.null(payload$post_warmup_draws_per_chain) &&
      length(default_post_warmup_draws_per_chain) == 1L &&
      is.finite(default_post_warmup_draws_per_chain)) {
    payload$post_warmup_draws_per_chain <- default_post_warmup_draws_per_chain
    default_reasons <- c(
      default_reasons,
      "diagnostic_counts_from_control_not_observed"
    )
  }

  specs <- list(
    rhat_max = list(positive = TRUE),
    ess_bulk_min = list(positive = TRUE),
    ess_tail_min = list(positive = TRUE),
    divergences = list(integerish = TRUE, nonnegative = TRUE),
    chains = list(integerish = TRUE, positive = TRUE),
    post_warmup_draws_per_chain = list(integerish = TRUE, positive = TRUE)
  )
  parsed <- lapply(names(specs), function(field) {
    args <- c(list(payload = payload, field = field), specs[[field]])
    do.call(pv_sampler_scalar_number, args)
  })
  names(parsed) <- names(specs)
  reasons <- c(
    default_reasons,
    unlist(lapply(parsed, `[[`, "reason"), use.names = FALSE)
  )

  chains <- parsed$chains$value
  ess_bulk_per_chain <- if (is.finite(parsed$ess_bulk_min$value) &&
                            is.finite(chains) && chains > 0L) {
    parsed$ess_bulk_min$value / chains
  } else {
    NA_real_
  }
  ess_tail_per_chain <- if (is.finite(parsed$ess_tail_min$value) &&
                            is.finite(chains) && chains > 0L) {
    parsed$ess_tail_min$value / chains
  } else {
    NA_real_
  }

  for (field in c("ess_bulk_per_chain_min", "ess_tail_per_chain_min")) {
    if (!is.null(payload[[field]])) {
      supplied <- pv_sampler_scalar_number(
        payload,
        field,
        positive = TRUE
      )
      expected <- if (identical(field, "ess_bulk_per_chain_min")) {
        ess_bulk_per_chain
      } else {
        ess_tail_per_chain
      }
      reasons <- c(reasons, supplied$reason)
      if (is.finite(supplied$value) && is.finite(expected) &&
          !isTRUE(all.equal(supplied$value, expected, tolerance = 1e-10))) {
        reasons <- c(reasons, "diagnostic_per_chain_mismatch")
      }
    }
  }
  reasons <- unique(reasons[nzchar(reasons)])

  out <- list(
    rhat_max = parsed$rhat_max$value,
    ess_bulk_min = parsed$ess_bulk_min$value,
    ess_tail_min = parsed$ess_tail_min$value,
    ess_bulk_per_chain_min = ess_bulk_per_chain,
    ess_tail_per_chain_min = ess_tail_per_chain,
    divergences = parsed$divergences$value,
    chains = parsed$chains$value,
    post_warmup_draws_per_chain = parsed$post_warmup_draws_per_chain$value,
    diagnostic_source = diagnostic_source,
    diagnostic_complete = length(reasons) == 0L,
    diagnostic_reason_codes = reasons
  )
  pv_validate_sampler_diagnostics(out)
  out
}

pv_validate_sampler_diagnostics <- function(diagnostics) {
  pv_assert_named_list(diagnostics, "sampler diagnostics")
  exact_fields <- c(
    pv_sampler_diagnostic_required_fields(),
    "diagnostic_reason_codes"
  )
  missing <- setdiff(exact_fields, names(diagnostics))
  if (length(missing) > 0L) {
    pv_abort(sprintf(
      "Sampler diagnostics are missing required field(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }
  extra <- setdiff(names(diagnostics), exact_fields)
  if (length(extra) > 0L) {
    pv_abort(sprintf(
      "Sampler diagnostics contain unsupported field(s): %s.",
      paste(extra, collapse = ", ")
    ))
  }
  pv_assert_scalar_string(
    diagnostics$diagnostic_source,
    "diagnostic_source"
  )
  if (!diagnostics$diagnostic_source %in%
      pv_sampler_diagnostic_allowed_sources()) {
    pv_abort("Sampler `diagnostic_source` is not recognized.")
  }
  if (!is.logical(diagnostics$diagnostic_complete) ||
      length(diagnostics$diagnostic_complete) != 1L ||
      is.na(diagnostics$diagnostic_complete)) {
    pv_abort("Sampler `diagnostic_complete` must be a non-missing logical scalar.")
  }
  pv_validate_character_field(
    diagnostics$diagnostic_reason_codes,
    "diagnostic_reason_codes"
  )
  if (length(diagnostics$diagnostic_reason_codes) > 0L &&
      any(!vapply(
        diagnostics$diagnostic_reason_codes,
        pv_sampler_diagnostic_reason_allowed,
        logical(1)
      ))) {
    pv_abort("Sampler diagnostic reason code is not recognized.")
  }

  numeric_fields <- c(
    "rhat_max", "ess_bulk_min", "ess_tail_min",
    "ess_bulk_per_chain_min", "ess_tail_per_chain_min",
    "divergences", "chains", "post_warmup_draws_per_chain"
  )
  for (field in numeric_fields) {
    value <- diagnostics[[field]]
    if (!is.numeric(value) || length(value) != 1L) {
      pv_abort(sprintf("Sampler `%s` must be a numeric scalar.", field))
    }
  }
  finite_positive <- c(
    "rhat_max", "ess_bulk_min", "ess_tail_min",
    "ess_bulk_per_chain_min", "ess_tail_per_chain_min", "chains",
    "post_warmup_draws_per_chain"
  )
  for (field in finite_positive) {
    value <- diagnostics[[field]]
    if (!is.na(value) && (!is.finite(value) || value <= 0)) {
      pv_abort(sprintf(
        "Sampler `%s` must be missing or a positive finite scalar.",
        field
      ))
    }
  }
  for (field in c("divergences", "chains", "post_warmup_draws_per_chain")) {
    value <- diagnostics[[field]]
    if (!is.na(value) &&
        (!is.finite(value) || value != floor(value) || value < 0 ||
         (field != "divergences" && value == 0))) {
      pv_abort(sprintf(
        "Sampler `%s` must be missing or a valid integer count.",
        field
      ))
    }
  }
  if (is.finite(diagnostics$chains) && diagnostics$chains > 0 &&
      is.finite(diagnostics$ess_bulk_min) &&
      is.finite(diagnostics$ess_bulk_per_chain_min) &&
      !isTRUE(all.equal(
        diagnostics$ess_bulk_per_chain_min,
        diagnostics$ess_bulk_min / diagnostics$chains,
        tolerance = 1e-10
      ))) {
    pv_abort("Sampler bulk per-chain ESS minimum is inconsistent.")
  }
  if (is.finite(diagnostics$chains) && diagnostics$chains > 0 &&
      is.finite(diagnostics$ess_tail_min) &&
      is.finite(diagnostics$ess_tail_per_chain_min) &&
      !isTRUE(all.equal(
        diagnostics$ess_tail_per_chain_min,
        diagnostics$ess_tail_min / diagnostics$chains,
        tolerance = 1e-10
      ))) {
    pv_abort("Sampler tail per-chain ESS minimum is inconsistent.")
  }
  if (isTRUE(diagnostics$diagnostic_complete)) {
    if (length(diagnostics$diagnostic_reason_codes) > 0L ||
        any(!is.finite(unlist(diagnostics[numeric_fields], use.names = FALSE))) ||
        diagnostics$rhat_max <= 0 || diagnostics$ess_bulk_min <= 0 ||
        diagnostics$ess_tail_min <= 0 || diagnostics$divergences < 0 ||
        diagnostics$chains <= 0 ||
        diagnostics$post_warmup_draws_per_chain <= 0 ||
        any(unlist(diagnostics[c(
          "divergences", "chains", "post_warmup_draws_per_chain"
        )], use.names = FALSE) != floor(unlist(diagnostics[c(
          "divergences", "chains", "post_warmup_draws_per_chain"
        )], use.names = FALSE)))) {
      pv_abort("Complete sampler diagnostics contain invalid values or reason codes.")
    }
    if (!isTRUE(all.equal(
      diagnostics$ess_bulk_per_chain_min,
      diagnostics$ess_bulk_min / diagnostics$chains,
      tolerance = 1e-10
    )) || !isTRUE(all.equal(
      diagnostics$ess_tail_per_chain_min,
      diagnostics$ess_tail_min / diagnostics$chains,
      tolerance = 1e-10
    ))) {
      pv_abort("Sampler per-chain ESS minima must equal total minima divided by chains.")
    }
  } else if (length(diagnostics$diagnostic_reason_codes) == 0L) {
    pv_abort("Incomplete sampler diagnostics require reason codes.")
  }
  invisible(diagnostics)
}

pv_stack_sampler_diagnostics <- function(
  fit,
  bundled_backend,
  diagnose_function,
  chains,
  iter,
  warmup,
  bundled_sampler_function = pv_backend_brms_sampler_diagnostics
) {
  post_warmup <- as.integer(iter - warmup)
  if (isTRUE(bundled_backend)) {
    sampler <- bundled_sampler_function(fit)
    if (is.null(diagnose_function)) {
      return(list(sampler = sampler))
    }
    custom_failed <- FALSE
    raw <- tryCatch(
      diagnose_function(fit),
      error = function(e) {
        custom_failed <<- TRUE
        NULL
      }
    )
    if (custom_failed) {
      return(list(
        custom_diagnostic_reason_codes = "diagnostic_extraction_failed",
        sampler = sampler
      ))
    }
    if (is.list(raw) && !is.null(raw$sampler)) {
      return(list(
        custom_sampler_override_ignored = TRUE,
        sampler = sampler
      ))
    }
    return(list(sampler = sampler))
  }
  if (is.null(diagnose_function)) {
    return(list(
      sampler = pv_sampler_diagnostics_incomplete(
        "injected_diagnose_function_absent",
        "diagnostic_extractor_not_supplied",
        chains = chains,
        post_warmup_draws_per_chain = post_warmup
      )
    ))
  }

  failed <- NULL
  raw <- tryCatch(
    diagnose_function(fit),
    error = function(e) {
      failed <<- TRUE
      NULL
    }
  )
  if (isTRUE(failed)) {
    return(list(
      sampler = pv_sampler_diagnostics_incomplete(
        "injected_diagnose_function",
        "diagnostic_extraction_failed",
        chains = chains,
        post_warmup_draws_per_chain = post_warmup
      )
    ))
  }

  payload <- if (is.list(raw) && !is.null(raw$sampler)) raw$sampler else raw
  sampler <- pv_sampler_diagnostics_normalize(
    payload,
    diagnostic_source = "injected_diagnose_function",
    default_chains = chains,
    default_post_warmup_draws_per_chain = post_warmup
  )
  list(sampler = sampler)
}

pv_sampler_gate_thresholds <- function() {
  list(
    rhat_warning_above = 1.01,
    rhat_blocked_above = 1.05,
    ess_total_blocked_below = 100,
    ess_per_chain_warning_below = 100,
    ess_per_chain_blocked_below = 25,
    divergences_required = 0L
  )
}

pv_sampler_gate_reason_codes <- function() {
  c(
    "sampler_diagnostics_incomplete",
    "sampler_configuration_blocked",
    "sampler_rhat_warning",
    "sampler_rhat_blocked",
    "sampler_ess_bulk_warning",
    "sampler_ess_bulk_blocked",
    "sampler_ess_tail_warning",
    "sampler_ess_tail_blocked",
    "sampler_divergences_blocked"
  )
}

pv_sampler_gate <- function(diagnostics, expected_chains,
                            expected_post_warmup_draws_per_chain) {
  pv_validate_sampler_diagnostics(diagnostics)
  expected_chains <- pv_sampler_scalar_number(
    list(value = expected_chains),
    "value",
    integerish = TRUE,
    positive = TRUE
  )$value
  expected_post_warmup_draws_per_chain <- pv_sampler_scalar_number(
    list(value = expected_post_warmup_draws_per_chain),
    "value",
    integerish = TRUE,
    positive = TRUE
  )$value
  if (!is.finite(expected_chains) ||
      !is.finite(expected_post_warmup_draws_per_chain)) {
    pv_abort("Sampler gate expected draw configuration is invalid.")
  }

  thresholds <- pv_sampler_gate_thresholds()
  reason_codes <- character()
  warnings <- character()
  bands <- character()
  add_finding <- function(band, reason_code, warning) {
    bands <<- c(bands, band)
    reason_codes <<- c(reason_codes, reason_code)
    warnings <<- c(warnings, warning)
  }

  if (!isTRUE(diagnostics$diagnostic_complete)) {
    add_finding(
      "blocked",
      "sampler_diagnostics_incomplete",
      "Sampler diagnostics are incomplete; reportable stack_direct output is blocked."
    )
  } else {
    if (!identical(diagnostics$chains, as.integer(expected_chains)) ||
        !identical(
          diagnostics$post_warmup_draws_per_chain,
          as.integer(expected_post_warmup_draws_per_chain)
        )) {
      add_finding(
        "blocked",
        "sampler_configuration_blocked",
        "Observed sampler chain/draw counts do not match the requested fit configuration."
      )
    }
    if (diagnostics$rhat_max > thresholds$rhat_blocked_above) {
      add_finding(
        "blocked",
        "sampler_rhat_blocked",
        sprintf(
          "Sampler maximum R-hat %.6g exceeds the block threshold %.6g.",
          diagnostics$rhat_max,
          thresholds$rhat_blocked_above
        )
      )
    } else if (diagnostics$rhat_max > thresholds$rhat_warning_above) {
      add_finding(
        "warning",
        "sampler_rhat_warning",
        sprintf(
          "Sampler maximum R-hat %.6g exceeds the warning threshold %.6g.",
          diagnostics$rhat_max,
          thresholds$rhat_warning_above
        )
      )
    }

    for (prefix in c("bulk", "tail")) {
      total <- diagnostics[[paste0("ess_", prefix, "_min")]]
      per_chain <- diagnostics[[paste0(
        "ess_", prefix, "_per_chain_min"
      )]]
      if (total < thresholds$ess_total_blocked_below ||
          per_chain < thresholds$ess_per_chain_blocked_below) {
        add_finding(
          "blocked",
          paste0("sampler_ess_", prefix, "_blocked"),
          sprintf(
            "Sampler minimum %s ESS %.6g (%.6g per chain) violates the total/per-chain block floors %.6g/%.6g.",
            prefix,
            total,
            per_chain,
            thresholds$ess_total_blocked_below,
            thresholds$ess_per_chain_blocked_below
          )
        )
      } else if (per_chain < thresholds$ess_per_chain_warning_below) {
        add_finding(
          "warning",
          paste0("sampler_ess_", prefix, "_warning"),
          sprintf(
            "Sampler minimum %s ESS per chain %.6g is below the warning threshold %.6g.",
            prefix,
            per_chain,
            thresholds$ess_per_chain_warning_below
          )
        )
      }
    }
    if (diagnostics$divergences > thresholds$divergences_required) {
      add_finding(
        "blocked",
        "sampler_divergences_blocked",
        sprintf(
          "Sampler recorded %d divergence(s); reportable output requires zero.",
          diagnostics$divergences
        )
      )
    }
  }

  status <- if ("blocked" %in% bands) {
    "blocked"
  } else if ("warning" %in% bands) {
    "warning"
  } else {
    "ok"
  }
  out <- list(
    diagnostic_role = "reportability_gate",
    status = status,
    reason_codes = unique(reason_codes),
    warnings = unique(warnings),
    thresholds = thresholds
  )
  pv_validate_sampler_gate(out)
  out
}

pv_validate_sampler_gate <- function(gate) {
  pv_assert_named_list(gate, "sampler gate")
  exact <- c(
    "diagnostic_role", "status", "reason_codes", "warnings", "thresholds"
  )
  if (!setequal(names(gate), exact) || length(gate) != length(exact)) {
    pv_abort("Sampler gate must contain the exact frozen status fields.")
  }
  if (!identical(gate$diagnostic_role, "reportability_gate") ||
      !is.character(gate$status) || length(gate$status) != 1L ||
      is.na(gate$status) || !gate$status %in% c("ok", "warning", "blocked")) {
    pv_abort("Sampler gate role or status is invalid.")
  }
  pv_validate_character_field(gate$reason_codes, "sampler gate reason_codes")
  pv_validate_character_field(gate$warnings, "sampler gate warnings")
  if (any(!gate$reason_codes %in% pv_sampler_gate_reason_codes())) {
    pv_abort("Sampler gate contains an unknown reason code.")
  }
  if (!identical(gate$thresholds, pv_sampler_gate_thresholds())) {
    pv_abort("Sampler gate thresholds do not match the frozen policy.")
  }
  if (identical(gate$status, "ok") &&
      (length(gate$reason_codes) > 0L || length(gate$warnings) > 0L)) {
    pv_abort("OK sampler gate must not contain reasons or warnings.")
  }
  if (!identical(gate$status, "ok") &&
      (length(gate$reason_codes) == 0L || length(gate$warnings) == 0L)) {
    pv_abort("Non-OK sampler gate requires reasons and warnings.")
  }
  invisible(gate)
}
