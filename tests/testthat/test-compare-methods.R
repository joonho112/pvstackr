compare_reference_draws <- function() {
  list(
    PV1 = cbind(
      b_Intercept = c(0.8, 0.9, 1.0, 1.1),
      b_x = c(0.2, 0.3, 0.4, 0.5),
      sigma = c(0.9, 1.0, 1.1, 1.2)
    ),
    PV2 = cbind(
      b_Intercept = c(1.0, 1.1, 1.2, 1.3),
      b_x = c(0.4, 0.5, 0.6, 0.7),
      sigma = c(1.0, 1.1, 1.2, 1.3)
    ),
    PV3 = cbind(
      b_Intercept = c(1.2, 1.3, 1.4, 1.5),
      b_x = c(0.6, 0.7, 0.8, 0.9),
      sigma = c(1.1, 1.2, 1.3, 1.4)
    )
  )
}

compare_stacked_draws <- function() {
  cbind(
    b_Intercept = c(0.9, 1.0, 1.1, 1.2, 1.3, 1.4),
    b_x = c(0.25, 0.35, 0.45, 0.55, 0.65, 0.75),
    sigma = c(0.9, 1.0, 1.1, 1.0, 1.2, 1.1)
  )
}

compare_weights <- function() {
  cbind(
    PV1 = c(1, 1, 1, 1, 1, 1),
    PV2 = c(1, 2, 1, 2, 1, 2),
    PV3 = c(2, 1, 2, 1, 2, 1)
  )
}

compare_control <- function(method, ...) {
  pv_control(method = method, backend = "injected", iter = 10L, warmup = 5L, ...)
}

compare_fit_reference <- function() {
  pv_fit_reference(
    per_pv_draws = compare_reference_draws(),
    control = compare_control("per_pv")
  )
}

compare_fit_psis <- function(pareto_k = c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4), fallback = "block", ...) {
  call_fit <- function() pv_fit_stack_psis(
    stacked_draws = compare_stacked_draws(),
    psis_weights = compare_weights(),
    pareto_k = pareto_k,
    psis_producer = "testthat_psis_fixture",
    psis_producer_version = "1.0.0",
    fallback = fallback,
    control = compare_control("stack_psis"),
    ...
  )
  if (identical(fallback, "warn")) suppressWarnings(call_fit()) else call_fit()
}

compare_fit_psis_one_term <- function() {
  pv_fit_stack_psis(
    stacked_draws = compare_stacked_draws()[, c("b_Intercept", "sigma"), drop = FALSE],
    psis_weights = compare_weights(),
    pareto_k = c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4),
    psis_producer = "testthat_psis_fixture",
    psis_producer_version = "1.0.0",
    control = compare_control("stack_psis")
  )
}

compare_direct_data <- function() {
  data.frame(
    id = seq_len(10),
    x = c(-2, -1.5, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3),
    PV1 = c(2.1, 2.4, 3.0, 3.2, 3.6, 4.0, 4.2, 4.8, 5.1, 5.5),
    PV2 = c(2.0, 2.5, 2.9, 3.4, 3.5, 4.1, 4.4, 4.7, 5.0, 5.8),
    W = c(1.0, 1.1, 0.9, 1.2, 1.0, 1.3, 0.8, 1.4, 1.1, 0.95),
    RW1 = c(0.8, 1.3, 1.0, 1.1, 0.9, 1.5, 0.9, 1.2, 1.2, 1.0),
    RW2 = c(1.2, 0.9, 0.8, 1.4, 1.1, 1.0, 0.7, 1.6, 1.0, 0.9),
    RW3 = c(1.1, 1.0, 1.1, 1.0, 1.3, 1.2, 0.8, 1.3, 1.4, 1.1),
    RW4 = c(0.9, 1.2, 0.95, 1.3, 1.0, 1.4, 1.1, 1.1, 0.9, 1.2)
  )
}

compare_fit_direct <- function(...) {
  data <- compare_direct_data()
  target <- pv_brr_target(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    ...
  )
  fit_function <- function(...) {
    beta <- target$beta[target$fe_names]
    se <- sqrt(diag(target$T_MI))[target$fe_names]
    e1 <- c(-0.03, -0.01, 0.02, 0.04, -0.02, 0.01, 0.03, -0.04)
    e2 <- c(0.02, -0.03, 0.01, -0.04, 0.03, -0.01, 0.04, -0.02)
    list(draws = cbind(
      b_Intercept = beta[["b_Intercept"]] + se[["b_Intercept"]] * e1,
      b_x = beta[["b_x"]] + se[["b_x"]] * e2,
      sigma = rep(1, 8)
    ))
  }
  pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = pv_control(
      method = "stack_direct",
      backend = "injected",
      iter = 10L,
      warmup = 5L,
      return_draws = FALSE
    ),
    fit_function = fit_function,
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(
      chains = 4L,
      post_warmup = 5L
    )
  )
}

test_that("pv_compare_methods builds a comparison table and diagnostics", {
  per_pv <- compare_fit_reference()
  psis <- compare_fit_psis()

  comparison <- pv_compare_methods(
    reference = per_pv,
    psis = psis,
    timings = c(reference = 3.2, psis = 1.1),
    include_fits = TRUE
  )

  expect_s3_class(comparison, "pvstackr_method_comparison")
  expect_equal(comparison$reference_method, "reference")
  expect_equal(comparison$method_labels, c("reference", "psis"))
  expect_equal(nrow(comparison$table), 4L)
  expect_equal(unique(comparison$table$term), c("b_Intercept", "b_x"))
  expect_equal(unique(comparison$table$method), c("per_pv", "stack_psis"))
  expect_equal(comparison$table$std.error, comparison$table$se)
  expect_true(all(c(
    "df_method", "df_complete", "conf_level", "interval_role",
    "coverage_claim_allowed"
  ) %in% names(comparison$table)))
  expect_identical(comparison$estimate_table, comparison$table)
  expect_identical(comparison$agreement, comparison$diagnostics$agreement)
  expect_identical(comparison$timing, comparison$diagnostics$timing)
  expect_identical(comparison$diagnostic_table, comparison$diagnostics$method_diagnostics)
  ref_rows <- comparison$table$method_label == "reference"
  psis_rows <- comparison$table$method_label == "psis"
  expect_equal(comparison$table$estimate_diff[ref_rows], c(0, 0), tolerance = 1e-14)
  expect_equal(comparison$table$se_ratio[ref_rows], c(1, 1), tolerance = 1e-14)
  expect_equal(comparison$table$agreement_band[ref_rows], c("close", "close"))
  expect_equal(comparison$table$df_method[ref_rows], per_pv$estimates$df_method)
  expect_equal(comparison$table$df_complete[ref_rows], per_pv$estimates$df_complete)
  expect_equal(comparison$table$conf_level[ref_rows], per_pv$estimates$conf_level)
  expect_equal(comparison$table$interval_role[ref_rows], per_pv$estimates$interval_role)
  expect_equal(comparison$table$coverage_claim_allowed[ref_rows], per_pv$estimates$coverage_claim_allowed)
  expect_equal(comparison$table$df_method[psis_rows], psis$estimates$df_method)
  expect_equal(comparison$table$interval_role[psis_rows], psis$estimates$interval_role)
  expect_equal(comparison$table$coverage_claim_allowed[psis_rows], psis$estimates$coverage_claim_allowed)
  expect_true(all(c(
    "df_method", "interval_role", "coverage_claim_allowed",
    "n_descriptive_intervals", "target_source", "target_hash",
    "pooling_source", "pooling_hash", "center", "shared_target_hash",
    "shared_pooling_hash", "shared_external_target",
    "shares_reference_target", "shares_reference_pooling", "weight_method",
    "weight_ess_iid_min", "weight_ess_fraction_min",
    "max_normalized_weight_max"
  ) %in% names(comparison$diagnostic_table)))
  expect_equal(comparison$diagnostic_table$df_method, c("classic", "classic"))
  expect_equal(
    comparison$diagnostic_table$interval_role,
    c("reference_classic_rubin", "psis_classic_rubin")
  )
  expect_equal(comparison$diagnostic_table$coverage_claim_allowed, c(FALSE, FALSE))
  expect_equal(comparison$diagnostic_table$n_descriptive_intervals, c(2L, 2L))
  expect_equal(
    comparison$diagnostic_table$target_source,
    c("per_pv_rubin_draws", "stack_psis_rubin_pooling")
  )
  expect_equal(
    comparison$diagnostic_table$pooling_source,
    c("per_pv_rubin_draws", "stack_psis_rubin_pooling")
  )
  expect_false(any(comparison$diagnostic_table$shared_target_hash))
  expect_false(any(comparison$diagnostic_table$shared_pooling_hash))
  expect_false(any(comparison$diagnostic_table$shared_external_target))
  expect_false(any(comparison$diagnostic_table$shares_reference_target))
  expect_false(any(comparison$diagnostic_table$shares_reference_pooling))
  expect_true(all(is.na(comparison$diagnostic_table$center)))
  expect_identical(
    comparison$diagnostic_table$weight_method,
    c(NA_character_, "caller_declared_external_psis")
  )
  expect_identical(
    comparison$diagnostic_table$psis_source,
    c(NA_character_, "supplied_psis_weights")
  )
  expect_identical(
    comparison$diagnostic_table$pareto_k_source,
    c(NA_character_, "supplied")
  )
  expect_identical(
    comparison$diagnostic_table$psis_producer,
    c(NA_character_, "testthat_psis_fixture")
  )
  expect_identical(
    comparison$diagnostic_table$psis_producer_version,
    c(NA_character_, "1.0.0")
  )
  psis_rows <- comparison$table$method == "stack_psis"
  expect_identical(
    comparison$table$weight_method[psis_rows],
    rep("caller_declared_external_psis", sum(psis_rows))
  )
  expect_identical(
    comparison$table$psis_source[psis_rows],
    rep("supplied_psis_weights", sum(psis_rows))
  )
  expect_identical(
    comparison$table$psis_producer[psis_rows],
    rep("testthat_psis_fixture", sum(psis_rows))
  )
  expect_equal(
    comparison$diagnostic_table$weight_ess_iid_min,
    c(NA_real_, 5.4),
    tolerance = 1e-14
  )
  expect_equal(
    comparison$diagnostic_table$weight_ess_fraction_min,
    c(NA_real_, 0.9),
    tolerance = 1e-14
  )
  expect_equal(
    comparison$diagnostic_table$max_normalized_weight_max,
    c(NA_real_, 2 / 9),
    tolerance = 1e-14
  )
  expect_false(comparison$diagnostics$target_overlap$shared_target_hash)
  expect_false(comparison$diagnostics$target_overlap$shared_pooling_hash)
  expect_false(comparison$diagnostics$target_overlap$shared_external_target)
  expect_false(comparison$diagnostics$target_overlap$independence_caveat_required)
  expect_match(comparison$diagnostics$target_overlap$independence_caveat, "independent corroboration", fixed = TRUE)
  expect_equal(comparison$diagnostics$timing$elapsed_seconds, c(3.2, 1.1))
  expect_equal(comparison$diagnostics$timing$n_fits, c(3L, 1L))
  expect_s3_class(comparison$fits$reference, "pvstackr_fit")
  expect_identical(
    comparison$source_fit_projection$psis$psis_source,
    "supplied_psis_weights"
  )
  expect_equal(
    comparison$source_fit_projection$psis$weight_ess_iid_min,
    5.4,
    tolerance = 1e-14
  )
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(comparison))
})

test_that("pv_compare_methods uses implemented method IDs and stable auto labels", {
  per_pv <- compare_fit_reference()
  psis <- compare_fit_psis()
  direct <- compare_fit_direct()

  comparison <- pv_compare_methods(per_pv, psis, direct)

  expect_equal(comparison$method_labels, c("per_pv", "stack_psis", "stack_direct"))
  expect_equal(comparison$methods, c(per_pv = "per_pv", stack_psis = "stack_psis", stack_direct = "stack_direct"))
  expect_setequal(unique(comparison$table$method), c("per_pv", "stack_psis", "stack_direct"))
  expect_equal(nrow(comparison$table), length(comparison$method_labels) * length(unique(comparison$table$term)))
  expect_equal(
    comparison$diagnostic_table$target_source,
    c("per_pv_rubin_draws", "stack_psis_rubin_pooling", "external_brr_fay_rubin")
  )
  expect_equal(comparison$diagnostic_table$center, c(NA_character_, NA_character_, "target"))
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(comparison))

  duplicates <- pv_compare_methods(per_pv, compare_fit_reference())
  expect_equal(duplicates$method_labels, c("per_pv", "per_pv_1"))
  expect_true(all(duplicates$diagnostic_table$shared_target_hash))
  expect_true(all(duplicates$diagnostic_table$shared_pooling_hash))
  expect_false(any(duplicates$diagnostic_table$shared_external_target))
  expect_equal(duplicates$diagnostic_table$shares_reference_target, c(FALSE, TRUE))
  expect_equal(duplicates$diagnostic_table$shares_reference_pooling, c(FALSE, TRUE))
  expect_true(duplicates$diagnostics$target_overlap$shared_target_hash)
  expect_true(duplicates$diagnostics$target_overlap$shared_pooling_hash)
  expect_false(duplicates$diagnostics$target_overlap$shared_external_target)
  expect_equal(duplicates$diagnostics$target_overlap$shared_target_sources, "per_pv_rubin_draws")
  expect_true(duplicates$diagnostics$target_overlap$independence_caveat_required)
})

test_that("pv_compare_methods preserves mixed interval coverage semantics", {
  per_pv <- pv_fit_reference(
    per_pv_draws = compare_reference_draws(),
    control = compare_control("per_pv"),
    df_method = "barnard_rubin",
    df_complete = c(b_x = 40, b_Intercept = 30)
  )
  psis <- compare_fit_psis(
    df_method = "barnard_rubin",
    df_complete = c(b_x = 45, b_Intercept = 35)
  )
  direct <- compare_fit_direct(
    df_method = "barnard_rubin",
    df_complete = c(b_x = 50, b_Intercept = 40)
  )

  comparison <- pv_compare_methods(per_pv = per_pv, psis = psis, direct = direct)

  expect_equal(
    comparison$diagnostic_table$interval_role,
    c("reference_barnard_rubin", "psis_barnard_rubin", "coverage_barnard_rubin")
  )
  expect_equal(
    comparison$diagnostic_table$coverage_claim_allowed,
    c(FALSE, FALSE, TRUE)
  )
  expect_equal(
    comparison$diagnostic_table$n_descriptive_intervals,
    c(2L, 2L, 0L)
  )
  expect_equal(pvstackr:::pv_interval_note(comparison$estimate_table), "some intervals are descriptive rather than coverage-claimable.")
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(comparison))
})

test_that("pv_compare_methods flags shared external target provenance", {
  direct_a <- compare_fit_direct()
  direct_b <- compare_fit_direct()

  comparison <- pv_compare_methods(direct_a = direct_a, direct_b = direct_b)

  expect_equal(comparison$diagnostic_table$target_source, rep("external_brr_fay_rubin", 2L))
  expect_equal(comparison$diagnostic_table$target_hash, rep(direct_a$target$target_hash, 2L))
  expect_true(all(comparison$diagnostic_table$shared_target_hash))
  expect_true(all(comparison$diagnostic_table$shared_external_target))
  expect_false(any(comparison$diagnostic_table$shared_pooling_hash))
  expect_equal(comparison$diagnostic_table$shares_reference_target, c(FALSE, TRUE))
  expect_false(any(comparison$diagnostic_table$shares_reference_pooling))
  expect_equal(comparison$diagnostic_table$center, rep("target", 2L))
  expect_true(comparison$diagnostics$target_overlap$shared_target_hash)
  expect_true(comparison$diagnostics$target_overlap$shared_external_target)
  expect_false(comparison$diagnostics$target_overlap$shared_pooling_hash)
  expect_true(comparison$diagnostics$target_overlap$shares_reference_target)
  expect_false(comparison$diagnostics$target_overlap$shares_reference_pooling)
  expect_equal(comparison$diagnostics$target_overlap$shared_target_sources, "external_brr_fay_rubin")
  expect_true(comparison$diagnostics$target_overlap$independence_caveat_required)
  expect_match(
    comparison$diagnostics$target_overlap$independence_caveat,
    "independent corroboration",
    fixed = TRUE
  )
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(comparison))
})

test_that("pv_compare_methods records blocked methods without dropping terms", {
  per_pv <- compare_fit_reference()
  blocked <- compare_fit_psis(c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.4))

  comparison <- pv_compare_methods(list(per_pv = per_pv, psis = blocked))

  expect_equal(comparison$diagnostics$blocked_methods, "psis")
  psis_rows <- comparison$table$method_label == "psis"
  expect_equal(sum(psis_rows), 2L)
  expect_true(all(is.na(comparison$table$estimate[psis_rows])))
  expect_true(all(is.na(comparison$table$estimate_diff[psis_rows])))
  expect_true(all(is.na(comparison$table$df_method[psis_rows])))
  expect_true(all(is.na(comparison$table$df_complete[psis_rows])))
  expect_true(all(is.na(comparison$table$conf_level[psis_rows])))
  expect_true(all(is.na(comparison$table$interval_role[psis_rows])))
  expect_true(all(is.na(comparison$table$coverage_claim_allowed[psis_rows])))
  expect_equal(comparison$table$agreement_band[psis_rows], rep("not_available", 2L))
  expect_equal(comparison$table$reason_codes[psis_rows], rep("psis_k_too_high", 2L))
  psis_diagnostics <- comparison$diagnostic_table[comparison$diagnostic_table$method_label == "psis", ]
  expect_true(is.na(psis_diagnostics$df_method))
  expect_true(is.na(psis_diagnostics$interval_role))
  expect_true(is.na(psis_diagnostics$coverage_claim_allowed))
  expect_equal(psis_diagnostics$n_descriptive_intervals, 0L)
  expect_true(is.na(psis_diagnostics$target_source))
  expect_true(is.na(psis_diagnostics$target_hash))
  expect_true(is.na(psis_diagnostics$pooling_source))
  expect_true(is.na(psis_diagnostics$pooling_hash))
  expect_false(psis_diagnostics$shared_target_hash)
  expect_false(psis_diagnostics$shared_pooling_hash)
  expect_false(psis_diagnostics$shared_external_target)
  expect_false(psis_diagnostics$shares_reference_target)
  expect_false(psis_diagnostics$shares_reference_pooling)
  expect_true(is.na(psis_diagnostics$center))
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(comparison))
})

test_that("legacy PSIS comparisons and summaries cannot expose derived numerics", {
  reference <- compare_fit_reference()
  psis <- compare_fit_psis()
  current <- pv_compare_methods(
    reference = reference,
    psis = psis,
    include_fits = TRUE
  )

  expect_identical(current$schema_version, "0.2.0")
  expect_identical(
    current$validation$schema_version,
    "pvstackr_method_comparison_validation_v1"
  )
  expect_identical(
    names(current$source_fit_validation),
    current$method_labels
  )
  expect_identical(
    names(current$source_fit_projection),
    current$method_labels
  )
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(current))
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(
    unserialize(serialize(current, NULL, version = 2L))
  ))

  without_fits <- pv_compare_methods(reference = reference, psis = psis)
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(
    unserialize(serialize(without_fits, NULL, version = 2L))
  ))

  legacy <- current
  legacy$source_fit_validation <- NULL
  legacy$validation <- NULL
  legacy$schema_version <- "0.1.0"
  psis_diag <- legacy$diagnostic_table$method_label == "psis"
  psis_rows <- legacy$estimate_table$method_label == "psis"
  legacy$diagnostic_table$status[psis_diag] <- "warning"
  legacy$diagnostic_table$reason_codes[psis_diag] <- "psis_k_too_high"
  legacy$diagnostic_table$psis_status[psis_diag] <- "warning"
  legacy$diagnostic_table$pareto_k_max[psis_diag] <- 0.91
  legacy$diagnostics$method_diagnostics <- legacy$diagnostic_table
  legacy$diagnostics$statuses[["psis"]] <- "warning"
  legacy$diagnostics$warning_methods <- "psis"
  legacy$table$status[psis_rows] <- "warning"
  legacy$table$reason_codes[psis_rows] <- "psis_k_too_high"
  legacy$estimate_table <- legacy$table

  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(legacy),
    "inspection-only"
  )
  expect_error(get_estimates(legacy), "inspection-only")
  expect_error(get_diagnostics(legacy), "inspection-only")
  expect_error(print(legacy), "inspection-only")
  expect_error(summary(legacy), "inspection-only")

  fit_summary <- summary(psis)
  expect_identical(fit_summary$summary_schema_version, "0.2.0")
  expect_false(fit_summary$source_reportability_fit$control$return_draws)
  expect_null(
    fit_summary$source_reportability_fit$diagnostics$weighted$proposal_draws
  )
  expect_null(fit_summary$source_reportability_fit$diagnostics$weighted$weights)
  expect_output(print(fit_summary), "pvstackr fit summary")
  expect_output(
    print(unserialize(serialize(fit_summary, NULL, version = 2L))),
    "pvstackr fit summary"
  )
  legacy_fit_summary <- fit_summary
  legacy_fit_summary$summary_schema_version <- NULL
  legacy_fit_summary$source_validation <- NULL
  legacy_fit_summary$source_reportability_fit <- NULL
  legacy_fit_summary$validation <- NULL
  legacy_fit_summary$status <- "warning"
  legacy_fit_summary$diagnostics$psis$status <- "warning"
  legacy_fit_summary$diagnostics$psis$pareto_k <-
    c(PV1 = 0.2, PV2 = 0.91, PV3 = 0.4)
  legacy_fit_summary$diagnostics$psis$pareto_k_max <- 0.91
  expect_error(print(legacy_fit_summary), "inspection-only")
  malformed_legacy_fit_summary <- legacy_fit_summary
  malformed_legacy_fit_summary$method <- c("stack_psis", "other")
  expect_error(
    print(malformed_legacy_fit_summary),
    "schema is unrecognized"
  )

  comparison_summary <- summary(current)
  expect_identical(comparison_summary$summary_schema_version, "0.2.0")
  expect_null(
    current$source_fit_reportability$psis$diagnostics$weighted$proposal_draws
  )
  expect_null(
    current$source_fit_reportability$psis$diagnostics$weighted$weights
  )
  expect_output(
    print(comparison_summary),
    "pvstackr method comparison summary"
  )
  expect_output(
    print(unserialize(serialize(comparison_summary, NULL, version = 2L))),
    "pvstackr method comparison summary"
  )
  legacy_comparison_summary <- comparison_summary
  legacy_comparison_summary$summary_schema_version <- NULL
  legacy_comparison_summary$source_validation <- NULL
  legacy_comparison_summary$source_reportability_comparison <- NULL
  legacy_comparison_summary$validation <- NULL
  expect_error(print(legacy_comparison_summary), "inspection-only")
  malformed_legacy_comparison_summary <- legacy_comparison_summary
  malformed_legacy_comparison_summary$methods <- as.list(
    malformed_legacy_comparison_summary$methods
  )
  expect_error(
    print(malformed_legacy_comparison_summary),
    "schema is unrecognized"
  )

  missing_marker <- current
  missing_marker$validation <- NULL
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(missing_marker),
    "inspection-only"
  )

  bad_stamp <- current
  bad_stamp$validation$stamp <- paste0("sha256:", strrep("f", 64L))
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad_stamp),
    "validation stamp"
  )

  bad_summary_stamp <- fit_summary
  bad_summary_stamp$validation$stamp <- paste0("sha256:", strrep("f", 64L))
  expect_error(print(bad_summary_stamp), "validation")

  restamped_estimate_summary <- fit_summary
  restamped_estimate_summary$estimates$estimate[[1L]] <- 987654321
  restamped_estimate_summary <- pvstackr:::pv_summary_issue_validation_stamp(
    restamped_estimate_summary,
    "fit"
  )
  expect_error(print(restamped_estimate_summary), "source reportability")

  restamped_warning_summary <- fit_summary
  restamped_warning_summary$status <- "warning"
  restamped_warning_summary$diagnostics$psis$status <- "warning"
  restamped_warning_summary$diagnostics$psis$pareto_k <-
    c(PV1 = 0.2, PV2 = 0.91, PV3 = 0.4)
  restamped_warning_summary$diagnostics$psis$pareto_k_max <- 0.91
  restamped_warning_summary <- pvstackr:::pv_summary_issue_validation_stamp(
    restamped_warning_summary,
    "fit"
  )
  expect_error(print(restamped_warning_summary), "source reportability")

  blocked_fit_summary <- summary(compare_fit_psis(
    c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.4)
  ))
  promoted_blocked_summary <- blocked_fit_summary
  visible_fields <- setdiff(
    names(fit_summary),
    c("source_validation", "source_reportability_fit", "validation")
  )
  promoted_blocked_summary[visible_fields] <- fit_summary[visible_fields]
  promoted_blocked_summary <- pvstackr:::pv_summary_issue_validation_stamp(
    promoted_blocked_summary,
    "fit"
  )
  expect_error(print(promoted_blocked_summary), "source reportability")

  restamped_comparison_summary <- comparison_summary
  restamped_comparison_summary$estimate_table$estimate[
    restamped_comparison_summary$estimate_table$method == "stack_psis"
  ] <- 123456789
  restamped_comparison_summary <- pvstackr:::pv_summary_issue_validation_stamp(
    restamped_comparison_summary,
    "comparison"
  )
  expect_error(
    print(restamped_comparison_summary),
    "validated source comparison"
  )
})

test_that("blocked PSIS comparison metadata remain redacted after restamping", {
  comparison <- pv_compare_methods(
    reference = compare_fit_reference(),
    psis = compare_fit_psis(c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.4))
  )
  psis_diag <- comparison$diagnostic_table$method_label == "psis"
  comparison$diagnostic_table$pooling_source[psis_diag] <-
    "legacy_psis_pooling"
  comparison$diagnostic_table$pooling_hash[psis_diag] <-
    "legacy_pooling_hash"
  comparison$diagnostics$method_diagnostics <- comparison$diagnostic_table
  comparison <- pvstackr:::pv_comparison_issue_validation_stamp(comparison)
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(comparison),
    "retains reportable or pooling metadata"
  )
})

test_that("comparison derivation remains bound to source fit projections", {
  reference <- compare_fit_reference()
  safe_psis <- compare_fit_psis()
  blocked_psis <- compare_fit_psis(c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.4))

  for (include_fits in c(FALSE, TRUE)) {
    safe <- pv_compare_methods(
      reference = reference,
      psis = safe_psis,
      include_fits = include_fits
    )
    blocked <- pv_compare_methods(
      reference = reference,
      psis = blocked_psis,
      include_fits = include_fits
    )
    substituted <- blocked
    derived_fields <- c(
      "table", "estimate_table", "diagnostics", "diagnostic_table",
      "agreement", "timing"
    )
    substituted[derived_fields] <- safe[derived_fields]
    substituted <- pvstackr:::pv_comparison_issue_validation_stamp(substituted)
    expect_error(
      pvstackr:::validate_pvstackr_method_comparison(substituted),
      "source projections",
      info = paste("include_fits", include_fits)
    )
    expect_error(
      get_estimates(substituted),
      "source projections",
      info = paste("include_fits", include_fits)
    )
  }
})

test_that("pv_compare_methods tolerates hollow blocked stack_direct fits", {
  per_pv <- compare_fit_reference()
  blocked_direct <- pvstackr:::new_pvstackr_fit(
    "stack_direct",
    status = "blocked",
    reason_codes = "preflight_failed",
    control = pv_control(method = "stack_direct", return_draws = FALSE)
  )

  comparison <- pv_compare_methods(list(per_pv = per_pv, direct = blocked_direct))

  direct_rows <- comparison$table$method_label == "direct"
  expect_equal(sum(direct_rows), 2L)
  expect_true(all(is.na(comparison$table$estimate[direct_rows])))
  expect_equal(comparison$table$agreement_band[direct_rows], rep("not_available", 2L))
  expect_equal(comparison$table$reason_codes[direct_rows], rep("preflight_failed", 2L))
  direct_diagnostics <- comparison$diagnostic_table[comparison$diagnostic_table$method_label == "direct", ]
  expect_equal(direct_diagnostics$status, "blocked")
  expect_true(is.na(direct_diagnostics$target_source))
  expect_true(is.na(direct_diagnostics$target_hash))
  expect_true(is.na(direct_diagnostics$pooling_source))
  expect_true(is.na(direct_diagnostics$pooling_hash))
  expect_equal(direct_diagnostics$center, "target")
  expect_false(direct_diagnostics$shared_external_target)
  expect_false(direct_diagnostics$shares_reference_target)
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(comparison))
})

test_that("pv_compare_methods aligns missing terms and computes agreement diagnostics", {
  per_pv <- compare_fit_reference()
  one_term <- compare_fit_psis_one_term()
  blocked <- compare_fit_psis(c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.4), fallback = "warn")

  comparison <- pv_compare_methods(
    fits = list(reference = per_pv, one_term = one_term, blocked = blocked),
    timings = list(one_term = 0.8, stack_psis = 1.1)
  )

  missing <- comparison$table$method_label == "one_term" & comparison$table$term == "b_x"
  expect_true(any(missing))
  expect_true(is.na(comparison$table$estimate[missing]))
  expect_true(is.na(comparison$table$estimate_diff[missing]))
  expect_equal(
    comparison$agreement$n_available[comparison$agreement$method_label == "one_term"],
    1L
  )
  expect_equal(comparison$diagnostics$warning_methods, character())
  expect_equal(comparison$diagnostics$blocked_methods, "blocked")
  expect_equal(
    comparison$timing$elapsed_seconds,
    c(NA_real_, 0.8, 1.1)
  )

  row <- comparison$table[comparison$table$method_label == "one_term" & comparison$table$term == "b_Intercept", ]
  expect_equal(row$estimate_diff, row$estimate - row$reference_estimate)
  expect_equal(row$se_ratio, row$se / row$reference_se)
  expect_equal(row$abs_z_diff, abs(row$estimate_diff) / sqrt(row$se^2 + row$reference_se^2))
  expect_equal(
    comparison$agreement$max_abs_z_diff[comparison$agreement$method_label == "one_term"],
    row$abs_z_diff
  )
  expect_invisible(pvstackr:::validate_pvstackr_method_comparison(comparison))
})

test_that("pv_compare_methods validates references and malformed inputs", {
  per_pv <- compare_fit_reference()
  psis <- compare_fit_psis()
  blocked <- compare_fit_psis(c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.4))

  explicit <- pv_compare_methods(
    fits = list(reference = per_pv, psis = psis),
    reference_method = "per_pv"
  )
  expect_equal(explicit$reference_method, "reference")

  expect_error(
    pv_compare_methods(per_pv),
    "at least two"
  )
  expect_error(
    pv_compare_methods(per_pv, blocked, reference_method = "stack_psis"),
    "blocked"
  )
  expect_error(
    pv_compare_methods(per_pv, list()),
    "pvstackr_fit"
  )
  expect_error(
    pv_compare_methods(per_pv, psis, timings = c(1.0, 2.0)),
    "named"
  )
  expect_error(
    pv_compare_methods(per_pv, psis, timings = c(per_pv = -1.0)),
    "non-negative"
  )

  bad <- explicit
  bad$table <- data.frame()
  bad$estimate_table <- bad$table
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "non-empty"
  )

  bad <- explicit
  bad$table$method[[1L]] <- "pipeline_a"
  bad$estimate_table <- bad$table
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "implemented method IDs"
  )

  bad <- explicit
  bad$table$std.error[[1L]] <- bad$table$std.error[[1L]] + 1
  bad$estimate_table <- bad$table
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "std.error"
  )

  bad <- explicit
  bad$table$coverage_claim_allowed <- as.character(bad$table$coverage_claim_allowed)
  bad$estimate_table <- bad$table
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "coverage_claim_allowed"
  )

  bad <- explicit
  bad$diagnostic_table$n_descriptive_intervals[[1L]] <- -1
  bad$diagnostics$method_diagnostics <- bad$diagnostic_table
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "n_descriptive_intervals"
  )

  bad <- explicit
  bad$diagnostic_table$shared_target_hash[[1L]] <- TRUE
  bad$diagnostics$method_diagnostics <- bad$diagnostic_table
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "shared-provenance"
  )

  bad <- explicit
  bad$diagnostics$target_overlap$shared_target_hash <- TRUE
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "target-overlap diagnostics"
  )

  bad <- explicit
  bad$diagnostics$agreement <- bad$diagnostics$agreement[0, ]
  bad$agreement <- bad$diagnostics$agreement
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "non-empty agreement"
  )

  bad <- explicit
  bad$table <- bad$table[-1L, ]
  bad$estimate_table <- bad$table
  expect_error(
    pvstackr:::validate_pvstackr_method_comparison(bad),
    "exactly one row"
  )
})

test_that("provenance-blocked PSIS fits remain non-reportable in comparisons", {
  per_pv <- compare_fit_reference()
  unverified <- pvstackr::pv_fit_stack_psis(
    stacked_draws = compare_stacked_draws(),
    psis_weights = compare_weights(),
    pareto_k = c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4),
    control = compare_control("stack_psis")
  )
  unsmoothed <- pvstackr::pv_fit_stack_psis(
    stacked_draws = compare_stacked_draws(),
    log_ratios = log(compare_weights()),
    pareto_k = c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4),
    control = compare_control("stack_psis")
  )

  for (fit in list(unverified = unverified, unsmoothed = unsmoothed)) {
    comparison <- pv_compare_methods(
      reference = per_pv,
      diagnostic = fit,
      reference_method = "reference"
    )
    row <- comparison$diagnostics$method_diagnostics[
      comparison$diagnostics$method_diagnostics$method_label == "diagnostic",
      ,
      drop = FALSE
    ]
    expect_equal(row$status, "blocked")
    expect_equal(row$n_available, 0L)
    expect_true(row$psis_status %in% c("provenance_incomplete", "unsmoothed"))
    expected_method <- if (identical(fit, unverified)) {
      "unspecified_external"
    } else {
      "self_normalized_raw_importance"
    }
    expect_identical(row$weight_method, expected_method)
    expect_equal(row$weight_ess_iid_min, 5.4, tolerance = 1e-14)
    expect_equal(row$weight_ess_fraction_min, 0.9, tolerance = 1e-14)
    expect_equal(row$max_normalized_weight_max, 2 / 9, tolerance = 1e-14)
    expect_true(all(is.na(
      comparison$table$estimate[
        comparison$table$method_label == "diagnostic"
      ]
    )))
    expect_invisible(
      pvstackr:::validate_pvstackr_method_comparison(comparison)
    )
  }
})
