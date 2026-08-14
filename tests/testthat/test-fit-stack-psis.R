stack_psis_fixture_data <- function() {
  data.frame(
    id = seq_len(5),
    x = c(-1, -0.5, 0, 0.5, 1),
    PV1 = c(1.0, 1.2, 1.5, 1.8, 2.0),
    PV2 = c(1.1, 1.3, 1.6, 1.7, 2.1),
    PV3 = c(0.9, 1.4, 1.5, 1.9, 2.2),
    W = c(1.0, 1.1, 0.9, 1.2, 1.0)
  )
}

stack_psis_draws_fixture <- function() {
  cbind(
    b_Intercept = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3),
    b_x = c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45),
    sigma = c(0.9, 1.0, 1.1, 1.0, 1.2, 1.1),
    lp__ = seq(-4, -3.5, length.out = 6)
  )
}

stack_psis_weights_fixture <- function() {
  cbind(
    PV1 = c(1, 1, 1, 1, 1, 1),
    PV2 = c(1, 2, 1, 2, 1, 2),
    PV3 = c(2, 1, 2, 1, 2, 1)
  )
}

stack_psis_control <- function(...) {
  pv_control(
    method = "stack_psis",
    backend = "injected",
    iter = 10L,
    warmup = 5L,
    chains = 2L,
    seed = 20260607L,
    ...
  )
}

# Existing reportable fixtures explicitly attest that their supplied or
# injected weights came from an external PSIS producer. Tests of missing
# provenance call pvstackr::pv_fit_stack_psis() directly.
pv_fit_stack_psis <- function(...) {
  args <- list(...)
  external_route <- !is.null(args$psis_weights) || !is.null(args$psis_function)
  if (external_route && !"psis_producer" %in% names(args)) {
    args$psis_producer <- "testthat_psis_fixture"
    args$psis_producer_version <- "1.0.0"
  }
  call_fit <- function() do.call(pvstackr::pv_fit_stack_psis, args)
  if (identical(args$fallback, "warn")) suppressWarnings(call_fit()) else call_fit()
}

stack_psis_oracle <- function(draws, weights) {
  weights <- sweep(weights, 2L, colSums(weights), FUN = "/")
  summary <- pvstackr:::pv_stack_psis_summarize(draws, weights)
  pool <- pvstackr:::rubin_pool_matrix(summary$beta, summary$U, orientation = "rows_pv")
  list(summary = summary, pool = pool)
}

test_that("pv_fit_stack_psis pools PSIS-weighted stacked draws when Pareto-k passes", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.4, PV3 = 0.6)
  oracle <- stack_psis_oracle(draws, weights)

  fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    control = stack_psis_control(return_draws = TRUE)
  )

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$method, "stack_psis")
  expect_equal(fit$status, "ok")
  expect_equal(fit$reason_codes, character())
  expect_null(fit$target)
  expect_null(fit$ccc)
  expect_null(fit$draws)
  expect_equal(fit$diagnostics$psis$status, "ok")
  expect_equal(fit$diagnostics$psis$source, "supplied_psis_weights")
  expect_equal(fit$diagnostics$psis$weight_method, "caller_declared_external_psis")
  expect_equal(fit$diagnostics$psis$producer, "testthat_psis_fixture")
  expect_equal(
    fit$diagnostics$psis$weight_diagnostic_authority,
    "retained_weights_recomputed"
  )
  expect_equal(fit$diagnostics$psis$n_draws, 6L)
  expect_equal(
    fit$diagnostics$psis$weight_ess_iid,
    c(PV1 = 6, PV2 = 5.4, PV3 = 5.4),
    tolerance = 1e-14
  )
  expect_equal(
    fit$diagnostics$psis$weight_ess_fraction,
    c(PV1 = 1, PV2 = 0.9, PV3 = 0.9),
    tolerance = 1e-14
  )
  expect_equal(
    fit$diagnostics$psis$max_normalized_weight,
    c(PV1 = 1 / 6, PV2 = 2 / 9, PV3 = 2 / 9),
    tolerance = 1e-14
  )
  expect_equal(fit$diagnostics$psis$pareto_k, pareto_k)
  expect_equal(fit$diagnostics$pooling$pooling_source, "stack_psis_rubin_pooling")
  expect_equal(
    fit$diagnostics$weighted$proposal_draws,
    draws[, c("b_Intercept", "b_x"), drop = FALSE],
    tolerance = 0
  )
  expect_equal(unname(colSums(fit$diagnostics$weighted$weights)), rep(1, 3), tolerance = 1e-14)
  expect_equal(fit$diagnostics$pooling$beta, oracle$pool$beta, tolerance = 1e-14)
  expect_equal(fit$diagnostics$pooling$T_MI, oracle$pool$T_MI, tolerance = 1e-14)
  expect_equal(fit$estimates$term, names(oracle$pool$beta))
  expect_equal(fit$estimates$estimate, unname(oracle$pool$beta), tolerance = 1e-14)
  expect_equal(fit$estimates$se, unname(oracle$pool$se), tolerance = 1e-14)
  expect_equal(fit$estimates$std.error, fit$estimates$se)
  expect_equal(fit$estimates$interval_role, rep("psis_classic_rubin", 2L))
  expect_false(any(fit$estimates$coverage_claim_allowed))
  expect_equal(fit$estimates$target_source, rep("stack_psis_rubin_pooling", 2L))
  expect_equal(fit$estimates$pooling_source, rep("stack_psis_rubin_pooling", 2L))
  expect_equal(fit$estimates$psis_status, rep("ok", 2L))
  expect_equal(fit$estimates$pareto_k_max, rep(max(pareto_k), 2L))
  expect_equal(fit$estimates$pooling_hash, rep(fit$diagnostics$pooling$pooling_hash, 2L))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  restamped_coverage <- fit
  restamped_coverage$estimates$coverage_claim_allowed <- TRUE
  restamped_coverage$estimates$interval_role <- "coverage_interval"
  restamped_coverage <- pvstackr:::pv_fit_issue_validation_stamp(
    restamped_coverage
  )
  expect_error(
    pvstackr:::validate_pvstackr_fit(restamped_coverage),
    "canonical pooled projection"
  )

  redacted <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    control = stack_psis_control(return_draws = FALSE)
  )
  expect_null(redacted$diagnostics$weighted$proposal_draws)
  expect_null(redacted$diagnostics$weighted$weights)
  expect_equal(
    redacted$diagnostics$psis$weight_diagnostic_authority,
    "owned_stamp_bounded_projection"
  )
  expect_identical(
    redacted$diagnostics$pooling$pooling_hash,
    fit$diagnostics$pooling$pooling_hash
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(redacted))

  injected <- redacted
  injected$diagnostics$weighted$proposal_draws <-
    fit$diagnostics$weighted$proposal_draws
  injected$diagnostics$weighted$weights <- fit$diagnostics$weighted$weights
  expect_error(pvstackr:::validate_pvstackr_fit(injected), "return_draws")

  removed <- fit
  removed$diagnostics$weighted["proposal_draws"] <- list(NULL)
  removed$diagnostics$weighted["weights"] <- list(NULL)
  expect_error(pvstackr:::validate_pvstackr_fit(removed), "return_draws")

  changed_draw <- fit
  changed_draw$diagnostics$weighted$proposal_draws[1L, 1L] <-
    changed_draw$diagnostics$weighted$proposal_draws[1L, 1L] + 1
  expect_error(
    pvstackr:::validate_pvstackr_fit(changed_draw),
    "reproduce weighted summaries"
  )

  changed_weight <- fit
  changed_weight$diagnostics$weighted$weights[, 1L] <-
    changed_weight$diagnostics$weighted$weights[, 1L] * 2
  expect_error(
    pvstackr:::validate_pvstackr_fit(changed_weight),
    "canonical PV-aligned matrix"
  )

  jointly_reordered <- fit
  reverse_rows <- nrow(jointly_reordered$diagnostics$weighted$proposal_draws):1L
  jointly_reordered$diagnostics$weighted$proposal_draws <-
    jointly_reordered$diagnostics$weighted$proposal_draws[reverse_rows, , drop = FALSE]
  jointly_reordered$diagnostics$weighted$weights <-
    jointly_reordered$diagnostics$weighted$weights[reverse_rows, , drop = FALSE]
  rownames(jointly_reordered$diagnostics$weighted$proposal_draws) <- NULL
  rownames(jointly_reordered$diagnostics$weighted$weights) <- NULL
  expect_error(
    pvstackr:::validate_pvstackr_fit(jointly_reordered),
    "validation stamp"
  )

  injected_summary <- redacted
  injected_summary$diagnostics$pooling$U_bar <- matrix(1, 2L, 2L)
  expect_error(pvstackr:::validate_pvstackr_fit(injected_summary))

  injected_summary <- redacted
  injected_summary$diagnostics$weighted$param_map$private_draws <-
    matrix(1, 2L, 2L)
  expect_error(pvstackr:::validate_pvstackr_fit(injected_summary))

  injected_summary <- redacted
  injected_summary$diagnostics$weighted$param_map$vc_names <-
    rep("PRIVATE_FALSE_ACCESSOR_SENTINEL", 100L)
  expect_error(pvstackr:::validate_pvstackr_fit(injected_summary))

  injected_summary <- redacted
  injected_summary$diagnostics$weighted$param_map$dropped_names <- "lp__"
  expect_error(pvstackr:::validate_pvstackr_fit(injected_summary))

  injected_summary <- redacted
  injected_summary$diagnostics$weighted$param_map$map_source <- "explicit"
  expect_error(pvstackr:::validate_pvstackr_fit(injected_summary))

  injected_summary <- redacted
  attr(injected_summary$diagnostics$psis$pareto_k, "private_draws") <-
    matrix(1, 2L, 2L)
  expect_error(pvstackr:::validate_pvstackr_fit(injected_summary))
})

test_that("stack_psis keep_data removes raw inputs and PSIS source labels", {
  data <- stack_psis_fixture_data()
  data$UNUSED_PSIS_PRIVATE_COLUMN <- rep(
    "PSIS_RAW_VALUE_SENTINEL",
    nrow(data)
  )
  formula_env <- new.env(parent = baseenv())
  formula_env$private_state <- "PSIS_FORMULA_ENV_SENTINEL"
  formula <- OUTCOME ~ x
  environment(formula) <- formula_env
  weights <- stack_psis_weights_fixture()
  rownames(weights) <- paste0("PSIS_WEIGHT_ROW_SENTINEL_", seq_len(nrow(weights)))
  attr(weights, "private_marker") <- "PSIS_WEIGHT_ATTR_SENTINEL"
  args <- list(
    data = data,
    formula = formula,
    pv_cols = c("PV1", "PV2", "PV3"),
    stacked_draws = stack_psis_draws_fixture(),
    psis_weights = weights,
    pareto_k = c(PV1 = 0.2, PV2 = 0.4, PV3 = 0.6),
    weight_col = "W",
    id_cols = "id"
  )
  redacted <- do.call(
    pv_fit_stack_psis,
    c(args, list(control = stack_psis_control(keep_data = FALSE)))
  )
  retained <- do.call(
    pv_fit_stack_psis,
    c(args, list(control = stack_psis_control(keep_data = TRUE)))
  )
  expect_null(redacted$design$data)
  expect_identical(environment(redacted$design$formula), baseenv())
  expect_identical(retained$design$data, data)
  expect_identical(environment(retained$design$formula), baseenv())
  expect_null(rownames(redacted$diagnostics$weighted$weights))
  expect_null(attr(redacted$diagnostics$weighted$weights, "private_marker"))
  serialized <- rawToChar(serialize(redacted, NULL, ascii = TRUE))
  for (marker in c(
    "UNUSED_PSIS_PRIVATE_COLUMN", "PSIS_RAW_VALUE_SENTINEL",
    "PSIS_FORMULA_ENV_SENTINEL", "PSIS_WEIGHT_ROW_SENTINEL",
    "PSIS_WEIGHT_ATTR_SENTINEL"
  )) {
    expect_false(grepl(marker, serialized, fixed = TRUE), info = marker)
  }
  expect_lt(
    length(serialize(redacted, NULL, xdr = TRUE)),
    length(serialize(retained, NULL, xdr = TRUE))
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(redacted))

  control_flip <- redacted
  control_flip$control$keep_data <- TRUE
  expect_error(pvstackr:::validate_pvstackr_fit(control_flip))
})

test_that("supplied stack_fit is projected without retention regression or payload injection", {
  data <- stack_psis_fixture_data()
  fit_function <- function(formula, data, ...) {
    list(draws = stack_psis_draws_fixture(), data = data)
  }
  stack_fit <- pvstackr:::pv_stack_fit(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2", "PV3"),
    weight_col = "W",
    control = pv_control(
      method = "stack_direct",
      backend = "injected",
      iter = 10L,
      warmup = 5L,
      chains = 2L
    ),
    fit_function = fit_function,
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )
  fit <- pv_fit_stack_psis(
    stack_fit = stack_fit,
    psis_weights = stack_psis_weights_fixture(),
    pareto_k = c(PV1 = 0.2, PV2 = 0.4, PV3 = 0.6),
    control = stack_psis_control()
  )
  expect_s3_class(fit, "pvstackr_fit")
  expect_null(fit$stack_fit$prepared_data)
  expect_null(fit$stack_fit$fit)
  expect_null(fit$stack_fit$stacked_draws)
  expect_false(fit$stack_fit$control$return_draws)
  expect_identical(environment(fit$stack_fit$formula), baseenv())
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  secret <- data.frame(secret = "EXTERNAL_STACK_FIT_SECRET")
  raw_secret <- serialize(secret, NULL)
  row_secret <- rep("EXTERNAL_STACK_FIT_SECRET", 5L)
  mutations <- list(
    root_field = function(x) { x$private_payload <- secret; x },
    root_attribute = function(x) { attr(x, "private_payload") <- secret; x },
    provenance = function(x) { x$provenance$private_payload <- secret; x },
    meta = function(x) { x$meta$private_payload <- secret; x },
    weight_summary = function(x) { x$weight_summary$private_payload <- secret; x },
    param_map = function(x) { x$param_map$private_payload <- secret; x },
    draws_attribute = function(x) {
      attr(x$stacked_draws, "private_payload") <- secret
      x
    },
    warnings_attribute = function(x) {
      attr(x$warnings, "private_payload") <- secret
      x
    },
    diagnostics_attribute = function(x) {
      attr(x$diagnostics, "private_payload") <- secret
      x
    },
    sampler_attribute = function(x) {
      attr(x$diagnostics$sampler, "private_payload") <- secret
      x
    },
    sampler_scalar_attribute = function(x) {
      attr(x$diagnostics$sampler$rhat_max, "private_payload") <- secret
      x
    },
    sampler_scalar_names = function(x) {
      names(x$diagnostics$sampler$rhat_max) <- "EXTERNAL_STACK_FIT_SECRET"
      x
    },
    engine_nested_extra = function(x) {
      x$provenance$engine$private_payload <- secret
      x
    },
    package_versions_extra = function(x) {
      x$provenance$engine$package_versions$private_payload <- secret
      x
    },
    cmdstan_state_extra = function(x) {
      x$provenance$engine$cmdstan_state$private_payload <- secret
      x
    },
    cache_nested_extra = function(x) {
      x$provenance$cache$private_payload <- secret
      x$meta$cache$private_payload <- secret
      x
    },
    vc_policy_extra = function(x) {
      x$meta$vc_policy$private_payload <- secret
      x
    },
    prior_policy_extra = function(x) {
      x$meta$prior_policy$private_payload <- secret
      x$meta$prior_diagnostic$private_payload <- secret
      x
    },
    weight_scalar_attribute = function(x) {
      attr(x$weight_summary$M, "private_payload") <- secret
      x
    },
    weight_scalar_names = function(x) {
      names(x$weight_summary$M) <- "EXTERNAL_STACK_FIT_SECRET"
      x
    },
    param_scalar_attribute = function(x) {
      attr(x$param_map$original_fe_idx, "private_payload") <- secret
      x
    },
    param_scalar_names = function(x) {
      names(x$param_map$map_source) <- "EXTERNAL_STACK_FIT_SECRET"
      x
    },
    psi_attribute = function(x) {
      attr(x$psi_hat_fe, "private_payload") <- secret
      x
    },
    formula_string_attribute = function(x) {
      attr(x$formula_string, "private_payload") <- secret
      x
    },
    raw_atomic_payload = function(x) {
      x$meta$vc_policy$status <- raw_secret
      x
    },
    vc_policy_vector_payload = function(x) {
      x$meta$vc_policy$status <- row_secret
      x
    },
    prior_policy_vector_payload = function(x) {
      x$meta$prior_policy$identity_scope <- row_secret
      x$meta$prior_diagnostic <- x$meta$prior_policy
      x
    },
    package_version_vector_payload = function(x) {
      x$provenance$engine$package_versions <- stats::setNames(
        rep(list(row_secret), 5L),
        c("brms", "posterior", "cmdstanr", "cmdstan", "rstan")
      )
      x
    },
    cmdstan_state_vector_payload = function(x) {
      x$provenance$engine$cmdstan_state$namespace_available <- row_secret
      x
    },
    dropped_columns_vector_payload = function(x) {
      x$param_map$dropped_names <- row_secret
      x$meta$dropped_draw_columns <- row_secret
      x$warnings <- c(
        x$meta$prior_policy$warning,
        pvstackr:::pv_stack_param_drop_warning(
          row_secret,
          x$param_map$map_source
        )
      )
      x
    },
    hash_columns_vector_payload = function(x) {
      x$weight_summary$long_data_hash_columns <- row_secret
      x$provenance$long_data_hash_columns <- row_secret
      x
    }
  )
  for (name in names(mutations)) {
    expect_error(
      pv_fit_stack_psis(
        stack_fit = mutations[[name]](stack_fit),
        psis_weights = stack_psis_weights_fixture(),
        pareto_k = c(PV1 = 0.2, PV2 = 0.4, PV3 = 0.6),
        control = stack_psis_control()
      ),
      info = name
    )
  }

  legacy <- stack_fit
  legacy$schema_version <- "0.1.0"
  legacy$provenance$schema_version <- "0.1.0"
  expect_error(
    pv_fit_stack_psis(
      stack_fit = legacy,
      psis_weights = stack_psis_weights_fixture(),
      pareto_k = c(PV1 = 0.2, PV2 = 0.4, PV3 = 0.6),
      control = stack_psis_control()
    ),
    "inspection-only"
  )
})

test_that("pv_fit_stack_psis Barnard-Rubin metadata remains diagnostic", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.4, PV3 = 0.6)

  fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    df_method = "barnard_rubin",
    df_complete = c(b_x = 35, b_Intercept = 25),
    control = stack_psis_control()
  )

  expect_equal(fit$diagnostics$pooling$df_method, "barnard_rubin")
  expect_equal(fit$diagnostics$pooling$df_complete, c(b_Intercept = 25, b_x = 35), tolerance = 0)
  expect_equal(fit$estimates$df_method, rep("barnard_rubin", 2L))
  expect_equal(fit$estimates$df_complete, c(25, 35), tolerance = 0)
  expect_equal(fit$estimates$interval_role, rep("psis_barnard_rubin", 2L))
  expect_false(any(fit$estimates$coverage_claim_allowed))
  expect_null(get_target(fit))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("self-normalized log ratios are labeled unsmoothed and blocked", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.4, PV3 = 0.6)

  fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    log_ratios = log(weights),
    pareto_k = pareto_k,
    control = stack_psis_control()
  )

  expect_equal(fit$status, "blocked")
  expect_equal(fit$reason_codes, "psis_smoothing_not_applied")
  expect_equal(fit$diagnostics$psis$status, "unsmoothed")
  expect_equal(fit$diagnostics$psis$source, "self_normalized_log_ratios")
  expect_equal(
    fit$diagnostics$psis$weight_method,
    "self_normalized_raw_importance"
  )
  expect_equal(fit$diagnostics$psis$pareto_k_source, "supplied")
  expect_equal(fit$diagnostics$psis$pareto_k, pareto_k)
  expect_equal(fit$diagnostics$psis$weight_ess_iid, c(PV1 = 6, PV2 = 5.4, PV3 = 5.4))
  expect_false(any(c("pooling", "weighted") %in% names(fit$diagnostics)))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_stack_psis can use an injected PSIS function", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3)

  fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    log_ratios = log(weights),
    psis_function = function(log_ratios) {
      list(weights = exp(log_ratios), pareto_k = pareto_k)
    },
    control = stack_psis_control()
  )

  expect_equal(fit$status, "ok")
  expect_equal(fit$diagnostics$psis$source, "injected_psis_function")
  expect_equal(fit$diagnostics$psis$pareto_k_source, "injected_function_output")
  expect_equal(fit$diagnostics$psis$weight_method, "caller_declared_external_psis")
  expect_equal(fit$diagnostics$psis$pareto_k, pareto_k)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("external PSIS provenance is explicit and fail-closed when absent", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3)

  unverified <- pvstackr::pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    control = stack_psis_control(return_draws = TRUE)
  )
  expect_equal(unverified$status, "blocked")
  expect_equal(
    unverified$reason_codes,
    "psis_weight_provenance_incomplete"
  )
  expect_equal(
    unverified$diagnostics$psis$status,
    "provenance_incomplete"
  )
  expect_equal(
    unverified$diagnostics$psis$weight_method,
    "unspecified_external"
  )
  expect_identical(unverified$diagnostics$psis$producer, NA_character_)
  expect_equal(unverified$diagnostics$psis$weight_ess_iid, c(PV1 = 6, PV2 = 5.4, PV3 = 5.4))
  expect_false(any(c("pooling", "weighted") %in% names(unverified$diagnostics)))
  expect_false(unverified$control$return_draws)
  expect_invisible(pvstackr:::validate_pvstackr_fit(unverified))

  fake_function <- pvstackr::pv_fit_stack_psis(
    stacked_draws = draws,
    log_ratios = log(weights),
    psis_function = function(log_ratios) {
      list(weights = exp(log_ratios), pareto_k = pareto_k)
    },
    control = stack_psis_control()
  )
  expect_equal(fake_function$status, "blocked")
  expect_equal(
    fake_function$diagnostics$psis$source,
    "injected_psis_function"
  )
  expect_equal(
    fake_function$diagnostics$psis$weight_method,
    "unspecified_external"
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(fake_function))

  for (blank in c(" ", "  ")) {
    expect_error(
      pvstackr::pv_fit_stack_psis(
        stacked_draws = draws,
        psis_weights = weights,
        pareto_k = pareto_k,
        psis_producer = blank,
        psis_producer_version = "1.0.0",
        control = stack_psis_control()
      ),
      "bounded scalar string"
    )
    expect_error(
      pvstackr::pv_fit_stack_psis(
        stacked_draws = draws,
        psis_weights = weights,
        pareto_k = pareto_k,
        psis_producer = "producer",
        psis_producer_version = blank,
        control = stack_psis_control()
      ),
      "bounded scalar string"
    )
  }

  calls <- 0L
  must_not_run <- function(log_ratios) {
    calls <<- calls + 1L
    stop("injected function must not run for malformed log ratios")
  }
  malformed_log_ratios <- list(
    matrix(Inf, nrow = 6L, ncol = 3L),
    matrix(0, nrow = 2L, ncol = 3L),
    matrix(0, nrow = 6L, ncol = 2L),
    new.env(parent = emptyenv())
  )
  for (candidate in malformed_log_ratios) {
    expect_error(
      pvstackr::pv_fit_stack_psis(
        stacked_draws = draws,
        pv_cols = c("PV1", "PV2", "PV3"),
        log_ratios = candidate,
        psis_function = must_not_run,
        psis_producer = "producer",
        psis_producer_version = "1.0.0",
        control = stack_psis_control()
      ),
      "log_ratios"
    )
  }
  expect_identical(calls, 0L)

  expect_error(
    pvstackr::pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = pareto_k,
      psis_producer = "loo",
      control = stack_psis_control()
    ),
    "supplied together"
  )
  expect_error(
    pvstackr::pv_fit_stack_psis(
      stacked_draws = draws,
      log_ratios = log(weights),
      pareto_k = pareto_k,
      psis_producer = "loo",
      psis_producer_version = "2.9.0",
      control = stack_psis_control()
    ),
    "cannot be relabeled"
  )
  expect_error(
    pvstackr::pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = pareto_k,
      log_ratios = log(weights),
      psis_producer = "loo",
      psis_producer_version = "2.9.0",
      control = stack_psis_control()
    ),
    "exactly one weight route"
  )
})

test_that("weight concentration diagnostics are scale-invariant and tamper-checked", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3)

  base_fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    control = stack_psis_control(return_draws = FALSE)
  )
  scaled_fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = sweep(weights, 2L, c(2, 9, 17), FUN = "*"),
    pareto_k = pareto_k,
    control = stack_psis_control(return_draws = FALSE)
  )
  fields <- c(
    "weight_ess_iid", "weight_ess_fraction", "max_normalized_weight"
  )
  expect_identical(
    base_fit$diagnostics$psis[fields],
    scaled_fit$diagnostics$psis[fields]
  )

  zero_weights <- cbind(
    PV1 = c(1, 0, 0, 0, 1, 0),
    PV2 = c(1, 0, 1, 0, 1, 0),
    PV3 = c(1, 1, 0, 1, 0, 1)
  )
  zero_fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = zero_weights,
    pareto_k = pareto_k,
    control = stack_psis_control(return_draws = TRUE)
  )
  expect_equal(
    zero_fit$diagnostics$psis$weight_ess_iid,
    c(PV1 = 2, PV2 = 3, PV3 = 4),
    tolerance = 1e-14
  )
  expect_equal(
    zero_fit$diagnostics$psis$max_normalized_weight,
    c(PV1 = 0.5, PV2 = 1 / 3, PV3 = 0.25),
    tolerance = 1e-14
  )

  tampered <- base_fit
  tampered$diagnostics$psis$weight_ess_iid[[1L]] <- 5
  tampered$diagnostics$psis$weight_ess_fraction[[1L]] <- 5 / 6
  tampered <- pvstackr:::pv_fit_issue_validation_stamp(tampered)
  expect_error(
    pvstackr:::validate_pvstackr_fit(tampered),
    "canonical bounds or algebra"
  )

  hidden <- base_fit
  attr(hidden$diagnostics$psis$weight_ess_iid, "private") <- "secret"
  hidden <- pvstackr:::pv_fit_issue_validation_stamp(hidden)
  expect_error(
    pvstackr:::validate_pvstackr_fit(hidden),
    "bare named finite vectors"
  )
})

test_that("named weight, Pareto-k, and PV columns align by label", {
  draws <- stack_psis_draws_fixture()
  aligned_weights <- cbind(
    PV1 = rep(1, 6),
    PV2 = c(1, 0, 1, 0, 0, 0),
    PV3 = c(1, 1, 1, 0, 0, 0)
  )
  reverse <- c("PV3", "PV2", "PV1")
  fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    pv_cols = c("PV1", "PV2", "PV3"),
    psis_weights = aligned_weights[, reverse, drop = FALSE],
    pareto_k = c(0.3, 0.2, 0.1),
    control = stack_psis_control(return_draws = TRUE)
  )
  expect_equal(
    fit$diagnostics$psis$weight_ess_iid,
    c(PV1 = 6, PV2 = 2, PV3 = 3),
    tolerance = 1e-14
  )
  expect_equal(
    fit$diagnostics$psis$pareto_k,
    c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3),
    tolerance = 0
  )
  oracle <- pvstackr:::pv_stack_psis_summarize(draws, aligned_weights)
  expect_equal(
    fit$diagnostics$weighted$beta,
    oracle$beta,
    tolerance = 1e-14
  )

  bad_names <- aligned_weights
  colnames(bad_names)[[3L]] <- "PV_OTHER"
  expect_error(
    pv_fit_stack_psis(
      stacked_draws = draws,
      pv_cols = c("PV1", "PV2", "PV3"),
      psis_weights = bad_names,
      pareto_k = c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3),
      control = stack_psis_control()
    ),
    "align exactly"
  )
})

test_that("pv_fit_stack_psis fails closed for every Pareto-k failure", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.3)

  blocked <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    control = stack_psis_control(psis_k_threshold = 0.7)
  )

  expect_equal(blocked$status, "blocked")
  expect_equal(blocked$reason_codes, "psis_k_too_high")
  expect_equal(nrow(blocked$estimates), 0L)
  expect_equal(blocked$diagnostics$psis$status, "failed")
  expect_equal(blocked$diagnostics$psis$bad_pv_cols, "PV2")
  expect_null(blocked$stack_fit)
  expect_null(blocked$design)
  expect_setequal(names(blocked$diagnostics), c("psis", "redaction"))
  expect_equal(blocked$diagnostics$redaction$status, "withheld")
  expect_false(any(c("pooling", "weighted") %in% names(blocked$diagnostics)))
  expect_invisible(pvstackr:::validate_pvstackr_fit(blocked))

  warn_requested <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    fallback = "warn",
    control = stack_psis_control(psis_k_threshold = 0.7)
  )

  expect_equal(warn_requested$status, "blocked")
  expect_equal(warn_requested$reason_codes, "psis_k_too_high")
  expect_equal(nrow(warn_requested$estimates), 0L)
  expect_equal(warn_requested$diagnostics$psis$fallback_requested, "warn")
  expect_equal(warn_requested$diagnostics$psis$fallback_effective, "block")
  expect_invisible(pvstackr:::validate_pvstackr_fit(warn_requested))

  expect_warning(
    pvstackr::pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = pareto_k,
      fallback = "warn",
      control = stack_psis_control(psis_k_threshold = 0.7)
    ),
    "deprecated and behaves as"
  )

  boundary <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = c(PV1 = 0.2, PV2 = 0.7, PV3 = 0.3),
    control = stack_psis_control(psis_k_threshold = 0.7)
  )
  expect_equal(boundary$status, "blocked")
  expect_equal(boundary$reason_codes, "psis_k_too_high")
  expect_equal(boundary$diagnostics$psis$bad_pv_cols, "PV2")

  for (pv in names(pareto_k)) {
    every_pv_boundary <- pareto_k
    every_pv_boundary[] <- 0.2
    every_pv_boundary[[pv]] <- 0.7
    boundary_fit <- pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = every_pv_boundary,
      fallback = "warn",
      control = stack_psis_control(psis_k_threshold = 0.7)
    )
    expect_equal(boundary_fit$status, "blocked", info = pv)
    expect_equal(boundary_fit$diagnostics$psis$bad_pv_cols, pv, info = pv)
    expect_equal(nrow(boundary_fit$estimates), 0L, info = pv)
  }

  not_evaluated <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = c(PV1 = 0.2, PV2 = NA_real_, PV3 = Inf),
    fallback = "warn",
    control = stack_psis_control(psis_k_threshold = 0.7)
  )
  expect_equal(not_evaluated$status, "blocked")
  expect_equal(not_evaluated$reason_codes, "psis_k_not_evaluated")
  expect_equal(not_evaluated$diagnostics$psis$status, "not_evaluated")
  expect_equal(not_evaluated$diagnostics$psis$bad_pv_cols, c("PV2", "PV3"))
  expect_equal(nrow(not_evaluated$estimates), 0L)
  expect_setequal(names(not_evaluated$diagnostics), c("psis", "redaction"))
  expect_invisible(pvstackr:::validate_pvstackr_fit(not_evaluated))
})

test_that("blocked stack_psis objects recursively redact numeric result payloads", {
  blocked <- pv_fit_stack_psis(
    data = stack_psis_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2", "PV3"),
    fit_function = function(...) list(
      draws = stack_psis_draws_fixture(),
      private_data = "BLOCKED_PSIS_BACKEND_SENTINEL"
    ),
    draws_function = function(fit) fit$draws,
    psis_weights = stack_psis_weights_fixture(),
    pareto_k = c(PV1 = 0.2, PV2 = 0.7, PV3 = 0.3),
    fallback = "warn",
    control = stack_psis_control(
      return_draws = TRUE,
      keep_data = TRUE,
      keep_backend_fit = TRUE
    )
  )

  expect_equal(blocked$status, "blocked")
  expect_true(all(vapply(
    blocked$control[c("return_draws", "keep_data", "keep_backend_fit", "keep_log_lik")],
    identical,
    logical(1),
    FALSE
  )))
  expect_null(blocked$design)
  expect_null(blocked$stack_fit)
  expect_null(blocked$draws)
  expect_null(blocked$ccc)
  expect_equal(nrow(blocked$estimates), 0L)

  recursive_names <- function(x) {
    if (!is.list(x)) return(character())
    present <- !vapply(x, is.null, logical(1))
    c(names(x)[present], unlist(lapply(x[present], recursive_names), use.names = FALSE))
  }
  forbidden <- c(
    "pooling", "weighted", "beta", "U", "U_bar", "B", "T_MI",
    "lambda", "df", "df_complete", "weights", "draws",
    "proposal_draws", "stacked_draws",
    "backend_fit", "log_lik"
  )
  payload_names <- recursive_names(blocked)
  # The redaction manifest names withheld fields as character values; it does
  # not recreate them as keys in the object.
  expect_length(intersect(payload_names, forbidden), 0L)
  expect_false(grepl(
    "BLOCKED_PSIS_BACKEND_SENTINEL",
    rawToChar(serialize(blocked, NULL, ascii = TRUE)),
    fixed = TRUE
  ))
  expect_invisible(pvstackr:::validate_pvstackr_fit(blocked))
})

test_that("stack_psis validators preserve stricter thresholds and the full PV universe", {
  reportable <- pv_fit_stack_psis(
    stacked_draws = stack_psis_draws_fixture(),
    psis_weights = stack_psis_weights_fixture(),
    pareto_k = c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4),
    control = stack_psis_control(psis_k_threshold = 0.5)
  )

  loosened <- reportable
  loosened$diagnostics$psis$threshold <- 0.7
  loosened$estimates$psis_k_threshold <- 0.7
  expect_error(
    pvstackr:::validate_pvstackr_fit(loosened),
    "selected fit control threshold"
  )

  missing_pv <- reportable
  missing_pv$diagnostics$psis$pareto_k <- missing_pv$diagnostics$psis$pareto_k[-3L]
  missing_pv$diagnostics$psis$pareto_k_max <- max(missing_pv$diagnostics$psis$pareto_k)
  expect_error(
    pvstackr:::validate_pvstackr_fit(missing_pv),
    "every weighted PV"
  )

  renamed_pv <- reportable
  names(renamed_pv$diagnostics$psis$pareto_k)[[3L]] <- "PV_OTHER"
  expect_error(
    pvstackr:::validate_pvstackr_fit(renamed_pv),
    "every weighted PV"
  )

  legacy_warning <- reportable
  legacy_warning$status <- "warning"
  legacy_warning$reason_codes <- "psis_k_too_high"
  legacy_warning$warnings <- "legacy warning payload"
  expect_error(
    pvstackr:::validate_pvstackr_fit(legacy_warning),
    "legacy unsafe"
  )
})

test_that("blocked stack_psis validator enforces slim and coherent diagnostics", {
  blocked <- pv_fit_stack_psis(
    stacked_draws = stack_psis_draws_fixture(),
    psis_weights = stack_psis_weights_fixture(),
    pareto_k = c(PV1 = 0.2, PV2 = 0.7, PV3 = 0.3),
    fallback = "warn",
    control = stack_psis_control(psis_k_threshold = 0.7)
  )

  with_data <- blocked
  with_data$design <- new_pvstackr_design(
    data = stack_psis_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2", "PV3")
  )
  expect_error(
    pvstackr:::validate_pvstackr_fit(with_data),
    "keep_data = FALSE|redact data"
  )

  aliased <- blocked
  aliased$diagnostics$hidden_results <- list(beta = 1, T_MI = diag(1))
  expect_error(
    pvstackr:::validate_pvstackr_fit(aliased),
    "exact method-specific retention envelope|only slim PSIS"
  )

  wrong_reason <- blocked
  wrong_reason$diagnostics$psis$reason_code <- "psis_k_not_evaluated"
  expect_error(pvstackr:::validate_pvstackr_fit(wrong_reason), "internally inconsistent")

  wrong_bad_pv <- blocked
  wrong_bad_pv$diagnostics$psis$bad_pv_cols <- "PV1"
  expect_error(pvstackr:::validate_pvstackr_fit(wrong_bad_pv), "internally inconsistent")

  wrong_fallback <- blocked
  wrong_fallback$diagnostics$psis$fallback_effective <- "warn"
  expect_error(pvstackr:::validate_pvstackr_fit(wrong_fallback), "fail-closed fallback")

  nested_psis <- blocked
  nested_psis$diagnostics$psis$hidden <- list(beta = 1, T_MI = diag(1))
  expect_error(
    pvstackr:::validate_pvstackr_fit(nested_psis),
    "exact ordered current schema"
  )

  reordered <- blocked
  reordered$diagnostics$psis <-
    reordered$diagnostics$psis[rev(names(reordered$diagnostics$psis))]
  reordered <- pvstackr:::pv_fit_issue_validation_stamp(reordered)
  expect_error(
    pvstackr:::validate_pvstackr_fit(reordered),
    "exact ordered current schema"
  )

  impossible_ess <- blocked
  impossible_ess$diagnostics$psis$weight_ess_iid[[1L]] <- 5
  impossible_ess$diagnostics$psis$weight_ess_fraction[[1L]] <- 5 / 6
  impossible_ess <- pvstackr:::pv_fit_issue_validation_stamp(impossible_ess)
  expect_error(
    pvstackr:::validate_pvstackr_fit(impossible_ess),
    "canonical bounds or algebra"
  )

  nested_redaction <- blocked
  nested_redaction$diagnostics$redaction$payload <- list(weights = 1)
  expect_error(pvstackr:::validate_pvstackr_fit(nested_redaction), "redaction record")

  nested_source <- blocked
  nested_source$diagnostics$psis$source <- list(beta = 1, T_MI = diag(1))
  expect_error(
    pvstackr:::validate_pvstackr_fit(nested_source),
    "weight provenance"
  )

  nested_max <- blocked
  nested_max$diagnostics$psis$pareto_k_max <- list(beta = 0.7)
  expect_error(pvstackr:::validate_pvstackr_fit(nested_max), "field types or source")

  missing_pv <- blocked
  missing_pv$diagnostics$psis$pareto_k <- missing_pv$diagnostics$psis$pareto_k[-1L]
  expect_error(pvstackr:::validate_pvstackr_fit(missing_pv), "declared PV universe")

  hidden_attr <- blocked
  attr(hidden_attr$diagnostics$psis$pareto_k, "private_draws") <-
    matrix(1, 10L, 10L)
  hidden_attr <- pvstackr:::pv_fit_issue_validation_stamp(hidden_attr)
  expect_error(
    pvstackr:::validate_pvstackr_fit(hidden_attr),
    "hidden attributes"
  )

  named_scalar <- blocked
  names(named_scalar$diagnostics$psis$pareto_k_max) <- "private"
  named_scalar <- pvstackr:::pv_fit_issue_validation_stamp(named_scalar)
  expect_error(
    pvstackr:::validate_pvstackr_fit(named_scalar),
    "hidden attributes"
  )

  private_warning <- blocked
  private_warning$warnings <- paste(rep("PRIVATE_DATA", 100L), collapse = "")
  private_warning <- pvstackr:::pv_fit_issue_validation_stamp(private_warning)
  expect_error(
    pvstackr:::validate_pvstackr_fit(private_warning),
    "warnings must exactly reproduce"
  )

  private_provenance <- blocked
  private_provenance$provenance$stacked_source <-
    paste(rep("PRIVATE_DATA", 100L), collapse = "")
  private_provenance <- pvstackr:::pv_fit_issue_validation_stamp(private_provenance)
  expect_error(
    pvstackr:::validate_pvstackr_fit(private_provenance),
    "provenance must exactly reproduce"
  )

  private_redaction <- blocked
  private_redaction$diagnostics$redaction$withheld <-
    paste(rep("PRIVATE_DATA", 100L), collapse = "")
  private_redaction <- pvstackr:::pv_fit_issue_validation_stamp(private_redaction)
  expect_error(
    pvstackr:::validate_pvstackr_fit(private_redaction),
    "redaction record"
  )
})

test_that("legacy PSIS fits migrate only to a slim inspection-only object", {
  current <- pv_fit_stack_psis(
    stacked_draws = stack_psis_draws_fixture(),
    psis_weights = stack_psis_weights_fixture(),
    pareto_k = c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4),
    control = stack_psis_control(return_draws = TRUE)
  )
  current_blocked <- pv_fit_stack_psis(
    stacked_draws = stack_psis_draws_fixture(),
    psis_weights = stack_psis_weights_fixture(),
    pareto_k = c(PV1 = 0.2, PV2 = 0.7, PV3 = 0.4),
    control = stack_psis_control()
  )
  expect_identical(pv_migrate_legacy_psis_fit(current), current)
  expect_identical(pv_migrate_legacy_psis_fit(current_blocked), current_blocked)

  legacy_base <- current
  legacy_base$validation <- NULL
  legacy_base$diagnostics$weighted$private_payload <-
    "PRIVATE_LEGACY_PSIS_SENTINEL"

  warning_high <- legacy_base
  warning_high$status <- "warning"
  warning_high$reason_codes <- "psis_k_too_high"
  warning_high$warnings <- "legacy warning"
  warning_high$diagnostics$psis$pareto_k <-
    c(PV1 = 0.2, PV2 = 0.91, PV3 = 0.4)
  warning_high$diagnostics$psis$pareto_k_max <- 0.91
  warning_high$diagnostics$psis$status <- "warning"
  warning_high$diagnostics$psis$fallback_requested <- "warn"
  warning_high$diagnostics$psis$fallback_effective <- "warn"

  blocked_with_results <- warning_high
  blocked_with_results$status <- "blocked"
  blocked_with_results$diagnostics$psis$status <- "failed"

  not_evaluated <- blocked_with_results
  not_evaluated$reason_codes <- "psis_k_not_evaluated"
  not_evaluated$diagnostics$psis$pareto_k <-
    c(PV1 = 0.2, PV2 = NA_real_, PV3 = Inf)
  not_evaluated$diagnostics$psis$pareto_k_max <- NA_real_
  not_evaluated$diagnostics$psis$status <- "not_evaluated"

  ok_at_boundary <- legacy_base
  ok_at_boundary$diagnostics$psis$pareto_k <-
    c(PV1 = 0.2, PV2 = 0.7, PV3 = 0.4)
  ok_at_boundary$diagnostics$psis$pareto_k_max <- 0.7

  warn_fallback_safe_k <- legacy_base
  warn_fallback_safe_k$diagnostics$psis$fallback_requested <- "warn"
  warn_fallback_safe_k$diagnostics$psis$fallback_effective <- "warn"

  loose_legacy_threshold <- legacy_base
  loose_legacy_threshold$diagnostics$psis$threshold <- 0.9
  loose_legacy_threshold$diagnostics$psis$pareto_k <-
    c(PV1 = 0.2, PV2 = 0.8, PV3 = 0.4)
  loose_legacy_threshold$diagnostics$psis$pareto_k_max <- 0.8

  stricter_threshold_boundary <- legacy_base
  stricter_threshold_boundary$diagnostics$psis$threshold <- 0.5
  stricter_threshold_boundary$diagnostics$psis$pareto_k <-
    c(PV1 = 0.2, PV2 = 0.5, PV3 = 0.4)
  stricter_threshold_boundary$diagnostics$psis$pareto_k_max <- 0.5

  stricter_control_boundary <- legacy_base
  stricter_control_boundary$control$psis_k_threshold <- 0.5
  stricter_control_boundary$diagnostics$psis$pareto_k <-
    c(PV1 = 0.2, PV2 = 0.5, PV3 = 0.4)
  stricter_control_boundary$diagnostics$psis$pareto_k_max <- 0.5

  legacy_fallback_alias <- legacy_base
  legacy_fallback_alias$diagnostics$psis$fallback_requested <- NULL
  legacy_fallback_alias$diagnostics$psis$fallback_effective <- NULL
  legacy_fallback_alias$diagnostics$psis$fallback <- "warn"

  legacy_cases <- list(
    warning_high = warning_high,
    blocked_with_results = blocked_with_results,
    not_evaluated = not_evaluated,
    ok_at_boundary = ok_at_boundary,
    warn_fallback_safe_k = warn_fallback_safe_k,
    loose_legacy_threshold = loose_legacy_threshold,
    stricter_threshold_boundary = stricter_threshold_boundary,
    stricter_control_boundary = stricter_control_boundary,
    legacy_fallback_alias = legacy_fallback_alias
  )
  expected_evidence <- c(
    warning_high = "failed",
    blocked_with_results = "failed",
    not_evaluated = "not_evaluated",
    ok_at_boundary = "failed",
    warn_fallback_safe_k = "legacy_unsafe",
    loose_legacy_threshold = "failed",
    stricter_threshold_boundary = "failed",
    stricter_control_boundary = "failed",
    legacy_fallback_alias = "legacy_unsafe"
  )

  for (name in names(legacy_cases)) {
    source <- legacy_cases[[name]]
    source_before <- serialize(source, NULL, version = 2L)
    inspection <- pv_migrate_legacy_psis_fit(source)

    expect_s3_class(
      inspection,
      "pvstackr_legacy_psis_inspection"
    )
    expect_true(inspection$inspection_only, info = name)
    expect_false(inspection$reportable, info = name)
    expect_identical(
      inspection$diagnostics$psis$evidence_status,
      unname(expected_evidence[[name]]),
      info = name
    )
    expect_identical(serialize(source, NULL, version = 2L), source_before)
    expect_identical(pv_migrate_legacy_psis_fit(inspection), inspection)
    expect_invisible(pvstackr:::pv_validate_legacy_psis_inspection(
      unserialize(serialize(inspection, NULL, version = 2L))
    ))
    if (identical(name, "legacy_fallback_alias")) {
      expect_identical(
        inspection$diagnostics$psis$fallback_requested,
        "warn"
      )
    }
    expect_error(get_estimates(source), info = name)
    expect_error(get_draws(source), info = name)
    expect_error(get_estimates(inspection), "[Ii]nspection-only", info = name)
    expect_error(get_draws(inspection), "[Ii]nspection-only", info = name)
    expect_null(get_target(inspection), info = name)
    expect_identical(
      get_diagnostics(inspection),
      inspection$diagnostics,
      info = name
    )
    expect_output(print(inspection), "withheld", info = name)
    expect_output(print(summary(inspection)), "withheld", info = name)
    expect_false(grepl(
      "PRIVATE_LEGACY_PSIS_SENTINEL",
      rawToChar(serialize(inspection, NULL, ascii = TRUE)),
      fixed = TRUE
    ))
    expect_invisible(
      pvstackr:::pv_validate_legacy_psis_inspection(inspection)
    )
  }

  hidden <- pv_migrate_legacy_psis_fit(warning_high)
  attr(hidden$diagnostics$psis$pareto_k, "private") <- matrix(1:4, 2L)
  expect_error(get_diagnostics(hidden), "Pareto-k evidence")

  injected <- pv_migrate_legacy_psis_fit(warning_high)
  injected$estimates <- warning_high$estimates
  expect_error(get_diagnostics(injected), "fields, order")

  inspection_summary <- summary(pv_migrate_legacy_psis_fit(warning_high))
  inspection_summary$evidence_status <- "legacy_unsafe"
  expect_error(
    print(inspection_summary),
    "validated source inspection"
  )
})

test_that("pv_fit_stack_psis permits stricter but never looser Pareto-k gates", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.5, PV3 = 0.6)

  strict <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    control = stack_psis_control(psis_k_threshold = 0.5)
  )
  expect_equal(strict$status, "blocked")
  expect_equal(strict$diagnostics$psis$bad_pv_cols, c("PV2", "PV3"))

  expect_error(
    stack_psis_control(psis_k_threshold = 0.8),
    "`psis_k_threshold`"
  )
})

test_that("pv_fit dispatches stack_psis through the injected stacked-fit route", {
  data <- stack_psis_fixture_data()
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4)
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  record$n_rows <- integer()

  fake_fit <- function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    record$n <- record$n + 1L
    record$n_rows <- c(record$n_rows, nrow(data))
    list(draws = draws)
  }

  fit <- pv_fit(
    data = data,
    formula = OUTCOME ~ x,
    method = "stack_psis",
    pv_cols = c("PV1", "PV2", "PV3"),
    weight_col = "W",
    control = stack_psis_control(),
    fit_function = fake_fit,
    draws_function = function(fit) fit$draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    psis_producer = "testthat_psis_fixture",
    psis_producer_version = "1.0.0"
  )

  expect_equal(record$n, 1L)
  expect_equal(record$n_rows, nrow(data) * 3L)
  expect_s3_class(fit$stack_fit, "pvstackr_stack_fit")
  expect_s3_class(fit$design, "pvstackr_design")
  expect_equal(fit$design$roles$method, "stack_psis")
  expect_equal(fit$status, "ok")
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit dispatches stack_psis rejects group terms before fitting", {
  data <- transform(
    stack_psis_fixture_data(),
    school = factor(c("A", "A", "B", "B", "C"))
  )
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fake_fit <- function(...) {
    record$n <- record$n + 1L
    list(draws = draws)
  }

  formulas <- list(
    OUTCOME ~ x + (1 | school),
    OUTCOME ~ x + (1 || school)
  )
  for (formula in formulas) {
    expect_error(
      pv_fit(
        data = data,
        formula = formula,
        method = "stack_psis",
        pv_cols = c("PV1", "PV2", "PV3"),
        weight_col = "W",
        control = stack_psis_control(),
        fit_function = fake_fit,
        draws_function = function(fit) fit$draws,
        psis_weights = weights,
        pareto_k = pareto_k
      ),
      "Random-effect/group terms.*`stack_psis`"
    )
  }
  expect_identical(record$n, 0L)
})

test_that("pv_fit_stack_psis guards group terms on formula-bearing precomputed paths", {
  data <- transform(
    stack_psis_fixture_data(),
    school = factor(c("A", "A", "B", "B", "C"))
  )
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4)

  expect_error(
    pv_fit_stack_psis(
      data = data,
      formula = OUTCOME ~ x + (1 | school),
      pv_cols = c("PV1", "PV2", "PV3"),
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = pareto_k,
      control = stack_psis_control()
    ),
    "Random-effect/group terms.*`stack_psis`"
  )
  expect_error(
    pv_fit_stack_psis(
      data = data,
      formula = OUTCOME ~ x + (1 || school),
      pv_cols = c("PV1", "PV2", "PV3"),
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = pareto_k,
      control = stack_psis_control()
    ),
    "Random-effect/group terms.*`stack_psis`"
  )

  fake_fit <- function(...) list(draws = draws)
  grouped_stack_fit <- pvstackr:::pv_stack_fit(
    data = data,
    formula = OUTCOME ~ x + (1 | school),
    pv_cols = c("PV1", "PV2", "PV3"),
    weight_col = "W",
    control = stack_psis_control(),
    fit_function = fake_fit,
    draws_function = function(fit) fit$draws
  )

  expect_error(
    pv_fit_stack_psis(
      stack_fit = grouped_stack_fit,
      psis_weights = weights,
      pareto_k = pareto_k,
      control = stack_psis_control()
    ),
    "Random-effect/group terms.*`stack_psis`"
  )

  fit_without_formula <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    control = stack_psis_control()
  )
  expect_equal(fit_without_formula$status, "ok")
})

test_that("pv_fit_stack_psis allows fixed-effect logical OR formulas", {
  data <- stack_psis_fixture_data()
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4)
  formulas <- list(
    OUTCOME ~ x + I((x > 0) | (id > 2)),
    OUTCOME ~ x + ifelse((x > 0) | (id > 2), 1, 0)
  )

  for (formula in formulas) {
    fit <- pv_fit_stack_psis(
      data = data,
      formula = formula,
      pv_cols = c("PV1", "PV2", "PV3"),
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = pareto_k,
      control = stack_psis_control()
    )

    expect_equal(fit$status, "ok")
    expect_s3_class(fit$design, "pvstackr_design")
    expect_equal(fit$design$roles$method, "stack_psis")
  }
})

test_that("pv_fit_stack_psis rejects malformed PSIS inputs and hollow objects", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()

  expect_error(
    pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights,
      control = stack_psis_control()
    ),
    "pareto_k"
  )
  expect_error(
    pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights[, 1:2],
      pareto_k = c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3),
      control = stack_psis_control()
    ),
    "align"
  )
  expect_error(
    pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights[-1L, ],
      pareto_k = c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3),
      control = stack_psis_control()
    ),
    "complete stacked draw rows"
  )
  expect_error(
    pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3),
      fit_function = function(...) list(),
      control = stack_psis_control()
    ),
    "exactly one"
  )
  expect_error(
    pv_fit_stack_psis(
      stacked_draws = draws,
      psis_weights = weights,
      pareto_k = c(PV1 = 0.1, PV2 = 0.2, PV3 = 0.3),
      control = pv_control(method = "per_pv")
    ),
    "control\\$method"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_psis", control = stack_psis_control()),
    "PSIS"
  )
})

test_that("pv_weighted_mean_cov matches the direct weighted-covariance formula", {
  set.seed(42)
  draws <- matrix(rnorm(400), nrow = 100, ncol = 4,
                  dimnames = list(NULL, paste0("b_", letters[1:4])))
  w <- rexp(100)
  w <- w / sum(w)
  out <- pv_weighted_mean_cov(draws, w)
  mu <- colSums(draws * w)
  cen <- sweep(draws, 2L, mu, FUN = "-")
  denom <- 1 - sum(w^2)
  ref <- matrix(0, 4, 4, dimnames = list(colnames(draws), colnames(draws)))
  for (s in seq_len(nrow(draws))) {
    ref <- ref + w[s] * tcrossprod(cen[s, ])
  }
  ref <- ref / denom
  expect_equal(out$mean, mu)
  expect_equal(out$cov, ref, tolerance = 1e-12)
  # uniform weights must reproduce the ordinary covariance up to the
  # reliability-weights denominator, not an S-dependent inflation
  wu <- rep(1 / 100, 100)
  outu <- pv_weighted_mean_cov(draws, wu)
  expect_equal(outu$cov, stats::cov(draws) * (99 / 100) / (1 - sum(wu^2)),
               tolerance = 1e-12)
})
