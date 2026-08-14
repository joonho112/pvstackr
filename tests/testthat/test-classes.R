classes_fixture_data <- function() {
  data.frame(
    id = seq_len(8),
    school = factor(c("A", "A", "B", "B", "C", "C", "D", "D")),
    x = c(-2, -1, -0.5, 0, 0.5, 1, 1.5, 2),
    PV1 = c(2.1, 2.4, 3.0, 3.2, 3.6, 4.0, 4.2, 4.8),
    PV2 = c(2.0, 2.5, 2.9, 3.4, 3.5, 4.1, 4.4, 4.7),
    W = c(1.0, 1.1, 0.9, 1.2, 1.0, 1.3, 0.8, 1.4),
    RW1 = c(0.8, 1.3, 1.0, 1.1, 0.9, 1.5, 0.9, 1.2),
    RW2 = c(1.2, 0.9, 0.8, 1.4, 1.1, 1.0, 0.7, 1.6),
    RW3 = c(1.1, 1.0, 1.1, 1.0, 1.3, 1.2, 0.8, 1.3),
    RW4 = c(0.9, 1.2, 0.95, 1.3, 1.0, 1.4, 1.1, 1.1)
  )
}

classes_design_fixture <- function() {
  pvstackr:::new_pvstackr_design(
    data = classes_fixture_data(),
    formula = OUTCOME ~ x + (1 | school),
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    roles = list(outcome_placeholder = "OUTCOME", method = "stack_direct"),
    provenance = list(source = "test")
  )
}

classes_direct_design_fixture <- function() {
  pvstackr:::new_pvstackr_design(
    data = classes_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    roles = list(outcome_placeholder = "OUTCOME", method = "stack_direct"),
    provenance = list(source = "test")
  )
}

classes_brr_target_fixture <- function() {
  pv_brr_target(
    data = classes_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )
}

classes_fake_stack_fit <- function(record) {
  function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    record$n <- record$n + 1L
    draws <- cbind(
      b_Intercept = seq(0.1, 0.8, length.out = 8),
      b_x = seq(1.2, 0.5, length.out = 8),
      sigma = seq(0.8, 1.1, length.out = 8),
      lp__ = seq(-4, -3, length.out = 8)
    )
    list(draws = draws, log_lik = matrix(0, nrow = 8, ncol = nrow(data)))
  }
}

classes_stack_fit_fixture <- function(formula = OUTCOME ~ x + (1 | school)) {
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  resolved_model_bundle <- if (pvstackr:::pv_formula_has_random_effect_bar(
    pvstackr:::pv_formula_rhs_checked(formula)
  )) {
    NULL
  } else {
    pvstackr:::pv_binding_resolve_model_bundle(classes_fixture_data(), formula)
  }
  resolved_binding_manifest <- if (is.null(resolved_model_bundle)) {
    NULL
  } else {
    classes_brr_target_fixture()$binding_manifest
  }
  pvstackr:::pv_stack_fit(
    data = classes_fixture_data(),
    formula = formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    control = pv_control(backend = "injected", keep_log_lik = TRUE),
    fit_function = classes_fake_stack_fit(record),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(
      chains = 4L,
      post_warmup = 1000L
    ),
    extract_log_lik = TRUE,
    log_lik_function = function(fit) fit$log_lik,
    resolved_model_bundle = resolved_model_bundle,
    resolved_binding_manifest = resolved_binding_manifest
  )
}

classes_ccc_fixture <- function(
  target = classes_brr_target_fixture(),
  binding_proof = NULL
) {
  beta <- target$beta
  draws <- cbind(
    b_Intercept = beta[["b_Intercept"]] + c(-0.03, -0.01, 0.01, 0.03),
    b_x = beta[["b_x"]] + c(-0.02, 0.02, -0.01, 0.01),
    sigma = c(0.8, 0.9, 1.0, 1.1)
  )
  pvstackr:::ccc_calibrate(
    draws,
    target,
    binding_proof = binding_proof
  )
}

test_that("pvstackr_design constructor records schema and stable hashes", {
  design <- classes_design_fixture()

  expect_s3_class(design, "pvstackr_design")
  expect_equal(design$M, 2L)
  expect_equal(design$R, 4L)
  expect_equal(design$n, nrow(classes_fixture_data()))
  expect_equal(design$pv_cols, c("PV1", "PV2"))
  expect_equal(design$rep_weight_cols, paste0("RW", 1:4))
  expect_equal(design$id_cols, "id")
  expect_equal(design$formula_string, "OUTCOME ~ x + (1 | school)")
  expect_equal(design$data_manifest$n, nrow(classes_fixture_data()))
  expect_equal(design$data_manifest$pv_cols, c("PV1", "PV2"))
  expect_match(design$row_support_hash, "^[0-9a-f]{8}$")
  expect_match(design$pv_value_hash, "^[0-9a-f]{8}$")
  expect_match(design$weight_design_hash, "^[0-9a-f]{8}$")
  expect_match(design$design_hash, "^[0-9a-f]{8}$")
  expect_equal(design$provenance$function_name, "new_pvstackr_design")
  expect_equal(design$provenance$schema_version, design$schema_version)
  expect_invisible(pvstackr:::validate_pvstackr_design(design))
})

test_that("pvstackr_design validator rejects malformed schema and stale hashes", {
  data <- classes_fixture_data()

  expect_error(
    pvstackr:::new_pvstackr_design(data, PV1 ~ x, c("PV1", "PV2")),
    "OUTCOME"
  )
  expect_error(
    pvstackr:::new_pvstackr_design(data, OUTCOME ~ x + weights(W), c("PV1", "PV2")),
    "weights"
  )
  expect_s3_class(
    pvstackr:::new_pvstackr_design(data, OUTCOME ~ x + I((x > 0) | (id > 2)), c("PV1", "PV2")),
    "pvstackr_design"
  )

  bad <- data
  bad$id[2] <- bad$id[1]
  expect_error(
    pvstackr:::new_pvstackr_design(bad, OUTCOME ~ x, c("PV1", "PV2"), id_cols = "id"),
    "unique rows"
  )

  bad <- data
  bad$W[1] <- 0
  expect_error(
    pvstackr:::new_pvstackr_design(bad, OUTCOME ~ x, c("PV1", "PV2"), weight_col = "W"),
    "strictly positive"
  )

  expect_error(
    pvstackr:::new_pvstackr_design(data, OUTCOME ~ x, c("PV1", "PV2"), rep_weight_cols = "RW1"),
    "at least two"
  )

  design <- classes_design_fixture()
  design$data$PV1[1] <- design$data$PV1[1] + 1
  expect_error(pvstackr:::validate_pvstackr_design(design), "hashes")
})

test_that("schema validators accept current target, stack-fit, and CCC objects", {
  target <- classes_brr_target_fixture()
  stack_fit <- classes_stack_fit_fixture()
  ccc <- classes_ccc_fixture()

  expect_invisible(pvstackr:::validate_pvstackr_brr_target(target))
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(stack_fit))
  expect_identical(ccc$schema_version, "0.1.0")
  expect_false("binding_proof" %in% names(ccc))
  expect_invisible(pvstackr:::validate_pvstackr_ccc(ccc))
})

test_that("schema-0.2 CCC validates proof and link before draw payload", {
  target <- classes_brr_target_fixture()
  preflight <- pvstackr:::pv_stack_direct_preflight(
    data = classes_fixture_data(),
    formula = OUTCOME ~ x,
    target = target
  )
  ccc <- classes_ccc_fixture(
    target = target,
    binding_proof = preflight$binding_proof
  )

  bad_proof <- ccc
  bad_proof$draws_calibrated <- "not-a-draw-payload"
  bad_proof$binding_proof$current_manifest_hash <-
    paste0("sha256:", strrep("0", 64L))
  proof_error <- tryCatch(
    pvstackr:::validate_pvstackr_ccc(bad_proof),
    pvstackr_binding_error = identity
  )
  expect_s3_class(proof_error, "pvstackr_binding_error")
  expect_identical(proof_error$code, "PV_BIND_E005")

  bad_link <- ccc
  bad_link$draws_calibrated <- "not-a-draw-payload"
  bad_link$target_hash <- "malformed"
  link_error <- tryCatch(
    pvstackr:::validate_pvstackr_ccc(bad_link),
    pvstackr_binding_error = identity
  )
  expect_s3_class(link_error, "pvstackr_binding_error")
  expect_identical(link_error$code, "PV_BIND_E090")
})

test_that("schema validators reject corrupted target, stack-fit, and CCC objects", {
  target <- classes_brr_target_fixture()
  expect_s3_class(pvstackr:::new_pvstackr_brr_target(target), "pvstackr_brr_target")

  bad_target <- target
  bad_target$target_source <- "raw posterior covariance"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "external_brr_fay_rubin")

  bad_target <- target
  bad_target$M <- 3L
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "per_pv")

  bad_target <- target
  rownames(bad_target$U_bar)[1] <- "wrong"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "U_bar")

  bad_target <- target
  bad_target$target_hash <- "not-a-hash"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "8-character")

  bad_target <- target
  bad_target$policy$target_repair <- "forbidden_by_default"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "target repair")

  bad_target <- target
  bad_target$df_method <- "unsupported"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "df_method")

  bad_target <- target
  bad_target$df_complete <- setNames(rep(20, length(target$fe_names)), target$fe_names)
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "classic")

  bad_target <- target
  bad_target$coverage_claim_allowed <- TRUE
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "coverage_claim_allowed")

  bad_target <- target
  bad_target$interval_role <- "coverage_barnard_rubin"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "interval_role")

  bad_target <- target
  bad_target$policy$df_method <- "barnard_rubin"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "policy df metadata")

  stack_fit <- classes_stack_fit_fixture()
  bad_stack <- stack_fit
  bad_stack$psi_hat_fe <- unname(bad_stack$psi_hat_fe)
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "psi_hat_fe")

  bad_stack <- stack_fit
  bad_stack$psi_hat_fe[1] <- bad_stack$psi_hat_fe[1] + 1
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "draw means")

  bad_stack <- stack_fit
  bad_stack$provenance$long_data_hash <- "00000000"
  expect_error(
    pvstackr:::validate_pvstackr_stack_fit(bad_stack),
    "provenance mirrors|provenance hash"
  )

  bad_stack <- stack_fit
  bad_stack["log_lik"] <- list(NULL)
  expect_error(
    pvstackr:::validate_pvstackr_stack_fit(bad_stack),
    "provenance mirrors|log-likelihood retention"
  )

  bad_stack <- stack_fit
  bad_stack$weight_summary$n_long <- bad_stack$weight_summary$n_long + 1
  expect_error(
    pvstackr:::validate_pvstackr_stack_fit(bad_stack),
    "weight summary|long-data dimensions"
  )

  bad_stack <- stack_fit
  bad_stack$meta$prior_policy$non_flat_prior_warning <- TRUE
  expect_error(
    pvstackr:::validate_pvstackr_stack_fit(bad_stack),
    "canonical fields|legacy prior diagnostic"
  )

  bad_stack <- stack_fit
  bad_stack$meta$prior_diagnostic$warn_nonflat_prior <- TRUE
  expect_error(
    pvstackr:::validate_pvstackr_stack_fit(bad_stack),
    "canonical fields|legacy prior diagnostic"
  )

  bad_stack <- stack_fit
  bad_stack$meta$prior_policy$explicit_prior_warning <- TRUE
  bad_stack$meta$prior_policy$reason_code <- "wrong_prior_code"
  expect_error(
    pvstackr:::validate_pvstackr_stack_fit(bad_stack),
    "prior policy|explicit_prior_warning"
  )

  bad_stack <- stack_fit
  bad_stack$meta$engine_id <- "wrong-engine"
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "backend metadata")

  bad_stack <- stack_fit
  bad_stack$provenance$cache$policy <- "bundled_brms_managed"
  bad_stack$meta$cache <- bad_stack$provenance$cache
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "cache provenance")

  bad_stack <- stack_fit
  bad_stack$provenance["engine"] <- list(NULL)
  bad_stack$provenance["cache"] <- list(NULL)
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "require both")

  bad_stack <- stack_fit
  bad_stack$provenance$engine$requested_backend <- "brms"
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "requested backend")

  bad_stack <- stack_fit
  bad_stack$meta$backend_selection_reason <- "wrong-reason"
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "backend metadata")

  bad_stack <- stack_fit
  bad_stack$provenance$cache$enabled <- FALSE
  bad_stack$provenance$cache$policy <- "disabled"
  bad_stack$meta$cache <- bad_stack$provenance$cache
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "inconsistent tuple")

  legacy_stack <- stack_fit
  legacy_stack$provenance$engine <- NULL
  legacy_stack$provenance$cache <- NULL
  legacy_stack$meta$adapter_source <- NULL
  legacy_stack$meta$resolved_backend <- NULL
  legacy_stack$meta$backend_selection_reason <- NULL
  legacy_stack$meta$cache <- NULL
  legacy_stack$schema_version <- "0.1.0"
  legacy_stack$provenance$schema_version <- "0.1.0"
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(legacy_stack))

  bundled_stack <- stack_fit
  bundled_stack$control$backend <- "brms"
  bundled_stack$provenance$backend <- "brms"
  bundled_stack$provenance$engine <- list(
    adapter_source = "bundled",
    adapter_id = "bundled_brms",
    requested_backend = "brms",
    resolved_backend = "cmdstanr",
    engine_id = "bundled_brms_cmdstanr",
    selection_policy = "cmdstanr_when_namespace_and_cmdstan_configured_else_rstan",
    selection_reason = "configured_cmdstan_selected",
    package_versions = list(),
    cmdstan_state = list(
      namespace_available = TRUE,
      package_version = "0.8.0",
      cmdstan_configured = TRUE,
      cmdstan_version = "2.38.0",
      cmdstan_path_basename = "cmdstan-2.38.0",
      state_reason = "cmdstanr_namespace_and_cmdstan_configured",
      toolchain_checked = FALSE
    )
  )
  bundled_stack$meta$engine_id <- "bundled_brms_cmdstanr"
  bundled_stack$meta$fit_engine <- "bundled_brms"
  bundled_stack$meta$adapter_source <- "bundled"
  bundled_stack$meta$resolved_backend <- "cmdstanr"
  bundled_stack$meta$backend_selection_reason <- "configured_cmdstan_selected"
  bundled_stack$provenance$cache <- list(
    enabled = TRUE,
    policy = "bundled_brms_managed",
    cache_stem = "fixture",
    file_refit = "on_change",
    directory_created = FALSE,
    writable = TRUE
  )
  bundled_stack$meta$cache <- bundled_stack$provenance$cache
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(bundled_stack))

  bad_stack <- bundled_stack
  bad_stack$provenance$engine$cmdstan_state$cmdstan_configured <- FALSE
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "configured CmdStan")

  bad_stack <- bundled_stack
  bad_stack$provenance$engine$selection_reason <- "cmdstanr_namespace_absent_rstan_selected"
  bad_stack$meta$backend_selection_reason <- "cmdstanr_namespace_absent_rstan_selected"
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "configured CmdStan")

  bad_stack <- bundled_stack
  bad_stack$provenance$engine$cmdstan_state$cmdstan_version <- NA_character_
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "configured CmdStan")

  rstan_stack <- bundled_stack
  rstan_stack$provenance$engine$resolved_backend <- "rstan"
  rstan_stack$provenance$engine$engine_id <- "bundled_brms_rstan"
  rstan_stack$provenance$engine$selection_reason <- "cmdstanr_namespace_absent_rstan_selected"
  rstan_stack$provenance$engine$cmdstan_state <- list(
    namespace_available = FALSE,
    package_version = NA_character_,
    cmdstan_configured = FALSE,
    cmdstan_version = NA_character_,
    cmdstan_path_basename = NA_character_,
    state_reason = "cmdstanr_namespace_unavailable",
    toolchain_checked = FALSE
  )
  rstan_stack$meta$engine_id <- "bundled_brms_rstan"
  rstan_stack$meta$resolved_backend <- "rstan"
  rstan_stack$meta$backend_selection_reason <- "cmdstanr_namespace_absent_rstan_selected"
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(rstan_stack))

  bad_stack <- rstan_stack
  bad_stack$provenance$engine$cmdstan_state$state_reason <- "wrong-state"
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "namespace absence")

  bad_stack <- stack_fit
  bad_stack$provenance$cache$policy <- c("disabled", "injected_adapter_managed")
  bad_stack$meta$cache <- bad_stack$provenance$cache
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "policy")

  bad_stack <- stack_fit
  bad_stack$diagnostics$sampler$rhat_max <- NULL
  expect_error(
    pvstackr:::validate_pvstackr_stack_fit(bad_stack),
    "canonical fields|missing required"
  )

  ccc <- classes_ccc_fixture()
  bad_ccc <- ccc
  bad_ccc$target_source <- "mock"
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "external_brr_fay_rubin")

  bad_ccc <- ccc
  bad_ccc$draws_fe_cal[1, 1] <- bad_ccc$draws_fe_cal[1, 1] + 1
  expect_error(
    pvstackr:::validate_pvstackr_ccc(bad_ccc),
    "retained draw matrices|draws_fe_cal"
  )

  bad_ccc <- ccc
  bad_ccc$A_full[3, 3] <- 2
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "A_full")

  bad_ccc <- ccc
  bad_ccc$flags$target_repaired <- TRUE
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "target repair")

  bad_ccc <- ccc
  bad_ccc$policy$target_repair <- "allowed"
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "target repair")

  bad_ccc <- ccc
  bad_ccc$control$allow_target_nearpd <- TRUE
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "not currently supported")

  bad_ccc <- ccc
  bad_ccc$flags$nearpd_target <- TRUE
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "nearPD target")

  bad_ccc <- ccc
  bad_ccc$ccc_status <- "blocked"
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "ccc_status")
})

classes_fit_estimates <- function(ccc, target = classes_brr_target_fixture()) {
  pvstackr:::pv_fit_direct_estimates(ccc, target, conf_level = 0.95)
}

classes_reportable_fit_fixture <- function(status = "ok", reason_codes = character(), warnings = character()) {
  target <- classes_brr_target_fixture()
  formula <- stats::as.formula("OUTCOME ~ x", env = baseenv())
  preflight <- pvstackr:::pv_stack_direct_preflight(
    data = classes_fixture_data(),
    formula = formula,
    target = target
  )
  design <- pvstackr:::new_pvstackr_design(
    data = classes_fixture_data(),
    formula = formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    roles = list(outcome_placeholder = "OUTCOME", method = "stack_direct"),
    provenance = list(
      source = "test",
      target_hash = target$target_hash,
      target_manifest_hash = preflight$binding_proof$target_manifest_hash,
      target_content_hash = target$target_content$target_content_hash,
      binding_verification_policy = preflight$binding_proof$verification_policy
    )
  )
  design$formula <- stats::as.formula(design$formula_string, env = baseenv())
  pvstackr:::validate_pvstackr_design(design)
  stack_fit <- classes_stack_fit_fixture(formula)
  control <- stack_fit$control
  control$method <- "stack_direct"
  control$keep_data <- TRUE
  control$keep_log_lik <- TRUE
  stack_fit <- pvstackr:::pv_stack_fit_composite_projection(
    stack_fit,
    control,
    canonicalize_formula = FALSE
  )
  ccc <- classes_ccc_fixture(
    target = target,
    binding_proof = preflight$binding_proof
  )
  retained_draws <- ccc$draws_fe_cal
  nested_control <- control
  nested_control$return_draws <- FALSE
  stack_fit <- pvstackr:::pv_stack_fit_composite_projection(
    stack_fit,
    nested_control,
    canonicalize_formula = FALSE
  )
  ccc <- pvstackr:::pv_ccc_draw_projection(ccc, FALSE)
  sampler <- stack_fit$diagnostics$sampler
  sampler_gate <- pvstackr:::pv_sampler_gate(
    sampler,
    expected_chains = stack_fit$control$chains,
    expected_post_warmup_draws_per_chain =
      stack_fit$control$iter - stack_fit$control$warmup
  )
  pvstackr:::new_pvstackr_fit(
    method = "stack_direct",
    design = design,
    target = target,
    stack_fit = stack_fit,
    ccc = ccc,
    estimates = classes_fit_estimates(ccc, target),
    draws = retained_draws,
    diagnostics = list(
      preflight = preflight,
      sampler = sampler,
      sampler_gate = sampler_gate,
      stack_fit = stack_fit$diagnostics,
      stack_fit_warnings = stack_fit$warnings,
      ccc = ccc$diagnostics
    ),
    status = status,
    control = control,
    reason_codes = reason_codes,
    provenance = list(
      source = "test",
      target_hash = target$target_hash,
      target_manifest_hash = preflight$binding_proof$target_manifest_hash,
      target_content_hash = target$target_content$target_content_hash,
      binding_verification_policy = preflight$binding_proof$verification_policy
    ),
    warnings = warnings
  )
}

test_that("pvstackr_fit constructor assembles method objects with fixed-effect draws only", {
  fit <- classes_reportable_fit_fixture()

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$method, "stack_direct")
  expect_equal(fit$status, "ok")
  expect_equal(fit$provenance$function_name, "new_pvstackr_fit")
  expect_equal(colnames(fit$draws), c("b_Intercept", "b_x"))
  expect_identical(fit$ccc$schema_version, "0.2.0")
  expect_identical(
    fit$ccc$binding_proof,
    fit$diagnostics$preflight$binding_proof
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  forged <- fit
  forged$diagnostics$sampler$rhat_max <- 1.02
  forged$stack_fit$diagnostics$sampler$rhat_max <- 1.02
  forged$diagnostics$stack_fit <- forged$stack_fit$diagnostics
  expect_error(
    pvstackr:::validate_pvstackr_fit(forged),
    "recomputed"
  )
})

test_that("fit validation tiers bind the exact current owned payload", {
  fit <- classes_reportable_fit_fixture()

  expect_identical(
    names(fit$validation),
    c(
      "schema_version", "policy_id", "canonicalizer_id",
      "fast_path_eligible", "stamp"
    )
  )
  expect_true(fit$validation$fast_path_eligible)
  expect_match(fit$validation$stamp, "^sha256:[0-9a-f]{64}$")
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit, tier = "deep"))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit, tier = "cheap"))

  leaf_attribute <- fit
  attr(leaf_attribute$estimates$estimate, "private_payload") <- "TAMPER"
  expect_error(
    pvstackr:::validate_pvstackr_fit(leaf_attribute, tier = "cheap"),
    "validation stamp"
  )
  expect_error(
    pvstackr:::validate_pvstackr_fit(leaf_attribute, tier = "deep"),
    "validation stamp"
  )

  same_moments <- fit
  same_moments$draws <- same_moments$draws[nrow(same_moments$draws):1L, , drop = FALSE]
  rownames(same_moments$draws) <- NULL
  expect_equal(colMeans(same_moments$draws), colMeans(fit$draws))
  expect_equal(stats::cov(same_moments$draws), stats::cov(fit$draws))
  expect_error(
    pvstackr:::validate_pvstackr_fit(same_moments, tier = "cheap"),
    "validation stamp"
  )

  stamp_tamper <- fit
  stamp_tamper$validation$stamp <- paste0("sha256:", strrep("f", 64L))
  expect_error(
    pvstackr:::validate_pvstackr_fit(stamp_tamper, tier = "cheap"),
    "validation stamp"
  )

  round_trip <- unserialize(serialize(fit, NULL, version = 3L))
  expect_identical(round_trip$validation, fit$validation)
  expect_invisible(pvstackr:::validate_pvstackr_fit(round_trip, tier = "cheap"))

  for (version in c(2L, 3L)) {
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(fit, path, version = version)
    persisted <- readRDS(path)
    expect_identical(persisted$validation, fit$validation, info = paste("RDS", version))
    expect_invisible(
      pvstackr:::validate_pvstackr_fit(persisted, tier = "cheap"),
      label = paste("RDS", version)
    )
  }
})

test_that("reportable stack_direct validator recomputes proof and enforces provenance links", {
  fit <- classes_reportable_fit_fixture()
  mutations <- list(
    preflight_field = function(x) { x$diagnostics$preflight$rhs_string <- "I(x)"; x },
    proof = function(x) {
      x$diagnostics$preflight$binding_proof$current_manifest_hash <-
        paste0("sha256:", strrep("0", 64L))
      x
    },
    ccc_proof = function(x) {
      x$ccc$binding_proof$current_manifest_hash <-
        paste0("sha256:", strrep("0", 64L))
      x
    },
    ccc_proof_missing = function(x) {
      x$ccc$binding_proof <- NULL
      x
    },
    ccc_target_content_link = function(x) {
      x$ccc$target_hash <- paste0("sha256:", strrep("3", 64L))
      x
    },
    design_manifest = function(x) {
      x$design$provenance$target_manifest_hash <- paste0("sha256:", strrep("1", 64L))
      x
    },
    fit_content = function(x) {
      x$provenance$target_content_hash <- paste0("sha256:", strrep("2", 64L))
      x
    },
    verification_policy = function(x) {
      x$provenance$binding_verification_policy <- "unchecked"
      x
    },
    materialized_weight_link = function(x) {
      x$stack_fit$weight_summary$model_matrix_bundle_hash <-
        paste0("sha256:", strrep("4", 64L))
      x
    },
    materialized_meta_link = function(x) {
      x$stack_fit$meta$model_matrix_materialized <- FALSE
      x
    },
    materialized_provenance_link = function(x) {
      x$stack_fit$provenance$model_matrix_bundle_hash <-
        paste0("sha256:", strrep("5", 64L))
      x
    },
    materialized_consistent_bundle_forgery = function(x) {
      forged_hash <- paste0("sha256:", strrep("6", 64L))
      x$stack_fit$weight_summary$model_matrix_bundle_hash <- forged_hash
      x$stack_fit$meta$model_matrix_bundle_hash <- forged_hash
      x$stack_fit$provenance$model_matrix_bundle_hash <- forged_hash
      x
    },
    materialized_formula_string = function(x) {
      x$stack_fit$formula_string <- paste0(x$stack_fit$formula_string, " ")
      x
    },
    materialized_provenance_formula_string = function(x) {
      x$stack_fit$provenance$formula_string <-
        paste0(x$stack_fit$provenance$formula_string, " ")
      x
    },
    materialized_hash_columns = function(x) {
      x$stack_fit$provenance$long_data_hash_columns <-
        rev(x$stack_fit$provenance$long_data_hash_columns)
      x
    },
    materialized_consistent_long_hash_forgery = function(x) {
      x$stack_fit$weight_summary$long_data_hash <- "deadbeef"
      x$stack_fit$meta$long_data_hash <- "deadbeef"
      x$stack_fit$provenance$long_data_hash <- "deadbeef"
      x
    }
  )
  for (mutation in mutations) {
    expect_error(pvstackr:::validate_pvstackr_fit(mutation(fit)))
  }
})

test_that("current stack_direct status tuple rejects reordered or tampered metadata", {
  fit <- classes_reportable_fit_fixture()
  prior_policy <- pvstackr:::pv_stack_prior_policy("fake_prior")
  fit$stack_fit$meta$prior_policy <- prior_policy
  fit$stack_fit$meta$prior_diagnostic <- prior_policy
  fit$stack_fit$warnings <- c(
    prior_policy$warning,
    pvstackr:::pv_stack_param_drop_warning(
      fit$stack_fit$param_map$dropped_names,
      fit$stack_fit$param_map$map_source
    )
  )
  fit$diagnostics$sampler$rhat_max <- 1.02
  fit$stack_fit$diagnostics$sampler$rhat_max <- 1.02
  fit$diagnostics$stack_fit <- fit$stack_fit$diagnostics
  fit$diagnostics$stack_fit_warnings <- fit$stack_fit$warnings
  fit$diagnostics$sampler_gate <- pvstackr:::pv_sampler_gate(
    fit$diagnostics$sampler,
    expected_chains = fit$stack_fit$control$chains,
    expected_post_warmup_draws_per_chain =
      fit$stack_fit$control$iter - fit$stack_fit$control$warmup
  )
  expected <- pvstackr:::pv_fit_direct_status(
    fit$stack_fit,
    fit$ccc,
    fit$diagnostics$sampler_gate
  )
  fit$status <- expected$status
  fit$reason_codes <- expected$reason_codes
  fit$warnings <- expected$warnings
  fit <- pvstackr:::pv_fit_issue_validation_stamp(fit)
  expect_identical(
    fit$reason_codes,
    c("sampler_rhat_warning", "explicit_prior_warning")
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  reordered_reasons <- fit
  reordered_reasons$reason_codes <- rev(reordered_reasons$reason_codes)
  expect_error(
    pvstackr:::validate_pvstackr_fit(reordered_reasons),
    "exactly equal"
  )

  reordered_warnings <- fit
  reordered_warnings$warnings <- rev(reordered_warnings$warnings)
  expect_error(
    pvstackr:::validate_pvstackr_fit(reordered_warnings),
    "exactly equal"
  )

  tampered <- fit
  tampered$reason_codes[[1L]] <- "sampler_rhat_blocked"
  expect_error(
    pvstackr:::validate_pvstackr_fit(tampered),
    "exactly equal"
  )
})

test_that("pvstackr_fit constructor enforces stack_direct reportable status invariants", {
  target <- classes_brr_target_fixture()
  preflight <- pvstackr:::pv_stack_direct_preflight(
    classes_fixture_data(),
    OUTCOME ~ x,
    target
  )
  design <- pvstackr:::new_pvstackr_design(
    data = classes_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    roles = list(outcome_placeholder = "OUTCOME", method = "stack_direct"),
    provenance = list(
      target_hash = target$target_hash,
      target_manifest_hash = preflight$binding_proof$target_manifest_hash,
      target_content_hash = target$target_content$target_content_hash,
      binding_verification_policy = preflight$binding_proof$verification_policy
    )
  )
  fit_provenance <- list(
    target_hash = target$target_hash,
    target_manifest_hash = preflight$binding_proof$target_manifest_hash,
    target_content_hash = target$target_content$target_content_hash,
    binding_verification_policy = preflight$binding_proof$verification_policy
  )
  stack_fit <- classes_stack_fit_fixture(OUTCOME ~ x)
  fit_control <- stack_fit$control
  fit_control$method <- "stack_direct"
  fit_control$keep_data <- TRUE
  fit_control$return_draws <- FALSE
  stack_fit <- pvstackr:::pv_stack_fit_composite_projection(
    stack_fit,
    fit_control,
    canonicalize_formula = FALSE
  )
  ccc <- classes_ccc_fixture(
    target = target,
    binding_proof = preflight$binding_proof
  )
  estimates <- classes_fit_estimates(ccc, target)
  ccc <- pvstackr:::pv_ccc_draw_projection(ccc, FALSE)
  sampler <- stack_fit$diagnostics$sampler
  sampler_gate <- pvstackr:::pv_sampler_gate(
    sampler,
    expected_chains = stack_fit$control$chains,
    expected_post_warmup_draws_per_chain =
      stack_fit$control$iter - stack_fit$control$warmup
  )
  diagnostics <- list(
    preflight = preflight,
    sampler = sampler,
    sampler_gate = sampler_gate,
    stack_fit = stack_fit$diagnostics,
    stack_fit_warnings = stack_fit$warnings,
    ccc = ccc$diagnostics
  )

  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct"),
    "Reportable stack_direct fit requires"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = design,
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = data.frame(),
      diagnostics = diagnostics,
      provenance = fit_provenance,
      control = fit_control
    ),
    "non-empty estimates"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = design,
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = data.frame(term = rev(names(ccc$psi_hat)), estimate = 1:2),
      diagnostics = diagnostics,
      provenance = fit_provenance,
      control = fit_control
    ),
    "estimate terms"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = classes_design_fixture(),
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = estimates,
      diagnostics = diagnostics,
      provenance = fit_provenance,
      control = fit_control
    ),
    class = "pvstackr_binding_error"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = design,
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = estimates,
      status = "blocked",
      reason_codes = "ccc_conditioning_red",
      control = fit_control
    ),
    "Blocked stack_direct fit must not include reportable estimates"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = design,
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = estimates,
      diagnostics = diagnostics,
      provenance = fit_provenance,
      control = within(fit_control, center <- "posterior")
    ),
    "Composite stack_direct|control\\$center = \"target\""
  )

  fit <- classes_reportable_fit_fixture()
  fit["ccc"] <- list(NULL)
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "Reportable stack_direct fit requires"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates <- data.frame(term = "b_x", estimate = 1)
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "estimate terms"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates$target_hash <- "00000000"
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "estimate target metadata"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates$target_source <- "mock"
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "estimate target metadata"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates$parameter_scope <- "variance_component"
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "parameter_scope"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates$estimate[1] <- fit$estimates$estimate[1] + 1
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "external target"
  )

  fit <- classes_reportable_fit_fixture()
  fit$ccc$center <- "posterior"
  fit$ccc$control$center <- "posterior"
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "target-centered CCC"
  )

  fit <- classes_reportable_fit_fixture()
  fit$ccc$diagnostics$center_status <- "blocked"
  fit$ccc$diagnostics$center_reason_code <- "center_separation_red"
  fit$ccc$diagnostics$center_separation$band <- "red"
  fit$ccc$diagnostics$center_separation$reason_code <- "center_separation_red"
  fit$ccc$flags$center_separation_blocked <- TRUE
  fit$ccc$warnings <- "Center separation diagnostic exceeded the block threshold."
  fit$diagnostics$ccc <- fit$ccc$diagnostics
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "blocked center separation"
  )

  fit <- classes_reportable_fit_fixture()
  fit$ccc$diagnostics$center_status <- "warning"
  fit$ccc$diagnostics$center_reason_code <- "center_separation_yellow"
  fit$ccc$diagnostics$center_separation$band <- "yellow"
  fit$ccc$diagnostics$center_separation$reason_code <- "center_separation_yellow"
  fit$ccc$flags$center_separation_warning <- TRUE
  fit$ccc$warnings <- "Center separation diagnostic exceeded the warning threshold."
  fit$diagnostics$ccc <- fit$ccc$diagnostics
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "status = \"warning\""
  )

  fit$status <- "warning"
  fit$reason_codes <- "center_separation_yellow"
  fit$warnings <- "Center separation diagnostic exceeded the warning threshold."
  fit <- pvstackr:::pv_fit_issue_validation_stamp(fit)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  fit <- classes_reportable_fit_fixture()
  fit$ccc$diagnostics$conditioning_status <- "blocked"
  fit$ccc$diagnostics$conditioning_reason_code <- "ccc_conditioning_red"
  fit$ccc$diagnostics$conditioning_band <- "red"
  fit$ccc$diagnostics$conditioning$band <- "red"
  fit$ccc$diagnostics$conditioning$reason_code <- "ccc_conditioning_red"
  fit$ccc$flags$conditioning_blocked <- TRUE
  fit$ccc$flags$kappa_a_blocked <- TRUE
  fit$ccc$warnings <- "CCC conditioning diagnostic exceeded the block threshold."
  fit$diagnostics$ccc <- fit$ccc$diagnostics
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "blocked CCC conditioning"
  )

  fit <- classes_reportable_fit_fixture()
  fit$ccc$diagnostics$conditioning_status <- "warning"
  fit$ccc$diagnostics$conditioning_reason_code <- "ccc_conditioning_yellow"
  fit$ccc$diagnostics$conditioning_band <- "yellow"
  fit$ccc$diagnostics$conditioning$band <- "yellow"
  fit$ccc$diagnostics$conditioning$reason_code <- "ccc_conditioning_yellow"
  fit$ccc$flags$conditioning_warning <- TRUE
  fit$ccc$flags$kappa_a_warning <- TRUE
  fit$ccc$warnings <- "CCC conditioning diagnostic exceeded the warning threshold."
  fit$diagnostics$ccc <- fit$ccc$diagnostics
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "status = \"warning\""
  )

  fit$status <- "warning"
  fit$reason_codes <- "ccc_conditioning_yellow"
  fit$warnings <- "CCC conditioning diagnostic exceeded the warning threshold."
  fit <- pvstackr:::pv_fit_issue_validation_stamp(fit)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  ccc <- classes_reportable_fit_fixture()$ccc
  ccc$flags$conditioning_warning <- TRUE
  expect_error(
    pvstackr:::validate_pvstackr_ccc(ccc),
    "conditioning flags"
  )
})

test_that("pvstackr_fit status metadata distinguishes ok, warning, and blocked objects", {
  expect_error(
    classes_reportable_fit_fixture(status = "ok", reason_codes = "unexpected"),
    "must not include `reason_codes`"
  )
  expect_error(
    classes_reportable_fit_fixture(status = "ok", warnings = "unexpected"),
    "must not include `warnings`"
  )
  expect_error(
    classes_reportable_fit_fixture(status = "warning"),
    "reason_codes"
  )
  expect_error(
    classes_reportable_fit_fixture(status = "warning", reason_codes = "center_separation_yellow"),
    "warnings"
  )

  warning_fit <- classes_reportable_fit_fixture()
  prior_policy <- pvstackr:::pv_stack_prior_policy("fake_prior")
  warning_fit$stack_fit$meta$prior_policy <- prior_policy
  warning_fit$stack_fit$meta$prior_diagnostic <- prior_policy
  warning_fit$stack_fit$warnings <- c(
    prior_policy$warning,
    pvstackr:::pv_stack_param_drop_warning(
      warning_fit$stack_fit$param_map$dropped_names,
      warning_fit$stack_fit$param_map$map_source
    )
  )
  warning_fit$diagnostics$stack_fit_warnings <- warning_fit$stack_fit$warnings
  expected <- pvstackr:::pv_fit_direct_status(
    warning_fit$stack_fit,
    warning_fit$ccc,
    warning_fit$diagnostics$sampler_gate
  )
  warning_fit$status <- expected$status
  warning_fit$reason_codes <- expected$reason_codes
  warning_fit$warnings <- expected$warnings
  warning_fit <- pvstackr:::pv_fit_issue_validation_stamp(warning_fit)
  expect_s3_class(warning_fit, "pvstackr_fit")
  expect_equal(warning_fit$status, "warning")
  expect_equal(warning_fit$reason_codes, "explicit_prior_warning")
  expect_invisible(pvstackr:::validate_pvstackr_fit(warning_fit))

  expect_s3_class(
    pvstackr:::new_pvstackr_fit(
      "stack_direct",
      status = "blocked",
      reason_codes = "preflight_failed",
      control = pv_control(method = "stack_direct", return_draws = FALSE)
    ),
    "pvstackr_fit"
  )
})

test_that("pvstackr_fit constructor enforces method, status, control, and draw contracts", {
  expect_error(
    pvstackr:::new_pvstackr_fit("pipeline_c"),
    "internal manuscript"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("per_pv"),
    "pvstackr_reference_pool"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_psis"),
    "PSIS"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", control = pv_control(method = "per_pv")),
    "control\\$method"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", status = "blocked"),
    "reason_codes"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      "stack_psis",
      status = "blocked",
      reason_codes = "psis_k_too_high",
      control = pv_control(method = "stack_psis", return_draws = FALSE)
    ),
    "PSIS diagnostics"
  )
  blocked <- pvstackr:::new_pvstackr_fit(
    "stack_direct",
    status = "blocked",
    reason_codes = "preflight_failed",
    control = pv_control(method = "stack_direct", return_draws = FALSE)
  )
  expect_s3_class(blocked, "pvstackr_fit")
  expect_equal(blocked$reason_codes, "preflight_failed")

  blocked$estimates <- data.frame(term = "b_x", estimate = 1)
  expect_error(
    pvstackr:::validate_pvstackr_fit(blocked),
    "Blocked stack_direct fit must not include reportable estimates"
  )

  expect_error(
    pvstackr:::new_pvstackr_fit(
      "stack_direct",
      diagnostics = list(
        ccc = list(
          harmless_name = matrix(1, 4L, 4L),
          private_payload = list(data = "PRIVATE_BLOCKED_PAYLOAD")
        )
      ),
      status = "blocked",
      reason_codes = "preflight_failed",
      control = pvstackr:::pv_fit_blocked_control(
        pv_control(method = "stack_direct")
      )
    ),
    "exact current slim variant"
  )

  bad_draws <- cbind(b_Intercept = c(1, 2, 3), sigma = c(0.8, 0.9, 1.0))
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", draws = bad_draws),
    "fixed-effect columns only"
  )

  bad_estimates <- data.frame(term = "sigma", estimate = 1)
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", estimates = bad_estimates),
    "fixed-effect terms"
  )

  target <- classes_brr_target_fixture()
  ccc <- classes_ccc_fixture()
  ccc$target_hash <- "00000000"
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", target = target, ccc = ccc),
    "same `target_hash`"
  )
})
