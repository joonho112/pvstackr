skip_if_brms_numeric_fixture_disabled <- function() {
  skip_if_not(
    identical(Sys.getenv("PVSTACKR_RUN_NUMERIC_FIXTURE_TESTS"), "true"),
    paste(
      "brms numeric fixture acceptance is dev-only; set",
      "PVSTACKR_RUN_NUMERIC_FIXTURE_TESTS=true"
    )
  )
}

test_that("deterministic four-chain diagnostics and CCC satisfy numeric acceptance", {
  skip_if_brms_numeric_fixture_disabled()
  brms_numeric_require_real_posterior()
  manifest <- brms_numeric_fixture_manifest()
  bundle <- brms_numeric_target()
  target <- bundle$target

  expect_identical(
    unname(manifest[["Fixture-Schema"]]),
    "pvstackr-brms-numeric-fixture/1"
  )
  expect_identical(
    unname(manifest[["Generator"]]),
    "seeded-whitened-gaussian-v1"
  )
  expect_identical(
    unname(manifest[["RNG-Kind"]]),
    "Mersenne-Twister,Inversion,Rejection"
  )
  expect_identical(target$target_hash, unname(manifest[["Target-Hash"]]))
  expect_identical(
    unname(manifest[["Acceptance-Role"]]),
    "deterministic_numeric_contract"
  )
  expect_identical(
    unname(manifest[["Fixture-Origin"]]),
    "synthetic_not_brmsfit"
  )
  expect_identical(unname(manifest[["Empirical-Backend-Accuracy"]]), "false")
  expect_identical(unname(manifest[["Bundled-Sampling-Tested-Here"]]), "false")
  expect_identical(unname(manifest[["Live-Sampler-Quality-Validated"]]), "false")
  expect_identical(unname(manifest[["Model-Recovery-Validated"]]), "false")
  expect_identical(unname(manifest[["Coverage-Validated"]]), "false")
  expect_identical(unname(manifest[["Real-Data-Evidence"]]), "false")
  expect_identical(
    unname(manifest[["Sampler-Reference"]]),
    "posterior_canonical_reference_parity"
  )

  fixture <- brms_numeric_generate(target, manifest)
  expect_identical(dim(fixture$draws_matrix), c(4000L, 4L))
  expect_identical(
    colnames(fixture$draws_matrix),
    strsplit(unname(manifest[["Parameter-Order"]]), ",", fixed = TRUE)[[1L]]
  )
  expect_identical(dim(fixture$draws_array), c(1000L, 4L, 4L))
  for (variable in seq_len(dim(fixture$draws_array)[[3L]])) {
    expect_identical(
      as.double(fixture$draws_array[, , variable]),
      as.double(fixture$draws_matrix[, variable])
    )
  }

  fixture_hash <- brms_numeric_fixture_hash(fixture, target$target_hash)
  expect_identical(fixture_hash, unname(manifest[["Draws-SHA256"]]))
  mutated_fixture <- fixture
  mutated_fixture$draws_matrix[[1L]] <-
    mutated_fixture$draws_matrix[[1L]] + 1e-8
  expect_false(identical(
    brms_numeric_fixture_hash(mutated_fixture, target$target_hash),
    fixture_hash
  ))
  expect_match(fixture_hash, "^sha256:[0-9a-f]{64}$")
  set.seed(91827L)
  rng_kind_before <- RNGkind()
  rng_seed_before <- .Random.seed
  repeated_fixture <- brms_numeric_generate(target, manifest)
  expect_identical(RNGkind(), rng_kind_before)
  expect_identical(.Random.seed, rng_seed_before)
  expect_identical(repeated_fixture$draws_matrix, fixture$draws_matrix)
  expect_identical(
    brms_numeric_fixture_hash(repeated_fixture, target$target_hash),
    fixture_hash
  )

  fe_draws <- fixture$draws_matrix[, target$fe_names, drop = FALSE]
  raw_covariance <- stats::cov(fe_draws)
  expect_equal(raw_covariance, fixture$raw_covariance, tolerance = 1e-12)
  expect_true(max(abs(raw_covariance[row(raw_covariance) != col(raw_covariance)])) < 1e-12)
  expect_gt(max(abs(target$T_MI[row(target$T_MI) != col(target$T_MI)])), 1)
  raw_gap <- brms_numeric_frobenius_relative(raw_covariance, target$T_MI)
  expect_gt(raw_gap, as.numeric(manifest[["CCC-Raw-Gap-Lower"]]))
  expect_lt(raw_gap, 0.55)

  reference <- brms_numeric_posterior_reference(fixture$draws_array)
  rhat_lower <- as.numeric(manifest[["Rhat-Lower"]])
  rhat_upper <- as.numeric(manifest[["Rhat-Upper"]])
  ess_lower <- as.numeric(manifest[["ESS-Lower"]])
  ess_upper <- as.numeric(manifest[["ESS-Upper"]])
  expect_true(
    reference$rhat_max >= rhat_lower && reference$rhat_max <= rhat_upper
  )
  expect_true(
    reference$ess_bulk_min >= ess_lower && reference$ess_bulk_min <= ess_upper
  )
  expect_true(
    reference$ess_tail_min >= ess_lower && reference$ess_tail_min <= ess_upper
  )
  expect_identical(reference$divergences, 0L)
  expect_identical(reference$chains, 4L)
  expect_identical(reference$post_warmup_draws_per_chain, 1000L)

  nuts_ok <- brms_numeric_nuts(
    fixture$chains,
    fixture$post_warmup,
    divergences = 0L
  )
  extracted <- pvstackr:::pv_backend_brms_sampler_diagnostics(
    fit = list(),
    draws_array_function = function(fit) {
      posterior::as_draws_array(fixture$draws_array)
    },
    nuts_function = function(fit) nuts_ok
  )
  expect_identical(extracted$diagnostic_complete, TRUE)
  expect_identical(extracted$diagnostic_reason_codes, character())
  expect_identical(
    extracted$diagnostic_source,
    "bundled_brms_posterior_and_nuts"
  )
  expect_equal(extracted$rhat_max, reference$rhat_max, tolerance = 0)
  expect_equal(extracted$ess_bulk_min, reference$ess_bulk_min, tolerance = 0)
  expect_equal(extracted$ess_tail_min, reference$ess_tail_min, tolerance = 0)
  expect_equal(
    extracted$ess_bulk_per_chain_min,
    reference$ess_bulk_min / 4,
    tolerance = 0
  )
  expect_equal(
    extracted$ess_tail_per_chain_min,
    reference$ess_tail_min / 4,
    tolerance = 0
  )
  expect_identical(extracted$divergences, 0L)
  expect_identical(extracted$chains, 4L)
  expect_identical(extracted$post_warmup_draws_per_chain, 1000L)
  gate <- pvstackr:::pv_sampler_gate(extracted, 4L, 1000L)
  expect_identical(gate$diagnostic_role, "reportability_gate")
  expect_identical(gate$status, "ok")
  expect_identical(gate$reason_codes, character())
  expect_identical(gate$warnings, character())

  bad_array <- fixture$draws_array
  bad_array[, 1L, "b_x"] <-
    bad_array[, 1L, "b_x"] + 2 * sqrt(target$T_MI["b_x", "b_x"])
  bad_chain <- pvstackr:::pv_backend_brms_sampler_diagnostics(
    fit = list(),
    draws_array_function = function(fit) posterior::as_draws_array(bad_array),
    nuts_function = function(fit) nuts_ok
  )
  bad_chain_gate <- pvstackr:::pv_sampler_gate(bad_chain, 4L, 1000L)
  expect_gt(bad_chain$rhat_max, 1.05)
  expect_lt(bad_chain$ess_bulk_min, 100)
  expect_identical(bad_chain_gate$status, "blocked")
  expect_true("sampler_rhat_blocked" %in% bad_chain_gate$reason_codes)
  expect_true("sampler_ess_bulk_blocked" %in% bad_chain_gate$reason_codes)

  nuts_bad <- brms_numeric_nuts(
    fixture$chains,
    fixture$post_warmup,
    divergences = 1L
  )
  bad_divergence <- pvstackr:::pv_backend_brms_sampler_diagnostics(
    fit = list(),
    draws_array_function = function(fit) {
      posterior::as_draws_array(fixture$draws_array)
    },
    nuts_function = function(fit) nuts_bad
  )
  divergence_gate <- pvstackr:::pv_sampler_gate(bad_divergence, 4L, 1000L)
  expect_identical(bad_divergence$divergences, 1L)
  expect_identical(divergence_gate$status, "blocked")
  expect_true("sampler_divergences_blocked" %in% divergence_gate$reason_codes)

  extraction_failed <- pvstackr:::pv_backend_brms_sampler_diagnostics(
    fit = list(),
    draws_array_function = function(fit) {
      posterior::as_draws_array(fixture$draws_array)
    },
    summary_function = function(draws) stop("forced numeric-fixture failure"),
    nuts_function = function(fit) nuts_ok
  )
  failed_gate <- pvstackr:::pv_sampler_gate(extraction_failed, 4L, 1000L)
  expect_identical(extraction_failed$diagnostic_complete, FALSE)
  expect_identical(
    extraction_failed$diagnostic_reason_codes,
    "diagnostic_extraction_failed"
  )
  expect_identical(failed_gate$status, "blocked")
  expect_true("sampler_diagnostics_incomplete" %in% failed_gate$reason_codes)
  mismatch_gate <- pvstackr:::pv_sampler_gate(extracted, 4L, 999L)
  expect_identical(mismatch_gate$status, "blocked")
  expect_true("sampler_configuration_blocked" %in% mismatch_gate$reason_codes)

  analysis_formula <- stats::as.formula(
    "OUTCOME ~ x + female",
    env = baseenv()
  )
  resolved <- pvstackr:::pv_stack_direct_preflight(
    data = bundle$data,
    formula = analysis_formula,
    target = target,
    family = stats::gaussian(),
    return_model_bundle = TRUE
  )
  ccc <- pvstackr:::ccc_calibrate(
    draws = fixture$draws_matrix,
    target = target,
    param_map = list(fe_names = target$fe_names, vc_names = "sigma"),
    center = "target",
    binding_proof = resolved$preflight$binding_proof
  )
  expect_s3_class(ccc, "pvstackr_ccc")
  calibrated_fe <- ccc$draws_calibrated[, target$fe_names, drop = FALSE]
  mean_error <- max(abs(colMeans(calibrated_fe) - target$beta))
  covariance_error <- brms_numeric_frobenius_relative(
    stats::cov(calibrated_fe),
    target$T_MI
  )
  diagonal_error <- max(
    abs(diag(stats::cov(calibrated_fe) - target$T_MI)) /
      pmax(abs(diag(target$T_MI)), .Machine$double.eps)
  )
  ccc_tolerance <- as.numeric(manifest[["CCC-Relative-Tolerance"]])
  expect_lt(mean_error, 1e-10)
  expect_lt(covariance_error, ccc_tolerance)
  expect_lt(diagonal_error, ccc_tolerance)
  expect_identical(
    ccc$draws_calibrated[, "sigma"],
    fixture$draws_matrix[, "sigma"]
  )

  lower_raw <- t(chol(raw_covariance))
  lower_target <- t(chol(target$T_MI))
  manual_A <- lower_target %*% solve(lower_raw)
  manual_calibrated <- sweep(
    sweep(fe_draws, 2L, colMeans(fe_draws), FUN = "-") %*% t(manual_A),
    2L,
    target$beta,
    FUN = "+"
  )
  colnames(manual_calibrated) <- target$fe_names
  expect_equal(unname(ccc$A), unname(manual_A), tolerance = 1e-12)
  expect_equal(calibrated_fe, manual_calibrated, tolerance = 1e-12)

  expect_identical(ccc$diagnostics$center_status, "ok")
  expect_identical(ccc$diagnostics$conditioning_status, "ok")
  expect_true(
    ccc$diagnostics$kappa_A >= as.numeric(manifest[["CCC-Kappa-Lower"]]) &&
      ccc$diagnostics$kappa_A <= as.numeric(manifest[["CCC-Kappa-Upper"]])
  )
  expect_lt(ccc$diagnostics$rho1, ccc_tolerance)
  expect_lt(ccc$diagnostics$rho2, ccc_tolerance)
  expect_lt(ccc$diagnostics$empirical_fro_rel, ccc_tolerance)

  accepted_shift_relative <- c(
    b_Intercept = 0.005,
    b_x = -0.004,
    b_female = 0.003
  )
  accepted_shift_draws <- fixture$draws_matrix
  accepted_shift_draws[, target$fe_names] <- sweep(
    accepted_shift_draws[, target$fe_names, drop = FALSE],
    2L,
    accepted_shift_relative * target$se[target$fe_names],
    FUN = "+"
  )
  accepted_shift_ccc <- pvstackr:::ccc_calibrate(
    draws = accepted_shift_draws,
    target = target,
    param_map = list(fe_names = target$fe_names, vc_names = "sigma"),
    center = "target",
    binding_proof = resolved$preflight$binding_proof
  )
  expect_identical(accepted_shift_ccc$diagnostics$center_status, "ok")
  expect_equal(
    accepted_shift_ccc$diagnostics$delta_c_by_term,
    abs(accepted_shift_relative),
    tolerance = 1e-10
  )
  expect_equal(accepted_shift_ccc$diagnostics$delta_c_max, 0.005, tolerance = 1e-10)
  expect_gt(
    max(abs(colMeans(accepted_shift_draws[, target$fe_names]) - target$beta)),
    1e-3
  )
  expect_lt(
    max(abs(colMeans(accepted_shift_ccc$draws_fe_cal) - target$beta)),
    1e-10
  )

  wrong_A_draws <- sweep(
    sweep(calibrated_fe, 2L, target$beta, FUN = "-") * 0.9,
    2L,
    target$beta,
    FUN = "+"
  )
  expect_gt(
    brms_numeric_frobenius_relative(
      stats::cov(wrong_A_draws),
      target$T_MI
    ),
    0.1
  )
  shifted_draws <- fixture$draws_matrix
  shifted_draws[, target$fe_names] <- sweep(
    shifted_draws[, target$fe_names, drop = FALSE],
    2L,
    0.1 * target$se[target$fe_names],
    FUN = "+"
  )
  shifted_ccc <- pvstackr:::ccc_calibrate(
    draws = shifted_draws,
    target = target,
    param_map = list(fe_names = target$fe_names, vc_names = "sigma"),
    center = "target",
    binding_proof = resolved$preflight$binding_proof
  )
  expect_identical(shifted_ccc$diagnostics$center_status, "blocked")
  expect_gt(shifted_ccc$diagnostics$delta_c_max, 0.05)

  control <- pv_control(
    method = "stack_direct",
    backend = "injected",
    chains = fixture$chains,
    iter = fixture$post_warmup + 200L,
    warmup = 200L,
    cores = 1L,
    seed = fixture$seed,
    return_draws = TRUE,
    keep_data = FALSE,
    keep_backend_fit = FALSE
  )
  fit <- pv_fit_direct(
    data = bundle$data,
    formula = analysis_formula,
    target = target,
    control = control,
    family = stats::gaussian(),
    fit_function = brms_numeric_fit_function(fixture$draws_matrix),
    draws_function = function(fit, ...) fit$draws,
    diagnose_function = brms_numeric_diagnose_function(reference),
    cache_dir = tempdir(),
    cache_stem = "brms-numeric-fixture"
  )
  expect_s3_class(fit, "pvstackr_fit")
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
  expect_identical(fit$status, "ok")
  expect_identical(fit$reason_codes, character())
  expect_identical(fit$diagnostics$sampler$diagnostic_complete, TRUE)
  expect_identical(
    fit$diagnostics$sampler$diagnostic_source,
    "injected_diagnose_function"
  )
  expect_identical(fit$diagnostics$sampler_gate$status, "ok")
  expect_identical(fit$stack_fit$meta$adapter_source, "injected")
  expect_identical(fit$stack_fit$meta$engine_id, "injected_fit_function")
  expect_identical(dim(fit$draws), c(4000L, 3L))
  expect_lt(max(abs(colMeans(fit$draws) - target$beta)), 1e-10)
  expect_lt(
    brms_numeric_frobenius_relative(stats::cov(fit$draws), target$T_MI),
    ccc_tolerance
  )
  expect_equal(fit$estimates$estimate, unname(target$beta), tolerance = 1e-12)
  expect_equal(fit$estimates$se, unname(target$se), tolerance = 1e-12)
  expect_identical(fit$provenance$target_hash, target$target_hash)
  expect_identical(fit$provenance$ccc_target_hash, target$target_hash)
  expect_identical(
    fit$provenance$target_manifest_hash,
    target$binding_manifest$manifest_hash
  )
  shifted_status <- pvstackr:::pv_fit_direct_status(
    fit$stack_fit,
    shifted_ccc,
    fit$diagnostics$sampler_gate
  )
  expect_identical(shifted_status$status, "blocked")
  expect_true("center_separation_red" %in% shifted_status$reason_codes)

  acceptance_evidence <- list(
    diagnostic_role = "numeric_acceptance",
    numeric_contract_acceptance = TRUE,
    empirical_backend_accuracy = FALSE,
    fixture_origin = "synthetic_not_brmsfit",
    bundled_sampling_tested_here = FALSE,
    live_sampler_quality_validated = FALSE,
    model_recovery_validated = FALSE,
    coverage_validated = FALSE,
    real_data_evidence = FALSE,
    sampler_reference = "posterior_canonical_reference_parity",
    separate_execution_evidence = "step3.2-live-bundled-brms",
    full_path_adapter_source = fit$stack_fit$meta$adapter_source,
    fixture_hash = fixture_hash,
    target_hash = target$target_hash,
    posterior_version = as.character(utils::packageVersion("posterior")),
    sampler_status = gate$status,
    ccc_center_status = ccc$diagnostics$center_status,
    ccc_conditioning_status = ccc$diagnostics$conditioning_status
  )
  expect_identical(acceptance_evidence$diagnostic_role, "numeric_acceptance")
  expect_identical(acceptance_evidence$numeric_contract_acceptance, TRUE)
  expect_identical(acceptance_evidence$empirical_backend_accuracy, FALSE)
  expect_identical(
    acceptance_evidence$fixture_origin,
    "synthetic_not_brmsfit"
  )
  expect_identical(acceptance_evidence$bundled_sampling_tested_here, FALSE)
  expect_identical(acceptance_evidence$live_sampler_quality_validated, FALSE)
  expect_identical(acceptance_evidence$model_recovery_validated, FALSE)
  expect_identical(acceptance_evidence$coverage_validated, FALSE)
  expect_identical(acceptance_evidence$real_data_evidence, FALSE)
  expect_identical(
    acceptance_evidence$sampler_reference,
    "posterior_canonical_reference_parity"
  )
  expect_identical(
    acceptance_evidence$separate_execution_evidence,
    "step3.2-live-bundled-brms"
  )
  expect_identical(acceptance_evidence$full_path_adapter_source, "injected")
  expect_match(acceptance_evidence$fixture_hash, "^sha256:[0-9a-f]{64}$")
  expect_match(acceptance_evidence$target_hash, "^sha256:[0-9a-f]{64}$")
  expect_match(acceptance_evidence$posterior_version, "^[0-9]+[.]")
  expect_identical(acceptance_evidence$sampler_status, "ok")
  expect_identical(acceptance_evidence$ccc_center_status, "ok")
  expect_identical(acceptance_evidence$ccc_conditioning_status, "ok")
})
