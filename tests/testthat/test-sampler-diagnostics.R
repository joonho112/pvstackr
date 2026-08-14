sampler_complete_payload <- function() {
  list(
    rhat_max = 1.008,
    ess_bulk_min = 480,
    ess_tail_min = 360,
    divergences = 0,
    chains = 4L,
    post_warmup_draws_per_chain = 500L
  )
}

test_that("sampler diagnostics normalize the frozen schema", {
  out <- pvstackr:::pv_sampler_diagnostics_normalize(
    sampler_complete_payload(),
    diagnostic_source = "injected_diagnose_function"
  )

  expect_true(all(
    pvstackr:::pv_sampler_diagnostic_required_fields() %in% names(out)
  ))
  expect_equal(out$rhat_max, 1.008)
  expect_equal(out$ess_bulk_min, 480)
  expect_equal(out$ess_tail_min, 360)
  expect_equal(out$ess_bulk_per_chain_min, 120)
  expect_equal(out$ess_tail_per_chain_min, 90)
  expect_identical(out$divergences, 0L)
  expect_identical(out$chains, 4L)
  expect_identical(out$post_warmup_draws_per_chain, 500L)
  expect_identical(out$diagnostic_source, "injected_diagnose_function")
  expect_true(out$diagnostic_complete)
  expect_identical(out$diagnostic_reason_codes, character())
  expect_invisible(pvstackr:::pv_validate_sampler_diagnostics(out))
})

test_that("sampler diagnostics record malformed payload reasons", {
  malformed <- pvstackr:::pv_sampler_diagnostics_normalize(
    list(
      rhat_max = Inf,
      ess_bulk_min = 200,
      ess_tail_min = 150,
      ess_bulk_per_chain_min = 999,
      divergences = -1,
      chains = 2L
    ),
    diagnostic_source = "injected_diagnose_function",
    default_post_warmup_draws_per_chain = 100L
  )

  expect_false(malformed$diagnostic_complete)
  expect_true(any(grepl("rhat_max", malformed$diagnostic_reason_codes)))
  expect_true(any(grepl("divergences", malformed$diagnostic_reason_codes)))
  expect_true("diagnostic_per_chain_mismatch" %in%
                malformed$diagnostic_reason_codes)
  expect_invisible(pvstackr:::pv_validate_sampler_diagnostics(malformed))

  not_list <- pvstackr:::pv_sampler_diagnostics_normalize(
    "bad",
    diagnostic_source = "injected_diagnose_function",
    default_chains = 2L,
    default_post_warmup_draws_per_chain = 100L
  )
  expect_false(not_list$diagnostic_complete)
  expect_identical(
    not_list$diagnostic_reason_codes,
    "diagnostic_payload_not_list"
  )
})

test_that("injected diagnose functions retain only normalized sampler data", {
  with_sampler <- pvstackr:::pv_stack_sampler_diagnostics(
    fit = list(),
    bundled_backend = FALSE,
    diagnose_function = function(fit) c(
      list(n_long = 20L),
      sampler_complete_payload()
    ),
    chains = 4L,
    iter = 600L,
    warmup = 100L
  )
  expect_identical(names(with_sampler), "sampler")
  expect_true(with_sampler$sampler$diagnostic_complete)
  expect_identical(
    with_sampler$sampler$diagnostic_source,
    "injected_diagnose_function"
  )

  failed <- pvstackr:::pv_stack_sampler_diagnostics(
    fit = list(),
    bundled_backend = FALSE,
    diagnose_function = function(fit) stop("backend diagnostic failure"),
    chains = 2L,
    iter = 200L,
    warmup = 50L
  )
  expect_false(failed$sampler$diagnostic_complete)
  expect_identical(
    failed$sampler$diagnostic_reason_codes,
    "diagnostic_extraction_failed"
  )
  expect_identical(failed$sampler$chains, 2L)
  expect_identical(failed$sampler$post_warmup_draws_per_chain, 150L)

  absent <- pvstackr:::pv_stack_sampler_diagnostics(
    fit = list(),
    bundled_backend = FALSE,
    diagnose_function = NULL,
    chains = 2L,
    iter = 200L,
    warmup = 50L
  )
  expect_identical(
    absent$sampler$diagnostic_source,
    "injected_diagnose_function_absent"
  )
  expect_identical(
    absent$sampler$diagnostic_reason_codes,
    "diagnostic_extractor_not_supplied"
  )

  assumed_counts <- sampler_complete_payload()
  assumed_counts$chains <- NULL
  assumed_counts$post_warmup_draws_per_chain <- NULL
  assumed <- pvstackr:::pv_sampler_diagnostics_normalize(
    assumed_counts,
    "injected_diagnose_function",
    default_chains = 4L,
    default_post_warmup_draws_per_chain = 500L
  )
  expect_false(assumed$diagnostic_complete)
  expect_true("diagnostic_counts_from_control_not_observed" %in%
                assumed$diagnostic_reason_codes)
})

test_that("bundled sampler extraction cannot be overridden by a custom hook", {
  authoritative <- pvstackr:::pv_sampler_diagnostics_normalize(
    sampler_complete_payload(),
    "bundled_brms_posterior_and_nuts"
  )
  fake_override <- sampler_complete_payload()
  fake_override$rhat_max <- 9

  out <- pvstackr:::pv_stack_sampler_diagnostics(
    fit = list(),
    bundled_backend = TRUE,
    diagnose_function = function(fit) list(
      n_long = 20L,
      sampler = fake_override
    ),
    chains = 4L,
    iter = 600L,
    warmup = 100L,
    bundled_sampler_function = function(fit) authoritative
  )

  expect_false("n_long" %in% names(out))
  expect_true(out$custom_sampler_override_ignored)
  expect_identical(out$sampler, authoritative)
  expect_equal(out$sampler$rhat_max, 1.008)
  expect_identical(
    out$sampler$diagnostic_source,
    "bundled_brms_posterior_and_nuts"
  )
})

test_that("bundled brms extractor combines posterior summaries and NUTS parameters", {
  fake_draws <- array(
    0,
    dim = c(100L, 2L, 3L),
    dimnames = list(NULL, NULL, c("b_Intercept", "b_x", "sigma"))
  )
  fake_summary <- data.frame(
    variable = c("b_Intercept", "b_x", "sigma"),
    rhat = c(1.003, 1.009, 1.005),
    ess_bulk = c(300, 220, 250),
    ess_tail = c(260, 180, 210)
  )
  fake_nuts <- data.frame(
    Parameter = rep("divergent__", 200),
    Value = c(rep(0, 199), 1)
  )

  out <- pvstackr:::pv_backend_brms_sampler_diagnostics(
    fit = list(),
    draws_array_function = function(fit) fake_draws,
    summary_function = function(draws) fake_summary,
    nuts_function = function(fit) fake_nuts
  )
  expect_true(out$diagnostic_complete)
  expect_equal(out$rhat_max, 1.009)
  expect_equal(out$ess_bulk_min, 220)
  expect_equal(out$ess_tail_min, 180)
  expect_equal(out$ess_bulk_per_chain_min, 110)
  expect_equal(out$ess_tail_per_chain_min, 90)
  expect_identical(out$divergences, 1L)
  expect_identical(out$chains, 2L)
  expect_identical(out$post_warmup_draws_per_chain, 100L)
  expect_identical(
    out$diagnostic_source,
    "bundled_brms_posterior_and_nuts"
  )

  failed <- pvstackr:::pv_backend_brms_sampler_diagnostics(
    fit = list(),
    draws_array_function = function(fit) fake_draws,
    summary_function = function(draws) stop("summary failure"),
    nuts_function = function(fit) fake_nuts
  )
  expect_false(failed$diagnostic_complete)
  expect_identical(
    failed$diagnostic_reason_codes,
    "diagnostic_extraction_failed"
  )

  missing_divergence <- pvstackr:::pv_backend_brms_sampler_diagnostics(
    fit = list(),
    draws_array_function = function(fit) fake_draws,
    summary_function = function(draws) fake_summary,
    nuts_function = function(fit) data.frame(
      Parameter = "treedepth__",
      Value = 5
    )
  )
  expect_false(missing_divergence$diagnostic_complete)
  expect_true(any(grepl(
    "divergences",
    missing_divergence$diagnostic_reason_codes
  )))
})

test_that("sampler validator rejects schema and per-chain tampering", {
  valid <- pvstackr:::pv_sampler_diagnostics_normalize(
    sampler_complete_payload(),
    "injected_diagnose_function"
  )

  bad <- valid
  bad$rhat_max <- NULL
  expect_error(
    pvstackr:::pv_validate_sampler_diagnostics(bad),
    "missing required"
  )

  bad <- valid
  bad$ess_bulk_per_chain_min <- bad$ess_bulk_per_chain_min + 1
  expect_error(
    pvstackr:::pv_validate_sampler_diagnostics(bad),
    "per-chain ESS"
  )

  bad <- valid
  bad$diagnostic_complete <- FALSE
  expect_error(
    pvstackr:::pv_validate_sampler_diagnostics(bad),
    "require reason codes"
  )

  bad <- pvstackr:::pv_sampler_diagnostics_incomplete(
    "injected_diagnose_function",
    "diagnostic_extraction_failed",
    chains = 2L,
    post_warmup_draws_per_chain = 100L
  )
  bad$chains <- -2L
  expect_error(
    pvstackr:::pv_validate_sampler_diagnostics(bad),
    "positive finite|valid integer count"
  )

  bad <- valid
  bad$unexpected <- TRUE
  expect_error(
    pvstackr:::pv_validate_sampler_diagnostics(bad),
    "unsupported field"
  )

  bad <- valid
  bad$diagnostic_source <- "caller_claimed_source"
  expect_error(
    pvstackr:::pv_validate_sampler_diagnostics(bad),
    "not recognized"
  )

  bad <- pvstackr:::pv_sampler_diagnostics_incomplete(
    "injected_diagnose_function",
    "diagnostic_extraction_failed",
    chains = 2L,
    post_warmup_draws_per_chain = 100L
  )
  bad$diagnostic_reason_codes <- "caller_claimed_reason"
  expect_error(
    pvstackr:::pv_validate_sampler_diagnostics(bad),
    "reason code is not recognized"
  )
})

sampler_gate_fixture <- function(chains = 4L, post_warmup = 100L, ...) {
  payload <- test_sampler_payload(
    chains = chains,
    post_warmup = post_warmup,
    ess_bulk_min = 400,
    ess_tail_min = 400,
    ...
  )
  diagnostics <- pvstackr:::pv_sampler_diagnostics_normalize(
    payload,
    "injected_diagnose_function"
  )
  pvstackr:::pv_sampler_gate(
    diagnostics,
    expected_chains = chains,
    expected_post_warmup_draws_per_chain = post_warmup
  )
}

test_that("sampler gate implements exact frozen boundaries", {
  expect_identical(sampler_gate_fixture(rhat_max = 1.01)$status, "ok")

  rhat_105 <- sampler_gate_fixture(rhat_max = 1.05)
  expect_identical(rhat_105$status, "warning")
  expect_identical(rhat_105$reason_codes, "sampler_rhat_warning")

  rhat_red <- sampler_gate_fixture(rhat_max = 1.050001)
  expect_identical(rhat_red$status, "blocked")
  expect_identical(rhat_red$reason_codes, "sampler_rhat_blocked")

  total_100_per_chain_25 <- sampler_gate_fixture(
    ess_bulk_min = 100,
    ess_tail_min = 100
  )
  expect_identical(total_100_per_chain_25$status, "warning")
  expect_equal(total_100_per_chain_25$reason_codes, c(
    "sampler_ess_bulk_warning",
    "sampler_ess_tail_warning"
  ))

  expect_identical(
    sampler_gate_fixture(ess_bulk_min = 400, ess_tail_min = 400)$status,
    "ok"
  )
  expect_identical(
    sampler_gate_fixture(ess_bulk_min = 99)$status,
    "blocked"
  )
  expect_identical(
    sampler_gate_fixture(chains = 5L, ess_tail_min = 100)$status,
    "blocked"
  )
  expect_identical(sampler_gate_fixture(divergences = 0L)$status, "ok")
  expect_identical(sampler_gate_fixture(divergences = 1L)$status, "blocked")
})

test_that("sampler gate keeps deterministic multi-finding order and configuration", {
  multi <- sampler_gate_fixture(
    rhat_max = 1.06,
    ess_bulk_min = 200,
    ess_tail_min = 80,
    divergences = 1L
  )
  expect_identical(multi$status, "blocked")
  expect_equal(multi$reason_codes, c(
    "sampler_rhat_blocked",
    "sampler_ess_bulk_warning",
    "sampler_ess_tail_blocked",
    "sampler_divergences_blocked"
  ))

  diagnostics <- pvstackr:::pv_sampler_diagnostics_normalize(
    test_sampler_payload(chains = 2L, post_warmup = 5L),
    "injected_diagnose_function"
  )
  mismatch <- pvstackr:::pv_sampler_gate(
    diagnostics,
    expected_chains = 4L,
    expected_post_warmup_draws_per_chain = 5L
  )
  expect_identical(mismatch$status, "blocked")
  expect_true("sampler_configuration_blocked" %in% mismatch$reason_codes)
})
