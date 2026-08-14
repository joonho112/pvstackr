article_output_bundle <- function(df_method = "barnard_rubin") {
  df_method <- match.arg(df_method, c("barnard_rubin", "classic"))
  source <- pisa_tiny_parity_load()
  design <- pv_design(
    data = source$data,
    formula = OUTCOME ~ x + female,
    pv_suffix = "READ",
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID"
  )
  target_args <- list(
    data = source$data,
    formula = design$formula,
    pv_cols = design$pv_cols,
    weight_col = design$weight_col,
    rep_weight_cols = design$rep_weight_cols,
    fay_k = design$fay_k,
    id_cols = design$id_cols,
    df_method = df_method
  )
  if (identical(df_method, "barnard_rubin")) {
    target_args$df_complete <- 80
  }
  target <- do.call(pv_brr_target, target_args)
  fit <- pv_fit_direct(
    data = source$data,
    formula = design$formula,
    target = target,
    control = pv_control(
      method = "stack_direct",
      backend = "injected",
      chains = 2L,
      iter = 12L,
      warmup = 6L,
      cores = 1L,
      seed = 20260713L,
      return_draws = TRUE,
      keep_data = FALSE,
      keep_backend_fit = FALSE
    ),
    family = stats::gaussian(),
    fit_function = pisa_tiny_parity_fit_function(target),
    draws_function = function(fit, ...) fit$draws,
    diagnose_function = test_sampler_diagnose_function(
      chains = 2L,
      post_warmup = 6L
    ),
    cache_dir = tempdir(),
    cache_stem = paste0("article-output-", df_method)
  )
  list(data = source$data, design = design, target = target, fit = fit)
}

test_that("article-shaped design target and fit console output is frozen", {
  bundle <- article_output_bundle()

  design_output <- capture.output(design_returned <- print(bundle$design))
  expect_identical(design_returned, bundle$design)
  expect_match(
    design_output[[8L]],
    "^  design hash: [0-9a-f]{8}$"
  )
  expect_identical(
    sub("^  design hash: ", "", design_output[[8L]]),
    bundle$design$design_hash
  )
  normalized_design_output <- design_output
  normalized_design_output[[8L]] <- "  design hash: <HASH8>"
  expect_identical(normalized_design_output, c(
    "pvstackr design",
    "  rows: 12",
    "  formula: OUTCOME ~ x + female",
    "  plausible values: 2",
    "  final weight: W_FSTUWT",
    "  replicate weights: 4",
    "  fay_k: 0.5",
    "  design hash: <HASH8>"
  ))

  target_output <- capture.output(target_returned <- print(bundle$target))
  expect_identical(target_returned, bundle$target)
  expect_identical(target_output, c(
    "pvstackr BRR-Fay target",
    "  fixed effects: 3",
    "  plausible values: 2",
    "  replicate weights: 4",
    "  fay_k: 0.5",
    "  df method: barnard_rubin",
    "  interval role: coverage_barnard_rubin",
    "  source: external_brr_fay_rubin"
  ))

  fit_output <- capture.output(fit_returned <- print(bundle$fit))
  expect_identical(fit_returned, bundle$fit)
  expect_identical(fit_output, c(
    "pvstackr fit",
    "  method: stack_direct",
    "  status: ok",
    "  fixed effects: 3",
    "  target: external_brr_fay_rubin",
    "  draws: 12 x 3",
    paste(
      "  diagnostics: preflight, sampler, sampler_gate, stack_fit,",
      "stack_fit_warnings, ccc"
    )
  ))
})

test_that("Barnard-Rubin estimate rows use t intervals and exact metadata", {
  bundle <- article_output_bundle()
  estimates <- get_estimates(bundle$fit)
  expected_columns <- c(
    "term", "estimate", "se", "std.error", "df", "df_method",
    "df_complete", "conf_level", "conf_low", "conf_high", "conf.low",
    "conf.high", "interval_role", "coverage_claim_allowed",
    "parameter_scope", "target_source", "target_hash"
  )
  expect_identical(names(estimates), expected_columns)
  expect_identical(dim(estimates), c(3L, 17L))
  expect_identical(estimates$term, bundle$target$fe_names)
  expect_equal(
    estimates$estimate,
    unname(bundle$target$beta[bundle$target$fe_names]),
    tolerance = 0
  )
  expect_equal(
    estimates$se,
    unname(sqrt(diag(bundle$target$T_MI))[bundle$target$fe_names]),
    tolerance = 0
  )
  expect_identical(estimates$std.error, estimates$se)
  expect_equal(
    estimates$df,
    unname(bundle$target$df[bundle$target$fe_names]),
    tolerance = 0
  )
  expect_identical(estimates$df_method, rep("barnard_rubin", 3L))
  expect_identical(estimates$df_complete, rep(80, 3L))
  expect_identical(estimates$conf_level, rep(0.95, 3L))

  target_beta <- unname(bundle$target$beta[bundle$target$fe_names])
  target_se <- unname(sqrt(diag(bundle$target$T_MI))[bundle$target$fe_names])
  target_df <- unname(bundle$target$df[bundle$target$fe_names])
  t_critical <- stats::qt(0.975, df = target_df)
  expected_low <- target_beta - t_critical * target_se
  expected_high <- target_beta + t_critical * target_se
  expect_equal(estimates$conf_low, expected_low, tolerance = 1e-12)
  expect_equal(estimates$conf_high, expected_high, tolerance = 1e-12)
  expect_identical(estimates$conf.low, estimates$conf_low)
  expect_identical(estimates$conf.high, estimates$conf_high)

  normal_critical <- stats::qnorm(0.975)
  expect_true(all(t_critical > normal_critical))
  expect_gt(min(t_critical - normal_critical), 1e-6)
  expect_true(all(estimates$conf_low < estimates$estimate))
  expect_true(all(estimates$conf_high > estimates$estimate))

  expect_identical(
    estimates$interval_role,
    rep("coverage_barnard_rubin", 3L)
  )
  expect_identical(estimates$coverage_claim_allowed, rep(TRUE, 3L))
  expect_identical(estimates$parameter_scope, rep("fixed_effect", 3L))
  expect_identical(
    estimates$target_source,
    rep("external_brr_fay_rubin", 3L)
  )
  expect_identical(
    estimates$target_hash,
    rep(bundle$target$target_hash, 3L)
  )
  expect_match(unique(estimates$target_hash), "^sha256:[0-9a-f]{64}$")
})

test_that("article fit accessors expose exact reportable shapes", {
  bundle <- article_output_bundle()
  fit <- bundle$fit

  estimates <- get_estimates(fit)
  target <- get_target(fit)
  draws <- get_draws(fit)
  diagnostics <- get_diagnostics(fit)
  expect_identical(dim(estimates), c(3L, 17L))
  expect_s3_class(target, "pvstackr_brr_target")
  expect_identical(target$target_hash, fit$provenance$target_hash)
  expect_identical(dim(draws), c(12L, 3L))
  expect_identical(colnames(draws), bundle$target$fe_names)
  expect_false("sigma" %in% colnames(draws))

  expect_identical(names(diagnostics), c(
    "preflight", "sampler", "sampler_gate", "stack_fit",
    "stack_fit_warnings", "ccc"
  ))
  expect_identical(diagnostics$sampler$diagnostic_complete, TRUE)
  expect_identical(diagnostics$sampler_gate$diagnostic_role, "reportability_gate")
  expect_identical(diagnostics$sampler_gate$status, "ok")
  expect_identical(diagnostics$ccc$center_status, "ok")
  expect_identical(diagnostics$ccc$conditioning_status, "ok")
  expect_null(fit$stack_fit$stacked_draws)
  expect_null(fit$ccc$draws_calibrated)
  expect_null(fit$ccc$draws_fe_cal)

})

test_that("coverage metadata fails closed to descriptive under classic Rubin df", {
  bundle <- article_output_bundle("classic")
  estimates <- get_estimates(bundle$fit)

  expect_identical(estimates$df_method, rep("classic", 3L))
  expect_true(all(is.na(estimates$df_complete)))
  expect_identical(
    estimates$interval_role,
    rep("descriptive_classic_rubin", 3L)
  )
  expect_identical(estimates$coverage_claim_allowed, rep(FALSE, 3L))
  expect_identical(bundle$fit$status, "ok")
  printed <- capture.output(print(bundle$fit))
  expect_true(any(grepl(
    "interval note: intervals are descriptive rather than coverage-claimable.",
    printed,
    fixed = TRUE
  )))
  summarized <- summary(bundle$fit)
  expect_identical(
    summarized$interval_note,
    "intervals are descriptive rather than coverage-claimable."
  )
})
