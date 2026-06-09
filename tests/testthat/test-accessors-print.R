accessor_reference_draws <- function() {
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

accessor_control <- function(method, ...) {
  pv_control(method = method, backend = "injected", iter = 10L, warmup = 5L, ...)
}

accessor_fit_reference <- function() {
  pv_fit_reference(
    per_pv_draws = accessor_reference_draws(),
    control = accessor_control("per_pv")
  )
}

accessor_stacked_draws <- function() {
  cbind(
    b_Intercept = c(0.9, 1.0, 1.1, 1.2, 1.3, 1.4),
    b_x = c(0.25, 0.35, 0.45, 0.55, 0.65, 0.75),
    sigma = c(0.9, 1.0, 1.1, 1.0, 1.2, 1.1)
  )
}

accessor_weights <- function() {
  cbind(
    PV1 = c(1, 1, 1, 1, 1, 1),
    PV2 = c(1, 2, 1, 2, 1, 2),
    PV3 = c(2, 1, 2, 1, 2, 1)
  )
}

accessor_fit_psis <- function(pareto_k = c(PV1 = 0.2, PV2 = 0.3, PV3 = 0.4), fallback = "block") {
  pv_fit_stack_psis(
    stacked_draws = accessor_stacked_draws(),
    psis_weights = accessor_weights(),
    pareto_k = pareto_k,
    fallback = fallback,
    control = accessor_control("stack_psis")
  )
}

accessor_direct_data <- function() {
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

accessor_fit_direct <- function(return_draws = TRUE) {
  data <- accessor_direct_data()
  target <- pv_brr_target(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )
  fit_function <- function(...) {
    beta <- target$beta[target$fe_names]
    se <- sqrt(diag(target$T_MI))[target$fe_names]
    list(draws = cbind(
      b_Intercept = beta[["b_Intercept"]] + se[["b_Intercept"]] * c(-0.03, -0.01, 0.02, 0.04, -0.02, 0.01, 0.03, -0.04),
      b_x = beta[["b_x"]] + se[["b_x"]] * c(0.02, -0.03, 0.01, -0.04, 0.03, -0.01, 0.04, -0.02),
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
      return_draws = return_draws
    ),
    fit_function = fit_function,
    draws_function = function(fit) fit$draws
  )
}

test_that("pvstackr fit accessors return validated components", {
  fit <- accessor_fit_reference()

  expect_equal(get_estimates(fit), fit$estimates)
  expect_equal(get_target(fit), fit$target)
  expect_null(get_draws(fit))
  expect_equal(get_diagnostics(fit), fit$diagnostics)
  expect_true(is.list(get_diagnostics(fit)$reference$per_pv_draws))

  direct <- accessor_fit_direct()
  expect_equal(get_draws(direct), direct$draws)
  expect_s3_class(get_target(direct), "pvstackr_brr_target")

  no_draws <- accessor_fit_direct(return_draws = FALSE)
  expect_null(get_draws(no_draws))

  blocked_psis <- accessor_fit_psis(c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.4))
  expect_null(get_target(blocked_psis))
  expect_null(get_draws(blocked_psis))
  expect_equal(nrow(get_estimates(blocked_psis)), 0L)
  expect_equal(get_diagnostics(blocked_psis)$psis$status, "failed")
})

test_that("pvstackr method comparison accessors return comparison components", {
  reference <- accessor_fit_reference()
  psis <- accessor_fit_psis()
  comparison <- pv_compare_methods(reference = reference, psis = psis)

  expect_equal(get_estimates(comparison), comparison$estimate_table)
  expect_equal(get_diagnostics(comparison), comparison$diagnostics)
  expect_error(get_target(comparison), "No `get_target\\(\\)` method")
  expect_error(get_draws(comparison), "No `get_draws\\(\\)` method")
})

test_that("accessor defaults reject unsupported objects", {
  expect_error(get_estimates(list()), "No `get_estimates\\(\\)` method")
  expect_error(get_target(list()), "No `get_target\\(\\)` method")
  expect_error(get_draws(list()), "No `get_draws\\(\\)` method")
  expect_error(get_diagnostics(list()), "No `get_diagnostics\\(\\)` method")
})

test_that("pvstackr fit print and summary output are stable", {
  fit <- accessor_fit_reference()

  printed <- paste(capture.output(returned <- print(fit)), collapse = "\n")
  expect_identical(returned, fit)
  expect_match(printed, "pvstackr fit", fixed = TRUE)
  expect_match(printed, "method: per_pv", fixed = TRUE)
  expect_match(printed, "status: ok", fixed = TRUE)
  expect_match(printed, "fixed effects: 2", fixed = TRUE)
  expect_match(printed, "target: per_pv_rubin_draws", fixed = TRUE)
  expect_match(printed, "draws: not retained", fixed = TRUE)
  expect_match(printed, "diagnostics: reference, pooling", fixed = TRUE)
  expect_match(
    printed,
    "interval note: intervals are descriptive rather than coverage-claimable.",
    fixed = TRUE
  )

  summarized <- summary(fit)
  expect_s3_class(summarized, "summary.pvstackr_fit")
  expect_equal(summarized$method, "per_pv")
  expect_equal(summarized$n_terms, 2L)
  expect_equal(summarized$target_source, "per_pv_rubin_draws")
  expect_equal(summarized$diagnostic_keys, c("reference", "pooling"))
  expect_equal(summarized$draw_dim, c(0L, 0L))
  expect_equal(
    summarized$interval_note,
    "intervals are descriptive rather than coverage-claimable."
  )
  expect_equal(summarized$estimates, fit$estimates)

  summary_printed <- paste(capture.output(print(summarized)), collapse = "\n")
  expect_match(summary_printed, "pvstackr fit summary", fixed = TRUE)
  expect_match(summary_printed, "method: per_pv", fixed = TRUE)
  expect_match(summary_printed, "interval note:", fixed = TRUE)
  expect_match(summary_printed, "b_Intercept", fixed = TRUE)
  expect_match(summary_printed, "b_x", fixed = TRUE)
  expect_false(grepl("interval_role", summary_printed, fixed = TRUE))
})

test_that("pvstackr method comparison print and summary output are stable", {
  reference <- accessor_fit_reference()
  warning <- accessor_fit_psis(c(PV1 = 0.2, PV2 = 0.9, PV3 = 0.4), fallback = "warn")
  comparison <- pv_compare_methods(reference = reference, warning = warning)

  printed <- paste(capture.output(returned <- print(comparison)), collapse = "\n")
  expect_identical(returned, comparison)
  expect_match(printed, "pvstackr method comparison", fixed = TRUE)
  expect_match(printed, "reference: reference", fixed = TRUE)
  expect_match(printed, "reference=per_pv", fixed = TRUE)
  expect_match(printed, "warning=stack_psis", fixed = TRUE)
  expect_match(printed, "warnings: warning", fixed = TRUE)
  expect_match(printed, "provenance note:", fixed = TRUE)
  expect_match(
    printed,
    "interval note: intervals are descriptive rather than coverage-claimable.",
    fixed = TRUE
  )

  summarized <- summary(comparison)
  expect_s3_class(summarized, "summary.pvstackr_method_comparison")
  expect_equal(summarized$reference_method, "reference")
  expect_equal(summarized$n_methods, 2L)
  expect_equal(summarized$n_terms, 2L)
  expect_equal(summarized$warning_methods, "warning")
  expect_equal(
    summarized$interval_note,
    "intervals are descriptive rather than coverage-claimable."
  )
  expect_match(summarized$provenance_note, "independent corroboration", fixed = TRUE)
  expect_equal(summarized$estimate_table, comparison$estimate_table)

  summary_printed <- paste(capture.output(print(summarized)), collapse = "\n")
  expect_match(summary_printed, "pvstackr method comparison summary", fixed = TRUE)
  expect_match(summary_printed, "reference: reference", fixed = TRUE)
  expect_match(summary_printed, "warnings: warning", fixed = TRUE)
  expect_match(summary_printed, "provenance note:", fixed = TRUE)
  expect_match(summary_printed, "interval note:", fixed = TRUE)
  expect_match(summary_printed, "max_abs_z_diff", fixed = TRUE)
  expect_false(grepl("interval_role", summary_printed, fixed = TRUE))
})
