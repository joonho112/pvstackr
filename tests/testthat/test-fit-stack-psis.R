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
  expect_equal(fit$diagnostics$psis$source, "supplied_weights")
  expect_equal(fit$diagnostics$psis$pareto_k, pareto_k)
  expect_equal(fit$diagnostics$pooling$pooling_source, "stack_psis_rubin_pooling")
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

test_that("pv_fit_stack_psis records self-normalized log-ratio source", {
  draws <- stack_psis_draws_fixture()
  weights <- stack_psis_weights_fixture()
  pareto_k <- c(PV1 = 0.2, PV2 = 0.4, PV3 = 0.6)

  fit <- pv_fit_stack_psis(
    stacked_draws = draws,
    log_ratios = log(weights),
    pareto_k = pareto_k,
    control = stack_psis_control()
  )

  expect_equal(fit$status, "ok")
  expect_equal(fit$diagnostics$psis$source, "log_ratios_self_normalized")
  expect_equal(fit$diagnostics$psis$pareto_k, pareto_k)
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
  expect_equal(fit$diagnostics$psis$source, "psis_function")
  expect_equal(fit$diagnostics$psis$pareto_k, pareto_k)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_stack_psis blocks or warns when Pareto-k fails", {
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
  expect_invisible(pvstackr:::validate_pvstackr_fit(blocked))

  warned <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = pareto_k,
    fallback = "warn",
    control = stack_psis_control(psis_k_threshold = 0.7)
  )

  expect_equal(warned$status, "warning")
  expect_equal(warned$reason_codes, "psis_k_too_high")
  expect_true(nrow(warned$estimates) > 0L)
  expect_equal(warned$estimates$psis_status, rep("failed", nrow(warned$estimates)))
  expect_invisible(pvstackr:::validate_pvstackr_fit(warned))

  boundary <- pv_fit_stack_psis(
    stacked_draws = draws,
    psis_weights = weights,
    pareto_k = c(PV1 = 0.2, PV2 = 0.7, PV3 = 0.3),
    control = stack_psis_control(psis_k_threshold = 0.7)
  )
  expect_equal(boundary$status, "blocked")
  expect_equal(boundary$reason_codes, "psis_k_too_high")
  expect_equal(boundary$diagnostics$psis$bad_pv_cols, "PV2")

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
  expect_invisible(pvstackr:::validate_pvstackr_fit(not_evaluated))
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
    pareto_k = pareto_k
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
