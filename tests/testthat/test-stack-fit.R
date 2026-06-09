stack_fixture_data <- function() {
  data.frame(
    id = 1:5,
    school = factor(c("A", "A", "B", "B", "C")),
    x = c(-1, -0.5, 0, 0.5, 1),
    z = c(0, 1, 0, 1, 0),
    W = c(10, 20, 30, 40, 50),
    PV1 = c(1.1, 2.1, 2.9, 4.2, 4.8),
    PV2 = c(0.9, 2.3, 3.2, 4.0, 5.1),
    PV3 = c(1.0, 2.0, 3.1, 4.1, 5.0)
  )
}

fake_stack_fit <- function(record) {
  function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    record$n <- record$n + 1L
    record$formula <- formula
    record$data <- data
    record$args <- list(
      family = family,
      prior = prior,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      backend = backend,
      file = file,
      file_refit = file_refit,
      extra = list(...)
    )
    y <- data$.pvstackr_y
    x <- data$x
    beta <- stats::coef(stats::lm(y ~ x))
    draws <- cbind(
      b_Intercept = beta[[1]] + seq(-0.02, 0.02, length.out = 8),
      b_x = beta[[2]] + seq(0.03, -0.03, length.out = 8),
      sigma = rep(stats::sd(stats::residuals(stats::lm(y ~ x))), 8),
      lp__ = seq(-10, -9, length.out = 8)
    )
    list(draws = draws, log_lik = matrix(0, nrow = 8, ncol = nrow(data)), data = data)
  }
}

fake_custom_stack_fit <- function(record) {
  function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    record$n <- record$n + 1L
    y <- data$.pvstackr_y
    x <- data$x
    beta <- stats::coef(stats::lm(y ~ x))
    draws <- cbind(
      theta_Intercept = beta[[1]] + seq(-0.02, 0.02, length.out = 8),
      theta_x = beta[[2]] + seq(0.03, -0.03, length.out = 8),
      tau = rep(stats::sd(stats::residuals(stats::lm(y ~ x))), 8),
      lp__ = seq(-10, -9, length.out = 8)
    )
    list(draws = draws, log_lik = matrix(0, nrow = 8, ncol = nrow(data)), data = data)
  }
}

test_that("stack preparation uses constant fractional weights without base weights", {
  data <- stack_fixture_data()
  out <- pvstackr:::pv_prepare_stack_data(data, OUTCOME ~ x + z, c("PV1", "PV2", "PV3"))

  expect_equal(nrow(out$data), nrow(data) * 3L)
  expect_equal(out$data$.pvstackr_weight, rep(1 / 3, nrow(out$data)))
  expect_true(all(pvstackr:::pv_stack_reserved_cols() %in% names(out$data)))
  expect_equal(unname(as.integer(table(out$data$.pvstackr_row))), rep(3L, nrow(data)))
  expect_equal(levels(out$data$.pvstackr_pv), c("PV1", "PV2", "PV3"))
  expect_equal(out$weight_summary$total_long_weight, nrow(data))
  expect_equal(unname(out$weight_summary$per_pv_weight_sum), rep(nrow(data) / 3, 3))
  expect_match(out$weight_summary$long_data_hash, "^[0-9a-f]{8}$")
  expect_match(out$formula_string, "weights\\(.pvstackr_weight\\)")
  expect_false(grepl("OUTCOME", out$formula_string, fixed = TRUE))
})

test_that("stack preparation normalizes base weights before fractional scaling", {
  data1 <- stack_fixture_data()
  data2 <- data1
  data2$W <- data2$W * 1000
  out1 <- pvstackr:::pv_prepare_stack_data(data1, OUTCOME ~ x, c("PV1", "PV2", "PV3"), weight_col = "W")
  out2 <- pvstackr:::pv_prepare_stack_data(data2, OUTCOME ~ x, c("PV1", "PV2", "PV3"), weight_col = "W")

  expect_equal(out1$data$.pvstackr_weight, out2$data$.pvstackr_weight)
  expect_equal(mean(out1$data$.pvstackr_weight), 1 / 3)
  expect_equal(sum(out1$data$.pvstackr_weight), nrow(data1))
  expect_equal(unname(out1$weight_summary$per_pv_weight_sum), rep(nrow(data1) / 3, 3))
})

test_that("stack preparation preserves PV order and formula structure", {
  data <- stack_fixture_data()
  out <- pvstackr:::pv_prepare_stack_data(data, OUTCOME ~ x + (1 | school), c("PV3", "PV2", "PV1"), weight_col = "W")

  expect_equal(levels(out$data$.pvstackr_pv), c("PV3", "PV2", "PV1"))
  expect_match(out$formula_string, "\\(1 \\| school\\)")

  logical_or <- pvstackr:::pv_prepare_stack_data(
    data,
    OUTCOME ~ x + I((z > 0) | (x > 0)),
    c("PV1", "PV2", "PV3"),
    weight_col = "W"
  )
  expect_match(logical_or$formula_string, "I\\(", fixed = FALSE)
})

test_that("stack preparation rejects invalid inputs", {
  data <- stack_fixture_data()

  expect_error(pvstackr:::pv_prepare_stack_data(data, y ~ x, c("PV1", "PV2")), "OUTCOME")
  expect_error(pvstackr:::pv_prepare_stack_data(data, OUTCOME | weights(W) ~ x, c("PV1", "PV2")), "OUTCOME")
  expect_error(pvstackr:::pv_prepare_stack_data(data, OUTCOME ~ x, c("PV1")), "at least two")
  expect_error(pvstackr:::pv_prepare_stack_data(data, OUTCOME ~ x, c("PV1", "PV1")), "unique")
  expect_error(pvstackr:::pv_prepare_stack_data(data, OUTCOME ~ x, c("PV1", "NOPE")), "not found")

  bad <- data
  bad$.pvstackr_y <- 1
  expect_error(pvstackr:::pv_prepare_stack_data(bad, OUTCOME ~ x, c("PV1", "PV2")), "reserved")

  bad <- data
  bad$W[1] <- NA_real_
  expect_error(pvstackr:::pv_prepare_stack_data(bad, OUTCOME ~ x, c("PV1", "PV2"), weight_col = "W"), "finite")

  bad <- data
  bad$PV1[1] <- Inf
  expect_error(pvstackr:::pv_prepare_stack_data(bad, OUTCOME ~ x, c("PV1", "PV2")), "finite and numeric")

  bad <- data
  bad$W[1] <- -1
  expect_error(pvstackr:::pv_prepare_stack_data(bad, OUTCOME ~ x, c("PV1", "PV2"), weight_col = "W"), "strictly positive")

  bad <- data
  bad$W[1] <- 0
  expect_error(pvstackr:::pv_prepare_stack_data(bad, OUTCOME ~ x, c("PV1", "PV2"), weight_col = "W"), "strictly positive")
})

test_that("stack parameter map reports dropped automatic draw columns", {
  draws <- matrix(seq_len(28), nrow = 4)
  colnames(draws) <- c(
    "b_Intercept",
    "b_x",
    "sigma",
    "sd_school__Intercept",
    "cor_school__Intercept__x",
    "r_school[A,Intercept]",
    "lp__"
  )

  map <- pvstackr:::pv_stack_param_map(draws)

  expect_equal(colnames(map$draws_selected), c("b_Intercept", "b_x", "sigma", "sd_school__Intercept"))
  expect_equal(map$fe_names, c("b_Intercept", "b_x"))
  expect_equal(map$vc_names, c("sigma", "sd_school__Intercept"))
  expect_equal(map$dropped_names, c("cor_school__Intercept__x", "r_school[A,Intercept]", "lp__"))
  expect_identical(map$map_source, "auto_regex")
  expect_match(map$warnings, "automatic regex")
  expect_match(map$warnings, "lp__")
})

test_that("stack parameter map supports explicit custom draw names", {
  draws <- matrix(seq_len(20), nrow = 5)
  colnames(draws) <- c("theta_Intercept", "theta_x", "tau", "lp__")

  map <- pvstackr:::pv_stack_param_map(
    draws,
    param_map = list(
      fe_names = c("theta_Intercept", "theta_x"),
      vc_names = "tau"
    )
  )

  expect_equal(colnames(map$draws_selected), c("theta_Intercept", "theta_x", "tau"))
  expect_equal(map$fe_names, c("theta_Intercept", "theta_x"))
  expect_equal(map$vc_names, "tau")
  expect_equal(map$dropped_names, "lp__")
  expect_identical(map$map_source, "explicit")
  expect_match(map$warnings, "explicit")

  map_no_vc <- pvstackr:::pv_stack_param_map(
    draws,
    param_map = list(fe_names = c("theta_Intercept", "theta_x"), vc_names = character())
  )
  expect_equal(colnames(map_no_vc$draws_selected), c("theta_Intercept", "theta_x"))
  expect_equal(map_no_vc$vc_names, character())
  expect_equal(map_no_vc$dropped_names, c("tau", "lp__"))

  distributional <- matrix(seq_len(25), nrow = 5)
  colnames(distributional) <- c("theta_Intercept", "theta_x", "b_sigma_x", "tau", "lp__")
  map_distributional <- pvstackr:::pv_stack_param_map(
    distributional,
    param_map = list(fe_names = c("theta_Intercept", "theta_x"), vc_names = character())
  )
  expect_equal(map_distributional$fe_names, c("theta_Intercept", "theta_x"))
  expect_equal(map_distributional$vc_names, character())
  expect_equal(map_distributional$dropped_names, c("b_sigma_x", "tau", "lp__"))
  expect_identical(map_distributional$map_source, "explicit")
})

test_that("stack parameter map rejects malformed explicit maps", {
  draws <- matrix(seq_len(20), nrow = 5)
  colnames(draws) <- c("theta_Intercept", "theta_x", "tau", "lp__")

  expect_error(pvstackr:::pv_stack_param_map(draws, param_map = "bad"), "param_map")
  expect_error(pvstackr:::pv_stack_param_map(draws, param_map = list(vc_names = "tau")), "fe_idx")
  expect_error(
    pvstackr:::pv_stack_param_map(draws, param_map = list(fe_names = "missing")),
    "not found"
  )
  expect_error(
    pvstackr:::pv_stack_param_map(draws, param_map = list(fe_idx = c(1L, 1L))),
    "unique draw columns"
  )
  expect_error(
    pvstackr:::pv_stack_param_map(draws, param_map = list(fe_names = "theta_Intercept", vc_names = "theta_Intercept")),
    "must not overlap"
  )
  expect_error(
    pvstackr:::pv_stack_param_map(draws, param_map = list(fe_idx = 1L, fe_names = "theta_Intercept")),
    "must not supply both"
  )
})

test_that("stack fit performs exactly one injected fit and extracts selected draws", {
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  control <- pv_control(
    iter = 10L,
    warmup = 5L,
    chains = 2L,
    seed = 20260514L,
    backend = "injected",
    keep_backend_fit = TRUE,
    keep_log_lik = TRUE,
    keep_data = TRUE
  )

  out <- pvstackr:::pv_stack_fit(
    data = stack_fixture_data(),
    formula = OUTCOME ~ x + z + (1 | school),
    pv_cols = c("PV1", "PV2", "PV3"),
    weight_col = "W",
    control = control,
    family = "fake_gaussian",
    prior = "fake_prior",
    fit_function = fake_stack_fit(record),
    draws_function = function(fit) fit$draws,
    diagnose_function = function(fit) list(n_long = nrow(fit$data)),
    extract_log_lik = TRUE,
    log_lik_function = function(fit) fit$log_lik,
    cache_dir = tempdir(),
    cache_stem = "single-long-fit",
    additional_args = list(sample_prior = "no")
  )

  expect_s3_class(out, "pvstackr_stack_fit")
  expect_equal(record$n, 1L)
  expect_match(paste(deparse(record$formula), collapse = ""), "weights\\(.pvstackr_weight\\)")
  expect_match(paste(deparse(record$formula), collapse = ""), "\\(1 \\| school\\)")
  expect_equal(nrow(record$data), nrow(stack_fixture_data()) * 3L)
  expect_equal(record$args$chains, 2L)
  expect_equal(record$args$iter, 10L)
  expect_equal(record$args$warmup, 5L)
  expect_equal(record$args$backend, "injected")
  expect_match(record$args$file, "single-long-fit$")
  expect_equal(record$args$extra$sample_prior, "no")
  expect_equal(colnames(out$stacked_draws), c("b_Intercept", "b_x", "sigma"))
  expect_false("lp__" %in% colnames(out$stacked_draws))
  expect_equal(out$param_map$fe_names, c("b_Intercept", "b_x"))
  expect_equal(out$param_map$vc_names, "sigma")
  expect_equal(out$param_map$dropped_names, "lp__")
  expect_identical(out$param_map$map_source, "auto_regex")
  expect_equal(dim(out$log_lik), c(8L, 15L))
  expect_equal(out$diagnostics$n_long, 15L)
  expect_equal(out$meta$n_fits, 1L)
  expect_true(out$meta$log_lik_extracted)
  expect_true(out$meta$log_lik_retained)
  expect_false(out$meta$vc_confirmatory_reporting_allowed)
  expect_true(out$meta$prior_policy$explicit_prior_supplied)
  expect_true(out$meta$prior_policy$explicit_prior_warning)
  expect_true(out$meta$prior_diagnostic$warn_explicit_prior)
  expect_identical(out$meta$prior_diagnostic$reason_code, "explicit_prior_warning")
  expect_false("non_flat_prior_warning" %in% names(out$meta$prior_policy))
  expect_false("warn_nonflat_prior" %in% names(out$meta$prior_policy))
  expect_identical(out$meta$param_map_source, "auto_regex")
  expect_equal(out$meta$dropped_draw_columns, "lp__")
  expect_true(any(grepl("Explicit priors", out$warnings)))
  expect_true(any(grepl("lp__", out$warnings)))
  expect_equal(out$provenance$function_name, "pv_stack_fit")
  expect_equal(out$provenance$pv_cols, c("PV1", "PV2", "PV3"))
  expect_equal(out$provenance$weight_col, "W")
  expect_match(out$provenance$long_data_hash, "^[0-9a-f]{8}$")
  expect_true(is.list(out$fit))
  expect_true(is.data.frame(out$prepared_data))
})

test_that("stack fit accepts explicit param_map for custom fixed-effect names", {
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  out <- pvstackr:::pv_stack_fit(
    data = stack_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2", "PV3"),
    control = pv_control(backend = "injected"),
    fit_function = fake_custom_stack_fit(record),
    draws_function = function(fit) fit$draws,
    param_map = list(
      fe_names = c("theta_Intercept", "theta_x"),
      vc_names = "tau"
    )
  )

  expect_equal(record$n, 1L)
  expect_equal(colnames(out$stacked_draws), c("theta_Intercept", "theta_x", "tau"))
  expect_equal(out$param_map$fe_names, c("theta_Intercept", "theta_x"))
  expect_equal(out$param_map$vc_names, "tau")
  expect_equal(out$param_map$dropped_names, "lp__")
  expect_identical(out$param_map$map_source, "explicit")
  expect_equal(names(out$psi_hat_fe), c("theta_Intercept", "theta_x"))
  expect_true(any(grepl("explicit", out$warnings)))
  expect_true(any(grepl("lp__", out$warnings)))
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(out))
})

test_that("stack fit keeps heavy fields light by default", {
  out <- pvstackr:::pv_stack_fit(
    data = stack_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    control = pv_control(backend = "injected"),
    fit_function = fake_stack_fit(new.env()),
    draws_function = function(fit) fit$draws,
    extract_log_lik = TRUE,
    log_lik_function = function(fit) fit$log_lik
  )

  expect_null(out$fit)
  expect_null(out$prepared_data)
  expect_null(out$log_lik)
  expect_true(out$meta$log_lik_extracted)
  expect_false(out$meta$log_lik_retained)
  expect_false(out$meta$prior_policy$explicit_prior_supplied)
  expect_false(out$meta$prior_diagnostic$warn_explicit_prior)
  expect_true(is.na(out$meta$prior_diagnostic$reason_code))
  expect_equal(out$param_map$dropped_names, "lp__")
  expect_true(any(grepl("lp__", out$warnings)))
  expect_false(out$provenance$backend_fit_retained)
  expect_false(out$provenance$log_lik_retained)
})

test_that("stack fit validates injected hooks and protected arguments", {
  data <- stack_fixture_data()
  prepared <- pvstackr:::pv_prepare_stack_data(data, OUTCOME ~ x, c("PV1", "PV2"))

  expect_error(
    pvstackr:::pv_stack_fit(data, OUTCOME ~ x, c("PV1", "PV2"), control = pv_control()),
    "fit_function"
  )
  expect_error(
    pvstackr:::pv_stack_fit(data, OUTCOME ~ x, c("PV1", "PV2"), fit_function = fake_stack_fit(new.env()), control = pv_control()),
    "draws_function"
  )
  expect_error(
    pvstackr:::pv_stack_build_fit_args(prepared, additional_args = list(data = data.frame(x = 1))),
    "protected"
  )
  expect_error(
    pvstackr:::pv_stack_build_fit_args(prepared, additional_args = list(a = 1, a = 2)),
    "unique"
  )
})
