fit_direct_fixture_data <- function() {
  data.frame(
    id = seq_len(10),
    school = factor(rep(c("A", "B"), each = 5)),
    x = c(-2, -1.5, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3),
    z = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    PV1 = c(2.1, 2.4, 3.0, 3.2, 3.6, 4.0, 4.2, 4.8, 5.1, 5.5),
    PV2 = c(2.0, 2.5, 2.9, 3.4, 3.5, 4.1, 4.4, 4.7, 5.0, 5.8),
    W = c(1.0, 1.1, 0.9, 1.2, 1.0, 1.3, 0.8, 1.4, 1.1, 0.95),
    RW1 = c(0.8, 1.3, 1.0, 1.1, 0.9, 1.5, 0.9, 1.2, 1.2, 1.0),
    RW2 = c(1.2, 0.9, 0.8, 1.4, 1.1, 1.0, 0.7, 1.6, 1.0, 0.9),
    RW3 = c(1.1, 1.0, 1.1, 1.0, 1.3, 1.2, 0.8, 1.3, 1.4, 1.1),
    RW4 = c(0.9, 1.2, 0.95, 1.3, 1.0, 1.4, 1.1, 1.1, 0.9, 1.2)
  )
}

fit_direct_target <- function(formula = OUTCOME ~ x, data = fit_direct_fixture_data(), ...) {
  pv_brr_target(
    data = data,
    formula = formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    ...
  )
}

fake_fit_direct_fit <- function(record, target, center_shift_se = 0) {
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
    beta <- target$beta[target$fe_names]
    se <- sqrt(diag(target$T_MI))[target$fe_names]
    e1 <- c(-0.03, -0.01, 0.02, 0.04, -0.02, 0.01, 0.03, -0.04)
    e2 <- c(0.02, -0.03, 0.01, -0.04, 0.03, -0.01, 0.04, -0.02)
    draws <- cbind(
      b_Intercept = beta[["b_Intercept"]] + se[["b_Intercept"]] * (center_shift_se + e1),
      b_x = beta[["b_x"]] + se[["b_x"]] * (center_shift_se + e2),
      sigma = rep(1, 8),
      lp__ = seq(-10, -9, length.out = 8)
    )
    list(draws = draws, log_lik = matrix(0, nrow = 8, ncol = nrow(data)), data = data)
  }
}

fit_direct_control <- function(...) {
  pv_control(
    method = "stack_direct",
    backend = "injected",
    iter = 10L,
    warmup = 5L,
    chains = 2L,
    seed = 20260607L,
    ...
  )
}

fit_direct_conditioning_ccc <- function(scale) {
  raw_center <- c(b_Intercept = 0, b_x = 0)
  raw_centered <- rbind(
    c(sqrt(6), 0),
    c(-sqrt(6), 0),
    c(0, sqrt(13.5)),
    c(0, -sqrt(13.5))
  )
  colnames(raw_centered) <- names(raw_center)
  target <- structure(
    list(
      beta = raw_center,
      T_MI = diag(c(4 * scale^2, 9)),
      fe_names = names(raw_center),
      target_source = "external_brr_fay_rubin",
      target_hash = "fit-direct-conditioning"
    ),
    class = c("pvstackr_brr_target", "list")
  )
  dimnames(target$T_MI) <- list(names(raw_center), names(raw_center))
  pvstackr:::ccc_calibrate(raw_centered, target)
}

test_that("pv_fit_direct returns a calibrated reportable stack_direct fit", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(keep_data = TRUE, keep_log_lik = TRUE),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    diagnose_function = function(fit) list(n_long = nrow(fit$data)),
    extract_log_lik = TRUE,
    log_lik_function = function(fit) fit$log_lik,
    cache_dir = tempdir(),
    cache_stem = "fit-direct",
    additional_args = list(sample_prior = "no")
  )

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(record$n, 1L)
  expect_equal(fit$method, "stack_direct")
  expect_equal(fit$status, "ok")
  expect_equal(fit$reason_codes, character())
  expect_equal(fit$warnings, character())
  expect_equal(fit$target$target_hash, target$target_hash)
  expect_equal(fit$ccc$target_hash, target$target_hash)
  expect_equal(fit$ccc$target_source, "external_brr_fay_rubin")
  expect_equal(fit$stack_fit$param_map$fe_names, target$fe_names)
  expect_equal(fit$stack_fit$param_map$dropped_names, "lp__")
  expect_true(any(grepl("lp__", fit$stack_fit$warnings)))
  expect_equal(colnames(fit$draws), target$fe_names)
  expect_equal(fit$draws, fit$ccc$draws_fe_cal)
  expect_equal(colMeans(fit$draws), fit$estimates$estimate, ignore_attr = TRUE)
  expect_equal(fit$ccc$Sigma_cal_emp, target$T_MI, tolerance = 1e-8)
  expect_equal(fit$estimates$term, target$fe_names)
  expect_equal(nrow(fit$estimates), length(target$fe_names))
  expect_equal(fit$estimates$estimate, unname(target$beta[target$fe_names]))
  expect_equal(fit$estimates$se, sqrt(diag(target$T_MI))[target$fe_names], ignore_attr = TRUE)
  expect_equal(fit$estimates$std.error, fit$estimates$se)
  expect_equal(fit$estimates$df, unname(target$df[target$fe_names]))
  expect_equal(fit$estimates$df_method, rep("classic", length(target$fe_names)))
  expect_true(all(is.na(fit$estimates$df_complete)))
  expect_equal(fit$estimates$interval_role, rep(target$interval_role, length(target$fe_names)))
  expect_equal(fit$estimates$coverage_claim_allowed, rep(target$coverage_claim_allowed, length(target$fe_names)))
  expect_false(any(fit$estimates$coverage_claim_allowed))
  expect_equal(fit$estimates$parameter_scope, rep("fixed_effect", length(target$fe_names)))
  expect_equal(fit$estimates$target_source, rep("external_brr_fay_rubin", length(target$fe_names)))
  expect_equal(fit$estimates$target_hash, rep(target$target_hash, length(target$fe_names)))
  expect_equal(fit$estimates$conf_low, fit$estimates$conf.low)
  expect_equal(fit$estimates$conf_high, fit$estimates$conf.high)
  expect_equal(fit$diagnostics$preflight$target_hash, target$target_hash)
  expect_equal(fit$diagnostics$preflight$target_source, "external_brr_fay_rubin")
  expect_equal(fit$diagnostics$stack_fit$n_long, nrow(data) * 2L)
  expect_equal(fit$design$roles$method, "stack_direct")
  expect_equal(fit$design$provenance$target_hash, target$target_hash)
  expect_equal(fit$provenance$wrapper_function, "pv_fit_direct")
  expect_equal(fit$provenance$target_hash, target$target_hash)
  expect_equal(fit$provenance$ccc_target_hash, target$target_hash)
  expect_equal(fit$provenance$stack_fit_long_data_hash, fit$stack_fit$weight_summary$long_data_hash)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_direct inherits Barnard-Rubin interval policy from the external target", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(
    OUTCOME ~ x,
    data = data,
    df_method = "barnard_rubin",
    df_complete = c(b_x = 30, b_Intercept = 20)
  )
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(return_draws = FALSE),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws
  )

  expect_equal(fit$estimates$df_method, rep("barnard_rubin", length(target$fe_names)))
  expect_equal(fit$estimates$df_complete, c(20, 30), tolerance = 0)
  expect_equal(fit$estimates$interval_role, rep("coverage_barnard_rubin", length(target$fe_names)))
  expect_true(all(fit$estimates$coverage_claim_allowed))
  expect_equal(fit$estimates$df, unname(target$df[target$fe_names]), tolerance = 1e-12)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_direct honors return_draws retention policy", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(return_draws = FALSE),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws
  )

  expect_equal(fit$status, "ok")
  expect_null(fit$draws)
  expect_equal(fit$estimates$term, target$fe_names)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

# Mimic the object an `as_draws_matrix()` draws_function would return without
# taking a dependency on the posterior package (deliberately absent per the
# light-check contract): a matrix carrying the draws_matrix S3 class plus
# draw-id row names, so is.matrix() is TRUE but identical() differs from a
# plain matrix on class -- the exact shape that tripped the strict CCC check.
as_draws_matrix_like <- function(x) {
  structure(
    x,
    dimnames = list(as.character(seq_len(nrow(x))), colnames(x)),
    class = c("draws_matrix", "draws", "matrix", "array")
  )
}

test_that("pv_fit_direct accepts a draws_function returning a draws_matrix-classed object", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  plain_fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws
  )

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) as_draws_matrix_like(fit$draws)
  )

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$status, "ok")
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
  expect_identical(class(fit$draws), c("matrix", "array"))
  expect_identical(class(fit$ccc$draws_calibrated), c("matrix", "array"))
  expect_equal(fit$draws, plain_fit$draws, tolerance = 0)
  expect_equal(fit$estimates, plain_fit$estimates, tolerance = 0)
})

test_that("pv_fit_direct preflight rejects incompatible calls before fitting", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x + z,
      target = target,
      control = fit_direct_control(),
      fit_function = fake_fit_direct_fit(record, target),
      draws_function = function(fit) fit$draws
    ),
    "formula RHS"
  )
  expect_equal(record$n, 0L)

  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x + (1 | school),
      target = target,
      control = fit_direct_control(),
      fit_function = fake_fit_direct_fit(record, target),
      draws_function = function(fit) fit$draws
    ),
    "Random-effect/group terms"
  )
  expect_equal(record$n, 0L)
})

test_that("pv_fit_direct promotes explicit prior diagnostics to warning status", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    prior = "fake_prior",
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws
  )

  expect_equal(fit$status, "warning")
  expect_equal(fit$reason_codes, "explicit_prior_warning")
  expect_true(any(grepl("Explicit priors", fit$warnings)))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_direct_status normalizes missing explicit-prior reason codes", {
  ccc <- list(diagnostics = list(center_status = "ok", conditioning_status = "ok"), warnings = character())
  missing_reason_policy <- list(
    explicit_prior_warning = TRUE,
    reason_code = NA_character_,
    warning = "Explicit priors were supplied."
  )
  empty_reason_policy <- list(
    explicit_prior_warning = TRUE,
    reason_code = "",
    warning = "Explicit priors were supplied."
  )

  for (policy in list(missing_reason_policy, empty_reason_policy)) {
    out <- pvstackr:::pv_fit_direct_status(
      stack_fit = list(meta = list(prior_policy = policy)),
      ccc = ccc
    )
    expect_equal(out$status, "warning")
    expect_equal(out$reason_codes, "explicit_prior_warning")
    expect_equal(out$warnings, "Explicit priors were supplied.")
  }
})

test_that("pv_fit_direct promotes yellow center separation diagnostics", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target, center_shift_se = 0.02),
    draws_function = function(fit) fit$draws
  )

  expect_equal(fit$status, "warning")
  expect_equal(fit$reason_codes, "center_separation_yellow")
  expect_true(any(grepl("Center separation", fit$warnings)))
  expect_equal(fit$ccc$diagnostics$center_status, "warning")
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_direct blocks red center separation diagnostics", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target, center_shift_se = 0.08),
    draws_function = function(fit) fit$draws
  )

  expect_equal(fit$status, "blocked")
  expect_equal(fit$reason_codes, "center_separation_red")
  expect_true(any(grepl("Center separation", fit$warnings)))
  expect_equal(fit$ccc$diagnostics$center_status, "blocked")
  expect_equal(nrow(fit$estimates), 0L)
  expect_null(fit$draws)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_direct_status propagates CCC conditioning diagnostics", {
  warning_ccc <- fit_direct_conditioning_ccc(1e6)
  warning_status <- pvstackr:::pv_fit_direct_status(
    stack_fit = list(meta = list()),
    ccc = warning_ccc
  )

  expect_equal(warning_ccc$diagnostics$center_status, "ok")
  expect_equal(warning_ccc$diagnostics$conditioning_status, "warning")
  expect_equal(warning_status$status, "warning")
  expect_equal(warning_status$reason_codes, "ccc_conditioning_yellow")
  expect_true(any(grepl("CCC conditioning diagnostic", warning_status$warnings)))

  blocked_ccc <- fit_direct_conditioning_ccc(1e8)
  blocked_status <- pvstackr:::pv_fit_direct_status(
    stack_fit = list(meta = list()),
    ccc = blocked_ccc
  )

  expect_equal(blocked_ccc$diagnostics$center_status, "ok")
  expect_equal(blocked_ccc$diagnostics$conditioning_status, "blocked")
  expect_equal(blocked_status$status, "blocked")
  expect_equal(blocked_status$reason_codes, "ccc_conditioning_red")
  expect_true(any(grepl("reportable stack_direct estimates are blocked", blocked_status$warnings)))
})

test_that("pv_fit_direct enforces stack_direct controls", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)

  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x,
      target = target,
      control = pv_control(method = "per_pv"),
      fit_function = fake_fit_direct_fit(new.env(parent = emptyenv()), target),
      draws_function = function(fit) fit$draws
    ),
    "control\\$method"
  )

  posterior_record <- new.env(parent = emptyenv())
  posterior_record$n <- 0L
  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x,
      target = target,
      control = fit_direct_control(center = "posterior"),
      fit_function = fake_fit_direct_fit(posterior_record, target),
      draws_function = function(fit) fit$draws
    ),
    "control\\$center = \"target\""
  )
  expect_equal(posterior_record$n, 0L)

  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x,
      target = target,
      control = fit_direct_control(),
      draws_function = function(fit) fit$draws
    ),
    "fit_function"
  )
  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x,
      target = target,
      control = fit_direct_control(),
      fit_function = fake_fit_direct_fit(new.env(parent = emptyenv()), target)
    ),
    "draws_function"
  )
})

test_that("pv_fit dispatches stack_direct and guards unimplemented methods", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    control = fit_direct_control()
  )

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$method, "stack_direct")
  expect_equal(record$n, 1L)

  expect_error(
    pv_fit(data = data, formula = OUTCOME ~ x, control = fit_direct_control()),
    "`target` is required"
  )

  expect_error(
    pv_fit(
      data = data,
      formula = OUTCOME ~ x,
      target = target,
      method = "stack_direct",
      control = pv_control(method = "per_pv")
    ),
    "control\\$method"
  )
})

test_that("pv_fit defaults to stack_direct control when control is omitted", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws
  )

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$method, "stack_direct")
  expect_equal(fit$control$method, "stack_direct")
  expect_equal(fit$control$backend, "none")
  expect_equal(record$n, 1L)
  expect_equal(record$args$backend, "none")
  expect_equal(fit$estimates$target_hash, rep(target$target_hash, length(target$fe_names)))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})
