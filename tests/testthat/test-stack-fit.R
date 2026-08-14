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

test_that("materialized stack data evaluates factor interactions and offsets without user closures", {
  data <- stack_fixture_data()
  data$g <- factor(rep(c("a", "b"), length.out = nrow(data)))
  formula <- OUTCOME ~ x * g + offset(z)
  bundle <- pvstackr:::pv_binding_resolve_model_bundle(data, formula)
  prepared <- pvstackr:::pv_prepare_stack_data(
    data = data,
    formula = formula,
    pv_cols = c("PV1", "PV2", "PV3"),
    weight_col = "W",
    model_bundle = bundle
  )

  terms <- stats::delete.response(stats::terms(prepared$formula))
  frame <- stats::model.frame(
    terms,
    data = prepared$data,
    na.action = stats::na.fail
  )
  backend_matrix <- stats::model.matrix(terms, data = frame)
  backend_offset <- stats::model.offset(frame)
  row_index <- rep(seq_len(nrow(data)), times = 3L)

  expect_identical(environment(prepared$formula), asNamespace("stats"))
  expected_matrix <- unname(bundle$model_matrix[row_index, , drop = FALSE])
  expect_identical(dim(backend_matrix), dim(expected_matrix))
  expect_equal(
    as.vector(backend_matrix),
    as.vector(expected_matrix),
    tolerance = 0
  )
  expect_identical(
    unname(backend_offset),
    unname(bundle$offset[row_index])
  )
  expect_identical(prepared$weight_summary$model_matrix_materialized, TRUE)
  expect_identical(
    prepared$weight_summary$model_matrix_bundle_hash,
    bundle$bundle_hash
  )
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
  expect_equal(
    map$dropped_names,
    c(
      "cor_school__Intercept__x",
      "<redacted_draw_column_002>",
      "lp__"
    )
  )
  expect_identical(map$map_source, "auto_regex")
  expect_match(map$warnings, "automatic regex")
  expect_match(map$warnings, "lp__")
  expect_false(grepl("r_school[", map$warnings, fixed = TRUE))

  privacy_draws <- cbind(draws, seq_len(nrow(draws)))
  colnames(privacy_draws)[ncol(privacy_draws)] <-
    "EXTERNAL_RAW_DATA_SENTINEL"
  privacy_map <- pvstackr:::pv_stack_param_map(privacy_draws)
  expect_true("<redacted_draw_column_004>" %in% privacy_map$dropped_names)
  expect_false(any(grepl(
    "EXTERNAL_RAW_DATA_SENTINEL",
    c(privacy_map$dropped_names, privacy_map$warnings),
    fixed = TRUE
  )))
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
  level_draws <- cbind(draws, `r_school[CONFIDENTIAL_ID,Intercept]` = 1:5)
  expect_error(
    pvstackr:::pv_stack_param_map(
      level_draws,
      param_map = list(
        fe_names = c("theta_Intercept", "theta_x"),
        vc_names = "r_school[CONFIDENTIAL_ID,Intercept]"
      )
    ),
    "level-specific `r_` draws"
  )
})

test_that("stack cache specification supports managed cache and explicit opt-out", {
  disabled <- pvstackr:::pv_stack_cache_spec(cache_dir = NULL, package_managed = TRUE)
  expect_false(disabled$enabled)
  expect_null(disabled$file)
  expect_identical(disabled$file_refit, "never")
  expect_identical(
    pvstackr:::pv_stack_cache_provenance(disabled)$policy,
    "disabled"
  )

  parent <- tempfile("pvstackr-cache-parent-")
  managed_dir <- file.path(parent, "nested", "cache")
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  managed <- pvstackr:::pv_stack_cache_spec(
    cache_dir = managed_dir,
    cache_stem = "managed-fit",
    package_managed = TRUE
  )
  expect_true(dir.exists(managed_dir))
  expect_true(managed$directory_created)
  expect_true(managed$writable)
  expect_identical(managed$file_refit, "on_change")
  expect_match(managed$file, "managed-fit$")
  expect_length(list.files(managed_dir, all.files = TRUE, no.. = TRUE), 0L)
  expect_identical(
    pvstackr:::pv_stack_cache_provenance(managed)$policy,
    "bundled_brms_managed"
  )
})

test_that("managed stack cache rejects unsafe stems and file paths", {
  expect_error(
    pvstackr:::pv_stack_cache_spec(tempdir(), "../escape", package_managed = TRUE),
    "safe file stem"
  )

  not_a_dir <- tempfile("pvstackr-cache-file-")
  on.exit(unlink(not_a_dir, force = TRUE), add = TRUE)
  expect_true(file.create(not_a_dir))
  expect_error(
    pvstackr:::pv_stack_cache_spec(not_a_dir, package_managed = TRUE),
    "not a directory"
  )

  collision_dir <- tempfile("pvstackr-cache-collision-")
  on.exit(unlink(collision_dir, recursive = TRUE, force = TRUE), add = TRUE)
  expect_true(dir.create(collision_dir))
  expect_true(dir.create(file.path(collision_dir, "collision.rds")))
  expect_error(
    pvstackr:::pv_stack_cache_spec(
      collision_dir,
      cache_stem = "collision",
      package_managed = TRUE
    ),
    "collides with an existing directory"
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
  expect_false("n_long" %in% names(out$diagnostics))
  expect_false(out$diagnostics$sampler$diagnostic_complete)
  expect_identical(
    out$diagnostics$sampler$diagnostic_source,
    "injected_diagnose_function"
  )
  expect_true(any(grepl(
    "diagnostic_missing_rhat_max",
    out$diagnostics$sampler$diagnostic_reason_codes,
    fixed = TRUE
  )))
  expect_equal(out$meta$n_fits, 1L)
  expect_identical(out$meta$engine_id, "injected_fit_function")
  expect_identical(out$meta$fit_engine, "injected_fit_function")
  expect_identical(out$meta$adapter_source, "injected")
  expect_identical(out$meta$resolved_backend, "injected")
  expect_identical(out$meta$cache$policy, "injected_adapter_managed")
  expect_identical(out$meta$cache$file_refit, "on_change")
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
  expect_identical(out$provenance$engine$requested_backend, "injected")
  expect_identical(out$provenance$engine$resolved_backend, "injected")
  expect_identical(out$provenance$cache, out$meta$cache)
  expect_true(is.list(out$fit))
  expect_true(is.data.frame(out$prepared_data))
})

test_that("low-level stack fit honors return_draws without losing summaries", {
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  out <- pvstackr:::pv_stack_fit(
    data = stack_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2", "PV3"),
    weight_col = "W",
    control = pv_control(
      backend = "injected",
      return_draws = FALSE
    ),
    fit_function = fake_stack_fit(record),
    draws_function = function(fit) fit$draws
  )

  expect_equal(record$n, 1L)
  expect_null(out$stacked_draws)
  expect_false(out$control$return_draws)
  expect_equal(names(out$psi_hat_fe), c("b_Intercept", "b_x"))
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(out))

  injected <- out
  injected$stacked_draws <- matrix(
    0,
    nrow = 2L,
    ncol = 3L,
    dimnames = list(NULL, c("b_Intercept", "b_x", "sigma"))
  )
  expect_error(
    pvstackr:::validate_pvstackr_stack_fit(injected),
    "draw retention"
  )
})

test_that("injected stack fit preserves backend label and cache opt-out", {
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  out <- pvstackr:::pv_stack_fit(
    data = stack_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    control = pv_control(backend = "cmdstanr"),
    fit_function = fake_stack_fit(record),
    draws_function = function(fit) fit$draws,
    cache_dir = NULL
  )

  expect_identical(record$args$backend, "cmdstanr")
  expect_null(record$args$file)
  expect_identical(record$args$file_refit, "never")
  expect_identical(out$provenance$engine$adapter_source, "injected")
  expect_identical(out$provenance$engine$requested_backend, "cmdstanr")
  expect_identical(out$provenance$engine$resolved_backend, "injected")
  expect_false(out$provenance$cache$enabled)
  expect_identical(out$provenance$cache$policy, "disabled")
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

test_that("materialized stack fit maps explicit custom backend names to bound names", {
  data <- stack_fixture_data()
  formula <- OUTCOME ~ x
  bundle <- pvstackr:::pv_binding_resolve_model_bundle(data, formula)
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  out <- pvstackr:::pv_stack_fit(
    data = data,
    formula = formula,
    pv_cols = c("PV1", "PV2", "PV3"),
    control = pv_control(backend = "injected"),
    fit_function = fake_custom_stack_fit(record),
    draws_function = function(fit) fit$draws,
    param_map = list(
      fe_names = c("theta_Intercept", "theta_x"),
      vc_names = "tau"
    ),
    resolved_model_bundle = bundle
  )

  expect_identical(
    colnames(out$stacked_draws),
    c("b_Intercept", "b_x", "tau")
  )
  expect_identical(out$param_map$fe_names, c("b_Intercept", "b_x"))
  expect_identical(out$param_map$vc_names, "tau")
  expect_identical(out$param_map$map_source, "explicit")
  expect_identical(out$param_map$original_fe_idx, c(1L, 2L))
  expect_identical(out$meta$model_matrix_bundle_hash, bundle$bundle_hash)
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(out))

  bad_record <- new.env(parent = emptyenv())
  bad_record$n <- 0L
  expect_error(
    pvstackr:::pv_stack_fit(
      data = data,
      formula = formula,
      pv_cols = c("PV1", "PV2", "PV3"),
      control = pv_control(backend = "injected"),
      fit_function = fake_custom_stack_fit(bad_record),
      draws_function = function(fit) fit$draws,
      param_map = list(fe_names = "theta_Intercept", vc_names = "tau"),
      resolved_model_bundle = bundle
    ),
    "complete bound design width"
  )
})

test_that("materialized prior policy permits only provably invariant tables", {
  with_intercept <- stats::setNames(
    c("b_Intercept", "b_x"),
    c("b_pvstackrMM001", "b_pvstackrMM002")
  )
  no_intercept <- stats::setNames("b_x", "b_pvstackrMM001")
  sigma_prior <- data.frame(
    prior = "normal(0, 1)", class = "sigma", coef = "",
    stringsAsFactors = FALSE
  )
  global_b <- data.frame(
    prior = "normal(0, 1)", class = "b", coef = "",
    stringsAsFactors = FALSE
  )
  coefficient_b <- data.frame(
    prior = "normal(0, 1)", class = "b", coef = "x",
    stringsAsFactors = FALSE
  )
  intercept_prior <- data.frame(
    prior = "normal(0, 1)", class = "Intercept", coef = "",
    stringsAsFactors = FALSE
  )

  expect_identical(
    pvstackr:::pv_stack_materialized_prior(sigma_prior, with_intercept),
    sigma_prior
  )
  expect_identical(
    pvstackr:::pv_stack_materialized_prior(global_b, no_intercept),
    global_b
  )
  # A population-level prior is expanded onto the non-intercept columns rather
  # than refused: `0 + ...` on the materialized matrix would otherwise let it
  # reach the intercept column, which the original formula excludes.
  expanded <- pvstackr:::pv_stack_materialized_prior(global_b, with_intercept)
  expect_identical(expanded$class, "b")
  expect_identical(expanded$coef, "pvstackrMM002")
  expect_identical(expanded$prior, "normal(0, 1)")
  expect_false("pvstackrMM001" %in% expanded$coef)

  expect_error(
    pvstackr:::pv_stack_materialized_prior(coefficient_b, with_intercept),
    "cannot be preserved exactly"
  )
  expect_error(
    pvstackr:::pv_stack_materialized_prior(intercept_prior, with_intercept),
    "cannot be preserved exactly"
  )
  expect_error(
    pvstackr:::pv_stack_materialized_prior("opaque-prior", with_intercept),
    "opaque"
  )
})

test_that("population-level prior expansion preserves the original prior scope", {
  design <- stats::setNames(
    c("b_Intercept", "b_escs", "b_female"),
    c("b_pvstackrMM001", "b_pvstackrMM002", "b_pvstackrMM003")
  )
  global_b <- data.frame(
    prior = "normal(0, 25)", class = "b", coef = "",
    stringsAsFactors = FALSE
  )

  expanded <- pvstackr:::pv_stack_materialized_prior(global_b, design)

  # One row per slope, none for the intercept column, same distribution.
  expect_identical(nrow(expanded), 2L)
  expect_setequal(expanded$coef, c("pvstackrMM002", "pvstackrMM003"))
  expect_true(all(expanded$prior == "normal(0, 25)"))
  expect_true(all(expanded$class == "b"))

  # Other classes ride along untouched.
  mixed <- rbind(
    global_b,
    data.frame(prior = "student_t(3, 0, 10)", class = "sigma", coef = "",
               stringsAsFactors = FALSE)
  )
  out <- pvstackr:::pv_stack_materialized_prior(mixed, design)
  expect_identical(nrow(out), 3L)
  expect_identical(out$prior[out$class == "sigma"], "student_t(3, 0, 10)")
  expect_setequal(out$coef[out$class == "b"], c("pvstackrMM002", "pvstackrMM003"))

  # A scoped prior is still refused: expansion cannot preserve group scope.
  scoped <- cbind(global_b, group = "school", stringsAsFactors = FALSE)
  expect_error(
    pvstackr:::pv_stack_materialized_prior(scoped, design),
    "cannot be preserved exactly"
  )

  # An intercept-only design has nothing to carry a population-level prior.
  expect_error(
    pvstackr:::pv_stack_materialized_prior(
      global_b,
      stats::setNames("b_Intercept", "b_pvstackrMM001")
    ),
    "no non-intercept coefficient"
  )
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
  expect_false(out$diagnostics$sampler$diagnostic_complete)
  expect_identical(
    out$diagnostics$sampler$diagnostic_reason_codes,
    "diagnostic_extractor_not_supplied"
  )
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
