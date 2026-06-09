test_that("pv_control returns a validated control object with conservative defaults", {
  control <- pv_control()

  expect_s3_class(control, "pvstackr_control")
  expect_identical(control$method, "stack_direct")
  expect_identical(control$backend, "none")
  expect_identical(control$chains, 4L)
  expect_identical(control$iter, 2000L)
  expect_identical(control$warmup, 1000L)
  expect_identical(control$cores, 1L)
  expect_null(control$seed)
  expect_identical(control$conf_level, 0.95)
  expect_identical(control$psis_k_threshold, 0.7)
  expect_identical(control$center, "target")
  expect_false(control$allow_target_nearpd)
  expect_true(control$return_draws)
  expect_false(control$keep_data)
  expect_false(control$keep_backend_fit)
  expect_false(control$keep_log_lik)
  expect_false(control$verbose)
  expect_output(print(control), "target repair: unsupported")
})

test_that("pv_control accepts only canonical public method identifiers", {
  for (method in c("stack_direct", "stack_psis", "per_pv")) {
    expect_identical(pv_control(method = method)$method, method)
  }

  for (method in c("pipeline_a", "pipeline_b", "pipeline_c_direct", "c_direct", "Pipeline A")) {
    expect_error(pv_control(method = method), "internal manuscript/pipeline label")
  }

  for (method in c("direct", "stack", "psis", "reference", "bad_method")) {
    expect_error(pv_control(method = method), "must be one of")
  }
})

test_that("pv_control validates scalar sampler controls", {
  expect_identical(pv_control(iter = 1201L)$warmup, 600L)
  expect_identical(pv_control(iter = 1201L, warmup = 500L)$warmup, 500L)
  expect_identical(pv_control(seed = 123)$seed, 123L)

  expect_error(pv_control(chains = 0L), "`chains`")
  expect_error(pv_control(chains = c(1L, 2L)), "`chains`")
  expect_error(pv_control(iter = 1L), "`iter`")
  expect_error(pv_control(iter = 100L, warmup = 100L), "`warmup` must be smaller")
  expect_error(pv_control(warmup = -1L), "`warmup`")
  expect_error(pv_control(cores = 0L), "`cores`")
  expect_error(pv_control(seed = NA), "`seed`")
  expect_error(pv_control(seed = 1.5), "`seed`")
})

test_that("pv_control validates probabilities and choices", {
  expect_identical(pv_control(backend = "injected")$backend, "injected")
  expect_identical(pv_control(backend = "brms")$backend, "brms")
  expect_identical(pv_control(backend = "cmdstanr")$backend, "cmdstanr")
  expect_identical(pv_control(center = "posterior")$center, "posterior")

  expect_error(pv_control(backend = "stan"), "`backend`")
  expect_error(pv_control(center = "mean"), "`center`")
  expect_error(pv_control(conf_level = 0), "`conf_level`")
  expect_error(pv_control(conf_level = 1), "`conf_level`")
  expect_error(pv_control(conf_level = c(0.9, 0.95)), "`conf_level`")
  expect_error(pv_control(psis_k_threshold = 0), "`psis_k_threshold`")
  expect_error(pv_control(psis_k_threshold = Inf), "`psis_k_threshold`")
})

test_that("pv_control validates retention and verbosity flags strictly", {
  expect_false(pv_control(return_draws = FALSE)$return_draws)
  expect_true(pv_control(keep_data = TRUE)$keep_data)
  expect_true(pv_control(keep_backend_fit = TRUE)$keep_backend_fit)
  expect_true(pv_control(keep_log_lik = TRUE)$keep_log_lik)
  expect_true(pv_control(verbose = TRUE)$verbose)

  expect_error(pv_control(allow_target_nearpd = c(TRUE, FALSE)), "`allow_target_nearpd`")
  expect_error(pv_control(return_draws = NA), "`return_draws`")
  expect_error(pv_control(keep_data = 1), "`keep_data`")
  expect_error(pv_control(keep_backend_fit = "yes"), "`keep_backend_fit`")
  expect_error(pv_control(keep_log_lik = NULL), "`keep_log_lik`")
  expect_error(pv_control(verbose = c(TRUE, FALSE)), "`verbose`")
})

test_that("pv_control reserves unsupported target repair explicitly", {
  expect_error(
    pv_control(allow_target_nearpd = TRUE),
    "not currently supported"
  )
})

test_that("pv_validate_control revalidates corrupted control values", {
  control <- pv_control()

  bad <- control
  bad$chains <- 0L
  expect_error(pvstackr:::pv_validate_control(bad), "`chains`")

  bad <- control
  bad$method <- "bad"
  expect_error(pvstackr:::pv_validate_control(bad), "`method`")

  bad <- control
  bad$conf_level <- 2
  expect_error(pvstackr:::pv_validate_control(bad), "`conf_level`")

  bad <- control
  bad$psis_k_threshold <- 0
  expect_error(pvstackr:::pv_validate_control(bad), "`psis_k_threshold`")

  bad <- control
  bad$warmup <- bad$iter
  expect_error(pvstackr:::pv_validate_control(bad), "`warmup` must be smaller")

  bad <- control
  bad$allow_target_nearpd <- TRUE
  expect_error(pvstackr:::pv_validate_control(bad), "not currently supported")

  bad <- control
  bad$return_draws <- NA
  expect_error(pvstackr:::pv_validate_control(bad), "`return_draws`")
})

test_that("internal validators catch malformed columns and weights", {
  data <- data.frame(
    pv1 = c(1, 2, 3),
    pv2 = c(2, 3, 4),
    w = c(1, 1, 2),
    x = letters[1:3]
  )

  expect_identical(pvstackr:::pv_validate_pv_columns(data, c("pv1", "pv2")), c("pv1", "pv2"))
  expect_identical(pvstackr:::pv_validate_weight_vector(data$w, "w", nrow(data)), data$w)
  expect_identical(pvstackr:::pv_validate_fay_k(0.5), 0.5)

  expect_error(pvstackr:::pv_validate_pv_columns(data, "pv1"), "at least two")
  expect_error(pvstackr:::pv_validate_pv_columns(data, c("pv1", "missing")), "not found")
  expect_error(pvstackr:::pv_validate_pv_columns(data, c("pv1", "pv1")), "unique")
  expect_error(pvstackr:::pv_validate_pv_columns(data, c("pv1", "x")), "finite and numeric")
  expect_error(pvstackr:::pv_validate_weight_vector(c(1, 0, 1), "w", 3), "strictly positive")
  expect_error(pvstackr:::pv_validate_weight_vector(c(1, NA, 1), "w", 3), "finite numeric")
  expect_error(pvstackr:::pv_validate_fay_k(-0.1), "`fay_k`")
  expect_error(pvstackr:::pv_validate_fay_k(1), "`fay_k`")
})

test_that("internal target validation blocks unsafe covariance targets", {
  target <- diag(c(1, 2))
  rownames(target) <- colnames(target) <- c("alpha", "beta")

  expect_identical(
    pvstackr:::pv_validate_target_matrix(target, c("alpha", "beta")),
    target
  )
  expect_error(
    pvstackr:::pv_validate_target_matrix(target, c("alpha", "beta"), allow_nearpd = TRUE),
    "not currently supported"
  )
  expect_error(
    pvstackr:::pv_validate_target_matrix(target, c("alpha", "beta"), allow_nearpd = NA),
    "`allow_nearpd`"
  )

  unnamed <- diag(2)
  expect_error(pvstackr:::pv_validate_target_matrix(unnamed, c("alpha", "beta")), "row and column names")

  misordered <- target
  rownames(misordered) <- colnames(misordered) <- c("beta", "alpha")
  expect_error(pvstackr:::pv_validate_target_matrix(misordered, c("alpha", "beta")), "exactly match")

  nonsymmetric <- target
  nonsymmetric[1, 2] <- 0.5
  expect_error(pvstackr:::pv_validate_target_matrix(nonsymmetric, c("alpha", "beta")), "symmetric")

  non_pd <- target
  non_pd[2, 2] <- -1
  expect_error(pvstackr:::pv_validate_target_matrix(non_pd, c("alpha", "beta")), "positive definite")
  expect_error(
    pvstackr:::pv_validate_target_matrix(non_pd, c("alpha", "beta"), allow_nearpd = TRUE),
    "not currently supported"
  )

  mock <- target
  attr(mock, "target_source") <- "mock fallback target"
  expect_error(pvstackr:::pv_validate_target_matrix(mock, c("alpha", "beta")), "forbidden")

  expect_error(
    pvstackr:::pv_validate_target_matrix(2 * target, c("alpha", "beta"), raw_fe_cov = target),
    "proportional"
  )
})
