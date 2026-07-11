skip_if_backend_smoke_disabled <- function() {
  skip_if_not(
    identical(Sys.getenv("PVSTACKR_RUN_BACKEND_TESTS"), "true"),
    "optional backend smoke checks are dev-only; set PVSTACKR_RUN_BACKEND_TESTS=true"
  )
}

backend_smoke_heavy_packages <- function() {
  c("brms", "cmdstanr", "posterior", "loo", "EdSurvey", "haven", "lme4", "rstan", "StanHeaders")
}

backend_smoke_fit_function <- function(target, record, draws = NULL) {
  draws_template <- draws
  function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    record$called <- (if (is.null(record$called)) 0L else record$called) + 1L
    record$backend <- backend
    record$family_class <- class(family)
    record$prior_class <- class(prior)
    record$chains <- chains
    record$iter <- iter
    record$warmup <- warmup
    record$seed <- seed
    record$extra <- list(...)
    out_draws <- draws_template
    if (is.null(out_draws)) {
      out_draws <- pisa_tiny_parity_fit_function(target)(
        formula = formula,
        data = data,
        family = family,
        prior = prior,
        chains = chains,
        iter = iter,
        warmup = warmup,
        cores = cores,
        seed = seed,
        backend = backend,
        file = file,
        file_refit = file_refit
      )$draws
    }
    record$draws_class <- class(out_draws)
    list(draws = out_draws, backend = backend, family = family, prior = prior)
  }
}

backend_smoke_direct_fit <- function(backend = "injected", family = NULL, prior = NULL, draws = NULL, ...) {
  bundle <- pisa_tiny_parity_load()
  record <- new.env(parent = emptyenv())
  fit <- pv_fit_direct(
    data = bundle$data,
    formula = OUTCOME ~ x + female,
    target = bundle$cached$target,
    control = pv_control(
      method = "stack_direct",
      backend = backend,
      iter = 20L,
      warmup = 10L,
      chains = 2L,
      seed = 20260607L,
      return_draws = FALSE
    ),
    family = family,
    prior = prior,
    fit_function = backend_smoke_fit_function(bundle$cached$target, record, draws = draws),
    draws_function = function(fit) fit$draws,
    cache_dir = tempdir(),
    cache_stem = paste0("optional-backend-", backend),
    additional_args = list(...)
  )
  list(fit = fit, record = record)
}

test_that("injected backend smoke path remains available without optional packages", {
  before <- loadedNamespaces()
  out <- backend_smoke_direct_fit(backend = "injected", sample_prior = "no")

  expect_s3_class(out$fit, "pvstackr_fit")
  expect_equal(out$fit$status, "ok")
  expect_equal(out$fit$method, "stack_direct")
  expect_equal(out$record$called, 1L)
  expect_equal(out$record$backend, "injected")
  expect_equal(out$record$extra$sample_prior, "no")
  expect_null(get_draws(out$fit))
  expect_equal(intersect(setdiff(loadedNamespaces(), before), backend_smoke_heavy_packages()), character())
})

test_that("brms backend smoke uses injected adapter boundary without sampling", {
  skip_if_backend_smoke_disabled()
  skip_if_not_installed("brms")

  family <- brms::brmsfamily("gaussian")
  out <- backend_smoke_direct_fit(backend = "brms", family = family, sample_prior = "no")

  expect_s3_class(out$fit, "pvstackr_fit")
  expect_equal(out$fit$status, "ok")
  expect_equal(out$record$backend, "brms")
  expect_true("brmsfamily" %in% out$record$family_class)
  expect_equal(out$record$extra$sample_prior, "no")
  expect_equal(out$fit$stack_fit$meta$fit_engine, "injected_fit_function")
})

test_that("cmdstanr backend smoke uses injected adapter boundary without CmdStan", {
  skip_if_backend_smoke_disabled()
  skip_if_not_installed("cmdstanr")

  loadNamespace("cmdstanr")
  out <- backend_smoke_direct_fit(backend = "cmdstanr", cmdstan_model = "not-used")

  expect_s3_class(out$fit, "pvstackr_fit")
  expect_equal(out$fit$status, "ok")
  expect_equal(out$record$backend, "cmdstanr")
  expect_equal(out$record$extra$cmdstan_model, "not-used")
  expect_equal(out$fit$stack_fit$meta$fit_engine, "injected_fit_function")
})

test_that("posterior draws_matrix output is accepted by the injected draw extractor", {
  skip_if_backend_smoke_disabled()
  skip_if_not_installed("posterior")

  bundle <- pisa_tiny_parity_load()
  raw_draws <- pisa_tiny_parity_fit_function(bundle$cached$target)(
    formula = OUTCOME ~ x + female,
    data = bundle$data,
    family = NULL,
    prior = NULL,
    chains = 2L,
    iter = 20L,
    warmup = 10L,
    cores = 1L,
    seed = 20260607L,
    backend = "injected",
    file = tempfile("posterior-smoke-"),
    file_refit = "never"
  )$draws
  posterior_draws <- posterior::as_draws_matrix(raw_draws)
  record <- new.env(parent = emptyenv())
  fit <- pvstackr:::pv_stack_fit(
    data = bundle$data,
    formula = OUTCOME ~ x + female,
    pv_cols = bundle$cached$target$pv_cols,
    weight_col = bundle$cached$target$weight_col,
    control = pv_control(backend = "injected", iter = 20L, warmup = 10L, chains = 2L, seed = 20260607L),
    fit_function = backend_smoke_fit_function(bundle$cached$target, record, draws = posterior_draws),
    draws_function = function(fit) fit$draws
  )

  expect_s3_class(fit, "pvstackr_stack_fit")
  expect_true("draws_matrix" %in% record$draws_class)
  expect_equal(fit$param_map$fe_names, bundle$cached$target$fe_names)
  expect_equal(fit$param_map$vc_names, "sigma")
  expect_equal(colnames(fit$stacked_draws), c(bundle$cached$target$fe_names, "sigma"))
  expect_equal(fit$meta$fit_engine, "injected_fit_function")
})

test_that("posterior draws_matrix output survives the full stack_direct CCC path", {
  skip_if_backend_smoke_disabled()
  skip_if_not_installed("posterior")

  bundle <- pisa_tiny_parity_load()
  raw_draws <- pisa_tiny_parity_fit_function(bundle$cached$target)(
    formula = OUTCOME ~ x + female,
    data = bundle$data,
    family = NULL,
    prior = NULL,
    chains = 2L,
    iter = 20L,
    warmup = 10L,
    cores = 1L,
    seed = 20260607L,
    backend = "injected",
    file = tempfile("posterior-ccc-"),
    file_refit = "never"
  )$draws

  # Baseline plain-matrix run and the posterior::draws_matrix run must agree.
  # The draws_matrix source used to fail CCC validation: draws_calibrated stayed
  # a draws_matrix while draws_fe_cal became a plain matrix built by matrix
  # algebra, so the strict identical() block check compared class, not values.
  plain <- backend_smoke_direct_fit(backend = "injected", draws = raw_draws)
  posterior_run <- backend_smoke_direct_fit(
    backend = "injected",
    draws = posterior::as_draws_matrix(raw_draws)
  )

  expect_true("draws_matrix" %in% posterior_run$record$draws_class)
  expect_s3_class(posterior_run$fit, "pvstackr_fit")
  expect_equal(posterior_run$fit$status, "ok")
  expect_invisible(pvstackr:::validate_pvstackr_fit(posterior_run$fit))
  expect_identical(class(posterior_run$fit$ccc$draws_calibrated), c("matrix", "array"))
  expect_identical(class(posterior_run$fit$ccc$draws_fe_cal), c("matrix", "array"))
  expect_equal(posterior_run$fit$ccc$draws_calibrated, plain$fit$ccc$draws_calibrated, tolerance = 0)
  expect_equal(posterior_run$fit$estimates, plain$fit$estimates, tolerance = 0)
})
