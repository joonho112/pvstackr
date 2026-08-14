skip_if_live_brms_smoke_disabled <- function() {
  skip_if_not(
    identical(Sys.getenv("PVSTACKR_RUN_LIVE_BACKEND_TESTS"), "true"),
    paste(
      "live bundled-brms smoke is dev-only; set",
      "PVSTACKR_RUN_LIVE_BACKEND_TESTS=true"
    )
  )
}

live_brms_require_real_package <- function(package) {
  desc <- suppressWarnings(utils::packageDescription(package))
  if (length(desc) == 1L && is.na(desc)) {
    stop(sprintf("Explicit live brms smoke requires `%s`.", package))
  }
  if (identical(desc$Title, "Sentinel Optional Package") ||
      identical(desc$Version, "0.0.0")) {
    stop(sprintf(
      "Explicit live brms smoke requires the real `%s` package, not a sentinel.",
      package
    ))
  }
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Explicit live brms smoke could not load `%s`.", package))
  }
  invisible(TRUE)
}

live_brms_toolchain_record <- function() {
  for (package in c("brms", "posterior", "cmdstanr")) {
    live_brms_require_real_package(package)
  }
  path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(error) stop(
      "Explicit live brms smoke requires configured CmdStan: ",
      conditionMessage(error)
    )
  )
  version <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
  if (!is.character(path) || length(path) != 1L || !dir.exists(path) ||
      length(version) != 1L || is.na(version)) {
    stop("Explicit live brms smoke requires a valid CmdStan path and version.")
  }
  tryCatch(
    cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE),
    error = function(error) stop(
      "Explicit live brms smoke requires a working CmdStan toolchain: ",
      conditionMessage(error)
    )
  )
  list(
    cmdstan_path_basename = basename(path),
    cmdstan_version = as.character(version),
    toolchain_check = "passed"
  )
}

live_brms_fixture <- function() {
  x <- seq(-1, 1, length.out = 8L)
  out <- data.frame(
    id = seq_along(x),
    x = x,
    PV1 = 1 + 0.5 * x + c(-0.10, 0.05, 0.02, -0.03, 0.04, -0.02, 0.03, 0.01),
    PV2 = 1.1 + 0.45 * x + c(0.02, -0.04, 0.03, 0.01, -0.02, 0.05, -0.03, 0.02),
    w = c(1, 1.1, 0.9, 1.2, 1, 0.8, 1.15, 0.95)
  )
  out$rw1 <- out$w * c(0.8, 1.2, 0.9, 1.1, 0.85, 1.15, 0.95, 1.05)
  out$rw2 <- out$w * c(1.2, 0.8, 1.1, 0.9, 1.15, 0.85, 1.05, 0.95)
  out$rw3 <- out$w * c(0.9, 1.1, 0.8, 1.2, 0.95, 1.05, 0.85, 1.15)
  out$rw4 <- out$w * c(1.1, 0.9, 1.2, 0.8, 1.05, 0.95, 1.15, 0.85)
  out
}

live_brms_capture_warnings <- function(code) {
  warnings <- character()
  value <- withCallingHandlers(
    code,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

live_brms_fit_once <- function(
  data,
  formula,
  cache_dir,
  model_bundle,
  binding_manifest
) {
  live_brms_capture_warnings(
    pvstackr:::pv_stack_fit(
      data = data,
      formula = formula,
      pv_cols = c("PV1", "PV2"),
      weight_col = "w",
      family = stats::gaussian(),
      control = pv_control(
        backend = "brms",
        chains = 1L,
        iter = 40L,
        warmup = 20L,
        cores = 1L,
        seed = 20260713L,
        return_draws = TRUE,
        keep_data = FALSE,
        keep_backend_fit = FALSE,
        verbose = FALSE
      ),
      cache_dir = cache_dir,
      cache_stem = "tiny-live-smoke",
      additional_args = list(refresh = 0L, silent = 2L),
      resolved_model_bundle = model_bundle,
      resolved_binding_manifest = binding_manifest
    )
  )
}

live_brms_execution_evidence <- function(fit, toolchain, warnings) {
  list(
    diagnostic_role = "execution_smoke",
    statistical_acceptance = FALSE,
    sampler_thresholds_applied = FALSE,
    observed_chains = fit$diagnostics$sampler$chains,
    observed_post_warmup_draws_per_chain =
      fit$diagnostics$sampler$post_warmup_draws_per_chain,
    toolchain_check = toolchain$toolchain_check,
    warning_count = as.integer(length(warnings))
  )
}

test_that("bundled brms compiles, samples, caches, and extracts a bounded draw shape", {
  skip_if_live_brms_smoke_disabled()
  toolchain <- live_brms_toolchain_record()
  data <- live_brms_fixture()
  analysis_formula <- stats::as.formula("OUTCOME ~ x", env = baseenv())
  target <- pv_brr_target(
    data = data,
    formula = analysis_formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "w",
    rep_weight_cols = paste0("rw", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )
  resolved <- pvstackr:::pv_stack_direct_preflight(
    data = data,
    formula = analysis_formula,
    target = target,
    family = stats::gaussian(),
    return_model_bundle = TRUE
  )
  cache_dir <- tempfile("pvstackr-live-brms-")
  expect_false(dir.exists(cache_dir))

  first <- live_brms_fit_once(
    data,
    analysis_formula,
    cache_dir,
    resolved$model_bundle,
    target$binding_manifest
  )
  fit <- first$value
  allowed_warning <- paste(
    "diverg|ESS has been capped|effective sample size|maximum treedepth|R-hat|Rhat",
    collapse = "|"
  )
  unknown_warnings <- first$warnings[
    !grepl(allowed_warning, first$warnings, ignore.case = TRUE)
  ]
  expect_length(unknown_warnings, 0L)

  expect_s3_class(fit, "pvstackr_stack_fit")
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(fit))
  expect_identical(fit$meta$engine_id, "bundled_brms_cmdstanr")
  expect_identical(fit$meta$fit_engine, "bundled_brms")
  expect_identical(fit$meta$adapter_source, "bundled")
  expect_identical(fit$meta$resolved_backend, "cmdstanr")
  expect_identical(fit$meta$backend_selection_reason, "configured_cmdstan_selected")
  expect_identical(fit$meta$cache$enabled, TRUE)
  expect_identical(fit$meta$cache$policy, "bundled_brms_managed")
  expect_identical(fit$meta$cache$cache_stem, "tiny-live-smoke")
  expect_identical(fit$meta$cache$file_refit, "on_change")
  expect_identical(fit$meta$cache$directory_created, TRUE)
  expect_identical(fit$meta$cache$writable, TRUE)
  expect_identical(fit$meta$model_matrix_materialized, TRUE)
  expect_identical(
    fit$meta$model_matrix_bundle_hash,
    target$binding_manifest$model_bundle_hash
  )
  expect_null(fit$fit)
  expect_null(fit$prepared_data)
  expect_match(fit$formula_string, "weights(.pvstackr_weight)", fixed = TRUE)
  expect_match(
    fit$formula_string,
    "0 + pvstackrMM001 + pvstackrMM002",
    fixed = TRUE
  )
  expect_identical(fit$weight_summary$M, 2L)
  expect_identical(fit$weight_summary$n_original, 8L)
  expect_identical(fit$weight_summary$n_long, 16L)
  expect_identical(fit$weight_summary$fractional_weight, 0.5)
  expect_identical(fit$weight_summary$weight_source, "w")
  expect_identical(fit$weight_summary$weight_col, "w")
  expect_equal(fit$weight_summary$mean_long_weight, 0.5, tolerance = 1e-14)
  expect_equal(fit$weight_summary$total_long_weight, 8, tolerance = 1e-14)
  expect_equal(
    unname(fit$weight_summary$per_pv_weight_sum),
    c(4, 4),
    tolerance = 1e-14
  )

  expect_identical(dim(fit$stacked_draws), c(20L, 3L))
  expect_identical(colnames(fit$stacked_draws), c("b_Intercept", "b_x", "sigma"))
  expect_true(all(is.finite(fit$stacked_draws)))
  expect_identical(fit$param_map$fe_names, c("b_Intercept", "b_x"))
  expect_identical(fit$param_map$vc_names, "sigma")
  expect_equal(
    fit$psi_hat_fe,
    colMeans(fit$stacked_draws[, c("b_Intercept", "b_x"), drop = FALSE]),
    tolerance = 0
  )
  expect_identical(
    fit$diagnostics$sampler$diagnostic_source,
    "bundled_brms_posterior_and_nuts"
  )
  expect_false(
    "diagnostic_extraction_failed" %in%
      fit$diagnostics$sampler$diagnostic_reason_codes
  )
  expect_identical(fit$diagnostics$sampler$chains, 1L)
  expect_identical(fit$diagnostics$sampler$post_warmup_draws_per_chain, 20L)

  engine <- fit$provenance$engine
  expect_identical(engine$adapter_source, "bundled")
  expect_identical(engine$requested_backend, "brms")
  expect_identical(engine$resolved_backend, "cmdstanr")
  expect_identical(engine$package_versions$brms, as.character(packageVersion("brms")))
  expect_identical(
    engine$package_versions$posterior,
    as.character(packageVersion("posterior"))
  )
  expect_identical(
    engine$package_versions$cmdstanr,
    as.character(packageVersion("cmdstanr"))
  )
  expect_identical(engine$package_versions$cmdstan, toolchain$cmdstan_version)
  expect_identical(
    engine$cmdstan_state$cmdstan_path_basename,
    toolchain$cmdstan_path_basename
  )
  retained_provenance_text <- paste(
    capture.output(str(list(meta = fit$meta, provenance = fit$provenance))),
    collapse = "\n"
  )
  expect_false(grepl(cache_dir, retained_provenance_text, fixed = TRUE))

  cache_file <- file.path(cache_dir, "tiny-live-smoke.rds")
  expect_true(file.exists(cache_file))
  backend_fit <- readRDS(cache_file)
  expect_s3_class(backend_fit, "brmsfit")
  expect_identical(backend_fit$family$family, "gaussian")
  expect_identical(backend_fit$family$link, "identity")
  expect_match(
    paste(deparse(backend_fit$formula$formula), collapse = ""),
    "weights(.pvstackr_weight)",
    fixed = TRUE
  )
  cache_hash <- digest::digest(cache_file, algo = "sha256", file = TRUE)
  cache_size <- file.info(cache_file)$size
  sentinel_mtime <- Sys.time() - 3600
  expect_true(Sys.setFileTime(cache_file, sentinel_mtime))
  cache_mtime <- file.info(cache_file)$mtime
  second <- live_brms_fit_once(
    data,
    analysis_formula,
    cache_dir,
    resolved$model_bundle,
    target$binding_manifest
  )
  second_unknown_warnings <- second$warnings[
    !grepl(allowed_warning, second$warnings, ignore.case = TRUE)
  ]
  expect_length(second_unknown_warnings, 0L)
  expect_identical(second$value$meta$cache$directory_created, FALSE)
  expect_identical(second$value$stacked_draws, fit$stacked_draws)
  expect_identical(
    digest::digest(cache_file, algo = "sha256", file = TRUE),
    cache_hash
  )
  expect_identical(file.info(cache_file)$mtime, cache_mtime)
  expect_identical(file.info(cache_file)$size, cache_size)

  # This is execution-contract evidence only. With one chain and 20 retained
  # draws, sampler values are recorded but never used for numerical acceptance.
  smoke_evidence <- live_brms_execution_evidence(fit, toolchain, first$warnings)
  expect_identical(smoke_evidence$diagnostic_role, "execution_smoke")
  expect_identical(smoke_evidence$statistical_acceptance, FALSE)
  expect_identical(smoke_evidence$sampler_thresholds_applied, FALSE)
  expect_identical(smoke_evidence$observed_chains, 1L)
  expect_identical(smoke_evidence$observed_post_warmup_draws_per_chain, 20L)
  expect_identical(smoke_evidence$toolchain_check, "passed")
})
