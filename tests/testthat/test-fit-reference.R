reference_fixture_data <- function() {
  data.frame(
    id = seq_len(6),
    school = factor(rep(c("A", "B", "C"), each = 2)),
    x = c(-1.5, -0.5, 0, 0.5, 1, 1.5),
    PV1 = c(1.0, 1.4, 1.8, 2.1, 2.5, 2.9),
    PV2 = c(1.1, 1.3, 1.9, 2.0, 2.7, 3.0),
    PV3 = c(0.9, 1.5, 1.7, 2.2, 2.6, 3.1),
    W = c(1.0, 1.1, 0.9, 1.2, 1.0, 1.3)
  )
}

reference_draw_matrix <- function(center) {
  cbind(
    b_Intercept = center[[1]] + c(-0.20, -0.05, 0.05, 0.20),
    b_x = center[[2]] + c(-0.10, 0.15, -0.15, 0.10),
    sigma = c(0.9, 1.0, 1.1, 1.2),
    lp__ = c(-5.0, -4.8, -4.7, -4.5)
  )
}

reference_draws_fixture <- function() {
  list(
    PV1 = reference_draw_matrix(c(1.0, 0.3)),
    PV2 = reference_draw_matrix(c(1.2, 0.5)),
    PV3 = reference_draw_matrix(c(1.4, 0.7))
  )
}

reference_pool_oracle <- function(draws, ...) {
  fe_draws <- lapply(draws, function(x) x[, c("b_Intercept", "b_x"), drop = FALSE])
  beta <- do.call(rbind, lapply(fe_draws, colMeans))
  U <- lapply(fe_draws, stats::cov)
  for (m in seq_along(U)) {
    dimnames(U[[m]]) <- list(colnames(beta), colnames(beta))
  }
  pvstackr:::rubin_pool_matrix(beta, U, orientation = "rows_pv", ...)
}

reference_control <- function(...) {
  pv_control(
    method = "per_pv",
    backend = "injected",
    iter = 10L,
    warmup = 5L,
    chains = 2L,
    seed = 20260607L,
    ...
  )
}

test_that("pv_fit_reference pools precomputed per-PV draws with Rubin rules", {
  draws <- reference_draws_fixture()
  oracle <- reference_pool_oracle(draws)

  fit <- pv_fit_reference(
    per_pv_draws = draws,
    control = reference_control(return_draws = TRUE)
  )

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$method, "per_pv")
  expect_equal(fit$status, "ok")
  expect_null(fit$stack_fit)
  expect_null(fit$ccc)
  expect_null(fit$draws)
  expect_s3_class(fit$target, "pvstackr_reference_pool")
  expect_equal(fit$target$target_source, "per_pv_rubin_draws")
  expect_equal(fit$target$beta, oracle$beta, tolerance = 1e-14)
  expect_equal(fit$target$T_MI, oracle$T_MI, tolerance = 1e-14)
  expect_equal(fit$target$df, oracle$df, tolerance = 1e-14)
  expect_equal(fit$estimates$term, names(oracle$beta))
  expect_equal(fit$estimates$estimate, unname(oracle$beta), tolerance = 1e-14)
  expect_equal(fit$estimates$se, unname(oracle$se), tolerance = 1e-14)
  expect_equal(fit$estimates$std.error, fit$estimates$se)
  expect_equal(fit$estimates$df, unname(oracle$df), tolerance = 1e-14)
  expect_equal(fit$estimates$conf_low, unname(oracle$ci_low), tolerance = 1e-14)
  expect_equal(fit$estimates$conf_high, unname(oracle$ci_high), tolerance = 1e-14)
  expect_equal(fit$estimates$conf.low, fit$estimates$conf_low)
  expect_equal(fit$estimates$conf.high, fit$estimates$conf_high)
  expect_equal(fit$estimates$interval_role, rep("reference_classic_rubin", 2L))
  expect_false(any(fit$estimates$coverage_claim_allowed))
  expect_equal(fit$estimates$parameter_scope, rep("fixed_effect", 2L))
  expect_equal(fit$estimates$target_source, rep("per_pv_rubin_draws", 2L))
  expect_equal(fit$estimates$pooling_source, rep("per_pv_rubin_draws", 2L))
  expect_equal(fit$estimates$pooling_hash, rep(fit$target$target_hash, 2L))
  expect_equal(names(fit$diagnostics$reference$per_pv_draws), names(draws))
  expect_true(all(vapply(
    fit$diagnostics$reference$per_pv_draws,
    function(x) identical(colnames(x), c("b_Intercept", "b_x")),
    logical(1)
  )))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  redacted <- pv_fit_reference(
    per_pv_draws = draws,
    control = reference_control(return_draws = FALSE)
  )
  expect_null(redacted$diagnostics$reference$per_pv_draws)
  expect_invisible(pvstackr:::validate_pvstackr_fit(redacted))

  injected <- redacted
  injected$diagnostics$reference$per_pv_draws <-
    fit$diagnostics$reference$per_pv_draws
  expect_error(pvstackr:::validate_pvstackr_fit(injected), "return_draws")

  removed <- fit
  removed$diagnostics$reference["per_pv_draws"] <- list(NULL)
  expect_error(pvstackr:::validate_pvstackr_fit(removed), "return_draws")

  changed <- fit
  changed$diagnostics$reference$per_pv_draws$PV1[1L, 1L] <-
    changed$diagnostics$reference$per_pv_draws$PV1[1L, 1L] + 1
  expect_error(
    pvstackr:::validate_pvstackr_fit(changed),
    "reproduce their pooled mean/covariance"
  )

  reordered <- fit
  reordered$diagnostics$reference$per_pv_draws <- lapply(
    reordered$diagnostics$reference$per_pv_draws,
    function(x) {
      out <- x[nrow(x):1L, , drop = FALSE]
      rownames(out) <- NULL
      out
    }
  )
  expect_error(
    pvstackr:::validate_pvstackr_fit(reordered),
    "validation stamp"
  )

  injected_summary <- redacted
  injected_summary$diagnostics$pooling$U_bar <- matrix(1, 2L, 2L)
  expect_error(
    pvstackr:::validate_pvstackr_fit(injected_summary),
    "exactly mirror"
  )

  injected_summary <- redacted
  injected_summary$diagnostics$reference$map_sources <- matrix(1, 2L, 2L)
  expect_error(
    pvstackr:::validate_pvstackr_fit(injected_summary),
    "exactly mirror"
  )
})

test_that("per_pv reportable payload cannot be relabeled as a blocked fit", {
  fit <- pv_fit_reference(
    per_pv_draws = reference_draws_fixture(),
    control = reference_control(return_draws = TRUE)
  )
  relabeled <- fit
  relabeled$status <- "blocked"
  relabeled$reason_codes <- "backend_failed"
  relabeled <- pvstackr:::pv_fit_issue_validation_stamp(relabeled)
  expect_error(
    pvstackr:::validate_pvstackr_fit(relabeled),
    "canonical empty estimates"
  )

  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "per_pv",
      status = "blocked",
      reason_codes = "backend_failed",
      control = pvstackr:::pv_fit_blocked_control(reference_control())
    ),
    "no typed blocked-object path"
  )
})

test_that("per_pv keep_data uses a data-free snapshot and canonical formula", {
  data <- reference_fixture_data()
  data$UNUSED_REFERENCE_PRIVATE_COLUMN <- rep(
    "REFERENCE_RAW_VALUE_SENTINEL",
    nrow(data)
  )
  formula_env <- new.env(parent = baseenv())
  formula_env$private_state <- "REFERENCE_FORMULA_ENV_SENTINEL"
  formula <- OUTCOME ~ x
  environment(formula) <- formula_env
  draws <- reference_draws_fixture()

  redacted <- pv_fit_reference(
    data = data,
    formula = formula,
    pv_cols = names(draws),
    per_pv_draws = draws,
    control = reference_control(keep_data = FALSE),
    weight_col = "W",
    id_cols = "id"
  )
  retained <- pv_fit_reference(
    data = data,
    formula = formula,
    pv_cols = names(draws),
    per_pv_draws = draws,
    control = reference_control(keep_data = TRUE),
    weight_col = "W",
    id_cols = "id"
  )
  expect_null(redacted$design$data)
  expect_identical(environment(redacted$design$formula), baseenv())
  expect_identical(retained$design$data, data)
  expect_identical(environment(retained$design$formula), baseenv())
  serialized <- rawToChar(serialize(redacted, NULL, ascii = TRUE))
  for (marker in c(
    "UNUSED_REFERENCE_PRIVATE_COLUMN", "REFERENCE_RAW_VALUE_SENTINEL",
    "REFERENCE_FORMULA_ENV_SENTINEL"
  )) {
    expect_false(grepl(marker, serialized, fixed = TRUE), info = marker)
  }
  expect_lt(
    length(serialize(redacted, NULL, xdr = TRUE)),
    length(serialize(retained, NULL, xdr = TRUE))
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(redacted))

  self_rehashed <- redacted
  self_rehashed$design$data_manifest$column_classes[["x"]] <- "integer"
  self_rehashed$design$design_hash <- pvstackr:::pv_hash_payload(c(
    self_rehashed$design$data_manifest,
    list(
      row_support_hash = self_rehashed$design$row_support_hash,
      pv_value_hash = self_rehashed$design$pv_value_hash,
      weight_design_hash = self_rehashed$design$weight_design_hash
    )
  ))
  expect_invisible(pvstackr:::validate_pvstackr_design(self_rehashed$design))
  expect_error(
    pvstackr:::validate_pvstackr_fit(self_rehashed),
    "validation stamp"
  )

  formula_environment_injection <- redacted
  private_environment <- new.env(parent = baseenv())
  private_environment$private_payload <- "PRIVATE_FORMULA_ENVIRONMENT_PAYLOAD"
  environment(formula_environment_injection$design$formula) <- private_environment
  expect_error(
    pvstackr:::validate_pvstackr_fit(
      formula_environment_injection,
      tier = "cheap"
    ),
    "formula environment"
  )

  control_flip <- redacted
  control_flip$control$keep_data <- TRUE
  expect_error(pvstackr:::validate_pvstackr_fit(control_flip))
  backend_injection <- redacted
  backend_injection$diagnostics$reference["backend_fits"] <- list(
    list(list(data = data))
  )
  expect_error(pvstackr:::validate_pvstackr_fit(backend_injection))
})

test_that("pv_fit_reference does not expose an inert allow_m1 API", {
  expect_false("allow_m1" %in% names(formals(pv_fit_reference)))

  draws <- reference_draws_fixture()
  expect_error(
    pv_fit_reference(
      per_pv_draws = draws[1],
      control = reference_control()
    ),
    "at least two"
  )
})

test_that("pv_fit_reference normalizes Barnard-Rubin df_complete before pooling", {
  draws <- reference_draws_fixture()
  named_df <- c(b_x = 60, b_Intercept = 40)
  oracle <- reference_pool_oracle(
    draws,
    df_method = "barnard_rubin",
    df_complete = c(b_Intercept = 40, b_x = 60)
  )

  fit <- pv_fit_reference(
    per_pv_draws = draws,
    control = reference_control(),
    df_method = "barnard_rubin",
    df_complete = named_df
  )

  expect_equal(fit$target$df_method, "barnard_rubin")
  expect_equal(fit$target$df_complete, c(b_Intercept = 40, b_x = 60), tolerance = 0)
  expect_equal(fit$target$df, oracle$df, tolerance = 1e-14)
  expect_equal(fit$estimates$df_complete, c(40, 60), tolerance = 0)
  expect_equal(fit$estimates$interval_role, rep("reference_barnard_rubin", 2L))
  expect_false(any(fit$estimates$coverage_claim_allowed))
  expect_true(all(fit$target$df <= fit$target$df_classic))
  expect_equal(fit$diagnostics$pooling$df_complete, fit$target$df_complete)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_reference accepts scalar and named df_complete and hashes normalized policy", {
  draws <- reference_draws_fixture()

  scalar_fit <- pv_fit_reference(
    per_pv_draws = draws,
    control = reference_control(),
    df_method = "barnard_rubin",
    df_complete = 40
  )
  named_fit <- pv_fit_reference(
    per_pv_draws = draws,
    control = reference_control(),
    df_method = "barnard_rubin",
    df_complete = c(b_x = 40, b_Intercept = 40)
  )
  canonical_fit <- pv_fit_reference(
    per_pv_draws = draws,
    control = reference_control(),
    df_method = "barnard_rubin",
    df_complete = c(b_Intercept = 40, b_x = 40)
  )
  reversed_fit <- pv_fit_reference(
    per_pv_draws = draws,
    control = reference_control(),
    df_method = "barnard_rubin",
    df_complete = c(b_x = 60, b_Intercept = 40)
  )
  reordered_canonical_fit <- pv_fit_reference(
    per_pv_draws = draws,
    control = reference_control(),
    df_method = "barnard_rubin",
    df_complete = c(b_Intercept = 40, b_x = 60)
  )

  expect_equal(scalar_fit$target$df_complete, c(b_Intercept = 40, b_x = 40), tolerance = 0)
  expect_equal(named_fit$target$df_complete, scalar_fit$target$df_complete, tolerance = 0)
  expect_equal(named_fit$target$target_hash, scalar_fit$target$target_hash)
  expect_equal(named_fit$diagnostics$pooling$target_hash, scalar_fit$diagnostics$pooling$target_hash)
  expect_equal(canonical_fit$target$target_hash, scalar_fit$target$target_hash)
  expect_equal(reversed_fit$target$df_complete, c(b_Intercept = 40, b_x = 60), tolerance = 0)
  expect_equal(reversed_fit$target$target_hash, reordered_canonical_fit$target$target_hash)
  expect_equal(reversed_fit$target$df, reordered_canonical_fit$target$df, tolerance = 0)
})

test_that("pv_fit_reference reports actionable df_complete errors", {
  draws <- reference_draws_fixture()

  expect_error(
    pv_fit_reference(
      per_pv_draws = draws,
      control = reference_control(),
      df_method = "classic",
      df_complete = 40
    ),
    'only supported when `df_method = "barnard_rubin"`'
  )
  expect_error(
    pv_fit_reference(
      per_pv_draws = draws,
      control = reference_control(),
      df_method = "barnard_rubin"
    ),
    'required when `df_method = "barnard_rubin"`'
  )
  expect_error(
    pv_fit_reference(
      per_pv_draws = draws,
      control = reference_control(),
      df_method = "barnard_rubin",
      df_complete = c(40, 60)
    ),
    "scalar or a named numeric vector aligned to `fe_names`"
  )
  expect_error(
    pv_fit_reference(
      per_pv_draws = draws,
      control = reference_control(),
      df_method = "barnard_rubin",
      df_complete = c(b_Intercept = 40, missing = 60)
    ),
    "scalar or a named numeric vector aligned to `fe_names`"
  )
  expect_error(
    pv_fit_reference(
      per_pv_draws = draws,
      control = reference_control(),
      df_method = "barnard_rubin",
      df_complete = c(b_Intercept = 40, b_Intercept = 60)
    ),
    "scalar or a named numeric vector aligned to `fe_names`"
  )
  expect_error(
    pv_fit_reference(
      per_pv_draws = draws,
      control = reference_control(),
      df_method = "barnard_rubin",
      df_complete = c(b_Intercept = 40, b_x = 0)
    ),
    "`df_complete` must be positive"
  )
  expect_error(
    pv_fit_reference(
      per_pv_draws = draws,
      control = reference_control(),
      df_method = "barnard_rubin",
      df_complete = c(b_Intercept = 40, b_x = NA_real_)
    ),
    "`df_complete` must be positive"
  )
})

test_that("pv_fit_reference runs one injected fit per plausible value in PV order", {
  data <- reference_fixture_data()
  draws <- reference_draws_fixture()
  record <- new.env(parent = emptyenv())
  record$formulas <- character()
  record$n_rows <- integer()
  record$files <- character()
  record$extra_names <- list()

  fake_fit <- function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    lhs <- all.vars(formula)[[1L]]
    extra <- list(...)
    record$formulas <- c(record$formulas, paste(deparse(formula, width.cutoff = 500L), collapse = ""))
    record$n_rows <- c(record$n_rows, nrow(data))
    record$files <- c(record$files, file)
    record$extra_names[[length(record$extra_names) + 1L]] <- names(extra)
    list(draws = draws[[lhs]], lhs = lhs)
  }

  fit <- pv_fit_reference(
    data = data,
    formula = OUTCOME ~ x + (1 | school),
    pv_cols = c("PV1", "PV2", "PV3"),
    control = reference_control(keep_backend_fit = TRUE, keep_data = TRUE),
    fit_function = fake_fit,
    draws_function = function(fit) fit$draws,
    weight_col = "W",
    cache_dir = tempdir(),
    cache_stem = "reference-fit",
    additional_args = list(adapter_weight_col = "W")
  )

  expect_equal(record$formulas, c(
    "PV1 ~ x + (1 | school)",
    "PV2 ~ x + (1 | school)",
    "PV3 ~ x + (1 | school)"
  ))
  expect_equal(record$n_rows, rep(nrow(data), 3L))
  expect_equal(length(unique(record$files)), 3L)
  expect_s3_class(fit$design, "pvstackr_design")
  expect_equal(fit$design$roles$method, "per_pv")
  expect_equal(fit$design$weight_col, "W")
  expect_true(all(vapply(record$extra_names, function(x) {
    "adapter_weight_col" %in% x &&
      !"weight_col" %in% x &&
      !"weights" %in% x
  }, logical(1))))
  expect_equal(names(fit$diagnostics$reference$backend_fits), c("PV1", "PV2", "PV3"))
  expect_false(fit$validation$fast_path_eligible)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit, tier = "cheap"))
  opaque_tamper <- fit
  opaque_tamper$estimates$estimate[[1L]] <-
    opaque_tamper$estimates$estimate[[1L]] + 1
  opaque_tamper <- pvstackr:::pv_fit_issue_validation_stamp(opaque_tamper)
  expect_error(
    pvstackr:::validate_pvstackr_fit(opaque_tamper, tier = "cheap"),
    "must match Rubin pooled centers"
  )
})

test_that("pv_fit dispatches per_pv through pv_fit_reference", {
  data <- reference_fixture_data()
  draws <- reference_draws_fixture()
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fake_fit <- function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    lhs <- all.vars(formula)[[1L]]
    record$n <- record$n + 1L
    list(draws = draws[[lhs]])
  }

  fit <- pv_fit(
    data = data,
    formula = OUTCOME ~ x,
    method = "per_pv",
    pv_cols = c("PV1", "PV2", "PV3"),
    control = reference_control(),
    fit_function = fake_fit,
    draws_function = function(fit) fit$draws
  )

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$method, "per_pv")
  expect_equal(record$n, 3L)
  expect_error(
    pv_fit(
      data = data,
      formula = OUTCOME ~ x,
      method = "per_pv",
      pv_cols = c("PV1", "PV2", "PV3"),
      control = pv_control(method = "stack_direct"),
      fit_function = fake_fit,
      draws_function = function(fit) fit$draws
    ),
    "control\\$method"
  )
})

test_that("pv_fit_reference rejects malformed or misaligned per-PV draws", {
  draws <- reference_draws_fixture()

  bad <- draws
  bad$PV2 <- bad$PV2[, c("b_x", "b_Intercept", "sigma", "lp__")]
  expect_error(
    pv_fit_reference(per_pv_draws = bad, control = reference_control()),
    "fixed-effect draw columns"
  )

  bad <- draws
  colnames(bad$PV2)[1] <- ""
  expect_error(
    pv_fit_reference(per_pv_draws = bad, control = reference_control()),
    "column names"
  )

  bad <- draws
  bad$PV2[1, 1] <- Inf
  expect_error(
    pv_fit_reference(per_pv_draws = bad, control = reference_control()),
    "finite"
  )

  expect_error(
    pv_fit_reference(per_pv_draws = draws[1], control = reference_control()),
    "at least two"
  )
  expect_error(
    pv_fit_reference(control = reference_control()),
    "exactly one"
  )
  expect_error(
    pv_fit_reference(
      per_pv_draws = draws,
      fit_function = function(...) list(),
      draws_function = function(fit) fit$draws,
      control = reference_control()
    ),
    "exactly one"
  )
  expect_error(
    pv_fit_reference(per_pv_draws = draws, control = pv_control(method = "stack_direct")),
    "control\\$method"
  )
})

test_that("pvstackr_fit enforces per_pv reportable status invariants", {
  fit <- pv_fit_reference(
    per_pv_draws = reference_draws_fixture(),
    control = reference_control()
  )

  expect_error(
    pvstackr:::new_pvstackr_fit("per_pv", control = reference_control()),
    "pvstackr_reference_pool"
  )

  bad <- fit
  bad$estimates <- data.frame()
  expect_error(
    pvstackr:::validate_pvstackr_fit(bad),
    "non-empty estimates"
  )

  bad <- fit
  bad$estimates$target_hash <- "00000000"
  expect_error(
    pvstackr:::validate_pvstackr_fit(bad),
    "pooling metadata"
  )

  bad <- fit
  bad$estimates$estimate[1] <- bad$estimates$estimate[1] + 1
  expect_error(
    pvstackr:::validate_pvstackr_fit(bad),
    "Rubin pooled centers"
  )

  bad <- fit
  bad$draws <- cbind(b_Intercept = c(1, 2, 3), b_x = c(1, 2, 3))
  expect_error(
    pvstackr:::validate_pvstackr_fit(bad),
    "top-level"
  )
})
