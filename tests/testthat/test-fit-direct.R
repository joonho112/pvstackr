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

fit_direct_legacy_target <- function(target) {
  legacy_fields <- setdiff(
    pvstackr:::pv_brr_target_v02_fields(),
    c("conf_level", "binding_manifest", "target_content")
  )
  legacy <- target[legacy_fields]
  if (identical(legacy$df_method, "classic")) {
    legacy$df_complete <- stats::setNames(
      rep(NA_real_, length(legacy$fe_names)),
      legacy$fe_names
    )
  }
  legacy$design_hash <- target$binding_manifest$legacy_hashes$design_hash
  legacy$target_hash <- target$binding_manifest$legacy_hashes$target_hash
  legacy$schema_version <- "0.1.0"
  legacy$provenance <- list(
    function_name = "pv_brr_target",
    assembled_at = "2026-07-12 12:00:00",
    package = "pvstackr"
  )
  class(legacy) <- c("pvstackr_brr_target", "list")
  pvstackr:::validate_pvstackr_brr_target(legacy)
  legacy
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
    diagnose_function = test_sampler_diagnose_function(
      extra = list(n_long = nrow(data) * 2L)
    ),
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
  expect_identical(fit$ccc$schema_version, "0.2.0")
  expect_identical(
    fit$ccc$binding_proof,
    fit$diagnostics$preflight$binding_proof
  )
  expect_equal(fit$stack_fit$param_map$fe_names, target$fe_names)
  expect_equal(fit$stack_fit$param_map$dropped_names, "lp__")
  expect_true(any(grepl("lp__", fit$stack_fit$warnings)))
  expect_equal(colnames(fit$draws), target$fe_names)
  expect_null(fit$stack_fit$stacked_draws)
  expect_null(fit$ccc$draws_calibrated)
  expect_null(fit$ccc$draws_fe_cal)
  expect_false(fit$ccc$provenance$draws_retained)
  expect_equal(colMeans(fit$draws), fit$ccc$psi_hat, tolerance = 1e-12)
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
  expect_identical(
    fit$diagnostics$preflight$binding_proof$target_manifest_hash,
    target$binding_manifest$manifest_hash
  )
  expect_false("n_long" %in% names(fit$diagnostics$stack_fit))
  expect_equal(fit$design$roles$method, "stack_direct")
  expect_equal(fit$design$provenance$target_hash, target$target_hash)
  expect_identical(
    fit$design$provenance$target_manifest_hash,
    target$binding_manifest$manifest_hash
  )
  expect_identical(
    fit$design$provenance$target_content_hash,
    target$target_content$target_content_hash
  )
  expect_equal(fit$provenance$wrapper_function, "pv_fit_direct")
  expect_equal(fit$provenance$target_hash, target$target_hash)
  expect_identical(
    fit$provenance$target_manifest_hash,
    target$binding_manifest$manifest_hash
  )
  expect_identical(
    fit$provenance$target_content_hash,
    target$target_content$target_content_hash
  )
  expect_identical(
    fit$provenance$binding_verification_policy,
    pvstackr:::pv_binding_proof_policy()
  )
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
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
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
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )

  expect_equal(fit$status, "ok")
  expect_null(fit$draws)
  expect_null(fit$stack_fit$stacked_draws)
  expect_false(fit$stack_fit$control$return_draws)
  expect_null(fit$ccc$draws_calibrated)
  expect_null(fit$ccc$draws_fe_cal)
  expect_false(fit$ccc$provenance$draws_retained)
  expect_equal(fit$estimates$term, target$fe_names)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  injected_top <- fit
  injected_top$draws <- matrix(
    0,
    nrow = 2L,
    ncol = length(target$fe_names),
    dimnames = list(NULL, target$fe_names)
  )
  expect_error(pvstackr:::validate_pvstackr_fit(injected_top), "return_draws")

  injected_stack <- fit
  injected_stack$stack_fit$stacked_draws <- injected_top$draws
  expect_error(pvstackr:::validate_pvstackr_fit(injected_stack), "draw retention")

  injected_ccc <- fit
  injected_ccc$ccc$draws_fe_cal <- injected_top$draws
  expect_error(pvstackr:::validate_pvstackr_fit(injected_ccc), "retention flag")

  injected_diagnostics <- fit
  injected_diagnostics$diagnostics$ccc$private_draws <- injected_top$draws
  expect_error(
    pvstackr:::validate_pvstackr_fit(injected_diagnostics),
    "nested diagnostics must exactly mirror"
  )

  injected_diagnostics <- fit
  injected_diagnostics$diagnostics$stack_fit$private_draws <-
    injected_top$draws
  expect_error(
    pvstackr:::validate_pvstackr_fit(injected_diagnostics),
    "nested diagnostics must exactly mirror"
  )
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
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) as_draws_matrix_like(fit$draws),
    diagnose_function = test_sampler_diagnose_function()
  )

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$status, "ok")
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
  expect_identical(class(fit$draws), c("matrix", "array"))
  expect_null(fit$ccc$draws_calibrated)
  expect_null(fit$ccc$draws_fe_cal)
  expect_equal(fit$draws, plain_fit$draws, tolerance = 0)
  expect_equal(fit$estimates, plain_fit$estimates, tolerance = 0)
})

test_that("pv_fit_direct preflight rejects incompatible calls before fitting", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  mismatch <- tryCatch(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x + z,
      target = target,
      control = fit_direct_control(),
      fit_function = fake_fit_direct_fit(record, target),
      draws_function = function(fit) fit$draws,
      diagnose_function = test_sampler_diagnose_function()
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(mismatch$code, "PV_BIND_E030")
  expect_equal(record$n, 0L)

  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x + (1 | school),
      target = target,
      control = fit_direct_control(),
      fit_function = fake_fit_direct_fit(record, target),
      draws_function = function(fit) fit$draws,
      diagnose_function = test_sampler_diagnose_function()
    ),
    "Random-effect/group terms"
  )
  expect_equal(record$n, 0L)
})

test_that("pv_fit_direct binding preflight ignores unused columns and runs backend once", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  changed <- data
  changed$unused_private <- paste0("unused_", rev(seq_len(nrow(changed))))
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = changed,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )
  expect_identical(record$n, 1L)
  expect_s3_class(fit, "pvstackr_fit")
  expect_identical(
    fit$diagnostics$preflight$binding_proof$target_manifest_hash,
    target$binding_manifest$manifest_hash
  )
})

test_that("pv_fit_direct binding mismatches fail before every backend callback", {
  data <- fit_direct_fixture_data()
  formula <- OUTCOME ~ x
  target <- fit_direct_target(formula, data = data)

  reordered <- data[nrow(data):1L, , drop = FALSE]
  changed_row <- data
  changed_row$id[[1L]] <- 1001L
  changed_pv <- data
  changed_pv$PV1[[1L]] <- changed_pv$PV1[[1L]] + 1
  changed_predictor <- data
  changed_predictor$x[[1L]] <- changed_predictor$x[[1L]] + 0.25
  changed_weight <- data
  changed_weight$W[[1L]] <- changed_weight$W[[1L]] + 0.25
  tampered_manifest <- target
  tampered_manifest$binding_manifest$components$pv$M <- 3L
  tampered_content <- target
  tampered_content$target_content$derived$beta[[1L]] <-
    tampered_content$target_content$derived$beta[[1L]] + 0.1
  relabelled_schema <- target
  relabelled_schema$schema_version <- "0.1.0"
  unbound <- target
  unbound$binding_manifest <- NULL
  unbound$target_content <- NULL
  legacy <- fit_direct_legacy_target(target)
  rehash_target_component <- function(target, component_name) {
    hash_name <- paste0(component_name, "_hash")
    manifest <- target$binding_manifest
    manifest$component_hashes[[hash_name]] <-
      pvstackr:::pv_binding_component_hash(
        manifest$components[[component_name]],
        component_name
      )
    manifest$manifest_hash <- pvstackr:::pv_binding_hash_payload(
      pvstackr:::pv_binding_manifest_hash_payload(manifest),
      "manifest"
    )
    target$binding_manifest <- manifest
    target$design_hash <- manifest$manifest_hash
    target$target_content$manifest_hash <- manifest$manifest_hash
    target$target_content$target_content_hash <-
      pvstackr:::pv_binding_hash_payload(
        pvstackr:::pv_binding_target_content_hash_payload(
          target$target_content
        ),
        "target_content"
      )
    target$target_hash <- target$target_content$target_content_hash
    target
  }
  incompatible_support <- target
  incompatible_support$binding_manifest$components$family_link$response_support_id <-
    "positive"
  incompatible_support <- rehash_target_component(
    incompatible_support,
    "family_link"
  )
  incompatible_scope <- target
  incompatible_scope$binding_manifest$components$estimand$parameter_scope <-
    "variance_component"
  incompatible_scope <- rehash_target_component(
    incompatible_scope,
    "estimand"
  )

  cases <- list(
    row_order = list(data = reordered, formula = formula, target = target, family = NULL, code = "PV_BIND_E012"),
    row_membership = list(data = changed_row, formula = formula, target = target, family = NULL, code = "PV_BIND_E011"),
    pv_value = list(data = changed_pv, formula = formula, target = target, family = NULL, code = "PV_BIND_E021"),
    predictor_value = list(data = changed_predictor, formula = formula, target = target, family = NULL, code = "PV_BIND_E031"),
    weight_value = list(data = changed_weight, formula = formula, target = target, family = NULL, code = "PV_BIND_E051"),
    formula = list(data = data, formula = OUTCOME ~ I(x), target = target, family = NULL, code = "PV_BIND_E040"),
    family = list(data = data, formula = formula, target = target, family = stats::binomial(), code = "PV_BIND_E060"),
    link = list(data = data, formula = formula, target = target, family = stats::gaussian(link = "log"), code = "PV_BIND_E061"),
    target_support = list(data = data, formula = formula, target = incompatible_support, family = NULL, code = "PV_BIND_E060"),
    target_scope = list(data = data, formula = formula, target = incompatible_scope, family = NULL, code = "PV_BIND_E070"),
    target_manifest = list(data = data, formula = formula, target = tampered_manifest, family = NULL, code = "PV_BIND_E005"),
    target_content = list(data = data, formula = formula, target = tampered_content, family = NULL, code = "PV_BIND_E090"),
    schema = list(data = data, formula = formula, target = relabelled_schema, family = NULL, code = "PV_BIND_E080"),
    unbound = list(data = data, formula = formula, target = unbound, family = NULL, code = "PV_BIND_E080"),
    legacy = list(data = data, formula = formula, target = legacy, family = NULL, code = "PV_BIND_E080")
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    callbacks <- new.env(parent = emptyenv())
    callbacks$fit <- 0L
    callbacks$draws <- 0L
    callbacks$diagnose <- 0L
    callbacks$loglik <- 0L
    error <- tryCatch(
      pv_fit_direct(
        data = case$data,
        formula = case$formula,
        target = case$target,
        family = case$family,
        control = fit_direct_control(keep_log_lik = TRUE),
        fit_function = function(...) { callbacks$fit <- callbacks$fit + 1L; list() },
        draws_function = function(fit) { callbacks$draws <- callbacks$draws + 1L; matrix(0) },
        diagnose_function = function(fit) { callbacks$diagnose <- callbacks$diagnose + 1L; list() },
        log_lik_function = function(fit) { callbacks$loglik <- callbacks$loglik + 1L; matrix(0) },
        extract_log_lik = TRUE
      ),
      pvstackr_binding_error = identity
    )
    expect_s3_class(error, "pvstackr_binding_error")
    expect_identical(error$code, case$code, info = name)
    expect_identical(
      c(callbacks$fit, callbacks$draws, callbacks$diagnose, callbacks$loglik),
      rep(0L, 4L),
      info = name
    )
  }
})

test_that("pv_fit_direct materializes verified stateful RHS without backend evaluation", {
  data <- fit_direct_fixture_data()
  state <- new.env(parent = baseenv())
  state$count <- 0L
  state$pulse_on_third <- FALSE
  state$transform_x <- function(x) {
    state$count <- state$count + 1L
    if (isTRUE(state$pulse_on_third) && state$count == 3L) 2 * x else x
  }
  formula <- stats::as.formula("OUTCOME ~ transform_x(x)", env = state)
  target <- fit_direct_target(formula, data = data)
  expected_bundle <- pvstackr:::pv_binding_resolve_model_bundle(data, formula)

  record <- new.env(parent = emptyenv())
  record$n <- 0L
  state$count <- 0L
  state$pulse_on_third <- TRUE
  fit_function <- function(formula, data, family, prior, chains, iter, warmup,
                           cores, seed, backend, file, file_refit, ...) {
    record$n <- record$n + 1L
    record$formula <- formula
    record$data <- data
    record$count_before_matrix <- state$count
    rhs_terms <- stats::delete.response(stats::terms(formula))
    record$backend_matrix <- stats::model.matrix(rhs_terms, data = data)
    record$count_after_matrix <- state$count

    beta <- target$beta[target$fe_names]
    se <- sqrt(diag(target$T_MI))[target$fe_names]
    e1 <- c(-0.03, -0.01, 0.02, 0.04, -0.02, 0.01, 0.03, -0.04)
    e2 <- c(0.02, -0.03, 0.01, -0.04, 0.03, -0.01, 0.04, -0.02)
    draws <- cbind(
      beta[[1L]] + se[[1L]] * e1,
      beta[[2L]] + se[[2L]] * e2,
      sigma = rep(1, 8),
      lp__ = seq(-10, -9, length.out = 8)
    )
    colnames(draws)[1:2] <- c("b_pvstackrMM001", "b_pvstackrMM002")
    record$returned_draws <- draws
    list(draws = draws)
  }

  fit <- pv_fit_direct(
    data = data,
    formula = formula,
    target = target,
    control = fit_direct_control(),
    fit_function = fit_function,
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )

  internal_cols <- c("pvstackrMM001", "pvstackrMM002")
  expected_long_x <- expected_bundle$model_matrix[
    rep(seq_len(nrow(data)), times = length(target$pv_cols)),
    ,
    drop = FALSE
  ]
  expect_s3_class(fit, "pvstackr_fit")
  expect_identical(fit$status, "ok")
  expect_identical(record$n, 1L)
  expect_identical(record$count_before_matrix, 1L)
  expect_identical(record$count_after_matrix, 1L)
  expect_identical(state$count, 1L)
  expect_identical(environment(record$formula), asNamespace("stats"))
  expect_false(grepl("transform_x", paste(deparse(record$formula), collapse = ""), fixed = TRUE))
  expect_true(all(internal_cols %in% names(record$data)))
  expect_equal(
    unname(as.matrix(record$data[internal_cols])),
    unname(expected_long_x),
    tolerance = 0
  )
  expect_equal(
    unname(record$backend_matrix),
    unname(as.matrix(record$data[internal_cols])),
    tolerance = 0,
    ignore_attr = TRUE
  )
  expect_identical(
    dim(record$backend_matrix),
    dim(as.matrix(record$data[internal_cols]))
  )
  expect_identical(
    colnames(record$returned_draws)[1:2],
    c("b_pvstackrMM001", "b_pvstackrMM002")
  )
  expect_identical(
    fit$stack_fit$param_map$fe_names,
    target$fe_names
  )
  expect_false(any(grepl("^b_pvstackrMM", colnames(fit$ccc$draws_calibrated))))
  expect_identical(
    fit$ccc$binding_proof,
    fit$diagnostics$preflight$binding_proof
  )
  expect_invisible(pvstackr:::pv_binding_proof_validate(
    fit$ccc$binding_proof,
    target_manifest = target$binding_manifest
  ))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
  expect_s3_class(get_estimates(fit), "data.frame")
  expect_identical(state$count, 1L)
})

test_that("pv_fit_direct materializes the verified offset with its model matrix", {
  data <- fit_direct_fixture_data()
  formula <- OUTCOME ~ x + offset(z)
  target <- fit_direct_target(formula, data = data)
  bundle <- pvstackr:::pv_binding_resolve_model_bundle(data, formula)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = data,
    formula = formula,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )

  row_index <- rep(seq_len(nrow(data)), times = length(target$pv_cols))
  internal_cols <- sprintf("pvstackrMM%03d", seq_len(ncol(bundle$model_matrix)))
  backend_terms <- stats::delete.response(stats::terms(record$formula))
  backend_frame <- stats::model.frame(
    backend_terms,
    data = record$data,
    na.action = stats::na.fail
  )
  backend_matrix <- stats::model.matrix(backend_terms, data = backend_frame)
  backend_offset <- stats::model.offset(backend_frame)
  expect_identical(record$n, 1L)
  expect_identical(environment(record$formula), asNamespace("stats"))
  expect_match(
    paste(deparse(record$formula), collapse = ""),
    "offset\\(pvstackrOffset\\)"
  )
  expect_equal(
    unname(as.matrix(record$data[internal_cols])),
    unname(bundle$model_matrix[row_index, , drop = FALSE]),
    tolerance = 0
  )
  expect_identical(
    unname(record$data$pvstackrOffset),
    unname(bundle$offset[row_index])
  )
  expected_backend_matrix <- unname(as.matrix(record$data[internal_cols]))
  expect_identical(dim(backend_matrix), dim(expected_backend_matrix))
  expect_equal(
    as.vector(backend_matrix),
    as.vector(expected_backend_matrix),
    tolerance = 0
  )
  expect_identical(
    unname(backend_offset),
    unname(bundle$offset[row_index])
  )
  expect_identical(fit$stack_fit$meta$model_matrix_materialized, TRUE)
  expect_identical(
    fit$stack_fit$meta$model_matrix_bundle_hash,
    bundle$bundle_hash
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_direct replaces forged Gaussian closures with canonical family", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  forged_family <- stats::gaussian(link = "identity")
  forged_family$linkfun <- function(mu) rep(-999, length(mu))
  forged_family$linkinv <- function(eta) rep(-999, length(eta))
  forged_family$variance <- function(mu) rep(-999, length(mu))
  forged_family$secret_payload <- "CALLER_PRIVATE_FAMILY_SECRET"
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    family = forged_family,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )

  backend_family <- record$args$family
  expect_identical(record$n, 1L)
  expect_s3_class(backend_family, "family")
  expect_identical(backend_family$family, "gaussian")
  expect_identical(backend_family$link, "identity")
  expect_identical(backend_family$linkfun(2), 2)
  expect_false(identical(backend_family$linkfun, forged_family$linkfun))
  serialized_family <- rawToChar(
    serialize(backend_family, NULL, ascii = TRUE)
  )
  expect_false(grepl(
    forged_family$secret_payload,
    serialized_family,
    fixed = TRUE
  ))
  expect_false(grepl("secret_payload", serialized_family, fixed = TRUE))
  expect_identical(fit$status, "ok")
  expect_identical(
    fit$ccc$binding_proof,
    fit$diagnostics$preflight$binding_proof
  )
  expect_invisible(pvstackr:::pv_binding_proof_validate(
    fit$ccc$binding_proof,
    target_manifest = target$binding_manifest
  ))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("keep_data enforces a recursive data-free direct-fit snapshot", {
  data <- fit_direct_fixture_data()
  data$UNUSED_PRIVATE_COLUMN_NAME <- rep(
    "RAW_ANALYSIS_VALUE_SENTINEL",
    nrow(data)
  )
  formula_env <- new.env(parent = baseenv())
  formula_env$private_formula_state <- "FORMULA_ENV_SENTINEL"
  formula <- OUTCOME ~ x
  environment(formula) <- formula_env
  target <- fit_direct_target(formula, data = data)

  build <- function(keep_data) {
    record <- new.env(parent = emptyenv())
    record$n <- 0L
    pv_fit_direct(
      data = data,
      formula = formula,
      target = target,
      control = fit_direct_control(
        keep_data = keep_data,
        keep_log_lik = TRUE
      ),
      fit_function = fake_fit_direct_fit(record, target),
      draws_function = function(fit) fit$draws,
      diagnose_function = test_sampler_diagnose_function(
        extra = list(
          private_nested_data = data,
          private_diagnostic_secret = "DIAGNOSTIC_SENTINEL"
        )
      ),
      extract_log_lik = TRUE,
      log_lik_function = function(fit) {
        out <- fit$log_lik
        attr(out, "private_marker") <- "LOG_LIK_ATTR_SENTINEL"
        out
      }
    )
  }

  redacted <- build(FALSE)
  retained <- build(TRUE)
  expect_null(redacted$design$data)
  expect_identical(redacted$design$schema_version, "0.2.0")
  expect_identical(environment(redacted$design$formula), baseenv())
  expect_identical(
    environment(redacted$diagnostics$preflight$formula),
    baseenv()
  )
  expect_null(redacted$stack_fit$prepared_data)
  expect_null(redacted$stack_fit$fit)
  expect_identical(names(redacted$stack_fit$diagnostics), "sampler")
  expect_null(attributes(redacted$stack_fit$log_lik)$private_marker)
  serialized_redacted <- rawToChar(
    serialize(redacted, NULL, ascii = TRUE)
  )
  for (marker in c(
    "UNUSED_PRIVATE_COLUMN_NAME", "RAW_ANALYSIS_VALUE_SENTINEL",
    "FORMULA_ENV_SENTINEL", "DIAGNOSTIC_SENTINEL",
    "LOG_LIK_ATTR_SENTINEL", "private_nested_data"
  )) {
    expect_false(grepl(marker, serialized_redacted, fixed = TRUE), info = marker)
  }

  expect_identical(retained$design$data, data)
  expect_identical(environment(retained$design$formula), baseenv())
  expect_identical(
    environment(retained$diagnostics$preflight$formula),
    baseenv()
  )
  expect_null(retained$stack_fit$prepared_data)
  expect_lt(
    length(serialize(redacted, NULL, xdr = TRUE)),
    length(serialize(retained, NULL, xdr = TRUE))
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(redacted))
  expect_invisible(pvstackr:::validate_pvstackr_fit(retained))

  mutations <- list(
    keep_data_flip = function(x) { x$control$keep_data <- TRUE; x },
    snapshot_role = function(x) { x$design$roles$method <- "per_pv"; x },
    snapshot_source = function(x) { x$design$provenance$source <- "forged"; x },
    snapshot_hash = function(x) {
      forged <- paste0("sha256:", strrep("d", 64L))
      x$design$row_support_hash <- forged
      x$design$row_support$hash <- forged
      x$design$data_manifest$row_hash <- forged
      x
    },
    snapshot_formula_environment = function(x) {
      environment(x$design$formula) <- globalenv()
      x
    },
    nested_prepared_data = function(x) {
      x$stack_fit["prepared_data"] <- list(data)
      x
    }
  )
  for (name in names(mutations)) {
    expect_error(
      pvstackr:::validate_pvstackr_fit(mutations[[name]](redacted)),
      info = name
    )
  }
})

test_that("opaque backend retention cannot bypass keep_data before callbacks", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  callbacks <- new.env(parent = emptyenv())
  callbacks$n <- 0L
  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x,
      target = target,
      control = fit_direct_control(
        keep_data = FALSE,
        keep_backend_fit = TRUE
      ),
      fit_function = function(...) {
        callbacks$n <- callbacks$n + 1L
        list()
      },
      draws_function = function(fit) matrix(0)
    ),
    "requires `keep_data = TRUE`"
  )
  expect_identical(callbacks$n, 0L)
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
    prior = data.frame(
      prior = "normal(0, 1)",
      class = "sigma",
      coef = "",
      stringsAsFactors = FALSE
    ),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )

  expect_equal(fit$status, "warning")
  expect_equal(fit$reason_codes, "explicit_prior_warning")
  expect_true(any(grepl("Explicit priors", fit$warnings)))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pv_fit_direct rejects coefficient-specific priors before backend callbacks", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  callbacks <- new.env(parent = emptyenv())
  callbacks$fit <- 0L
  callbacks$draws <- 0L
  callbacks$diagnose <- 0L

  expect_error(
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x,
      target = target,
      control = fit_direct_control(),
      prior = data.frame(
        prior = "normal(0, 1)",
        class = "b",
        coef = "x",
        stringsAsFactors = FALSE
      ),
      fit_function = function(...) {
        callbacks$fit <- callbacks$fit + 1L
        list()
      },
      draws_function = function(fit) {
        callbacks$draws <- callbacks$draws + 1L
        matrix(0)
      },
      diagnose_function = function(fit) {
        callbacks$diagnose <- callbacks$diagnose + 1L
        list()
      }
    ),
    "cannot be preserved exactly"
  )
  expect_identical(
    c(callbacks$fit, callbacks$draws, callbacks$diagnose),
    rep(0L, 3L)
  )
})

test_that("pv_fit_direct promotes sampler warnings and blocks before CCC", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L

  warning_fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(rhat_max = 1.02)
  )
  expect_identical(warning_fit$status, "warning")
  expect_identical(warning_fit$reason_codes, "sampler_rhat_warning")
  expect_true(nrow(warning_fit$estimates) > 0L)
  expect_identical(warning_fit$ccc$schema_version, "0.2.0")
  expect_identical(
    warning_fit$ccc$binding_proof,
    warning_fit$diagnostics$preflight$binding_proof
  )
  expect_identical(
    warning_fit$diagnostics$sampler_gate$status,
    "warning"
  )

  blocked <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(divergences = 1L)
  )
  expect_identical(blocked$status, "blocked")
  expect_identical(blocked$reason_codes, "sampler_divergences_blocked")
  expect_false(blocked$control$return_draws)
  expect_s3_class(blocked$target, "pvstackr_brr_target")
  expect_true(blocked$provenance$independent_target_retained)
  expect_null(blocked$design)
  expect_null(blocked$stack_fit)
  expect_null(blocked$ccc)
  expect_null(blocked$draws)
  expect_equal(nrow(blocked$estimates), 0L)
  expect_setequal(
    names(blocked$diagnostics),
    c("preflight", "sampler", "sampler_gate", "redaction")
  )
  expect_identical(
    blocked$diagnostics$redaction,
    pvstackr:::pv_fit_direct_blocked_redaction("sampler_gate")
  )
  expect_invisible(pvstackr:::validate_pvstackr_fit(blocked))

  hidden_attribute <- blocked
  attr(
    hidden_attribute$diagnostics$sampler$rhat_max,
    "private_draws"
  ) <- matrix(1, 10L, 10L)
  hidden_attribute <- pvstackr:::pv_fit_issue_validation_stamp(hidden_attribute)
  expect_error(
    pvstackr:::validate_pvstackr_fit(hidden_attribute),
    "hidden attributes"
  )

  incomplete <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws
  )
  expect_identical(incomplete$status, "blocked")
  expect_identical(
    incomplete$reason_codes,
    "sampler_diagnostics_incomplete"
  )
  expect_s3_class(incomplete$target, "pvstackr_brr_target")
})

test_that("schema-0.2 sampler snapshot preserves exact links and rejects raw extras", {
  data <- fit_direct_fixture_data()
  private <- new.env(parent = baseenv())
  private$secret_raw <- "PRIVATE_V02_FORMULA_ENVIRONMENT"
  formula <- stats::as.formula("OUTCOME ~ x", env = private)
  target <- fit_direct_target(formula, data = data, conf_level = 0.9)

  expect_identical(environment(target$formula), baseenv())
  snapshot <- pvstackr:::pv_fit_direct_independent_target(target)
  expect_identical(snapshot, target)
  expect_identical(pvstackr:::pv_fit_direct_independent_target(snapshot), snapshot)
  expect_identical(snapshot$conf_level, 0.9)
  expect_identical(snapshot$binding_manifest, target$binding_manifest)
  expect_identical(snapshot$target_content, target$target_content)
  expect_identical(snapshot$design_hash, snapshot$binding_manifest$manifest_hash)
  expect_identical(snapshot$target_hash, snapshot$target_content$target_content_hash)
  expect_false(grepl(
    private$secret_raw,
    rawToChar(serialize(snapshot, NULL, ascii = TRUE)),
    fixed = TRUE
  ))

  record <- new.env(parent = emptyenv())
  record$n <- 0L
  blocked <- pv_fit_direct(
    data = data,
    formula = formula,
    target = target,
    control = fit_direct_control(),
    fit_function = fake_fit_direct_fit(record, target),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(divergences = 1L)
  )
  expect_identical(blocked$target, target)
  expect_identical(record$n, 1L)
  expect_identical(
    blocked$diagnostics$preflight$binding_proof$target_manifest_hash,
    target$binding_manifest$manifest_hash
  )
  expect_identical(
    blocked$provenance$target_manifest_hash,
    target$binding_manifest$manifest_hash
  )
  expect_identical(
    blocked$provenance$target_content_hash,
    target$target_content$target_content_hash
  )
  expect_identical(
    blocked$provenance$binding_verification_policy,
    pvstackr:::pv_binding_proof_policy()
  )
  blocked_mutations <- list(
    proof = function(x) {
      x$diagnostics$preflight$binding_proof$current_manifest_hash <-
        paste0("sha256:", strrep("0", 64L))
      x
    },
    manifest_link = function(x) {
      x$provenance$target_manifest_hash <- paste0("sha256:", strrep("1", 64L))
      x
    },
    content_link = function(x) {
      x$provenance$target_content_hash <- paste0("sha256:", strrep("2", 64L))
      x
    },
    policy_link = function(x) {
      x$provenance$binding_verification_policy <- "unchecked"
      x
    }
  )
  for (mutation in blocked_mutations) {
    expect_error(pvstackr:::validate_pvstackr_fit(mutation(blocked)))
  }

  raw_mutations <- list(
    root = function(x) { x$raw_data <- data; x },
    policy = function(x) { x$policy$raw_data <- data; x },
    provenance = function(x) { x$provenance$raw_data <- data; x },
    per_pv = function(x) { x$per_pv[[1L]]$raw_data <- data; x },
    warnings = function(x) { x$warnings <- private$secret_raw; x }
  )
  for (name in names(raw_mutations)) {
    record$n <- 0L
    candidate <- raw_mutations[[name]](target)
    expect_error(
      pv_fit_direct(
        data = data,
        formula = formula,
        target = candidate,
        control = fit_direct_control(),
        fit_function = fake_fit_direct_fit(record, candidate),
        draws_function = function(fit) fit$draws,
        diagnose_function = test_sampler_diagnose_function(divergences = 1L)
      ),
      info = name
    )
    expect_identical(record$n, 0L, info = name)
  }
})

test_that("schema-0.1 target remains inspectable but public direct fit fails E080", {
  data <- fit_direct_fixture_data()
  formula <- OUTCOME ~ x
  target <- fit_direct_legacy_target(fit_direct_target(formula, data = data))
  expect_invisible(pvstackr:::validate_pvstackr_brr_target(target))
  inspected <- pvstackr:::pv_fit_direct_independent_target(target)
  expect_identical(inspected$schema_version, "0.1.0")

  callbacks <- new.env(parent = emptyenv())
  callbacks$fit <- 0L
  callbacks$draws <- 0L
  callbacks$diagnose <- 0L
  callbacks$loglik <- 0L
  error <- tryCatch(
    pv_fit_direct(
      data = data,
      formula = formula,
      target = target,
      control = fit_direct_control(),
      fit_function = function(...) { callbacks$fit <- callbacks$fit + 1L; list() },
      draws_function = function(fit) { callbacks$draws <- callbacks$draws + 1L; matrix(0) },
      diagnose_function = function(fit) { callbacks$diagnose <- callbacks$diagnose + 1L; list() },
      log_lik_function = function(fit) { callbacks$loglik <- callbacks$loglik + 1L; matrix(0) },
      extract_log_lik = TRUE
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(error$code, "PV_BIND_E080")
  expect_identical(
    unlist(as.list(callbacks), use.names = FALSE),
    rep(0L, 4L)
  )
})

test_that("sampler-blocked target snapshot recomputes pooling and BRR topology", {
  data <- fit_direct_fixture_data()
  formula <- OUTCOME ~ x
  target <- fit_direct_target(formula, data = data)

  record <- new.env(parent = emptyenv())
  blocked_attempt <- function(candidate) {
    record$n <- 0L
    pv_fit_direct(
      data = data,
      formula = formula,
      target = candidate,
      control = fit_direct_control(),
      fit_function = fake_fit_direct_fit(record, candidate),
      draws_function = function(fit) fit$draws,
      diagnose_function = test_sampler_diagnose_function(divergences = 1L)
    )
  }

  valid <- blocked_attempt(target)
  expect_identical(valid$target$beta, target$beta)
  expect_identical(valid$target$U_bar, target$U_bar)
  expect_identical(valid$target$T_MI, target$T_MI)
  expect_identical(valid$target$df, target$df)
  expect_invisible(pvstackr:::validate_pvstackr_fit(valid))

  shifted_beta <- target
  shifted_beta$beta[[1L]] <- shifted_beta$beta[[1L]] + 0.1
  expect_error(blocked_attempt(shifted_beta))
  expect_identical(record$n, 0L)

  inconsistent_u_bar <- target
  inconsistent_u_bar$U_bar[[1L, 1L]] <-
    inconsistent_u_bar$U_bar[[1L, 1L]] + 0.1
  expect_error(blocked_attempt(inconsistent_u_bar))
  expect_identical(record$n, 0L)

  wrong_R <- target
  wrong_R$R <- 3L
  wrong_R$fay_variance_multiplier <-
    1 / (wrong_R$R * (1 - wrong_R$fay_k)^2)
  wrong_R$per_pv <- lapply(wrong_R$per_pv, function(item) {
    item$R <- 3L
    item$fay_variance_multiplier <-
      1 / (item$R * (1 - item$fay_k)^2)
    item
  })
  expect_error(blocked_attempt(wrong_R))
  expect_identical(record$n, 0L)

  wrong_replicate_beta <- target
  wrong_replicate_beta$per_pv[[1L]]$replicate_beta <-
    wrong_replicate_beta$per_pv[[1L]]$replicate_beta[, -1L, drop = FALSE]
  expect_error(blocked_attempt(wrong_replicate_beta))
  expect_identical(record$n, 0L)
})

test_that("sampler warning followed by CCC block retains only slim gate evidence", {
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
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(rhat_max = 1.02)
  )

  expect_identical(fit$status, "blocked")
  expect_identical(
    fit$reason_codes,
    c("sampler_rhat_warning", "center_separation_red")
  )
  expect_identical(fit$diagnostics$sampler_gate$status, "warning")
  expect_null(fit$design)
  expect_null(fit$stack_fit)
  expect_null(fit$ccc)
  expect_s3_class(fit$target, "pvstackr_brr_target")
  expect_true(isTRUE(fit$provenance$independent_target_retained))
  expect_identical(fit$diagnostics$ccc$source, "ccc_reportability_gate")
  expect_identical(fit$diagnostics$ccc$center$status, "blocked")
  expect_identical(fit$diagnostics$ccc$conditioning$status, "ok")
  expect_identical(
    fit$diagnostics$redaction,
    pvstackr:::pv_fit_direct_blocked_redaction("ccc_reportability_gate")
  )
  expect_equal(nrow(fit$estimates), 0L)
  expect_null(fit$draws)
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  stale_preflight_proof <- fit
  stale_preflight_proof$diagnostics$preflight$binding_proof$target_manifest_hash <-
    paste0("sha256:", strrep("0", 64L))
  stale_preflight_proof$diagnostics$preflight$binding_proof$current_manifest_hash <-
    paste0("sha256:", strrep("0", 64L))
  expect_error(
    pvstackr:::validate_pvstackr_fit(stale_preflight_proof),
    class = "pvstackr_binding_error"
  )

  stale_ccc_target_link <- fit
  stale_ccc_target_link$diagnostics$ccc$target_hash <-
    paste0("sha256:", strrep("3", 64L))
  expect_error(
    pvstackr:::validate_pvstackr_fit(stale_ccc_target_link),
    "slim root schema"
  )

  missing_independent_target <- fit
  missing_independent_target["target"] <- list(NULL)
  expect_error(
    pvstackr:::validate_pvstackr_fit(missing_independent_target),
    "canonical independent external target"
  )

  stale_fit_link <- fit
  stale_fit_link$provenance$target_content_hash <-
    paste0("sha256:", strrep("2", 64L))
  expect_error(
    pvstackr:::validate_pvstackr_fit(stale_fit_link),
    "provenance"
  )
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
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
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
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function()
  )

  expect_equal(fit$status, "blocked")
  expect_false(fit$control$return_draws)
  expect_equal(fit$reason_codes, "center_separation_red")
  expect_true(any(grepl("Center separation", fit$warnings)))
  expect_null(fit$design)
  expect_null(fit$stack_fit)
  expect_null(fit$ccc)
  expect_s3_class(fit$target, "pvstackr_brr_target")
  expect_true(fit$provenance$independent_target_retained)
  expect_equal(fit$diagnostics$ccc$center$status, "blocked")
  expect_equal(nrow(fit$estimates), 0L)
  expect_null(fit$draws)
  expect_false(any(vapply(
    fit$control[c("return_draws", "keep_data", "keep_backend_fit", "keep_log_lik")],
    isTRUE,
    logical(1)
  )))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  injected <- fit
  injected$draws <- matrix(
    0,
    nrow = 2L,
    ncol = length(target$fe_names),
    dimnames = list(NULL, target$fe_names)
  )
  expect_error(
    pvstackr:::validate_pvstackr_fit(injected),
    "redact component"
  )
})

test_that("CCC-blocked direct overrides every heavy retention request", {
  data <- fit_direct_fixture_data()
  data$PRIVATE_BLOCKED_DATA <- rep("BLOCKED_RAW_DATA_SENTINEL", nrow(data))
  target <- fit_direct_target(OUTCOME ~ x, data = data)
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  inner_fit <- fake_fit_direct_fit(record, target, center_shift_se = 0.08)
  fit <- pv_fit_direct(
    data = data,
    formula = OUTCOME ~ x,
    target = target,
    control = fit_direct_control(
      return_draws = TRUE,
      keep_data = TRUE,
      keep_backend_fit = TRUE,
      keep_log_lik = TRUE
    ),
    fit_function = function(...) {
      out <- inner_fit(...)
      out$private_data <- "BLOCKED_BACKEND_DATA_SENTINEL"
      out$log_lik <- matrix(987654321, nrow(out$draws), 20L)
      out
    },
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(),
    log_lik_function = function(fit) fit$log_lik,
    extract_log_lik = TRUE
  )

  expect_identical(fit$status, "blocked")
  expect_true(all(vapply(
    fit$control[c("return_draws", "keep_data", "keep_backend_fit", "keep_log_lik")],
    identical,
    logical(1),
    FALSE
  )))
  expect_null(fit$design)
  expect_null(fit$stack_fit)
  expect_null(fit$ccc)
  expect_null(fit$draws)
  expect_s3_class(fit$target, "pvstackr_brr_target")
  serialized <- rawToChar(serialize(fit, NULL, ascii = TRUE))
  expect_false(grepl("BLOCKED_RAW_DATA_SENTINEL", serialized, fixed = TRUE))
  expect_false(grepl("BLOCKED_BACKEND_DATA_SENTINEL", serialized, fixed = TRUE))
  expect_false(grepl("987654321", serialized, fixed = TRUE))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("public conditioning-only and joint CCC blocks use the slim schema", {
  data <- fit_direct_fixture_data()
  target <- fit_direct_target(OUTCOME ~ x, data = data)

  run_fit <- function(center_shift_se) {
    record <- new.env(parent = emptyenv())
    record$n <- 0L
    inner_fit <- fake_fit_direct_fit(
      record,
      target,
      center_shift_se = center_shift_se
    )
    pv_fit_direct(
      data = data,
      formula = OUTCOME ~ x,
      target = target,
      control = fit_direct_control(),
      fit_function = function(...) {
        out <- inner_fit(...)
        e2 <- c(0.02, -0.03, 0.01, -0.04, 0.03, -0.01, 0.04, -0.02)
        shifted_center <- target$beta[["b_x"]] +
          sqrt(target$T_MI["b_x", "b_x"]) * center_shift_se
        out$draws[, "b_x"] <- shifted_center + 1e-12 * e2
        out
      },
      draws_function = function(fit) fit$draws,
      diagnose_function = test_sampler_diagnose_function()
    )
  }

  conditioning_only <- run_fit(0)
  expect_identical(conditioning_only$status, "blocked")
  expect_identical(conditioning_only$reason_codes, "ccc_conditioning_red")
  expect_identical(conditioning_only$diagnostics$ccc$center$status, "ok")
  expect_identical(
    conditioning_only$diagnostics$ccc$conditioning$status,
    "blocked"
  )
  expect_null(conditioning_only$ccc)
  expect_invisible(pvstackr:::validate_pvstackr_fit(conditioning_only))

  joint <- run_fit(0.08)
  expect_identical(joint$status, "blocked")
  expect_identical(
    joint$reason_codes,
    c("center_separation_red", "ccc_conditioning_red")
  )
  expect_identical(joint$diagnostics$ccc$center$status, "blocked")
  expect_identical(joint$diagnostics$ccc$conditioning$status, "blocked")
  expect_null(joint$ccc)
  expect_invisible(pvstackr:::validate_pvstackr_fit(joint))

  hidden_attribute <- conditioning_only
  attr(
    hidden_attribute$diagnostics$ccc$conditioning$kappa_A,
    "private_payload"
  ) <- matrix(1:4, 2L)
  hidden_attribute <- pvstackr:::pv_fit_issue_validation_stamp(hidden_attribute)
  expect_error(
    pvstackr:::validate_pvstackr_fit(hidden_attribute),
    "canonical nonnegative scalar"
  )

  named_reason <- conditioning_only
  names(named_reason$diagnostics$ccc$conditioning$reason_code) <- "private"
  named_reason <- pvstackr:::pv_fit_issue_validation_stamp(named_reason)
  expect_error(
    pvstackr:::validate_pvstackr_fit(named_reason),
    "does not reproduce"
  )

  wrong_independence <- conditioning_only
  wrong_independence$provenance$independent_target_retained <- FALSE
  wrong_independence <- pvstackr:::pv_fit_issue_validation_stamp(wrong_independence)
  expect_error(
    pvstackr:::validate_pvstackr_fit(wrong_independence),
    "independent external target"
  )

  reordered_reasons <- joint
  reordered_reasons$reason_codes <- rev(reordered_reasons$reason_codes)
  reordered_reasons <- pvstackr:::pv_fit_issue_validation_stamp(reordered_reasons)
  expect_error(
    pvstackr:::validate_pvstackr_fit(reordered_reasons),
    "status"
  )
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
      draws_function = function(fit) fit$draws,
      diagnose_function = test_sampler_diagnose_function()
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
      draws_function = function(fit) fit$draws,
      diagnose_function = test_sampler_diagnose_function()
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
      draws_function = function(fit) fit$draws,
      diagnose_function = test_sampler_diagnose_function()
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
    diagnose_function = test_sampler_diagnose_function(),
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
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(
      chains = 4L,
      post_warmup = 1000L
    )
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
