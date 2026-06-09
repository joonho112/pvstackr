classes_fixture_data <- function() {
  data.frame(
    id = seq_len(8),
    school = factor(c("A", "A", "B", "B", "C", "C", "D", "D")),
    x = c(-2, -1, -0.5, 0, 0.5, 1, 1.5, 2),
    PV1 = c(2.1, 2.4, 3.0, 3.2, 3.6, 4.0, 4.2, 4.8),
    PV2 = c(2.0, 2.5, 2.9, 3.4, 3.5, 4.1, 4.4, 4.7),
    W = c(1.0, 1.1, 0.9, 1.2, 1.0, 1.3, 0.8, 1.4),
    RW1 = c(0.8, 1.3, 1.0, 1.1, 0.9, 1.5, 0.9, 1.2),
    RW2 = c(1.2, 0.9, 0.8, 1.4, 1.1, 1.0, 0.7, 1.6),
    RW3 = c(1.1, 1.0, 1.1, 1.0, 1.3, 1.2, 0.8, 1.3),
    RW4 = c(0.9, 1.2, 0.95, 1.3, 1.0, 1.4, 1.1, 1.1)
  )
}

classes_design_fixture <- function() {
  pvstackr:::new_pvstackr_design(
    data = classes_fixture_data(),
    formula = OUTCOME ~ x + (1 | school),
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    roles = list(outcome_placeholder = "OUTCOME"),
    provenance = list(source = "test")
  )
}

classes_direct_design_fixture <- function() {
  pvstackr:::new_pvstackr_design(
    data = classes_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    roles = list(outcome_placeholder = "OUTCOME"),
    provenance = list(source = "test")
  )
}

classes_brr_target_fixture <- function() {
  pv_brr_target(
    data = classes_fixture_data(),
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )
}

classes_fake_stack_fit <- function(record) {
  function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    record$n <- record$n + 1L
    draws <- cbind(
      b_Intercept = seq(0.1, 0.8, length.out = 8),
      b_x = seq(1.2, 0.5, length.out = 8),
      sigma = seq(0.8, 1.1, length.out = 8),
      lp__ = seq(-4, -3, length.out = 8)
    )
    list(draws = draws, log_lik = matrix(0, nrow = 8, ncol = nrow(data)))
  }
}

classes_stack_fit_fixture <- function(formula = OUTCOME ~ x + (1 | school)) {
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  pvstackr:::pv_stack_fit(
    data = classes_fixture_data(),
    formula = formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    control = pv_control(backend = "injected", keep_log_lik = TRUE),
    fit_function = classes_fake_stack_fit(record),
    draws_function = function(fit) fit$draws,
    extract_log_lik = TRUE,
    log_lik_function = function(fit) fit$log_lik
  )
}

classes_ccc_fixture <- function() {
  target <- classes_brr_target_fixture()
  beta <- target$beta
  draws <- cbind(
    b_Intercept = beta[["b_Intercept"]] + c(-0.03, -0.01, 0.01, 0.03),
    b_x = beta[["b_x"]] + c(-0.02, 0.02, -0.01, 0.01),
    sigma = c(0.8, 0.9, 1.0, 1.1)
  )
  pvstackr:::ccc_calibrate(draws, target)
}

test_that("pvstackr_design constructor records schema and stable hashes", {
  design <- classes_design_fixture()

  expect_s3_class(design, "pvstackr_design")
  expect_equal(design$M, 2L)
  expect_equal(design$R, 4L)
  expect_equal(design$n, nrow(classes_fixture_data()))
  expect_equal(design$pv_cols, c("PV1", "PV2"))
  expect_equal(design$rep_weight_cols, paste0("RW", 1:4))
  expect_equal(design$id_cols, "id")
  expect_equal(design$formula_string, "OUTCOME ~ x + (1 | school)")
  expect_equal(design$data_manifest$n, nrow(classes_fixture_data()))
  expect_equal(design$data_manifest$pv_cols, c("PV1", "PV2"))
  expect_match(design$row_support_hash, "^[0-9a-f]{8}$")
  expect_match(design$pv_value_hash, "^[0-9a-f]{8}$")
  expect_match(design$weight_design_hash, "^[0-9a-f]{8}$")
  expect_match(design$design_hash, "^[0-9a-f]{8}$")
  expect_equal(design$provenance$function_name, "new_pvstackr_design")
  expect_equal(design$provenance$schema_version, design$schema_version)
  expect_invisible(pvstackr:::validate_pvstackr_design(design))
})

test_that("pvstackr_design validator rejects malformed schema and stale hashes", {
  data <- classes_fixture_data()

  expect_error(
    pvstackr:::new_pvstackr_design(data, PV1 ~ x, c("PV1", "PV2")),
    "OUTCOME"
  )
  expect_error(
    pvstackr:::new_pvstackr_design(data, OUTCOME ~ x + weights(W), c("PV1", "PV2")),
    "weights"
  )
  expect_s3_class(
    pvstackr:::new_pvstackr_design(data, OUTCOME ~ x + I((x > 0) | (id > 2)), c("PV1", "PV2")),
    "pvstackr_design"
  )

  bad <- data
  bad$id[2] <- bad$id[1]
  expect_error(
    pvstackr:::new_pvstackr_design(bad, OUTCOME ~ x, c("PV1", "PV2"), id_cols = "id"),
    "unique rows"
  )

  bad <- data
  bad$W[1] <- 0
  expect_error(
    pvstackr:::new_pvstackr_design(bad, OUTCOME ~ x, c("PV1", "PV2"), weight_col = "W"),
    "strictly positive"
  )

  expect_error(
    pvstackr:::new_pvstackr_design(data, OUTCOME ~ x, c("PV1", "PV2"), rep_weight_cols = "RW1"),
    "at least two"
  )

  design <- classes_design_fixture()
  design$data$PV1[1] <- design$data$PV1[1] + 1
  expect_error(pvstackr:::validate_pvstackr_design(design), "hashes")
})

test_that("schema validators accept current target, stack-fit, and CCC objects", {
  target <- classes_brr_target_fixture()
  stack_fit <- classes_stack_fit_fixture()
  ccc <- classes_ccc_fixture()

  expect_invisible(pvstackr:::validate_pvstackr_brr_target(target))
  expect_invisible(pvstackr:::validate_pvstackr_stack_fit(stack_fit))
  expect_invisible(pvstackr:::validate_pvstackr_ccc(ccc))
})

test_that("schema validators reject corrupted target, stack-fit, and CCC objects", {
  target <- classes_brr_target_fixture()
  expect_s3_class(pvstackr:::new_pvstackr_brr_target(target), "pvstackr_brr_target")

  bad_target <- target
  bad_target$target_source <- "raw posterior covariance"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "external_brr_fay_rubin")

  bad_target <- target
  bad_target$M <- 3L
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "per_pv")

  bad_target <- target
  rownames(bad_target$U_bar)[1] <- "wrong"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "U_bar")

  bad_target <- target
  bad_target$target_hash <- "not-a-hash"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "8-character")

  bad_target <- target
  bad_target$policy$target_repair <- "forbidden_by_default"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "target repair")

  bad_target <- target
  bad_target$df_method <- "unsupported"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "df_method")

  bad_target <- target
  bad_target$df_complete <- setNames(rep(20, length(target$fe_names)), target$fe_names)
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "classic")

  bad_target <- target
  bad_target$coverage_claim_allowed <- TRUE
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "coverage_claim_allowed")

  bad_target <- target
  bad_target$interval_role <- "coverage_barnard_rubin"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "interval_role")

  bad_target <- target
  bad_target$policy$df_method <- "barnard_rubin"
  expect_error(pvstackr:::validate_pvstackr_brr_target(bad_target), "policy df metadata")

  stack_fit <- classes_stack_fit_fixture()
  bad_stack <- stack_fit
  bad_stack$psi_hat_fe <- unname(bad_stack$psi_hat_fe)
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "psi_hat_fe")

  bad_stack <- stack_fit
  bad_stack$psi_hat_fe[1] <- bad_stack$psi_hat_fe[1] + 1
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "draw means")

  bad_stack <- stack_fit
  bad_stack$provenance$long_data_hash <- "00000000"
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "provenance hash")

  bad_stack <- stack_fit
  bad_stack["log_lik"] <- list(NULL)
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "log-likelihood retention")

  bad_stack <- stack_fit
  bad_stack$weight_summary$n_long <- bad_stack$weight_summary$n_long + 1
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "long-data dimensions")

  bad_stack <- stack_fit
  bad_stack$meta$prior_policy$non_flat_prior_warning <- TRUE
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "legacy prior diagnostic")

  bad_stack <- stack_fit
  bad_stack$meta$prior_diagnostic$warn_nonflat_prior <- TRUE
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "legacy prior diagnostic")

  bad_stack <- stack_fit
  bad_stack$meta$prior_policy$explicit_prior_warning <- TRUE
  bad_stack$meta$prior_policy$reason_code <- "wrong_prior_code"
  expect_error(pvstackr:::validate_pvstackr_stack_fit(bad_stack), "explicit_prior_warning")

  ccc <- classes_ccc_fixture()
  bad_ccc <- ccc
  bad_ccc$target_source <- "mock"
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "external_brr_fay_rubin")

  bad_ccc <- ccc
  bad_ccc$draws_fe_cal[1, 1] <- bad_ccc$draws_fe_cal[1, 1] + 1
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "draws_fe_cal")

  bad_ccc <- ccc
  bad_ccc$A_full[3, 3] <- 2
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "A_full")

  bad_ccc <- ccc
  bad_ccc$flags$target_repaired <- TRUE
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "target repair")

  bad_ccc <- ccc
  bad_ccc$policy$target_repair <- "allowed"
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "target repair")

  bad_ccc <- ccc
  bad_ccc$control$allow_target_nearpd <- TRUE
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "not currently supported")

  bad_ccc <- ccc
  bad_ccc$flags$nearpd_target <- TRUE
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "nearPD target")

  bad_ccc <- ccc
  bad_ccc$ccc_status <- "blocked"
  expect_error(pvstackr:::validate_pvstackr_ccc(bad_ccc), "ccc_status")
})

classes_fit_estimates <- function(ccc, target = classes_brr_target_fixture()) {
  pvstackr:::pv_fit_direct_estimates(ccc, target, conf_level = 0.95)
}

classes_reportable_fit_fixture <- function(status = "ok", reason_codes = character(), warnings = character()) {
  design <- classes_direct_design_fixture()
  target <- classes_brr_target_fixture()
  stack_fit <- classes_stack_fit_fixture(OUTCOME ~ x)
  ccc <- classes_ccc_fixture()
  pvstackr:::new_pvstackr_fit(
    method = "stack_direct",
    design = design,
    target = target,
    stack_fit = stack_fit,
    ccc = ccc,
    estimates = classes_fit_estimates(ccc, target),
    draws = ccc$draws_fe_cal,
    diagnostics = list(ccc = ccc$diagnostics),
    status = status,
    control = pv_control(method = "stack_direct"),
    reason_codes = reason_codes,
    provenance = list(source = "test"),
    warnings = warnings
  )
}

test_that("pvstackr_fit constructor assembles method objects with fixed-effect draws only", {
  fit <- classes_reportable_fit_fixture()

  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$method, "stack_direct")
  expect_equal(fit$status, "ok")
  expect_equal(fit$provenance$function_name, "new_pvstackr_fit")
  expect_equal(colnames(fit$draws), c("b_Intercept", "b_x"))
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))
})

test_that("pvstackr_fit constructor enforces stack_direct reportable status invariants", {
  design <- classes_direct_design_fixture()
  target <- classes_brr_target_fixture()
  stack_fit <- classes_stack_fit_fixture(OUTCOME ~ x)
  ccc <- classes_ccc_fixture()
  estimates <- classes_fit_estimates(ccc, target)

  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct"),
    "Reportable stack_direct fit requires"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = design,
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = data.frame(),
      control = pv_control(method = "stack_direct")
    ),
    "non-empty estimates"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = design,
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = data.frame(term = rev(names(ccc$psi_hat)), estimate = 1:2),
      control = pv_control(method = "stack_direct")
    ),
    "estimate terms"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = classes_design_fixture(),
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = estimates,
      control = pv_control(method = "stack_direct")
    ),
    "Random-effect/group terms"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = design,
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = estimates,
      status = "blocked",
      reason_codes = "ccc_conditioning_red",
      control = pv_control(method = "stack_direct")
    ),
    "Blocked stack_direct fit must not include reportable estimates"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      method = "stack_direct",
      design = design,
      target = target,
      stack_fit = stack_fit,
      ccc = ccc,
      estimates = estimates,
      control = pv_control(method = "stack_direct", center = "posterior")
    ),
    "control\\$center = \"target\""
  )

  fit <- classes_reportable_fit_fixture()
  fit["ccc"] <- list(NULL)
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "Reportable stack_direct fit requires"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates <- data.frame(term = "b_x", estimate = 1)
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "estimate terms"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates$target_hash <- "00000000"
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "estimate target metadata"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates$target_source <- "mock"
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "estimate target metadata"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates$parameter_scope <- "variance_component"
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "parameter_scope"
  )

  fit <- classes_reportable_fit_fixture()
  fit$estimates$estimate[1] <- fit$estimates$estimate[1] + 1
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "external target"
  )

  fit <- classes_reportable_fit_fixture()
  fit$ccc$center <- "posterior"
  fit$ccc$control$center <- "posterior"
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "target-centered CCC"
  )

  fit <- classes_reportable_fit_fixture()
  fit$ccc$diagnostics$center_status <- "blocked"
  fit$ccc$diagnostics$center_reason_code <- "center_separation_red"
  fit$ccc$diagnostics$center_separation$band <- "red"
  fit$ccc$diagnostics$center_separation$reason_code <- "center_separation_red"
  fit$ccc$flags$center_separation_blocked <- TRUE
  fit$ccc$warnings <- "Center separation diagnostic exceeded the block threshold."
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "blocked center separation"
  )

  fit <- classes_reportable_fit_fixture()
  fit$ccc$diagnostics$center_status <- "warning"
  fit$ccc$diagnostics$center_reason_code <- "center_separation_yellow"
  fit$ccc$diagnostics$center_separation$band <- "yellow"
  fit$ccc$diagnostics$center_separation$reason_code <- "center_separation_yellow"
  fit$ccc$flags$center_separation_warning <- TRUE
  fit$ccc$warnings <- "Center separation diagnostic exceeded the warning threshold."
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "status = \"warning\""
  )

  fit$status <- "warning"
  fit$reason_codes <- "center_separation_yellow"
  fit$warnings <- "Center separation diagnostic exceeded the warning threshold."
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  fit <- classes_reportable_fit_fixture()
  fit$ccc$diagnostics$conditioning_status <- "blocked"
  fit$ccc$diagnostics$conditioning_reason_code <- "ccc_conditioning_red"
  fit$ccc$diagnostics$conditioning_band <- "red"
  fit$ccc$diagnostics$conditioning$band <- "red"
  fit$ccc$diagnostics$conditioning$reason_code <- "ccc_conditioning_red"
  fit$ccc$flags$conditioning_blocked <- TRUE
  fit$ccc$flags$kappa_a_blocked <- TRUE
  fit$ccc$warnings <- "CCC conditioning diagnostic exceeded the block threshold."
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "blocked CCC conditioning"
  )

  fit <- classes_reportable_fit_fixture()
  fit$ccc$diagnostics$conditioning_status <- "warning"
  fit$ccc$diagnostics$conditioning_reason_code <- "ccc_conditioning_yellow"
  fit$ccc$diagnostics$conditioning_band <- "yellow"
  fit$ccc$diagnostics$conditioning$band <- "yellow"
  fit$ccc$diagnostics$conditioning$reason_code <- "ccc_conditioning_yellow"
  fit$ccc$flags$conditioning_warning <- TRUE
  fit$ccc$flags$kappa_a_warning <- TRUE
  fit$ccc$warnings <- "CCC conditioning diagnostic exceeded the warning threshold."
  expect_error(
    pvstackr:::validate_pvstackr_fit(fit),
    "status = \"warning\""
  )

  fit$status <- "warning"
  fit$reason_codes <- "ccc_conditioning_yellow"
  fit$warnings <- "CCC conditioning diagnostic exceeded the warning threshold."
  expect_invisible(pvstackr:::validate_pvstackr_fit(fit))

  ccc <- classes_reportable_fit_fixture()$ccc
  ccc$flags$conditioning_warning <- TRUE
  expect_error(
    pvstackr:::validate_pvstackr_ccc(ccc),
    "conditioning flags"
  )
})

test_that("pvstackr_fit status metadata distinguishes ok, warning, and blocked objects", {
  expect_error(
    classes_reportable_fit_fixture(status = "ok", reason_codes = "unexpected"),
    "must not include `reason_codes`"
  )
  expect_error(
    classes_reportable_fit_fixture(status = "ok", warnings = "unexpected"),
    "must not include `warnings`"
  )
  expect_error(
    classes_reportable_fit_fixture(status = "warning"),
    "reason_codes"
  )
  expect_error(
    classes_reportable_fit_fixture(status = "warning", reason_codes = "center_separation_yellow"),
    "warnings"
  )

  warning_fit <- classes_reportable_fit_fixture(
    status = "warning",
    reason_codes = "center_separation_yellow",
    warnings = "Center separation diagnostic exceeded the advisory threshold."
  )
  expect_s3_class(warning_fit, "pvstackr_fit")
  expect_equal(warning_fit$status, "warning")
  expect_equal(warning_fit$reason_codes, "center_separation_yellow")
  expect_invisible(pvstackr:::validate_pvstackr_fit(warning_fit))

  expect_s3_class(
    pvstackr:::new_pvstackr_fit(
      "stack_direct",
      status = "blocked",
      reason_codes = "preflight_failed",
      control = pv_control(method = "stack_direct")
    ),
    "pvstackr_fit"
  )
})

test_that("pvstackr_fit constructor enforces method, status, control, and draw contracts", {
  expect_error(
    pvstackr:::new_pvstackr_fit("pipeline_c"),
    "internal manuscript"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("per_pv"),
    "pvstackr_reference_pool"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_psis"),
    "PSIS"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", control = pv_control(method = "per_pv")),
    "control\\$method"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", status = "blocked"),
    "reason_codes"
  )
  expect_error(
    pvstackr:::new_pvstackr_fit(
      "stack_psis",
      status = "blocked",
      reason_codes = "psis_k_too_high",
      control = pv_control(method = "stack_psis")
    ),
    "PSIS diagnostics"
  )
  blocked <- pvstackr:::new_pvstackr_fit(
    "stack_direct",
    status = "blocked",
    reason_codes = "preflight_failed",
    control = pv_control(method = "stack_direct")
  )
  expect_s3_class(blocked, "pvstackr_fit")
  expect_equal(blocked$reason_codes, "preflight_failed")

  blocked$estimates <- data.frame(term = "b_x", estimate = 1)
  expect_error(
    pvstackr:::validate_pvstackr_fit(blocked),
    "Blocked stack_direct fit must not include reportable estimates"
  )

  bad_draws <- cbind(b_Intercept = c(1, 2, 3), sigma = c(0.8, 0.9, 1.0))
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", draws = bad_draws),
    "fixed-effect columns only"
  )

  bad_estimates <- data.frame(term = "sigma", estimate = 1)
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", estimates = bad_estimates),
    "fixed-effect terms"
  )

  target <- classes_brr_target_fixture()
  ccc <- classes_ccc_fixture()
  ccc$target_hash <- "00000000"
  expect_error(
    pvstackr:::new_pvstackr_fit("stack_direct", target = target, ccc = ccc),
    "same `target_hash`"
  )
})
