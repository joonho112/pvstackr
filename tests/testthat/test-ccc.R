ccc_target_fixture <- function(beta = c(b_Intercept = 5, b_x = -3),
                               Sigma = diag(c(1, 16))) {
  dimnames(Sigma) <- list(names(beta), names(beta))
  structure(
    list(
      beta = beta,
      T_MI = Sigma,
      fe_names = names(beta),
      target_source = "external_brr_fay_rubin",
      target_hash = "fixture-external"
    ),
    class = c("pvstackr_brr_target", "list")
  )
}

ccc_draw_fixture <- function(include_vc = TRUE) {
  raw_center <- c(b_Intercept = 10, b_x = 20)
  raw_centered <- rbind(
    c(sqrt(6), 0),
    c(-sqrt(6), 0),
    c(0, sqrt(13.5)),
    c(0, -sqrt(13.5))
  )
  colnames(raw_centered) <- names(raw_center)
  fe <- sweep(raw_centered, 2L, raw_center, FUN = "+")
  if (!include_vc) {
    return(fe)
  }
  cbind(fe, sigma = c(1.1, 1.2, 1.3, 1.4))
}

ccc_bound_fixture <- function() {
  data <- data.frame(
    id = seq_len(10),
    x = c(-2, -1.5, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3),
    PV1 = c(2.1, 2.4, 3.0, 3.2, 3.6, 4.0, 4.2, 4.8, 5.1, 5.5),
    PV2 = c(2.0, 2.5, 2.9, 3.4, 3.5, 4.1, 4.4, 4.7, 5.0, 5.8),
    W = c(1.0, 1.1, 0.9, 1.2, 1.0, 1.3, 0.8, 1.4, 1.1, 0.95),
    RW1 = c(0.8, 1.3, 1.0, 1.1, 0.9, 1.5, 0.9, 1.2, 1.2, 1.0),
    RW2 = c(1.2, 0.9, 0.8, 1.4, 1.1, 1.0, 0.7, 1.6, 1.0, 0.9),
    RW3 = c(1.1, 1.0, 1.1, 1.0, 1.3, 1.2, 0.8, 1.3, 1.4, 1.1),
    RW4 = c(0.9, 1.2, 0.95, 1.3, 1.0, 1.4, 1.1, 1.1, 0.9, 1.2)
  )
  target <- pv_brr_target(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )
  preflight <- pvstackr:::pv_stack_direct_preflight(
    data = data,
    formula = OUTCOME ~ x,
    target = target
  )
  list(
    data = data,
    target = target,
    preflight = preflight,
    draws = ccc_draw_fixture()
  )
}

test_that("CCC calibration matches target mean and covariance exactly for fixed effects", {
  draws <- ccc_draw_fixture()
  target <- ccc_target_fixture()
  out <- pvstackr:::ccc_calibrate(draws, target)

  expect_s3_class(out, "pvstackr_ccc")
  expect_identical(out$ccc_status, "ok")
  expect_identical(out$schema_version, "0.1.0")
  expect_false("binding_proof" %in% names(out))
  expect_invisible(pvstackr:::validate_pvstackr_ccc(out))
  expect_equal(out$psi_hat, target$beta, tolerance = 1e-14)
  expect_equal(colMeans(out$draws_fe_cal), target$beta, tolerance = 1e-14)
  expect_equal(unname(out$Sigma_raw), diag(c(4, 9)), tolerance = 1e-14)
  expect_equal(out$Sigma_cal_emp_raw, target$T_MI, tolerance = 1e-14)
  expect_equal(out$Sigma_cal_emp, target$T_MI, tolerance = 1e-14)
  expect_false(out$flags$target_repaired)
  expect_identical(out$control$center, "target")
  expect_identical(out$provenance$function_name, "ccc_calibrate")
  expect_lt(out$diagnostics$rho1, 1e-12)
  expect_lt(out$diagnostics$rho2, 1e-12)
  expect_lt(out$diagnostics$empirical_fro_rel, 1e-12)
  expect_equal(
    out$diagnostics$delta_c_rel,
    sqrt(mean(((target$beta - colMeans(draws[, target$fe_names])) / sqrt(diag(target$T_MI)))^2)),
    tolerance = 1e-14
  )
  expect_equal(
    out$diagnostics$delta_c_max,
    max(abs((target$beta - colMeans(draws[, target$fe_names])) / sqrt(diag(target$T_MI)))),
    tolerance = 1e-14
  )
  expect_equal(out$psi_target, target$beta, tolerance = 0)
  expect_equal(out$diagnostics$a_matrix_fro_rel, sqrt(((0.5 - 1)^2 + ((4 / 3) - 1)^2) / 2), tolerance = 1e-14)
  expect_equal(out$diagnostics$center_separation$raw_center, colMeans(draws[, target$fe_names]), tolerance = 1e-14)
  expect_equal(out$diagnostics$center_separation$target_center, target$beta, tolerance = 0)
  expect_identical(out$diagnostics$center_separation$gate_metric, "delta_c_max")
  expect_equal(out$diagnostics$center_separation$gate_value, out$diagnostics$delta_c_max, tolerance = 0)
  expect_identical(out$diagnostics$center_separation$band, "red")
  expect_identical(out$diagnostics$center_separation$reason_code, "center_separation_red")
  expect_null(out$diagnostics$delta_c_rel_old)
})

test_that("proof-bearing CCC is schema 0.2 and preserves the exact verified proof", {
  fixture <- ccc_bound_fixture()
  out <- pvstackr:::ccc_calibrate(
    fixture$draws,
    fixture$target,
    binding_proof = fixture$preflight$binding_proof
  )

  expect_identical(out$schema_version, "0.2.0")
  expect_identical(names(out), pvstackr:::pv_ccc_schema02_fields())
  expect_identical(out$binding_proof, fixture$preflight$binding_proof)
  expect_identical(out$target_hash, fixture$target$target_content$target_content_hash)
  expect_invisible(pvstackr:::validate_pvstackr_ccc(out))
})

test_that("CCC validates a supplied proof before inspecting calibration payload", {
  fixture <- ccc_bound_fixture()
  stale <- fixture$preflight$binding_proof
  stale$current_manifest_hash <- paste0("sha256:", strrep("0", 64))

  error <- tryCatch(
    pvstackr:::ccc_calibrate(
      draws = "not-a-draw-payload",
      target = fixture$target,
      binding_proof = stale
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(error, "pvstackr_binding_error")
  expect_identical(error$code, "PV_BIND_E005")

  legacy_error <- tryCatch(
    pvstackr:::ccc_calibrate(
      draws = "not-a-draw-payload",
      target = ccc_target_fixture(),
      binding_proof = fixture$preflight$binding_proof
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(legacy_error, "pvstackr_binding_error")
  expect_identical(legacy_error$code, "PV_BIND_E080")

  stale_target <- fixture$target
  stale_target$target_hash <- paste0("sha256:", strrep("1", 64))
  target_error <- tryCatch(
    pvstackr:::ccc_calibrate(
      draws = "not-a-draw-payload",
      target = stale_target,
      binding_proof = fixture$preflight$binding_proof
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(target_error, "pvstackr_binding_error")
  expect_identical(target_error$code, "PV_BIND_E090")

  secret <- "PRIVATE_ATOMIC_TARGET_CONTENT_MUST_NOT_LEAK"
  atomic_target <- fixture$target
  atomic_target$target_content <- secret
  atomic_error <- tryCatch(
    pvstackr:::ccc_calibrate(
      draws = "not-a-draw-payload",
      target = atomic_target,
      binding_proof = fixture$preflight$binding_proof
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(atomic_error, "pvstackr_binding_error")
  expect_identical(atomic_error$code, "PV_BIND_E090")
  expect_false(grepl(secret, conditionMessage(atomic_error), fixed = TRUE))
  expect_false(grepl(
    secret,
    rawToChar(serialize(atomic_error, NULL, ascii = TRUE)),
    fixed = TRUE
  ))
})

test_that("CCC computes the expected diagonal calibration matrix", {
  draws <- ccc_draw_fixture()
  target <- ccc_target_fixture()
  out <- pvstackr:::ccc_calibrate(draws, target)
  A_hand <- diag(c(0.5, 4 / 3))

  expect_equal(unname(out$A), A_hand, tolerance = 1e-14)
  expect_equal(unname(out$A_full[1:2, 1:2]), A_hand, tolerance = 1e-14)
  expect_equal(unname(out$A_full[3, 3]), 1, tolerance = 0)
})

test_that("CCC preserves non-fixed-effect draw columns verbatim", {
  draws <- ccc_draw_fixture(include_vc = TRUE)
  target <- ccc_target_fixture()
  out <- pvstackr:::ccc_calibrate(
    draws,
    target,
    param_map = list(fe_names = c("b_Intercept", "b_x"), vc_names = "sigma")
  )

  expect_equal(out$draws_calibrated[, "sigma"], draws[, "sigma"], tolerance = 0)
  expect_true(out$policy$vc_passthrough)
  expect_identical(out$param_map$vc_names, "sigma")
})

test_that("CCC accepts numeric fixed-effect and passthrough column maps", {
  base_draws <- ccc_draw_fixture(include_vc = TRUE)
  draws <- cbind(
    sigma = base_draws[, "sigma"],
    b_x = base_draws[, "b_x"],
    b_Intercept = base_draws[, "b_Intercept"]
  )
  target <- ccc_target_fixture()
  out <- pvstackr:::ccc_calibrate(
    draws,
    target,
    param_map = list(fe_idx = c(3L, 2L), vc_idx = 1L)
  )

  expect_identical(out$param_map$fe_names, target$fe_names)
  expect_identical(out$param_map$vc_names, "sigma")
  expect_equal(out$draws_fe_cal, out$draws_calibrated[, target$fe_names], tolerance = 0)
  expect_equal(colMeans(out$draws_fe_cal), target$beta, tolerance = 1e-14)
  expect_equal(out$Sigma_cal_emp_raw, target$T_MI, tolerance = 1e-14)
  expect_equal(out$draws_calibrated[, "sigma"], draws[, "sigma"], tolerance = 0)
  expect_equal(out$A_full[1L, -1L], c(b_x = 0, b_Intercept = 0), tolerance = 0)
  expect_equal(out$A_full[-1L, 1L], c(b_x = 0, b_Intercept = 0), tolerance = 0)
  expect_equal(out$A_full[1L, 1L], 1, tolerance = 0)
})

test_that("CCC supports off-diagonal target covariance", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  Sigma <- matrix(c(1, 0.4, 0.4, 2), 2, 2)
  target <- ccc_target_fixture(Sigma = Sigma)
  out <- pvstackr:::ccc_calibrate(draws, target)

  algebraic <- pvstackr:::pv_symmetrize(out$A %*% out$Sigma_raw %*% t(out$A))
  expect_equal(algebraic, target$T_MI, tolerance = 1e-12)
  expect_equal(out$Sigma_cal_emp_raw, target$T_MI, tolerance = 1e-12)
  expect_lt(max(abs(out$diagnostics$matrix_residual)), 1e-12)
})

test_that("CCC can retain posterior center while matching target covariance", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  target <- ccc_target_fixture()
  out <- pvstackr:::ccc_calibrate(draws, target, center = "posterior")

  expect_equal(out$psi_hat, colMeans(draws), tolerance = 1e-14)
  expect_equal(out$psi_target, target$beta, tolerance = 0)
  expect_equal(colMeans(out$draws_fe_cal), colMeans(draws), tolerance = 1e-14)
  expect_equal(out$Sigma_cal_emp_raw, target$T_MI, tolerance = 1e-14)
  expect_equal(
    out$diagnostics$delta_c_rel,
    sqrt(mean(((target$beta - colMeans(draws)) / sqrt(diag(target$T_MI)))^2)),
    tolerance = 1e-14
  )
  expect_equal(
    out$diagnostics$delta_c_max,
    max(abs((target$beta - colMeans(draws)) / sqrt(diag(target$T_MI)))),
    tolerance = 1e-14
  )
})

test_that("CCC calibration is equivariant to common fixed-effect translations", {
  draws <- ccc_draw_fixture(include_vc = TRUE)
  target <- ccc_target_fixture()
  base <- pvstackr:::ccc_calibrate(draws, target)

  shift <- c(b_Intercept = 100, b_x = -50)
  shifted_draws <- draws
  shifted_draws[, names(shift)] <- sweep(shifted_draws[, names(shift), drop = FALSE], 2L, shift, FUN = "+")
  shifted_target <- target
  shifted_target$beta <- target$beta + shift

  shifted <- pvstackr:::ccc_calibrate(shifted_draws, shifted_target)

  expect_equal(shifted$A, base$A, tolerance = 1e-14)
  expect_equal(shifted$diagnostics$delta_c_rel, base$diagnostics$delta_c_rel, tolerance = 1e-14)
  expect_equal(
    sweep(shifted$draws_fe_cal, 2L, shift, FUN = "-"),
    base$draws_fe_cal,
    tolerance = 1e-14
  )
  expect_equal(shifted$draws_calibrated[, "sigma"], base$draws_calibrated[, "sigma"], tolerance = 0)
})

test_that("CCC center diagnostics flag warning and blocked separations", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  raw_center <- colMeans(draws)
  zero_target <- ccc_target_fixture(
    beta = raw_center,
    Sigma = diag(c(1, 16))
  )
  zero_out <- pvstackr:::ccc_calibrate(draws, zero_target)
  expect_equal(zero_out$diagnostics$delta_c_rel, 0, tolerance = 1e-14)
  expect_identical(zero_out$diagnostics$center_status, "ok")
  expect_identical(zero_out$diagnostics$center_separation$band, "green")
  expect_true(is.na(zero_out$diagnostics$center_separation$reason_code))

  small_shift <- c(b_Intercept = 0.02, b_x = 0)
  warning_target <- ccc_target_fixture(
    beta = raw_center + small_shift,
    Sigma = diag(c(1, 16))
  )
  warning_out <- pvstackr:::ccc_calibrate(draws, warning_target)

  expect_identical(warning_out$diagnostics$center_status, "warning")
  expect_identical(warning_out$diagnostics$center_reason_code, "center_separation_yellow")
  expect_identical(warning_out$diagnostics$center_separation$band, "yellow")
  expect_true(warning_out$flags$center_separation_warning)
  expect_false(warning_out$flags$center_separation_blocked)
  expect_match(warning_out$warnings, "warning threshold")

  blocked_target <- ccc_target_fixture(
    beta = raw_center + c(b_Intercept = 0.1, b_x = 0),
    Sigma = diag(c(1, 16))
  )
  blocked_out <- pvstackr:::ccc_calibrate(draws, blocked_target)

  expect_identical(blocked_out$diagnostics$center_status, "blocked")
  expect_identical(blocked_out$diagnostics$center_reason_code, "center_separation_red")
  expect_identical(blocked_out$diagnostics$center_separation$band, "red")
  expect_false(blocked_out$flags$center_separation_warning)
  expect_true(blocked_out$flags$center_separation_blocked)
  expect_match(blocked_out$warnings, "block threshold")
})

test_that("CCC center diagnostics block when any term exceeds the red max threshold", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  raw_center <- colMeans(draws)

  target <- ccc_target_fixture(
    beta = raw_center + c(b_Intercept = 0.06, b_x = 0),
    Sigma = diag(c(1, 16))
  )
  out <- pvstackr:::ccc_calibrate(draws, target)

  expect_equal(out$diagnostics$delta_c_by_term, c(b_Intercept = 0.06, b_x = 0), tolerance = 1e-14)
  expect_equal(out$diagnostics$delta_c_rel, sqrt(0.06^2 / 2), tolerance = 1e-14)
  expect_equal(out$diagnostics$delta_c_max, 0.06, tolerance = 1e-14)
  expect_lt(out$diagnostics$delta_c_rel, out$diagnostics$center_threshold_block)
  expect_gt(out$diagnostics$delta_c_max, out$diagnostics$center_threshold_block)
  expect_identical(out$diagnostics$center_separation$gate_metric, "delta_c_max")
  expect_equal(out$diagnostics$center_separation$gate_value, out$diagnostics$delta_c_max, tolerance = 0)
  expect_identical(out$diagnostics$center_status, "blocked")
  expect_identical(out$diagnostics$center_reason_code, "center_separation_red")
  expect_identical(out$diagnostics$center_separation$band, "red")
  expect_identical(out$diagnostics$center_separation$reason_code, "center_separation_red")
  expect_false(out$flags$center_separation_warning)
  expect_true(out$flags$center_separation_blocked)
  expect_match(out$warnings, "delta_c_max")
  expect_match(out$warnings, "block threshold")
})

test_that("CCC center diagnostics use inclusive warning and block bands", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  raw_center <- colMeans(draws)

  warning_target <- ccc_target_fixture(
    beta = raw_center + c(b_Intercept = 0.0100001, b_x = 0.0400004),
    Sigma = diag(c(1, 16))
  )
  warning_out <- pvstackr:::ccc_calibrate(draws, warning_target)
  expect_gte(warning_out$diagnostics$delta_c_max, 0.01)
  expect_identical(warning_out$diagnostics$center_status, "warning")
  expect_identical(warning_out$diagnostics$center_separation$band, "yellow")
  expect_identical(warning_out$diagnostics$center_reason_code, "center_separation_yellow")

  blocked_target <- ccc_target_fixture(
    beta = raw_center + c(b_Intercept = 0.0500001, b_x = 0.2000004),
    Sigma = diag(c(1, 16))
  )
  blocked_out <- pvstackr:::ccc_calibrate(draws, blocked_target)
  expect_gte(blocked_out$diagnostics$delta_c_max, 0.05)
  expect_identical(blocked_out$diagnostics$center_status, "blocked")
  expect_identical(blocked_out$diagnostics$center_separation$band, "red")
  expect_identical(blocked_out$diagnostics$center_reason_code, "center_separation_red")
})

test_that("CCC conditioning diagnostics warn and block on kappa_A", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  raw_center <- colMeans(draws)

  warning_target <- ccc_target_fixture(
    beta = raw_center,
    Sigma = diag(c(4 * 1e12, 9))
  )
  warning_out <- pvstackr:::ccc_calibrate(draws, warning_target)

  expect_equal(unname(warning_out$A), diag(c(1e6, 1)), tolerance = 1e-10)
  expect_equal(warning_out$diagnostics$kappa_A, 1e6, tolerance = 1e-10)
  expect_identical(warning_out$diagnostics$center_status, "ok")
  expect_identical(warning_out$diagnostics$conditioning_status, "warning")
  expect_identical(warning_out$diagnostics$conditioning_reason_code, "ccc_conditioning_yellow")
  expect_identical(warning_out$diagnostics$conditioning$gate_metric, "kappa_A")
  expect_equal(warning_out$diagnostics$conditioning$gate_value, warning_out$diagnostics$kappa_A, tolerance = 0)
  expect_equal(warning_out$diagnostics$kappa_A_threshold_warn, 1e6, tolerance = 0)
  expect_equal(warning_out$diagnostics$kappa_A_threshold_block, 1e8, tolerance = 0)
  expect_true(warning_out$flags$conditioning_warning)
  expect_false(warning_out$flags$conditioning_blocked)
  expect_true(warning_out$flags$kappa_a_warning)
  expect_false(warning_out$flags$kappa_a_blocked)
  expect_match(warning_out$warnings, "CCC conditioning diagnostic")
  expect_match(warning_out$warnings, "warning threshold")
  expect_invisible(pvstackr:::validate_pvstackr_ccc(warning_out))

  blocked_target <- ccc_target_fixture(
    beta = raw_center,
    Sigma = diag(c(4 * 1e16, 9))
  )
  blocked_out <- pvstackr:::ccc_calibrate(draws, blocked_target)

  expect_equal(unname(blocked_out$A), diag(c(1e8, 1)), tolerance = 1e-8)
  expect_equal(blocked_out$diagnostics$kappa_A, 1e8, tolerance = 1e-8)
  expect_identical(blocked_out$diagnostics$center_status, "ok")
  expect_identical(blocked_out$diagnostics$conditioning_status, "blocked")
  expect_identical(blocked_out$diagnostics$conditioning_reason_code, "ccc_conditioning_red")
  expect_identical(blocked_out$diagnostics$conditioning$band, "red")
  expect_false(blocked_out$flags$conditioning_warning)
  expect_true(blocked_out$flags$conditioning_blocked)
  expect_false(blocked_out$flags$kappa_a_warning)
  expect_true(blocked_out$flags$kappa_a_blocked)
  expect_match(blocked_out$warnings, "CCC conditioning diagnostic")
  expect_match(blocked_out$warnings, "block threshold")
  expect_invisible(pvstackr:::validate_pvstackr_ccc(blocked_out))
})

test_that("CCC blocks reusing an already calibrated cloud as a self target", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  raw_cov <- stats::cov(draws)
  target <- ccc_target_fixture(beta = colMeans(draws), Sigma = raw_cov + diag(c(0.01, 0.02)))

  out1 <- pvstackr:::ccc_calibrate(draws, target, center = "posterior")
  expect_error(
    pvstackr:::ccc_calibrate(out1$draws_fe_cal, target, center = "posterior"),
    "raw fixed-effect covariance|proportional"
  )
})

test_that("CCC rejects self, mock, raw, fallback, and missing target provenance", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  target <- ccc_target_fixture()

  bad <- target
  bad$target_source <- "diagnostic_self_target"
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "forbidden")

  bad <- target
  bad$target_source <- "mock fallback target"
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "forbidden")

  bad <- target
  bad$target_source <- "raw posterior covariance"
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "forbidden")

  bad <- target
  bad$target_source <- NULL
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "missing required")

  bad <- target
  bad$target_source <- "external_schema_fixture"
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "external_brr_fay_rubin")

  bad <- target
  attr(bad, "source") <- "mock self target"
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "forbidden")
})

test_that("CCC validates fixed-effect alignment and target definiteness", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  target <- ccc_target_fixture()

  bad <- target
  names(bad$beta) <- c("wrong", "b_x")
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "missing fixed-effect")

  bad <- target
  colnames(bad$T_MI) <- rownames(bad$T_MI) <- c("b_x", "b_Intercept")
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "exactly match")

  bad <- target
  bad$T_MI[2, 2] <- -1
  expect_error(pvstackr:::ccc_calibrate(draws, bad), "positive definite")

  bad_draws <- draws
  colnames(bad_draws) <- NULL
  expect_error(pvstackr:::ccc_calibrate(bad_draws, target), "column names")

  bad_draws <- cbind(draws, draws[, 1, drop = FALSE])
  colnames(bad_draws)[3] <- "b_x"
  expect_error(pvstackr:::ccc_calibrate(bad_draws, target), "unique column names")
})

test_that("CCC rejects singular raw fixed-effect covariance", {
  draws <- cbind(
    b_Intercept = c(1, 2, 3, 4),
    b_x = c(2, 4, 6, 8)
  )
  target <- ccc_target_fixture()
  expect_error(pvstackr:::ccc_calibrate(draws, target), "Sigma_raw")
})

test_that("CCC blocks raw and proportional raw covariance targets", {
  draws <- ccc_draw_fixture(include_vc = FALSE)
  raw_cov <- stats::cov(draws)

  raw_target <- ccc_target_fixture(beta = colMeans(draws), Sigma = raw_cov)
  expect_error(pvstackr:::ccc_calibrate(draws, raw_target), "raw fixed-effect covariance")

  proportional_target <- ccc_target_fixture(beta = colMeans(draws), Sigma = 2 * raw_cov)
  expect_error(pvstackr:::ccc_calibrate(draws, proportional_target), "proportional")
})

test_that("CCC supports a one-dimensional fixed-effect block", {
  draws <- matrix(c(1, 2, 3, 4), ncol = 1)
  colnames(draws) <- "b_Intercept"
  target <- ccc_target_fixture(
    beta = c(b_Intercept = 10),
    Sigma = matrix(9, 1, 1, dimnames = list("b_Intercept", "b_Intercept"))
  )
  out <- pvstackr:::ccc_calibrate(draws, target)
  expect_equal(unname(colMeans(out$draws_fe_cal)), 10, tolerance = 1e-14)
  expect_equal(as.numeric(stats::cov(out$draws_fe_cal)), 9, tolerance = 1e-14)
})

test_that("ccc_as_draw_matrix strips class/attributes and row names to a plain matrix", {
  base <- ccc_draw_fixture(include_vc = TRUE)
  # Mimic an as_draws_matrix() result: a matrix carrying an extra S3 subclass,
  # draw-id row names, and an incidental attribute.
  classed <- structure(
    base,
    dimnames = list(as.character(seq_len(nrow(base))), colnames(base)),
    class = c("draws_matrix", "draws", "matrix", "array"),
    extra_attr = "provenance"
  )
  expect_true(is.matrix(classed))

  out <- pvstackr:::ccc_as_draw_matrix(classed)
  expect_identical(class(out), c("matrix", "array"))
  expect_null(rownames(out))
  expect_identical(colnames(out), colnames(base))
  expect_null(attr(out, "extra_attr"))
  expect_identical(out, `dimnames<-`(base, list(NULL, colnames(base))))
})

test_that("CCC calibrates and validates a draws_matrix-classed input identically to a plain matrix", {
  base <- ccc_draw_fixture(include_vc = TRUE)
  target <- ccc_target_fixture()
  classed <- structure(
    base,
    dimnames = list(as.character(seq_len(nrow(base))), colnames(base)),
    class = c("draws_matrix", "draws", "matrix", "array")
  )

  plain_out <- pvstackr:::ccc_calibrate(base, target)
  classed_out <- pvstackr:::ccc_calibrate(classed, target)

  # The classed source used to fail validation because draws_calibrated stayed a
  # draws_matrix while draws_fe_cal was a plain matrix built by matrix algebra,
  # so the strict identical() block check compared class, not values.
  expect_invisible(pvstackr:::validate_pvstackr_ccc(classed_out))
  expect_identical(class(classed_out$draws_calibrated), c("matrix", "array"))
  expect_identical(classed_out$draws_calibrated, plain_out$draws_calibrated)
  expect_identical(classed_out$draws_fe_cal, plain_out$draws_fe_cal)
})
