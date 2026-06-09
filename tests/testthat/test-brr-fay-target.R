brr_fixture_data <- function() {
  data.frame(
    id = seq_len(10),
    x = c(-2, -1.5, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3),
    PV1READ = c(2.1, 2.4, 3.0, 3.2, 3.6, 4.0, 4.2, 4.8, 5.1, 5.5),
    PV2READ = c(2.0, 2.5, 2.9, 3.4, 3.5, 4.1, 4.4, 4.7, 5.0, 5.8),
    W_FSTUWT = c(1.0, 1.1, 0.9, 1.2, 1.0, 1.3, 0.8, 1.4, 1.1, 0.95),
    W_FSTURWT1 = c(0.8, 1.3, 1.0, 1.1, 0.9, 1.5, 0.9, 1.2, 1.2, 1.0),
    W_FSTURWT2 = c(1.2, 0.9, 0.8, 1.4, 1.1, 1.0, 0.7, 1.6, 1.0, 0.9),
    W_FSTURWT3 = c(1.1, 1.0, 1.1, 1.0, 1.3, 1.2, 0.8, 1.3, 1.4, 1.1),
    W_FSTURWT4 = c(0.9, 1.2, 0.95, 1.3, 1.0, 1.4, 1.1, 1.1, 0.9, 1.2)
  )
}

hand_wls_beta <- function(data, formula, weights) {
  frame <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
  x <- stats::model.matrix(stats::terms(formula), data = frame)
  y <- stats::model.response(frame)
  beta <- stats::coef(stats::lm.wfit(x, y, weights))
  names(beta) <- c("b_Intercept", "b_x")
  beta
}

hand_one_pv <- function(data, pv_col, fay_k = 0.5) {
  rep_cols <- paste0("W_FSTURWT", 1:4)
  formula <- stats::as.formula(paste(pv_col, "~ x"))
  beta0 <- hand_wls_beta(data, formula, data$W_FSTUWT)
  diff <- vapply(rep_cols, function(col) {
    hand_wls_beta(data, formula, data[[col]]) - beta0
  }, numeric(length(beta0)))
  multiplier <- 1 / (length(rep_cols) * (1 - fay_k)^2)
  U <- multiplier * tcrossprod(diff)
  dimnames(U) <- list(names(beta0), names(beta0))
  list(beta = beta0, U = U)
}

survey_oracle_fe_names <- function(x) {
  out <- paste0("b_", x)
  out[out == "b_(Intercept)"] <- "b_Intercept"
  out
}

survey_brr_fay_oracle_one_pv <- function(data, pv_col, fay_k = 0.5) {
  rep_cols <- paste0("W_FSTURWT", 1:4)
  design <- survey::svrepdesign(
    weights = stats::as.formula("~ W_FSTUWT"),
    repweights = data[, rep_cols, drop = FALSE],
    data = data,
    type = "Fay",
    rho = fay_k,
    combined.weights = TRUE,
    mse = TRUE
  )
  fit <- survey::svyglm(
    stats::as.formula(paste(pv_col, "~ x")),
    design = design,
    family = stats::gaussian(),
    rescale = FALSE
  )
  beta <- stats::coef(fit)
  U <- stats::vcov(fit)
  attr(U, "means") <- NULL
  names(beta) <- survey_oracle_fe_names(names(beta))
  dimnames(U) <- list(survey_oracle_fe_names(rownames(U)), survey_oracle_fe_names(colnames(U)))
  list(beta = beta, U = U)
}

survey_oracle_rubin_target <- function(pieces) {
  beta_rows <- do.call(rbind, lapply(pieces, `[[`, "beta"))
  fe_names <- names(pieces[[1L]]$beta)
  colnames(beta_rows) <- fe_names
  beta <- colMeans(beta_rows)
  names(beta) <- fe_names
  U_bar <- Reduce(`+`, lapply(pieces, `[[`, "U")) / length(pieces)
  centered <- sweep(t(beta_rows), 1L, beta, FUN = "-")
  B <- centered %*% t(centered) / (length(pieces) - 1L)
  T_MI <- U_bar + (1 + 1 / length(pieces)) * B
  dimnames(U_bar) <- dimnames(B) <- dimnames(T_MI) <- list(fe_names, fe_names)
  list(beta = beta, U_bar = U_bar, B = B, T_MI = T_MI)
}

skip_survey_oracle_tests <- function() {
  skip_if_not(
    identical(Sys.getenv("PVSTACKR_RUN_ORACLE_TESTS"), "true"),
    "survey oracle checks are dev-only; set PVSTACKR_RUN_ORACLE_TESTS=true"
  )
  desc <- suppressWarnings(utils::packageDescription("survey"))
  skip_if(
    length(desc) == 1L && is.na(desc),
    "survey oracle checks require survey"
  )
  skip_if(
    identical(desc$Title, "Sentinel Optional Package") ||
      identical(desc$Version, "0.0.0"),
    "survey oracle checks require the real survey package, not the light-path sentinel"
  )
  skip_if_not_installed("survey")
}

test_that("PISA-style column detection uses natural numeric order", {
  data <- brr_fixture_data()
  data$PV10READ <- data$PV2READ + 1
  data$PV3READ <- data$PV2READ + 0.1
  data$PV4READ <- data$PV2READ + 0.2
  data$PV5READ <- data$PV2READ + 0.3
  data$PV6READ <- data$PV2READ + 0.4
  data$PV7READ <- data$PV2READ + 0.5
  data$PV8READ <- data$PV2READ + 0.6
  data$PV9READ <- data$PV2READ + 0.7

  expect_identical(
    detect_pisa_pv_columns(data, suffix = "READ", expected_M = 10),
    paste0("PV", 1:10, "READ")
  )
  expect_identical(
    detect_pisa_brr_replicate_weights(data, expected_R = 4),
    paste0("W_FSTURWT", 1:4)
  )
  expect_error(detect_pisa_pv_columns(data, suffix = "MATH"), "No plausible-value")
  expect_error(detect_pisa_brr_replicate_weights(data, expected_R = 5), "Expected 5")

  gap <- data[, setdiff(names(data), "W_FSTURWT3")]
  expect_error(detect_pisa_brr_replicate_weights(gap), "contiguous")
})

test_that("PISA-style plausible value detection guides subject suffix use", {
  math <- data.frame(
    PV1MATH = 1:4,
    PV2MATH = 2:5,
    PV3MATH = 3:6
  )

  expect_identical(
    detect_pisa_pv_columns(math, suffix = "MATH", expected_M = 3L),
    paste0("PV", 1:3, "MATH")
  )
  expect_error(
    detect_pisa_pv_columns(math),
    "pv_suffix = \"MATH\""
  )
  expect_error(
    detect_pisa_pv_columns(math),
    "suffix = \"MATH\""
  )

  bare <- data.frame(PV1 = 1:4, PV2 = 2:5)
  expect_warning(
    bare_cols <- detect_pisa_pv_columns(bare, expected_M = 2L),
    NA
  )
  expect_identical(bare_cols, c("PV1", "PV2"))

  mixed <- bare
  mixed$PV1MATH <- 1:4
  mixed$PV2MATH <- 2:5
  expect_warning(
    mixed_cols <- detect_pisa_pv_columns(mixed, expected_M = 2L),
    "subject-suffixed"
  )
  expect_identical(mixed_cols, c("PV1", "PV2"))
})

test_that("pv_brr_target assembles a Rubin/BRR-Fay target from WLS pieces", {
  data <- brr_fixture_data()
  target <- pv_brr_target(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    conf_level = 0.95
  )

  expect_s3_class(target, "pvstackr_brr_target")
  expect_identical(target$target_source, "external_brr_fay_rubin")
  expect_identical(target$engine, "lm")
  expect_identical(target$fe_names, c("b_Intercept", "b_x"))
  expect_equal(target$formula_string, "OUTCOME ~ x")
  expect_equal(target$rhs_string, "x")
  expect_identical(target$M, 2L)
  expect_identical(target$R, 4L)
  expect_equal(target$fay_k, 0.5)
  expect_equal(target$fay_variance_multiplier, 1)
  expect_match(target$design_hash, "^[0-9a-f]{8}$")
  expect_match(target$target_hash, "^[0-9a-f]{8}$")
  expect_true(target$policy$fixed_effects_only)
  expect_identical(target$policy$replicate_weight_role, "external_design_variance_only")
  expect_identical(target$policy$target_repair, "forbidden")
  expect_identical(target$df_method, "classic")
  expect_equal(target$df, target$df_classic, tolerance = 1e-12)
  expect_true(all(is.na(target$df_complete)))
  expect_identical(names(target$df_complete), target$fe_names)
  expect_identical(target$interval_role, "descriptive_classic_rubin")
  expect_false(target$coverage_claim_allowed)
  expect_identical(target$policy$df_method, target$df_method)
  expect_identical(target$policy$interval_role, target$interval_role)
  expect_identical(target$policy$coverage_claim_allowed, target$coverage_claim_allowed)

  pv1 <- hand_one_pv(data, "PV1READ")
  pv2 <- hand_one_pv(data, "PV2READ")
  beta_rows <- rbind(pv1$beta, pv2$beta)
  pooled <- pvstackr:::rubin_pool_matrix(beta_rows, list(pv1$U, pv2$U), orientation = "rows_pv")

  expect_equal(target$beta, pooled$beta, tolerance = 1e-12)
  expect_equal(target$U_bar, pooled$U_bar, tolerance = 1e-12)
  expect_equal(target$B, pooled$B, tolerance = 1e-12)
  expect_equal(target$T_MI, pooled$T_MI, tolerance = 1e-12)
  expect_equal(target$df, pooled$df, tolerance = 1e-12)
})

test_that("pv_brr_target agrees with an external survey Fay replicate-weight oracle", {
  skip_survey_oracle_tests()

  data <- brr_fixture_data()
  target <- pv_brr_target(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    conf_level = 0.95
  )
  expect_equal(
    target$fay_variance_multiplier,
    1 / (target$R * (1 - target$fay_k)^2),
    tolerance = 0
  )

  pieces <- lapply(target$pv_cols, function(pv_col) {
    survey_brr_fay_oracle_one_pv(data, pv_col, fay_k = target$fay_k)
  })
  oracle <- survey_oracle_rubin_target(pieces)
  oracle_tol <- 1e-10

  for (m in seq_along(pieces)) {
    expect_equal(target$per_pv[[m]]$beta, pieces[[m]]$beta, tolerance = oracle_tol)
    expect_equal(target$per_pv[[m]]$U, pieces[[m]]$U, tolerance = oracle_tol)
  }
  expect_equal(target$beta, oracle$beta, tolerance = oracle_tol)
  expect_equal(target$U_bar, oracle$U_bar, tolerance = oracle_tol)
  expect_equal(target$B, oracle$B, tolerance = oracle_tol)
  expect_equal(target$T_MI, oracle$T_MI, tolerance = oracle_tol)
  expect_lt(sqrt(sum((target$T_MI - oracle$T_MI)^2)), oracle_tol)
})

test_that("pv_brr_target per-PV payloads reconstruct replicate and pooled variance pieces", {
  data <- brr_fixture_data()
  target <- pv_brr_target(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )

  for (m in seq_along(target$per_pv)) {
    item <- target$per_pv[[m]]
    expect_identical(item$pv_col, target$pv_cols[[m]])
    expect_equal(dimnames(item$replicate_diff), list(target$fe_names, target$rep_weight_cols))
    expect_equal(
      item$replicate_beta,
      sweep(item$replicate_diff, 1L, item$beta, FUN = "+"),
      tolerance = 1e-12
    )
    expect_equal(
      item$U,
      target$fay_variance_multiplier * tcrossprod(item$replicate_diff),
      tolerance = 1e-12
    )
    expect_equal(item$fay_variance_multiplier, target$fay_variance_multiplier, tolerance = 0)
  }

  beta_rows <- do.call(rbind, lapply(target$per_pv, `[[`, "beta"))
  colnames(beta_rows) <- target$fe_names
  centered <- sweep(t(beta_rows), 1L, target$beta, FUN = "-")
  B <- pvstackr:::pv_symmetrize(centered %*% t(centered) / (target$M - 1L))
  dimnames(B) <- list(target$fe_names, target$fe_names)
  U_bar <- pvstackr:::pv_symmetrize(Reduce(`+`, lapply(target$per_pv, `[[`, "U")) / target$M)
  dimnames(U_bar) <- list(target$fe_names, target$fe_names)

  expect_equal(target$U_bar, U_bar, tolerance = 1e-12)
  expect_equal(target$B, B, tolerance = 1e-12)
  expect_equal(target$T_MI, U_bar + (1 + 1 / target$M) * B, tolerance = 1e-12)
})

test_that("pv_brr_target exposes Barnard-Rubin df metadata only with complete-data df", {
  data <- brr_fixture_data()
  args <- list(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.5
  )

  expect_error(
    do.call(pv_brr_target, c(args, list(df_method = "barnard_rubin"))),
    "df_complete"
  )
  expect_error(
    do.call(pv_brr_target, c(args, list(df_method = "barnard_rubin", df_complete = 0))),
    "positive"
  )
  expect_error(
    do.call(pv_brr_target, c(args, list(
      df_method = "barnard_rubin",
      df_complete = c(b_x = 30, missing = 20)
    ))),
    "aligned"
  )
  expect_error(
    do.call(pv_brr_target, c(args, list(df_method = "classic", df_complete = 20))),
    "only supported"
  )

  target <- do.call(pv_brr_target, c(args, list(
    df_method = "barnard_rubin",
    df_complete = c(b_x = 30, b_Intercept = 20)
  )))

  pv1 <- hand_one_pv(data, "PV1READ")
  pv2 <- hand_one_pv(data, "PV2READ")
  beta_rows <- rbind(pv1$beta, pv2$beta)
  pooled <- pvstackr:::rubin_pool_matrix(
    beta_rows,
    list(pv1$U, pv2$U),
    orientation = "rows_pv",
    df_method = "barnard_rubin",
    df_complete = c(b_Intercept = 20, b_x = 30)
  )

  expect_identical(target$df_method, "barnard_rubin")
  expect_equal(target$df_complete, c(b_Intercept = 20, b_x = 30), tolerance = 0)
  expect_equal(target$df, pooled$df, tolerance = 1e-12)
  expect_equal(target$df_classic, pooled$df_classic, tolerance = 1e-12)
  expect_true(all(target$df <= target$df_classic))
  expect_identical(target$interval_role, "coverage_barnard_rubin")
  expect_true(target$coverage_claim_allowed)
  expect_identical(target$policy$df_method, target$df_method)
  expect_identical(target$policy$interval_role, target$interval_role)
  expect_identical(target$policy$coverage_claim_allowed, target$coverage_claim_allowed)
})

test_that("Fay k controls the BRR design variance multiplier", {
  data <- brr_fixture_data()
  target_fay05 <- pv_brr_target(
    data,
    OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.5
  )
  target_fay0 <- pv_brr_target(
    data,
    OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0
  )

  expect_equal(target_fay05$beta, target_fay0$beta, tolerance = 1e-12)
  expect_equal(target_fay05$U_bar, 4 * target_fay0$U_bar, tolerance = 1e-12)
  expect_equal(target_fay05$B, target_fay0$B, tolerance = 1e-12)
})

test_that("pv_brr_target is equivariant to common plausible-value outcome shifts", {
  data <- brr_fixture_data()
  shifted <- data
  shifted[c("PV1READ", "PV2READ")] <- shifted[c("PV1READ", "PV2READ")] + 25
  args <- list(
    formula = OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )

  target <- do.call(pv_brr_target, c(list(data = data), args))
  shifted_target <- do.call(pv_brr_target, c(list(data = shifted), args))

  expect_equal(
    shifted_target$beta - target$beta,
    c(b_Intercept = 25, b_x = 0),
    tolerance = 1e-12
  )
  expect_equal(shifted_target$U_bar, target$U_bar, tolerance = 1e-12)
  expect_equal(shifted_target$B, target$B, tolerance = 1e-12)
  expect_equal(shifted_target$T_MI, target$T_MI, tolerance = 1e-12)
  expect_equal(shifted_target$df, target$df, tolerance = 1e-12)
  expect_identical(shifted_target$design_hash, target$design_hash)
  expect_false(identical(shifted_target$target_hash, target$target_hash))
})

test_that("pv_brr_target supports the explicit one-plausible-value policy", {
  data <- brr_fixture_data()
  target <- pv_brr_target(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = "PV1READ",
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.5,
    id_cols = "id",
    allow_m1 = TRUE
  )

  zero <- 0 * target$U_bar
  expect_identical(target$M, 1L)
  expect_length(target$per_pv, 1L)
  expect_equal(target$B, zero, tolerance = 0)
  expect_equal(target$U_bar, target$per_pv[[1L]]$U, tolerance = 1e-12)
  expect_equal(target$T_MI, target$U_bar, tolerance = 1e-12)
  expect_equal(target$lambda, c(b_Intercept = 0, b_x = 0), tolerance = 1e-14)
  expect_equal(target$riv, c(b_Intercept = 0, b_x = 0), tolerance = 1e-14)
  expect_true(all(is.infinite(target$df)))
})

test_that("pv_brr_target hashes are deterministic for identical target inputs", {
  data <- brr_fixture_data()
  args <- list(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )

  first <- do.call(pv_brr_target, args)
  second <- do.call(pv_brr_target, args)
  fay0_args <- args
  fay0_args$fay_k <- 0
  fay0 <- do.call(pv_brr_target, fay0_args)

  expect_identical(second$design_hash, first$design_hash)
  expect_identical(second$target_hash, first$target_hash)
  expect_equal(second$beta, first$beta, tolerance = 0)
  expect_equal(second$U_bar, first$U_bar, tolerance = 0)
  expect_equal(second$B, first$B, tolerance = 0)
  expect_equal(second$T_MI, first$T_MI, tolerance = 0)
  expect_equal(fay0$beta, first$beta, tolerance = 1e-12)
  expect_false(identical(fay0$target_hash, first$target_hash))
})

test_that("pv_brr_target validates malformed weights, columns, IDs, formula, and engine", {
  data <- brr_fixture_data()
  rep_cols <- paste0("W_FSTURWT", 1:4)

  expect_error(
    pv_brr_target(data, OUTCOME ~ x, pv_cols = "PV1READ", weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "at least two"
  )
  expect_error(
    pv_brr_target(data, OUTCOME ~ x, pv_cols = c("PV1READ", "missing"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "not found"
  )
  expect_error(
    pv_brr_target(data, OUTCOME ~ x, pv_cols = c("PV1READ", "PV1READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "unique"
  )
  expect_error(
    pv_brr_target(data, OUTCOME ~ x, pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = "W_FSTURWT1"),
    "at least two"
  )
  expect_error(
    pv_brr_target(data, OUTCOME ~ x, pv_cols = c("PV1READ", "PV2READ"), weight_col = "missing", rep_weight_cols = rep_cols),
    "not found"
  )
  expect_error(
    pv_brr_target(data, PV1READ ~ x, pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "OUTCOME"
  )
  expect_error(
    pv_brr_target(data, OUTCOME ~ x + weights(W_FSTUWT), pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "weights"
  )
  data$school <- rep(1:2, length.out = nrow(data))
  expect_error(
    pv_brr_target(data, OUTCOME ~ x + (1 | school), pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "Random-effect"
  )
  expect_error(
    pv_brr_target(data, OUTCOME ~ x, pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols, engine = "lme4"),
    "Only `engine"
  )

  bad <- data
  bad$W_FSTUWT[1] <- 0
  expect_error(
    pv_brr_target(bad, OUTCOME ~ x, pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "strictly positive"
  )
  bad <- data
  bad$W_FSTURWT2[2] <- Inf
  expect_error(
    pv_brr_target(bad, OUTCOME ~ x, pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "finite numeric"
  )
  bad <- data
  bad$id[2] <- bad$id[1]
  expect_error(
    pv_brr_target(bad, OUTCOME ~ x, pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols, id_cols = "id"),
    "unique rows"
  )
})

test_that("pv_brr_target formula guards allow logical OR inside I but reject model bars", {
  data <- brr_fixture_data()
  data$a <- c(-1, 1, -1, 1, -1, 1, -1, 1, -1, 1)
  data$b <- c(1, 1, -1, -1, 1, 1, -1, -1, 1, -1)
  data$school <- rep(1:2, length.out = nrow(data))
  rep_cols <- paste0("W_FSTURWT", 1:4)

  target <- pv_brr_target(
    data,
    OUTCOME ~ x + I((a > 0) | (b < 0)),
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = rep_cols
  )

  expect_s3_class(target, "pvstackr_brr_target")
  expect_true(any(grepl("I", target$fe_names, fixed = TRUE)))

  expect_error(
    pv_brr_target(data, OUTCOME ~ x + (1 | school), pv_cols = c("PV1READ", "PV2READ"), weight_col = "W_FSTUWT", rep_weight_cols = rep_cols),
    "Random-effect"
  )
})
