test_that("Rubin scalar pooling matches hand calculations", {
  out <- pvstackr:::rubin_pool_scalar(q = c(1, 2, 3), u = c(0.25, 0.25, 0.25))

  expect_equal(out$q_bar, 2, tolerance = 1e-14)
  expect_equal(out$u_bar, 0.25, tolerance = 1e-14)
  expect_equal(out$b, 1, tolerance = 1e-14)
  expect_equal(out$total_var, 0.25 + (4 / 3), tolerance = 1e-14)
  expect_equal(out$se, sqrt(19 / 12), tolerance = 1e-14)

  lambda <- (4 / 3) / (19 / 12)
  expect_equal(out$df, 2 / lambda^2, tolerance = 1e-14)
  expect_identical(out$df_method, "classic")
})

test_that("Rubin scalar pooling preserves large between-imputation variance", {
  out <- pvstackr:::rubin_pool_scalar(q = c(1, 3, 8), u = c(0.2, 0.4, 0.6))

  expect_equal(out$q_bar, 4, tolerance = 1e-14)
  expect_equal(out$u_bar, 0.4, tolerance = 1e-14)
  expect_equal(out$b, 13, tolerance = 1e-14)
  expect_equal(out$total_var, 0.4 + (4 / 3) * 13, tolerance = 1e-14)
})

test_that("Rubin pooling reduces to within variance when between variance is zero", {
  scalar <- pvstackr:::rubin_pool_scalar(q = c(2, 2, 2), u = c(0.5, 0.75, 1.0))
  expect_equal(scalar$q_bar, 2, tolerance = 0)
  expect_equal(scalar$u_bar, 0.75, tolerance = 1e-14)
  expect_equal(scalar$b, 0, tolerance = 1e-14)
  expect_equal(scalar$total_var, scalar$u_bar, tolerance = 1e-14)
  expect_equal(scalar$se, sqrt(0.75), tolerance = 1e-14)
  expect_gt(scalar$df, 1e12)

  beta <- rbind(
    PV1 = c(a = 1, b = 3),
    PV2 = c(a = 1, b = 3),
    PV3 = c(a = 1, b = 3)
  )
  U <- list(
    diag(c(0.5, 2.0)),
    diag(c(0.7, 2.2)),
    diag(c(0.9, 2.4))
  )
  out <- pvstackr:::rubin_pool_matrix(beta, U, orientation = "rows_pv")
  U_bar <- diag(c(0.7, 2.2))
  dimnames(U_bar) <- list(c("a", "b"), c("a", "b"))

  expect_equal(out$B, 0 * U_bar, tolerance = 1e-14)
  expect_equal(out$T_MI, U_bar, tolerance = 1e-14)
  expect_equal(out$lambda, c(a = 0, b = 0), tolerance = 1e-14)
  expect_equal(out$riv, c(a = 0, b = 0), tolerance = 1e-14)
  expect_true(all(out$df > 1e12))
})

test_that("Rubin matrix pooling matches independent closed-form oracle", {
  beta <- rbind(
    PV1 = c(theta1 = 1, theta2 = 2),
    PV2 = c(theta1 = 3, theta2 = 5),
    PV3 = c(theta1 = 5, theta2 = 8)
  )
  U <- list(
    matrix(c(1, 0.1, 0.1, 2), 2, 2),
    matrix(c(2, 0.2, 0.2, 4), 2, 2),
    matrix(c(3, 0.3, 0.3, 6), 2, 2)
  )

  out <- pvstackr:::rubin_pool_matrix(beta, U, orientation = "rows_pv")

  U_bar <- matrix(c(2, 0.2, 0.2, 4), 2, 2)
  B <- matrix(c(4, 6, 6, 9), 2, 2)
  T_MI <- U_bar + (4 / 3) * B
  dimnames(U_bar) <- dimnames(B) <- dimnames(T_MI) <- list(colnames(beta), colnames(beta))
  rho <- diag((4 / 3) * B) / diag(T_MI)
  df <- 2 / rho^2

  expect_equal(out$beta, c(theta1 = 3, theta2 = 5), tolerance = 0)
  expect_equal(out$U_bar, U_bar, tolerance = 1e-14)
  expect_equal(out$B, B, tolerance = 1e-14)
  expect_equal(out$T_MI, T_MI, tolerance = 1e-14)
  expect_equal(out$total_var, T_MI, tolerance = 1e-14)
  expect_equal(out$lambda, setNames(rho, colnames(beta)), tolerance = 1e-14)
  expect_equal(out$fmi, setNames(rho, colnames(beta)), tolerance = 1e-14)
  expect_equal(out$rho, setNames(rho, colnames(beta)), tolerance = 1e-14)
  expect_equal(out$riv, setNames(diag((4 / 3) * B) / diag(U_bar), colnames(beta)), tolerance = 1e-14)
  expect_equal(out$df, setNames(df, colnames(beta)), tolerance = 1e-14)
  expect_equal(out$se, setNames(sqrt(diag(T_MI)), colnames(beta)), tolerance = 1e-14)
  expect_false(isTRUE(all.equal(T_MI, U_bar + (2 / 3) * B, tolerance = 1e-14)))
})

test_that("Rubin matrix pooling is invariant to paired plausible-value order", {
  beta <- rbind(
    PV1 = c(theta1 = 1, theta2 = 2),
    PV2 = c(theta1 = 3, theta2 = 5),
    PV3 = c(theta1 = 5, theta2 = 8)
  )
  U <- list(
    matrix(c(1, 0.1, 0.1, 2), 2, 2),
    matrix(c(2, 0.2, 0.2, 4), 2, 2),
    matrix(c(3, 0.3, 0.3, 6), 2, 2)
  )
  perm <- c(3L, 1L, 2L)

  original <- pvstackr:::rubin_pool_matrix(beta, U, orientation = "rows_pv")
  reordered <- pvstackr:::rubin_pool_matrix(beta[perm, , drop = FALSE], U[perm], orientation = "rows_pv")

  for (field in c("beta", "U_bar", "B", "T_MI", "lambda", "riv", "df", "ci_low", "ci_high")) {
    expect_equal(reordered[[field]], original[[field]], tolerance = 1e-14)
  }
})

test_that("Rubin matrix pooling exposes named df metadata", {
  beta <- rbind(
    PV1 = c(a = 1, b = 2),
    PV2 = c(a = 2, b = 4),
    PV3 = c(a = 4, b = 7)
  )
  U <- list(
    diag(c(0.5, 1.0)),
    diag(c(0.6, 1.1)),
    diag(c(0.7, 1.2))
  )

  classic <- pvstackr:::rubin_pool_matrix(beta, U, orientation = "rows_pv")
  adjusted <- pvstackr:::rubin_pool_matrix(
    beta,
    U,
    orientation = "rows_pv",
    df_method = "barnard_rubin",
    df_complete = c(b = 40, a = 20)
  )
  adjusted_oracle <- pvstackr:::rubin_pool_matrix(
    beta,
    U,
    orientation = "rows_pv",
    df_method = "barnard_rubin",
    df_complete = c(a = 20, b = 40)
  )

  expect_identical(classic$df_method, "classic")
  expect_true(all(is.na(classic$df_complete)))
  expect_identical(names(classic$df_complete), colnames(beta))
  expect_identical(adjusted$df_method, "barnard_rubin")
  expect_equal(adjusted$df_complete, c(a = 20, b = 40), tolerance = 0)
  expect_equal(adjusted$df, adjusted_oracle$df, tolerance = 1e-14)
  expect_true(all(adjusted$df <= adjusted$df_classic))
  expect_error(
    pvstackr:::rubin_pool_matrix(
      beta,
      U,
      orientation = "rows_pv",
      df_method = "barnard_rubin",
      df_complete = c(a = 20, missing = 40)
    ),
    "parameter names"
  )
})

test_that("Rubin matrix pooling handles mixed finite and infinite complete-data df", {
  beta <- rbind(
    PV1 = c(a = 1, b = 4),
    PV2 = c(a = 2, b = 6),
    PV3 = c(a = 4, b = 9)
  )
  U <- list(
    diag(c(0.4, 1.0)),
    diag(c(0.5, 1.2)),
    diag(c(0.6, 1.4))
  )

  classic <- pvstackr:::rubin_pool_matrix(beta, U, orientation = "rows_pv")
  adjusted <- pvstackr:::rubin_pool_matrix(
    beta,
    U,
    orientation = "rows_pv",
    df_method = "barnard_rubin",
    df_complete = c(a = Inf, b = 25)
  )

  expect_identical(adjusted$df_method, "barnard_rubin")
  expect_equal(adjusted$df_complete, c(a = Inf, b = 25), tolerance = 0)
  expect_equal(adjusted$df["a"], classic$df["a"], tolerance = 1e-14)
  expect_lt(adjusted$df["b"], classic$df["b"])
  expect_equal(adjusted$df_classic, classic$df_classic, tolerance = 1e-14)
})

test_that("Rubin matrix pooling accepts rows_pv and cols_pv orientations", {
  beta <- rbind(
    PV1 = c(a = 0, b = 2),
    PV2 = c(a = 1, b = 3),
    PV3 = c(a = 2, b = 4)
  )
  U <- replicate(3, diag(c(0.5, 1.5)), simplify = FALSE)

  rows <- pvstackr:::rubin_pool_matrix(beta, U, orientation = "rows_pv")
  cols <- pvstackr:::rubin_pool_matrix(t(beta), U, orientation = "cols_pv")

  expect_equal(cols$beta, rows$beta, tolerance = 1e-14)
  expect_equal(cols$T_MI, rows$T_MI, tolerance = 1e-14)
  expect_identical(cols$orientation, "cols_pv")
})

test_that("Rubin pooling accepts list and array forms of U", {
  beta <- rbind(
    c(a = 1, b = 2),
    c(a = 2, b = 3),
    c(a = 3, b = 4)
  )
  U_list <- replicate(3, matrix(c(1, 0.2, 0.2, 2), 2, 2), simplify = FALSE)
  U_array <- array(unlist(U_list), dim = c(2, 2, 3))

  list_out <- pvstackr:::rubin_pool_matrix(beta, U_list, orientation = "rows_pv")
  array_out <- pvstackr:::rubin_pool_matrix(beta, U_array, orientation = "rows_pv")

  expect_equal(array_out$U_bar, list_out$U_bar, tolerance = 1e-14)
  expect_equal(array_out$T_MI, list_out$T_MI, tolerance = 1e-14)
})

test_that("Barnard-Rubin degrees of freedom are available but not default", {
  classic <- pvstackr:::rubin_pool_scalar(q = c(1, 2, 3), u = c(0.25, 0.25, 0.25))
  adjusted <- pvstackr:::rubin_pool_scalar(
    q = c(1, 2, 3),
    u = c(0.25, 0.25, 0.25),
    df_method = "barnard_rubin",
    df_complete = 20
  )

  rho <- (4 / 3) / (19 / 12)
  nu_old <- 2 / rho^2
  nu_obs <- ((20 + 1) / (20 + 3)) * 20 * (1 - rho)
  hand <- 1 / (1 / nu_old + 1 / nu_obs)

  expect_identical(classic$df_method, "classic")
  expect_identical(adjusted$df_method, "barnard_rubin")
  expect_equal(classic$df, nu_old, tolerance = 1e-14)
  expect_equal(adjusted$df, hand, tolerance = 1e-14)
  expect_lt(adjusted$df, classic$df)
  expect_error(
    pvstackr:::rubin_pool_scalar(c(1, 2, 3), c(0.25, 0.25, 0.25), df_method = "barnard_rubin"),
    "df_complete"
  )
  expect_equal(
    pvstackr:::rubin_pool_scalar(
      c(1, 2, 3),
      c(0.25, 0.25, 0.25),
      df_method = "barnard_rubin",
      df_complete = Inf
    )$df,
    classic$df,
    tolerance = 1e-14
  )
})

test_that("Rubin pooling validates dimensions and M=1 policy", {
  beta <- rbind(
    c(a = 1, b = 2),
    c(a = 2, b = 3),
    c(a = 3, b = 4)
  )
  U <- replicate(3, diag(2), simplify = FALSE)

  expect_error(
    pvstackr:::rubin_pool_matrix(beta[1, , drop = FALSE], U[1], orientation = "rows_pv"),
    "M >= 2"
  )
  one <- pvstackr:::rubin_pool_matrix(
    beta[1, , drop = FALSE],
    U[1],
    orientation = "rows_pv",
    allow_m1 = TRUE
  )
  expect_equal(unname(one$B), matrix(0, 2, 2), tolerance = 0)
  expect_equal(unname(one$T_MI), diag(2), tolerance = 1e-14)
  expect_true(all(is.infinite(one$df)))

  expect_error(
    pvstackr:::rubin_pool_matrix(beta, U[1:2], orientation = "rows_pv"),
    "list length"
  )
  expect_error(
    pvstackr:::rubin_pool_matrix(beta, U, orientation = "rows_pv", conf_level = 1),
    "`conf_level`"
  )
  bad_beta <- beta
  bad_beta[1, 1] <- Inf
  expect_error(
    pvstackr:::rubin_pool_matrix(bad_beta, U, orientation = "rows_pv"),
    "`beta`"
  )
  bad_U <- U
  bad_U[[1]] <- matrix(1, 3, 3)
  expect_error(
    pvstackr:::rubin_pool_matrix(beta, bad_U, orientation = "rows_pv"),
    "p x p"
  )
  bad_U <- U
  bad_U[[1]][1, 1] <- NA_real_
  expect_error(
    pvstackr:::rubin_pool_matrix(beta, bad_U, orientation = "rows_pv"),
    "finite numeric"
  )
  bad_U <- U
  bad_U[[1]][1, 2] <- 0.5
  expect_error(
    pvstackr:::rubin_pool_matrix(beta, bad_U, orientation = "rows_pv"),
    "symmetric"
  )
  bad_U <- U
  bad_U[[1]][1, 1] <- -1
  expect_error(
    pvstackr:::rubin_pool_matrix(beta, bad_U, orientation = "rows_pv"),
    "positive semidefinite"
  )
  named_bad <- U
  rownames(named_bad[[1]]) <- c("a", "b")
  colnames(named_bad[[1]]) <- c("b", "a")
  expect_error(
    pvstackr:::rubin_pool_matrix(beta, named_bad, orientation = "rows_pv"),
    "row and column names"
  )
})

test_that("Rubin scalar pooling validates q and u", {
  expect_error(pvstackr:::rubin_pool_scalar(c(1, 2), c(0.1)), "equal positive length")
  expect_error(pvstackr:::rubin_pool_scalar(c(1, NA), c(0.1, 0.2)), "finite numeric")
  expect_error(pvstackr:::rubin_pool_scalar(c(1, 2), c(0.1, -0.2)), "non-negative")
})
