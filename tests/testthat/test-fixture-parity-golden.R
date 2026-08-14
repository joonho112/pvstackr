reference_rubin_fixture <- function() {
  beta <- matrix(
    c(
      0.10, 0.15, 0.05, 0.12,
      0.40, 0.35, 0.50, 0.45,
      -0.20, -0.10, -0.25, -0.15
    ),
    nrow = 4L,
    ncol = 3L,
    dimnames = list(NULL, c("b_Intercept", "b_x", "b_z"))
  )
  u <- list(
    matrix(c(0.040, 0.004, 0.002, 0.004, 0.030, 0.003, 0.002, 0.003, 0.020), 3L, 3L),
    matrix(c(0.050, 0.006, 0.001, 0.006, 0.025, 0.004, 0.001, 0.004, 0.018), 3L, 3L),
    matrix(c(0.045, 0.005, 0.003, 0.005, 0.035, 0.002, 0.003, 0.002, 0.022), 3L, 3L),
    matrix(c(0.042, 0.003, 0.002, 0.003, 0.028, 0.003, 0.002, 0.003, 0.019), 3L, 3L)
  )
  names <- c("b_Intercept", "b_x", "b_z")

  list(
    provenance = list(
      source_fixture = "reference Rubin-metrics golden fixture",
      case_id = "rubin_metrics_inline"
    ),
    beta = beta,
    U = u,
    expected = list(
      beta = c(b_Intercept = 0.105, b_x = 0.425, b_z = -0.175),
      U_bar = matrix(c(0.04425, 0.0045, 0.002, 0.0045, 0.0295, 0.003, 0.002, 0.003, 0.01975), 3L, 3L),
      B = matrix(
        c(
          0.00176666666666667, -0.00233333333333333, 0.00266666666666667,
          -0.00233333333333333, 0.00416666666666667, -0.00333333333333333,
          0.00266666666666667, -0.00333333333333333, 0.00416666666666667
        ),
        3L,
        3L,
        dimnames = list(names, names)
      ),
      T_MI = matrix(
        c(
          0.0464583333333333, 0.00158333333333333, 0.00533333333333333,
          0.00158333333333333, 0.0347083333333333, -0.00116666666666667,
          0.00533333333333333, -0.00116666666666667, 0.0249583333333333
        ),
        3L,
        3L,
        dimnames = list(names, names)
      ),
      se = c(b_Intercept = 0.215541952606293, b_x = 0.186301726597832, b_z = 0.157982066492793),
      df = c(b_Intercept = 1327.75898896404, b_x = 133.226688, b_z = 68.889792),
      ci_low = c(b_Intercept = -0.317839912589382, b_x = 0.0565081575992401, b_z = -0.490174594330593),
      ci_high = c(b_Intercept = 0.527839912589382, b_x = 0.79349184240076, b_z = 0.140174594330593),
      M = 4L,
      p = 3L
    )
  )
}

reference_stack_fixture <- function() {
  list(
    provenance = list(
      source_fixture = "reference stacked single-long-fit golden fixture",
      case_id = "stack_single_long_fit"
    ),
    data = data.frame(
      id = 1:5,
      school = c(1, 1, 2, 2, 3),
      x = c(-1, -0.5, 0, 0.5, 1),
      z = c(0, 1, 0, 1, 0),
      W = c(10, 20, 30, 40, 50),
      PV1 = c(1.1, 2.1, 2.9, 4.2, 4.8),
      PV2 = c(0.9, 2.3, 3.2, 4.0, 5.1),
      PV3 = c(1.0, 2.0, 3.1, 4.1, 5.0)
    ),
    draws = matrix(
      c(
        3.03333333333333, 3.03904761904762, 3.04476190476190, 3.05047619047619,
        3.05619047619048, 3.06190476190476, 3.06761904761905, 3.07333333333333,
        2.01000000000000, 2.00142857142857, 1.99285714285714, 1.98428571428571,
        1.97571428571428, 1.96714285714286, 1.95857142857143, 1.95000000000000,
        0.129375790772831, 0.129375790772831, 0.129375790772831, 0.129375790772831,
        0.129375790772831, 0.129375790772831, 0.129375790772831, 0.129375790772831
      ),
      nrow = 8L,
      ncol = 3L,
      dimnames = list(NULL, c("b_Intercept", "b_x", "sigma"))
    ),
    expected_psi = c(b_Intercept = 3.05333333333333, b_x = 1.98)
  )
}

reference_brr_scope_data <- function() {
  data.frame(
    school = factor(rep(1:2, each = 3L)),
    x_within = c(-1.0, -0.5, 0.0, 0.2, 0.7, 1.0),
    x_between = rep(c(0.1, 0.3), each = 3L),
    PV1 = 1:6,
    PV2 = 2:7,
    PV3 = 3:8,
    W_FSTUWT = rep(1, 6),
    W_FSTURWT1 = c(1.2, 0.8, 1.2, 0.8, 1.2, 0.8),
    W_FSTURWT2 = c(0.8, 1.2, 0.8, 1.2, 0.8, 1.2),
    W_FSTURWT3 = c(1.1, 0.9, 1.1, 0.9, 1.1, 0.9),
    W_FSTURWT4 = c(0.9, 1.1, 0.9, 1.1, 0.9, 1.1)
  )
}

test_that("rubin_pool_matrix matches the reference Rubin golden values", {
  fixture <- reference_rubin_fixture()
  pool <- pvstackr:::rubin_pool_matrix(
    beta = fixture$beta,
    U = fixture$U,
    orientation = "rows_pv"
  )
  expected <- fixture$expected

  expect_equal(pool$beta, expected$beta, tolerance = 1e-12)
  expect_equal(unname(pool$U_bar), expected$U_bar, tolerance = 1e-12)
  expect_equal(pool$B, expected$B, tolerance = 1e-12)
  expect_equal(pool$T_MI, expected$T_MI, tolerance = 1e-12)
  expect_equal(pool$se, expected$se, tolerance = 1e-12)
  expect_equal(pool$df, expected$df, tolerance = 1e-8)
  expect_equal(pool$ci_low, expected$ci_low, tolerance = 1e-12)
  expect_equal(pool$ci_high, expected$ci_high, tolerance = 1e-12)
  expect_equal(pool$M, expected$M)
  expect_equal(pool$p, expected$p)
  expect_equal(pool$df_method, "classic")
})

test_that("pv_stack_fit preserves the single-long-fit topology contract", {
  fixture <- reference_stack_fixture()
  record <- new.env(parent = emptyenv())
  fit_function <- function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    record$formula <- formula
    record$data <- data
    record$chains <- chains
    record$iter <- iter
    record$warmup <- warmup
    record$seed <- seed
    list(draws = fixture$draws)
  }

  fit <- pvstackr:::pv_stack_fit(
    data = fixture$data,
    formula = OUTCOME ~ x + z,
    pv_cols = c("PV1", "PV2", "PV3"),
    weight_col = "W",
    control = pv_control(iter = 10L, warmup = 5L, chains = 2L, seed = 20260514L),
    fit_function = fit_function,
    draws_function = function(fit) fit$draws
  )

  expect_equal(nrow(record$data), 15L)
  expect_equal(record$chains, 2L)
  expect_equal(record$iter, 10L)
  expect_equal(record$warmup, 5L)
  expect_equal(record$seed, 20260514L)
  expect_equal(fit$stacked_draws, fixture$draws, tolerance = 1e-12)
  expect_equal(fit$psi_hat_fe, fixture$expected_psi, tolerance = 1e-12)
  expect_equal(fit$weight_summary$topology, "single_long_fit")
  expect_equal(fit$weight_summary$M, 3L)
  expect_equal(fit$weight_summary$n_original, 5L)
  expect_equal(fit$weight_summary$n_long, 15L)
  expect_equal(fit$weight_summary$total_long_weight, 5, tolerance = 1e-12)
  expect_equal(
    fit$weight_summary$per_pv_weight_sum,
    setNames(rep(5 / 3, 3L), c("PV1", "PV2", "PV3")),
    tolerance = 1e-12
  )
  expect_equal(fit$param_map$fe_names, c("b_Intercept", "b_x"))
  expect_equal(fit$param_map$vc_names, "sigma")
  expect_match(fit$weight_summary$long_data_hash, "^[0-9a-f]{8}$")
})

test_that("BRR-Fay mixed-model fixture remains outside current WLS target scope", {
  data <- reference_brr_scope_data()
  expect_error(
    pv_brr_target(
      data = data,
      formula = OUTCOME ~ x_within + x_between + (1 | school),
      pv_cols = c("PV1", "PV2", "PV3"),
      weight_col = "W_FSTUWT",
      rep_weight_cols = paste0("W_FSTURWT", 1:4),
      fay_k = 0.5
    ),
    "Random-effect formula terms"
  )
})

test_that("pisa_tiny shipped artifacts match golden sidecar hashes", {
  bundle <- pisa_tiny_parity_load()
  golden <- pisa_tiny_parity_golden()

  expect_equal(unname(tools::sha256sum(bundle$paths$csv)), unname(golden[["CSV-SHA256"]]))
  expect_equal(unname(tools::sha256sum(bundle$paths$rds)), unname(golden[["RDS-SHA256"]]))
  expect_equal(bundle$cached$data_hash, golden[["Data-Hash"]])
  expect_equal(bundle$cached$row_support_hash, golden[["Row-Support-Hash"]])
  expect_equal(bundle$cached$pv_value_hash, golden[["PV-Value-Hash"]])
  expect_equal(bundle$cached$weight_design_hash, golden[["Weight-Design-Hash"]])
  expect_equal(bundle$cached$design_hash, golden[["Design-Hash"]])
  expect_equal(bundle$cached$target_hash, golden[["Target-Hash"]])
  expect_equal(bundle$cached$long_data_hash, golden[["Long-Data-Hash"]])
  expect_equal(bundle$cached$draws_hash, golden[["Draws-Hash"]])
  expect_equal(bundle$cached$fit_estimates_hash, golden[["Estimates-Hash"]])
  expect_equal(bundle$cached$canonical_fit_hash, golden[["Canonical-Fit-Hash"]])
  expect_equal(bundle$manifest[["Real-PISA-Data"]], "false")
  expect_equal(golden[["Real-PISA-Data"]], "false")
})

test_that("live package rebuild reproduces cached tiny parity signature", {
  skip_if_not_hash_tests()

  live <- pisa_tiny_parity_live_fit()
  cached <- pisa_tiny_parity_load()$cached
  golden <- pisa_tiny_parity_golden()

  expect_equal(live$design$design_hash, golden[["Design-Hash"]])
  expect_equal(live$target$target_hash, golden[["Target-Hash"]])
  expect_equal(live$fit$target$target_hash, cached$fit$target$target_hash)
  expect_equal(live$fit$stack_fit$weight_summary$long_data_hash, golden[["Long-Data-Hash"]])
  expect_equal(pvstackr:::pv_hash_payload(live$fit$stack_fit$stacked_draws), golden[["Draws-Hash"]])
  expect_equal(pvstackr:::pv_hash_payload(live$fit$estimates), golden[["Estimates-Hash"]])
  expect_equal(pisa_tiny_parity_canonical_hash(live$fit, live$target$target_hash), golden[["Canonical-Fit-Hash"]])
})

test_that("cached tiny fit public contract matches compact golden signature", {
  skip_if_not_hash_tests()

  cached <- pisa_tiny_parity_load()$cached
  fit <- cached$fit
  estimates <- round_numeric_columns(get_estimates(fit)[c(
    "term", "estimate", "se", "df", "conf_low", "conf_high", "target_hash"
  )])
  expected_estimates <- data.frame(
    term = c("b_Intercept", "b_x", "b_female"),
    estimate = c(457.894088, 46.883361, 2.143702),
    se = c(1.287312, 0.371793, 3.555031),
    df = c(1.021194, 1.402308, 1.013730),
    conf_low = c(442.318045, 44.414570, -41.606872),
    conf_high = c(473.470132, 49.352151, 45.894276),
    target_hash = rep(
      "sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa",
      3L
    )
  )

  expect_equal(fit$method, "stack_direct")
  expect_equal(fit$status, "ok")
  expect_equal(fit$schema_version, "0.1.0")
  expect_equal(get_target(fit)$target_source, "external_brr_fay_rubin")
  expect_false(fit$control$return_draws)
  expect_null(get_draws(fit))
  expect_equal(estimates, expected_estimates)
})

test_that("fixture parity tests do not replay external paper-repo code", {
  files <- testthat::test_path(c("helper-fixture-parity.R", "test-fixture-parity-golden.R"))
  text <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  local_home_marker <- paste0("/", "Users", "/")
  external_fixture_marker <- paste0("codebase-v2", "/", "inst", "/", "fixtures")
  source_call <- paste0("source", "[(]")

  expect_false(grepl(local_home_marker, text, fixed = TRUE))
  expect_false(grepl(external_fixture_marker, text, fixed = TRUE))
  expect_false(grepl(source_call, text))
})
