skip_if_not_hash_tests <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("PVSTACKR_RUN_HASH_TESTS"), "true"),
    "recomputed content-hash parity is environment-specific; set PVSTACKR_RUN_HASH_TESTS=true (dev/reference environment only)"
  )
}

pisa_tiny_parity_path <- function(...) {
  path <- system.file(..., package = "pvstackr")
  if (nzchar(path)) {
    return(path)
  }
  file.path(getwd(), "inst", ...)
}

pisa_tiny_parity_load <- function() {
  csv_path <- pisa_tiny_parity_path("extdata", "pisa_tiny.csv")
  manifest_path <- pisa_tiny_parity_path("extdata", "pisa_tiny_manifest.dcf")
  rds_path <- pisa_tiny_parity_path("extdata", "examples", "pisa_tiny_stack_direct.rds")

  list(
    paths = list(csv = csv_path, manifest = manifest_path, rds = rds_path),
    data = utils::read.csv(csv_path, stringsAsFactors = FALSE),
    manifest = read.dcf(manifest_path)[1L, ],
    cached = readRDS(rds_path)
  )
}

pisa_tiny_parity_golden <- function() {
  read.dcf(testthat::test_path("fixtures", "pisa_tiny_golden.dcf"))[1L, ]
}

pisa_tiny_parity_fit_function <- function(target) {
  function(formula, data, family, prior, chains, iter, warmup, cores, seed, backend, file, file_refit, ...) {
    beta <- target$beta[target$fe_names]
    se <- sqrt(diag(target$T_MI))[target$fe_names]
    e1 <- c(-0.06, -0.03, 0.00, 0.03, 0.06, -0.04, -0.01, 0.02, 0.05, -0.05, 0.01, 0.02)
    e2 <- c(0.04, -0.02, 0.05, -0.01, 0.03, -0.04, 0.02, -0.05, 0.01, 0.00, 0.03, -0.06)
    e3 <- c(-0.03, 0.02, -0.04, 0.01, 0.05, -0.02, 0.04, -0.01, 0.00, 0.03, -0.05, 0.00)
    list(
      draws = cbind(
        b_Intercept = beta[["b_Intercept"]] + se[["b_Intercept"]] * e1,
        b_x = beta[["b_x"]] + se[["b_x"]] * e2,
        b_female = beta[["b_female"]] + se[["b_female"]] * e3,
        sigma = rep(1, length(e1))
      )
    )
  }
}

pisa_tiny_parity_canonical_hash <- function(fit, target_hash) {
  pvstackr:::pv_hash_payload(list(
    method = fit$method,
    status = fit$status,
    terms = fit$estimates$term,
    estimates = fit$estimates,
    target_hash = target_hash,
    long_data_hash = fit$stack_fit$weight_summary$long_data_hash
  ))
}

pisa_tiny_parity_live_fit <- function() {
  bundle <- pisa_tiny_parity_load()
  design <- pv_design(
    data = bundle$data,
    formula = OUTCOME ~ x + female,
    pv_suffix = "READ",
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID"
  )
  target <- pv_brr_target(
    data = bundle$data,
    formula = OUTCOME ~ x + female,
    pv_cols = design$pv_cols,
    weight_col = design$weight_col,
    rep_weight_cols = design$rep_weight_cols,
    fay_k = design$fay_k,
    id_cols = design$id_cols
  )
  fit <- pv_fit_direct(
    data = bundle$data,
    formula = OUTCOME ~ x + female,
    target = target,
    control = pv_control(
      method = "stack_direct",
      backend = "injected",
      iter = 20L,
      warmup = 10L,
      chains = 2L,
      seed = 20260607L,
      return_draws = FALSE
    ),
    fit_function = pisa_tiny_parity_fit_function(target),
    draws_function = function(fit) fit$draws,
    diagnose_function = test_sampler_diagnose_function(
      chains = 2L,
      post_warmup = 10L
    ),
    cache_dir = tempdir(),
    cache_stem = "pisa-tiny-stack-direct"
  )

  list(data = bundle$data, design = design, target = target, fit = fit)
}

round_numeric_columns <- function(data, digits = 6L) {
  numeric <- vapply(data, is.numeric, logical(1))
  data[numeric] <- lapply(data[numeric], round, digits = digits)
  data
}
