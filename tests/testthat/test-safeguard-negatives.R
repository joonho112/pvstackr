safeguard_direct_bundle <- function() {
  source <- pisa_tiny_parity_load()
  design <- pv_design(
    data = source$data,
    formula = OUTCOME ~ x + female,
    pv_suffix = "READ",
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID"
  )
  target <- pv_brr_target(
    data = source$data,
    formula = design$formula,
    pv_cols = design$pv_cols,
    weight_col = design$weight_col,
    rep_weight_cols = design$rep_weight_cols,
    fay_k = design$fay_k,
    id_cols = design$id_cols,
    df_method = "barnard_rubin",
    df_complete = 80
  )
  list(data = source$data, design = design, target = target)
}

safeguard_direct_draws <- function(target, center_shift_se = 0,
                                   near_singular = FALSE) {
  n <- 12L
  z <- cbind(
    c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, -2.5, 2.5, 0),
    c(1, -1, 2, -2, 0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 0.25, -0.25),
    c(-1.8, 0.2, 1.1, -0.7, 2.2, -2.1, 0.9, -0.3, 1.7, -1.4, 0.5, -0.3)
  )
  z <- sweep(z, 2L, colMeans(z), FUN = "-")
  whitened <- z %*% solve(chol(stats::cov(z)))
  raw_covariance <- diag(diag(target$T_MI) * c(0.7, 1.2, 0.9))
  fe <- whitened %*% chol(raw_covariance)
  shifted_center <- target$beta + center_shift_se * target$se
  fe <- sweep(fe, 2L, shifted_center, FUN = "+")
  colnames(fe) <- target$fe_names
  if (isTRUE(near_singular)) {
    epsilon <- seq(-1, 1, length.out = n)
    epsilon <- epsilon - mean(epsilon)
    fe[, "b_female"] <-
      shifted_center[["b_female"]] + 1e-12 * epsilon
  }
  cbind(fe, sigma = rep(1, n))
}

safeguard_direct_fit <- function(center_shift_se = 0,
                                 near_singular = FALSE,
  prior = NULL) {
  bundle <- safeguard_direct_bundle()
  callbacks <- new.env(parent = emptyenv())
  callbacks$count <- 0L
  fit <- pv_fit(
    data = bundle$data,
    formula = bundle$design$formula,
    target = bundle$target,
    method = "stack_direct",
    control = pv_control(
      method = "stack_direct",
      backend = "injected",
      chains = 2L,
      iter = 12L,
      warmup = 6L,
      cores = 1L,
      seed = 20260713L,
      return_draws = TRUE,
      keep_data = FALSE,
      keep_backend_fit = FALSE
    ),
    prior = prior,
    fit_function = function(...) {
      callbacks$count <- callbacks$count + 1L
      list(draws = safeguard_direct_draws(
        bundle$target,
        center_shift_se = center_shift_se,
        near_singular = near_singular
      ))
    },
    draws_function = function(fit, ...) fit$draws,
    diagnose_function = test_sampler_diagnose_function(
      chains = 2L,
      post_warmup = 6L
    ),
    cache_dir = tempdir(),
    cache_stem = "safeguard-negative"
  )
  list(bundle = bundle, fit = fit, callbacks = callbacks)
}

safeguard_psis_draws <- function() {
  cbind(
    b_Intercept = seq(0.8, 1.9, length.out = 12L),
    b_x = seq(-0.4, 0.7, length.out = 12L),
    sigma = seq(0.8, 1.2, length.out = 12L)
  )
}

safeguard_psis_weights <- function() {
  cbind(
    PV1 = rep(1, 12L),
    PV2 = rep(c(1, 2), 6L),
    PV3 = rep(c(2, 1, 1), 4L)
  )
}

safeguard_psis_data <- function() {
  data.frame(
    id = seq_len(5L),
    x = seq(-1, 1, length.out = 5L),
    PV1 = 1:5,
    PV2 = 1.1 + 1:5,
    PV3 = 0.9 + 1:5,
    W = 1
  )
}

safeguard_psis_fit <- function(pareto_k, fallback = "warn") {
  suppressWarnings(pv_fit(
    data = safeguard_psis_data(),
    formula = OUTCOME ~ x,
    method = "stack_psis",
    pv_cols = c("PV1", "PV2", "PV3"),
    weight_col = "W",
    stacked_draws = safeguard_psis_draws(),
    psis_weights = safeguard_psis_weights(),
    pareto_k = pareto_k,
    psis_producer = "safeguard_fixture",
    psis_producer_version = "1.0.0",
    fallback = fallback,
    control = pv_control(
      method = "stack_psis",
      backend = "injected",
      chains = 2L,
      iter = 12L,
      warmup = 6L,
      return_draws = TRUE,
      keep_data = FALSE,
      keep_backend_fit = FALSE
    )
  ))
}

test_that("prior safeguards distinguish warning and pre-backend rejection", {
  invariant_prior <- data.frame(
    prior = "normal(0, 1)",
    class = "sigma",
    coef = "",
    stringsAsFactors = FALSE
  )
  warning_result <- safeguard_direct_fit(prior = invariant_prior)
  warning_fit <- warning_result$fit
  expect_identical(warning_fit$status, "warning")
  expect_identical(warning_fit$reason_codes, "explicit_prior_warning")
  expect_true(any(grepl("Explicit priors", warning_fit$warnings, fixed = TRUE)))
  expect_gt(nrow(get_estimates(warning_fit)), 0L)
  expect_identical(warning_result$callbacks$count, 1L)
  expect_invisible(pvstackr:::validate_pvstackr_fit(warning_fit))

  bundle <- safeguard_direct_bundle()
  # A population-level `b` prior is not listed here: it is expanded onto the
  # slope columns instead of refused (see the prior-scope test below). What
  # stays invalid is anything whose scope expansion would not be exact.
  invalid_priors <- list(
    coefficient_b = data.frame(
      prior = "normal(0, 1)",
      class = "b",
      coef = "x",
      stringsAsFactors = FALSE
    )
  )
  for (name in names(invalid_priors)) {
    callbacks <- new.env(parent = emptyenv())
    callbacks$count <- 0L
    error <- tryCatch(
      pv_fit(
        data = bundle$data,
        formula = bundle$design$formula,
        target = bundle$target,
        method = "stack_direct",
        control = pv_control(
          method = "stack_direct",
          backend = "injected",
          chains = 2L,
          iter = 12L,
          warmup = 6L
        ),
        prior = invalid_priors[[name]],
        fit_function = function(...) {
          callbacks$count <- callbacks$count + 1L
          list()
        },
        draws_function = function(fit, ...) fit$draws,
        diagnose_function = test_sampler_diagnose_function(
          chains = 2L,
          post_warmup = 6L
        )
      ),
      error = identity
    )
    expect_s3_class(error, "simpleError")
    expect_match(
      conditionMessage(error),
      "cannot be preserved exactly after stack_direct model-matrix materialization",
      fixed = TRUE,
      info = name
    )
    expect_identical(callbacks$count, 0L, info = name)
  }

  center_callbacks <- new.env(parent = emptyenv())
  center_callbacks$count <- 0L
  center_error <- tryCatch(
    pv_fit(
      data = bundle$data,
      formula = bundle$design$formula,
      target = bundle$target,
      method = "stack_direct",
      control = pv_control(
        method = "stack_direct",
        backend = "injected",
        center = "posterior",
        chains = 2L,
        iter = 12L,
        warmup = 6L
      ),
      fit_function = function(...) {
        center_callbacks$count <- center_callbacks$count + 1L
        list()
      },
      draws_function = function(fit, ...) fit$draws,
      diagnose_function = test_sampler_diagnose_function(
        chains = 2L,
        post_warmup = 6L
      )
    ),
    error = identity
  )
  expect_s3_class(center_error, "simpleError")
  expect_match(
    conditionMessage(center_error),
    "control$center = \"target\"",
    fixed = TRUE
  )
  expect_identical(center_callbacks$count, 0L)
})

test_that("center and conditioning gates preserve warning or redact on block", {
  yellow <- safeguard_direct_fit(center_shift_se = 0.02)$fit
  expect_identical(yellow$status, "warning")
  expect_identical(yellow$reason_codes, "center_separation_yellow")
  expect_identical(yellow$ccc$diagnostics$center_status, "warning")
  expect_gt(nrow(yellow$estimates), 0L)
  expect_false(is.null(yellow$draws))

  red <- safeguard_direct_fit(center_shift_se = 0.08)$fit
  expect_identical(red$status, "blocked")
  expect_identical(red$reason_codes, "center_separation_red")
  expect_identical(red$diagnostics$ccc$center$status, "blocked")
  expect_identical(red$diagnostics$ccc$conditioning$status, "ok")
  expect_true(any(grepl("block threshold", red$warnings, fixed = TRUE)))
  expect_null(red$design)
  expect_null(red$stack_fit)
  expect_null(red$ccc)
  expect_null(red$draws)
  expect_identical(nrow(red$estimates), 0L)
  expect_s3_class(red$target, "pvstackr_brr_target")
  expect_true(isTRUE(red$provenance$independent_target_retained))
  expect_true(all(vapply(
    red$control[c("return_draws", "keep_data", "keep_backend_fit", "keep_log_lik")],
    identical,
    logical(1),
    FALSE
  )))
  expect_setequal(
    names(red$diagnostics),
    c("preflight", "sampler", "sampler_gate", "ccc", "redaction")
  )
  expect_identical(
    red$diagnostics$redaction$policy,
    "generic_blocked_fail_closed"
  )
  expect_identical(red$diagnostics$redaction$source, "ccc_reportability_gate")
  expect_invisible(pvstackr:::validate_pvstackr_fit(red))

  conditioning <- safeguard_direct_fit(near_singular = TRUE)$fit
  expect_identical(conditioning$status, "blocked")
  expect_identical(conditioning$reason_codes, "ccc_conditioning_red")
  expect_identical(conditioning$diagnostics$ccc$center$status, "ok")
  expect_identical(
    conditioning$diagnostics$ccc$conditioning$status,
    "blocked"
  )
  expect_gt(conditioning$diagnostics$ccc$conditioning$kappa_A, 1e8)
  expect_true(any(grepl(
    "reportable stack_direct estimates are blocked",
    conditioning$warnings,
    fixed = TRUE
  )))
  expect_null(conditioning$design)
  expect_null(conditioning$stack_fit)
  expect_null(conditioning$ccc)
  expect_null(conditioning$draws)
  expect_identical(nrow(conditioning$estimates), 0L)
  expect_invisible(pvstackr:::validate_pvstackr_fit(conditioning))
})

test_that("group binding tamper and family mismatches fail before callbacks", {
  bundle <- safeguard_direct_bundle()
  grouped_data <- transform(
    bundle$data,
    school = factor(rep(c("A", "B", "C"), length.out = nrow(bundle$data)))
  )
  changed_predictor <- bundle$data
  changed_predictor$x[[1L]] <- changed_predictor$x[[1L]] + 0.1
  tampered_target <- bundle$target
  tampered_target$target_content$derived$beta[[1L]] <-
    tampered_target$target_content$derived$beta[[1L]] + 0.1
  cases <- list(
    group = list(
      data = grouped_data,
      formula = OUTCOME ~ x + female + (1 | school),
      target = bundle$target,
      family = NULL,
      class = "simpleError",
      code = NULL,
      message = "Random-effect/group terms"
    ),
    predictor_tamper = list(
      data = changed_predictor,
      formula = bundle$design$formula,
      target = bundle$target,
      family = NULL,
      class = "pvstackr_binding_predictor_values_mismatch",
      code = "PV_BIND_E031",
      message = "PREDICTOR_VALUES_MISMATCH"
    ),
    target_content_tamper = list(
      data = bundle$data,
      formula = bundle$design$formula,
      target = tampered_target,
      family = NULL,
      class = "pvstackr_binding_target_content_hash_stale",
      code = "PV_BIND_E090",
      message = "TARGET_CONTENT_HASH_STALE"
    ),
    family = list(
      data = bundle$data,
      formula = bundle$design$formula,
      target = bundle$target,
      family = stats::binomial(),
      class = "pvstackr_binding_family_mismatch",
      code = "PV_BIND_E060",
      message = "FAMILY_MISMATCH"
    ),
    link = list(
      data = bundle$data,
      formula = bundle$design$formula,
      target = bundle$target,
      family = stats::gaussian(link = "log"),
      class = "pvstackr_binding_link_mismatch",
      code = "PV_BIND_E061",
      message = "LINK_MISMATCH"
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    callbacks <- new.env(parent = emptyenv())
    callbacks$count <- 0L
    error <- tryCatch(
      pv_fit(
        data = case$data,
        formula = case$formula,
        target = case$target,
        method = "stack_direct",
        family = case$family,
        control = pv_control(
          method = "stack_direct",
          backend = "injected",
          chains = 2L,
          iter = 12L,
          warmup = 6L
        ),
        fit_function = function(...) {
          callbacks$count <- callbacks$count + 1L
          list()
        },
        draws_function = function(fit, ...) fit$draws,
        diagnose_function = test_sampler_diagnose_function(
          chains = 2L,
          post_warmup = 6L
        )
      ),
      error = identity
    )
    expect_s3_class(error, case$class)
    expect_match(conditionMessage(error), case$message, fixed = TRUE, info = name)
    if (!is.null(case$code)) {
      expect_s3_class(error, "pvstackr_binding_error")
      expect_identical(error$code, case$code, info = name)
    }
    expect_identical(callbacks$count, 0L, info = name)
  }
})

test_that("PSIS boundary and unevaluated diagnostics are immutable blocks", {
  boundary <- safeguard_psis_fit(c(PV1 = 0.2, PV2 = 0.7, PV3 = 0.3))
  expect_identical(boundary$status, "blocked")
  expect_identical(boundary$reason_codes, "psis_k_too_high")
  expect_identical(boundary$diagnostics$psis$status, "failed")
  expect_identical(boundary$diagnostics$psis$bad_pv_cols, "PV2")
  expect_identical(boundary$diagnostics$psis$fallback_requested, "warn")
  expect_identical(boundary$diagnostics$psis$fallback_effective, "block")
  expect_true(any(grepl("met or exceeded", boundary$warnings, fixed = TRUE)))
  expect_setequal(names(boundary$diagnostics), c("psis", "redaction"))
  expect_identical(boundary$diagnostics$redaction$status, "withheld")
  expect_identical(
    boundary$diagnostics$redaction$policy,
    "immutable_psis_fail_closed"
  )
  expect_false(any(c("pooling", "weighted") %in% names(boundary$diagnostics)))
  expect_identical(nrow(boundary$estimates), 0L)
  expect_null(boundary$design)
  expect_null(boundary$stack_fit)
  expect_null(boundary$ccc)
  expect_null(boundary$draws)
  expect_true(all(vapply(
    boundary$control[c("return_draws", "keep_data", "keep_backend_fit", "keep_log_lik")],
    identical,
    logical(1),
    FALSE
  )))
  expect_invisible(pvstackr:::validate_pvstackr_fit(boundary))

  not_evaluated <- safeguard_psis_fit(c(PV1 = 0.2, PV2 = NA_real_, PV3 = Inf))
  expect_identical(not_evaluated$status, "blocked")
  expect_identical(not_evaluated$reason_codes, "psis_k_not_evaluated")
  expect_identical(not_evaluated$diagnostics$psis$status, "not_evaluated")
  expect_true(any(grepl("not fully evaluated", not_evaluated$warnings, fixed = TRUE)))
  expect_identical(
    not_evaluated$diagnostics$psis$bad_pv_cols,
    c("PV2", "PV3")
  )
  expect_identical(nrow(not_evaluated$estimates), 0L)
  expect_setequal(names(not_evaluated$diagnostics), c("psis", "redaction"))
  expect_invisible(pvstackr:::validate_pvstackr_fit(not_evaluated))
})

test_that("a population-level prior reaches the backend scoped to the slopes", {
  bundle <- safeguard_direct_bundle()
  seen <- new.env(parent = emptyenv())
  global_b <- data.frame(
    prior = "normal(0, 25)", class = "b", coef = "",
    stringsAsFactors = FALSE
  )

  fit <- pv_fit(
    data = bundle$data,
    formula = bundle$design$formula,
    target = bundle$target,
    method = "stack_direct",
    control = pv_control(
      method = "stack_direct", backend = "injected", chains = 2L,
      iter = 12L, warmup = 6L, cores = 1L, seed = 20260713L,
      return_draws = TRUE
    ),
    prior = global_b,
    fit_function = function(formula, data, family, prior, ...) {
      seen$prior <- prior
      seen$formula <- paste(deparse(formula), collapse = " ")
      list(draws = safeguard_direct_draws(bundle$target, center_shift_se = 0))
    },
    draws_function = function(fit, ...) fit$draws,
    diagnose_function = test_sampler_diagnose_function(chains = 2L, post_warmup = 6L),
    cache_dir = tempdir(),
    cache_stem = "safeguard-prior-scope"
  )

  # The materialized design fits `0 + ...`, so the intercept is an ordinary
  # column. The prior must arrive attached to the slope columns only, which is
  # what the original formula's population-level prior meant.
  intercept_column <- "pvstackrMM001"
  expect_match(seen$formula, paste0("0 \\+ ", intercept_column))
  expect_true(all(seen$prior$class == "b"))
  expect_false(intercept_column %in% seen$prior$coef)
  expect_setequal(seen$prior$coef, c("pvstackrMM002", "pvstackrMM003"))
  expect_true(all(seen$prior$prior == "normal(0, 25)"))

  # The prior policy still flags the fit rather than passing it silently.
  expect_identical(fit$status, "warning")
  expect_identical(fit$reason_codes, "explicit_prior_warning")
})
