stack_direct_fixture_data <- function() {
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

stack_direct_target <- function(formula = OUTCOME ~ x, data = stack_direct_fixture_data()) {
  pv_brr_target(
    data = data,
    formula = formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = paste0("RW", 1:4),
    fay_k = 0.5,
    id_cols = "id"
  )
}

test_that("stack_direct preflight accepts compatible fixed-effect target", {
  data <- stack_direct_fixture_data()
  target <- stack_direct_target(OUTCOME ~ x)

  out <- pvstackr:::pv_stack_direct_preflight(
    data = data,
    formula = OUTCOME ~ x,
    target = target
  )

  expect_s3_class(out, "pvstackr_stack_direct_preflight")
  expect_identical(out$rhs_string, "x")
  expect_identical(out$fe_names, target$fe_names)
  expect_identical(out$target_hash, target$target_hash)
  expect_identical(out$policy$fixed_effects_only, TRUE)
  expect_identical(names(out$binding_proof), pvstackr:::pv_binding_proof_fields())
  expect_identical(
    out$binding_proof$target_manifest_hash,
    target$binding_manifest$manifest_hash
  )
  expect_identical(
    out$binding_proof$current_manifest_hash,
    target$binding_manifest$manifest_hash
  )
  expect_invisible(pvstackr:::pv_binding_proof_validate(
    out$binding_proof,
    target_manifest = target$binding_manifest
  ))
  expect_false("current_manifest" %in% names(out))
})

test_that("stack_direct preflight rejects target formula mismatch before fitting", {
  data <- stack_direct_fixture_data()
  target <- stack_direct_target(OUTCOME ~ x)
  record <- new.env(parent = emptyenv())
  record$n <- 0L
  fit_function <- function(...) {
    record$n <- record$n + 1L
    list()
  }

  first_error <- tryCatch(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ x + z,
      target = target
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(first_error$code, "PV_BIND_E030")
  expect_true("PV_BIND_E040" %in% first_error$all_codes)
  second_error <- tryCatch(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ z + x,
      target = stack_direct_target(OUTCOME ~ x + z, data = data)
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(second_error$code, "PV_BIND_E040")
  expect_identical(record$n, 0L)
  expect_true(is.function(fit_function))
})

test_that("stack_direct preflight rejects group terms for BRR target compatibility", {
  data <- stack_direct_fixture_data()
  target <- stack_direct_target(OUTCOME ~ x)

  expect_error(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ x + (1 | school),
      target = target
    ),
    "Random-effect/group terms"
  )
  expect_error(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ x + (1 || school),
      target = target
    ),
    "Random-effect/group terms"
  )
})

test_that("stack_direct preflight allows logical OR when target and formula agree", {
  data <- stack_direct_fixture_data()
  formula <- OUTCOME ~ x + I((z > 0) | (x > 0))
  target <- stack_direct_target(formula)

  out <- pvstackr:::pv_stack_direct_preflight(
    data = data,
    formula = formula,
    target = target
  )

  expect_identical(out$rhs_string, target$rhs_string)
  expect_identical(out$fe_names, target$fe_names)
})

test_that("stack_direct preflight rejects stale targets without formula metadata", {
  data <- stack_direct_fixture_data()
  target <- stack_direct_target(OUTCOME ~ x)
  target$rhs_string <- NULL

  expect_error(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ x,
      target = target
    ),
    "missing required field"
  )
})

test_that("stack_direct preflight rejects fixed-effect name mismatches", {
  data <- stack_direct_fixture_data()
  data$g <- factor(rep(c("A", "B"), each = 5))
  target <- stack_direct_target(OUTCOME ~ x + g, data = data)
  data$g <- factor(rep("A", nrow(data)), levels = "A")

  error <- tryCatch(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ x + g,
      target = target
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(error$code, "PV_BIND_E081")
})

test_that("stack_direct preflight resolves a stable stateful RHS once and returns its exact bundle", {
  data <- stack_direct_fixture_data()
  state <- new.env(parent = baseenv())
  state$count <- 0L
  state$scale <- 1
  state$transform_x <- function(x) {
    state$count <- state$count + 1L
    state$scale * x
  }
  formula <- stats::as.formula("OUTCOME ~ transform_x(x)", env = state)
  target <- stack_direct_target(formula, data = data)

  state$count <- 0L
  expected_bundle <- pvstackr:::pv_binding_resolve_model_bundle(data, formula)
  expect_identical(state$count, 1L)

  state$count <- 0L
  resolved <- pvstackr:::pv_stack_direct_preflight(
    data,
    formula,
    target,
    return_model_bundle = TRUE
  )
  expect_identical(names(resolved), c("preflight", "model_bundle"))
  expect_s3_class(resolved$preflight, "pvstackr_stack_direct_preflight")
  expect_identical(resolved$model_bundle, expected_bundle)
  expect_identical(state$count, 1L)

  state$count <- 0L
  state$scale <- 2
  drift <- tryCatch(
    pvstackr:::pv_stack_direct_preflight(data, formula, target),
    pvstackr_binding_error = identity
  )
  expect_identical(drift$code, "PV_BIND_E042")
  expect_identical(state$count, 1L)
})

test_that("stack_direct preflight normalizes Gaussian family and rejects incompatible family", {
  data <- stack_direct_fixture_data()
  target <- stack_direct_target(OUTCOME ~ x, data = data)
  expect_s3_class(
    pvstackr:::pv_stack_direct_preflight(
      data,
      OUTCOME ~ x,
      target,
      family = stats::gaussian(link = "identity")
    ),
    "pvstackr_stack_direct_preflight"
  )
  error <- tryCatch(
    pvstackr:::pv_stack_direct_preflight(
      data,
      OUTCOME ~ x,
      target,
      family = stats::binomial()
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(error$code, "PV_BIND_E060")
})
