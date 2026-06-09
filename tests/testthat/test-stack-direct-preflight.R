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

  expect_error(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ x + z,
      target = target
    ),
    "same order"
  )
  expect_error(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ z + x,
      target = stack_direct_target(OUTCOME ~ x + z, data = data)
    ),
    "Fit RHS: `z \\+ x`; target RHS: `x \\+ z`"
  )
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

  expect_error(
    pvstackr:::pv_stack_direct_preflight(
      data = data,
      formula = OUTCOME ~ x + g,
      target = target
    ),
    "factor levels"
  )
})
