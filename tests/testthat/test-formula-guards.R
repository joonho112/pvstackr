test_that("formula guards detect embedded weights calls without text false positives", {
  expect_true(pvstackr:::pv_formula_has_weights_call(OUTCOME ~ x + weights(W)))
  expect_true(pvstackr:::pv_formula_has_weights_call(OUTCOME | weights(W) ~ x))
  expect_true(pvstackr:::pv_formula_has_weights_call(OUTCOME ~ x + I(weights(W) > 0)))
  expect_true(pvstackr:::pv_formula_has_weights_call(OUTCOME ~ x + foo(weights(W))))
  expect_false(pvstackr:::pv_formula_has_weights_call(OUTCOME ~ x + weights))
  expect_false(pvstackr:::pv_formula_has_weights_call(OUTCOME ~ x + I("weights(W)")))
})

test_that("formula guards preserve OUTCOME LHS contract before RHS inspection", {
  expect_error(pvstackr:::pv_stack_formula_rhs(OUTCOME | weights(W) ~ x), "OUTCOME")
  expect_error(pvstackr:::pv_validate_design_formula(OUTCOME | weights(W) ~ x), "OUTCOME")
  expect_error(pvstackr:::pv_formula_for_pv(OUTCOME | weights(W) ~ x, "PV1"), "OUTCOME")
  expect_error(pvstackr:::pv_formula_for_pv(OUTCOME | school ~ x, "PV1"), "OUTCOME")
})

test_that("formula guards distinguish model bars from fixed-effect logical OR", {
  expect_true(pvstackr:::pv_formula_has_random_effect_bar(OUTCOME ~ x + (1 | school)))
  expect_true(pvstackr:::pv_formula_has_random_effect_bar(OUTCOME ~ x + (1 || school)))
  expect_true(pvstackr:::pv_formula_has_random_effect_bar(OUTCOME ~ x + ((1 | school))))
  expect_false(pvstackr:::pv_formula_has_random_effect_bar(OUTCOME ~ x + I((a > 0) | (b < 0))))
  expect_false(pvstackr:::pv_formula_has_random_effect_bar(OUTCOME ~ x + ifelse(a | b, 1, 0)))

  expect_silent(pvstackr:::pv_formula_for_pv(OUTCOME ~ x + I((a > 0) | (b < 0)), "PV1"))
  expect_silent(pvstackr:::pv_formula_for_pv(OUTCOME ~ x + ifelse(a | b, 1, 0), "PV1"))
  expect_error(pvstackr:::pv_formula_for_pv(OUTCOME ~ x + (1 | school), "PV1"), "Random-effect")
  expect_error(pvstackr:::pv_formula_for_pv(OUTCOME ~ x + (1 || school), "PV1"), "Random-effect")
})
