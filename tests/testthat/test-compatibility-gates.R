compatibility_error_code <- function(expr) {
  error <- tryCatch(
    force(expr),
    pvstackr_binding_error = identity
  )
  expect_s3_class(error, "pvstackr_binding_error")
  error$code
}

test_that("family compatibility normalizes supported inputs to one safe tuple", {
  expected <- pvstackr:::pv_binding_family_link_projection(
    "gaussian",
    "identity",
    "real",
    "estimated"
  )
  brms_like <- structure(
    list(
      family = "gaussian",
      link = "identity",
      backend_fit = "PRIVATE_BACKEND_PAYLOAD",
      raw_data = data.frame(secret = 1)
    ),
    class = c("brmsfamily", "list")
  )
  inputs <- list(
    NULL,
    "gaussian",
    stats::setNames("gaussian", "PRIVATE_SCALAR_ATTRIBUTE"),
    stats::gaussian(link = "identity"),
    stats::gaussian(identity),
    brms_like
  )

  for (input in inputs) {
    observed <- pvstackr:::pv_compatibility_family_link(input)
    expect_identical(observed, expected)
    expect_identical(
      names(observed),
      c("family_id", "link_id", "response_support_id", "dispersion_role")
    )
    serialized <- rawToChar(serialize(observed, NULL, ascii = TRUE))
    expect_false(grepl("PRIVATE_BACKEND_PAYLOAD", serialized, fixed = TRUE))
    expect_false(grepl("raw_data", serialized, fixed = TRUE))
  }
})

test_that("canonical backend family does not retain caller family state", {
  secret <- "CALLER_PRIVATE_FAMILY_SECRET"
  forged <- stats::gaussian(link = "identity")
  forged$secret_payload <- secret
  forged$linkfun <- local({
    private_state <- secret
    function(mu) rep(-999, length(mu))
  })

  canonical <- pvstackr:::pv_compatibility_stack_direct_family(forged)
  expect_s3_class(canonical, "family")
  expect_identical(canonical$family, "gaussian")
  expect_identical(canonical$link, "identity")
  expect_identical(canonical$linkfun(2), 2)
  serialized <- rawToChar(serialize(canonical, NULL, ascii = TRUE))
  expect_false(grepl(secret, serialized, fixed = TRUE))
  expect_false(grepl("secret_payload", serialized, fixed = TRUE))
})

test_that("family compatibility uses E060 and E061 for semantic failures", {
  non_gaussian <- list(
    "poisson",
    stats::poisson(link = "identity"),
    list(family = "binomial", link = "identity")
  )
  for (input in non_gaussian) {
    expect_identical(
      compatibility_error_code(
        pvstackr:::pv_compatibility_family_link(input)
      ),
      "PV_BIND_E060"
    )
  }

  non_identity <- list(
    stats::gaussian(link = "log"),
    list(family = "gaussian", link = "inverse")
  )
  for (input in non_identity) {
    expect_identical(
      compatibility_error_code(
        pvstackr:::pv_compatibility_family_link(input)
      ),
      "PV_BIND_E061"
    )
  }
})

test_that("family compatibility rejects malformed metadata with E081", {
  duplicate_names <- list("gaussian", "poisson", "identity")
  names(duplicate_names) <- c("family", "family", "link")
  malformed <- list(
    c("gaussian", "poisson"),
    NA_character_,
    1,
    list(family = "gaussian"),
    list(family = c("gaussian", "poisson"), link = "identity"),
    list(family = "gaussian", link = NA_character_),
    duplicate_names,
    list("gaussian", "identity")
  )
  for (input in malformed) {
    expect_identical(
      compatibility_error_code(
        pvstackr:::pv_compatibility_family_link(input)
      ),
      "PV_BIND_E081"
    )
  }
})

test_that("BRR-Fay estimand compatibility delegates exact classic and coverage tuples", {
  fe_names <- c("b_Intercept", "b_x")
  cases <- list(
    list(
      interval_role = "descriptive_classic_rubin",
      coverage_claim_allowed = FALSE
    ),
    list(
      interval_role = "coverage_barnard_rubin",
      coverage_claim_allowed = TRUE
    )
  )
  for (case in cases) {
    contrast <- pvstackr:::pv_binding_estimand_contrast_projection(fe_names)
    expected <- pvstackr:::pv_binding_estimand_projection(
      fe_names = fe_names,
      estimand_contrast = contrast,
      interval_role = case$interval_role,
      coverage_claim_allowed = case$coverage_claim_allowed
    )
    observed <- pvstackr:::pv_compatibility_brr_fay_estimand(
      fe_names = fe_names,
      interval_role = case$interval_role,
      coverage_claim_allowed = case$coverage_claim_allowed
    )
    expect_identical(observed, expected)
    expect_identical(
      names(observed),
      c(
        "estimand_id", "target_source", "target_engine_id",
        "parameter_scope", "fe_names", "interval_role",
        "coverage_claim_allowed", "estimand_contrast_hash"
      )
    )
  }

  named_fe <- stats::setNames(
    fe_names,
    c("PRIVATE_FE_ATTRIBUTE_1", "PRIVATE_FE_ATTRIBUTE_2")
  )
  normalized <- pvstackr:::pv_compatibility_brr_fay_estimand(
    fe_names = named_fe,
    interval_role = "descriptive_classic_rubin",
    coverage_claim_allowed = FALSE
  )
  expect_null(names(normalized$fe_names))
  serialized <- rawToChar(serialize(normalized, NULL, ascii = TRUE))
  expect_false(grepl("PRIVATE_FE_ATTRIBUTE", serialized, fixed = TRUE))
})

test_that("BRR-Fay estimand compatibility rejects semantic mismatches with E070", {
  base_args <- list(
    fe_names = c("b_Intercept", "b_x"),
    interval_role = "descriptive_classic_rubin",
    coverage_claim_allowed = FALSE
  )
  incompatible <- list(
    c(base_args, list(estimand_id = "other_estimand")),
    c(base_args, list(target_source = "stacked_posterior")),
    c(base_args, list(target_engine_id = "lm")),
    c(base_args, list(parameter_scope = "variance_component")),
    within(base_args, interval_role <- "coverage_barnard_rubin"),
    within(base_args, coverage_claim_allowed <- TRUE)
  )
  for (args in incompatible) {
    expect_identical(
      compatibility_error_code(do.call(
        pvstackr:::pv_compatibility_brr_fay_estimand,
        args
      )),
      "PV_BIND_E070"
    )
  }
})

test_that("BRR-Fay estimand compatibility rejects malformed fields with E081", {
  base_args <- list(
    fe_names = c("b_Intercept", "b_x"),
    interval_role = "descriptive_classic_rubin",
    coverage_claim_allowed = FALSE
  )
  malformed <- list(
    c(base_args, list(target_source = NA_character_)),
    c(base_args, list(target_engine_id = c("lm", "other"))),
    c(base_args, list(parameter_scope = NA_character_)),
    within(base_args, interval_role <- NA_character_),
    within(base_args, coverage_claim_allowed <- NA),
    within(base_args, coverage_claim_allowed <- c(TRUE, FALSE)),
    within(base_args, fe_names <- c("b_x", "b_x"))
  )
  for (args in malformed) {
    expect_identical(
      compatibility_error_code(do.call(
        pvstackr:::pv_compatibility_brr_fay_estimand,
        args
      )),
      "PV_BIND_E081"
    )
  }
})
