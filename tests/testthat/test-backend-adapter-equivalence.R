adapter_equivalence_fakes <- function() {
  list(
    draws = array(
      seq_len(2L * 100L * 3L) / 97,
      dim = c(100L, 2L, 3L),
      dimnames = list(NULL, NULL, c("b_Intercept", "b_x", "sigma"))
    ),
    summary = data.frame(
      variable = c("b_Intercept", "b_x", "sigma"),
      rhat = c(1.003, 1.009, 1.005),
      ess_bulk = c(300, 220, 250),
      ess_tail = c(260, 180, 210)
    ),
    nuts = data.frame(
      Parameter = rep("divergent__", 200),
      Value = c(rep(0, 199), 1)
    )
  )
}

# The bundled adapter driven through its documented injection points, so the
# equivalence can be checked without brms, cmdstanr, or a live sampler.
adapter_equivalence_diagnose <- function() {
  fakes <- adapter_equivalence_fakes()
  function(fit) {
    pv_backend_brms_sampler_diagnostics(
      fit,
      draws_array_function = function(fit) fakes$draws,
      summary_function = function(draws) fakes$summary,
      nuts_function = function(fit) fakes$nuts
    )
  }
}

test_that("the bundled brms adapter triple is part of the public surface", {
  for (name in c(
    "pv_backend_brms_fit_function",
    "pv_backend_brms_draws_function",
    "pv_backend_brms_sampler_diagnostics"
  )) {
    expect_true(is.function(getExportedValue("pvstackr", name)), info = name)
  }

  # The article documents a twelve-argument fit contract; keep it pinned so the
  # published listing cannot drift away from the shipped adapter.
  expect_identical(
    names(formals(pv_backend_brms_fit_function)),
    c("formula", "data", "family", "prior", "chains", "iter", "warmup",
      "cores", "seed", "backend", "file", "file_refit", "...")
  )
  expect_identical(names(formals(pv_backend_brms_draws_function)), c("fit", "..."))
  expect_identical(
    names(formals(pv_backend_brms_sampler_diagnostics)),
    c("fit", "draws_array_function", "summary_function", "nuts_function")
  )
})

test_that("bundled and injected routes agree once the diagnose function is supplied", {
  diagnose <- adapter_equivalence_diagnose()

  bundled <- pvstackr:::pv_stack_sampler_diagnostics(
    fit = list(),
    bundled_backend = TRUE,
    diagnose_function = NULL,
    chains = 2L,
    iter = 200L,
    warmup = 100L,
    bundled_sampler_function = diagnose
  )
  injected <- pvstackr:::pv_stack_sampler_diagnostics(
    fit = list(),
    bundled_backend = FALSE,
    diagnose_function = diagnose,
    chains = 2L,
    iter = 200L,
    warmup = 100L
  )

  expect_true(bundled$sampler$diagnostic_complete)
  expect_true(injected$sampler$diagnostic_complete)

  # Everything the gate reads on is identical; only the provenance label, which
  # is meant to differ, separates the two routes.
  numeric_fields <- c(
    "rhat_max", "ess_bulk_min", "ess_tail_min", "ess_bulk_per_chain_min",
    "ess_tail_per_chain_min", "divergences", "chains",
    "post_warmup_draws_per_chain", "diagnostic_complete"
  )
  expect_identical(bundled$sampler[numeric_fields], injected$sampler[numeric_fields])
  expect_identical(bundled$sampler$diagnostic_source, "bundled_brms_posterior_and_nuts")
  expect_identical(injected$sampler$diagnostic_source, "injected_diagnose_function")
})

test_that("an injected pair without a diagnose function reports incomplete evidence", {
  injected <- pvstackr:::pv_stack_sampler_diagnostics(
    fit = list(),
    bundled_backend = FALSE,
    diagnose_function = NULL,
    chains = 2L,
    iter = 200L,
    warmup = 100L
  )

  expect_false(injected$sampler$diagnostic_complete)
  expect_identical(
    injected$sampler$diagnostic_source,
    "injected_diagnose_function_absent"
  )
  expect_identical(
    injected$sampler$diagnostic_reason_codes,
    "diagnostic_extractor_not_supplied"
  )
})

test_that("an injected fit function takes precedence over a bundled backend request", {
  # Supplying `fit_function` alongside `pv_control(backend = "brms")` resolves to
  # the injected route. The choice is silent at the call site, so the engine
  # record is the only place it is visible; keep that record pinned.
  spec <- pvstackr:::pv_backend_injected_engine_spec("brms")

  expect_identical(spec$requested_backend, "brms")
  expect_identical(spec$resolved_backend, "injected")
  expect_identical(spec$selection_reason, "fit_function_supplied")
  expect_identical(spec$selection_policy, "caller_supplied_fit_function")
})
