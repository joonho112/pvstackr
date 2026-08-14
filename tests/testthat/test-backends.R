test_that("stack fit without fit_function aborts unless the brms backend is selected", {
  dat <- data.frame(
    PV1 = rnorm(6), PV2 = rnorm(6),
    x = rnorm(6), w = rep(1, 6)
  )
  expect_error(
    pv_stack_fit(dat, OUTCOME ~ x, pv_cols = c("PV1", "PV2"),
                 weight_col = "w", control = pv_control(backend = "none")),
    "bundled brms backend"
  )
})

test_that("bundled brms draws function returns a plain named base matrix", {
  skip_if_not_installed("posterior")
  fake <- posterior::as_draws_matrix(posterior::example_draws())
  colnames(fake) <- c("b_Intercept", "b_x", "sigma", paste0("junk", seq_len(ncol(fake) - 3L)))
  out <- pv_backend_brms_draws_function(fake)
  expect_true(is.matrix(out) && is.numeric(out))
  expect_identical(class(out), c("matrix", "array"))
  expect_identical(colnames(out), c("b_Intercept", "b_x", "sigma"))
  expect_null(rownames(out))
  expect_equal(out[, "b_x"], as.numeric(fake[, "b_x"]))
})

cmdstan_state_fixture <- function(namespace_available = TRUE,
                                  cmdstan_configured = TRUE) {
  list(
    namespace_available = namespace_available,
    package_version = if (namespace_available) "0.8.0" else NA_character_,
    cmdstan_configured = cmdstan_configured,
    cmdstan_version = if (cmdstan_configured) "2.38.0" else NA_character_,
    cmdstan_path_basename = if (cmdstan_configured) "cmdstan-2.38.0" else NA_character_,
    state_reason = if (cmdstan_configured) {
      "cmdstanr_namespace_and_cmdstan_configured"
    } else {
      "cmdstanr_namespace_unavailable"
    },
    toolchain_checked = FALSE
  )
}

test_that("brms engine resolution is deterministic from CmdStan state", {
  cmdstanr_resolution <- pvstackr:::pv_backend_resolve_brms_engine(
    cmdstan_state_fixture()
  )
  expect_identical(cmdstanr_resolution$resolved_backend, "cmdstanr")
  expect_identical(cmdstanr_resolution$selection_reason, "configured_cmdstan_selected")

  rstan_resolution <- pvstackr:::pv_backend_resolve_brms_engine(
    cmdstan_state_fixture(
      namespace_available = FALSE,
      cmdstan_configured = FALSE
    )
  )
  expect_identical(rstan_resolution$resolved_backend, "rstan")
  expect_identical(
    rstan_resolution$selection_reason,
    "cmdstanr_namespace_absent_rstan_selected"
  )
})

test_that("brms engine resolution rejects ambiguous or incomplete state", {
  expect_error(
    pvstackr:::pv_backend_resolve_brms_engine(
      cmdstan_state_fixture(cmdstan_configured = FALSE)
    ),
    "does not silently retry"
  )

  incomplete <- cmdstan_state_fixture()
  incomplete$cmdstan_version <- NULL
  expect_error(
    pvstackr:::pv_backend_resolve_brms_engine(incomplete),
    "state is incomplete"
  )

  inconsistent <- cmdstan_state_fixture(namespace_available = FALSE)
  expect_error(
    pvstackr:::pv_backend_resolve_brms_engine(inconsistent),
    "cannot be configured"
  )
})

test_that("injected engine provenance does not resolve an optional backend", {
  spec <- pvstackr:::pv_backend_injected_engine_spec("cmdstanr")

  expect_identical(spec$adapter_source, "injected")
  expect_identical(spec$adapter_id, "injected_fit_function")
  expect_identical(spec$requested_backend, "cmdstanr")
  expect_identical(spec$resolved_backend, "injected")
  expect_identical(spec$cmdstan_state$state_reason, "not_evaluated_for_injected_adapter")
})
