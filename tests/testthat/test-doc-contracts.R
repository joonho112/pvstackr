doc_contract_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
}

doc_contract_read <- function(path) {
  root <- doc_contract_root()
  full <- file.path(root, path)
  skip_if_not(file.exists(full), paste("source contract file absent:", path))
  paste(readLines(full, warn = FALSE), collapse = "\n")
}

doc_contract_flat <- function(text) {
  gsub("[[:space:]]+", " ", text)
}

test_that("help pages keep the interval and source vocabulary", {
  object_contracts <- doc_contract_read(file.path("R", "object-contracts.R"))
  fit_direct <- doc_contract_read(file.path("R", "fit-direct.R"))
  fit_reference <- doc_contract_read(file.path("R", "fit-reference.R"))
  fit_psis <- doc_contract_read(file.path("R", "fit-stack-psis.R"))
  accessors <- doc_contract_read(file.path("R", "accessors-print.R"))
  compare <- doc_contract_read(file.path("R", "compare-methods.R"))

  for (role in c(
    "descriptive_classic_rubin",
    "coverage_barnard_rubin",
    "reference_classic_rubin",
    "reference_barnard_rubin",
    "psis_classic_rubin",
    "psis_barnard_rubin"
  )) {
    expect_match(object_contracts, role, fixed = TRUE)
  }
  expect_match(object_contracts, "coverage_claim_allowed = FALSE", fixed = TRUE)
  expect_match(object_contracts, "coverage_claim_allowed = TRUE", fixed = TRUE)
  expect_match(object_contracts, "`\"external_brr_fay_rubin\"` (`stack_direct`): the target that", fixed = TRUE)
  expect_match(object_contracts, "self_normalized_log_ratios", fixed = TRUE)
  expect_match(object_contracts, "weight_ess_iid", fixed = TRUE)

  expect_match(fit_direct, "are copied from `target`", fixed = TRUE)
  expect_match(fit_direct, "not residual degrees of freedom", fixed = TRUE)
  expect_match(fit_direct, "estimated from the stacked fit", fixed = TRUE)
  expect_match(fit_direct, "requires `control$center = \"target\"`", fixed = TRUE)
  expect_match(fit_direct, "would keep the", fixed = TRUE)
  expect_match(fit_direct, "calibrate only their covariance", fixed = TRUE)
  expect_match(fit_direct, "b_sigma_*", fixed = TRUE)
  expect_match(fit_reference, "not passed as", fixed = TRUE)
  expect_match(fit_reference, "survey weights are not passed as an argument", fixed = TRUE)
  expect_match(fit_reference, "b_sigma_*", fixed = TRUE)
  expect_match(fit_psis, "b_sigma_*", fixed = TRUE)
  expect_match(accessors, "pisa_tiny_stack_direct.rds", fixed = TRUE)
  expect_match(accessors, "returns the table of fixed-effect estimates", fixed = TRUE)
  expect_match(accessors, "has no target object", fixed = TRUE)
  expect_match(object_contracts, "used by no fitting function", fixed = TRUE)
  expect_match(object_contracts, "fits must use `control$center = \"target\"`", fixed = TRUE)
  expect_match(object_contracts, "not automatically", fixed = TRUE)
  expect_match(object_contracts, "additional_args", fixed = TRUE)
  expect_match(fit_reference, "posterior covariance of the draws, not a design-based", fixed = TRUE)
  expect_match(fit_reference, "coverage_claim_allowed = FALSE", fixed = TRUE)
  expect_match(fit_psis, "does not run Pareto smoothing", fixed = TRUE)
  expect_match(fit_psis, "how the weights were made, but it does not check them", fixed = TRUE)
  expect_match(fit_psis, "owned_stamp_bounded_projection", fixed = TRUE)
  expect_match(accessors, "but that is only a label", fixed = TRUE)
  expect_match(compare, "only when all available", fixed = TRUE)
  expect_match(compare, "has no target object", fixed = TRUE)
})

test_that("public CCC docs say delta_c_max sets the center status and delta_c_rel does not", {
  public_paths <- c(
    file.path("vignettes", "m4-ccc-calibration.Rmd"),
    file.path("vignettes", "a3-reading-results.Rmd"),
    file.path("R", "object-contracts.R")
  )
  paths <- c(
    public_paths,
    file.path("dev", "facts-and-notation.md"),
    file.path("dev", "method-track-outline.md")
  )

  all_text <- paste(vapply(paths, doc_contract_read, character(1)), collapse = "\n")
  public_flat <- doc_contract_flat(
    paste(vapply(public_paths, doc_contract_read, character(1)), collapse = "\n")
  )

  # The center status is set from delta_c_max (thresholds 0.01 and 0.05);
  # delta_c_rel is the root mean square of the same ratios and is only recorded.
  expect_match(public_flat, "delta_c_max", fixed = TRUE)
  expect_true(grepl(
    "center_status[^.]{0,40}delta_c_max|delta_c_max[^.]{0,200}0\\.01",
    public_flat,
    perl = TRUE
  ))
  expect_true(grepl(
    "delta_c_rel[^.]{0,180}(RMS|root-mean-square|root mean square|descriptive)",
    public_flat,
    perl = TRUE
  ))

  stale_fixed <- c(
    "RMS reportable gate",
    "the gate is the RMS one",
    "package currently gates",
    "package reportable gate uses the **RMS",
    "package gate: **RMS",
    "`delta_c_rel` is the reportable gate"
  )
  for (needle in stale_fixed) {
    expect_false(grepl(needle, all_text, fixed = TRUE))
  }
})

test_that("public coverage wording is scoped to provenance rather than guarantees", {
  paths <- c(
    "README.Rmd",
    "NEWS.md",
    file.path("vignettes", "a1-getting-started.Rmd"),
    file.path("vignettes", "a2-the-workflow.Rmd"),
    file.path("vignettes", "m2-brr-fay-target.Rmd")
  )

  all_text <- paste(vapply(paths, doc_contract_read, character(1)), collapse = "\n")
  flat <- doc_contract_flat(all_text)

  stale_fixed <- c(
    "coverage guarantees",
    "This target is what makes",
    "target is what makes",
    "Calibrating to the external",
    "makes `stack_direct` intervals coverage-claimable",
    "what makes a `stack_direct` fit coverage-claimable",
    "genuine coverage claim"
  )
  for (needle in stale_fixed) {
    expect_false(grepl(needle, flat, fixed = TRUE))
  }

  # The interval rule, stated in the README and in "Getting started".
  expect_match(flat, "only a `stack_direct` fit whose target uses Barnard", fixed = TRUE)
  expect_match(flat, "confidence intervals with nominal coverage", fixed = TRUE)
  expect_match(flat, "coverage-claimable only when backed", fixed = TRUE)
  expect_match(flat, "not on CCC arithmetic alone", fixed = TRUE)
})

test_that("public docs do not claim standard-error equivalence for package methods", {
  paths <- c(
    file.path("vignettes", "a1-getting-started.Rmd"),
    file.path("vignettes", "a5-real-pisa-guidance.Rmd"),
    file.path("dev", "facts-and-notation.md")
  )

  all_text <- paste(vapply(paths, doc_contract_read, character(1)), collapse = "\n")
  flat <- doc_contract_flat(all_text)

  stale_fixed <- c(
    "point estimates and standard errors match",
    "standard errors match the orthodox",
    "genuine coverage claim",
    "point estimate and standard error agreeing to about machine precision",
    "point + SE agreement"
  )
  for (needle in stale_fixed) {
    expect_false(grepl(needle, flat, fixed = TRUE))
  }
})

test_that("CCC docs describe plain Cholesky without automatic repair and never re-introduce the Tikhonov-regularization claim", {
  m4    <- doc_contract_read(file.path("vignettes", "m4-ccc-calibration.Rmd"))
  facts <- doc_contract_read(file.path("dev", "facts-and-notation.md"))
  spec  <- doc_contract_read(file.path("dev", "docs-spec.md"))
  mt    <- doc_contract_read(file.path("dev", "method-track-outline.md"))

  # FORBIDDEN. The package implements fail-fast plain Cholesky: R/ccc.R `ccc_chol()`
  # calls `chol(..., pivot = FALSE)` and aborts on a non-PD matrix, with policy
  # `target_repair = "forbidden"` (R/utils-validate.R rejects every near-PD-repair
  # path as reserved/unsupported). The false claim that production CCC
  # Tikhonov-/near-PD-regularizes Sigma_raw is a prior-lineage behaviour and must
  # never return to the docs. These are PRECISE multi-word fixed phrases; bare
  # "Tikhonov" is intentionally allowed because the honest docs use it in negations
  # ("no automatic Tikhonov or near-PD repair") and in the facts/method-track
  # lineage notes.
  all_text <- paste(c(m4, facts, spec, mt), collapse = "\n")
  flat     <- doc_contract_flat(all_text)
  forbidden <- c(
    "Tikhonov-regularized",
    "Cholesky-plus-Tikhonov",
    "Tikhonov regularization of",
    "deterministic algorithm (Tikhonov"
  )
  for (needle in forbidden) {
    expect_false(grepl(needle, flat, fixed = TRUE))
    expect_false(grepl(needle, all_text, fixed = TRUE))
  }

  # REQUIRED. The public M4 vignette says that pvstackr uses plain Cholesky
  # factorizations, stops with an error on a matrix that is not positive definite,
  # adds no ridge, and reports the kappa_A conditioning diagnostic, in both its
  # code-token and LaTeX forms.
  m4_flat <- doc_contract_flat(m4)
  expect_match(m4_flat, "pvstackr stops with an error", fixed = TRUE)
  expect_match(m4_flat, "positive definite", fixed = TRUE)
  expect_match(m4_flat, "automatic Tikhonov or near-PD repair", fixed = TRUE)
  expect_match(m4_flat, "plain Cholesky factorizations", fixed = TRUE)
  expect_match(m4_flat, "it adds no ridge", fixed = TRUE)
  expect_match(m4_flat, "kappa_A", fixed = TRUE)
  expect_match(m4_flat, "$\\kappa_A$", fixed = TRUE)

  # REQUIRED. The dev specs name the fail-fast algorithm, and the facts sheet
  # quarantines any Tikhonov wording to the prior lineage.
  expect_match(doc_contract_flat(facts), "fail-fast plain Cholesky", fixed = TRUE)
  expect_match(
    doc_contract_flat(facts),
    "any Tikhonov / near-PD wording belongs to the prior lineage, not the package",
    fixed = TRUE
  )
  expect_match(
    doc_contract_flat(spec),
    "deterministic fail-fast Cholesky algorithm",
    fixed = TRUE
  )
  expect_match(
    doc_contract_flat(mt),
    "no Tikhonov or near-PD repair",
    fixed = TRUE
  )
})

test_that("method-comparison interval cell is conditional, not a bare coverage claim", {
  a4    <- doc_contract_read(file.path("vignettes", "a4-comparing-methods.Rmd"))
  m5    <- doc_contract_read(file.path("vignettes", "m5-methods-and-coverage.Rmd"))
  facts <- doc_contract_read(file.path("dev", "facts-and-notation.md"))

  # FORBIDDEN. The bare bolded interval cell "**coverage-claimable**"
  # immediately followed by the next column delimiter ("coverage-claimable** |")
  # over-reads: a stack_direct interval is coverage-claimable ONLY against an
  # external Barnard-Rubin BRR-Fay target, and is DESCRIPTIVE under classic Rubin df
  # (the bundled fixture is the latter). This exact substring occurs ONLY in the
  # three method-comparison table rows -- every prose use of "coverage-claimable" is
  # followed by other words, never by " |" -- so after the fix it must be absent.
  for (txt in list(a4, m5, facts)) {
    expect_false(grepl("coverage-claimable** |", txt, fixed = TRUE))
  }

  # REQUIRED. Each file's interval cell must carry the conditional wording. In
  # these three files the phrase below occurs only in the stack_direct cell of
  # the method table, so this check tests the table cell.
  for (txt in list(a4, m5, facts)) {
    expect_match(txt, "coverage-claimable only with an external Barnard", fixed = TRUE)
  }
})

test_that("public PSIS comparison vignette keeps placeholder weights blocked", {
  a4 <- doc_contract_read(file.path("vignettes", "a4-comparing-methods.Rmd"))

  expect_match(a4, "pvstackr blocks the fit", fixed = TRUE)
  expect_match(a4, "provenance_incomplete", fixed = TRUE)
  expect_match(a4, "a producer name alone is not evidence that PSIS was run", fixed = TRUE)
  expect_false(grepl("synthetic_vignette_fixture", a4, fixed = TRUE))
  expect_false(grepl("psis_producer =", a4, fixed = TRUE))
  expect_false(grepl("status == \"ok\"", a4, fixed = TRUE))
})

test_that("M2 prints the Barnard-Rubin formula and does not evaluate it on the classic-df example", {
  m2 <- doc_contract_flat(doc_contract_read(file.path("vignettes", "m2-brr-fay-target.Rmd")))

  # The closed form that the code implements (R/rubin-pool.R:74-117) is printed in
  # the Method track, with the coefficient index ell.
  expect_match(m2, "nu_{\\text{BR},\\ell}", fixed = TRUE)
  expect_match(m2, "nu_{\\text{obs},\\ell}", fixed = TRUE)
  expect_match(m2, "(M-1)/\\lambda_\\ell^{2}", fixed = TRUE)

  # The formula is printed but not evaluated on the example target, which uses
  # classic degrees of freedom, so no Barnard-Rubin value is shown.
  expect_match(m2, "prints the formula but does", fixed = TRUE)
  expect_match(m2, "The example target uses classic degrees of freedom", fixed = TRUE)

  # It must not revert to declining to show the formula.
  expect_false(grepl("neither prints a formula nor", m2, fixed = TRUE))
})
