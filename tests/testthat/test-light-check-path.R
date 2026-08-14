light_check_pkg_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
}

skip_if_not_source_package <- function() {
  skip_if_not(
    file.exists(file.path(light_check_pkg_root(), "DESCRIPTION")),
    "light check path audit runs only from the package source tree"
  )
}

skip_if_no_dev_file <- function(...) {
  skip_if_not(
    file.exists(file.path(light_check_pkg_root(), ...)),
    "maintainer dev/ scripts are not distributed with the package"
  )
}

light_check_read <- function(...) {
  paste(readLines(file.path(light_check_pkg_root(), ...), warn = FALSE), collapse = "\n")
}

light_check_suggested_packages <- function(desc) {
  suggests <- if ("Suggests" %in% names(desc)) desc[["Suggests"]] else ""
  parts <- trimws(unlist(strsplit(suggests, ",")))
  parts <- sub("\\s*\\(.*\\)$", "", parts)
  parts[nzchar(parts)]
}

light_check_package_list <- function(field) {
  parts <- trimws(unlist(strsplit(field, ",")))
  parts <- sub("\\s*\\(.*\\)$", "", parts)
  parts[nzchar(parts)]
}

test_that("DESCRIPTION keeps default install dependencies light", {
  skip_if_not_source_package()
  desc <- read.dcf(file.path(light_check_pkg_root(), "DESCRIPTION"))[1L, ]
  depends <- if ("Depends" %in% names(desc)) desc[["Depends"]] else ""
  imports <- if ("Imports" %in% names(desc)) desc[["Imports"]] else ""
  linking_to <- if ("LinkingTo" %in% names(desc)) desc[["LinkingTo"]] else ""
  enhances <- if ("Enhances" %in% names(desc)) desc[["Enhances"]] else ""
  suggests <- light_check_suggested_packages(desc)
  heavy <- c(
    "brms", "cmdstanr", "posterior", "loo", "EdSurvey", "haven", "lme4",
    "rstan", "StanHeaders", "future", "furrr", "pkgdown", "quarto"
  )

  expect_match(depends, "R \\(>=", fixed = FALSE)
  expect_false(any(heavy %in% light_check_package_list(depends)))
  expect_setequal(light_check_package_list(imports), "digest")
  expect_equal(linking_to, "")
  expect_equal(enhances, "")
  expect_setequal(suggests, c("brms", "cmdstanr", "posterior", "knitr", "rmarkdown", "survey", "testthat"))
  # The bundled backend packages are permitted only in Suggests and remain
  # guarded by requireNamespace; everything else heavy stays out entirely.
  expect_false(any(setdiff(heavy, c("brms", "cmdstanr", "posterior")) %in% suggests))
  expect_identical(desc[["Additional_repositories"]], "https://stan-dev.r-universe.dev")
  system_requirements <- gsub("\\s+", " ", desc[["SystemRequirements"]])
  expect_match(system_requirements, "Optional CmdStan", fixed = TRUE)
  expect_match(
    system_requirements,
    "not required for package installation, loading, examples, or default tests",
    fixed = TRUE
  )
})

test_that("NAMESPACE has no runtime imports or heavy backend bindings", {
  skip_if_not_source_package()
  namespace <- readLines(file.path(light_check_pkg_root(), "NAMESPACE"), warn = FALSE)
  heavy <- c("brms", "cmdstanr", "posterior", "loo", "EdSurvey", "lme4", "rstan", "StanHeaders")

  expect_false(any(grepl("^import(From)?\\(", namespace)))
  for (pkg in heavy) {
    expect_false(any(grepl(pkg, namespace, fixed = TRUE)), info = pkg)
  }
})

test_that("runtime R code does not require optional backend packages", {
  skip_if_not_source_package()
  files <- list.files(file.path(light_check_pkg_root(), "R"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
  backends_file <- files[basename(files) == "backends.R"]
  core_files <- files[basename(files) != "backends.R"]
  text <- paste(unlist(lapply(core_files, readLines, warn = FALSE)), collapse = "\n")
  heavy <- c("brms", "cmdstanr", "posterior", "loo", "EdSurvey", "lme4", "rstan", "StanHeaders")
  heavy_alt <- paste(heavy, collapse = "|")

  # Core runtime code never touches optional backends.
  expect_false(grepl(paste0("\\b(", heavy_alt, ")::"), text))
  expect_false(grepl(paste0("\\b(", heavy_alt, "):::"), text))
  expect_false(grepl(paste0("library\\s*\\(\\s*['\"]?(", heavy_alt, ")"), text))
  expect_false(grepl(paste0("require\\s*\\(\\s*['\"]?(", heavy_alt, ")"), text))
  expect_false(grepl(paste0("requireNamespace\\s*\\(\\s*['\"](", heavy_alt, ")"), text))

  # The bundled adapters live in backends.R only, and every backend call there
  # is availability-guarded.
  expect_length(backends_file, 1L)
  backend_text <- paste(readLines(backends_file, warn = FALSE), collapse = "\n")
  expect_match(backend_text, "requireNamespace\\(\"brms\"", all = FALSE)
  expect_match(backend_text, "requireNamespace\\(\"cmdstanr\"", all = FALSE)
  expect_match(backend_text, "requireNamespace\\(\"posterior\"", all = FALSE)
  expect_false(grepl(paste0("library\\s*\\(\\s*['\"]?(", heavy_alt, ")"), backend_text))
})

test_that("default stack_direct path stays light without the bundled backend", {
  bundle <- pisa_tiny_parity_load()
  before <- loadedNamespaces()

  expect_error(
    pv_fit_direct(
      data = bundle$data,
      formula = OUTCOME ~ x + female,
      target = bundle$cached$target,
      control = pv_control(method = "stack_direct", backend = "none", iter = 20L, warmup = 10L, chains = 2L)
    ),
    "bundled brms backend"
  )

  heavy <- c("brms", "cmdstanr", "posterior", "loo", "EdSurvey", "lme4", "rstan", "StanHeaders")
  expect_equal(intersect(setdiff(loadedNamespaces(), before), heavy), character())
})

test_that("tiny fixture workflow does not load heavy optional backends", {
  before <- loadedNamespaces()
  live <- pisa_tiny_parity_live_fit()
  invisible(get_estimates(live$fit))
  invisible(get_target(live$fit))
  invisible(get_diagnostics(live$fit))
  invisible(get_draws(live$fit))

  heavy <- c("brms", "cmdstanr", "posterior", "loo", "EdSurvey", "lme4", "rstan", "StanHeaders")
  expect_equal(intersect(setdiff(loadedNamespaces(), before), heavy), character())
})

test_that("development check script preserves the light check policy", {
  skip_if_not_source_package()
  skip_if_no_dev_file("dev", "02_check.R")

  script <- light_check_read("dev", "02_check.R")

  expect_match(script, "--no-build-vignettes", fixed = TRUE)
  expect_match(script, "--no-manual", fixed = TRUE)
  expect_match(script, "_R_CHECK_FORCE_SUGGESTS_=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_RENDER_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RENDER_SITE=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_BACKEND_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_LIVE_BACKEND_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_NUMERIC_FIXTURE_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_ORACLE_TESTS=false", fixed = TRUE)
  expect_match(script, "check-build-hygiene.R", fixed = TRUE)
})

test_that("installed package-controlled payload stays light", {
  root <- system.file(package = "pvstackr")
  skip_if_not(nzchar(root), "installed package root is unavailable")
  source_inst <- file.path(light_check_pkg_root(), "inst")
  if (dir.exists(source_inst)) {
    skip_if(
      identical(normalizePath(root, mustWork = TRUE), normalizePath(source_inst, mustWork = TRUE)),
      "installed payload audit runs against the installed package, not pkgload's source inst root"
    )
  }

  forbidden_dirs <- c("dev", "log", "docs", "pkgdown", "cache", "data-cache", "pisa", "stan", "cmdstan", "brms")
  for (directory in forbidden_dirs) {
    expect_false(dir.exists(file.path(root, directory)), info = directory)
  }

  extdata <- list.files(file.path(root, "extdata"), recursive = TRUE)
  expect_setequal(extdata, c(
    "README.md",
    "examples/pisa_tiny_stack_direct.rds",
    "pisa_tiny.csv",
    "pisa_tiny_manifest.dcf"
  ))

  docs <- list.files(file.path(root, "doc"), recursive = TRUE)
  expect_true("a5-real-pisa-guidance.html" %in% docs)
  expect_true("a5-real-pisa-guidance.R" %in% docs)

  package_controlled <- c(
    list.files(file.path(root, "extdata"), recursive = TRUE, full.names = TRUE),
    list.files(file.path(root, "doc"), recursive = TRUE, full.names = TRUE),
    list.files(file.path(root, "tests"), recursive = TRUE, full.names = TRUE)
  )
  forbidden_extensions <- "[.](sav|por|sas7bdat|xpt|dta|parquet|duckdb|db|stan|cpp)$"
  expect_false(any(grepl(forbidden_extensions, package_controlled, ignore.case = TRUE)))
})

test_that("light verifier script exists and uses sentinel optional packages", {
  skip_if_not_source_package()
  skip_if_no_dev_file("dev", "verify-light-path.R")

  script <- light_check_read("dev", "verify-light-path.R")

  expect_match(script, "sentinel", fixed = TRUE)
  expect_match(script, "R CMD check", fixed = TRUE)
  expect_match(script, "--no-build-vignettes", fixed = TRUE)
  expect_match(script, "_R_CHECK_FORCE_SUGGESTS_=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_RENDER_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_BACKEND_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_LIVE_BACKEND_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_NUMERIC_FIXTURE_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_ORACLE_TESTS=false", fixed = TRUE)
  expect_match(script, "PISA_DATA_DIR", fixed = TRUE)
  for (pkg in c("brms", "cmdstanr", "posterior", "pkgdown", "quarto", "loo", "EdSurvey", "haven", "survey")) {
    expect_match(script, pkg, fixed = TRUE)
  }
})

test_that("live bundled-brms smoke remains opt-in and sentinel-fail-hard", {
  skip_if_not_source_package()
  live_test <- light_check_read("tests", "testthat", "test-live-brms-smoke.R")
  env_gate <- regexpr("PVSTACKR_RUN_LIVE_BACKEND_TESTS", live_test, fixed = TRUE)
  load_gate <- regexpr("requireNamespace(package", live_test, fixed = TRUE)

  expect_gt(env_gate, 0)
  expect_gt(load_gate, 0)
  expect_lt(env_gate, load_gate)
  expect_match(live_test, "Explicit live brms smoke requires the real", fixed = TRUE)
  expect_match(live_test, "check_cmdstan_toolchain", fixed = TRUE)
  expect_match(live_test, "diagnostic_role = \"execution_smoke\"", fixed = TRUE)
  expect_match(live_test, "statistical_acceptance = FALSE", fixed = TRUE)
  expect_false(grepl("library\\s*\\(\\s*(brms|cmdstanr|posterior)", live_test))
})

test_that("brms numeric fixture remains opt-in and sentinel-fail-hard", {
  skip_if_not_source_package()
  numeric_test <- light_check_read(
    "tests", "testthat", "test-brms-numeric-fixture.R"
  )
  helper <- light_check_read(
    "tests", "testthat", "helper-brms-numeric-fixture.R"
  )
  env_gate <- regexpr(
    "PVSTACKR_RUN_NUMERIC_FIXTURE_TESTS",
    numeric_test,
    fixed = TRUE
  )
  load_gate <- regexpr("requireNamespace(\"posterior\"", helper, fixed = TRUE)

  expect_gt(env_gate, 0)
  expect_gt(load_gate, 0)
  expect_match(
    helper,
    "requires the real posterior",
    fixed = TRUE
  )
  expect_match(helper, "package, not a sentinel", fixed = TRUE)
  expect_match(numeric_test, "diagnostic_role = \"numeric_acceptance\"", fixed = TRUE)
  expect_match(numeric_test, "numeric_contract_acceptance = TRUE", fixed = TRUE)
  expect_match(numeric_test, "empirical_backend_accuracy = FALSE", fixed = TRUE)
  expect_match(numeric_test, "fixture_origin = \"synthetic_not_brmsfit\"", fixed = TRUE)
  expect_match(numeric_test, "bundled_sampling_tested_here = FALSE", fixed = TRUE)
  expect_match(numeric_test, "live_sampler_quality_validated = FALSE", fixed = TRUE)
  expect_match(numeric_test, "model_recovery_validated = FALSE", fixed = TRUE)
  expect_match(numeric_test, "coverage_validated = FALSE", fixed = TRUE)
  expect_match(numeric_test, "real_data_evidence = FALSE", fixed = TRUE)
  expect_match(
    numeric_test,
    "sampler_reference = \"posterior_canonical_reference_parity\"",
    fixed = TRUE
  )
  expect_match(numeric_test, "diagnostic_extraction_failed", fixed = TRUE)
  expect_match(numeric_test, "sampler_divergences_blocked", fixed = TRUE)
  expect_match(numeric_test, "center_separation_red", fixed = TRUE)
  expect_false(grepl("library\\s*\\(\\s*posterior", numeric_test))
  expect_false(grepl("library\\s*\\(\\s*posterior", helper))
})

test_that("survey oracle checks remain opt-in and sentinel-safe", {
  skip_if_not_source_package()
  oracle_test <- light_check_read("tests", "testthat", "test-brr-fay-target.R")
  env_gate <- regexpr("PVSTACKR_RUN_ORACLE_TESTS", oracle_test, fixed = TRUE)
  load_gate <- regexpr("requireNamespace(\"survey\"", oracle_test, fixed = TRUE)

  expect_gt(env_gate, 0)
  expect_gt(load_gate, 0)
  expect_lt(env_gate, load_gate)
  expect_match(
    oracle_test,
    "Explicitly enabled survey oracle checks require the real survey",
    fixed = TRUE
  )
  expect_match(oracle_test, "survey::svrepdesign", fixed = TRUE)
  expect_match(oracle_test, "survey::svyglm", fixed = TRUE)
  expect_match(oracle_test, "combined.weights = TRUE", fixed = TRUE)
  expect_match(oracle_test, "mse = TRUE", fixed = TRUE)
  expect_match(oracle_test, "mse = FALSE", fixed = TRUE)
  expect_false(grepl("library\\s*\\(\\s*survey", oracle_test))
})

test_that("build hygiene blocks heavy backend and real-data directories", {
  skip_if_not_source_package()
  skip_if_no_dev_file("dev", "check-build-hygiene.R")

  script <- light_check_read("dev", "check-build-hygiene.R")

  for (directory in c("cache", "results", "data-cache", "cloud", "pisa", "stan", "cmdstan", "brms")) {
    expect_match(script, directory, fixed = TRUE)
  }
  for (extension in c("sav", "por", "sas7bdat", "xpt", "dta", "duckdb", "db")) {
    expect_match(script, extension, fixed = TRUE)
  }
})
