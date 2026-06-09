light_check_pkg_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
}

skip_if_not_source_package <- function() {
  skip_if_not(
    file.exists(file.path(light_check_pkg_root(), "DESCRIPTION")),
    "light check path audit runs only from the package source tree"
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
  expect_equal(imports, "")
  expect_equal(linking_to, "")
  expect_equal(enhances, "")
  expect_setequal(suggests, c("knitr", "rmarkdown", "survey", "testthat"))
  expect_false(any(heavy %in% suggests))
  expect_match(desc[["SystemRequirements"]], "Optional CmdStan", fixed = TRUE)
  expect_match(
    desc[["SystemRequirements"]],
    "not\\s+required for installation, loading, examples, or default tests"
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
  text <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  heavy <- c("brms", "cmdstanr", "posterior", "loo", "EdSurvey", "lme4", "rstan", "StanHeaders")
  heavy_alt <- paste(heavy, collapse = "|")

  expect_false(grepl(paste0("\\b(", heavy_alt, ")::"), text))
  expect_false(grepl(paste0("\\b(", heavy_alt, "):::"), text))
  expect_false(grepl(paste0("library\\s*\\(\\s*['\"]?(", heavy_alt, ")"), text))
  expect_false(grepl(paste0("require\\s*\\(\\s*['\"]?(", heavy_alt, ")"), text))
  expect_false(grepl(paste0("requireNamespace\\s*\\(\\s*['\"](", heavy_alt, ")"), text))
})

test_that("default stack_direct path stops at injected backend boundary", {
  bundle <- pisa_tiny_parity_load()
  before <- loadedNamespaces()

  expect_error(
    pv_fit_direct(
      data = bundle$data,
      formula = OUTCOME ~ x + female,
      target = bundle$cached$target,
      control = pv_control(method = "stack_direct", backend = "brms", iter = 20L, warmup = 10L, chains = 2L)
    ),
    "`fit_function` is required"
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
  skip_if_not(
    file.exists(file.path(light_check_pkg_root(), "dev", "02_check.R")),
    "dev/ tooling is not present in this tree (installed or public checkout)"
  )
  script <- light_check_read("dev", "02_check.R")

  expect_match(script, "--no-build-vignettes", fixed = TRUE)
  expect_match(script, "--no-manual", fixed = TRUE)
  expect_match(script, "_R_CHECK_FORCE_SUGGESTS_=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_RENDER_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RENDER_SITE=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_BACKEND_TESTS=false", fixed = TRUE)
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
  skip_if_not(
    file.exists(file.path(light_check_pkg_root(), "dev", "verify-light-path.R")),
    "dev/ tooling is not present in this tree (installed or public checkout)"
  )
  script <- light_check_read("dev", "verify-light-path.R")

  expect_match(script, "sentinel", fixed = TRUE)
  expect_match(script, "R CMD check", fixed = TRUE)
  expect_match(script, "--no-build-vignettes", fixed = TRUE)
  expect_match(script, "_R_CHECK_FORCE_SUGGESTS_=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_RENDER_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_BACKEND_TESTS=false", fixed = TRUE)
  expect_match(script, "PVSTACKR_RUN_ORACLE_TESTS=false", fixed = TRUE)
  expect_match(script, "PISA_DATA_DIR", fixed = TRUE)
  for (pkg in c("brms", "cmdstanr", "posterior", "pkgdown", "quarto", "loo", "EdSurvey", "haven", "survey")) {
    expect_match(script, pkg, fixed = TRUE)
  }
})

test_that("survey oracle checks remain opt-in and sentinel-safe", {
  skip_if_not_source_package()
  oracle_test <- light_check_read("tests", "testthat", "test-brr-fay-target.R")
  env_gate <- regexpr("PVSTACKR_RUN_ORACLE_TESTS", oracle_test, fixed = TRUE)
  install_gate <- regexpr("skip_if_not_installed\\(\"survey\"\\)", oracle_test)

  expect_gt(env_gate, 0)
  expect_gt(install_gate, 0)
  expect_lt(env_gate, install_gate)
  expect_match(oracle_test, "survey oracle checks require the real survey package", fixed = TRUE)
  expect_match(oracle_test, "survey::svrepdesign", fixed = TRUE)
  expect_match(oracle_test, "survey::svyglm", fixed = TRUE)
  expect_match(oracle_test, "combined.weights = TRUE", fixed = TRUE)
  expect_match(oracle_test, "mse = TRUE", fixed = TRUE)
  expect_false(grepl("library\\s*\\(\\s*survey", oracle_test))
})

test_that("build hygiene blocks heavy backend and real-data directories", {
  skip_if_not_source_package()
  skip_if_not(
    file.exists(file.path(light_check_pkg_root(), "dev", "check-build-hygiene.R")),
    "dev/ tooling is not present in this tree (installed or public checkout)"
  )
  script <- light_check_read("dev", "check-build-hygiene.R")

  for (directory in c("cache", "results", "data-cache", "cloud", "pisa", "stan", "cmdstan", "brms")) {
    expect_match(script, directory, fixed = TRUE)
  }
  for (extension in c("sav", "por", "sas7bdat", "xpt", "dta", "duckdb", "db")) {
    expect_match(script, extension, fixed = TRUE)
  }
})
