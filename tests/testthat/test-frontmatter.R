pkg_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
}

skip_if_not_source_tree <- function() {
  skip_if_not(
    file.exists(file.path(pkg_root(), "DESCRIPTION")),
    "source-tree frontmatter checks run only from package source"
  )
}

test_that("frontmatter files exist and record scope", {
  skip_if_not_source_tree()
  root_files <- c("LICENSE", "NEWS.md", "README.md")
  expect_true(all(file.exists(file.path(pkg_root(), root_files))))
  expect_true(file.exists(file.path(pkg_root(), "inst", "CITATION")))
  expect_true(file.exists(file.path(pkg_root(), ".gitignore")))
  expect_false(file.exists(file.path(pkg_root(), "README.html")))

  gitignore <- paste(readLines(file.path(pkg_root(), ".gitignore"), warn = FALSE), collapse = "\n")
  rbuildignore <- paste(readLines(file.path(pkg_root(), ".Rbuildignore"), warn = FALSE), collapse = "\n")
  expect_match(gitignore, ".DS_Store", fixed = TRUE)
  expect_match(gitignore, "README.html", fixed = TRUE)
  expect_match(rbuildignore, "^\\.gitignore$", fixed = TRUE)

  desc <- utils::packageDescription("pvstackr")
  expect_identical(desc$License, "MIT + file LICENSE")
  expect_identical(desc$Language, "en-US")

  readme <- paste(readLines(file.path(pkg_root(), "README.md"), warn = FALSE), collapse = "\n")
  readme_flat <- gsub("[[:space:]]+", " ", readme)
  expect_match(readme, "stack_direct", fixed = TRUE)
  expect_match(readme, "stack_psis", fixed = TRUE)
  expect_match(readme, "per_pv", fixed = TRUE)
  expect_match(readme, "Bayesian-backend fixed-effect calibration\\s+and reporting")
  expect_match(readme, "design-based external Rubin/BRR-Fay fixed-effect target for reported estimates and SEs", fixed = TRUE)
  expect_match(readme, "Per-PV Bayesian/backend reference", fixed = TRUE)
  expect_match(readme, "model-based within-PV covariance", fixed = TRUE)
  expect_match(readme, "supplied/precomputed or injected PSIS", fixed = TRUE)
  expect_match(readme_flat, "do not construct or calibrate to the external BRR-Fay target", fixed = TRUE)
  expect_match(readme, "Rubin/BRR-Fay fixed-effect target", fixed = TRUE)
  expect_match(readme, "calibrated reporting is scoped to fixed", fixed = TRUE)
  expect_match(readme, "not a full posterior for all model parameters", fixed = TRUE)
  expect_match(readme, "Under the default `control$center = \"target\"`", fixed = TRUE)
  expect_match(readme, "reported fixed-effect\\s+estimate and standard error are the external Rubin/BRR-Fay target\\s+values")
  expect_match(readme_flat, "calibrated draws and agreement diagnostics", fixed = TRUE)
  expect_match(readme, "Reportable `stack_direct` fits require `control$center = \"target\"`", fixed = TRUE)
  expect_match(readme, "posterior-centered CCC convention is reserved for diagnostic/exploratory", fixed = TRUE)
  expect_match(readme, "Five-Minute Synthetic Example", fixed = TRUE)
  expect_match(readme, "pisa_tiny.csv", fixed = TRUE)
  expect_match(readme, "pisa_tiny_stack_direct.rds", fixed = TRUE)
  expect_match(readme, "no real PISA records", fixed = TRUE)
  expect_match(readme, "not suitable for real inference", fixed = TRUE)
  expect_match(readme, "detect_pisa_pv_columns", fixed = TRUE)
  expect_match(readme, "pv_brr_target", fixed = TRUE)
  expect_match(readme, "get_target", fixed = TRUE)
  expect_match(readme, "get_diagnostics", fixed = TRUE)
  expect_match(readme, "get_draws", fixed = TRUE)
  expect_match(readme, "Conceptual Live Workflow", fixed = TRUE)
  expect_match(readme, "real-data/live-backend\\s+sketch")
  expect_match(readme, "vignette(\"a5-real-pisa-guidance\"", fixed = TRUE)
  expect_match(readme, "PV1MATH", fixed = TRUE)
  expect_match(readme, "pv_suffix = \"MATH\"", fixed = TRUE)
  expect_match(readme, "bare plausible-value names[[:space:]]+such[[:space:]]+as `PV1`")
  expect_match(readme, "N \\* M")
  expect_match(readme, "M \\* \\(R \\+ 1\\)")
  expect_match(readme, "no real PISA records", fixed = TRUE)
  expect_false(grepl("replicates PISA results", readme, fixed = TRUE))
  expect_false(grepl("real PISA benchmark", readme, fixed = TRUE))
  expect_false(grepl("per_pv` | Orthodox reference: one Bayesian fit per plausible value plus Rubin/BRR-Fay pooling", readme, fixed = TRUE))
  expect_false(grepl("per_pv.*Rubin/BRR-Fay pooling", readme))
})

test_that("README tiny workflow smoke code stays lightweight", {
  csv_path <- system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
  rds_path <- system.file("extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr")
  expect_true(nzchar(csv_path))
  expect_true(nzchar(rds_path))

  pisa_tiny <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  design <- pv_design(
    data = pisa_tiny,
    formula = OUTCOME ~ x + female,
    pv_suffix = "READ",
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID"
  )
  expect_s3_class(design, "pvstackr_design")
  expect_equal(design$design_hash, "4776c0f1")

  cached <- readRDS(rds_path)
  fit <- cached$fit
  expect_s3_class(fit, "pvstackr_fit")
  expect_equal(fit$method, "stack_direct")
  expect_equal(fit$status, "ok")
  expect_equal(fit$design$design_hash, design$design_hash)
  expect_equal(fit$target$target_hash, cached$target_hash)

  estimates <- get_estimates(fit)
  expect_equal(estimates$term, c("b_Intercept", "b_x", "b_female"))
  expect_true(all(c("estimate", "se", "df", "conf_low", "conf_high") %in% names(estimates)))
  expect_true(all(is.finite(estimates$estimate)))
  expect_true(all(is.finite(estimates$se)))

  diagnostics <- get_diagnostics(fit)
  expect_true(all(c("preflight", "stack_fit", "stack_fit_warnings", "ccc") %in% names(diagnostics)))

  expect_equal(get_target(fit)$target_hash, cached$target_hash)
  expect_null(get_draws(fit))
})

test_that("README.Rmd renders to github markdown when rmarkdown is available", {
  skip_if_not_source_tree()
  skip_if_not(
    identical(Sys.getenv("PVSTACKR_RUN_RENDER_TESTS"), "true"),
    "render checks are dev-only; set PVSTACKR_RUN_RENDER_TESTS=true"
  )
  testthat::skip_if_not_installed("rmarkdown")

  out_dir <- tempfile("pvstackr-readme-")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rendered <- rmarkdown::render(
    input = file.path(pkg_root(), "README.Rmd"),
    output_format = "github_document",
    output_file = "README.md",
    output_dir = out_dir,
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
  expect_true(file.exists(rendered))
  readme <- paste(readLines(rendered, warn = FALSE), collapse = "\n")
  expect_match(readme, "pvstackr design", fixed = TRUE)
  expect_match(readme, "pvstackr BRR-Fay target", fixed = TRUE)
  expect_match(readme, "pvstackr fit", fixed = TRUE)
  expect_match(readme, "b_Intercept", fixed = TRUE)
})

test_that("repository-only frontmatter exists in the local source tree", {
  skip_if_not_source_tree()
  repo_files <- c("LICENSE.md", "README.Rmd")
  missing <- repo_files[!file.exists(file.path(pkg_root(), repo_files))]
  skip_if(length(missing) > 0L, paste("repository-only files absent from built package:", paste(missing, collapse = ", ")))
  expect_true(all(file.exists(file.path(pkg_root(), repo_files))))
})

test_that("citation file is parseable", {
  skip_if_not_source_tree()
  citation_file <- file.path(pkg_root(), "inst", "CITATION")
  env <- new.env(parent = asNamespace("utils"))
  expect_error(sys.source(citation_file, envir = env), NA)

  citation_text <- paste(readLines(citation_file, warn = FALSE), collapse = "\n")
  expect_match(
    citation_text,
    "pvstackr: Stacked-Fit Calibration to Rubin/BRR-Fay Fixed-Effect Targets for Plausible Values",
    fixed = TRUE
  )
  expect_false(grepl("pvstackr: One-MCMC Bayesian Inference for Plausible Values", citation_text, fixed = TRUE))
})

test_that("repository-only frontmatter is excluded from package build", {
  skip_if_not_source_tree()
  ignore_file <- file.path(pkg_root(), ".Rbuildignore")
  skip_if_not(file.exists(ignore_file), ".Rbuildignore is not included in built package tests")
  ignore <- readLines(ignore_file, warn = FALSE)
  expect_true("^LICENSE\\.md$" %in% ignore)
  expect_true("^README\\.Rmd$" %in% ignore)
  expect_true("^README\\.html$" %in% ignore)
})
