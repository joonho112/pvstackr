test_that("DESCRIPTION records package identity and light dependency policy", {
  desc <- utils::packageDescription("pvstackr")
  expected_title <- "Stacked-Fit Calibration to Rubin/BRR-Fay Fixed-Effect Targets for Plausible Values"

  expect_identical(desc$Package, "pvstackr")
  expect_identical(gsub("\\s+", " ", desc$Title), expected_title)
  desc_text <- gsub("\\s+", " ", desc$Description)
  expect_match(desc_text, "Bayesian-backend fixed-effect calibration and reporting", fixed = TRUE)
  expect_match(desc_text, "one-stacked-fit workflows", fixed = TRUE)
  expect_match(desc_text, "Rubin/BRR-Fay fixed-effect", fixed = TRUE)
  expect_match(desc_text, "calibrating fixed-effect draws", fixed = TRUE)
  expect_false(grepl("expensive Bayesian sampler", desc_text, fixed = TRUE))
  expect_match(desc$URL, "https://github.com/joonho112/pvstackr", fixed = TRUE)
  expect_identical(desc$BugReports, "https://github.com/joonho112/pvstackr/issues")
  expect_match(desc$Suggests, "brms")
  expect_match(desc$Suggests, "cmdstanr")
  expect_match(desc$Suggests, "posterior")
  expect_match(desc$Suggests, "testthat")
  expect_match(desc$Suggests, "knitr")
  expect_match(desc$Suggests, "rmarkdown")
  expect_match(desc$Suggests, "survey")
  expect_identical(desc$Additional_repositories, "https://stan-dev.r-universe.dev")
  system_requirements <- gsub("\\s+", " ", desc$SystemRequirements)
  expect_match(system_requirements, "Optional CmdStan", fixed = TRUE)
  expect_match(system_requirements, "not required for package installation", fixed = TRUE)
  expect_identical(desc$VignetteBuilder, "knitr")

  imports <- desc$Imports
  if (is.null(imports)) imports <- ""
  expect_match(imports, "digest", fixed = TRUE)
  depends <- desc$Depends
  if (is.null(depends)) depends <- ""
  linking_to <- desc$LinkingTo
  if (is.null(linking_to)) linking_to <- ""
  enhances <- desc$Enhances
  if (is.null(enhances)) enhances <- ""
  heavy <- c("brms", "cmdstanr", "EdSurvey", "loo", "quarto", "future", "furrr", "survey")
  expect_false(any(vapply(heavy, grepl, logical(1), x = imports, fixed = TRUE)))
  expect_false(grepl("survey", depends, fixed = TRUE))
  expect_false(grepl("survey", linking_to, fixed = TRUE))
  expect_false(grepl("survey", enhances, fixed = TRUE))

  desc_path <- file.path(system.file(package = "pvstackr"), "DESCRIPTION")
  if (!file.exists(desc_path)) {
    desc_path <- test_path("../../DESCRIPTION")
  }
  desc_lines <- readLines(desc_path, warn = FALSE)
  expect_false(any(grepl("TODO@example.com", desc_lines, fixed = TRUE)))
  expect_true(any(grepl('person\\("JoonHo", "Lee"', desc_lines)))
  expect_true(any(grepl("jlee296@ua.edu", desc_lines, fixed = TRUE)))
  expect_true(any(grepl("ORCID = \"0009-0006-4019-8703\"", desc_lines, fixed = TRUE)))

  source_desc_path <- file.path(getwd(), "DESCRIPTION")
  if (file.exists(source_desc_path) && dir.exists(file.path(getwd(), "R"))) {
    source_desc_lines <- readLines(source_desc_path, warn = FALSE)
    expect_false(any(grepl("^Author:", source_desc_lines)))
	  expect_false(any(grepl("^Maintainer:", source_desc_lines)))
  }
})

test_that("README source describes survey as optional oracle infrastructure", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  readme_path <- file.path(root, "README.Rmd")
  skip_if_not(file.exists(readme_path), "README.Rmd absent from installed package")

  readme <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")
  flat <- gsub("[[:space:]]+", " ", readme)

  expect_false(grepl("\\b(relies on|depends on|requires)\\s+\\[?survey\\b", flat, perl = TRUE, ignore.case = TRUE))
  expect_match(flat, "runtime BRR-Fay target engine is dependency-free", fixed = TRUE)
  expect_true(grepl("optional development tests[^.]{0,160}survey", flat, perl = TRUE))
})
