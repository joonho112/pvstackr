fixture_path <- function(...) {
  path <- system.file(..., package = "pvstackr")
  if (nzchar(path)) {
    return(path)
  }
  file.path(getwd(), "inst", ...)
}

test_that("tiny PISA-shaped CSV fixture is packaged and license-clean", {
  path <- fixture_path("extdata", "pisa_tiny.csv")
  readme_path <- fixture_path("extdata", "README.md")
  expect_true(file.exists(path))
  expect_true(file.exists(readme_path))
  readme <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")
  expect_match(readme, "no real PISA records", fixed = TRUE)

  data <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_equal(nrow(data), 12L)
  expect_equal(
    names(data),
    c(
      "CNT", "CNTSCHID", "CNTSTUID", "x", "female",
      "PV1READ", "PV2READ", "W_FSTUWT", paste0("W_FSTURWT", 1:4)
    )
  )
  expect_equal(unique(data$CNT), "SYN")
  expect_false(anyNA(data))
  expect_true(all(data$W_FSTUWT > 0))
  expect_true(all(data[paste0("W_FSTURWT", 1:4)] > 0))
  expect_identical(detect_pisa_pv_columns(data, suffix = "READ"), c("PV1READ", "PV2READ"))
  expect_identical(detect_pisa_brr_replicate_weights(data), paste0("W_FSTURWT", 1:4))
})

test_that("cached tiny stack_direct fixture is fresh against packaged CSV", {
  skip_if_not_hash_tests()

  csv_path <- fixture_path("extdata", "pisa_tiny.csv")
  manifest_path <- fixture_path("extdata", "pisa_tiny_manifest.dcf")
  rds_path <- fixture_path("extdata", "examples", "pisa_tiny_stack_direct.rds")
  expect_true(file.exists(csv_path))
  expect_true(file.exists(manifest_path))
  expect_true(file.exists(rds_path))

  data <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  cached <- readRDS(rds_path)
  manifest <- read.dcf(manifest_path)[1L, ]
  required <- c(
    "fixture_version", "generated_by", "license", "data_hash", "design_hash",
    "row_support_hash", "pv_value_hash", "weight_design_hash", "target_hash",
    "long_data_hash", "draws_hash", "fit_estimates_hash", "canonical_fit_hash",
    "design", "target", "fit"
  )
  expect_true(all(required %in% names(cached)))
  expect_equal(cached$fixture_version, "0.1.0")
  expect_equal(cached$generated_by, "dev/build-fixtures.R")
  expect_match(cached$license, "synthetic", fixed = TRUE)
  expect_equal(cached$data_hash, pvstackr:::pv_hash_payload(data))
  expect_equal(unname(tools::sha256sum(csv_path)), unname(manifest[["CSV-SHA256"]]))
  expect_equal(unname(tools::sha256sum(rds_path)), unname(manifest[["RDS-SHA256"]]))
  expect_equal(manifest[["Real-PISA-Data"]], "false")
  expect_equal(manifest[["License-Clean"]], "true")
  license_note <- gsub("[[:space:]]+", " ", manifest[["Data-License-Note"]])
  expect_match(license_note, "no OECD/PISA records", fixed = TRUE)
  expect_equal(manifest[["Generated-At"]], "omitted_for_deterministic_fixture")
  expect_equal(manifest[["Data-Hash"]], cached$data_hash)

  design <- pv_design(
    data = data,
    formula = OUTCOME ~ x + female,
    pv_suffix = "READ",
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID"
  )
  expect_equal(cached$design_hash, design$design_hash)
  expect_equal(cached$row_support_hash, design$row_support_hash)
  expect_equal(cached$pv_value_hash, design$pv_value_hash)
  expect_equal(cached$weight_design_hash, design$weight_design_hash)
  expect_equal(cached$design$design_hash, design$design_hash)
  expect_equal(manifest[["Design-Hash"]], design$design_hash)
  expect_equal(manifest[["Row-Support-Hash"]], design$row_support_hash)
  expect_equal(manifest[["PV-Value-Hash"]], design$pv_value_hash)
  expect_equal(manifest[["Weight-Design-Hash"]], design$weight_design_hash)
  expect_invisible(pvstackr:::validate_pvstackr_design(cached$design))

  target <- pv_brr_target(
    data = data,
    formula = OUTCOME ~ x + female,
    pv_cols = design$pv_cols,
    weight_col = design$weight_col,
    rep_weight_cols = design$rep_weight_cols,
    fay_k = design$fay_k,
    id_cols = design$id_cols
  )
  expect_equal(cached$target_hash, target$target_hash)
  expect_equal(cached$target$target_hash, target$target_hash)
  expect_equal(manifest[["Target-Hash"]], target$target_hash)
  expect_equal(cached$target$beta, target$beta, tolerance = 1e-12)
  expect_equal(cached$target$T_MI, target$T_MI, tolerance = 1e-12)
  expect_invisible(pvstackr:::validate_pvstackr_brr_target(cached$target))

  expect_s3_class(cached$fit, "pvstackr_fit")
  expect_equal(cached$fit$method, "stack_direct")
  expect_equal(cached$fit$status, "ok")
  expect_equal(cached$fit$target$target_hash, target$target_hash)
  expect_equal(cached$long_data_hash, cached$fit$stack_fit$weight_summary$long_data_hash)
  expect_equal(cached$draws_hash, pvstackr:::pv_hash_payload(cached$fit$stack_fit$stacked_draws))
  expect_equal(cached$fit_estimates_hash, pvstackr:::pv_hash_payload(cached$fit$estimates))
  expect_equal(cached$canonical_fit_hash, pvstackr:::pv_hash_payload(list(
    method = cached$fit$method,
    status = cached$fit$status,
    terms = cached$fit$estimates$term,
    estimates = cached$fit$estimates,
    target_hash = target$target_hash,
    long_data_hash = cached$fit$stack_fit$weight_summary$long_data_hash
  )))
  expect_equal(manifest[["Long-Data-Hash"]], cached$long_data_hash)
  expect_equal(manifest[["Draws-Hash"]], cached$draws_hash)
  expect_equal(manifest[["Estimates-Hash"]], cached$fit_estimates_hash)
  expect_equal(manifest[["Canonical-Fit-Hash"]], cached$canonical_fit_hash)
  expect_equal(cached$fit$estimates$term, target$fe_names)
  expect_equal(cached$fit$estimates$target_hash, rep(target$target_hash, length(target$fe_names)))
  expect_invisible(pvstackr:::validate_pvstackr_fit(cached$fit))
})

test_that("tiny fixture files stay inside the explicit build-hygiene allowlist", {
  extdata_root <- fixture_path("extdata")
  expected <- c(
    "README.md",
    "examples/pisa_tiny_stack_direct.rds",
    "pisa_tiny.csv",
    "pisa_tiny_manifest.dcf"
  )
  expect_equal(sort(list.files(extdata_root, recursive = TRUE)), sort(expected))
  expect_equal(basename(fixture_path("extdata", "README.md")), "README.md")
  expect_equal(basename(fixture_path("extdata", "pisa_tiny.csv")), "pisa_tiny.csv")
  expect_equal(basename(fixture_path("extdata", "pisa_tiny_manifest.dcf")), "pisa_tiny_manifest.dcf")
  expect_equal(
    basename(fixture_path("extdata", "examples", "pisa_tiny_stack_direct.rds")),
    "pisa_tiny_stack_direct.rds"
  )
})
