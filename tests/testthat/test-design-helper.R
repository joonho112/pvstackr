design_helper_fixture_data <- function() {
  data.frame(
    CNTSTUID = sprintf("stu%02d", 1:10),
    school = factor(rep(c("A", "B"), each = 5)),
    x = c(-2, -1.5, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3),
    PV1READ = c(2.1, 2.4, 3.0, 3.2, 3.6, 4.0, 4.2, 4.8, 5.1, 5.5),
    PV2READ = c(2.0, 2.5, 2.9, 3.4, 3.5, 4.1, 4.4, 4.7, 5.0, 5.8),
    W_FSTUWT = c(1.0, 1.1, 0.9, 1.2, 1.0, 1.3, 0.8, 1.4, 1.1, 0.95),
    W_FSTURWT1 = c(0.8, 1.3, 1.0, 1.1, 0.9, 1.5, 0.9, 1.2, 1.2, 1.0),
    W_FSTURWT2 = c(1.2, 0.9, 0.8, 1.4, 1.1, 1.0, 0.7, 1.6, 1.0, 0.9),
    W_FSTURWT3 = c(1.1, 1.0, 1.1, 1.0, 1.3, 1.2, 0.8, 1.3, 1.4, 1.1),
    W_FSTURWT4 = c(0.9, 1.2, 0.95, 1.3, 1.0, 1.4, 1.1, 1.1, 0.9, 1.2)
  )
}

test_that("pv_design detects PISA-style columns and records hashes", {
  data <- design_helper_fixture_data()
  design <- pv_design(
    data = data,
    formula = OUTCOME ~ x + (1 | school),
    pv_suffix = "READ",
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID"
  )

  expect_s3_class(design, "pvstackr_design")
  expect_equal(design$pv_cols, c("PV1READ", "PV2READ"))
  expect_equal(design$weight_col, "W_FSTUWT")
  expect_equal(design$rep_weight_cols, paste0("W_FSTURWT", 1:4))
  expect_equal(design$fay_k, 0.5)
  expect_equal(design$M, 2L)
  expect_equal(design$R, 4L)
  expect_equal(design$n, nrow(data))
  expect_equal(design$id_cols, "CNTSTUID")
  expect_equal(design$row_support$type, "id_cols")
  expect_equal(design$row_support$id_cols, "CNTSTUID")
  expect_equal(design$row_support$n, nrow(data))
  expect_equal(design$row_support$hash, design$row_support_hash)
  expect_equal(design$roles$outcome_placeholder, "OUTCOME")
  expect_equal(design$roles$helper, "pv_design")
  expect_true(design$roles$column_detection$pv_cols_detected)
  expect_true(design$roles$column_detection$rep_weight_cols_detected)
  expect_equal(design$roles$column_detection$pv_suffix, "READ")
  expect_equal(design$provenance$wrapper_function, "pv_design")
  expect_true(design$provenance$pv_cols_detected)
  expect_true(design$provenance$rep_weight_cols_detected)
  expect_match(design$row_support_hash, "^[0-9a-f]{8}$")
  expect_match(design$pv_value_hash, "^[0-9a-f]{8}$")
  expect_match(design$weight_design_hash, "^[0-9a-f]{8}$")
  expect_match(design$design_hash, "^[0-9a-f]{8}$")
  expect_invisible(pvstackr:::validate_pvstackr_design(design))
})

test_that("pv_design supports explicit columns and passes components to targets", {
  data <- design_helper_fixture_data()
  design <- pv_design(
    data = data,
    formula = OUTCOME ~ x,
    pv_cols = c("PV1READ", "PV2READ"),
    weight_col = "W_FSTUWT",
    rep_weight_cols = paste0("W_FSTURWT", 1:4),
    fay_k = 0.3,
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID",
    roles = list(analysis = "reading")
  )

  expect_false(design$roles$column_detection$pv_cols_detected)
  expect_false(design$roles$column_detection$rep_weight_cols_detected)
  expect_equal(design$roles$analysis, "reading")
  expect_equal(design$fay_k, 0.3)

  target <- pv_brr_target(
    data = design$data,
    formula = design$formula,
    pv_cols = design$pv_cols,
    weight_col = design$weight_col,
    rep_weight_cols = design$rep_weight_cols,
    fay_k = design$fay_k,
    id_cols = design$id_cols
  )

  expect_s3_class(target, "pvstackr_brr_target")
  expect_equal(target$M, design$M)
  expect_equal(target$R, design$R)
  expect_equal(target$fay_k, design$fay_k)
  expect_equal(target$formula_string, design$formula_string)
})

test_that("pv_design guides subject-suffixed PISA plausible values", {
  data <- design_helper_fixture_data()
  names(data)[names(data) == "PV1READ"] <- "PV1MATH"
  names(data)[names(data) == "PV2READ"] <- "PV2MATH"

  expect_error(
    pv_design(data, OUTCOME ~ x),
    "pv_suffix = \"MATH\""
  )
  design <- pv_design(
    data,
    OUTCOME ~ x,
    pv_suffix = "MATH",
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID"
  )
  expect_equal(design$pv_cols, c("PV1MATH", "PV2MATH"))

  mixed <- data
  mixed$PV1 <- mixed$PV1MATH
  mixed$PV2 <- mixed$PV2MATH
  expect_warning(
    mixed_design <- pv_design(mixed, OUTCOME ~ x, expected_M = 2L),
    "subject-suffixed"
  )
  expect_equal(mixed_design$pv_cols, c("PV1", "PV2"))
})

test_that("pv_design hashes row support, PV values, and weight design separately", {
  data <- design_helper_fixture_data()
  base <- pv_design(data, OUTCOME ~ x, pv_suffix = "READ", id_cols = "CNTSTUID")
  repeat_base <- pv_design(data, OUTCOME ~ x, pv_suffix = "READ", id_cols = "CNTSTUID")

  expect_equal(repeat_base$row_support_hash, base$row_support_hash)
  expect_equal(repeat_base$pv_value_hash, base$pv_value_hash)
  expect_equal(repeat_base$weight_design_hash, base$weight_design_hash)
  expect_equal(repeat_base$design_hash, base$design_hash)

  changed_id <- data
  changed_id$CNTSTUID[1] <- "new-id"
  row_changed <- pv_design(changed_id, OUTCOME ~ x, pv_suffix = "READ", id_cols = "CNTSTUID")
  expect_false(identical(row_changed$row_support_hash, base$row_support_hash))
  expect_equal(row_changed$pv_value_hash, base$pv_value_hash)
  expect_equal(row_changed$weight_design_hash, base$weight_design_hash)

  changed_pv <- data
  changed_pv$PV1READ[1] <- changed_pv$PV1READ[1] + 0.1
  pv_changed <- pv_design(changed_pv, OUTCOME ~ x, pv_suffix = "READ", id_cols = "CNTSTUID")
  expect_equal(pv_changed$row_support_hash, base$row_support_hash)
  expect_false(identical(pv_changed$pv_value_hash, base$pv_value_hash))
  expect_equal(pv_changed$weight_design_hash, base$weight_design_hash)

  fay_changed <- pv_design(data, OUTCOME ~ x, pv_suffix = "READ", fay_k = 0.3, id_cols = "CNTSTUID")
  expect_equal(fay_changed$row_support_hash, base$row_support_hash)
  expect_equal(fay_changed$pv_value_hash, base$pv_value_hash)
  expect_false(identical(fay_changed$weight_design_hash, base$weight_design_hash))
  expect_false(identical(fay_changed$design_hash, base$design_hash))

  no_ids <- pv_design(data, OUTCOME ~ x, pv_suffix = "READ")
  expect_equal(no_ids$row_support$type, "row_number")
  expect_equal(no_ids$row_support$id_cols, character())
  expect_equal(no_ids$row_support$hash, no_ids$row_support_hash)
})

test_that("pv_design print output is compact and stable", {
  design <- pv_design(
    design_helper_fixture_data(),
    OUTCOME ~ x,
    pv_suffix = "READ",
    id_cols = "CNTSTUID"
  )

  printed <- paste(capture.output(returned <- print(design)), collapse = "\n")
  expect_identical(returned, design)
  expect_match(printed, "pvstackr design", fixed = TRUE)
  expect_match(printed, "rows: 10", fixed = TRUE)
  expect_match(printed, "formula: OUTCOME ~ x", fixed = TRUE)
  expect_match(printed, "plausible values: 2", fixed = TRUE)
  expect_match(printed, "final weight: W_FSTUWT", fixed = TRUE)
  expect_match(printed, "replicate weights: 4", fixed = TRUE)
  expect_match(printed, "fay_k: 0.5", fixed = TRUE)
})

test_that("pv_design rejects malformed PISA-shaped declarations", {
  data <- design_helper_fixture_data()

  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "MATH"),
    "No plausible-value"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x),
    "pv_suffix = \"READ\""
  )
  expect_error(
    pv_design(data, OUTCOME ~ x),
    "suffix = \"READ\""
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", expected_M = 3L),
    "Expected 3"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", expected_R = 5L),
    "Expected 5"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_cols = c("PV1READ", "missing"), rep_weight_cols = paste0("W_FSTURWT", 1:4)),
    "not found"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_cols = c("PV1READ", "PV1READ"), rep_weight_cols = paste0("W_FSTURWT", 1:4)),
    "unique"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", weight_col = "W_FSTURWT1"),
    "distinct"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", weight_col = "missing"),
    "not found"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", rep_weight_cols = "W_FSTURWT1"),
    "at least two"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", id_cols = "school"),
    "unique rows"
  )

  gap <- data[, setdiff(names(data), "PV2READ")]
  gap$PV3READ <- data$PV2READ + 0.1
  expect_error(
    pv_design(gap, OUTCOME ~ x, pv_suffix = "READ"),
    "contiguous"
  )

  duplicated_numeric_suffix <- data
  duplicated_numeric_suffix$PV01READ <- duplicated_numeric_suffix$PV1READ + 0.1
  expect_error(
    pv_design(duplicated_numeric_suffix, OUTCOME ~ x, pv_suffix = "READ"),
    "Duplicate plausible-value"
  )

  duplicate_names <- data.frame(
    id = 1:2,
    PV1READ = 1:2,
    PV1READ = 2:3,
    W_FSTUWT = c(1, 1),
    W_FSTURWT1 = c(1, 1),
    W_FSTURWT2 = c(1, 1),
    check.names = FALSE
  )
  expect_error(
    pv_design(duplicate_names, OUTCOME ~ id, pv_suffix = "READ", id_cols = "id"),
    "column names must be unique"
  )

  expect_s3_class(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", fay_k = 0, id_cols = "CNTSTUID"),
    "pvstackr_design"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", fay_k = 1, id_cols = "CNTSTUID"),
    "fay_k"
  )
  expect_error(
    pv_design(data, OUTCOME ~ x, pv_suffix = "READ", fay_k = NA_real_, id_cols = "CNTSTUID"),
    "fay_k"
  )
})
