test_that("binding errors expose a stable typed registry", {
  registry <- pvstackr:::pv_binding_error_registry()

  expect_length(registry, 28L)
  expect_identical(unname(registry[["PV_BIND_E005"]]), "MANIFEST_HASH_STALE")
  expect_identical(unname(registry[["PV_BIND_E081"]]), "BINDING_RECOMPUTE_FAILED")

  error <- tryCatch(
    pvstackr:::pv_binding_abort(
      "PV_BIND_E081",
      "Canonicalization failed.",
      component = "canonicalizer"
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(error, "pvstackr_binding_error")
  expect_identical(error$code, "PV_BIND_E081")
  expect_identical(error$condition_name, "BINDING_RECOMPUTE_FAILED")
  expect_identical(error$component, "canonicalizer")
})

test_that("TLV unsigned lengths are explicit eight-byte big-endian integers", {
  expect_identical(pvstackr:::pv_binding_u64be(0), as.raw(c(0, 0, 0, 0, 0, 0, 0, 0)))
  expect_identical(pvstackr:::pv_binding_u64be(1), as.raw(c(0, 0, 0, 0, 0, 0, 0, 1)))
  expect_identical(pvstackr:::pv_binding_u64be(255), as.raw(c(0, 0, 0, 0, 0, 0, 0, 255)))
  expect_identical(pvstackr:::pv_binding_u64be(256), as.raw(c(0, 0, 0, 0, 0, 0, 1, 0)))

  expect_error(pvstackr:::pv_binding_u64be(-1), class = "pvstackr_binding_error")
  expect_error(pvstackr:::pv_binding_u64be(1.5), class = "pvstackr_binding_error")
})

test_that("SHA-256 provider and domain separation match frozen literals", {
  expect_identical(
    pvstackr:::pv_binding_sha256_raw(charToRaw("abc")),
    paste0(
      "sha256:",
      "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )
  )

  expect_identical(
    pvstackr:::pv_binding_hash_payload("abc", "manifest"),
    paste0(
      "sha256:",
      "86da25093e43a20341d8e8bee8097344c9176cd2336306e172672c45382dae4b"
    )
  )
  expect_false(identical(
    pvstackr:::pv_binding_hash_payload("abc", "manifest"),
    pvstackr:::pv_binding_hash_payload("abc", "row")
  ))
  expect_error(
    pvstackr:::pv_binding_hash_payload("abc", "not_a_domain"),
    class = "pvstackr_binding_error"
  )
  expect_error(
    pvstackr:::pv_binding_sha256_raw("abc"),
    class = "pvstackr_binding_error"
  )
})

test_that("canonical bytes preserve type, boundaries, names, and order", {
  hash <- function(x) pvstackr:::pv_binding_hash_payload(x, "manifest")

  expect_false(identical(hash(NULL), hash(list())))
  expect_false(identical(hash(""), hash(NA_character_)))
  expect_false(identical(hash(TRUE), hash(1L)))
  expect_false(identical(hash(1L), hash(1)))
  expect_false(identical(hash(c("ab", "c")), hash(c("a", "bc"))))
  expect_false(identical(hash(c(a = 1L, b = 2L)), hash(c(b = 1L, a = 2L))))
  expect_false(identical(hash(list(a = 1L, b = 2L)), hash(list(b = 2L, a = 1L))))

  matrix_a <- matrix(1:4, nrow = 2L, dimnames = list(c("r1", "r2"), c("a", "b")))
  matrix_b <- matrix(1:4, nrow = 4L)
  matrix_c <- matrix_a
  dimnames(matrix_c)[[2L]] <- c("a", "c")
  expect_false(identical(hash(matrix_a), hash(matrix_b)))
  expect_false(identical(hash(matrix_a), hash(matrix_c)))
})

test_that("canonical bytes normalize negative zero and ignore display options", {
  hash <- function(x) pvstackr:::pv_binding_hash_payload(x, "manifest")

  expect_identical(hash(-0), hash(0))
  expect_identical(hash(c(-0, 1.5)), hash(c(0, 1.5)))

  old_options <- options(digits = 3, width = 20)
  on.exit(options(old_options), add = TRUE)
  compact <- hash(list(beta = c(pi, sqrt(2)), label = "alpha"))
  options(digits = 15, width = 200)
  expanded <- hash(list(beta = c(pi, sqrt(2)), label = "alpha"))
  expect_identical(compact, expanded)
})

test_that("class dispatch binds factor, ordered, date, timezone, and rejects unknown objects", {
  hash <- function(x) pvstackr:::pv_binding_hash_payload(x, "manifest")

  factor_a <- factor(c("a", "b"), levels = c("a", "b"))
  factor_b <- factor(c("a", "b"), levels = c("b", "a"))
  expect_false(identical(hash(factor_a), hash(factor_b)))
  expect_false(identical(hash(factor_a), hash(ordered(factor_a))))
  expect_false(identical(hash(as.Date("2026-01-01")), hash(as.double(as.Date("2026-01-01")))))

  utc <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  ny <- utc
  attr(ny, "tzone") <- "America/New_York"
  expect_false(identical(hash(utc), hash(ny)))

  expect_error(hash(structure(1, class = "mystery")), class = "pvstackr_binding_error")
  expect_error(hash(data.frame(x = 1)), class = "pvstackr_binding_error")
  expect_error(hash(list(a = 1, a = 2)), class = "pvstackr_binding_error")
})

test_that("formula AST excludes whitespace, source text, and environment identity", {
  env_a <- new.env(parent = baseenv())
  env_b <- new.env(parent = baseenv())
  formula_a <- stats::as.formula("OUTCOME~x+z", env = env_a)
  formula_b <- stats::as.formula("  OUTCOME  ~  x + z  ", env = env_b)
  data <- data.frame(x = 1:3, z = 4:6)

  ast_a <- pvstackr:::pv_binding_formula_ast(formula_a, data)
  ast_b <- pvstackr:::pv_binding_formula_ast(formula_b, data)

  expect_identical(ast_a, ast_b)
  expect_identical(ast_a$data_variables, c("x", "z"))
  expect_identical(
    pvstackr:::pv_binding_hash_payload(ast_a, "formula"),
    pvstackr:::pv_binding_hash_payload(ast_b, "formula")
  )
  expect_false(any(vapply(ast_a, is.environment, logical(1))))
  expect_false(any(vapply(ast_a, is.language, logical(1))))
})

test_that("formula AST preserves named and positional call argument order", {
  ast_positional <- pvstackr:::pv_binding_formula_ast(OUTCOME ~ transform_x(x, 2))
  ast_named <- pvstackr:::pv_binding_formula_ast(OUTCOME ~ transform_x(x, degree = 2))
  ast_reordered <- pvstackr:::pv_binding_formula_ast(OUTCOME ~ transform_x(degree = 2, x))

  hash <- function(ast) pvstackr:::pv_binding_hash_payload(ast, "formula")
  expect_false(identical(hash(ast_positional), hash(ast_named)))
  expect_false(identical(hash(ast_named), hash(ast_reordered)))

  named_call <- ast_named$rhs
  expect_identical(named_call$argument_names, c("", "degree"))
  reordered_call <- ast_reordered$rhs
  expect_identical(reordered_call$argument_names, c("degree", ""))
})

test_that("formula AST binds term order, intercept, and transformations", {
  hash_formula <- function(formula) {
    pvstackr:::pv_binding_hash_payload(
      pvstackr:::pv_binding_formula_ast(formula),
      "formula"
    )
  }

  expect_false(identical(hash_formula(OUTCOME ~ x + z), hash_formula(OUTCOME ~ z + x)))
  expect_false(identical(hash_formula(OUTCOME ~ x), hash_formula(OUTCOME ~ 0 + x)))
  expect_false(identical(hash_formula(OUTCOME ~ x), hash_formula(OUTCOME ~ I(x^2))))
  expect_false(identical(hash_formula(OUTCOME ~ x + z), hash_formula(OUTCOME ~ x:z)))
})

test_that("formula data validation rejects external variables but not function names", {
  data <- data.frame(x = 1:3)

  expect_identical(
    pvstackr:::pv_binding_formula_rhs_info(OUTCOME ~ custom_transform(x), data)$data_variables,
    "x"
  )
  expect_identical(
    pvstackr:::pv_binding_formula_rhs_info(OUTCOME ~ stats::poly(x, 2), data)$data_variables,
    "x"
  )

  error <- tryCatch(
    pvstackr:::pv_binding_formula_ast(OUTCOME ~ x + external_scale, data),
    pvstackr_binding_error = identity
  )
  expect_s3_class(error, "pvstackr_binding_error")
  expect_identical(error$code, "PV_BIND_E081")
  expect_identical(error$component, "formula")
  expect_match(conditionMessage(error), "external_scale", fixed = TRUE)

  expect_error(
    pvstackr:::pv_binding_formula_ast(y ~ x, data),
    class = "pvstackr_binding_error"
  )
})

binding_projection_fixture <- function() {
  data.frame(
    x = c(-2, -1, 0, 1, 2, 3),
    z = c(0.5, 1, 1.5, 2, 2.5, 3),
    g = factor(
      c("private_alpha", "private_beta", "private_alpha", "private_beta", "private_alpha", "private_beta"),
      levels = c("private_alpha", "private_beta")
    ),
    unused = letters[1:6]
  )
}

test_that("predictor projection binds class, factor level order, and ordered status", {
  data <- binding_projection_fixture()
  base <- pvstackr:::pv_binding_predictor_projection(data, OUTCOME ~ x + g)

  integer_data <- data
  integer_data$x <- as.integer(integer_data$x)
  integer_projection <- pvstackr:::pv_binding_predictor_projection(integer_data, OUTCOME ~ x + g)
  expect_false(identical(base$schema_hash, integer_projection$schema_hash))
  expect_identical(integer_projection$columns[[1L]]$type_id, "integer")

  releveled <- data
  releveled$g <- factor(releveled$g, levels = rev(levels(releveled$g)))
  relevel_projection <- pvstackr:::pv_binding_predictor_projection(releveled, OUTCOME ~ x + g)
  expect_false(identical(
    base$columns[[2L]]$ordered_level_labels_hash,
    relevel_projection$columns[[2L]]$ordered_level_labels_hash
  ))

  ordered_data <- data
  ordered_data$g <- ordered(ordered_data$g)
  ordered_projection <- pvstackr:::pv_binding_predictor_projection(ordered_data, OUTCOME ~ x + g)
  expect_identical(ordered_projection$columns[[2L]]$type_id, "ordered_factor")
  expect_true(ordered_projection$columns[[2L]]$ordered)
  expect_false(identical(base$schema_hash, ordered_projection$schema_hash))
})

test_that("model projection binds transformations, interactions, and offsets", {
  data <- binding_projection_fixture()
  linear <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ x + z)
  transformed <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ I(x^2) + z)
  interaction <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ x * z)
  offset <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ x + z + offset(z))

  expect_false(identical(linear$values_hash, transformed$values_hash))
  expect_false(identical(linear$ordered_term_labels_hash, transformed$ordered_term_labels_hash))
  expect_gt(interaction$column_count, linear$column_count)
  expect_false(identical(linear$values_hash, interaction$values_hash))
  expect_false(linear$has_offset)
  expect_true(offset$has_offset)
  expect_false(identical(linear$offset_values_hash, offset$offset_values_hash))

  shifted_offset <- data
  shifted_offset$z <- shifted_offset$z + 10
  shifted <- pvstackr:::pv_binding_model_matrix_projection(
    shifted_offset,
    OUTCOME ~ x + offset(z)
  )
  original <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ x + offset(z))
  expect_false(identical(original$offset_values_hash, shifted$offset_values_hash))
})

test_that("global contrast changes alter privacy-preserving contrast and matrix hashes", {
  data <- binding_projection_fixture()
  old_options <- options("contrasts")
  on.exit(options(old_options), add = TRUE)

  options(contrasts = c("contr.treatment", "contr.poly"))
  treatment <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ x + g)
  options(contrasts = c("contr.sum", "contr.poly"))
  sum_contrast <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ x + g)

  expect_identical(treatment$contrast_count, 1L)
  expect_identical(sum_contrast$contrast_count, 1L)
  expect_false(identical(
    treatment$contrasts[[1L]]$generator_label_hash,
    sum_contrast$contrasts[[1L]]$generator_label_hash
  ))
  expect_false(identical(
    treatment$contrasts[[1L]]$values_hash,
    sum_contrast$contrasts[[1L]]$values_hash
  ))
  expect_false(identical(treatment$values_hash, sum_contrast$values_hash))
})

test_that("predictor mutation changes projections while unused-column mutation is invariant", {
  data <- binding_projection_fixture()
  predictor <- pvstackr:::pv_binding_predictor_projection(data, OUTCOME ~ x + g)
  model <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ x + g)

  changed <- data
  changed$x[[2L]] <- changed$x[[2L]] + 0.25
  changed_predictor <- pvstackr:::pv_binding_predictor_projection(changed, OUTCOME ~ x + g)
  changed_model <- pvstackr:::pv_binding_model_matrix_projection(changed, OUTCOME ~ x + g)
  expect_identical(predictor$schema_hash, changed_predictor$schema_hash)
  expect_false(identical(predictor$ordered_values_hash, changed_predictor$ordered_values_hash))
  expect_false(identical(model$values_hash, changed_model$values_hash))

  unused_changed <- data
  unused_changed$unused <- rev(unused_changed$unused)
  unused_changed$another_unused <- paste0("respondent_", seq_len(nrow(unused_changed)))
  expect_identical(
    predictor,
    pvstackr:::pv_binding_predictor_projection(unused_changed, OUTCOME ~ x + g)
  )
  expect_identical(
    model,
    pvstackr:::pv_binding_model_matrix_projection(unused_changed, OUTCOME ~ x + g)
  )
})

test_that("stored projections redact raw level, xlevel, and contrast labels", {
  data <- binding_projection_fixture()
  predictor <- pvstackr:::pv_binding_predictor_projection(data, OUTCOME ~ x + g)
  model <- pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ x + g)
  forbidden <- c("private_alpha", "private_beta", "contr.treatment", "contr.sum")

  expect_false(any(unlist(predictor, use.names = FALSE) %in% forbidden))
  expect_false(any(unlist(model$xlevels, use.names = FALSE) %in% forbidden))
  expect_false(any(unlist(model$contrasts, use.names = FALSE) %in% forbidden))
  expect_false(any(unlist(model, use.names = FALSE) %in% forbidden))
  expect_false("fe_colnames" %in% names(model))
  expect_identical(model$fe_colname_count, model$column_count)
  expect_match(model$ordered_fe_colnames_hash, "^sha256:[0-9a-f]{64}$")
  expect_identical(model$xlevels[[1L]]$level_count, 2L)
  expect_match(model$xlevels[[1L]]$ordered_level_labels_hash, "^sha256:[0-9a-f]{64}$")
  expect_identical(model$contrasts[[1L]]$dimensions, c(2L, 1L))
})

test_that("missing, nonfinite, unsupported, and invalid transformed predictors fail closed", {
  data <- binding_projection_fixture()

  missing <- data
  missing$x[[1L]] <- NA_real_
  expect_error(
    pvstackr:::pv_binding_predictor_projection(missing, OUTCOME ~ x),
    class = "pvstackr_binding_error"
  )

  nonfinite <- data
  nonfinite$x[[1L]] <- Inf
  expect_error(
    pvstackr:::pv_binding_model_matrix_projection(nonfinite, OUTCOME ~ x),
    class = "pvstackr_binding_error"
  )

  unsupported <- data
  unsupported$list_predictor <- I(lapply(seq_len(nrow(data)), function(index) index))
  unsupported_error <- tryCatch(
    pvstackr:::pv_binding_predictor_projection(unsupported, OUTCOME ~ list_predictor),
    pvstackr_binding_error = identity
  )
  expect_s3_class(unsupported_error, "pvstackr_binding_error")
  expect_identical(unsupported_error$code, "PV_BIND_E081")
  expect_identical(unsupported_error$component, "predictor")

  zero_data <- data
  zero_data$x[[3L]] <- 0
  transformed_error <- tryCatch(
    pvstackr:::pv_binding_model_matrix_projection(zero_data, OUTCOME ~ I(1 / x)),
    pvstackr_binding_error = identity
  )
  expect_s3_class(transformed_error, "pvstackr_binding_error")
  expect_identical(transformed_error$code, "PV_BIND_E081")
  expect_identical(transformed_error$component, "model_matrix")

  factor_missing <- data
  factor_missing$g[[1L]] <- NA
  expect_error(
    pvstackr:::pv_binding_predictor_projection(factor_missing, OUTCOME ~ g),
    class = "pvstackr_binding_error"
  )

  character_data <- data
  character_data$character_predictor <- as.character(character_data$g)
  expect_error(
    pvstackr:::pv_binding_predictor_projection(character_data, OUTCOME ~ character_predictor),
    "explicit factor",
    class = "pvstackr_binding_error"
  )
  character_data$character_predictor[[1L]] <- NA_character_
  expect_error(
    pvstackr:::pv_binding_predictor_projection(character_data, OUTCOME ~ character_predictor),
    class = "pvstackr_binding_error"
  )

  date_data <- data
  date_data$when <- structure(rep(1, nrow(data)), class = "Date")
  date_data$when[[2L]] <- Inf
  expect_error(
    pvstackr:::pv_binding_predictor_projection(date_data, OUTCOME ~ when),
    class = "pvstackr_binding_error"
  )

  time_data <- data
  time_data$when <- structure(
    rep(1, nrow(data)),
    class = c("POSIXct", "POSIXt"),
    tzone = "UTC"
  )
  time_data$when[[2L]] <- Inf
  expect_error(
    pvstackr:::pv_binding_predictor_projection(time_data, OUTCOME ~ when),
    class = "pvstackr_binding_error"
  )

  warning_error <- tryCatch(
    pvstackr:::pv_binding_model_matrix_projection(data, OUTCOME ~ log(x)),
    pvstackr_binding_error = identity
  )
  expect_s3_class(warning_error, "pvstackr_binding_error")
  expect_identical(warning_error$code, "PV_BIND_E081")
  expect_identical(warning_error$component, "model_matrix")
})

binding_component_fixture <- function() {
  data.frame(
    respondent_id = paste0("respondent_secret_", 1:6),
    PV1 = c(101.25, 102.5, 103.75, 104.25, 105.5, 106.75),
    PV2 = c(100.75, 102.25, 103.5, 104.75, 105.25, 107.0),
    x = c(-2, -1, 0, 1, 2, 3),
    g = factor(
      c("category_secret_a", "category_secret_b", "category_secret_a", "category_secret_b", "category_secret_a", "category_secret_b"),
      levels = c("category_secret_a", "category_secret_b")
    ),
    W = c(1.05, 1.15, 0.95, 1.25, 1.35, 0.85),
    RW1 = c(0.95, 1.25, 1.05, 1.15, 1.45, 0.75),
    RW2 = c(1.15, 1.05, 0.85, 1.35, 1.25, 0.95),
    unused = paste0("unused_secret_", 1:6)
  )
}

test_that("row component separates order from membership for IDs and bound-row fallback", {
  data <- binding_component_fixture()
  declared <- pvstackr:::pv_binding_row_projection(data, id_cols = "respondent_id")
  reordered <- pvstackr:::pv_binding_row_projection(data[6:1, ], id_cols = "respondent_id")

  expect_identical(declared$identity_mode, "declared_id")
  expect_false(identical(declared$row_order_hash, reordered$row_order_hash))
  expect_identical(declared$row_membership_hash, reordered$row_membership_hash)
  expect_false(any(grepl("respondent_secret_", unlist(declared), fixed = TRUE)))
  expect_match(
    pvstackr:::pv_binding_component_hash(declared, "row"),
    "^sha256:[0-9a-f]{64}$"
  )

  bound_cols <- c("PV1", "PV2", "x", "g", "W", "RW1", "RW2")
  fallback <- pvstackr:::pv_binding_row_projection(data, bound_cols = bound_cols)
  fallback_reordered <- pvstackr:::pv_binding_row_projection(
    data[6:1, ],
    bound_cols = bound_cols
  )
  expect_identical(fallback$identity_mode, "bound_row_digest")
  expect_false(identical(fallback$row_order_hash, fallback_reordered$row_order_hash))
  expect_identical(fallback$row_membership_hash, fallback_reordered$row_membership_hash)

  character_id_changed <- data
  character_id_changed$respondent_id[[1L]] <- "replacement_secret_id"
  changed <- pvstackr:::pv_binding_row_projection(character_id_changed, id_cols = "respondent_id")
  expect_false(identical(declared$row_membership_hash, changed$row_membership_hash))

  duplicate <- data
  duplicate$respondent_id[[2L]] <- duplicate$respondent_id[[1L]]
  expect_error(
    pvstackr:::pv_binding_row_projection(duplicate, id_cols = "respondent_id"),
    class = "pvstackr_binding_error"
  )
  expect_error(
    pvstackr:::pv_binding_row_projection(data),
    class = "pvstackr_binding_error"
  )
})

test_that("PV component binds declared order, values, and keyed membership", {
  data <- binding_component_fixture()
  pv <- pvstackr:::pv_binding_pv_projection(
    data,
    c("PV1", "PV2"),
    id_cols = "respondent_id"
  )
  reordered_columns <- pvstackr:::pv_binding_pv_projection(
    data,
    c("PV2", "PV1"),
    id_cols = "respondent_id"
  )
  expect_false(identical(
    pvstackr:::pv_binding_component_hash(pv, "pv"),
    pvstackr:::pv_binding_component_hash(reordered_columns, "pv")
  ))

  changed <- data
  changed$PV2[[3L]] <- changed$PV2[[3L]] + 0.125
  changed_pv <- pvstackr:::pv_binding_pv_projection(
    changed,
    c("PV1", "PV2"),
    id_cols = "respondent_id"
  )
  expect_identical(pv$per_column_hashes[["PV1"]], changed_pv$per_column_hashes[["PV1"]])
  expect_false(identical(pv$per_column_hashes[["PV2"]], changed_pv$per_column_hashes[["PV2"]]))
  expect_false(identical(pv$membership_values_hash, changed_pv$membership_values_hash))

  reordered_rows <- pvstackr:::pv_binding_pv_projection(
    data[6:1, ],
    c("PV1", "PV2"),
    id_cols = "respondent_id"
  )
  expect_false(identical(pv$ordered_values_hash, reordered_rows$ordered_values_hash))
  expect_identical(pv$membership_values_hash, reordered_rows$membership_values_hash)
  expect_false(any(unlist(pv, use.names = FALSE) %in% as.character(data$PV1)))

  expect_error(
    pvstackr:::pv_binding_pv_projection(data, "PV1", id_cols = "respondent_id"),
    class = "pvstackr_binding_error"
  )
  invalid <- data
  invalid$PV1[[1L]] <- NA_real_
  expect_error(
    pvstackr:::pv_binding_pv_projection(invalid, c("PV1", "PV2"), id_cols = "respondent_id"),
    class = "pvstackr_binding_error"
  )
})

test_that("weight component binds base, replicate order and values, Fay, and transform IDs", {
  data <- binding_component_fixture()
  weight <- pvstackr:::pv_binding_weight_projection(
    data,
    "W",
    c("RW1", "RW2"),
    id_cols = "respondent_id"
  )
  changed_base <- data
  changed_base$W[[1L]] <- changed_base$W[[1L]] + 0.1
  base_projection <- pvstackr:::pv_binding_weight_projection(
    changed_base,
    "W",
    c("RW1", "RW2"),
    id_cols = "respondent_id"
  )
  expect_false(identical(weight$base_values_hash, base_projection$base_values_hash))
  expect_identical(weight$replicate_values_hash, base_projection$replicate_values_hash)

  changed_rep <- data
  changed_rep$RW2[[2L]] <- changed_rep$RW2[[2L]] + 0.1
  rep_projection <- pvstackr:::pv_binding_weight_projection(
    changed_rep,
    "W",
    c("RW1", "RW2"),
    id_cols = "respondent_id"
  )
  expect_false(identical(weight$replicate_values_hash, rep_projection$replicate_values_hash))

  reversed <- pvstackr:::pv_binding_weight_projection(
    data,
    "W",
    c("RW2", "RW1"),
    id_cols = "respondent_id"
  )
  fay_changed <- pvstackr:::pv_binding_weight_projection(
    data,
    "W",
    c("RW1", "RW2"),
    fay_k = 0.4,
    id_cols = "respondent_id"
  )
  hash <- function(x) pvstackr:::pv_binding_component_hash(x, "weight")
  expect_false(identical(hash(weight), hash(reversed)))
  expect_false(identical(hash(weight), hash(fay_changed)))

  invalid <- data
  invalid$W[[1L]] <- 0
  expect_error(
    pvstackr:::pv_binding_weight_projection(
      invalid, "W", c("RW1", "RW2"), id_cols = "respondent_id"
    ),
    class = "pvstackr_binding_error"
  )
  expect_error(
    pvstackr:::pv_binding_weight_projection(
      data, "W", c("RW1", "RW1"), id_cols = "respondent_id"
    ),
    class = "pvstackr_binding_error"
  )
  expect_error(
    pvstackr:::pv_binding_weight_projection(
      data,
      "W",
      c("RW1", "RW2"),
      id_cols = "respondent_id",
      stack_transform_id = "unknown_transform"
    ),
    class = "pvstackr_binding_error"
  )
})

test_that("factor-contrast component changes with coding but retains no category labels", {
  data <- binding_component_fixture()
  old_options <- options("contrasts")
  on.exit(options(old_options), add = TRUE)

  options(contrasts = c("contr.treatment", "contr.poly"))
  treatment <- pvstackr:::pv_binding_factor_contrast_component(data, OUTCOME ~ x + g)
  options(contrasts = c("contr.sum", "contr.poly"))
  sum_contrast <- pvstackr:::pv_binding_factor_contrast_component(data, OUTCOME ~ x + g)

  hash <- function(x) pvstackr:::pv_binding_component_hash(x, "factor_contrast")
  expect_false(identical(hash(treatment), hash(sum_contrast)))
  expect_false(identical(treatment$resolved_hash, sum_contrast$resolved_hash))
  expect_false(any(unlist(treatment, use.names = FALSE) %in% c(
    "category_secret_a", "category_secret_b", "contr.treatment", "contr.sum"
  )))

  releveled <- data
  releveled$g <- factor(releveled$g, levels = rev(levels(releveled$g)))
  relevel_projection <- pvstackr:::pv_binding_factor_contrast_component(
    releveled,
    OUTCOME ~ x + g
  )
  expect_false(identical(hash(treatment), hash(relevel_projection)))

  missing <- data
  missing$g[[1L]] <- NA
  expect_error(
    pvstackr:::pv_binding_factor_contrast_component(missing, OUTCOME ~ x + g),
    class = "pvstackr_binding_error"
  )
})

test_that("affine estimand contrasts bind L, names, order, and coverage inheritance", {
  fe_names <- c("b_Intercept", "b_x")
  identity_contrast <- pvstackr:::pv_binding_estimand_contrast_projection(fe_names)
  expect_identical(identity_contrast$matrix_dimensions, c(2L, 2L))
  expect_true(identity_contrast$coverage_inheritance)
  expect_false("L" %in% names(identity_contrast))

  L <- matrix(c(1, 0, 1, -1), nrow = 2L, byrow = TRUE)
  dimnames(L) <- list(c("level", "difference"), fe_names)
  affine <- pvstackr:::pv_binding_estimand_contrast_projection(
    fe_names,
    L = L,
    contrast_id = "level_difference_v1"
  )
  expect_false(identical(identity_contrast$matrix_hash, affine$matrix_hash))
  expect_false(identical(
    pvstackr:::pv_binding_component_hash(identity_contrast, "estimand_contrast"),
    pvstackr:::pv_binding_component_hash(affine, "estimand_contrast")
  ))

  no_inheritance <- pvstackr:::pv_binding_estimand_contrast_projection(
    fe_names,
    coverage_inheritance = FALSE
  )
  expect_false(identical(
    pvstackr:::pv_binding_component_hash(identity_contrast, "estimand_contrast"),
    pvstackr:::pv_binding_component_hash(no_inheritance, "estimand_contrast")
  ))

  bad_L <- L
  colnames(bad_L) <- rev(fe_names)
  expect_error(
    pvstackr:::pv_binding_estimand_contrast_projection(
      fe_names, L = bad_L, contrast_id = "bad"
    ),
    class = "pvstackr_binding_error"
  )
  expect_error(
    pvstackr:::pv_binding_estimand_contrast_projection(
      fe_names, transformation_class = "nonlinear"
    ),
    class = "pvstackr_binding_error"
  )
})

test_that("family-link and estimand components bind all semantic metadata", {
  family <- pvstackr:::pv_binding_family_link_projection("gaussian", "identity")
  changed_link <- pvstackr:::pv_binding_family_link_projection("gaussian", "log")
  expect_false(identical(
    pvstackr:::pv_binding_component_hash(family, "family_link"),
    pvstackr:::pv_binding_component_hash(changed_link, "family_link")
  ))
  expect_error(
    pvstackr:::pv_binding_family_link_projection("", "identity"),
    class = "pvstackr_binding_error"
  )

  fe_names <- c("b_Intercept", "b_x")
  contrast <- pvstackr:::pv_binding_estimand_contrast_projection(fe_names)
  estimand <- pvstackr:::pv_binding_estimand_projection(
    fe_names,
    contrast,
    interval_role = "native_brr_fay_t"
  )
  changed_source <- pvstackr:::pv_binding_estimand_projection(
    fe_names,
    contrast,
    interval_role = "native_brr_fay_t",
    target_source = "external_model_rubin"
  )
  expect_false(identical(
    pvstackr:::pv_binding_component_hash(estimand, "estimand"),
    pvstackr:::pv_binding_component_hash(changed_source, "estimand")
  ))

  reversed_names <- rev(fe_names)
  reversed_contrast <- pvstackr:::pv_binding_estimand_contrast_projection(reversed_names)
  reversed_estimand <- pvstackr:::pv_binding_estimand_projection(
    reversed_names,
    reversed_contrast,
    interval_role = "native_brr_fay_t"
  )
  expect_false(identical(
    pvstackr:::pv_binding_component_hash(estimand, "estimand"),
    pvstackr:::pv_binding_component_hash(reversed_estimand, "estimand")
  ))

  no_inheritance <- pvstackr:::pv_binding_estimand_contrast_projection(
    fe_names,
    coverage_inheritance = FALSE
  )
  expect_error(
    pvstackr:::pv_binding_estimand_projection(
      fe_names,
      no_inheritance,
      interval_role = "native_brr_fay_t",
      coverage_claim_allowed = TRUE
    ),
    class = "pvstackr_binding_error"
  )
  expect_error(
    pvstackr:::pv_binding_estimand_projection(
      reversed_names,
      contrast,
      interval_role = "native_brr_fay_t"
    ),
    class = "pvstackr_binding_error"
  )
})

binding_manifest_metadata <- function(
  data,
  formula,
  target_source = "external_brr_fay_rubin",
  interval_role = "native_brr_fay_t",
  coverage_claim_allowed = FALSE
) {
  list(
    estimand_id = "brr_fay_fixed_effects",
    target_source = target_source,
    target_engine_id = "lm_wls_brr_fay_v1",
    parameter_scope = "fixed_effect",
    fe_names = pvstackr:::pv_stack_direct_fe_names(data, formula),
    interval_role = interval_role,
    coverage_claim_allowed = coverage_claim_allowed
  )
}

binding_manifest_build <- function(
  data = binding_component_fixture(),
  formula = OUTCOME ~ x + g,
  pv_cols = c("PV1", "PV2"),
  fay_k = 0.5,
  id_cols = "respondent_id",
  family_link = pvstackr:::pv_binding_family_link_projection("gaussian", "identity"),
  estimand_contrast = NULL,
  estimand_metadata = NULL
) {
  if (is.null(estimand_metadata)) {
    estimand_metadata <- binding_manifest_metadata(data, formula)
  }
  pvstackr:::pv_binding_manifest_build(
    data = data,
    formula = formula,
    pv_cols = pv_cols,
    weight_col = "W",
    rep_weight_cols = c("RW1", "RW2"),
    fay_k = fay_k,
    id_cols = id_cols,
    family_link = family_link,
    estimand_contrast = estimand_contrast,
    estimand_metadata = estimand_metadata
  )
}

test_that("manifest builder fixes exact root, component, and hash schema order", {
  manifest <- binding_manifest_build()

  expect_identical(names(manifest), pvstackr:::pv_binding_manifest_root_fields())
  expect_identical(names(manifest$components), pvstackr:::pv_binding_component_names())
  expect_identical(
    names(manifest$component_hashes),
    pvstackr:::pv_binding_component_hash_names()
  )
  field_registry <- pvstackr:::pv_binding_component_field_registry()
  for (name in pvstackr:::pv_binding_component_names()) {
    expect_identical(names(manifest$components[[name]]), field_registry[[name]])
  }
  expect_match(manifest$manifest_hash, "^sha256:[0-9a-f]{64}$")
  expect_invisible(pvstackr:::pv_binding_manifest_validate(manifest))
  expect_identical(manifest, binding_manifest_build())

  comparison <- pvstackr:::pv_binding_manifest_compare(manifest, manifest)
  expect_true(comparison$ok)
  expect_length(comparison$all_codes, 0L)
  expect_true(is.na(comparison$primary_code))
  expect_length(comparison$mismatched_components, 0L)
})

test_that("4a resolved model bundle evaluates RHS once and binds raw data, formula, and exact X", {
  data <- binding_component_fixture()
  evaluation_env <- new.env(parent = baseenv())
  evaluation_env$counter <- new.env(parent = emptyenv())
  evaluation_env$counter$n <- 0L
  evaluation_env$offset <- stats::offset
  evalq(
    stateful_transform <- function(x) {
      counter$n <- counter$n + 1L
      x + counter$n / 100
    },
    evaluation_env
  )
  formula <- stats::as.formula(
    "OUTCOME ~ stateful_transform(x) + offset(x)",
    env = evaluation_env
  )
  bundle <- pvstackr:::pv_binding_resolve_model_bundle(data, formula)

  expect_identical(evaluation_env$counter$n, 1L)
  expect_invisible(
    pvstackr:::pv_binding_validate_model_bundle(bundle, data, formula)
  )
  expect_identical(evaluation_env$counter$n, 1L)
  expect_s3_class(bundle$frame, "data.frame")
  expect_true(is.matrix(bundle$model_matrix))
  expect_identical(bundle$offset, unname(as.double(data$x)))
  expect_identical(bundle$formula_environment, evaluation_env)
  expect_equal(
    unname(bundle$model_matrix[, 2L]),
    data$x + 0.01,
    tolerance = 0
  )

  fe_names <- paste0("b_", colnames(bundle$model_matrix))
  fe_names[fe_names == "b_(Intercept)"] <- "b_Intercept"
  metadata <- list(
    estimand_id = "brr_fay_fixed_effects",
    target_source = "external_brr_fay_rubin",
    target_engine_id = "lm_wls_brr_fay_v1",
    parameter_scope = "fixed_effect",
    fe_names = fe_names,
    interval_role = "native_brr_fay_t",
    coverage_claim_allowed = FALSE
  )
  build <- function(input_data = data, input_formula = formula,
                    input_bundle = bundle, input_metadata = metadata) {
    pvstackr:::pv_binding_manifest_build(
      data = input_data,
      formula = input_formula,
      pv_cols = c("PV1", "PV2"),
      weight_col = "W",
      rep_weight_cols = c("RW1", "RW2"),
      fay_k = 0.5,
      id_cols = "respondent_id",
      family_link = pvstackr:::pv_binding_family_link_projection(
        "gaussian",
        "identity"
      ),
      estimand_contrast = NULL,
      estimand_metadata = input_metadata,
      model_bundle = input_bundle
    )
  }
  manifest <- build()
  expect_identical(evaluation_env$counter$n, 1L)
  expect_identical(
    manifest$components$model_matrix$values_hash,
    pvstackr:::pv_binding_numeric_matrix_hash(
      bundle$model_matrix,
      "model_matrix",
      "resolved_model_matrix_values"
    )
  )
  expect_identical(
    manifest$components$model_matrix$offset_values_hash,
    bundle$offset_values_hash
  )

  projection <- pvstackr:::pv_binding_model_matrix_projection(
    data,
    formula,
    model_bundle = bundle
  )
  factor_component <- pvstackr:::pv_binding_factor_contrast_component(
    data,
    formula,
    model_bundle = bundle,
    model_projection = projection
  )
  expect_identical(evaluation_env$counter$n, 1L)
  expect_identical(factor_component$resolved_specs, projection$contrasts)

  collect_fields <- function(x) {
    if (!is.list(x)) {
      return(character())
    }
    c(names(x), unlist(lapply(x, collect_fields), use.names = FALSE))
  }
  expect_false(any(c(
    "model_bundle", "frame", "resolved_contrasts", "bundle_hash",
    "formula_environment"
  ) %in% collect_fields(manifest)))

  changed_data <- data
  changed_data$x[[1L]] <- changed_data$x[[1L]] + 0.5
  changed_error <- tryCatch(
    build(input_data = changed_data),
    pvstackr_binding_error = identity
  )
  expect_identical(changed_error$code, "PV_BIND_E081")

  reordered_error <- tryCatch(
    build(input_data = data[6:1, ]),
    pvstackr_binding_error = identity
  )
  expect_identical(reordered_error$code, "PV_BIND_E081")

  changed_formula <- stats::as.formula(
    "OUTCOME ~ stateful_transform(x) + g + offset(x)",
    env = evaluation_env
  )
  formula_error <- tryCatch(
    build(input_formula = changed_formula),
    pvstackr_binding_error = identity
  )
  expect_identical(formula_error$code, "PV_BIND_E081")

  forged <- bundle
  forged$model_matrix[[1L, 2L]] <- forged$model_matrix[[1L, 2L]] + 1
  forged_error <- tryCatch(
    build(input_bundle = forged),
    pvstackr_binding_error = identity
  )
  expect_identical(forged_error$code, "PV_BIND_E081")

  forged_frame <- bundle
  forged_frame$frame[[1L]][[1L]] <- forged_frame$frame[[1L]][[1L]] + 1
  frame_error <- tryCatch(
    build(input_bundle = forged_frame),
    pvstackr_binding_error = identity
  )
  expect_identical(frame_error$code, "PV_BIND_E081")

  forged_offset <- bundle
  forged_offset$offset[[1L]] <- forged_offset$offset[[1L]] + 1
  offset_error <- tryCatch(
    build(input_bundle = forged_offset),
    pvstackr_binding_error = identity
  )
  expect_identical(offset_error$code, "PV_BIND_E081")
  expect_identical(evaluation_env$counter$n, 1L)
})

test_that("4a model bundle rejects identical AST reuse across formula environments", {
  data <- binding_component_fixture()
  env_a <- new.env(parent = baseenv())
  env_a$count <- 0L
  evalq(
    transform_x <- function(x) {
      count <<- count + 1L
      x
    },
    env_a
  )
  formula_a <- stats::as.formula("OUTCOME ~ transform_x(x)", env = env_a)
  bundle <- pvstackr:::pv_binding_resolve_model_bundle(data, formula_a)
  expect_identical(env_a$count, 1L)

  env_b <- new.env(parent = baseenv())
  env_b$count <- 0L
  evalq(
    transform_x <- function(x) {
      count <<- count + 1L
      2 * x
    },
    env_b
  )
  formula_b <- stats::as.formula("OUTCOME ~ transform_x(x)", env = env_b)
  error <- tryCatch(
    pvstackr:::pv_binding_validate_model_bundle(bundle, data, formula_b),
    pvstackr_binding_error = identity
  )
  expect_identical(error$code, "PV_BIND_E081")
  expect_identical(env_a$count, 1L)
  expect_identical(env_b$count, 0L)
})

test_that("4a default manifest builder resolves one model bundle", {
  data <- binding_component_fixture()
  evaluation_env <- new.env(parent = baseenv())
  evaluation_env$count <- 0L
  evalq(
    one_call_transform <- function(x) {
      count <<- count + 1L
      x
    },
    evaluation_env
  )
  formula <- stats::as.formula(
    "OUTCOME ~ one_call_transform(x)",
    env = evaluation_env
  )
  metadata <- list(
    estimand_id = "brr_fay_fixed_effects",
    target_source = "external_brr_fay_rubin",
    target_engine_id = "lm_wls_brr_fay_v1",
    parameter_scope = "fixed_effect",
    fe_names = c("b_Intercept", "b_one_call_transform(x)"),
    interval_role = "native_brr_fay_t",
    coverage_claim_allowed = FALSE
  )
  manifest <- pvstackr:::pv_binding_manifest_build(
    data = data,
    formula = formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = c("RW1", "RW2"),
    fay_k = 0.5,
    id_cols = "respondent_id",
    family_link = pvstackr:::pv_binding_family_link_projection(
      "gaussian",
      "identity"
    ),
    estimand_metadata = metadata
  )
  expect_invisible(pvstackr:::pv_binding_manifest_validate(manifest))
  expect_identical(evaluation_env$count, 1L)
})

test_that("4c binding proof is exact, minimal, and linked to identical manifests", {
  manifest <- binding_manifest_build()
  rebuilt <- binding_manifest_build()
  proof <- pvstackr:::pv_binding_proof_build(manifest, rebuilt)

  expect_identical(names(proof), pvstackr:::pv_binding_proof_fields())
  expect_identical(proof$contract_id, "pvstackr_data_binding_v1")
  expect_identical(proof$manifest_schema_version, "0.2.0")
  expect_identical(proof$target_manifest_hash, manifest$manifest_hash)
  expect_identical(proof$current_manifest_hash, manifest$manifest_hash)
  expect_identical(
    proof$verification_policy,
    "recomputed_current_inputs_exact_component_match_v1"
  )
  expect_invisible(
    pvstackr:::pv_binding_proof_validate(proof, target_manifest = manifest)
  )
  expect_invisible(
    pvstackr:::pv_binding_manifest_assert_match(manifest, rebuilt)
  )

  fields <- unlist(lapply(proof, function(value) {
    if (is.character(value)) value else character()
  }), use.names = FALSE)
  expect_false(any(grepl("respondent|PV1|RW1|PRIVATE", fields, fixed = FALSE)))
})

test_that("4c mismatch conditions preserve safe registry-ordered all_codes", {
  data <- binding_component_fixture()
  expected <- binding_manifest_build(data = data)
  changed <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  changed$PV1[[1L]] <- changed$PV1[[1L]] + 0.25
  changed$x[[2L]] <- changed$x[[2L]] + 0.5
  changed$W[[3L]] <- changed$W[[3L]] + 0.75
  observed <- binding_manifest_build(data = changed)
  comparison <- pvstackr:::pv_binding_manifest_compare(expected, observed)
  expect_false(comparison$ok)
  expect_gt(length(comparison$all_codes), 1L)
  registry <- names(pvstackr:::pv_binding_error_registry())
  expect_identical(
    comparison$all_codes,
    registry[registry %in% comparison$all_codes]
  )
  expect_identical(comparison$primary_code, comparison$all_codes[[1L]])

  error <- tryCatch(
    pvstackr:::pv_binding_manifest_assert_match(expected, observed),
    pvstackr_binding_error = identity
  )
  expect_s3_class(error, "pvstackr_binding_error")
  expect_identical(error$code, comparison$primary_code)
  expect_identical(error$all_codes, comparison$all_codes)
  expect_identical(error$component, "manifest")
  prefix_pattern <- sprintf(
    "^sha256:[0-9a-f]{%d}$",
    pvstackr:::pv_binding_hash_prefix_hex_length()
  )
  expect_match(error$expected_hash, prefix_pattern)
  expect_match(error$observed_hash, prefix_pattern)
  serialized <- rawToChar(serialize(error, NULL, ascii = TRUE))
  expect_false(grepl(expected$manifest_hash, serialized, fixed = TRUE))
  expect_false(grepl(observed$manifest_hash, serialized, fixed = TRUE))
})

test_that("4c binding proof rejects reordered, extra, stale, and unlinked fields", {
  manifest <- binding_manifest_build()
  proof <- pvstackr:::pv_binding_proof_build(manifest, manifest)

  invalid <- list(
    reordered = proof[rev(names(proof))],
    extra = c(proof, list(extra = TRUE)),
    wrong_contract = within(proof, contract_id <- "other_contract"),
    wrong_schema = within(proof, manifest_schema_version <- "9.9.9"),
    wrong_policy = within(proof, verification_policy <- "unchecked"),
    mismatched_hash = within(
      proof,
      current_manifest_hash <- paste0("sha256:", strrep("a", 64L))
    )
  )
  attributed <- proof
  attr(attributed, "marker") <- "not-canonical"
  invalid$attributed <- attributed
  attributed_hash <- proof
  attr(attributed_hash$target_manifest_hash, "marker") <- "not-canonical"
  invalid$attributed_hash <- attributed_hash

  for (candidate in invalid) {
    expect_error(
      pvstackr:::pv_binding_proof_validate(candidate, target_manifest = manifest),
      class = "pvstackr_binding_error"
    )
  }

  other <- manifest
  other$components$pv$per_column_hashes[[1L]] <- paste0(
    "sha256:", strrep("b", 64L)
  )
  other$component_hashes$pv_hash <- pvstackr:::pv_binding_component_hash(
    other$components$pv,
    "pv"
  )
  other$manifest_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(other),
    "manifest"
  )
  unlinked <- proof
  unlinked$target_manifest_hash <- other$manifest_hash
  unlinked$current_manifest_hash <- other$manifest_hash
  expect_error(
    pvstackr:::pv_binding_proof_validate(unlinked, target_manifest = manifest),
    class = "pvstackr_binding_error"
  )
})

test_that("row fallback derives the exact ordered bound-role union and excludes unused data", {
  data <- binding_component_fixture()
  manifest <- binding_manifest_build(data = data, id_cols = NULL)
  expect_identical(
    manifest$components$row$id_schema$identity_columns,
    c("PV1", "PV2", "x", "g", "W", "RW1", "RW2")
  )
  expect_false("respondent_id" %in% manifest$components$row$id_schema$identity_columns)
  expect_false("unused" %in% manifest$components$row$id_schema$identity_columns)

  unused_changed <- data
  unused_changed$unused <- rev(unused_changed$unused)
  unused_changed$another_unused <- paste0("another_secret_", seq_len(nrow(data)))
  expect_identical(manifest, binding_manifest_build(data = unused_changed, id_cols = NULL))

  mutators <- list(
    pv = function(x) { x$PV1[[1L]] <- x$PV1[[1L]] + 0.125; x },
    predictor = function(x) { x$x[[1L]] <- x$x[[1L]] + 0.125; x },
    base_weight = function(x) { x$W[[1L]] <- x$W[[1L]] + 0.125; x },
    replicate_weight = function(x) { x$RW1[[1L]] <- x$RW1[[1L]] + 0.125; x }
  )
  for (mutate in mutators) {
    changed <- binding_manifest_build(data = mutate(data), id_cols = NULL)
    expect_false(identical(
      manifest$components$row$row_membership_hash,
      changed$components$row$row_membership_hash
    ))
  }
})

test_that("manifest comparator classifies pure reorder without redundant value codes", {
  data <- binding_component_fixture()
  expected <- binding_manifest_build(data = data)
  observed <- binding_manifest_build(data = data[6:1, ])
  comparison <- pvstackr:::pv_binding_manifest_compare(expected, observed)

  expect_false(comparison$ok)
  expect_identical(comparison$all_codes, "PV_BIND_E012")
  expect_identical(comparison$primary_code, "PV_BIND_E012")
  expect_true("row" %in% comparison$mismatched_components)
  expect_true("pv" %in% comparison$mismatched_components)
  expect_true("predictor" %in% comparison$mismatched_components)
})

test_that("manifest comparator classifies PV, predictor, and weight mutations", {
  data <- binding_component_fixture()
  expected <- binding_manifest_build(data = data)

  pv_changed <- data
  pv_changed$PV1[[1L]] <- pv_changed$PV1[[1L]] + 0.125
  pv_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(data = pv_changed)
  )
  expect_identical(pv_result$primary_code, "PV_BIND_E021")
  expect_identical(pv_result$all_codes, "PV_BIND_E021")

  pv_order_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(data = data, pv_cols = c("PV2", "PV1"))
  )
  expect_identical(pv_order_result$primary_code, "PV_BIND_E020")
  expect_true("PV_BIND_E021" %in% pv_order_result$all_codes)

  predictor_changed <- data
  predictor_changed$x[[2L]] <- predictor_changed$x[[2L]] + 0.125
  predictor_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(data = predictor_changed)
  )
  expect_identical(predictor_result$primary_code, "PV_BIND_E031")
  expect_true("PV_BIND_E042" %in% predictor_result$all_codes)

  weight_changed <- data
  weight_changed$RW2[[3L]] <- weight_changed$RW2[[3L]] + 0.125
  weight_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(data = weight_changed)
  )
  expect_identical(weight_result$primary_code, "PV_BIND_E051")

  fay_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(data = data, fay_k = 0.4)
  )
  expect_identical(fay_result$primary_code, "PV_BIND_E053")
})

test_that("manifest comparator classifies formula and factor-contrast mutations", {
  data <- binding_component_fixture()
  expected <- binding_manifest_build(data = data)
  reordered_formula <- OUTCOME ~ g + x
  formula_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(data = data, formula = reordered_formula)
  )
  expect_identical(formula_result$primary_code, "PV_BIND_E040")
  expect_true("PV_BIND_E041" %in% formula_result$all_codes)
  expect_true("PV_BIND_E071" %in% formula_result$all_codes)
  expect_true("PV_BIND_E072" %in% formula_result$all_codes)

  old_options <- options("contrasts")
  on.exit(options(old_options), add = TRUE)
  options(contrasts = c("contr.treatment", "contr.poly"))
  treatment <- binding_manifest_build(data = data)
  options(contrasts = c("contr.sum", "contr.poly"))
  sum_contrast <- binding_manifest_build(data = data)
  contrast_result <- pvstackr:::pv_binding_manifest_compare(treatment, sum_contrast)
  expect_true("PV_BIND_E043" %in% contrast_result$all_codes)
  expect_true("factor_contrast" %in% contrast_result$mismatched_components)
})

test_that("manifest comparator classifies family and estimand mutations", {
  data <- binding_component_fixture()
  expected <- binding_manifest_build(data = data)

  changed_family <- pvstackr:::pv_binding_family_link_projection("gaussian", "log")
  family_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(data = data, family_link = changed_family)
  )
  expect_identical(family_result$primary_code, "PV_BIND_E061")

  changed_metadata <- binding_manifest_metadata(
    data,
    OUTCOME ~ x + g,
    target_source = "external_model_rubin"
  )
  estimand_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(data = data, estimand_metadata = changed_metadata)
  )
  expect_identical(estimand_result$primary_code, "PV_BIND_E070")

  fe_names <- changed_metadata$fe_names
  L <- diag(length(fe_names))
  L[[2L, 2L]] <- 2
  dimnames(L) <- list(fe_names, fe_names)
  changed_contrast <- pvstackr:::pv_binding_estimand_contrast_projection(
    fe_names,
    L = L,
    contrast_id = "scaled_slope_v1"
  )
  contrast_result <- pvstackr:::pv_binding_manifest_compare(
    expected,
    binding_manifest_build(
      data = data,
      estimand_metadata = binding_manifest_metadata(data, OUTCOME ~ x + g),
      estimand_contrast = changed_contrast
    )
  )
  expect_identical(contrast_result$primary_code, "PV_BIND_E072")
})

test_that("manifest validator rejects root, self-hash, and component-hash tampering", {
  manifest <- binding_manifest_build()

  component_tamper <- manifest
  component_tamper$components$pv$M <- 99L
  component_error <- tryCatch(
    pvstackr:::pv_binding_manifest_validate(component_tamper),
    pvstackr_binding_error = identity
  )
  expect_identical(component_error$code, "PV_BIND_E005")

  self_tamper <- manifest
  self_tamper$manifest_hash <- paste0("sha256:", strrep("0", 64L))
  self_error <- tryCatch(
    pvstackr:::pv_binding_manifest_validate(self_tamper),
    pvstackr_binding_error = identity
  )
  expect_identical(self_error$code, "PV_BIND_E005")

  stored_hash_tamper <- manifest
  stored_hash_tamper$component_hashes$pv_hash <- paste0("sha256:", strrep("1", 64L))
  stored_hash_tamper$manifest_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(stored_hash_tamper),
    "manifest"
  )
  stored_error <- tryCatch(
    pvstackr:::pv_binding_manifest_validate(stored_hash_tamper),
    pvstackr_binding_error = identity
  )
  expect_identical(stored_error$code, "PV_BIND_E005")

  extra_root <- manifest
  extra_root$unknown <- TRUE
  root_error <- tryCatch(
    pvstackr:::pv_binding_manifest_validate(extra_root),
    pvstackr_binding_error = identity
  )
  expect_identical(root_error$code, "PV_BIND_E002")
})

test_that("manifest self-hash excludes legacy, timestamp, and registered migration warnings only", {
  manifest <- binding_manifest_build()
  migrated_a <- c(
    manifest,
    list(
      legacy_hashes = list(
        algorithm_id = pvstackr:::pv_binding_legacy_hash_algorithm_id(),
        design_hash = "1234abcd"
      ),
      created_at = "2026-01-01T00:00:00Z",
      migration = list(
        migration_from_schema = "0.1.0",
        migration_function = "pv_binding_revalidate_brr_target_v1",
        binding_revalidated = TRUE,
        inspection_only = FALSE,
        warnings = "PV_BIND_MIGRATION_LEGACY_HASH_ONLY"
      )
    )
  )
  migrated_a$manifest_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(migrated_a),
    "manifest"
  )
  expect_invisible(pvstackr:::pv_binding_manifest_validate(migrated_a))

  migrated_b <- migrated_a
  migrated_b$legacy_hashes$design_hash <- "ffffffff"
  migrated_b$created_at <- "2030-12-31T23:59:59Z"
  migrated_b$migration$warnings <- "PV_BIND_MIGRATION_SOURCE_METADATA_DROPPED"
  expect_identical(migrated_a$manifest_hash, pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(migrated_b),
    "manifest"
  ))

  migrated_c <- migrated_b
  migrated_c$migration$migration_from_schema <- "different_schema"
  expect_false(identical(migrated_a$manifest_hash, pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(migrated_c),
    "manifest"
  )))
})

test_that("exact estimand contrast and semantic metadata validators reject forgeries", {
  fe_names <- c("b_Intercept", "b_x")
  contrast <- pvstackr:::pv_binding_estimand_contrast_projection(fe_names)

  bad_dimensions <- contrast
  bad_dimensions$matrix_dimensions <- c(3L, 2L)
  expect_error(
    pvstackr:::pv_binding_validate_estimand_contrast_component(bad_dimensions),
    class = "pvstackr_binding_error"
  )

  bad_hash <- contrast
  bad_hash$matrix_hash <- "sha256:bad"
  expect_error(
    pvstackr:::pv_binding_validate_estimand_contrast_component(bad_hash),
    class = "pvstackr_binding_error"
  )

  bad_output <- contrast
  bad_output$output_count <- 3L
  expect_error(
    pvstackr:::pv_binding_validate_estimand_contrast_component(bad_output),
    class = "pvstackr_binding_error"
  )

  bad_class <- contrast
  bad_class$transformation_class <- "nonlinear"
  expect_error(
    pvstackr:::pv_binding_validate_estimand_contrast_component(bad_class),
    class = "pvstackr_binding_error"
  )

  expect_error(
    pvstackr:::pv_binding_family_link_projection(" gaussian", "identity"),
    class = "pvstackr_binding_error"
  )
  metadata <- binding_manifest_metadata(binding_component_fixture(), OUTCOME ~ x + g)
  metadata$interval_role <- " "
  expect_error(
    binding_manifest_build(estimand_metadata = metadata),
    class = "pvstackr_binding_error"
  )
})

binding_target_content_fixture <- function(
  df_method = "classic",
  df_complete = NULL,
  conf_level = 0.95
) {
  data <- binding_component_fixture()
  formula <- OUTCOME ~ x
  target <- pv_brr_target(
    data = data,
    formula = formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = c("RW1", "RW2"),
    fay_k = 0.5,
    id_cols = "respondent_id",
    conf_level = conf_level,
    df_method = df_method,
    df_complete = df_complete
  )
  estimand_metadata <- binding_manifest_metadata(
    data,
    formula,
    target_source = target$target_source,
    interval_role = target$interval_role,
    coverage_claim_allowed = target$coverage_claim_allowed
  )
  manifest <- binding_manifest_build(
    data = data,
    formula = formula,
    estimand_metadata = estimand_metadata
  )
  manifest <- c(
    manifest,
    list(model_bundle_hash = target$binding_manifest$model_bundle_hash)
  )
  manifest$manifest_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(manifest),
    "manifest"
  )
  pvstackr:::pv_binding_manifest_validate(manifest)
  stored_derived <- target[pvstackr:::pv_binding_target_derived_fields()]
  content <- pvstackr:::pv_binding_target_content_build(
    per_pv = target$per_pv,
    M = target$M,
    R = target$R,
    fay_k = target$fay_k,
    df_method = target$df_method,
    df_complete = if (identical(df_method, "classic")) NULL else target$df_complete,
    conf_level = conf_level,
    target_source = target$target_source,
    target_engine_id = "lm_wls_brr_fay_v1",
    target_policy = target$policy,
    manifest_hash = manifest$manifest_hash,
    stored_derived = stored_derived
  )
  pvstackr:::pv_binding_target_manifest_validate(content, manifest)
  list(
    data = data,
    manifest = manifest,
    target = target,
    stored_derived = stored_derived,
    content = content
  )
}

binding_schema01_target <- function(target) {
  legacy_fields <- setdiff(
    pvstackr:::pv_brr_target_v02_fields(),
    c("conf_level", "binding_manifest", "target_content")
  )
  legacy <- target[legacy_fields]
  if (identical(legacy$df_method, "classic")) {
    legacy$df_complete <- stats::setNames(
      rep(NA_real_, length(legacy$fe_names)),
      legacy$fe_names
    )
  }
  legacy$design_hash <- target$binding_manifest$legacy_hashes$design_hash
  legacy$target_hash <- target$binding_manifest$legacy_hashes$target_hash
  legacy$schema_version <- "0.1.0"
  legacy$provenance <- list(
    function_name = "pv_brr_target",
    assembled_at = "2026-07-12 12:00:00",
    package = "pvstackr"
  )
  class(legacy) <- c("pvstackr_brr_target", "list")
  pvstackr:::validate_pvstackr_brr_target(legacy)
  legacy
}

test_that("target content builder stores exact primitives and one canonical derived representation", {
  fixture <- binding_target_content_fixture()
  content <- fixture$content

  expect_identical(names(content), pvstackr:::pv_binding_target_content_fields())
  expect_identical(
    names(content$interval_policy),
    pvstackr:::pv_binding_target_interval_fields()
  )
  expect_identical(
    names(content$target_policy),
    pvstackr:::pv_binding_target_policy_fields()
  )
  expect_identical(
    names(content$derived),
    pvstackr:::pv_binding_target_derived_fields()
  )
  expect_true(all(vapply(
    content$per_pv,
    function(item) identical(names(item), pvstackr:::pv_binding_target_per_pv_fields()),
    logical(1)
  )))
  expect_match(content$target_content_hash, "^sha256:[0-9a-f]{64}$")
  expect_invisible(pvstackr:::pv_binding_target_content_validate(content))
  expect_identical(content, binding_target_content_fixture()$content)

  derived_names <- pvstackr:::pv_binding_target_derived_fields()
  expect_length(intersect(derived_names, names(content)), 0L)
  expect_identical(content$numeric_comparison_policy, list(
    policy_id = "absolute_relative_v1",
    absolute_tolerance = 1e-12,
    relative_tolerance = 1e-10,
    infinity_policy = "positive_only_for_df"
  ))

  swapped <- pvstackr:::pv_binding_target_content_build(
    per_pv = rev(fixture$target$per_pv),
    M = fixture$target$M,
    R = fixture$target$R,
    fay_k = fixture$target$fay_k,
    df_method = fixture$target$df_method,
    df_complete = NULL,
    conf_level = 0.95,
    target_source = fixture$target$target_source,
    target_engine_id = "lm_wls_brr_fay_v1",
    target_policy = fixture$target$policy,
    manifest_hash = fixture$manifest$manifest_hash
  )
  expect_false(identical(content$target_content_hash, swapped$target_content_hash))
  expect_identical(vapply(swapped$per_pv, `[[`, character(1), "pv_col"), c("PV2", "PV1"))
  swapped_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(swapped, fixture$manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(swapped_error$code, "PV_BIND_E020")

  changed_fay_per_pv <- fixture$target$per_pv
  changed_fay_per_pv <- lapply(changed_fay_per_pv, function(item) {
    item$fay_k <- 0.4
    item
  })
  changed_fay <- pvstackr:::pv_binding_target_content_build(
    per_pv = changed_fay_per_pv,
    M = fixture$target$M,
    R = fixture$target$R,
    fay_k = 0.4,
    df_method = fixture$target$df_method,
    df_complete = NULL,
    conf_level = 0.90,
    target_source = fixture$target$target_source,
    target_engine_id = "lm_wls_brr_fay_v1",
    target_policy = fixture$target$policy,
    manifest_hash = fixture$manifest$manifest_hash
  )
  expect_false(identical(content$target_content_hash, changed_fay$target_content_hash))
  expect_identical(changed_fay$fay_k, 0.4)
  expect_identical(changed_fay$interval_policy$conf_level, 0.90)
  fay_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(changed_fay, fixture$manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(fay_error$code, "PV_BIND_E053")
})

test_that("stored derived copies are checked by explicit absolute-relative tolerance then canonicalized", {
  fixture <- binding_target_content_fixture()
  stored_close <- fixture$stored_derived
  stored_close$beta[[1L]] <- stored_close$beta[[1L]] + 1e-11
  rebuilt <- pvstackr:::pv_binding_target_content_build(
    per_pv = fixture$target$per_pv,
    M = fixture$target$M,
    R = fixture$target$R,
    fay_k = fixture$target$fay_k,
    df_method = fixture$target$df_method,
    df_complete = NULL,
    conf_level = 0.95,
    target_source = fixture$target$target_source,
    target_engine_id = "lm_wls_brr_fay_v1",
    target_policy = fixture$target$policy,
    manifest_hash = fixture$manifest$manifest_hash,
    stored_derived = stored_close
  )
  expect_identical(rebuilt, fixture$content)

  stored_bad <- fixture$stored_derived
  stored_bad$beta[[1L]] <- stored_bad$beta[[1L]] + 1e-3
  error <- tryCatch(
    pvstackr:::pv_binding_target_content_build(
      per_pv = fixture$target$per_pv,
      M = fixture$target$M,
      R = fixture$target$R,
      fay_k = fixture$target$fay_k,
      df_method = fixture$target$df_method,
      df_complete = NULL,
      conf_level = 0.95,
      target_source = fixture$target$target_source,
      target_engine_id = "lm_wls_brr_fay_v1",
      target_policy = fixture$target$policy,
      manifest_hash = fixture$manifest$manifest_hash,
      stored_derived = stored_bad
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(error$code, "PV_BIND_E090")
})

test_that("target content supports positive infinite df and Barnard-Rubin infinite complete df", {
  fixture <- binding_target_content_fixture()
  one_pv <- fixture$target$per_pv[1L]
  classic_policy <- pvstackr:::pv_binding_target_policy_projection(
    df_method = "classic",
    interval_role = "descriptive_classic_rubin",
    coverage_claim_allowed = FALSE
  )
  one_content <- pvstackr:::pv_binding_target_content_build(
    per_pv = one_pv,
    M = 1L,
    R = fixture$target$R,
    fay_k = fixture$target$fay_k,
    df_method = "classic",
    df_complete = NULL,
    conf_level = 0.95,
    target_source = "external_brr_fay_rubin",
    target_engine_id = "lm_wls_brr_fay_v1",
    target_policy = classic_policy,
    manifest_hash = fixture$manifest$manifest_hash
  )
  expect_true(all(is.infinite(one_content$derived$df)))
  expect_true(all(one_content$derived$df > 0))
  expect_invisible(pvstackr:::pv_binding_target_content_validate(one_content))
  expect_error(
    pvstackr:::pv_binding_target_manifest_validate(one_content, fixture$manifest),
    class = "pvstackr_binding_error"
  )

  br_fixture <- binding_target_content_fixture(
    df_method = "barnard_rubin",
    df_complete = Inf
  )
  expect_true(all(is.infinite(br_fixture$content$interval_policy$df_complete)))
  expect_true(all(br_fixture$content$derived$df > 0))
  expect_invisible(pvstackr:::pv_binding_target_content_validate(br_fixture$content))
})

test_that("raw per-PV NA/nonfinite values and improper df policies fail with E090", {
  fixture <- binding_target_content_fixture()
  build_with <- function(per_pv = fixture$target$per_pv, df_method = "classic",
                         df_complete = NULL, policy = fixture$target$policy) {
    pvstackr:::pv_binding_target_content_build(
      per_pv = per_pv,
      M = length(per_pv),
      R = fixture$target$R,
      fay_k = fixture$target$fay_k,
      df_method = df_method,
      df_complete = df_complete,
      conf_level = 0.95,
      target_source = fixture$target$target_source,
      target_engine_id = "lm_wls_brr_fay_v1",
      target_policy = policy,
      manifest_hash = fixture$manifest$manifest_hash
    )
  }

  bad_beta <- fixture$target$per_pv
  bad_beta[[1L]]$beta[[1L]] <- NA_real_
  beta_error <- tryCatch(build_with(bad_beta), pvstackr_binding_error = identity)
  expect_identical(beta_error$code, "PV_BIND_E090")

  secret_value <- 987654.321
  bad_beta[[1L]]$beta[[1L]] <- secret_value
  bad_beta[[1L]]$beta[[2L]] <- Inf
  inf_error <- tryCatch(build_with(bad_beta), pvstackr_binding_error = identity)
  expect_identical(inf_error$code, "PV_BIND_E090")
  expect_false(grepl(as.character(secret_value), conditionMessage(inf_error), fixed = TRUE))

  bad_U <- fixture$target$per_pv
  bad_U[[1L]]$U[[1L, 1L]] <- Inf
  U_error <- tryCatch(build_with(bad_U), pvstackr_binding_error = identity)
  expect_identical(U_error$code, "PV_BIND_E090")

  expect_error(build_with(df_complete = 20), class = "pvstackr_binding_error")
  br_policy <- pvstackr:::pv_binding_target_policy_projection(
    df_method = "barnard_rubin",
    interval_role = "coverage_barnard_rubin",
    coverage_claim_allowed = TRUE
  )
  for (bad_df in list(NA_real_, 0, -Inf)) {
    error <- tryCatch(
      build_with(
        df_method = "barnard_rubin",
        df_complete = bad_df,
        policy = br_policy
      ),
      pvstackr_binding_error = identity
    )
    expect_identical(error$code, "PV_BIND_E090")
  }
})

test_that("target content self-verification rejects derived, primitive, policy, and manifest tampering", {
  fixture <- binding_target_content_fixture()
  content <- fixture$content

  mutate_and_expect_e090 <- function(mutate, recompute_hash = FALSE) {
    changed <- mutate(content)
    if (recompute_hash) {
      changed$target_content_hash <- pvstackr:::pv_binding_hash_payload(
        pvstackr:::pv_binding_target_content_hash_payload(changed),
        "target_content"
      )
    }
    error <- tryCatch(
      pvstackr:::pv_binding_target_content_validate(changed),
      pvstackr_binding_error = identity
    )
    expect_identical(error$code, "PV_BIND_E090")
  }

  mutate_and_expect_e090(function(x) { x$derived$beta[[1L]] <- x$derived$beta[[1L]] + 1; x })
  mutate_and_expect_e090(function(x) { x$derived$U_bar[[1L, 1L]] <- x$derived$U_bar[[1L, 1L]] + 1; x })
  mutate_and_expect_e090(function(x) { x$derived$T_MI[[1L, 1L]] <- x$derived$T_MI[[1L, 1L]] + 1; x })
  mutate_and_expect_e090(function(x) { x$derived$df[[1L]] <- x$derived$df[[1L]] + 1; x })
  mutate_and_expect_e090(function(x) { x$target_policy$target_repair <- "allowed"; x })
  mutate_and_expect_e090(function(x) { x$manifest_hash <- paste0("sha256:", strrep("0", 64)); x })

  mutate_and_expect_e090(
    function(x) {
      x$per_pv[[1L]]$U[[1L, 1L]] <- x$per_pv[[1L]]$U[[1L, 1L]] + 0.01
      x
    },
    recompute_hash = TRUE
  )
  mutate_and_expect_e090(
    function(x) {
      x$derived$beta[[1L]] <- x$derived$beta[[1L]] + 0.01
      x
    },
    recompute_hash = TRUE
  )
})

test_that("target content exact validators reject extra and reordered fields", {
  content <- binding_target_content_fixture()$content

  extra <- content
  extra$unexpected <- TRUE
  expect_error(
    pvstackr:::pv_binding_target_content_validate(extra),
    class = "pvstackr_binding_error"
  )

  reordered <- content[rev(names(content))]
  expect_error(
    pvstackr:::pv_binding_target_content_validate(reordered),
    class = "pvstackr_binding_error"
  )

  bad_per_pv <- content
  bad_per_pv$per_pv[[1L]] <- bad_per_pv$per_pv[[1L]][rev(names(bad_per_pv$per_pv[[1L]]))]
  bad_per_pv$target_content_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_target_content_hash_payload(bad_per_pv),
    "target_content"
  )
  error <- tryCatch(
    pvstackr:::pv_binding_target_content_validate(bad_per_pv),
    pvstackr_binding_error = identity
  )
  expect_identical(error$code, "PV_BIND_E090")
})

test_that("3c remediation cross-link rejects refreshed manifest, PV, R, and Fay mismatches", {
  fixture <- binding_target_content_fixture()
  content <- fixture$content

  refresh <- function(x) {
    x$target_content_hash <- pvstackr:::pv_binding_hash_payload(
      pvstackr:::pv_binding_target_content_hash_payload(x),
      "target_content"
    )
    expect_invisible(pvstackr:::pv_binding_target_content_validate(x))
    x
  }

  wrong_manifest <- content
  wrong_manifest$manifest_hash <- paste0("sha256:", strrep("a", 64L))
  wrong_manifest <- refresh(wrong_manifest)
  manifest_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(wrong_manifest, fixture$manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(manifest_error$code, "PV_BIND_E090")

  wrong_pv <- content
  wrong_pv$per_pv[[1L]]$pv_col <- "PV_WRONG"
  wrong_pv <- refresh(wrong_pv)
  pv_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(wrong_pv, fixture$manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(pv_error$code, "PV_BIND_E020")

  wrong_R <- content
  wrong_R$R <- wrong_R$R + 1L
  wrong_R <- refresh(wrong_R)
  R_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(wrong_R, fixture$manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(R_error$code, "PV_BIND_E053")

  wrong_fay <- content
  wrong_fay$fay_k <- 0.4
  wrong_fay <- refresh(wrong_fay)
  fay_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(wrong_fay, fixture$manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(fay_error$code, "PV_BIND_E053")
})

test_that("3c remediation cross-link enforces FE, family-link, source-engine, and interval semantics", {
  fixture <- binding_target_content_fixture()
  data <- fixture$data
  formula <- OUTCOME ~ x
  relink <- function(manifest) {
    content <- fixture$content
    content$manifest_hash <- manifest$manifest_hash
    content$target_content_hash <- pvstackr:::pv_binding_hash_payload(
      pvstackr:::pv_binding_target_content_hash_payload(content),
      "target_content"
    )
    expect_invisible(pvstackr:::pv_binding_target_content_validate(content))
    content
  }
  matching_metadata <- binding_manifest_metadata(
    data,
    formula,
    target_source = fixture$content$target_source,
    interval_role = fixture$content$interval_policy$interval_role,
    coverage_claim_allowed = fixture$content$interval_policy$coverage_claim_allowed
  )

  link_manifest <- binding_manifest_build(
    data = data,
    formula = formula,
    family_link = pvstackr:::pv_binding_family_link_projection("gaussian", "log"),
    estimand_metadata = matching_metadata
  )
  link_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(relink(link_manifest), link_manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(link_error$code, "PV_BIND_E061")

  support_manifest <- binding_manifest_build(
    data = data,
    formula = formula,
    family_link = pvstackr:::pv_binding_family_link_projection(
      "gaussian", "identity", "positive", "estimated"
    ),
    estimand_metadata = matching_metadata
  )
  support_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(
      relink(support_manifest),
      support_manifest
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(support_error$code, "PV_BIND_E060")

  dispersion_manifest <- binding_manifest_build(
    data = data,
    formula = formula,
    family_link = pvstackr:::pv_binding_family_link_projection(
      "gaussian", "identity", "real", "fixed"
    ),
    estimand_metadata = matching_metadata
  )
  dispersion_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(
      relink(dispersion_manifest),
      dispersion_manifest
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(dispersion_error$code, "PV_BIND_E060")

  source_metadata <- matching_metadata
  source_metadata$target_source <- "external_model_rubin"
  source_manifest <- binding_manifest_build(
    data = data,
    formula = formula,
    estimand_metadata = source_metadata
  )
  source_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(relink(source_manifest), source_manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(source_error$code, "PV_BIND_E070")

  interval_metadata <- matching_metadata
  interval_metadata$interval_role <- "different_interval_role"
  interval_manifest <- binding_manifest_build(
    data = data,
    formula = formula,
    estimand_metadata = interval_metadata
  )
  interval_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(relink(interval_manifest), interval_manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(interval_error$code, "PV_BIND_E070")

  expanded_formula <- OUTCOME ~ x + g
  expanded_metadata <- binding_manifest_metadata(
    data,
    expanded_formula,
    target_source = fixture$content$target_source,
    interval_role = fixture$content$interval_policy$interval_role,
    coverage_claim_allowed = fixture$content$interval_policy$coverage_claim_allowed
  )
  expanded_manifest <- binding_manifest_build(
    data = data,
    formula = expanded_formula,
    estimand_metadata = expanded_metadata
  )
  fe_error <- tryCatch(
    pvstackr:::pv_binding_target_manifest_validate(relink(expanded_manifest), expanded_manifest),
    pvstackr_binding_error = identity
  )
  expect_identical(fe_error$code, "PV_BIND_E071")

  rehash_contrast <- function(manifest) {
    manifest$components$estimand$estimand_contrast_hash <-
      pvstackr:::pv_binding_hash_payload(
        manifest$components$estimand_contrast,
        "estimand_contrast"
      )
    for (component_name in c("estimand_contrast", "estimand")) {
      manifest$component_hashes[[paste0(component_name, "_hash")]] <-
        pvstackr:::pv_binding_component_hash(
          manifest$components[[component_name]],
          component_name
        )
    }
    manifest$manifest_hash <- pvstackr:::pv_binding_hash_payload(
      pvstackr:::pv_binding_manifest_hash_payload(manifest),
      "manifest"
    )
    expect_invisible(pvstackr:::pv_binding_manifest_validate(manifest))
    manifest
  }
  contrast_mutators <- list(
    matrix_hash = function(x) {
      x$matrix_hash <- paste0("sha256:", strrep("a", 64L))
      x
    },
    output_names_hash = function(x) {
      x$ordered_output_names_hash <- paste0("sha256:", strrep("b", 64L))
      x
    },
    coverage_inheritance = function(x) {
      x$coverage_inheritance <- FALSE
      x
    }
  )
  for (name in names(contrast_mutators)) {
    contrast_manifest <- fixture$manifest
    contrast_manifest$components$estimand_contrast <-
      contrast_mutators[[name]](
        contrast_manifest$components$estimand_contrast
      )
    contrast_manifest <- rehash_contrast(contrast_manifest)
    contrast_error <- tryCatch(
      pvstackr:::pv_binding_target_manifest_validate(
        relink(contrast_manifest),
        contrast_manifest
      ),
      pvstackr_binding_error = identity
    )
    expect_identical(contrast_error$code, "PV_BIND_E072", info = name)
  }
})

test_that("3c remediation adapts current pv_brr_target objects to canonical linked content", {
  fixture <- binding_target_content_fixture()
  target_before <- serialize(fixture$target, NULL)
  adapted <- pvstackr:::pv_binding_target_content_from_brr_target(
    fixture$target,
    fixture$manifest,
    conf_level = 0.95,
    data = fixture$data,
    formula = fixture$target$formula
  )
  expect_identical(adapted, fixture$content)
  expect_identical(serialize(fixture$target, NULL), target_before)
  expect_null(adapted$interval_policy$df_complete)
  expect_identical(adapted$target_engine_id, "lm_wls_brr_fay_v1")
  expect_invisible(
    pvstackr:::pv_binding_target_manifest_validate(adapted, fixture$manifest)
  )

  target_with_level <- fixture$target
  target_with_level$conf_level <- 0.95
  expect_identical(
    pvstackr:::pv_binding_target_content_from_brr_target(
      target_with_level,
      fixture$manifest,
      data = fixture$data,
      formula = fixture$target$formula
    ),
    fixture$content
  )

  expect_identical(
    pvstackr:::pv_binding_target_content_from_brr_target(
      fixture$target,
      fixture$manifest,
      conf_level = 0.95
    ),
    fixture$content
  )

  expect_identical(
    pvstackr:::pv_binding_target_content_from_brr_target(
      fixture$target,
      fixture$manifest,
      data = fixture$data,
      formula = fixture$target$formula
    ),
    fixture$content
  )
  override_error <- tryCatch(
    pvstackr:::pv_binding_target_content_from_brr_target(
      fixture$target,
      fixture$manifest,
      conf_level = 0.9
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(override_error$code, "PV_BIND_E090")

  br_fixture <- binding_target_content_fixture(
    df_method = "barnard_rubin",
    df_complete = Inf
  )
  br_adapted <- pvstackr:::pv_binding_target_content_from_brr_target(
    br_fixture$target,
    br_fixture$manifest,
    conf_level = 0.95,
    data = br_fixture$data,
    formula = br_fixture$target$formula
  )
  expect_identical(br_adapted, br_fixture$content)
})

test_that("3c final legacy adapter rejects manifests or targets not rebuilt from one raw source", {
  fixture <- binding_target_content_fixture()
  target <- binding_schema01_target(fixture$target)
  formula <- target$formula
  source_before <- serialize(target, NULL)
  expect_e080 <- function(manifest, data, supplied_formula = formula) {
    error <- tryCatch(
      pvstackr:::pv_binding_target_content_from_brr_target(
        target,
        manifest,
        conf_level = 0.95,
        data = data,
        formula = supplied_formula
      ),
      pvstackr_binding_error = identity
    )
    expect_s3_class(error, "pvstackr_binding_error")
    expect_identical(error$code, "PV_BIND_E080")
  }
  metadata_for <- function(data, formula) {
    binding_manifest_metadata(
      data,
      formula,
      target_source = target$target_source,
      interval_role = target$interval_role,
      coverage_claim_allowed = target$coverage_claim_allowed
    )
  }

  missing_raw_error <- tryCatch(
    pvstackr:::pv_binding_target_content_from_brr_target(
      target,
      fixture$manifest,
      conf_level = 0.95
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(missing_raw_error$code, "PV_BIND_E080")
  missing_level_error <- tryCatch(
    pvstackr:::pv_binding_target_content_from_brr_target(
      target,
      fixture$manifest,
      data = fixture$data,
      formula = formula
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(missing_level_error$code, "PV_BIND_E080")

  changed_pv <- fixture$data
  changed_pv$PV1[[1L]] <- changed_pv$PV1[[1L]] + 10
  changed_pv_manifest <- binding_manifest_build(
    data = changed_pv,
    formula = formula,
    estimand_metadata = metadata_for(changed_pv, formula)
  )
  expect_e080(changed_pv_manifest, changed_pv)

  alternate_weights <- fixture$data
  alternate_weights$W_ALT <- alternate_weights$W
  alternate_weights$RW_ALT1 <- alternate_weights$RW1
  alternate_weights$RW_ALT2 <- alternate_weights$RW2
  alternate_weight_manifest <- pvstackr:::pv_binding_manifest_build(
    data = alternate_weights,
    formula = formula,
    pv_cols = target$pv_cols,
    weight_col = "W_ALT",
    rep_weight_cols = c("RW_ALT1", "RW_ALT2"),
    fay_k = target$fay_k,
    id_cols = target$id_cols,
    family_link = pvstackr:::pv_binding_family_link_projection("gaussian", "identity"),
    estimand_contrast = NULL,
    estimand_metadata = metadata_for(alternate_weights, formula)
  )
  expect_e080(alternate_weight_manifest, alternate_weights)

  offset_formula <- OUTCOME ~ x + offset(x)
  offset_manifest <- binding_manifest_build(
    data = fixture$data,
    formula = offset_formula,
    estimand_metadata = metadata_for(fixture$data, offset_formula)
  )
  expect_e080(offset_manifest, fixture$data, offset_formula)

  expect_identical(serialize(target, NULL), source_before)
})

test_that("3c remediation requires positive-definite T_MI and strictly positive SE", {
  fixture <- binding_target_content_fixture()
  expect_invisible(pvstackr:::pv_binding_target_content_validate(fixture$content))
  expect_lte(
    min(eigen(fixture$content$derived$B, symmetric = TRUE, only.values = TRUE)$values),
    pvstackr:::pv_binding_target_matrix_tolerance()
  )

  ill_conditioned <- fixture$content$derived
  fe_names <- names(ill_conditioned$beta)
  ill_conditioned$T_MI <- diag(c(1, 1e-14))
  dimnames(ill_conditioned$T_MI) <- list(fe_names, fe_names)
  ill_conditioned$se <- stats::setNames(
    sqrt(diag(ill_conditioned$T_MI)),
    fe_names
  )
  expect_invisible(pvstackr:::pv_binding_target_derived_validate(ill_conditioned))
  expect_silent(
    pvstackr:::pv_validate_target_matrix(ill_conditioned$T_MI, fe_names)
  )

  zero_primitives <- fixture$target$per_pv
  reference_beta <- zero_primitives[[1L]]$beta
  zero_U <- matrix(
    0,
    nrow = length(reference_beta),
    ncol = length(reference_beta),
    dimnames = list(names(reference_beta), names(reference_beta))
  )
  zero_primitives <- lapply(zero_primitives, function(item) {
    item$beta <- reference_beta
    item$U <- zero_U
    item
  })
  zero_error <- tryCatch(
    pvstackr:::pv_binding_target_content_build(
      per_pv = zero_primitives,
      M = fixture$target$M,
      R = fixture$target$R,
      fay_k = fixture$target$fay_k,
      df_method = "classic",
      df_complete = NULL,
      conf_level = 0.95,
      target_source = fixture$target$target_source,
      target_engine_id = "lm_wls_brr_fay_v1",
      target_policy = fixture$target$policy,
      manifest_hash = fixture$manifest$manifest_hash
    ),
    pvstackr_binding_error = identity
  )
  expect_identical(zero_error$code, "PV_BIND_E090")
  expect_match(conditionMessage(zero_error), "positive definite")
})

binding_manifest_rehash_component <- function(manifest, component_name) {
  hash_name <- paste0(component_name, "_hash")
  manifest$component_hashes[[hash_name]] <- pvstackr:::pv_binding_component_hash(
    manifest$components[[component_name]],
    component_name
  )
  manifest$manifest_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(manifest),
    "manifest"
  )
  manifest
}

test_that("3b remediation deep-validates all ten rehashed component schemas", {
  manifest <- binding_manifest_build()
  mutators <- list(
    row = function(x) { x$analysis_row_policy <- list("na.fail"); x },
    pv = function(x) { names(x$types) <- rev(names(x$types)); x },
    predictor = function(x) { x$columns[[1L]]$type_id <- "character"; x },
    formula = function(x) { x$intercept <- 2L; x },
    model_matrix = function(x) { x$assign <- x$assign[-1L]; x },
    weight = function(x) { x$stack_transform_id <- "forged_transform"; x },
    factor_contrast = function(x) { x$option_snapshot$count <- x$option_snapshot$count + 1L; x },
    estimand_contrast = function(x) { x$output_count <- x$output_count + 1L; x },
    family_link = function(x) { x$family_id <- " "; x },
    estimand = function(x) { x$coverage_claim_allowed <- NA; x }
  )

  for (component_name in names(mutators)) {
    forged <- manifest
    forged$components[[component_name]] <- mutators[[component_name]](
      forged$components[[component_name]]
    )
    forged <- binding_manifest_rehash_component(forged, component_name)
    error <- tryCatch(
      pvstackr:::pv_binding_manifest_validate(forged),
      pvstackr_binding_error = identity
    )
    expect_true(inherits(error, "pvstackr_binding_error"), info = component_name)
    expect_true(error$code %in% c("PV_BIND_E002", "PV_BIND_E081"), info = component_name)
  }

  stale <- manifest
  stale$components$row$analysis_row_policy <- list("na.fail")
  stale_error <- tryCatch(
    pvstackr:::pv_binding_manifest_validate(stale),
    pvstackr_binding_error = identity
  )
  expect_identical(stale_error$code, "PV_BIND_E005")
})

test_that("3b remediation comparator never reports ok for an unmapped component mismatch", {
  expected <- binding_manifest_build()
  observed <- expected
  forged_hash <- pvstackr:::pv_binding_hash_payload(
    "forged but schema-valid aggregate",
    "predictor"
  )
  observed$components$predictor$ordered_values_hash <- forged_hash
  observed$components$model_matrix$predictor_values_hash <- forged_hash
  observed <- binding_manifest_rehash_component(observed, "predictor")
  observed <- binding_manifest_rehash_component(observed, "model_matrix")
  expect_invisible(pvstackr:::pv_binding_manifest_validate(observed))

  comparison <- pvstackr:::pv_binding_manifest_compare(expected, observed)
  expect_false(comparison$ok)
  expect_identical(comparison$mismatched_components, c("predictor", "model_matrix"))
  expect_identical(comparison$all_codes, c("PV_BIND_E031", "PV_BIND_E042"))
  expect_identical(comparison$primary_code, "PV_BIND_E031")
  expect_identical(
    comparison$ok,
    length(comparison$mismatched_components) == 0L
  )
})

test_that("3b remediation comparator codes remain in registry order under multiple mutations", {
  data <- binding_component_fixture()
  expected <- binding_manifest_build(data = data)
  changed <- data
  changed$PV1[[1L]] <- changed$PV1[[1L]] + 0.25
  changed$x[[2L]] <- changed$x[[2L]] + 0.25
  changed$RW2[[3L]] <- changed$RW2[[3L]] + 0.25
  family <- pvstackr:::pv_binding_family_link_projection("gaussian", "log")
  metadata <- binding_manifest_metadata(
    changed,
    OUTCOME ~ x + g,
    target_source = "external_model_rubin"
  )
  observed <- binding_manifest_build(
    data = changed,
    family_link = family,
    estimand_metadata = metadata
  )
  comparison <- pvstackr:::pv_binding_manifest_compare(expected, observed)
  expect_identical(
    comparison$all_codes,
    c(
      "PV_BIND_E021", "PV_BIND_E031", "PV_BIND_E042", "PV_BIND_E051",
      "PV_BIND_E061", "PV_BIND_E070"
    )
  )
  expect_identical(comparison$primary_code, "PV_BIND_E021")
  expect_false(comparison$ok)
})

test_that("3b remediation keyed matrix hashes detect closure changes during row reorder", {
  data <- binding_component_fixture()
  env_identity <- new.env(parent = baseenv())
  env_identity$transform_x <- function(x) x
  formula_identity <- stats::as.formula("OUTCOME ~ transform_x(x)", env = env_identity)
  env_doubled <- new.env(parent = baseenv())
  env_doubled$transform_x <- function(x) 2 * x
  formula_doubled <- stats::as.formula("OUTCOME ~ transform_x(x)", env = env_doubled)

  expected <- binding_manifest_build(data = data, formula = formula_identity)
  observed <- binding_manifest_build(data = data[6:1, ], formula = formula_doubled)
  comparison <- pvstackr:::pv_binding_manifest_compare(expected, observed)
  expect_identical(comparison$all_codes, c("PV_BIND_E012", "PV_BIND_E042"))
  expect_false("PV_BIND_E040" %in% comparison$all_codes)

  expected_fallback <- binding_manifest_build(
    data = data,
    formula = formula_identity,
    id_cols = NULL
  )
  observed_fallback <- binding_manifest_build(
    data = data[6:1, ],
    formula = formula_doubled,
    id_cols = NULL
  )
  fallback_comparison <- pvstackr:::pv_binding_manifest_compare(
    expected_fallback,
    observed_fallback
  )
  expect_identical(fallback_comparison$all_codes, c("PV_BIND_E012", "PV_BIND_E042"))

  env_offset_a <- new.env(parent = asNamespace("stats"))
  env_offset_a$offset_fun <- function(x) x
  formula_offset_a <- stats::as.formula(
    "OUTCOME ~ x + offset(offset_fun(x))",
    env = env_offset_a
  )
  env_offset_b <- new.env(parent = asNamespace("stats"))
  env_offset_b$offset_fun <- function(x) 2 * x
  formula_offset_b <- stats::as.formula(
    "OUTCOME ~ x + offset(offset_fun(x))",
    env = env_offset_b
  )
  offset_result <- pvstackr:::pv_binding_manifest_compare(
    binding_manifest_build(data = data, formula = formula_offset_a),
    binding_manifest_build(data = data[6:1, ], formula = formula_offset_b)
  )
  expect_identical(offset_result$all_codes, c("PV_BIND_E012", "PV_BIND_E042"))
})

test_that("3b remediation optional root tail and excluded metadata are privacy-safe exact schemas", {
  manifest <- binding_manifest_build()
  migration <- list(
    migration_from_schema = "0.1.0",
    migration_function = "pv_binding_revalidate_brr_target_v1",
    binding_revalidated = TRUE,
    inspection_only = FALSE,
    warnings = "PV_BIND_MIGRATION_LEGACY_HASH_ONLY"
  )
  valid <- c(
    manifest,
    list(
      legacy_hashes = list(
        algorithm_id = pvstackr:::pv_binding_legacy_hash_algorithm_id(),
        design_hash = "1234abcd"
      ),
      created_at = "2026-07-12T12:00:00Z",
      migration = migration
    )
  )
  valid$manifest_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(valid),
    "manifest"
  )
  expect_invisible(pvstackr:::pv_binding_manifest_validate(valid))
  fractional <- valid
  fractional$created_at <- "2026-07-12T12:00:00.123456Z"
  expect_invisible(pvstackr:::pv_binding_manifest_validate(fractional))

  wrong_order <- c(
    manifest,
    list(
      migration = migration,
      created_at = "2026-07-12T12:00:00Z",
      legacy_hashes = list(
        algorithm_id = pvstackr:::pv_binding_legacy_hash_algorithm_id(),
        design_hash = "1234abcd"
      )
    )
  )
  expect_error(
    pvstackr:::pv_binding_manifest_validate(wrong_order),
    class = "pvstackr_binding_error"
  )

  privacy_mutators <- list(
    function(x) { x$legacy_hashes <- data.frame(id = 1:2); x },
    function(x) { x$created_at <- charToRaw("respondent payload"); x },
    function(x) { x$migration$warnings <- data.frame(id = 1:2); x },
    function(x) { x$migration$warnings <- new.env(parent = emptyenv()); x },
    function(x) { x$migration$warnings <- list(raw = 1:3); x },
    function(x) { x$migration <- new.env(parent = emptyenv()); x }
  )
  for (mutate in privacy_mutators) {
    expect_error(
      pvstackr:::pv_binding_manifest_validate(mutate(valid)),
      class = "pvstackr_binding_error"
    )
  }
})

test_that("3b final remediation rejects ignored attributes and optional payload channels", {
  manifest <- binding_manifest_build()
  valid <- c(
    manifest,
    list(
      legacy_hashes = list(
        algorithm_id = pvstackr:::pv_binding_legacy_hash_algorithm_id(),
        design_hash = "1234abcd",
        target_hash = "90abcdef"
      ),
      created_at = "2026-07-12T12:00:00Z",
      migration = list(
        migration_from_schema = "0.1.0",
        migration_function = "pv_binding_revalidate_brr_target_v1",
        binding_revalidated = TRUE,
        inspection_only = FALSE,
        warnings = "PV_BIND_MIGRATION_LEGACY_HASH_ONLY"
      )
    )
  )
  valid$manifest_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(valid),
    "manifest"
  )
  expect_invisible(pvstackr:::pv_binding_manifest_validate(valid))

  marker <- "benign_respondent_marker"
  forged <- list(
    root_attribute = function(x) {
      attr(x, "marker") <- marker
      x
    },
    component_attribute = function(x) {
      attr(x$components$pv, "marker") <- marker
      x
    },
    created_at_attribute = function(x) {
      attr(x$created_at, "marker") <- marker
      x
    },
    legacy_list_attribute = function(x) {
      attr(x$legacy_hashes, "marker") <- marker
      x
    },
    legacy_value_attribute = function(x) {
      attr(x$legacy_hashes$design_hash, "marker") <- marker
      x
    },
    migration_attribute = function(x) {
      attr(x$migration, "marker") <- marker
      x
    },
    migration_scalar_attribute = function(x) {
      attr(x$migration$migration_function, "marker") <- marker
      x
    },
    arbitrary_legacy_key = function(x) {
      x$legacy_hashes$respondent_payload <- "deadbeef"
      x
    },
    unsupported_legacy_algorithm = function(x) {
      x$legacy_hashes$algorithm_id <- "respondent_selected_algorithm"
      x
    },
    impossible_created_at = function(x) {
      x$created_at <- "2026-02-30T12:00:00Z"
      x
    },
    invalid_created_at_hour = function(x) {
      x$created_at <- "2026-07-12T99:00:00Z"
      x
    },
    excessive_created_at_precision = function(x) {
      x$created_at <- "2026-07-12T12:00:00.1234567Z"
      x
    },
    arbitrary_migration_schema = function(x) {
      x$migration$migration_from_schema <- "respondent_schema"
      x
    },
    arbitrary_migration_function = function(x) {
      x$migration$migration_function <- "respondent_migration_function"
      x
    },
    repeated_migration_schema = function(x) {
      x$migration$migration_from_schema <- rep("0.1.0", 1000L)
      x
    },
    repeated_migration_function = function(x) {
      x$migration$migration_function <- rep(
        "pv_binding_revalidate_brr_target_v1",
        1000L
      )
      x
    },
    long_warning_vector = function(x) {
      x$migration$warnings <- rep(
        "PV_BIND_MIGRATION_LEGACY_HASH_ONLY",
        1000L
      )
      x
    }
  )

  for (mutate in forged) {
    changed <- mutate(valid)
    validation_error <- tryCatch(
      pvstackr:::pv_binding_manifest_validate(changed),
      pvstackr_binding_error = identity
    )
    expect_s3_class(validation_error, "pvstackr_binding_error")
    expect_identical(validation_error$code, "PV_BIND_E002")
    expect_error(
      pvstackr:::pv_binding_manifest_compare(valid, changed),
      class = "pvstackr_binding_error"
    )
  }

  expect_error(
    pvstackr:::pv_binding_manifest_hash_payload(forged$root_attribute(valid)),
    class = "pvstackr_binding_error"
  )
  expect_error(
    pvstackr:::pv_binding_manifest_hash_payload(forged$created_at_attribute(valid)),
    class = "pvstackr_binding_error"
  )

  serialized <- paste(
    rawToChar(serialize(valid, NULL, ascii = TRUE), multiple = TRUE),
    collapse = ""
  )
  expect_false(grepl(marker, serialized, fixed = TRUE))
})

test_that("canonical TLV rejects list and atomic attributes outside exact names", {
  marked_list <- list(a = 1L)
  attr(marked_list, "marker") <- "ignored before remediation"
  expect_error(
    pvstackr:::pv_binding_hash_payload(marked_list, "manifest"),
    class = "pvstackr_binding_error"
  )

  marked_atomic <- c(a = 1L)
  attr(marked_atomic, "marker") <- "ignored before remediation"
  expect_error(
    pvstackr:::pv_binding_hash_payload(marked_atomic, "manifest"),
    class = "pvstackr_binding_error"
  )

  expect_silent(pvstackr:::pv_binding_hash_payload(list(a = c(x = 1L)), "manifest"))
})

test_that("3b remediation E005 condition retains only fixed SHA-256 prefixes", {
  manifest <- binding_manifest_build()
  full_observed <- paste0("sha256:", strrep("0", 64L))
  tampered <- manifest
  tampered$manifest_hash <- full_observed
  error <- tryCatch(
    pvstackr:::pv_binding_manifest_validate(tampered),
    pvstackr_binding_error = identity
  )

  prefix_pattern <- sprintf(
    "^sha256:[0-9a-f]{%d}$",
    pvstackr:::pv_binding_hash_prefix_hex_length()
  )
  expect_identical(error$code, "PV_BIND_E005")
  expect_match(error$expected_hash, prefix_pattern)
  expect_match(error$observed_hash, prefix_pattern)
  expect_lt(nchar(error$expected_hash), 71L)
  expect_lt(nchar(error$observed_hash), 71L)
  expect_false(grepl(full_observed, conditionMessage(error), fixed = TRUE))
  expect_false(identical(error$observed_hash, full_observed))
  condition_dump <- paste(capture.output(str(unclass(error))), collapse = "\n")
  expect_false(grepl(full_observed, condition_dump, fixed = TRUE))
  expect_identical(
    pvstackr:::pv_binding_sanitize_hash_metadata(data.frame(secret = 1)),
    "malformed"
  )
})

test_that("090 errata E4 stores reportable FE labels once, only in estimand", {
  manifest <- binding_manifest_build()
  fe_names <- manifest$components$estimand$fe_names
  flattened <- unname(unlist(manifest$components, recursive = TRUE, use.names = FALSE))

  expect_true(all(vapply(fe_names, function(name) sum(flattened == name) == 1L, logical(1))))
  expect_false(any(c("input_fe_names", "output_names") %in%
    names(manifest$components$estimand_contrast)))
  expect_identical(
    manifest$components$estimand_contrast$input_fe_count,
    as.integer(length(fe_names))
  )
  expect_match(
    manifest$components$estimand_contrast$ordered_input_fe_names_hash,
    "^sha256:[0-9a-f]{64}$"
  )

  reversed <- pvstackr:::pv_binding_estimand_contrast_projection(rev(fe_names))
  expect_false(identical(
    manifest$components$estimand_contrast$ordered_input_fe_names_hash,
    reversed$ordered_input_fe_names_hash
  ))
  output_names <- paste0("contrast_", seq_along(fe_names))
  identity_L <- diag(length(fe_names))
  dimnames(identity_L) <- list(output_names, fe_names)
  renamed_output <- pvstackr:::pv_binding_estimand_contrast_projection(
    fe_names,
    L = identity_L,
    output_names = output_names,
    contrast_id = "renamed_identity_v1"
  )
  expect_false(identical(
    manifest$components$estimand_contrast$ordered_output_names_hash,
    renamed_output$ordered_output_names_hash
  ))

  renamed_data <- binding_component_fixture()
  levels(renamed_data$g) <- c("renamed_category_a", "renamed_category_b")
  renamed_manifest <- binding_manifest_build(data = renamed_data)
  renamed_comparison <- pvstackr:::pv_binding_manifest_compare(
    manifest,
    renamed_manifest
  )
  expect_true("PV_BIND_E071" %in% renamed_comparison$all_codes)
  expect_true("PV_BIND_E072" %in% renamed_comparison$all_codes)
})

test_that("cosmetic survey metadata does not block or shift the binding manifest", {
  # PISA and its siblings arrive from SAS/SPSS with a `label` on nearly every
  # column and a `format.*` string on many. Those must neither refuse the fit
  # nor change the binding hash.
  plain <- data.frame(
    OUTCOME = c(1.5, 2.5, 3.5, 4.5),
    escs = c(-1, -0.5, 0.5, 1),
    female = c(0L, 1L, 0L, 1L)
  )
  labelled <- plain
  attr(labelled$escs, "label") <- "Index of economic, social and cultural status"
  attr(labelled$escs, "format.sas") <- "F8.2"
  attr(labelled$female, "label") <- "Female"

  bare <- pvstackr:::pv_binding_predictor_projection(plain, OUTCOME ~ escs + female)
  tagged <- pvstackr:::pv_binding_predictor_projection(labelled, OUTCOME ~ escs + female)

  expect_identical(tagged$schema_hash, bare$schema_hash)
  expect_identical(tagged$ordered_values_hash, bare$ordered_values_hash)
  expect_identical(
    vapply(tagged$columns, function(x) x$values_hash, character(1)),
    vapply(bare$columns, function(x) x$values_hash, character(1))
  )
})

test_that("plain column canonicalization keeps classed columns for their own branch", {
  strip <- pvstackr:::pv_binding_column_values

  # Cosmetic attributes go; names stay.
  tagged <- structure(c(a = 1, b = 2), label = "x", format.sas = "F8.2")
  expect_identical(strip(tagged), c(a = 1, b = 2))

  # Classed columns are handed on untouched, so the canonicaliser's typed
  # branches still see their class and still police their own attributes.
  f <- factor(c("a", "b"))
  expect_identical(strip(f), f)
  d <- as.Date("2026-08-14")
  expect_identical(strip(d), d)

  # A value change still moves the hash: stripping metadata is not stripping
  # content.
  expect_false(identical(
    pvstackr:::pv_binding_hash_payload(list(v = strip(c(1, 2))), "predictor"),
    pvstackr:::pv_binding_hash_payload(list(v = strip(c(1, 3))), "predictor")
  ))
})
