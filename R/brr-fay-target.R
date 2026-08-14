pv_escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

pv_pv_suffix_families <- function(data, prefix) {
  if (!is.data.frame(data) || is.null(names(data))) {
    return(data.frame(suffix = character(), example = character(), stringsAsFactors = FALSE))
  }
  pattern <- paste0("^", pv_escape_regex(prefix), "([0-9]+)([^0-9].*)$")
  hits <- grep(pattern, names(data), value = TRUE)
  if (length(hits) == 0L) {
    return(data.frame(suffix = character(), example = character(), stringsAsFactors = FALSE))
  }
  index <- as.integer(sub(pattern, "\\1", hits))
  suffix <- sub(pattern, "\\2", hits)
  suffixes <- sort(unique(suffix))
  examples <- vapply(suffixes, function(x) {
    family_hits <- hits[suffix == x]
    family_index <- index[suffix == x]
    family_hits[order(family_index)][[1L]]
  }, character(1))
  data.frame(suffix = suffixes, example = examples, stringsAsFactors = FALSE)
}

pv_format_pv_suffix_guidance <- function(data, prefix, suffix) {
  families <- pv_pv_suffix_families(data, prefix)
  if (nrow(families) == 0L) {
    return("")
  }
  shown <- families[seq_len(min(3L, nrow(families))), , drop = FALSE]
  family_text <- paste(
    sprintf("`%s` (suffix `%s`)", shown$example, shown$suffix),
    collapse = ", "
  )
  recommended <- families$suffix[[1L]]
  if (identical(suffix, "")) {
    return(paste(
      "Modern PISA files commonly use subject-suffixed plausible values,",
      sprintf("such as %s.", family_text),
      sprintf(
        "If calling `pv_design()` or `pv_brr_target()`, try `pv_suffix = \"%s\"`; if calling `detect_pisa_pv_columns()`, use `suffix = \"%s\"`; or pass explicit `pv_cols`.",
        recommended,
        recommended
      )
    ))
  }
  paste(
    sprintf("Detected plausible-value suffix family/families: %s.", family_text),
    "Check the subject suffix for this PISA cycle or pass explicit `pv_cols`."
  )
}

pv_warn_if_bare_pv_with_subject_families <- function(data, prefix, suffix) {
  if (!identical(suffix, "")) {
    return(invisible(FALSE))
  }
  families <- pv_pv_suffix_families(data, prefix)
  if (nrow(families) == 0L) {
    return(invisible(FALSE))
  }
  shown <- families[seq_len(min(3L, nrow(families))), , drop = FALSE]
  family_text <- paste(
    sprintf("`%s` (suffix `%s`)", shown$example, shown$suffix),
    collapse = ", "
  )
  warning(
    paste(
      "pvstackr: bare plausible-value columns were detected, but subject-suffixed",
      sprintf("PISA plausible-value families were also found: %s.", family_text),
      "Confirm bare columns are intended, set `pv_suffix`/`suffix` to the subject suffix, or pass explicit `pv_cols`."
    ),
    call. = FALSE
  )
  invisible(TRUE)
}

pv_hash_payload <- function(payload) {
  bytes <- as.integer(serialize(payload, NULL, ascii = FALSE, xdr = TRUE, version = 2))
  if (length(bytes) == 0L) {
    return("00000000")
  }

  # Adler-32 over canonical R serialization keeps the existing 8-hex contract
  # while hashing full payload bytes rather than truncated display text.
  modulus <- 65521
  a <- 1
  b <- 0
  chunk_size <- 5552L
  for (start in seq.int(1L, length(bytes), by = chunk_size)) {
    end <- min(start + chunk_size - 1L, length(bytes))
    chunk <- bytes[start:end]
    csum <- cumsum(chunk)
    b <- (b + length(chunk) * a + sum(csum)) %% modulus
    a <- (a + sum(chunk)) %% modulus
  }
  sprintf("%04x%04x", b, a)
}

pv_natural_prefixed_cols <- function(data, prefix, suffix, label) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  prefix <- pv_assert_scalar_string(prefix, paste0(label, "_prefix"))
  if (!is.character(suffix) || length(suffix) != 1L || is.na(suffix)) {
    pv_abort(sprintf("`%s_suffix` must be a character scalar.", label))
  }

  pattern <- paste0("^", pv_escape_regex(prefix), "([0-9]+)", pv_escape_regex(suffix), "$")
  hits <- grep(pattern, names(data), value = TRUE)
  if (length(hits) == 0L) {
    msg <- sprintf(
      "No %s columns detected with prefix `%s` and suffix `%s`.",
      label,
      prefix,
      suffix
    )
    if (identical(label, "plausible-value")) {
      guidance <- pv_format_pv_suffix_guidance(data, prefix, suffix)
      if (nzchar(guidance)) {
        msg <- paste(msg, guidance)
      }
    }
    pv_abort(msg)
  }

  index <- as.integer(sub(pattern, "\\1", hits))
  if (anyDuplicated(index)) {
    pv_abort(sprintf("Duplicate %s numeric suffixes detected.", label))
  }
  ordered <- sort(index)
  if (!identical(ordered, seq_len(length(ordered)))) {
    pv_abort(sprintf("%s numeric suffixes must be contiguous from 1.", label))
  }

  hits[order(index)]
}

#' Detect PISA-Style Plausible-Value Columns
#'
#' Modern PISA plausible values are often subject-suffixed, for example
#' `PV1MATH`, `PV2MATH`, or `PV1READ`, `PV2READ`. Pass the subject/domain
#' suffix explicitly, such as `suffix = "MATH"` here or
#' `pv_suffix = "MATH"` in [pv_design()] and [pv_brr_target()]. The default
#' `suffix = ""` is reserved for data whose plausible values are intentionally
#' named as bare `PV1`, `PV2`, and so on. Detection is anchored, so
#' `suffix = ""` does not match `PV1MATH`. Use `expected_M` to guard against
#' selecting the wrong subject or an incomplete plausible-value set.
#'
#' @param data Data frame containing plausible-value columns.
#' @param prefix Character scalar. Column prefix before the numeric
#'   plausible-value index. Default `"PV"`.
#' @param suffix Character scalar. Column suffix after the numeric
#'   plausible-value index, such as `"MATH"` for columns named `PV1MATH`,
#'   `PV2MATH`, and so on. Default `""` matches bare `PV1`, `PV2`, ... only;
#'   detection is anchored, so the default does not match subject-suffixed
#'   columns.
#' @param expected_M Optional integer scalar. Expected plausible-value count;
#'   if supplied, detection errors unless exactly `expected_M` columns are
#'   found. If `NULL` (default), the count is not checked.
#'
#' @returns Character vector of column names in natural numeric order.
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#' # Subject-suffixed plausible values: pass the subject suffix explicitly.
#' detect_pisa_pv_columns(pisa_tiny, suffix = "READ", expected_M = 2L)
#' @family pvstackr-detection
#' @seealso
#' [detect_pisa_brr_replicate_weights()]; build a design with [pv_design()] or a
#' target with [pv_brr_target()].
#' @export
detect_pisa_pv_columns <- function(data, prefix = "PV", suffix = "", expected_M = NULL) {
  cols <- pv_natural_prefixed_cols(data, prefix, suffix, "plausible-value")
  pv_warn_if_bare_pv_with_subject_families(data, prefix, suffix)
  if (!is.null(expected_M)) {
    expected_M <- pv_assert_scalar_number(expected_M, "expected_M", integer = TRUE, lower = 1)
    if (length(cols) != expected_M) {
      pv_abort(sprintf("Expected %d plausible-value columns, detected %d.", expected_M, length(cols)))
    }
  }
  cols
}

#' Detect PISA-Style BRR Replicate-Weight Columns
#'
#' PISA balanced-repeated-replication (BRR) replicate weights are named as a
#' prefix followed by a numeric replicate index, for example `W_FSTURWT1`,
#' `W_FSTURWT2`, and so on. This finds the `<prefix><r>` columns and returns
#' them in natural numeric order using the default prefix `"W_FSTURWT"`. It is
#' the replicate-weight companion to [detect_pisa_pv_columns()]. Use
#' `expected_R` to guard against an incomplete replicate set or the wrong
#' prefix.
#'
#' @param data Data frame containing replicate-weight columns.
#' @param prefix Character scalar. Column prefix before the numeric replicate
#'   index. Default `"W_FSTURWT"`.
#' @param expected_R Optional integer scalar. Expected replicate-weight count;
#'   if supplied, detection errors unless exactly `expected_R` columns are
#'   found. If `NULL` (default), the count is not checked.
#'
#' @returns Character vector of column names in natural numeric order.
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#' detect_pisa_brr_replicate_weights(pisa_tiny, expected_R = 4L)
#' @family pvstackr-detection
#' @seealso
#' [detect_pisa_pv_columns()]; build a design with [pv_design()] or a target
#' with [pv_brr_target()].
#' @export
detect_pisa_brr_replicate_weights <- function(data, prefix = "W_FSTURWT", expected_R = NULL) {
  cols <- pv_natural_prefixed_cols(data, prefix, "", "replicate-weight")
  if (!is.null(expected_R)) {
    expected_R <- pv_assert_scalar_number(expected_R, "expected_R", integer = TRUE, lower = 1)
    if (length(cols) != expected_R) {
      pv_abort(sprintf("Expected %d replicate-weight columns, detected %d.", expected_R, length(cols)))
    }
  }
  cols
}

pv_validate_id_columns <- function(data, id_cols) {
  if (is.null(id_cols) || length(id_cols) == 0L) {
    return(character())
  }
  id_cols <- pv_validate_columns(data, id_cols, "id_cols")
  row_key <- data[id_cols]
  if (anyNA(row_key)) {
    pv_abort("`id_cols` must not contain missing values.")
  }
  if (anyDuplicated(row_key)) {
    pv_abort("`id_cols` must jointly identify unique rows.")
  }
  id_cols
}

pv_brr_target_formula_rhs <- function(formula) {
  rhs <- pv_formula_rhs_checked(formula)
  if (pv_formula_has_weights_call(rhs)) {
    pv_abort("Do not embed `weights()` in `formula`; pass `weight_col` and `rep_weight_cols` explicitly.")
  }
  if (pv_formula_has_random_effect_bar(rhs)) {
    pv_abort("Random-effect formula terms are not supported by the base WLS BRR-Fay target engine yet.")
  }
  rhs
}

pv_formula_for_pv <- function(formula, pv_col) {
  rhs <- pv_brr_target_formula_rhs(formula)
  stats::as.formula(
    paste(pv_col, paste(deparse(rhs, width.cutoff = 500L), collapse = ""), sep = " ~ "),
    env = environment(formula)
  )
}

pv_wls_beta_from_bundle <- function(model_bundle, outcome, weights) {
  x <- model_bundle$model_matrix
  offset <- model_bundle$offset
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
      !is.numeric(outcome) || length(outcome) != nrow(x) ||
      any(!is.finite(outcome)) ||
      !(is.null(offset) ||
        (is.numeric(offset) && length(offset) == nrow(x) && all(is.finite(offset))))) {
    pv_abort("Resolved BRR-Fay WLS inputs must contain one finite outcome per model-matrix row and an aligned finite offset.")
  }
  weights <- pv_validate_weight_vector(weights, "weights", nrow(x))
  fit <- stats::lm.wfit(
    x = x,
    y = unname(as.double(outcome)),
    w = weights,
    offset = offset
  )
  beta <- stats::coef(fit)
  if (any(!is.finite(beta))) {
    pv_abort("Weighted least-squares fit produced non-finite fixed-effect estimates.")
  }
  names(beta) <- pv_normalize_fe_names(names(beta))
  beta
}

pv_normalize_fe_names <- function(names) {
  out <- paste0("b_", names)
  out[out == "b_(Intercept)"] <- "b_Intercept"
  out
}

pv_brr_fay_one_pv <- function(
  data,
  model_bundle,
  pv_col,
  weight_col,
  rep_weight_cols,
  fay_k
) {
  n <- nrow(data)
  pv_validate_columns(data, c(pv_col, weight_col, rep_weight_cols), "BRR-Fay input columns")
  if (!is.numeric(data[[pv_col]]) || any(!is.finite(data[[pv_col]]))) {
    pv_abort(sprintf("Plausible-value column `%s` must be finite and numeric.", pv_col))
  }
  base_w <- pv_validate_weight_vector(data[[weight_col]], weight_col, n)
  rep_w <- lapply(rep_weight_cols, function(col) {
    pv_validate_weight_vector(data[[col]], col, n)
  })
  beta0 <- pv_wls_beta_from_bundle(model_bundle, data[[pv_col]], base_w)

  p <- length(beta0)
  R <- length(rep_weight_cols)
  diff <- matrix(NA_real_, nrow = p, ncol = R)
  for (r in seq_len(R)) {
    beta_r <- pv_wls_beta_from_bundle(model_bundle, data[[pv_col]], rep_w[[r]])
    if (!identical(names(beta_r), names(beta0))) {
      pv_abort("Fixed-effect names differ across replicate fits.")
    }
    diff[, r] <- beta_r - beta0
  }

  multiplier <- 1 / (R * (1 - fay_k)^2)
  U <- pv_symmetrize(multiplier * tcrossprod(diff))
  dimnames(U) <- list(names(beta0), names(beta0))
  dimnames(diff) <- list(names(beta0), rep_weight_cols)
  replicate_beta <- sweep(diff, 1L, beta0, FUN = "+")

  list(
    beta = beta0,
    U = U,
    fe_names = names(beta0),
    pv_col = pv_col,
    R = R,
    fay_k = fay_k,
    fay_variance_multiplier = multiplier,
    replicate_beta = replicate_beta,
    replicate_diff = diff
  )
}

pv_validate_named_numeric <- function(x, name, names_expected, allow_infinite = FALSE) {
  if (!is.numeric(x) || is.null(names(x)) || !identical(names(x), names_expected)) {
    pv_abort(sprintf("BRR-Fay target `%s` must be a named numeric vector aligned to `fe_names`.", name))
  }
  bad_values <- if (allow_infinite) is.na(x) else !is.finite(x)
  if (any(bad_values)) {
    pv_abort(sprintf("BRR-Fay target `%s` contains invalid numeric values.", name))
  }
  x
}

pv_validate_aligned_symmetric_matrix <- function(x, name, fe_names) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
    pv_abort(sprintf("BRR-Fay target `%s` must be a finite numeric matrix.", name))
  }
  if (!identical(dim(x), c(length(fe_names), length(fe_names)))) {
    pv_abort(sprintf("BRR-Fay target `%s` dimensions must match `fe_names`.", name))
  }
  if (!identical(rownames(x), fe_names) || !identical(colnames(x), fe_names)) {
    pv_abort(sprintf("BRR-Fay target `%s` row and column names must match `fe_names`.", name))
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-10))) {
    pv_abort(sprintf("BRR-Fay target `%s` must be symmetric.", name))
  }
  x
}

pv_validate_hash_scalar <- function(x, name) {
  x <- pv_assert_scalar_string(x, name)
  if (!grepl("^[0-9a-f]{8}$", x)) {
    pv_abort(sprintf("BRR-Fay target `%s` must be an 8-character hexadecimal hash.", name))
  }
  x
}

pv_validate_df_complete_field <- function(x, fe_names, df_method) {
  if (!is.numeric(x) || is.null(names(x)) || !identical(names(x), fe_names)) {
    pv_abort("BRR-Fay target `df_complete` must be a named numeric vector aligned to `fe_names`.")
  }
  if (identical(df_method, "classic")) {
    if (!all(is.na(x))) {
      pv_abort("BRR-Fay target `df_complete` must be all NA when `df_method = \"classic\"`.")
    }
  } else if (identical(df_method, "barnard_rubin")) {
    if (any(is.na(x)) || any(x <= 0)) {
      pv_abort("BRR-Fay target `df_complete` must be positive when `df_method = \"barnard_rubin\"`.")
    }
  } else {
    pv_abort("BRR-Fay target `df_method` must be `classic` or `barnard_rubin`.")
  }
  x
}

pv_normalize_df_complete_input <- function(df_complete, fe_names, df_method) {
  df_method <- match.arg(df_method, c("classic", "barnard_rubin"))
  if (identical(df_method, "classic")) {
    if (!is.null(df_complete)) {
      pv_abort("`df_complete` is only supported when `df_method = \"barnard_rubin\"`.")
    }
    return(NULL)
  }

  if (is.null(df_complete)) {
    pv_abort("`df_complete` is required when `df_method = \"barnard_rubin\"`.")
  }
  if (!is.numeric(df_complete) || length(df_complete) < 1L ||
      any(is.na(df_complete)) || any(df_complete <= 0)) {
    pv_abort("`df_complete` must be positive.")
  }
  if (length(df_complete) == 1L) {
    return(stats::setNames(rep(df_complete, length(fe_names)), fe_names))
  }
  if (length(df_complete) != length(fe_names) || is.null(names(df_complete)) ||
      anyDuplicated(names(df_complete)) || !setequal(names(df_complete), fe_names)) {
    pv_abort("`df_complete` must be scalar or a named numeric vector aligned to `fe_names`.")
  }
  df_complete[fe_names]
}

pv_df_policy <- function(df_method) {
  df_method <- match.arg(df_method, c("classic", "barnard_rubin"))
  if (identical(df_method, "barnard_rubin")) {
    list(
      interval_role = "coverage_barnard_rubin",
      coverage_claim_allowed = TRUE
    )
  } else {
    list(
      interval_role = "descriptive_classic_rubin",
      coverage_claim_allowed = FALSE
    )
  }
}

pv_brr_target_estimand_metadata <- function(model_bundle, df_policy) {
  fe_names <- pv_normalize_fe_names(colnames(model_bundle$model_matrix))
  metadata <- list(
    estimand_id = "brr_fay_fixed_effects",
    target_source = "external_brr_fay_rubin",
    target_engine_id = "lm_wls_brr_fay_v1",
    parameter_scope = "fixed_effect",
    fe_names = fe_names,
    interval_role = df_policy$interval_role,
    coverage_claim_allowed = df_policy$coverage_claim_allowed
  )
  expected_fields <- c(
    "estimand_id", "target_source", "target_engine_id", "parameter_scope",
    "fe_names", "interval_role", "coverage_claim_allowed"
  )
  if (!identical(names(metadata), expected_fields)) {
    pv_abort("Internal BRR-Fay estimand metadata did not use its exact schema.")
  }
  metadata
}

pv_validate_brr_target_v01 <- function(target) {
  pv_assert_named_list(target, "target")
  required <- c(
    "beta", "U_bar", "B", "T_MI", "df", "df_classic", "df_method",
    "df_complete", "coverage_claim_allowed", "interval_role", "fe_names", "per_pv",
    "formula", "formula_string", "rhs_string", "M", "R", "fay_k", "design_hash",
    "target_source", "target_hash", "engine", "policy", "schema_version",
    "provenance", "warnings"
  )
  pv_required_fields(target, required, "BRR-Fay target")
  if (!is.numeric(target$beta) || any(!is.finite(target$beta)) ||
      is.null(names(target$beta)) || any(!nzchar(names(target$beta)))) {
    pv_abort("BRR-Fay target `beta` must be finite and named.")
  }
  if (!identical(target$fe_names, names(target$beta))) {
    pv_abort("BRR-Fay target `fe_names` must match `beta` names.")
  }
  fe_names <- pv_validate_unique_character(target$fe_names, "fe_names")
  if (!inherits(target$formula, "formula")) {
    pv_abort("BRR-Fay target `formula` must be a formula.")
  }
  rhs <- pv_brr_target_formula_rhs(target$formula)
  if (!identical(target$formula_string, pv_formula_string(target$formula))) {
    pv_abort("BRR-Fay target `formula_string` must match `formula`.")
  }
  if (!identical(target$rhs_string, pv_deparse_expr(rhs))) {
    pv_abort("BRR-Fay target `rhs_string` must match the formula RHS.")
  }
  pv_validate_named_numeric(target$df, "df", fe_names, allow_infinite = TRUE)
  pv_validate_named_numeric(target$df_classic, "df_classic", fe_names, allow_infinite = TRUE)
  target$df_method <- pv_assert_scalar_string(target$df_method, "df_method")
  if (!target$df_method %in% c("classic", "barnard_rubin")) {
    pv_abort("BRR-Fay target `df_method` must be `classic` or `barnard_rubin`.")
  }
  pv_validate_df_complete_field(target$df_complete, fe_names, target$df_method)
  expected_df_policy <- pv_df_policy(target$df_method)
  if (!identical(target$coverage_claim_allowed, expected_df_policy$coverage_claim_allowed)) {
    pv_abort("BRR-Fay target `coverage_claim_allowed` is inconsistent with `df_method`.")
  }
  if (!identical(target$interval_role, expected_df_policy$interval_role)) {
    pv_abort("BRR-Fay target `interval_role` is inconsistent with `df_method`.")
  }
  pv_validate_aligned_symmetric_matrix(target$U_bar, "U_bar", fe_names)
  pv_validate_aligned_symmetric_matrix(target$B, "B", fe_names)
  pv_validate_target_matrix(target$T_MI, target$fe_names)
  target$M <- pv_assert_scalar_number(target$M, "M", integer = TRUE, lower = 1)
  target$R <- pv_assert_scalar_number(target$R, "R", integer = TRUE, lower = 2)
  pv_validate_fay_k(target$fay_k)
  if (!is.list(target$per_pv) || length(target$per_pv) != target$M) {
    pv_abort("BRR-Fay target `per_pv` length must equal `M`.")
  }
  if (!identical(target$target_source, "external_brr_fay_rubin")) {
    pv_abort("BRR-Fay target source must be `external_brr_fay_rubin`.")
  }
  for (m in seq_along(target$per_pv)) {
    item <- target$per_pv[[m]]
    pv_assert_named_list(item, sprintf("per_pv[[%d]]", m))
    item_required <- c("beta", "U", "fe_names", "pv_col", "R", "fay_k")
    item_missing <- setdiff(item_required, names(item))
    if (length(item_missing) > 0L) {
      pv_abort(sprintf("BRR-Fay target `per_pv[[%d]]` is missing required field(s): %s.", m, paste(item_missing, collapse = ", ")))
    }
    if (!identical(item$fe_names, target$fe_names)) {
      pv_abort("Per-PV fixed-effect names must align with target `fe_names`.")
    }
    pv_validate_named_numeric(item$beta, sprintf("per_pv[[%d]]$beta", m), fe_names)
    pv_validate_aligned_symmetric_matrix(item$U, sprintf("per_pv[[%d]]$U", m), fe_names)
    if (!identical(item$R, target$R)) {
      pv_abort("Per-PV replicate-weight counts must align with target `R`.")
    }
    if (!identical(item$fay_k, target$fay_k)) {
      pv_abort("Per-PV Fay coefficients must align with target `fay_k`.")
    }
  }
  pv_validate_hash_scalar(target$design_hash, "design_hash")
  pv_validate_hash_scalar(target$target_hash, "target_hash")
  pv_assert_scalar_string(target$engine, "engine")
  pv_assert_named_list(target$policy, "policy")
  if (!isTRUE(target$policy$fixed_effects_only)) {
    pv_abort("BRR-Fay target policy must mark fixed effects as the only reportable scope.")
  }
  if (!identical(target$policy$target_repair, "forbidden")) {
    pv_abort("BRR-Fay target policy must mark target repair as forbidden.")
  }
  if (!identical(target$policy$df_method, target$df_method) ||
      !identical(target$policy$interval_role, target$interval_role) ||
      !identical(target$policy$coverage_claim_allowed, target$coverage_claim_allowed)) {
    pv_abort("BRR-Fay target policy df metadata (df_method, interval_role, coverage_claim_allowed) must match top-level target fields.")
  }
  pv_validate_schema_version(target$schema_version)
  pv_validate_named_list_field(target$provenance, "provenance")
  pv_validate_character_field(target$warnings, "warnings")
  invisible(target)
}

pv_brr_target_v02_fields <- function() {
  c(
    "beta", "beta_bar", "U_bar", "B", "T_MI", "total_var", "se",
    "df", "df_classic", "df_method", "df_complete", "conf_level",
    "interval_role", "coverage_claim_allowed", "lambda", "fmi", "riv",
    "fe_names", "per_pv", "formula", "formula_string", "rhs_string",
    "M", "R", "fay_k", "fay_variance_multiplier", "pv_cols",
    "weight_col", "rep_weight_cols", "id_cols", "design_hash",
    "target_source", "target_hash", "engine", "policy", "schema_version",
    "provenance", "binding_manifest", "target_content", "warnings"
  )
}

pv_brr_target_v02_per_pv_fields <- function() {
  c(
    "beta", "U", "fe_names", "pv_col", "R", "fay_k",
    "fay_variance_multiplier", "replicate_beta", "replicate_diff"
  )
}

pv_brr_target_v02_provenance_fields <- function() {
  c("function_name", "assembled_at", "package", "schema_version")
}

pv_brr_target_assert_exact_list <- function(x, fields, label) {
  if (!is.list(x) || !identical(names(x), fields) ||
      !identical(attributes(x), list(names = fields))) {
    pv_abort(sprintf("BRR-Fay target `%s` must use the exact ordered schema-0.2 fields.", label))
  }
  invisible(x)
}

pv_brr_target_require_identical <- function(expected, observed, label) {
  if (!identical(expected, observed)) {
    pv_abort(sprintf("BRR-Fay target `%s` must exactly match its canonical schema-0.2 cross-link.", label))
  }
  invisible(TRUE)
}

pv_brr_target_is_exact_matrix <- function(x, dimensions, dimnames_expected) {
  is.matrix(x) && is.numeric(x) && all(is.finite(x)) &&
    identical(dim(x), dimensions) && identical(dimnames(x), dimnames_expected) &&
    identical(
      attributes(x),
      list(dim = dimensions, dimnames = dimnames_expected)
    )
}

pv_validate_brr_target_v02 <- function(target) {
  fields <- pv_brr_target_v02_fields()
  root_attributes <- attributes(target)
  if (!is.list(target) || !identical(names(target), fields) ||
      !identical(names(root_attributes), c("names", "class")) ||
      !identical(root_attributes$names, fields) ||
      !identical(root_attributes$class, c("pvstackr_brr_target", "list"))) {
    pv_abort("BRR-Fay target is missing required field(s), has extras, or does not use exact schema-0.2 root order and class.")
  }
  if (!identical(target$schema_version, "0.2.0")) {
    pv_abort("BRR-Fay target schema version must be 0.2.0.")
  }
  canonical_formula <- tryCatch(
    stats::as.formula(target$formula_string, env = baseenv()),
    error = function(error) NULL
  )
  if (!inherits(target$formula, "formula") || is.null(canonical_formula) ||
      !identical(target$formula, canonical_formula)) {
    pv_abort("BRR-Fay target `formula` must be the exact base-environment projection of `formula_string`.")
  }
  if (!is.character(target$formula_string) || length(target$formula_string) != 1L ||
      is.na(target$formula_string) || !is.null(attributes(target$formula_string))) {
    pv_abort("BRR-Fay target `formula_string` must be an attribute-free scalar string.")
  }
  if (!identical(target$formula_string, pv_formula_string(target$formula))) {
    pv_abort("BRR-Fay target `formula_string` must be the canonical deparse of `formula`.")
  }
  rhs <- pv_brr_target_formula_rhs(target$formula)
  if (!identical(target$rhs_string, pv_deparse_expr(rhs))) {
    pv_abort("BRR-Fay target `rhs_string` must match the formula RHS.")
  }
  target$conf_level <- pv_assert_probability(target$conf_level, "conf_level")
  target$M <- pv_assert_scalar_number(target$M, "M", integer = TRUE, lower = 1)
  target$R <- pv_assert_scalar_number(target$R, "R", integer = TRUE, lower = 2)
  pv_validate_fay_k(target$fay_k)
  if (!identical(target$engine, "lm") ||
      !identical(target$target_source, "external_brr_fay_rubin")) {
    pv_abort("BRR-Fay target must use the canonical lm external_brr_fay_rubin engine.")
  }
  fe_names <- pv_validate_unique_character(target$fe_names, "fe_names")
  pv_brr_target_assert_exact_list(
    target$policy,
    pv_binding_target_policy_fields(),
    "policy"
  )
  if (!identical(target$policy$target_repair, "forbidden")) {
    pv_abort("BRR-Fay target policy must mark target repair as forbidden.")
  }
  if (!identical(target$policy$replicate_weight_role, "external_design_variance_only") ||
      !identical(target$policy$fixed_effects_only, TRUE)) {
    pv_abort("BRR-Fay target policy must preserve external replicate weights and fixed effects only.")
  }
  if (!identical(target$policy$df_method, target$df_method) ||
      !identical(target$policy$interval_role, target$interval_role) ||
      !identical(target$policy$coverage_claim_allowed, target$coverage_claim_allowed)) {
    pv_abort("BRR-Fay target policy df metadata (df_method, interval_role, coverage_claim_allowed) must match top-level target fields.")
  }
  if (identical(target$df_method, "classic") && !is.null(target$df_complete)) {
    pv_abort("BRR-Fay target classic df_method requires df_complete = NULL in schema 0.2.0.")
  }
  pv_brr_target_require_identical(
    target$policy,
    target$target_content$target_policy,
    "policy"
  )
  pv_brr_target_assert_exact_list(
    target$provenance,
    pv_brr_target_v02_provenance_fields(),
    "provenance"
  )
  if (!identical(target$provenance$function_name, "pv_brr_target") ||
      !identical(target$provenance$package, "pvstackr") ||
      !identical(target$provenance$schema_version, target$schema_version)) {
    pv_abort("BRR-Fay target provenance must use the canonical pv_brr_target schema identity.")
  }
  tryCatch(
    pv_binding_validate_optional_created_at(target$provenance$assembled_at),
    error = function(error) pv_abort("BRR-Fay target provenance timestamp must be canonical UTC ISO-8601.")
  )
  if (!identical(target$warnings, character())) {
    pv_abort("BRR-Fay target schema-0.2 warnings must be the exact empty registry subset.")
  }

  pv_binding_manifest_validate(target$binding_manifest)
  if (!"model_bundle_hash" %in% names(target$binding_manifest)) {
    pv_abort("BRR-Fay target schema-0.2 manifest must authenticate its resolved model bundle hash.")
  }
  pv_binding_target_manifest_validate(target$target_content, target$binding_manifest)
  if (!grepl("^sha256:[0-9a-f]{64}$", target$design_hash) ||
      !grepl("^sha256:[0-9a-f]{64}$", target$target_hash)) {
    pv_abort("BRR-Fay target schema-0.2 hashes must be SHA-256 cross-links; 8-character hashes are legacy-only.")
  }
  pv_brr_target_require_identical(
    target$binding_manifest$manifest_hash,
    target$design_hash,
    "design_hash"
  )
  pv_brr_target_require_identical(
    target$target_content$target_content_hash,
    target$target_hash,
    "target_hash"
  )

  components <- target$binding_manifest$components
  pv_brr_target_require_identical(components$pv$names, target$pv_cols, "pv_cols")
  pv_brr_target_require_identical(components$weight$base_name, target$weight_col, "weight_col")
  pv_brr_target_require_identical(
    components$weight$replicate_names,
    target$rep_weight_cols,
    "rep_weight_cols"
  )
  pv_brr_target_require_identical(components$row$id_cols, target$id_cols, "id_cols")
  pv_brr_target_require_identical(components$estimand$fe_names, fe_names, "fe_names")
  expected_rhs_hash <- pv_binding_hash_payload(
    pv_binding_formula_ast(target$formula)$rhs,
    "formula"
  )
  pv_brr_target_require_identical(
    components$formula$rhs_ast_hash,
    expected_rhs_hash,
    "formula RHS AST"
  )

  if (!is.list(target$per_pv) || length(target$per_pv) != target$M ||
      !is.null(attributes(target$per_pv))) {
    pv_abort("BRR-Fay target `per_pv` must be an unnamed exact list aligned with M.")
  }
  multiplier <- 1 / (target$R * (1 - target$fay_k)^2)
  pv_brr_target_require_identical(multiplier, target$fay_variance_multiplier, "fay_variance_multiplier")
  for (index in seq_along(target$per_pv)) {
    item <- target$per_pv[[index]]
    pv_brr_target_assert_exact_list(
      item,
      pv_brr_target_v02_per_pv_fields(),
      sprintf("per_pv[[%d]]", index)
    )
    if (!identical(item$fe_names, fe_names) ||
        !identical(item$pv_col, target$pv_cols[[index]]) ||
        !identical(item$R, target$R) ||
        !identical(item$fay_k, target$fay_k) ||
        !identical(item$fay_variance_multiplier, multiplier)) {
      pv_abort("BRR-Fay target per-PV identity fields must align exactly with top-level primitives.")
    }
    pv_validate_named_numeric(item$beta, sprintf("per_pv[[%d]]$beta", index), fe_names)
    pv_validate_aligned_symmetric_matrix(item$U, sprintf("per_pv[[%d]]$U", index), fe_names)
    replicate_dimensions <- c(length(fe_names), target$R)
    replicate_dimnames <- list(fe_names, target$rep_weight_cols)
    if (!pv_brr_target_is_exact_matrix(
      item$replicate_diff,
      replicate_dimensions,
      replicate_dimnames
    ) || !pv_brr_target_is_exact_matrix(
      item$replicate_beta,
      replicate_dimensions,
      replicate_dimnames
    )) {
      pv_abort("BRR-Fay target per-PV replicate differences must align with fixed effects and replicate weights.")
    }
    expected_U <- pv_symmetrize(multiplier * tcrossprod(item$replicate_diff))
    dimnames(expected_U) <- list(fe_names, fe_names)
    expected_replicate_beta <- sweep(item$replicate_diff, 1L, item$beta, FUN = "+")
    if (!isTRUE(all.equal(expected_U, item$U, tolerance = 1e-12)) ||
        !isTRUE(all.equal(expected_replicate_beta, item$replicate_beta, tolerance = 1e-12))) {
      pv_abort("BRR-Fay target per-PV replicate topology is not canonical.")
    }
    content_item <- target$target_content$per_pv[[index]]
    pv_brr_target_require_identical(as.integer(index), content_item$pv_index, "target_content pv_index")
    pv_brr_target_require_identical(item$pv_col, content_item$pv_col, "target_content pv_col")
    pv_brr_target_require_identical(item$beta, content_item$beta, "target_content beta")
    pv_brr_target_require_identical(item$U, content_item$U, "target_content U")
  }

  content <- target$target_content
  primitive_links <- list(
    M = content$M,
    R = content$R,
    fay_k = content$fay_k,
    fe_names = content$fe_names,
    df_method = content$interval_policy$df_method,
    df_complete = content$interval_policy$df_complete,
    conf_level = content$interval_policy$conf_level,
    interval_role = content$interval_policy$interval_role,
    coverage_claim_allowed = content$interval_policy$coverage_claim_allowed,
    target_source = content$target_source
  )
  for (field in names(primitive_links)) {
    pv_brr_target_require_identical(primitive_links[[field]], target[[field]], field)
  }
  for (field in pv_binding_target_derived_fields()) {
    pv_brr_target_require_identical(content$derived[[field]], target[[field]], field)
  }
  pv_brr_target_require_identical(target$beta, target$beta_bar, "beta_bar")

  beta_rows <- do.call(rbind, lapply(target$per_pv, `[[`, "beta"))
  colnames(beta_rows) <- fe_names
  pooled <- rubin_pool_matrix(
    beta = beta_rows,
    U = lapply(target$per_pv, `[[`, "U"),
    orientation = "rows_pv",
    conf_level = target$conf_level,
    allow_m1 = target$M == 1L,
    df_method = target$df_method,
    df_complete = target$df_complete
  )
  for (field in c("total_var", "lambda", "fmi", "riv")) {
    pv_brr_target_require_identical(pooled[[field]], target[[field]], field)
  }
  invisible(target)
}

pv_validate_brr_target <- function(target) {
  pv_assert_named_list(target, "target")
  schema_version <- target$schema_version
  if (identical(schema_version, "0.2.0")) {
    return(pv_validate_brr_target_v02(target))
  }
  if (identical(schema_version, "0.1.0")) {
    return(pv_validate_brr_target_v01(target))
  }
  pv_abort("BRR-Fay target `schema_version` must be 0.1.0 or 0.2.0.")
}

#' Assemble a Rubin/BRR-Fay Fixed-Effect Target
#'
#' `pv_brr_target()` assembles the external, design-based fixed-effect covariance
#' target that `stack_direct` calibrates against. For each plausible value it
#' forms a BRR-Fay sandwich covariance from the supplied replicate weights, then
#' Rubin-combines the per-plausible-value estimates and covariances across
#' plausible values into a single fixed-effect mean, total covariance, standard
#' errors, and degrees of freedom. The target engine is a dependency-free
#' weighted least-squares fit; in this package stage the target is
#' fixed-effect-only.
#'
#' @details
#' The construction follows the canonical equation set (per-plausible-value
#' BRR-Fay within-covariance, then Rubin combining).
#'
#' ## Per-plausible-value BRR-Fay sandwich
#' For plausible value \eqn{m}, the design-based within-covariance is the
#' replicate sandwich
#' \deqn{\hat U_m = a_d \sum_{r=1}^{R}
#'   (\hat\beta_m^{(r)} - \hat\beta_m)(\hat\beta_m^{(r)} - \hat\beta_m)^\top,
#'   \qquad a_d = \frac{1}{R\,(1 - k)^2},}
#' where \eqn{\hat\beta_m} is the full-sample weighted fixed-effect estimate,
#' \eqn{\hat\beta_m^{(r)}} is the estimate under replicate weight \eqn{r}, \eqn{R}
#' is the replicate count, and \eqn{k} is the Fay coefficient (`fay_k`). The
#' multiplier \eqn{a_d} is returned as `fay_variance_multiplier`.
#'
#' ## Rubin combining across plausible values
#' The per-plausible-value estimates and covariances are combined into the Rubin
#' mean \eqn{\bar\beta}, the within-imputation covariance \eqn{\bar U}, and the
#' between-imputation covariance \eqn{B},
#' \deqn{\bar\beta = \frac{1}{M}\sum_{m=1}^{M}\hat\beta_m, \qquad
#'   \bar U = \frac{1}{M}\sum_{m=1}^{M}\hat U_m, \qquad
#'   B = \frac{1}{M-1}\sum_{m=1}^{M}
#'   (\hat\beta_m - \bar\beta)(\hat\beta_m - \bar\beta)^\top,}
#' giving the Rubin total covariance --- the external target ---
#' \deqn{T_{\mathrm{MI}} = \bar U + \left(1 + \tfrac{1}{M}\right) B.}
#' Standard errors are \eqn{\sqrt{\mathrm{diag}(T_{\mathrm{MI}})}}. The fraction
#' of missing information \eqn{\gamma_k = (1 + 1/M)\,B_{kk}/T_{\mathrm{MI},kk}}
#' and the related Rubin quantities are returned as `fmi`, `riv`, and `lambda`.
#'
#' ## Degrees of freedom
#' The default `df_method = "classic"` uses the classic Rubin imputation degrees
#' of freedom and is descriptive only. `df_method = "barnard_rubin"` uses the
#' Barnard-Rubin small-sample degrees of freedom and requires `df_complete` (the
#' complete-data degrees of freedom); `df_classic` is retained alongside the
#' active `df` for reference.
#'
#' ## Computation
#' The target requires on the order of \eqn{M\,(R+1)} weighted fixed-effect fits:
#' one full-sample fit plus \eqn{R} replicate fits for each of the \eqn{M}
#' plausible values. This is a description of the computational work, not a
#' benchmarked performance claim.
#'
#' @section Interval metadata:
#' The target fixes the interval-policy fields that downstream fits carry on
#' reportable rows (consistent with [pvstackr_object_contracts]):
#' \describe{
#'   \item{`interval_role`}{`"descriptive_classic_rubin"` under the default
#'     `df_method = "classic"`, or `"coverage_barnard_rubin"` under
#'     `df_method = "barnard_rubin"` with `df_complete`.}
#'   \item{`coverage_claim_allowed`}{`FALSE` for the classic-df target;
#'     `TRUE` only under `"coverage_barnard_rubin"`. A coverage claim is actually
#'     realized downstream by a `stack_direct` fit backed by this external target;
#'     this object only declares the policy.}
#'   \item{`df_method`}{`"classic"` or `"barnard_rubin"`. Selecting
#'     `"barnard_rubin"` does not by itself make a row coverage-claimable absent
#'     the `stack_direct` target provenance.}
#'   \item{`df_complete`}{Complete-data degrees of freedom recorded when
#'     `df_method = "barnard_rubin"`; otherwise unset.}
#' }
#' This is a design-based external target. Coverage is reserved for `stack_direct`
#' rows under Barnard-Rubin df (see [pv_fit_direct()]); no speed or efficiency
#' claim attaches to assembling the target.
#'
#' @param data Data frame containing the outcome, predictors, plausible-value
#'   columns, weight column, and replicate-weight columns.
#' @param formula Two-sided formula with `OUTCOME` as the left-hand-side
#'   placeholder, for example `OUTCOME ~ x + female`. The right-hand side defines
#'   the fixed-effect design; group terms such as `(1 | school)` are not accepted
#'   in this fixed-effect-only stage.
#' @param pv_cols Character vector of plausible-value column names. If `NULL`,
#'   columns are detected with [detect_pisa_pv_columns()] using `pv_prefix`,
#'   `pv_suffix`, and `expected_M`.
#' @param weight_col Character scalar naming the full-sample weight column.
#'   Required; there is no default.
#' @param rep_weight_cols Character vector of BRR replicate-weight column names
#'   (at least two, distinct from `weight_col`). If `NULL`, columns are detected
#'   with [detect_pisa_brr_replicate_weights()] using `rep_weight_prefix` and
#'   `expected_R`.
#' @param fay_k Numeric scalar Fay coefficient \eqn{k} entering
#'   \eqn{a_d = 1/(R(1-k)^2)}. Must satisfy `0 <= fay_k < 1`. Default `0.5`
#'   (the PISA convention).
#' @param pv_prefix,pv_suffix Character scalars giving the prefix and suffix used
#'   when detecting plausible-value columns. Modern PISA files use
#'   subject-suffixed plausible values such as `PV1MATH`, so set, for example,
#'   `pv_suffix = "MATH"` rather than relying on the bare `pv_suffix = ""`
#'   default with `pv_prefix = "PV"`. Ignored when `pv_cols` is supplied.
#' @param rep_weight_prefix Character scalar prefix used when detecting replicate
#'   weights. Default `"W_FSTURWT"`. Ignored when `rep_weight_cols` is supplied.
#' @param expected_M Optional integer expected plausible-value count. If supplied,
#'   it is enforced against detected or supplied `pv_cols`. Default `NULL`.
#' @param expected_R Optional integer expected replicate-weight count. If
#'   supplied, it is enforced against detected or supplied `rep_weight_cols`.
#'   Default `NULL`.
#' @param id_cols Optional character vector of row-identifier columns. If
#'   supplied, they must jointly identify unique rows. Default `NULL`.
#' @param conf_level Numeric scalar interval level in `(0, 1)`, passed through to
#'   Rubin pooling. Default `0.95`.
#' @param allow_m1 Logical scalar; whether to allow a single plausible value
#'   (`M = 1`). Default `FALSE`.
#' @param df_method Degrees-of-freedom rule for Rubin pooling, one of
#'   `"classic"` or `"barnard_rubin"`. The default `"classic"` is descriptive
#'   only; `"barnard_rubin"` requires `df_complete`.
#' @param df_complete Complete-data degrees of freedom used when
#'   `df_method = "barnard_rubin"`; a positive numeric scalar or a per-term named
#'   numeric vector. Ignored under `"classic"`. Default `NULL`.
#' @param engine Character scalar target engine. Only `"lm"` (dependency-free
#'   weighted least squares) is implemented in this package stage. Default
#'   `"lm"`.
#' @param verbose Logical scalar; whether to emit per-plausible-value progress
#'   messages. Default `FALSE`.
#'
#' @returns A `pvstackr_brr_target` object: a list carrying the external
#'   fixed-effect target and its provenance, with class `pvstackr_brr_target`.
#'   Reportable fields include:
#'   \describe{
#'     \item{`beta_bar`, `U_bar`, `B`, `T_MI`}{The Rubin mean \eqn{\bar\beta},
#'       within-imputation covariance \eqn{\bar U}, between-imputation covariance
#'       \eqn{B}, and total (target) covariance \eqn{T_{\mathrm{MI}}} over the
#'       fixed-effect block.}
#'     \item{`se`}{Target standard errors,
#'       \eqn{\sqrt{\mathrm{diag}(T_{\mathrm{MI}})}}.}
#'     \item{`df`, `df_classic`, `df_method`, `df_complete`}{Active Rubin degrees
#'       of freedom, the classic-df reference, the df rule (`"classic"` or
#'       `"barnard_rubin"`), and the complete-data df when Barnard-Rubin is used.}
#'     \item{`interval_role`, `coverage_claim_allowed`}{Interval-policy fields for
#'       downstream fits (see the Interval metadata section and
#'       [pvstackr_object_contracts]).}
#'     \item{`lambda`, `fmi`, `riv`}{Per-term Rubin missing-information
#'       quantities: `lambda` (proportion of total variance attributable to the
#'       between-imputation component), `fmi` (fraction of missing information),
#'       and `riv` (relative increase in variance).}
#'     \item{`fe_names`, `M`, `R`, `fay_k`}{Fixed-effect names, number of
#'       plausible values \eqn{M}, replicate count \eqn{R}, and the Fay
#'       coefficient \eqn{k}.}
#'     \item{`fay_variance_multiplier`}{The BRR-Fay multiplier
#'       \eqn{a_d = 1/(R(1-k)^2)}.}
#'     \item{`pv_cols`, `weight_col`, `rep_weight_cols`, `id_cols`}{Resolved
#'       column names defining the plausible values, full-sample weight, replicate
#'       weights, and row identifiers.}
#'     \item{`target_source`, `target_hash`, `design_hash`}{Provenance:
#'       `target_source = "external_brr_fay_rubin"`, and stable content hashes of
#'       the target and its design manifest.}
#'   }
#'
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#' design <- pv_design(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
#' )
#' target <- pv_brr_target(
#'   pisa_tiny, OUTCOME ~ x + female,
#'   pv_cols = design$pv_cols, weight_col = design$weight_col,
#'   rep_weight_cols = design$rep_weight_cols, fay_k = design$fay_k,
#'   id_cols = design$id_cols
#' )
#' target
#'
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
#'
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#'
#' Judkins, D. R. (1990). Fay's method for variance estimation. *Journal of
#' Official Statistics*, 6(3), 223-239.
#'
#' @family pvstackr-target
#' @seealso
#' Build the design with [pv_design()]; detect columns with
#' [detect_pisa_pv_columns()] and [detect_pisa_brr_replicate_weights()]. Calibrate
#' against this target with [pv_fit_direct()] or via the dispatcher [pv_fit()].
#' @export
pv_brr_target <- function(
  data,
  formula,
  pv_cols = NULL,
  weight_col,
  rep_weight_cols = NULL,
  fay_k = 0.5,
  pv_prefix = "PV",
  pv_suffix = "",
  rep_weight_prefix = "W_FSTURWT",
  expected_M = NULL,
  expected_R = NULL,
  id_cols = NULL,
  conf_level = 0.95,
  allow_m1 = FALSE,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL,
  engine = "lm",
  verbose = FALSE
) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  engine <- pv_assert_scalar_string(engine, "engine")
  if (!identical(engine, "lm")) {
    pv_abort("Only `engine = \"lm\"` is implemented for BRR-Fay targets in this package stage.")
  }
  weight_col <- pv_assert_scalar_string(weight_col, "weight_col")
  fay_k <- pv_validate_fay_k(fay_k)
  conf_level <- pv_assert_probability(conf_level, "conf_level")
  allow_m1 <- pv_assert_scalar_logical(allow_m1, "allow_m1")
  df_method <- match.arg(df_method)
  verbose <- pv_assert_scalar_logical(verbose, "verbose")
  rhs <- pv_brr_target_formula_rhs(formula)
  formula_string <- pv_formula_string(formula)
  rhs_string <- pv_deparse_expr(rhs)

  if (is.null(pv_cols)) {
    pv_cols <- detect_pisa_pv_columns(data, prefix = pv_prefix, suffix = pv_suffix, expected_M = expected_M)
  } else if (!is.null(expected_M)) {
    expected_M <- pv_assert_scalar_number(expected_M, "expected_M", integer = TRUE, lower = 1)
    if (length(pv_cols) != expected_M) {
      pv_abort(sprintf("Expected %d plausible-value columns, received %d.", expected_M, length(pv_cols)))
    }
  }
  if (allow_m1) {
    pv_cols <- pv_validate_columns(data, pv_cols, "pv_cols")
    for (col in pv_cols) {
      if (!is.numeric(data[[col]]) || any(!is.finite(data[[col]]))) {
        pv_abort(sprintf("Plausible-value column `%s` must be finite and numeric.", col))
      }
    }
  } else {
    pv_cols <- pv_validate_pv_columns(data, pv_cols)
  }

  if (is.null(rep_weight_cols)) {
    rep_weight_cols <- detect_pisa_brr_replicate_weights(data, prefix = rep_weight_prefix, expected_R = expected_R)
  } else if (!is.null(expected_R)) {
    expected_R <- pv_assert_scalar_number(expected_R, "expected_R", integer = TRUE, lower = 1)
    if (length(rep_weight_cols) != expected_R) {
      pv_abort(sprintf("Expected %d replicate-weight columns, received %d.", expected_R, length(rep_weight_cols)))
    }
  }
  rep_weight_cols <- pv_validate_columns(data, rep_weight_cols, "rep_weight_cols")
  if (length(rep_weight_cols) < 2L) {
    pv_abort("`rep_weight_cols` must contain at least two BRR replicate weights.")
  }
  if (anyDuplicated(c(weight_col, rep_weight_cols))) {
    pv_abort("`weight_col` and `rep_weight_cols` must name distinct columns.")
  }
  pv_validate_columns(data, weight_col, "weight_col")
  id_cols <- pv_validate_id_columns(data, id_cols)
  model_bundle <- pv_binding_resolve_model_bundle(data, formula)

  per_pv <- vector("list", length(pv_cols))
  for (m in seq_along(pv_cols)) {
    if (verbose) {
      message(sprintf("pvstackr: BRR-Fay target PV %d/%d (%s)", m, length(pv_cols), pv_cols[[m]]))
    }
    per_pv[[m]] <- pv_brr_fay_one_pv(
      data = data,
      model_bundle = model_bundle,
      pv_col = pv_cols[[m]],
      weight_col = weight_col,
      rep_weight_cols = rep_weight_cols,
      fay_k = fay_k
    )
  }

  fe_names <- per_pv[[1L]]$fe_names
  if (!all(vapply(per_pv, function(x) identical(x$fe_names, fe_names), logical(1)))) {
    pv_abort("Fixed-effect names differ across plausible values.")
  }
  beta_rows <- do.call(rbind, lapply(per_pv, `[[`, "beta"))
  colnames(beta_rows) <- fe_names
  df_complete <- pv_normalize_df_complete_input(df_complete, fe_names, df_method)
  pooled <- rubin_pool_matrix(
    beta = beta_rows,
    U = lapply(per_pv, `[[`, "U"),
    orientation = "rows_pv",
    conf_level = conf_level,
    allow_m1 = allow_m1,
    df_method = df_method,
    df_complete = df_complete
  )
  df_policy <- pv_df_policy(pooled$df_method)

  design_manifest <- list(
    n = nrow(data),
    pv_cols = pv_cols,
    weight_col = weight_col,
    rep_weight_cols = rep_weight_cols,
    id_cols = id_cols,
    fay_k = fay_k,
    formula = formula_string,
    rhs = rhs_string,
    engine = engine
  )
  target_hash_payload <- list(
    beta = pooled$beta,
    U_bar = pooled$U_bar,
    B = pooled$B,
    T_MI = pooled$T_MI,
    df = pooled$df,
    df_classic = pooled$df_classic,
    df_method = pooled$df_method,
    df_complete = pooled$df_complete,
    interval_role = df_policy$interval_role,
    coverage_claim_allowed = df_policy$coverage_claim_allowed,
    design_hash = pv_hash_payload(design_manifest)
  )

  legacy_design_hash <- pv_hash_payload(design_manifest)
  legacy_target_hash <- pv_hash_payload(target_hash_payload)
  estimand_metadata <- pv_brr_target_estimand_metadata(model_bundle, df_policy)
  binding_manifest <- pv_binding_manifest_build(
    data = data,
    formula = formula,
    pv_cols = pv_cols,
    weight_col = weight_col,
    rep_weight_cols = rep_weight_cols,
    fay_k = fay_k,
    id_cols = id_cols,
    family_link = pv_binding_family_link_projection("gaussian", "identity"),
    estimand_contrast = NULL,
    estimand_metadata = estimand_metadata,
    model_bundle = model_bundle
  )
  binding_manifest <- c(
    binding_manifest,
    list(
      legacy_hashes = list(
        algorithm_id = pv_binding_legacy_hash_algorithm_id(),
        design_hash = legacy_design_hash,
        target_hash = legacy_target_hash
      ),
      model_bundle_hash = model_bundle$bundle_hash
    )
  )
  binding_manifest$manifest_hash <- pv_binding_hash_payload(
    pv_binding_manifest_hash_payload(binding_manifest),
    "manifest"
  )
  pv_binding_manifest_validate(binding_manifest)
  target_policy <- pv_binding_target_policy_projection(
    df_method = pooled$df_method,
    interval_role = df_policy$interval_role,
    coverage_claim_allowed = df_policy$coverage_claim_allowed
  )
  stored_derived <- list(
    beta = pooled$beta,
    U_bar = pooled$U_bar,
    B = pooled$B,
    T_MI = pooled$T_MI,
    se = pooled$se,
    df = pooled$df,
    df_classic = pooled$df_classic
  )
  target_content <- pv_binding_target_content_build(
    per_pv = per_pv,
    M = pooled$M,
    R = length(rep_weight_cols),
    fay_k = fay_k,
    df_method = pooled$df_method,
    df_complete = if (identical(pooled$df_method, "classic")) NULL else pooled$df_complete,
    conf_level = conf_level,
    target_source = "external_brr_fay_rubin",
    target_engine_id = "lm_wls_brr_fay_v1",
    target_policy = target_policy,
    manifest_hash = binding_manifest$manifest_hash,
    stored_derived = stored_derived
  )
  pv_binding_target_manifest_validate(target_content, binding_manifest)

  target <- list(
    beta = target_content$derived$beta,
    beta_bar = target_content$derived$beta,
    U_bar = target_content$derived$U_bar,
    B = target_content$derived$B,
    T_MI = target_content$derived$T_MI,
    total_var = pooled$total_var,
    se = target_content$derived$se,
    df = target_content$derived$df,
    df_classic = target_content$derived$df_classic,
    df_method = target_content$interval_policy$df_method,
    df_complete = target_content$interval_policy$df_complete,
    conf_level = target_content$interval_policy$conf_level,
    interval_role = target_content$interval_policy$interval_role,
    coverage_claim_allowed = target_content$interval_policy$coverage_claim_allowed,
    lambda = pooled$lambda,
    fmi = pooled$fmi,
    riv = pooled$riv,
    fe_names = fe_names,
    per_pv = per_pv,
    formula = stats::as.formula(formula_string, env = baseenv()),
    formula_string = formula_string,
    rhs_string = rhs_string,
    M = pooled$M,
    R = length(rep_weight_cols),
    fay_k = fay_k,
    fay_variance_multiplier = 1 / (length(rep_weight_cols) * (1 - fay_k)^2),
    pv_cols = pv_cols,
    weight_col = weight_col,
    rep_weight_cols = rep_weight_cols,
    id_cols = id_cols,
    design_hash = binding_manifest$manifest_hash,
    target_source = "external_brr_fay_rubin",
    target_hash = target_content$target_content_hash,
    engine = engine,
    policy = target_content$target_policy,
    schema_version = "0.2.0",
    provenance = list(
      function_name = "pv_brr_target",
      assembled_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      package = "pvstackr",
      schema_version = "0.2.0"
    ),
    binding_manifest = binding_manifest,
    target_content = target_content,
    warnings = character()
  )
  class(target) <- c("pvstackr_brr_target", "list")
  pv_validate_brr_target(target)
}

#' @rdname pv_brr_target
#' @param x A `pvstackr_brr_target` object from [pv_brr_target()].
#' @param ... Ignored.
#' @export
print.pvstackr_brr_target <- function(x, ...) {
  cat("pvstackr BRR-Fay target\n")
  cat("  fixed effects: ", length(x$fe_names), "\n", sep = "")
  cat("  plausible values: ", x$M, "\n", sep = "")
  cat("  replicate weights: ", x$R, "\n", sep = "")
  cat("  fay_k: ", x$fay_k, "\n", sep = "")
  cat("  df method: ", x$df_method, "\n", sep = "")
  cat("  interval role: ", x$interval_role, "\n", sep = "")
  cat("  source: ", x$target_source, "\n", sep = "")
  invisible(x)
}
