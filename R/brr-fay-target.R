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

pv_wls_beta <- function(data, formula, weights) {
  frame <- tryCatch(
    stats::model.frame(formula, data = data, na.action = stats::na.fail),
    error = function(e) pv_abort(sprintf("Model frame contains missing or invalid values: %s", conditionMessage(e)))
  )
  y <- stats::model.response(frame)
  x <- stats::model.matrix(stats::terms(formula), data = frame)
  weights <- pv_validate_weight_vector(weights, "weights", nrow(frame))
  fit <- stats::lm.wfit(x = x, y = y, w = weights)
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

pv_brr_fay_one_pv <- function(data, formula, pv_col, weight_col, rep_weight_cols, fay_k) {
  n <- nrow(data)
  pv_validate_columns(data, c(pv_col, weight_col, rep_weight_cols), "BRR-Fay input columns")
  if (!is.numeric(data[[pv_col]]) || any(!is.finite(data[[pv_col]]))) {
    pv_abort(sprintf("Plausible-value column `%s` must be finite and numeric.", pv_col))
  }
  base_w <- pv_validate_weight_vector(data[[weight_col]], weight_col, n)
  rep_w <- lapply(rep_weight_cols, function(col) {
    pv_validate_weight_vector(data[[col]], col, n)
  })
  formula_pv <- pv_formula_for_pv(formula, pv_col)
  beta0 <- pv_wls_beta(data, formula_pv, base_w)

  p <- length(beta0)
  R <- length(rep_weight_cols)
  diff <- matrix(NA_real_, nrow = p, ncol = R)
  for (r in seq_len(R)) {
    beta_r <- pv_wls_beta(data, formula_pv, rep_w[[r]])
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

pv_validate_brr_target <- function(target) {
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
    pv_abort("BRR-Fay target policy df metadata must match top-level target fields.")
  }
  pv_validate_schema_version(target$schema_version)
  pv_validate_named_list_field(target$provenance, "provenance")
  pv_validate_character_field(target$warnings, "warnings")
  invisible(target)
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

  per_pv <- vector("list", length(pv_cols))
  for (m in seq_along(pv_cols)) {
    if (verbose) {
      message(sprintf("pvstackr: BRR-Fay target PV %d/%d (%s)", m, length(pv_cols), pv_cols[[m]]))
    }
    per_pv[[m]] <- pv_brr_fay_one_pv(data, formula, pv_cols[[m]], weight_col, rep_weight_cols, fay_k)
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

  target <- list(
    beta = pooled$beta,
    beta_bar = pooled$beta_bar,
    U_bar = pooled$U_bar,
    B = pooled$B,
    T_MI = pooled$T_MI,
    total_var = pooled$total_var,
    se = pooled$se,
    df = pooled$df,
    df_classic = pooled$df_classic,
    df_method = pooled$df_method,
    df_complete = pooled$df_complete,
    interval_role = df_policy$interval_role,
    coverage_claim_allowed = df_policy$coverage_claim_allowed,
    lambda = pooled$lambda,
    fmi = pooled$fmi,
    riv = pooled$riv,
    fe_names = fe_names,
    per_pv = per_pv,
    formula = formula,
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
    design_hash = pv_hash_payload(design_manifest),
    target_source = "external_brr_fay_rubin",
    target_hash = pv_hash_payload(target_hash_payload),
    engine = engine,
    policy = list(
      replicate_weight_role = "external_design_variance_only",
      target_repair = "forbidden",
      fixed_effects_only = TRUE,
      df_method = pooled$df_method,
      interval_role = df_policy$interval_role,
      coverage_claim_allowed = df_policy$coverage_claim_allowed
    ),
    schema_version = "0.1.0",
    provenance = list(
      function_name = "pv_brr_target",
      assembled_at = as.character(Sys.time()),
      package = "pvstackr"
    ),
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
