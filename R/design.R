#' Declare the plausible-value and weight columns
#'
#' `pv_design()` records which columns of `data` hold the plausible values,
#' the final weight and the replicate weights, and the Fay coefficient used to
#' make the replicate weights. It checks these columns and returns them in one
#' object; it fits no model. By default it finds PISA-style column names with
#' [detect_pisa_pv_columns()] and [detect_pisa_brr_replicate_weights()]. No
#' other pvstackr function takes the design object: pass its columns, such as
#' `design$pv_cols`, on to [pv_brr_target()], as in the examples.
#'
#' @details
#' ## PISA data
#'
#' A PISA 2022 student file has ten plausible values per subject, such as
#' `PV1READ` to `PV10READ` for reading, the final weight `W_FSTUWT` and 80
#' replicate weights `W_FSTURWT1` to `W_FSTURWT80`, made with Fay's method
#' with coefficient 0.5. The defaults `weight_col = "W_FSTUWT"`,
#' `rep_weight_prefix = "W_FSTURWT"` and `fay_k = 0.5` match these, so only
#' the subject needs to be given, as the suffix of the plausible-value names:
#'
#' ```
#' design <- pv_design(pisa, OUTCOME ~ ESCS, pv_suffix = "READ",
#'                     expected_M = 10, expected_R = 80)
#' ```
#'
#' Here `pisa` is the student data frame. `expected_M` and `expected_R` are
#' optional: they stop the function when a different number of columns is
#' found, for example because the suffix is wrong.
#'
#' ## Other data
#'
#' To choose the columns yourself, give their names in `pv_cols` and
#' `rep_weight_cols`, the final weight in `weight_col` and the Fay coefficient
#' of your replicate weights in `fay_k`. `pv_prefix`, `pv_suffix` and
#' `rep_weight_prefix` are then not used.
#'
#' ## Checks
#'
#' `pv_design()` stops with an error when a named column is not in `data`,
#' when there are fewer than two plausible-value or replicate-weight columns,
#' when a plausible value is not a finite number, when a weight is not a
#' finite positive number, when `weight_col` is also listed as a replicate
#' weight, or when `id_cols` do not identify the rows. It checks only the form
#' of `formula`; the variables in it are checked by [pv_brr_target()].
#'
#' @param data A data frame with unique, non-empty column names that contains
#'   the plausible-value and weight columns.
#' @param formula A two-sided formula with the placeholder `OUTCOME` on the
#'   left-hand side, for example `OUTCOME ~ x + female`; `OUTCOME` stands for
#'   the plausible values. Do not put `weights()` in it; name the weights in
#'   `weight_col` and `rep_weight_cols`.
#' @param pv_cols A character vector with the names of at least two
#'   plausible-value columns, or `NULL` (default) to find them with
#'   [detect_pisa_pv_columns()] from `pv_prefix` and `pv_suffix`.
#' @param weight_col The name of the final weight column. Default
#'   `"W_FSTUWT"`, the PISA name.
#' @param rep_weight_cols A character vector with the names of at least two
#'   replicate-weight columns, or `NULL` (default) to find them with
#'   [detect_pisa_brr_replicate_weights()] from `rep_weight_prefix`. The
#'   weights must be positive, so plain BRR weights, which are zero for half
#'   of the sample, are not accepted.
#' @param fay_k The Fay coefficient used to make the replicate weights, a
#'   number with `0 <= fay_k < 1`. Default `0.5`, the PISA value.
#'   [pv_brr_target()] multiplies the replicate variance by
#'   `1 / (1 - fay_k)^2`.
#' @param pv_prefix,pv_suffix The text before and after the number in the
#'   names of the plausible-value columns, used only when `pv_cols` is `NULL`.
#'   For PISA data give the subject as the suffix, for example
#'   `pv_suffix = "READ"` for `PV1READ`, `PV2READ`, and so on. The defaults,
#'   `"PV"` and `""`, match only bare names such as `PV1` and `PV2`.
#' @param rep_weight_prefix The text before the number in the names of the
#'   replicate-weight columns, used only when `rep_weight_cols` is `NULL`.
#'   Default `"W_FSTURWT"`, which matches the PISA names `W_FSTURWT1`,
#'   `W_FSTURWT2`, and so on.
#' @param expected_M,expected_R The numbers of plausible-value and
#'   replicate-weight columns you expect, as whole numbers (10 and 80 for
#'   PISA 2022). Each is compared with the columns found or given in
#'   `pv_cols` and `rep_weight_cols`, and a different number stops the
#'   function with an error. `NULL` (default) skips the check.
#' @param id_cols A character vector of columns that together identify each
#'   row, such as a student ID; their combined values must be unique and not
#'   missing. `NULL` (default) identifies the rows by their position.
#' @param roles A named list of extra entries for the `roles` field of the
#'   result; an entry with the same name as one recorded by `pv_design()`
#'   replaces it. Default `list()`.
#'
#' @returns A `pvstackr_design` object, a list with these fields:
#'   \describe{
#'     \item{`data`}{The data frame as supplied.}
#'     \item{`n`, `M`, `R`}{The numbers of rows, plausible values and
#'       replicate weights.}
#'     \item{`formula`, `formula_string`}{The formula and its text.}
#'     \item{`pv_cols`, `weight_col`, `rep_weight_cols`}{The column names, as
#'       given or found.}
#'     \item{`fay_k`, `id_cols`}{The Fay coefficient and the row-identifier
#'       columns.}
#'     \item{`row_support`}{How the rows are identified: a list with `type`
#'       (`"id_cols"` or `"row_number"`), `id_cols`, `n` and `hash`.}
#'     \item{`design_hash`, `row_support_hash`, `pv_value_hash`,
#'       `weight_design_hash`}{Checksums of the whole design, the row
#'       identifiers, the plausible values, and the weights with `fay_k`. They
#'       do not cover the values of other columns, such as covariates.}
#'     \item{`roles`, `provenance`}{How the design was made: the `OUTCOME`
#'       placeholder, the prefixes, suffix and expected counts used to find
#'       the columns, and whether the columns were found or given.}
#'     \item{`data_manifest`, `created_at`, `schema_version`, `warnings`}{The
#'       names and classes of the columns of `data` with the declared columns,
#'       the creation time, the format version of the object, and warnings,
#'       which `pv_design()` leaves empty.}
#'   }
#'   `print()` shows the number of rows, the formula, `M`, the final weight,
#'   `R`, `fay_k` and `design_hash`, and returns the object invisibly. It
#'   first checks the object and stops with an error if the plausible values,
#'   the weights, the row identifiers or the declared columns were changed
#'   after the object was created.
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#' design <- pv_design(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
#' )
#' design
#' @family pvstackr-design
#' @seealso [detect_pisa_pv_columns()] and
#'   [detect_pisa_brr_replicate_weights()] for the column search;
#'   [pv_brr_target()] for the next step.
#' @export
pv_design <- function(
  data,
  formula,
  pv_cols = NULL,
  weight_col = "W_FSTUWT",
  rep_weight_cols = NULL,
  fay_k = 0.5,
  pv_prefix = "PV",
  pv_suffix = "",
  rep_weight_prefix = "W_FSTURWT",
  expected_M = NULL,
  expected_R = NULL,
  id_cols = NULL,
  roles = list()
) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  if (is.null(names(data)) || any(!nzchar(names(data)))) {
    pv_abort("`data` must have non-empty column names.")
  }
  if (anyDuplicated(names(data))) {
    pv_abort("`data` column names must be unique.")
  }
  detected_pv_cols <- is.null(pv_cols)
  detected_rep_weight_cols <- is.null(rep_weight_cols)
  fay_k <- pv_validate_fay_k(fay_k)
  roles <- pv_validate_named_list_field(roles, "roles")

  if (is.null(pv_cols)) {
    pv_cols <- detect_pisa_pv_columns(
      data,
      prefix = pv_prefix,
      suffix = pv_suffix,
      expected_M = expected_M
    )
  } else if (!is.null(expected_M)) {
    expected_M <- pv_assert_scalar_number(expected_M, "expected_M", integer = TRUE, lower = 1)
    if (length(pv_cols) != expected_M) {
      pv_abort(sprintf("Expected %d plausible-value columns, received %d.", expected_M, length(pv_cols)))
    }
  }

  if (is.null(rep_weight_cols)) {
    rep_weight_cols <- detect_pisa_brr_replicate_weights(
      data,
      prefix = rep_weight_prefix,
      expected_R = expected_R
    )
  } else if (!is.null(expected_R)) {
    expected_R <- pv_assert_scalar_number(expected_R, "expected_R", integer = TRUE, lower = 1)
    if (length(rep_weight_cols) != expected_R) {
      pv_abort(sprintf("Expected %d replicate-weight columns, received %d.", expected_R, length(rep_weight_cols)))
    }
  }

  default_roles <- list(
    outcome_placeholder = "OUTCOME",
    helper = "pv_design",
    column_detection = list(
      pv_cols_detected = detected_pv_cols,
      rep_weight_cols_detected = detected_rep_weight_cols,
      pv_prefix = pv_prefix,
      pv_suffix = pv_suffix,
      rep_weight_prefix = rep_weight_prefix,
      expected_M = expected_M,
      expected_R = expected_R
    )
  )
  roles <- utils::modifyList(default_roles, roles)

  new_pvstackr_design(
    data = data,
    formula = formula,
    pv_cols = pv_cols,
    weight_col = weight_col,
    rep_weight_cols = rep_weight_cols,
    fay_k = fay_k,
    id_cols = id_cols,
    roles = roles,
    provenance = list(
      wrapper_function = "pv_design",
      pv_cols_detected = detected_pv_cols,
      rep_weight_cols_detected = detected_rep_weight_cols
    )
  )
}

#' @rdname pv_design
#' @param x A `pvstackr_design` object from `pv_design()`.
#' @param ... Ignored.
#' @export
print.pvstackr_design <- function(x, ...) {
  validate_pvstackr_design(x)
  cat("pvstackr design\n")
  cat("  rows: ", x$n, "\n", sep = "")
  cat("  formula: ", x$formula_string, "\n", sep = "")
  cat("  plausible values: ", x$M, "\n", sep = "")
  cat("  final weight: ", x$weight_col %||% "none", "\n", sep = "")
  cat("  replicate weights: ", x$R, "\n", sep = "")
  cat("  fay_k: ", x$fay_k, "\n", sep = "")
  cat("  design hash: ", x$design_hash, "\n", sep = "")
  invisible(x)
}
