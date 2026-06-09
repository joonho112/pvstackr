#' Declare a PISA-Style Plausible-Value Design
#'
#' `pv_design()` records plausible-value columns, final weights, BRR replicate
#' weights, Fay coefficient, row-support metadata, and stable hashes without
#' running a model or assembling a target. It is the applied-user wrapper around
#' the internal `pvstackr_design` object.
#'
#' @param data Data frame containing plausible-value, weight, and replicate
#'   weight columns. Must have unique, non-empty column names.
#' @param formula Formula using `OUTCOME` as the plausible-value placeholder,
#'   for example `OUTCOME ~ x + female`. Do not embed `weights()`; pass weight
#'   columns explicitly.
#' @param pv_cols Optional character vector of plausible-value columns. If
#'   `NULL` (default), columns are detected with [detect_pisa_pv_columns()]
#'   using `pv_prefix`/`pv_suffix`.
#' @param weight_col Character scalar naming the final survey-weight column.
#'   Defaults to PISA's `"W_FSTUWT"`.
#' @param rep_weight_cols Optional character vector of BRR replicate-weight
#'   columns. If `NULL` (default), columns are detected with
#'   [detect_pisa_brr_replicate_weights()] using `rep_weight_prefix`.
#' @param fay_k Numeric scalar Fay coefficient used by the BRR replicate design.
#'   Must satisfy `0 <= fay_k < 1`. Default `0.5`.
#' @param pv_prefix,pv_suffix Character scalars giving the prefix and suffix
#'   used when detecting plausible values. Modern PISA files commonly use
#'   subject-suffixed plausible values such as `PV1MATH`; set
#'   `pv_suffix = "MATH"` or the relevant subject suffix instead of relying on
#'   the bare `PV1`, `PV2`, ... default. Defaults `"PV"` and `""`.
#' @param rep_weight_prefix Character scalar prefix used when detecting
#'   replicate weights. Default `"W_FSTURWT"`.
#' @param expected_M Optional integer scalar. Expected plausible-value count;
#'   checked against detected or supplied `pv_cols`. If `NULL` (default), the
#'   count is not checked.
#' @param expected_R Optional integer scalar. Expected replicate-weight count;
#'   checked against detected or supplied `rep_weight_cols`. If `NULL`
#'   (default), the count is not checked.
#' @param id_cols Optional character vector of columns that jointly identify
#'   unique rows. If `NULL` (default), row support is tracked by row number.
#' @param roles Optional named list of extra design roles, merged over the
#'   defaults recorded by the wrapper. Default `list()`.
#'
#' @returns A `pvstackr_design` object: a list recording the declared design and
#'   its provenance. User-facing fields include:
#'   \describe{
#'     \item{`n`, `M`, `R`}{Row count, number of plausible values, and number
#'       of BRR replicate weights.}
#'     \item{`formula`, `formula_string`}{The design formula (with the `OUTCOME`
#'       placeholder) and its string form.}
#'     \item{`pv_cols`, `weight_col`, `rep_weight_cols`}{Resolved
#'       plausible-value, final-weight, and replicate-weight column names.}
#'     \item{`fay_k`}{Fay coefficient for the BRR replicate design.}
#'     \item{`id_cols`, `row_support`}{Row-identifier columns and the
#'       row-support descriptor.}
#'     \item{`design_hash`, `row_support_hash`, `pv_value_hash`,
#'       `weight_design_hash`}{Stable content hashes of the design and its
#'       components.}
#'     \item{`roles`, `provenance`}{Design roles and wrapper provenance.}
#'   }
#'   The compact `print()` method shows `n`, `formula_string`, `M`,
#'   `weight_col`, `R`, `fay_k`, and `design_hash`.
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
#' @seealso
#' Detect inputs with [detect_pisa_pv_columns()] and
#' [detect_pisa_brr_replicate_weights()]; assemble the external target with
#' [pv_brr_target()]; fit a method with [pv_fit()].
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
#' @param x A `pvstackr_design` object.
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
