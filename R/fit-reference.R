pv_reference_formula_for_pv <- function(formula, pv_col) {
  rhs <- pv_formula_rhs_checked(formula)
  if (pv_formula_has_weights_call(rhs)) {
    pv_abort("Do not embed `weights()` in `formula`; pass weights through the backend or data explicitly.")
  }
  stats::as.formula(
    paste(pv_col, pv_deparse_expr(rhs), sep = " ~ "),
    env = environment(formula)
  )
}

pv_reference_draw_source <- function(per_pv_draws, fit_function, draws_function) {
  has_draws <- !is.null(per_pv_draws)
  has_fit <- !is.null(fit_function)
  if (identical(has_draws, has_fit)) {
    pv_abort("Supply exactly one reference source: `per_pv_draws` or `fit_function`.")
  }
  if (has_fit && is.null(draws_function)) {
    pv_abort("`draws_function` is required when using an injected `fit_function`.")
  }
  if (has_fit && !is.function(fit_function)) {
    pv_abort("`fit_function` must be a function.")
  }
  if (!is.null(draws_function) && !is.function(draws_function)) {
    pv_abort("`draws_function` must be a function.")
  }
  if (has_draws) "draws" else "injected_fit"
}

pv_reference_normalize_draw_names <- function(per_pv_draws, pv_cols = NULL) {
  if (!is.list(per_pv_draws) || length(per_pv_draws) < 2L) {
    pv_abort("`per_pv_draws` must be a list with at least two plausible-value draw matrices.")
  }
  if (is.null(pv_cols)) {
    names_in <- names(per_pv_draws)
    if (is.null(names_in) || any(!nzchar(names_in))) {
      names(per_pv_draws) <- paste0("PV", seq_along(per_pv_draws))
    }
    return(per_pv_draws)
  }

  pv_cols <- pv_validate_unique_character(pv_cols, "pv_cols", min_len = 2L)
  if (length(pv_cols) != length(per_pv_draws)) {
    pv_abort("`pv_cols` length must match `per_pv_draws` length.")
  }
  names_in <- names(per_pv_draws)
  if (is.null(names_in) || all(!nzchar(names_in))) {
    names(per_pv_draws) <- pv_cols
    return(per_pv_draws)
  }
  if (any(!nzchar(names_in)) || anyDuplicated(names_in) || !setequal(names_in, pv_cols)) {
    pv_abort("Named `per_pv_draws` must align with `pv_cols`.")
  }
  per_pv_draws[pv_cols]
}

pv_reference_summarize_draws <- function(per_pv_draws, pv_cols = NULL, param_map = NULL) {
  per_pv_draws <- pv_reference_normalize_draw_names(per_pv_draws, pv_cols)
  M <- length(per_pv_draws)
  pv_cols <- names(per_pv_draws)

  fe_draws <- vector("list", M)
  U <- vector("list", M)
  beta_rows <- NULL
  fe_names <- NULL
  draw_counts <- integer(M)
  dropped <- vector("list", M)
  map_sources <- character(M)

  for (m in seq_len(M)) {
    label <- pv_cols[[m]]
    map <- pv_stack_param_map(per_pv_draws[[m]], param_map = param_map)
    draws_m <- map$draws_selected[, map$fe_idx, drop = FALSE]
    if (nrow(draws_m) < 2L) {
      pv_abort(sprintf("Per-PV draw matrix `%s` must contain at least two posterior draws.", label))
    }
    if (is.null(fe_names)) {
      fe_names <- colnames(draws_m)
      beta_rows <- matrix(NA_real_, nrow = M, ncol = length(fe_names))
      colnames(beta_rows) <- fe_names
      rownames(beta_rows) <- pv_cols
    } else if (!identical(colnames(draws_m), fe_names)) {
      pv_abort("Per-PV fixed-effect draw columns must align in name and order.")
    }

    fe_draws[[m]] <- draws_m
    draw_counts[[m]] <- nrow(draws_m)
    beta_rows[m, ] <- colMeans(draws_m)
    U_m <- stats::cov(draws_m)
    dimnames(U_m) <- list(fe_names, fe_names)
    U[[m]] <- U_m
    dropped[[m]] <- map$dropped_names
    map_sources[[m]] <- map$map_source
  }

  names(fe_draws) <- names(U) <- names(dropped) <- names(map_sources) <- pv_cols
  list(
    pv_cols = pv_cols,
    fe_names = fe_names,
    beta = beta_rows,
    U = U,
    draws = fe_draws,
    draw_counts = stats::setNames(draw_counts, pv_cols),
    dropped_draw_columns = dropped,
    map_sources = map_sources
  )
}

pv_reference_fit_draws <- function(
  data,
  formula,
  pv_cols,
  control,
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  weight_col = NULL,
  cache_dir = "cache",
  cache_stem = "pvstackr-per-pv",
  additional_args = list()
) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame when using an injected `fit_function`.")
  }
  pv_cols <- pv_validate_pv_columns(data, pv_cols)
  if (!is.null(weight_col)) {
    weight_col <- pv_assert_scalar_string(weight_col, "weight_col")
    pv_validate_columns(data, weight_col, "weight_col")
    pv_validate_weight_vector(data[[weight_col]], weight_col, nrow(data))
  }

  fits <- vector("list", length(pv_cols))
  per_pv_draws <- vector("list", length(pv_cols))
  for (m in seq_along(pv_cols)) {
    pv_col <- pv_cols[[m]]
    formula_m <- pv_reference_formula_for_pv(formula, pv_col)
    fit_args <- pv_stack_build_fit_args(
      prepared = list(formula = formula_m, data = data),
      family = family,
      prior = prior,
      chains = control$chains,
      iter = control$iter,
      warmup = control$warmup,
      cores = control$cores,
      seed = control$seed,
      backend = control$backend,
      cache_dir = cache_dir,
      cache_stem = paste0(cache_stem, "-", pv_col),
      additional_args = additional_args
    )
    fit_m <- do.call(fit_function, fit_args)
    fits[[m]] <- fit_m
    per_pv_draws[[m]] <- draws_function(fit_m)
  }
  names(fits) <- names(per_pv_draws) <- pv_cols
  list(fits = fits, draws = per_pv_draws)
}

pv_reference_pool_hash <- function(summary, pool) {
  pv_hash_payload(list(
    pv_cols = summary$pv_cols,
    fe_names = summary$fe_names,
    beta = summary$beta,
    U = summary$U,
    df_method = pool$df_method,
    df_complete = pool$df_complete,
    conf_level = pool$conf_level,
    df = pool$df,
    df_classic = pool$df_classic,
    pooled_beta = pool$beta,
    T_MI = pool$T_MI
  ))
}

new_pvstackr_reference_pool <- function(summary, pool, pooling_hash) {
  fe_names <- names(pool$beta)
  target <- list(
    beta = pool$beta,
    beta_bar = pool$beta_bar,
    U_bar = pool$U_bar,
    B = pool$B,
    T_MI = pool$T_MI,
    total_var = pool$total_var,
    se = pool$se,
    df = pool$df,
    df_classic = pool$df_classic,
    df_method = pool$df_method,
    df_complete = pool$df_complete,
    lambda = pool$lambda,
    fmi = pool$fmi,
    riv = pool$riv,
    rho = pool$rho,
    ci_low = pool$ci_low,
    ci_high = pool$ci_high,
    conf_level = pool$conf_level,
    M = pool$M,
    p = pool$p,
    pv_cols = summary$pv_cols,
    fe_names = fe_names,
    per_pv = list(
      beta = summary$beta,
      U = summary$U,
      draw_counts = summary$draw_counts,
      dropped_draw_columns = summary$dropped_draw_columns,
      map_sources = summary$map_sources
    ),
    target_source = "per_pv_rubin_draws",
    target_hash = pooling_hash,
    engine = "rubin_pool_matrix",
    policy = list(
      fixed_effects_only = TRUE,
      target_repair = "not_applicable",
      reportable_parameter_scope = "fixed_effect"
    ),
    schema_version = pv_schema_version(),
    provenance = pv_provenance(
      "new_pvstackr_reference_pool",
      source = "pv_fit_reference",
      pooling_hash = pooling_hash
    )
  )
  class(target) <- c("pvstackr_reference_pool", "list")
  validate_pvstackr_reference_pool(target)
  target
}

validate_pvstackr_reference_pool <- function(target) {
  pv_assert_named_list(target, "reference_pool")
  required <- c(
    "beta", "U_bar", "B", "T_MI", "se", "df", "df_classic", "df_method",
    "df_complete", "lambda", "fmi", "riv", "rho", "ci_low", "ci_high",
    "conf_level", "M", "p", "pv_cols", "fe_names", "per_pv",
    "target_source", "target_hash", "engine", "policy", "schema_version",
    "provenance"
  )
  missing <- setdiff(required, names(target))
  if (length(missing) > 0L) {
    pv_abort(sprintf("Reference pool is missing required field(s): %s.", paste(missing, collapse = ", ")))
  }
  if (!identical(target$target_source, "per_pv_rubin_draws")) {
    pv_abort("Reference pool `target_source` must be `per_pv_rubin_draws`.")
  }
  pv_validate_hash_scalar(target$target_hash, "target_hash")
  fe_names <- pv_validate_unique_character(target$fe_names, "fe_names")
  pv_validate_named_numeric(target$beta, "beta", fe_names)
  pv_validate_named_numeric(target$se, "se", fe_names)
  pv_validate_named_numeric(target$df, "df", fe_names, allow_infinite = TRUE)
  pv_validate_named_numeric(target$df_classic, "df_classic", fe_names, allow_infinite = TRUE)
  if (!is.numeric(target$df_complete) ||
      is.null(names(target$df_complete)) ||
      !identical(names(target$df_complete), fe_names)) {
    pv_abort("Reference pool `df_complete` must be a named numeric vector aligned to `fe_names`.")
  }
  if (identical(target$df_method, "classic")) {
    if (!all(is.na(target$df_complete))) {
      pv_abort("Reference pool `df_complete` must be all NA when `df_method = \"classic\"`.")
    }
  } else if (identical(target$df_method, "barnard_rubin")) {
    if (any(is.na(target$df_complete)) || any(target$df_complete <= 0)) {
      pv_abort("Reference pool `df_complete` must be positive when `df_method = \"barnard_rubin\"`.")
    }
  } else {
    pv_abort("Reference pool `df_method` must be `classic` or `barnard_rubin`.")
  }
  pv_validate_aligned_symmetric_matrix(target$U_bar, "U_bar", fe_names)
  pv_validate_aligned_symmetric_matrix(target$B, "B", fe_names)
  pv_validate_aligned_symmetric_matrix(target$T_MI, "T_MI", fe_names)
  pv_cols <- pv_validate_unique_character(target$pv_cols, "pv_cols", min_len = 2L)
  if (!identical(target$M, length(pv_cols)) || !identical(target$p, length(fe_names))) {
    pv_abort("Reference pool dimensions must match `pv_cols` and `fe_names`.")
  }
  pv_assert_named_list(target$per_pv, "per_pv")
  per_pv_required <- c("beta", "U", "draw_counts", "dropped_draw_columns", "map_sources")
  per_pv_missing <- setdiff(per_pv_required, names(target$per_pv))
  if (length(per_pv_missing) > 0L) {
    pv_abort(sprintf("Reference pool `per_pv` is missing required field(s): %s.", paste(per_pv_missing, collapse = ", ")))
  }
  if (!is.matrix(target$per_pv$beta) ||
      !identical(rownames(target$per_pv$beta), pv_cols) ||
      !identical(colnames(target$per_pv$beta), fe_names)) {
    pv_abort("Reference pool `per_pv$beta` must be an M x p matrix aligned to PV and fixed-effect names.")
  }
  if (!is.list(target$per_pv$U) || length(target$per_pv$U) != length(pv_cols)) {
    pv_abort("Reference pool `per_pv$U` must contain one covariance matrix per plausible value.")
  }
  for (m in seq_along(target$per_pv$U)) {
    pv_validate_aligned_symmetric_matrix(target$per_pv$U[[m]], sprintf("per_pv$U[[%d]]", m), fe_names)
  }
  if (!is.numeric(target$per_pv$draw_counts) ||
      !identical(names(target$per_pv$draw_counts), pv_cols) ||
      any(target$per_pv$draw_counts < 2L)) {
    pv_abort("Reference pool `per_pv$draw_counts` must be named positive draw counts.")
  }
  pv_assert_named_list(target$policy, "policy")
  if (!identical(target$policy$fixed_effects_only, TRUE) ||
      !identical(target$policy$reportable_parameter_scope, "fixed_effect")) {
    pv_abort("Reference pool policy must be fixed-effect-only.")
  }
  pv_validate_schema_version(target$schema_version)
  pv_validate_named_list_field(target$provenance, "provenance")
  invisible(target)
}

pv_fit_reference_estimates <- function(target) {
  fe_names <- target$fe_names
  data.frame(
    term = fe_names,
    estimate = unname(target$beta[fe_names]),
    se = unname(target$se[fe_names]),
    std.error = unname(target$se[fe_names]),
    df = unname(target$df[fe_names]),
    df_method = target$df_method,
    df_complete = unname(target$df_complete[fe_names]),
    conf_level = target$conf_level,
    conf_low = unname(target$ci_low[fe_names]),
    conf_high = unname(target$ci_high[fe_names]),
    conf.low = unname(target$ci_low[fe_names]),
    conf.high = unname(target$ci_high[fe_names]),
    interval_role = if (identical(target$df_method, "barnard_rubin")) {
      "reference_barnard_rubin"
    } else {
      "reference_classic_rubin"
    },
    coverage_claim_allowed = FALSE,
    parameter_scope = "fixed_effect",
    target_source = target$target_source,
    target_hash = target$target_hash,
    pooling_source = target$target_source,
    pooling_hash = target$target_hash,
    stringsAsFactors = FALSE
  )
}

#' Combine one fit per plausible value with Rubin's rules
#'
#' `pv_fit_reference()` runs the `"per_pv"` method of [pv_fit()]. It takes
#' the posterior draws of one model fit per plausible value, either computed
#' elsewhere (`per_pv_draws`) or returned by your own fitting functions
#' (`fit_function` and `draws_function`, called once per plausible value),
#' and combines the fixed effects with Rubin's rules. The variance within
#' each plausible value is the covariance of that fit's posterior draws, so
#' the intervals are always descriptive.
#'
#' @details
#' ## Calculation
#'
#' For each plausible value \eqn{m}, pvstackr takes the fixed-effect columns
#' of the draws (see `param_map`) and computes their mean
#' \eqn{\hat\beta_m}{b_m} and their covariance matrix \eqn{U_m}. Rubin's
#' rules then combine the \eqn{M} results as in [pv_brr_target()]: the
#' estimate \eqn{\bar\beta}{beta_bar} is the mean of the
#' \eqn{\hat\beta_m}{b_m}, and the total covariance is
#' \eqn{T_{\mathrm{MI}} = \bar U + (1 + 1/M) B}{T_MI = U_bar + (1 + 1/M) B},
#' where \eqn{\bar U}{U_bar} is the mean of the \eqn{U_m} and \eqn{B} is the
#' covariance matrix of the \eqn{\hat\beta_m}{b_m}. The standard errors are
#' the square roots of the diagonal of \eqn{T_{\mathrm{MI}}}{T_MI}, and the
#' degrees of freedom follow `df_method`.
#'
#' ## What is reported
#'
#' As for every method, only the fixed effects are reported ([pv_fit()]
#' gives the scope and the full rule for the intervals). Each \eqn{U_m} is
#' the model-based posterior covariance of the draws, not a design-based
#' variance from replicate weights, so every row of the estimate table has
#' `interval_role = "reference_classic_rubin"` or
#' `"reference_barnard_rubin"` and `coverage_claim_allowed = FALSE`. The
#' replicate weights play no part in this method.
#'
#' ## Fitting functions and status
#'
#' pvstackr itself fits no model for this method, and there is no bundled
#' engine: `pv_control(backend = "brms")` does not supply one. With
#' `fit_function`, pvstackr replaces `OUTCOME` in `formula` by each
#' plausible-value column in turn (for example `PV1READ ~ x + female`),
#' calls `fit_function` once for each, and passes each result to
#' `draws_function`. It does not check the model further: random-effect
#' terms such as `(1 | school)`, `family` and `prior` reach `fit_function`
#' unchanged. The survey weights are not passed as an argument to
#' `fit_function`; a weighted fit must take them from `data` inside that
#' function or from `additional_args`.
#'
#' pvstackr collects no sampler diagnostics (R-hat, effective sample sizes,
#' divergent transitions) from these fits and has no other check for this
#' method, so the status of a `per_pv` fit is always `"ok"`. Check the
#' convergence of each fit yourself.
#'
#' @param data A data frame with the plausible-value columns and the
#'   variables in `formula`, or `NULL` (default). Required with
#'   `fit_function`, which receives it unchanged. With `per_pv_draws` it is
#'   optional: given together with `formula`, it is used only to record the
#'   design (the fit's `design`), and the plausible-value names must then be
#'   columns of `data`.
#' @param formula A two-sided formula with the placeholder `OUTCOME` on the
#'   left-hand side, such as `OUTCOME ~ x + female` (see [pv_fit()]), or
#'   `NULL` (default). Required with `fit_function`. A `weights()` term stops
#'   the function with an error; all other terms, random-effect terms
#'   included, are passed to `fit_function` as written.
#' @param pv_cols A character vector with the names of the plausible-value
#'   columns, at least two, or `NULL` (default). Required with
#'   `fit_function`: the columns of `data` that are fitted one at a time.
#'   With `per_pv_draws` it is optional and names the elements of the list;
#'   if the list has names, they must be the same. Without `pv_cols`, the
#'   names of `per_pv_draws` are used, or `PV1`, `PV2`, ... for an unnamed
#'   list.
#' @param per_pv_draws A list of draws computed elsewhere, one element per
#'   plausible value (at least two), or `NULL` (default). Each element is a
#'   numeric matrix or data frame with one row per posterior draw (at least
#'   two), finite values and unique column names. The fixed-effect columns
#'   (see `param_map`) must have the same names, in the same order, in every
#'   element. Give either `per_pv_draws` or `fit_function`, not both.
#' @param control A [pv_control()] object with `method = "per_pv"`. The
#'   default is `pv_control(method = "per_pv")`, and `NULL` gives the same.
#'   Its `conf_level` sets the interval level; its sampler settings and
#'   `backend` are only passed on to `fit_function`.
#' @param family,prior Passed unchanged to `fit_function`; pvstackr does not
#'   check them. Default `NULL`.
#' @param fit_function Your function that fits the model to one plausible
#'   value, or `NULL` (default) when you give `per_pv_draws`. It is called
#'   once per plausible value with the arguments `formula` (with `OUTCOME`
#'   replaced by that plausible-value column), `data`, `family`, `prior`,
#'   `chains`, `iter`, `warmup`, `cores`, `seed`, `backend`, `file` and
#'   `file_refit` (the settings from `control` and the cache file; see
#'   `cache_dir`), plus the elements of `additional_args`. There is no
#'   `weights` argument (see Details). It may return any object that
#'   `draws_function` accepts.
#' @param draws_function Your function that takes the object returned by
#'   `fit_function` and returns its draws, in the form described for
#'   `per_pv_draws`. Required with `fit_function`.
#' @param param_map `NULL` (default) or a named list that says which columns
#'   of the draws are the fixed effects, by name (`fe_names`) or by position
#'   (`fe_idx`), as for [pv_fit_direct()]. With `NULL`, the columns whose
#'   names start with `b_` are the fixed effects. Only the fixed effects are
#'   used, and their column names become the terms of the estimate table, so
#'   they must start with `b_`. Give a list to leave out other columns that
#'   start with `b_`, such as distributional parameters named `b_sigma_*`.
#' @param weight_col The name of the final weight column, or `NULL`
#'   (default). When `data` is used, pvstackr checks that the weights are
#'   positive and records the column in the fit's `design`, but it does not
#'   use the weights: they are not passed to `fit_function`.
#' @param rep_weight_cols,fay_k,id_cols Replicate-weight columns, the Fay
#'   coefficient and row-identifier columns, as in [pv_design()]. With `data`
#'   and `formula`, they are checked and recorded in the fit's `design`; the
#'   calculation does not use them. Defaults `NULL`, `0.5` and `NULL`.
#' @param df_method The rule for the degrees of freedom: `"classic"`
#'   (default) or `"barnard_rubin"`, which needs `df_complete`. Both rules are
#'   defined in [pv_brr_target()]. With either rule the intervals are
#'   descriptive.
#' @param df_complete For `df_method = "barnard_rubin"`, the complete-data
#'   degrees of freedom, which you state: the degrees of freedom that the
#'   analysis would have if the outcome were observed directly. Give one
#'   positive number for all fixed effects, or a vector with one value per
#'   fixed effect, named by the fixed-effect columns (such as `b_Intercept`);
#'   an unnamed vector of several values is an error. Default `NULL`; giving
#'   a value with `df_method = "classic"` is an error.
#' @param cache_dir,cache_stem The folder and the start of the file name that
#'   are passed to `fit_function`; defaults `"cache"` and
#'   `"pvstackr-per-pv"`. For each plausible value, `fit_function` receives
#'   `file = file.path(path.expand(cache_dir), paste0(cache_stem, "-", pv))`,
#'   where `pv` is the plausible-value column (such as
#'   `"cache/pvstackr-per-pv-PV1READ"`), and `file_refit = "on_change"`; with
#'   `cache_dir = NULL` it receives `file = NULL` and `file_refit = "never"`.
#'   They are meant for the arguments of the same names of the brms function
#'   `brm()`. pvstackr does not create the folder.
#' @param additional_args A named list of further arguments passed to
#'   `fit_function` for every plausible value, `list()` by default. They may
#'   not repeat the arguments that pvstackr sets (see `fit_function`).
#'
#' @returns A `pvstackr_fit` object with `method = "per_pv"`. Read it with
#'   [get_estimates()], [get_target()] and [get_diagnostics()];
#'   [get_draws()] returns `NULL` for this method, and
#'   [pvstackr_object_contracts] describes every part of the object. The
#'   estimate table has one row per fixed effect with the combined estimate,
#'   its standard error and its degrees of freedom; besides the columns of
#'   every method, it has `pooling_source` and `pooling_hash`.
#'   [get_target()] returns the combined result, a list of class
#'   `pvstackr_reference_pool` with `beta` (the estimates), `U_bar`, `B`,
#'   `T_MI`, `se`, `df`, `fmi` and related fields as in [pv_brr_target()],
#'   and `per_pv` (the mean `beta` and the covariance `U` of the draws of
#'   each plausible value). [get_diagnostics()] returns `reference` (where
#'   the draws came from, the plausible values, the numbers of draws and,
#'   with `control$return_draws = TRUE`, the fixed-effect draws of each
#'   plausible value in `per_pv_draws`) and `pooling` (the combined
#'   quantities).
#'
#' @examples
#' # Draws computed elsewhere: one matrix per plausible value, one row per
#' # draw, with the fixed-effect columns b_Intercept and b_x. Real draws
#' # would come from a model fitted to each plausible value; these are
#' # simulated to show the input.
#' set.seed(1)
#' make_draws <- function(n, b0, bx) {
#'   cbind(
#'     b_Intercept = rnorm(n, b0, 0.5),
#'     b_x         = rnorm(n, bx, 0.2)
#'   )
#' }
#' per_pv_draws <- list(
#'   PV1READ = make_draws(200, 1.0, 0.30),
#'   PV2READ = make_draws(200, 1.1, 0.28)
#' )
#' fit_ref <- pv_fit_reference(
#'   per_pv_draws = per_pv_draws,
#'   pv_cols      = c("PV1READ", "PV2READ"),
#'   control      = pv_control(method = "per_pv")
#' )
#' fit_ref                     # method, status "ok" and interval note
#' get_estimates(fit_ref)      # interval_role "reference_classic_rubin"
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
#'
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#' @family pvstackr-fitting
#' @seealso [pv_compare_methods()] to compare fits, and
#'   [pvstackr_object_contracts] for the parts of a fit.
#' @export
pv_fit_reference <- function(
  data = NULL,
  formula = NULL,
  pv_cols = NULL,
  per_pv_draws = NULL,
  control = pv_control(method = "per_pv"),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  param_map = NULL,
  weight_col = NULL,
  rep_weight_cols = NULL,
  fay_k = 0.5,
  id_cols = NULL,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL,
  cache_dir = "cache",
  cache_stem = "pvstackr-per-pv",
  additional_args = list()
) {
  df_method <- match.arg(df_method)
  control <- if (is.null(control)) pv_control(method = "per_pv") else pv_validate_control(control)
  if (!identical(control$method, "per_pv")) {
    pv_abort("`control$method` must be `per_pv` for `pv_fit_reference()`.")
  }
  pv_validate_fit_data_retention_control(control)
  source <- pv_reference_draw_source(per_pv_draws, fit_function, draws_function)

  fits <- NULL
  if (identical(source, "injected_fit")) {
    if (is.null(formula)) {
      pv_abort("`formula` is required when using an injected `fit_function`.")
    }
    fit_out <- pv_reference_fit_draws(
      data = data,
      formula = formula,
      pv_cols = pv_cols,
      control = control,
      family = family,
      prior = prior,
      fit_function = fit_function,
      draws_function = draws_function,
      weight_col = weight_col,
      cache_dir = cache_dir,
      cache_stem = cache_stem,
      additional_args = additional_args
    )
    fits <- fit_out$fits
    per_pv_draws <- fit_out$draws
    pv_cols <- names(per_pv_draws)
  }

  summary <- pv_reference_summarize_draws(
    per_pv_draws = per_pv_draws,
    pv_cols = pv_cols,
    param_map = param_map
  )
  df_complete <- pv_normalize_df_complete_input(df_complete, summary$fe_names, df_method)
  pool <- rubin_pool_matrix(
    beta = summary$beta,
    U = summary$U,
    orientation = "rows_pv",
    conf_level = control$conf_level,
    df_method = df_method,
    df_complete = df_complete
  )
  pooling_hash <- pv_reference_pool_hash(summary, pool)
  target <- new_pvstackr_reference_pool(summary, pool, pooling_hash)
  estimates <- pv_fit_reference_estimates(target)

  design <- NULL
  if (!is.null(data) && !is.null(formula) && !is.null(summary$pv_cols)) {
    design <- new_pvstackr_design(
      data = data,
      formula = formula,
      pv_cols = summary$pv_cols,
      weight_col = weight_col,
      rep_weight_cols = rep_weight_cols,
      fay_k = fay_k,
      id_cols = id_cols,
      roles = list(
        outcome_placeholder = "OUTCOME",
        method = "per_pv"
      ),
      provenance = list(
        source = "pv_fit_reference",
        target_hash = target$target_hash
      )
    )
    design <- pv_design_canonicalize_formula(design)
    if (!isTRUE(control$keep_data)) {
      design <- pv_design_data_free_snapshot(design)
    }
  }

  new_pvstackr_fit(
    method = "per_pv",
    design = design,
    target = target,
    stack_fit = NULL,
    ccc = NULL,
    estimates = estimates,
    draws = NULL,
    diagnostics = list(
      reference = list(
        source = source,
        topology = "one_fit_per_plausible_value",
        M = target$M,
        pv_cols = target$pv_cols,
        fe_names = target$fe_names,
        draw_counts = target$per_pv$draw_counts,
        dropped_draw_columns = target$per_pv$dropped_draw_columns,
        map_sources = target$per_pv$map_sources,
        backend_fits = if (isTRUE(control$keep_backend_fit)) fits else NULL,
        per_pv_draws = if (isTRUE(control$return_draws)) summary$draws else NULL
      ),
      pooling = list(
        beta = target$beta,
        U_bar = target$U_bar,
        B = target$B,
        T_MI = target$T_MI,
        lambda = target$lambda,
        df = target$df,
        df_classic = target$df_classic,
        df_method = target$df_method,
        df_complete = target$df_complete,
        target_hash = target$target_hash,
        target_source = target$target_source
      )
    ),
    status = "ok",
    control = control,
    reason_codes = character(),
    provenance = list(
      wrapper_function = "pv_fit_reference",
      source = source,
      target_hash = target$target_hash,
      pooling_hash = target$target_hash
    ),
    warnings = character()
  )
}
