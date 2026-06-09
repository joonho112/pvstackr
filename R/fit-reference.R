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

#' Fit the Per-PV Bayesian/Backend Reference Method
#'
#' `pv_fit_reference()` implements the `per_pv` reference path: one backend fit
#' or posterior-draw source per plausible value, fixed-effect summary extraction,
#' and Rubin pooling with model-based within-PV covariance from the selected
#' fixed-effect draws.
#'
#' @details
#' This reference method accepts either already-computed per-PV posterior draws
#' via `per_pv_draws`, or an injected `fit_function`/`draws_function` pair. The
#' injected route rewrites the `OUTCOME` placeholder to each plausible-value
#' column and calls the backend once per plausible value. Current reportable
#' output is fixed-effect-only, selected from draw columns named like `b_*` or by
#' an explicit `param_map`. The within-PV variance component is the posterior
#' draw covariance of the selected fixed-effect columns in each plausible value.
#' `weight_col`, `rep_weight_cols`, `fay_k`, and `id_cols` are recorded as design
#' metadata for this path; BRR/Fay replicate weights are not used to form that
#' within-PV covariance or an external target. These weight arguments are not
#' passed automatically to an injected backend. If a backend fit should use
#' weights, handle that explicitly inside `fit_function`, in the data supplied
#' to the backend, or through backend-specific `additional_args`.
#'
#' `per_pv` interval metadata is source-specific. Classic Rubin pooling is
#' labeled `interval_role = "reference_classic_rubin"`; Barnard-Rubin pooling is
#' labeled `interval_role = "reference_barnard_rubin"` and uses the explicit
#' `df_complete` supplied by the user. Both roles set
#' `coverage_claim_allowed = FALSE` because the within-PV covariance is
#' model-based posterior-draw covariance, not an external design replicate
#' variance.
#'
#' @param data Optional analysis data frame. Required for the injected fitting
#'   route.
#' @param formula Optional two-sided formula with `OUTCOME` on the left-hand
#'   side. Required for the injected fitting route.
#' @param pv_cols Plausible-value columns. Required for the injected fitting
#'   route; optional alignment labels for `per_pv_draws`.
#' @param per_pv_draws Optional list of per-PV posterior draw matrices.
#' @param control A [pv_control()] object with `method = "per_pv"`.
#' @param family Optional backend family object passed to the injected fit
#'   function.
#' @param prior Optional backend prior object passed to the injected fit
#'   function.
#' @param fit_function Injected backend fitting function.
#' @param draws_function Function extracting posterior draws from each injected
#'   backend fit.
#' @param param_map Optional explicit draw-column map passed to the fixed-effect
#'   extraction layer. Use this when backend draw columns do not follow the
#'   package's automatic `b_*` fixed-effect naming convention. Supply
#'   `fe_names` or `fe_idx` to identify fixed-effect columns, and optional
#'   `vc_names` or `vc_idx` for nuisance variance-component columns. Use
#'   `vc_names = character()` to drop all nuisance columns, including
#'   distributional names such as `b_sigma_*`.
#' @param weight_col Optional main weight column recorded in the design object
#'   and validated when the injected fitting route is used. It is not passed as
#'   an automatic backend argument.
#' @param rep_weight_cols Optional replicate-weight columns recorded in the
#'   design object. They are not used by `per_pv` pooling and are not passed as
#'   automatic backend arguments.
#' @param fay_k Fay coefficient recorded in the design object when replicate
#'   weights are supplied.
#' @param id_cols Optional row identifier columns recorded in the design object.
#' @param df_method Rubin degrees-of-freedom method. The default `"classic"` is
#'   descriptive only; `"barnard_rubin"` requires `df_complete`.
#' @param df_complete Complete-data degrees of freedom used when
#'   `df_method = "barnard_rubin"`. Supply a positive scalar for all fixed
#'   effects or a named numeric vector aligned by fixed-effect name. Unnamed
#'   length-p vectors are rejected to avoid positional ambiguity.
#' @param cache_dir,cache_stem Cache location metadata passed to injected fits.
#' @param additional_args Additional named arguments passed to each injected fit.
#'
#' @returns A `pvstackr_fit` object with `method = "per_pv"`.
#' @section Reportable scope and coverage:
#' In this package stage, reportable output is **fixed-effect-only**; variance
#' components are fit but not calibrated to the target. Coverage claims are
#' enabled **only** for `stack_direct` rows backed by the external Rubin/BRR-Fay
#' target (`interval_role = "coverage_barnard_rubin"`,
#' `coverage_claim_allowed = TRUE`). `per_pv` intervals are
#' **descriptive/reference** even with Barnard-Rubin degrees of freedom. "One
#' stacked fit" describes the computational topology, not a benchmarked speed
#' claim. The within-PV variance is the **model-based** posterior covariance of
#' the selected fixed-effect draws, not a BRR/Fay replicate variance.
#'
#' @examples
#' # Injected per-PV draws: one posterior draw matrix per plausible value, with
#' # fixed-effect columns named like `b_*` (auto-selected). Real fits would come
#' # from a backend; here we use small synthetic draw clouds to show the shape.
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
#' fit_ref                     # a per_pv pvstackr_fit
#' get_estimates(fit_ref)      # interval_role = "reference_*"; coverage_claim_allowed = FALSE
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
#'
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#' @family pvstackr-fitting
#' @seealso [pv_fit()], [pv_fit_direct()], [pv_fit_stack_psis()],
#'   [pv_control()]; [pv_compare_methods()], [get_estimates()],
#'   [pvstackr_object_contracts].
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
