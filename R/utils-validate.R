pv_abort <- function(message) {
  stop(sprintf("pvstackr: %s", message), call. = FALSE)
}

pv_allowed_methods <- function() {
  c("stack_direct", "stack_psis", "per_pv")
}

pv_forbidden_pipeline_methods <- function() {
  c(
    "pipeline_a",
    "pipeline_b",
    "pipeline_c",
    "pipeline_c_direct",
    "c_direct",
    "Pipeline A",
    "Pipeline B",
    "Pipeline C",
    "Pipeline C-Direct"
  )
}

pv_assert_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    pv_abort(sprintf("`%s` must be a non-missing logical scalar.", name))
  }
  x
}

pv_assert_scalar_string <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    pv_abort(sprintf("`%s` must be a non-empty character scalar.", name))
  }
  x
}

pv_assert_scalar_number <- function(
  x,
  name,
  integer = FALSE,
  lower = -Inf,
  upper = Inf,
  inclusive_lower = TRUE,
  inclusive_upper = TRUE,
  allow_null = FALSE
) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    pv_abort(sprintf("`%s` must be a finite numeric scalar.", name))
  }
  if (integer && x != floor(x)) {
    pv_abort(sprintf("`%s` must be an integer-valued scalar.", name))
  }
  lower_ok <- if (inclusive_lower) x >= lower else x > lower
  upper_ok <- if (inclusive_upper) x <= upper else x < upper
  if (!lower_ok || !upper_ok) {
    lower_op <- if (inclusive_lower) ">=" else ">"
    upper_op <- if (inclusive_upper) "<=" else "<"
    pv_abort(sprintf("`%s` must satisfy %s %s %s and %s %s %s.", name, name, lower_op, lower, name, upper_op, upper))
  }
  if (integer) {
    return(as.integer(x))
  }
  x
}

pv_assert_probability <- function(x, name) {
  pv_assert_scalar_number(
    x,
    name = name,
    lower = 0,
    upper = 1,
    inclusive_lower = FALSE,
    inclusive_upper = FALSE
  )
}

pv_assert_named_list <- function(x, name) {
  if (!is.list(x) || length(names(x)) != length(x) || any(!nzchar(names(x)))) {
    pv_abort(sprintf("`%s` must be a fully named list.", name))
  }
  x
}

pv_validate_unique_character <- function(x, name, min_len = 1L) {
  if (!is.character(x) || any(is.na(x)) || any(!nzchar(x))) {
    pv_abort(sprintf("`%s` must be a character vector with no missing or empty values.", name))
  }
  if (length(x) < min_len) {
    pv_abort(sprintf("`%s` must contain at least %d value(s).", name, min_len))
  }
  if (anyDuplicated(x)) {
    pv_abort(sprintf("`%s` must contain unique values.", name))
  }
  x
}

pv_validate_columns <- function(data, cols, label) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  cols <- pv_validate_unique_character(cols, label)
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0L) {
    pv_abort(sprintf("`%s` contains column(s) not found in `data`: %s.", label, paste(missing, collapse = ", ")))
  }
  cols
}

pv_validate_method <- function(method) {
  method <- pv_assert_scalar_string(method, "method")
  allowed <- pv_allowed_methods()
  forbidden <- pv_forbidden_pipeline_methods()

  if (method %in% forbidden || tolower(method) %in% tolower(forbidden)) {
    pv_abort(sprintf(
      "`%s` is an internal manuscript/pipeline label, not a public method. Use one of: %s.",
      method,
      paste(allowed, collapse = ", ")
    ))
  }

  if (!method %in% allowed) {
    pv_abort(sprintf("`method` must be one of: %s.", paste(allowed, collapse = ", ")))
  }

  method
}

pv_validate_backend <- function(backend) {
  backend <- pv_assert_scalar_string(backend, "backend")
  allowed <- c("none", "injected", "brms", "cmdstanr")
  if (!backend %in% allowed) {
    pv_abort(sprintf("`backend` must be one of: %s.", paste(allowed, collapse = ", ")))
  }
  backend
}

pv_validate_center <- function(center) {
  center <- pv_assert_scalar_string(center, "center")
  allowed <- c("target", "posterior")
  if (!center %in% allowed) {
    pv_abort(sprintf("`center` must be one of: %s.", paste(allowed, collapse = ", ")))
  }
  center
}

pv_expr_call_name <- function(expr) {
  if (!is.call(expr)) {
    return(NULL)
  }
  head <- expr[[1L]]
  if (is.symbol(head)) {
    return(as.character(head))
  }
  if (is.call(head) &&
      as.character(head[[1L]]) %in% c("::", ":::") &&
      length(head) == 3L &&
      is.symbol(head[[3L]])) {
    return(as.character(head[[3L]]))
  }
  NULL
}

pv_expr_is_call_named <- function(expr, name) {
  identical(pv_expr_call_name(expr), name)
}

pv_expr_contains_call <- function(expr, name) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  if (pv_expr_is_call_named(expr, name)) {
    return(TRUE)
  }
  args <- as.list(expr)[-1L]
  if (length(args) == 0L) {
    return(FALSE)
  }
  any(vapply(args, pv_expr_contains_call, logical(1), name = name))
}

pv_expr_has_bar <- function(expr, ignore_I = TRUE) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  op <- pv_expr_call_name(expr)
  if (ignore_I && identical(op, "I")) {
    return(FALSE)
  }
  if (op %in% c("|", "||")) {
    return(TRUE)
  }
  if (!op %in% c("~", "(", "+", "-", "*", "/", ":", "^")) {
    return(FALSE)
  }
  args <- as.list(expr)[-1L]
  if (length(args) == 0L) {
    return(FALSE)
  }
  any(vapply(args, pv_expr_has_bar, logical(1), ignore_I = ignore_I))
}

pv_formula_has_weights_call <- function(formula) {
  pv_expr_contains_call(formula, "weights")
}

pv_formula_has_random_effect_bar <- function(formula) {
  pv_expr_has_bar(formula, ignore_I = TRUE)
}

pv_formula_rhs_checked <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    pv_abort("`formula` must be a two-sided formula with `OUTCOME` on the left-hand side.")
  }
  if (!identical(formula[[2L]], as.name("OUTCOME"))) {
    pv_abort("`formula` must use `OUTCOME` as the left-hand-side placeholder.")
  }
  formula[[3L]]
}

pv_deparse_expr <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = "")
}

pv_formula_string <- function(formula) {
  pv_formula_rhs_checked(formula)
  pv_deparse_expr(formula)
}

pv_formula_rhs_string <- function(formula) {
  pv_deparse_expr(pv_formula_rhs_checked(formula))
}

pv_validate_target_repair_control <- function(allow_target_nearpd) {
  allow_target_nearpd <- pv_assert_scalar_logical(allow_target_nearpd, "allow_target_nearpd")
  if (allow_target_nearpd) {
    pv_abort("`allow_target_nearpd = TRUE` is reserved for future target repair and is not currently supported.")
  }
  allow_target_nearpd
}

pv_validate_fay_k <- function(fay_k) {
  pv_assert_scalar_number(
    fay_k,
    name = "fay_k",
    lower = 0,
    upper = 1,
    inclusive_lower = TRUE,
    inclusive_upper = FALSE
  )
}

pv_validate_weight_vector <- function(w, name, n, strictly_positive = TRUE) {
  if (!is.numeric(w) || length(w) != n || any(!is.finite(w))) {
    pv_abort(sprintf("`%s` must be a finite numeric vector of length %d.", name, n))
  }
  if (strictly_positive && any(w <= 0)) {
    pv_abort(sprintf("`%s` must contain strictly positive weights.", name))
  }
  if (!strictly_positive && any(w < 0)) {
    pv_abort(sprintf("`%s` must contain non-negative weights.", name))
  }
  w
}

pv_validate_pv_columns <- function(data, pv_cols) {
  pv_cols <- pv_validate_columns(data, pv_cols, "pv_cols")
  if (length(pv_cols) < 2L) {
    pv_abort("`pv_cols` must contain at least two plausible-value columns.")
  }

  for (col in pv_cols) {
    values <- data[[col]]
    if (!is.numeric(values) || any(!is.finite(values))) {
      pv_abort(sprintf("Plausible-value column `%s` must be finite and numeric.", col))
    }
  }

  pv_cols
}

pv_target_marker <- function(target) {
  values <- c(
    attr(target, "target_source", exact = TRUE),
    attr(target, "provenance", exact = TRUE),
    attr(target, "source", exact = TRUE)
  )
  values <- values[!vapply(values, is.null, logical(1))]
  if (length(values) == 0L) {
    return(character())
  }
  tolower(as.character(values))
}

pv_reject_forbidden_target_marker <- function(target) {
  marker <- pv_target_marker(target)
  if (length(marker) == 0L) {
    return(invisible(TRUE))
  }
  blocked <- c("mock", "sham", "fallback", "self", "raw posterior", "stacked posterior", "proportional")
  bad <- blocked[vapply(blocked, function(term) any(grepl(term, marker, fixed = TRUE)), logical(1))]
  if (length(bad) > 0L) {
    pv_abort(sprintf("Target covariance provenance is forbidden for confirmatory reporting: %s.", paste(unique(bad), collapse = ", ")))
  }
  invisible(TRUE)
}

pv_validate_target_matrix <- function(target, fe_names, raw_fe_cov = NULL, allow_nearpd = FALSE) {
  fe_names <- pv_validate_unique_character(fe_names, "fe_names")
  allow_nearpd <- pv_assert_scalar_logical(allow_nearpd, "allow_nearpd")
  if (allow_nearpd) {
    pv_abort("`allow_nearpd = TRUE` is reserved for future target repair and is not currently supported.")
  }
  if (!is.matrix(target) || !is.numeric(target) || length(dim(target)) != 2L) {
    pv_abort("`target` must be a numeric matrix.")
  }
  if (nrow(target) != ncol(target)) {
    pv_abort("`target` must be square.")
  }
  if (any(!is.finite(target))) {
    pv_abort("`target` must contain only finite numeric values.")
  }
  if (!identical(rownames(target), fe_names) || !identical(colnames(target), fe_names)) {
    pv_abort("`target` row and column names must exactly match `fe_names` in order.")
  }
  if (!isTRUE(all.equal(target, t(target), tolerance = 1e-10))) {
    pv_abort("`target` must be symmetric.")
  }
  pv_reject_forbidden_target_marker(target)

  eigenvalues <- eigen(target, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= 0)) {
    if (!allow_nearpd) {
      pv_abort("`target` must be positive definite; automatic repair is disabled by default.")
    }
  }

  if (!is.null(raw_fe_cov)) {
    if (!is.matrix(raw_fe_cov) || !identical(dim(raw_fe_cov), dim(target))) {
      pv_abort("`raw_fe_cov`, when supplied, must have the same dimensions as `target`.")
    }
    if (any(!is.finite(raw_fe_cov))) {
      pv_abort("`raw_fe_cov`, when supplied, must contain only finite numeric values.")
    }
    exact_scale <- max(1, max(abs(target)), max(abs(raw_fe_cov)))
    if (max(abs(target - raw_fe_cov)) <= 1e-10 * exact_scale) {
      pv_abort("`target` must not equal the raw fixed-effect covariance.")
    }
    ratio <- target / raw_fe_cov
    ratio <- ratio[is.finite(ratio) & raw_fe_cov != 0]
    if (length(ratio) > 1L && stats::sd(ratio) < 1e-10) {
      pv_abort("`target` must not be merely proportional to the raw fixed-effect covariance.")
    }
  }

  target
}

pv_validate_psis_k_threshold <- function(x, name = "psis_k_threshold") {
  pv_assert_scalar_number(
    x,
    name,
    lower = 0,
    upper = 0.7,
    inclusive_lower = FALSE,
    inclusive_upper = TRUE
  )
}

pv_validate_control <- function(control) {
  original_control <- control
  pv_assert_named_list(control, "control")
  required <- c(
    "method",
    "chains",
    "iter",
    "warmup",
    "cores",
    "seed",
    "backend",
    "conf_level",
    "psis_k_threshold",
    "center",
    "allow_target_nearpd",
    "return_draws",
    "keep_data",
    "keep_backend_fit",
    "keep_log_lik",
    "verbose"
  )
  root_attributes <- attributes(control)
  if (!identical(names(control), required) ||
      !identical(names(root_attributes), c("names", "class")) ||
      !identical(root_attributes$names, required) ||
      !identical(root_attributes$class, c("pvstackr_control", "list"))) {
    pv_abort("`control` must use the exact canonical field order, class, and attributes.")
  }
  control$method <- pv_validate_method(control$method)
  control$chains <- pv_assert_scalar_number(control$chains, "chains", integer = TRUE, lower = 1)
  control$iter <- pv_assert_scalar_number(control$iter, "iter", integer = TRUE, lower = 2)
  control$warmup <- pv_assert_scalar_number(control$warmup, "warmup", integer = TRUE, lower = 0)
  if (control$warmup >= control$iter) {
    pv_abort("`warmup` must be smaller than `iter`.")
  }
  control$cores <- pv_assert_scalar_number(control$cores, "cores", integer = TRUE, lower = 1)
  control["seed"] <- list(pv_assert_scalar_number(control$seed, "seed", integer = TRUE, lower = 0, allow_null = TRUE))
  control$backend <- pv_validate_backend(control$backend)
  control$conf_level <- pv_assert_probability(control$conf_level, "conf_level")
  control$psis_k_threshold <- pv_validate_psis_k_threshold(control$psis_k_threshold)
  control$center <- pv_validate_center(control$center)
  control$allow_target_nearpd <- pv_validate_target_repair_control(control$allow_target_nearpd)
  control$return_draws <- pv_assert_scalar_logical(control$return_draws, "return_draws")
  control$keep_data <- pv_assert_scalar_logical(control$keep_data, "keep_data")
  control$keep_backend_fit <- pv_assert_scalar_logical(control$keep_backend_fit, "keep_backend_fit")
  control$keep_log_lik <- pv_assert_scalar_logical(control$keep_log_lik, "keep_log_lik")
  control$verbose <- pv_assert_scalar_logical(control$verbose, "verbose")
  canonical <- list(
    method = unname(as.character(control$method)),
    chains = as.integer(control$chains),
    iter = as.integer(control$iter),
    warmup = as.integer(control$warmup),
    cores = as.integer(control$cores),
    seed = if (is.null(control$seed)) NULL else as.integer(control$seed),
    backend = unname(as.character(control$backend)),
    conf_level = as.numeric(control$conf_level),
    psis_k_threshold = as.numeric(control$psis_k_threshold),
    center = unname(as.character(control$center)),
    allow_target_nearpd = as.logical(control$allow_target_nearpd),
    return_draws = as.logical(control$return_draws),
    keep_data = as.logical(control$keep_data),
    keep_backend_fit = as.logical(control$keep_backend_fit),
    keep_log_lik = as.logical(control$keep_log_lik),
    verbose = as.logical(control$verbose)
  )
  class(canonical) <- c("pvstackr_control", "list")
  if (!identical(original_control, canonical)) {
    pv_abort("`control` fields must be bare canonical scalars without nested attributes.")
  }
  canonical
}
