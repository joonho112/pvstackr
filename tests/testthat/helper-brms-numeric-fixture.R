brms_numeric_fixture_manifest <- function() {
  path <- testthat::test_path("fixtures", "brms_numeric_fixture.dcf")
  read.dcf(path)[1L, ]
}

brms_numeric_require_real_posterior <- function() {
  desc <- suppressWarnings(utils::packageDescription("posterior"))
  if (length(desc) == 1L && is.na(desc)) {
    stop("Explicit numeric fixture acceptance requires `posterior`.")
  }
  if (identical(desc$Title, "Sentinel Optional Package") ||
      identical(desc$Version, "0.0.0")) {
    stop(
      paste(
        "Explicit numeric fixture acceptance requires the real posterior",
        "package, not a sentinel."
      )
    )
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Explicit numeric fixture acceptance could not load `posterior`.")
  }
  invisible(TRUE)
}

brms_numeric_target <- function() {
  bundle <- pisa_tiny_parity_load()
  design <- pv_design(
    data = bundle$data,
    formula = OUTCOME ~ x + female,
    pv_suffix = "READ",
    expected_M = 2L,
    expected_R = 4L,
    id_cols = "CNTSTUID"
  )
  target <- pv_brr_target(
    data = bundle$data,
    formula = OUTCOME ~ x + female,
    pv_cols = design$pv_cols,
    weight_col = design$weight_col,
    rep_weight_cols = design$rep_weight_cols,
    fay_k = design$fay_k,
    id_cols = design$id_cols
  )
  list(data = bundle$data, design = design, target = target)
}

brms_numeric_generate <- function(target, manifest) {
  seed <- as.integer(manifest[["Seed"]])
  chains <- as.integer(manifest[["Chains"]])
  post_warmup <- as.integer(manifest[["Post-Warmup-Per-Chain"]])
  n_draws <- chains * post_warmup
  fe_names <- target$fe_names
  parameter_names <- c(fe_names, "sigma")

  old_rng <- RNGkind()
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    do.call(RNGkind, as.list(old_rng))
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)

  z <- matrix(stats::rnorm(n_draws * length(fe_names)), ncol = length(fe_names))
  z <- sweep(z, 2L, colMeans(z), FUN = "-")
  whitened <- z %*% solve(chol(stats::cov(z)))
  raw_scales <- c(0.6, 1.4, 0.8)
  raw_covariance <- diag(diag(target$T_MI) * raw_scales)
  dimnames(raw_covariance) <- list(fe_names, fe_names)
  raw_fe <- whitened %*% chol(raw_covariance)
  raw_fe <- sweep(raw_fe, 2L, target$beta[fe_names], FUN = "+")
  colnames(raw_fe) <- fe_names

  sigma <- exp(0.25 * stats::rnorm(n_draws))
  draws_matrix <- cbind(raw_fe, sigma = sigma)
  draws_array <- array(
    as.double(draws_matrix),
    dim = c(post_warmup, chains, length(parameter_names)),
    dimnames = list(
      iteration = as.character(seq_len(post_warmup)),
      chain = as.character(seq_len(chains)),
      variable = parameter_names
    )
  )

  list(
    draws_matrix = draws_matrix,
    draws_array = draws_array,
    raw_covariance = raw_covariance,
    raw_scales = raw_scales,
    chains = chains,
    post_warmup = post_warmup,
    seed = seed,
    rng = c("Mersenne-Twister", "Inversion", "Rejection"),
    generator = unname(manifest[["Generator"]]),
    schema = unname(manifest[["Fixture-Schema"]]),
    claim_scope = list(
      acceptance_role = unname(manifest[["Acceptance-Role"]]),
      fixture_origin = unname(manifest[["Fixture-Origin"]]),
      empirical_backend_accuracy = unname(manifest[["Empirical-Backend-Accuracy"]]),
      bundled_sampling_tested_here = unname(manifest[["Bundled-Sampling-Tested-Here"]]),
      live_sampler_quality_validated = unname(manifest[["Live-Sampler-Quality-Validated"]]),
      model_recovery_validated = unname(manifest[["Model-Recovery-Validated"]]),
      coverage_validated = unname(manifest[["Coverage-Validated"]]),
      real_data_evidence = unname(manifest[["Real-Data-Evidence"]]),
      sampler_reference = unname(manifest[["Sampler-Reference"]])
    )
  )
}

brms_numeric_fixture_hash <- function(fixture, target_hash) {
  paste0(
    "sha256:",
    digest::digest(list(
      schema_version = fixture$schema,
      generator_id = fixture$generator,
      seed = fixture$seed,
      rng_kind = fixture$rng,
      dimensions = c(
        iterations_per_chain = fixture$post_warmup,
        chains = fixture$chains
      ),
      parameter_names = colnames(fixture$draws_matrix),
      target_hash = target_hash,
      claim_scope = fixture$claim_scope,
      raw_covariance_scale = unname(fixture$raw_scales),
      raw_covariance = unname(round(fixture$raw_covariance, 12L)),
      values = as.double(round(fixture$draws_matrix, 12L))
    ), algo = "sha256", serialize = TRUE, serializeVersion = 2)
  )
}

brms_numeric_posterior_reference <- function(draws_array) {
  summary <- posterior::summarise_draws(
    posterior::as_draws_array(draws_array),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  list(
    summary = summary,
    rhat_max = max(summary$rhat),
    ess_bulk_min = min(summary$ess_bulk),
    ess_tail_min = min(summary$ess_tail),
    divergences = 0L,
    chains = as.integer(dim(draws_array)[[2L]]),
    post_warmup_draws_per_chain = as.integer(dim(draws_array)[[1L]])
  )
}

brms_numeric_nuts <- function(chains, post_warmup, divergences = 0L) {
  total <- chains * post_warmup
  values <- integer(total)
  if (divergences > 0L) {
    values[seq_len(divergences)] <- 1L
  }
  data.frame(
    Parameter = rep("divergent__", total),
    Value = values
  )
}

brms_numeric_fit_function <- function(draws_matrix) {
  force(draws_matrix)
  function(formula, data, family, prior, chains, iter, warmup, cores, seed,
           backend, file, file_refit, ...) {
    list(draws = draws_matrix)
  }
}

brms_numeric_diagnose_function <- function(reference) {
  force(reference)
  function(fit) {
    list(sampler = list(
      rhat_max = reference$rhat_max,
      ess_bulk_min = reference$ess_bulk_min,
      ess_tail_min = reference$ess_tail_min,
      divergences = reference$divergences,
      chains = reference$chains,
      post_warmup_draws_per_chain = reference$post_warmup_draws_per_chain
    ))
  }
}

brms_numeric_frobenius_relative <- function(observed, expected) {
  sqrt(sum((observed - expected)^2)) / sqrt(sum(expected^2))
}
