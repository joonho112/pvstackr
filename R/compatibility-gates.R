pv_compatibility_family_link <- function(family = NULL) {
  if (is.null(family)) {
    family_id <- "gaussian"
    link_id <- "identity"
  } else if (is.character(family)) {
    family_id <- unname(pv_binding_validate_metadata_string(
      family,
      "family",
      "family_link"
    ))
    link_id <- "identity"
  } else if (is.list(family)) {
    input_names <- names(family)
    if (is.null(input_names) || anyNA(input_names) ||
        any(!nzchar(input_names)) || anyDuplicated(input_names) ||
        !all(c("family", "link") %in% input_names)) {
      pv_binding_abort(
        "PV_BIND_E081",
        "Family metadata must be a uniquely named list containing family and link.",
        "family_link"
      )
    }
    family_id <- unname(pv_binding_validate_metadata_string(
      family[["family"]],
      "family",
      "family_link"
    ))
    link_id <- unname(pv_binding_validate_metadata_string(
      family[["link"]],
      "link",
      "family_link"
    ))
  } else {
    pv_binding_abort(
      "PV_BIND_E081",
      "Family input must be NULL, one family name, or family/link metadata.",
      "family_link"
    )
  }

  if (!identical(family_id, "gaussian")) {
    pv_binding_abort(
      "PV_BIND_E060",
      "BRR-Fay reportability requires the Gaussian family.",
      "family_link"
    )
  }
  if (!identical(link_id, "identity")) {
    pv_binding_abort(
      "PV_BIND_E061",
      "BRR-Fay reportability requires the identity link.",
      "family_link"
    )
  }

  pv_binding_family_link_projection(
    family_id = "gaussian",
    link_id = "identity",
    response_support_id = "real",
    dispersion_role = "estimated"
  )
}

pv_compatibility_stack_direct_family <- function(family = NULL) {
  # Validate the caller-facing representation first, then discard every
  # executable closure or incidental field supplied by the caller. The backend
  # receives one package-owned Gaussian/identity implementation.
  pv_compatibility_family_link(family)
  pv_compatibility_canonical_gaussian()
}

pv_compatibility_canonical_gaussian <- function() {
  # Keep construction in a zero-argument frame. stats::gaussian() retains a
  # local closure environment, so constructing it in the caller-validation
  # frame would indirectly serialize the caller's family and private state.
  out <- stats::gaussian()
  out$linkfun(0)
  out$linkinv(0)
  out$variance(0)
  out
}

pv_compatibility_brr_fay_estimand <- function(
  fe_names,
  interval_role,
  coverage_claim_allowed,
  estimand_id = "brr_fay_fixed_effects",
  target_source = "external_brr_fay_rubin",
  target_engine_id = "lm_wls_brr_fay_v1",
  parameter_scope = "fixed_effect"
) {
  estimand_id <- unname(pv_binding_validate_metadata_string(
    estimand_id,
    "estimand_id",
    "estimand"
  ))
  target_source <- unname(pv_binding_validate_metadata_string(
    target_source,
    "target_source",
    "estimand"
  ))
  target_engine_id <- unname(pv_binding_validate_metadata_string(
    target_engine_id,
    "target_engine_id",
    "estimand"
  ))
  parameter_scope <- unname(pv_binding_validate_metadata_string(
    parameter_scope,
    "parameter_scope",
    "estimand"
  ))
  interval_role <- unname(pv_binding_validate_metadata_string(
    interval_role,
    "interval_role",
    "estimand"
  ))
  if (!is.logical(coverage_claim_allowed) ||
      length(coverage_claim_allowed) != 1L ||
      is.na(coverage_claim_allowed)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "coverage_claim_allowed must be one non-missing logical value.",
      "estimand"
    )
  }
  fe_names <- unname(pv_binding_validate_string_vector(fe_names, "fe_names"))

  expected_identity <- c(
    estimand_id = "brr_fay_fixed_effects",
    target_source = "external_brr_fay_rubin",
    target_engine_id = "lm_wls_brr_fay_v1",
    parameter_scope = "fixed_effect"
  )
  observed_identity <- c(
    estimand_id = estimand_id,
    target_source = target_source,
    target_engine_id = target_engine_id,
    parameter_scope = parameter_scope
  )
  valid_interval <-
    (identical(interval_role, "descriptive_classic_rubin") &&
      identical(coverage_claim_allowed, FALSE)) ||
    (identical(interval_role, "coverage_barnard_rubin") &&
      identical(coverage_claim_allowed, TRUE))
  if (!identical(observed_identity, expected_identity) || !valid_interval) {
    pv_binding_abort(
      "PV_BIND_E070",
      "BRR-Fay estimand identity, source, engine, scope, interval, or coverage is incompatible.",
      "estimand"
    )
  }

  estimand_contrast <- pv_binding_estimand_contrast_projection(fe_names)
  pv_binding_estimand_projection(
    fe_names = fe_names,
    estimand_contrast = estimand_contrast,
    interval_role = interval_role,
    estimand_id = estimand_id,
    target_source = target_source,
    target_engine_id = target_engine_id,
    parameter_scope = parameter_scope,
    coverage_claim_allowed = coverage_claim_allowed
  )
}
