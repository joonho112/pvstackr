pv_binding_revalidate_brr_target_v1 <- function(
  target,
  data = NULL,
  formula = NULL,
  conf_level = NULL
) {
  source_before <- serialize(target, NULL, version = 3L)
  tryCatch(
    validate_pvstackr_brr_target(target),
    error = function(error) {
      if (inherits(error, "pvstackr_binding_error")) {
        stop(error)
      }
      pv_binding_legacy_revalidation_abort(
        "The source BRR-Fay target is not a valid supported schema object."
      )
    }
  )

  if (identical(target$schema_version, "0.2.0")) {
    return(target)
  }
  if (!identical(target$schema_version, "0.1.0") ||
      is.null(data) || is.null(formula) ||
      !is.data.frame(data) || !inherits(formula, "formula")) {
    pv_binding_legacy_revalidation_abort(
      paste0(
        "Schema-0.1 BRR-Fay targets are inspection-only until complete ",
        "original data and formula inputs are supplied for revalidation."
      )
    )
  }
  retained_conf_level <- target$conf_level
  if (is.null(conf_level)) {
    conf_level <- retained_conf_level
  }
  if (is.null(conf_level)) {
    pv_binding_legacy_revalidation_abort(
      "Legacy target revalidation requires an explicit conf_level."
    )
  }
  conf_level <- tryCatch(
    pv_assert_probability(conf_level, "conf_level"),
    error = function(error) {
      pv_binding_legacy_revalidation_abort(
        "Legacy target revalidation requires a valid confidence level."
      )
    }
  )
  if (!is.null(retained_conf_level)) {
    retained_conf_level <- tryCatch(
      pv_assert_probability(retained_conf_level, "target$conf_level"),
      error = function(error) NULL
    )
    if (is.null(retained_conf_level) ||
        !identical(as.double(conf_level), as.double(retained_conf_level))) {
      pv_binding_legacy_revalidation_abort(
        "The supplied confidence level does not match the value retained on the legacy target."
      )
    }
    conf_level <- retained_conf_level
  }

  df_complete <- if (identical(target$df_method, "classic")) {
    NULL
  } else {
    target$df_complete
  }
  rebuilt <- tryCatch(
    pv_brr_target(
      data = data,
      formula = formula,
      pv_cols = target$pv_cols,
      weight_col = target$weight_col,
      rep_weight_cols = target$rep_weight_cols,
      fay_k = target$fay_k,
      id_cols = target$id_cols,
      conf_level = conf_level,
      allow_m1 = identical(target$M, 1L),
      df_method = target$df_method,
      df_complete = df_complete,
      engine = "lm",
      verbose = FALSE
    ),
    pvstackr_binding_error = function(error) stop(error),
    error = function(error) {
      pv_binding_abort(
        "PV_BIND_E081",
        "Legacy target raw-input reconstruction did not complete successfully.",
        "manifest"
      )
    }
  )

  proven_content <- pv_binding_target_content_from_brr_target(
    target = target,
    manifest = rebuilt$binding_manifest,
    conf_level = conf_level,
    data = data,
    formula = formula
  )
  if (!identical(proven_content, rebuilt$target_content)) {
    pv_binding_legacy_revalidation_abort(
      "Legacy target reconstruction did not produce one canonical target content object."
    )
  }

  migration <- list(
    migration_from_schema = "0.1.0",
    migration_function = "pv_binding_revalidate_brr_target_v1",
    binding_revalidated = TRUE,
    inspection_only = FALSE,
    warnings = c(
      "PV_BIND_MIGRATION_LEGACY_HASH_ONLY",
      "PV_BIND_MIGRATION_SOURCE_METADATA_DROPPED"
    )
  )
  manifest <- c(rebuilt$binding_manifest, list(migration = migration))
  manifest$manifest_hash <- pv_binding_hash_payload(
    pv_binding_manifest_hash_payload(manifest),
    "manifest"
  )
  pv_binding_manifest_validate(manifest)

  migrated_content <- pv_binding_target_content_build(
    per_pv = rebuilt$per_pv,
    M = rebuilt$M,
    R = rebuilt$R,
    fay_k = rebuilt$fay_k,
    df_method = rebuilt$df_method,
    df_complete = if (identical(rebuilt$df_method, "classic")) {
      NULL
    } else {
      rebuilt$df_complete
    },
    conf_level = rebuilt$conf_level,
    target_source = rebuilt$target_source,
    target_engine_id = "lm_wls_brr_fay_v1",
    target_policy = rebuilt$policy,
    manifest_hash = manifest$manifest_hash,
    stored_derived = rebuilt[pv_binding_target_derived_fields()]
  )
  pv_binding_target_manifest_validate(migrated_content, manifest)

  migrated <- rebuilt
  migrated$design_hash <- manifest$manifest_hash
  migrated$target_hash <- migrated_content$target_content_hash
  migrated$binding_manifest <- manifest
  migrated$target_content <- migrated_content
  validate_pvstackr_brr_target(migrated)
  if (!identical(serialize(target, NULL, version = 3L), source_before)) {
    pv_binding_legacy_revalidation_abort(
      "Legacy revalidation mutated its source target."
    )
  }
  migrated
}

#' Check a BRR-Fay target or rebuild it in the current format
#'
#' `pv_revalidate_brr_target()` checks a target from [pv_brr_target()] and
#' returns it in the current format. A target in the current format
#' (`schema_version` `"0.2.0"`) is returned unchanged; the function stops
#' with an error if that target was changed after it was created. A target
#' in the format of pvstackr 0.1.x (`schema_version` `"0.1.0"`) is computed
#' again from `data` and `formula`, and the new target is returned if it
#' matches the old one. In pvstackr 0.2.x no target saved by pvstackr 0.1.x
#' matches, because pvstackr 0.2.0 computes the design checksum
#' `design_hash` differently (see the package NEWS); compute such a target
#' again with [pv_brr_target()]. The target that you pass is not modified.
#'
#' @details
#' The old target is computed again with [pv_brr_target()], using the
#' columns and settings stored on it. The new target must match the old one:
#' the formula, the settings (such as the numbers of plausible values and
#' replicate weights, the Fay coefficient and the degrees-of-freedom rule),
#' the fixed-effect names and `design_hash` exactly, and the estimates,
#' covariance matrices, standard errors and degrees of freedom up to a small
#' numerical tolerance (absolute 1e-12 plus relative 1e-10), because the last
#' digits of least-squares results can differ between computers. Otherwise,
#' and for a target in the old format without `data` and `formula`, the
#' function stops with an error. The new target records the upgrade in
#' `binding_manifest$migration`.
#'
#' @param target A BRR-Fay target (class `pvstackr_brr_target`).
#' @param data For a target in the old format, the data frame from which it
#'   was built; default `NULL`. Not used for a target in the current format.
#' @param formula For a target in the old format, the formula from which it
#'   was built; it must be the formula stored on the target. Default `NULL`;
#'   not used for a target in the current format.
#' @param conf_level For a target in the old format, the confidence level to
#'   store in the new target, a number between 0 and 1. It is required when
#'   the old target does not store one (targets saved by pvstackr 0.1.x do
#'   not); when it does, give `NULL` (default) or the same value. Not used
#'   for a target in the current format.
#'
#' @returns A target (class `pvstackr_brr_target`) in the current format: the
#'   input itself when it is already in that format, otherwise the new
#'   target.
#' @export
pv_revalidate_brr_target <- function(
  target,
  data = NULL,
  formula = NULL,
  conf_level = NULL
) {
  pv_binding_revalidate_brr_target_v1(
    target = target,
    data = data,
    formula = formula,
    conf_level = conf_level
  )
}
