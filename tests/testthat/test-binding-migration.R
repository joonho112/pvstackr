binding_legacy_from_current <- function(current) {
  legacy_fields <- setdiff(
    names(current),
    c("conf_level", "binding_manifest", "target_content")
  )
  legacy <- current[legacy_fields]
  if (identical(legacy$df_method, "classic")) {
    legacy$df_complete <- stats::setNames(
      rep(NA_real_, length(legacy$fe_names)),
      legacy$fe_names
    )
  }
  legacy$design_hash <- current$binding_manifest$legacy_hashes$design_hash
  legacy$target_hash <- current$binding_manifest$legacy_hashes$target_hash
  legacy$schema_version <- "0.1.0"
  legacy$provenance <- list(
    function_name = "pv_brr_target",
    package = "pvstackr",
    schema_version = "0.1.0"
  )
  class(legacy) <- c("pvstackr_brr_target", "list")
  pvstackr:::validate_pvstackr_brr_target(legacy)
  legacy
}

binding_migration_fixture <- function() {
  data <- data.frame(
    id = 1:8,
    x = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
    PV1 = c(1.2, 1.8, 2.4, 2.9, 3.6, 4.1, 4.8, 5.3),
    PV2 = c(1.0, 2.0, 2.2, 3.1, 3.4, 4.3, 4.6, 5.5),
    W = c(1, 1.2, 0.9, 1.1, 1, 0.8, 1.3, 1),
    RW1 = c(1.5, 0.8, 1.2, 0.7, 1.4, 0.9, 1.1, 0.6),
    RW2 = c(0.7, 1.4, 0.8, 1.3, 0.6, 1.2, 0.9, 1.5),
    stringsAsFactors = FALSE
  )
  formula <- OUTCOME ~ x
  current <- pv_brr_target(
    data = data,
    formula = formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = c("RW1", "RW2"),
    fay_k = 0.5,
    id_cols = "id",
    conf_level = 0.9,
    verbose = FALSE
  )

  legacy <- binding_legacy_from_current(current)
  list(data = data, formula = formula, current = current, legacy = legacy)
}

binding_inspection_target <- function(target) {
  manifest <- target$binding_manifest
  manifest$migration <- list(
    migration_from_schema = "0.1.0",
    migration_function = "pv_binding_revalidate_brr_target_v1",
    binding_revalidated = FALSE,
    inspection_only = TRUE,
    warnings = c(
      "PV_BIND_MIGRATION_LEGACY_HASH_ONLY",
      "PV_BIND_MIGRATION_SOURCE_METADATA_DROPPED"
    )
  )
  manifest$manifest_hash <- pvstackr:::pv_binding_hash_payload(
    pvstackr:::pv_binding_manifest_hash_payload(manifest),
    "manifest"
  )
  content <- pvstackr:::pv_binding_target_content_build(
    per_pv = target$per_pv,
    M = target$M,
    R = target$R,
    fay_k = target$fay_k,
    df_method = target$df_method,
    df_complete = if (identical(target$df_method, "classic")) NULL else target$df_complete,
    conf_level = target$conf_level,
    target_source = target$target_source,
    target_engine_id = target$target_content$target_engine_id,
    target_policy = target$policy,
    manifest_hash = manifest$manifest_hash,
    stored_derived = target[pvstackr:::pv_binding_target_derived_fields()]
  )
  target$binding_manifest <- manifest
  target$target_content <- content
  target$design_hash <- manifest$manifest_hash
  target$target_hash <- content$target_content_hash
  pvstackr:::validate_pvstackr_brr_target(target)
  target
}

test_that("legacy BRR target revalidation rebuilds one fresh bound target", {
  fixture <- binding_migration_fixture()
  source_before <- serialize(fixture$legacy, NULL, version = 3L)
  migration_data <- fixture$data
  migration_data$unused_private_note <- "respondent-secret-not-for-target"
  migrated <- pv_revalidate_brr_target(
    target = fixture$legacy,
    data = migration_data,
    formula = fixture$formula,
    conf_level = 0.9
  )

  expect_identical(serialize(fixture$legacy, NULL, version = 3L), source_before)
  expect_identical(migrated$schema_version, "0.2.0")
  expect_invisible(pvstackr:::validate_pvstackr_brr_target(migrated))
  expect_identical(
    migrated$binding_manifest$migration,
    list(
      migration_from_schema = "0.1.0",
      migration_function = "pv_binding_revalidate_brr_target_v1",
      binding_revalidated = TRUE,
      inspection_only = FALSE,
      warnings = c(
        "PV_BIND_MIGRATION_LEGACY_HASH_ONLY",
        "PV_BIND_MIGRATION_SOURCE_METADATA_DROPPED"
      )
    )
  )
  expect_identical(
    migrated$binding_manifest$legacy_hashes$design_hash,
    fixture$legacy$design_hash
  )
  expect_identical(
    migrated$binding_manifest$legacy_hashes$target_hash,
    fixture$legacy$target_hash
  )
  expect_identical(
    migrated$design_hash,
    migrated$binding_manifest$manifest_hash
  )
  expect_identical(
    migrated$target_hash,
    migrated$target_content$target_content_hash
  )
  expect_false(identical(migrated$design_hash, fixture$current$design_hash))
  migrated_again <- pv_revalidate_brr_target(
    target = fixture$legacy,
    data = migration_data,
    formula = fixture$formula,
    conf_level = 0.9
  )
  expect_identical(migrated_again$design_hash, migrated$design_hash)
  expect_identical(migrated_again$target_hash, migrated$target_hash)
  serialized_migrated <- serialize(migrated, NULL, version = 3L)
  secret_raw <- charToRaw("respondent-secret-not-for-target")
  candidate_starts <- seq_len(length(serialized_migrated) - length(secret_raw) + 1L)
  contains_secret <- any(vapply(
    candidate_starts,
    function(start) identical(
      serialized_migrated[start:(start + length(secret_raw) - 1L)],
      secret_raw
    ),
    logical(1)
  ))
  expect_false(contains_secret)
  expect_identical(
    pvstackr:::pv_stack_direct_preflight(
      migration_data,
      fixture$formula,
      migrated
    )$binding_proof$target_manifest_hash,
    migrated$binding_manifest$manifest_hash
  )
})

test_that("legacy targets remain inspection-only without complete raw inputs", {
  fixture <- binding_migration_fixture()
  expect_invisible(print(fixture$legacy))

  missing_data <- tryCatch(
    pv_revalidate_brr_target(
      fixture$legacy,
      formula = fixture$formula,
      conf_level = 0.9
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(missing_data, "pvstackr_binding_error")
  expect_identical(missing_data$code, "PV_BIND_E080")

  missing_level <- tryCatch(
    pv_revalidate_brr_target(
      fixture$legacy,
      data = fixture$data,
      formula = fixture$formula
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(missing_level, "pvstackr_binding_error")
  expect_identical(missing_level$code, "PV_BIND_E080")
})

test_that("self-consistent inspection-only schema-0.2 targets cannot become reportable", {
  fixture <- binding_migration_fixture()
  migrated <- pv_revalidate_brr_target(
    fixture$legacy,
    data = fixture$data,
    formula = fixture$formula,
    conf_level = 0.9
  )
  inspection <- binding_inspection_target(migrated)
  expect_invisible(pvstackr:::validate_pvstackr_brr_target(inspection))
  expect_output(print(inspection), "pvstackr BRR-Fay target")

  preflight_error <- tryCatch(
    pvstackr:::pv_stack_direct_preflight(
      fixture$data,
      fixture$formula,
      inspection
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(preflight_error, "pvstackr_binding_error")
  expect_identical(preflight_error$code, "PV_BIND_E080")

  callbacks <- new.env(parent = emptyenv())
  callbacks$fit <- 0L
  callbacks$draws <- 0L
  callbacks$diagnose <- 0L
  callbacks$loglik <- 0L
  fit_error <- tryCatch(
    pv_fit_direct(
      data = fixture$data,
      formula = fixture$formula,
      target = inspection,
      fit_function = function(...) { callbacks$fit <- callbacks$fit + 1L; list() },
      draws_function = function(fit) { callbacks$draws <- callbacks$draws + 1L; matrix(0) },
      diagnose_function = function(fit) { callbacks$diagnose <- callbacks$diagnose + 1L; list() },
      log_lik_function = function(fit) { callbacks$loglik <- callbacks$loglik + 1L; matrix(0) },
      extract_log_lik = TRUE
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(fit_error, "pvstackr_binding_error")
  expect_identical(fit_error$code, "PV_BIND_E080")
  expect_identical(
    c(callbacks$fit, callbacks$draws, callbacks$diagnose, callbacks$loglik),
    rep(0L, 4L)
  )

  forged_proof <- list(
    contract_id = inspection$binding_manifest$contract_id,
    manifest_schema_version = inspection$binding_manifest$manifest_schema_version,
    target_manifest_hash = inspection$binding_manifest$manifest_hash,
    current_manifest_hash = inspection$binding_manifest$manifest_hash,
    verification_policy = pvstackr:::pv_binding_proof_policy()
  )
  ccc_error <- tryCatch(
    pvstackr:::ccc_calibrate(
      draws = "invalid-draw-secret",
      target = inspection,
      binding_proof = forged_proof
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(ccc_error, "pvstackr_binding_error")
  expect_identical(ccc_error$code, "PV_BIND_E080")
  expect_false(grepl("invalid-draw-secret", conditionMessage(ccc_error), fixed = TRUE))
})

test_that("retained legacy confidence level cannot be overridden", {
  fixture <- binding_migration_fixture()
  retained <- fixture$legacy
  retained$conf_level <- 0.9
  expect_invisible(pvstackr:::validate_pvstackr_brr_target(retained))

  defaulted <- pv_revalidate_brr_target(
    retained,
    data = fixture$data,
    formula = fixture$formula
  )
  matched <- pv_revalidate_brr_target(
    retained,
    data = fixture$data,
    formula = fixture$formula,
    conf_level = 0.9
  )
  expect_identical(defaulted$target_hash, matched$target_hash)
  source_before <- serialize(retained, NULL, version = 3L)
  mismatch <- tryCatch(
    pv_revalidate_brr_target(
      retained,
      data = fixture$data,
      formula = fixture$formula,
      conf_level = 0.95
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(mismatch, "pvstackr_binding_error")
  expect_identical(mismatch$code, "PV_BIND_E080")
  expect_identical(serialize(retained, NULL, version = 3L), source_before)

  helper_mismatch <- tryCatch(
    pvstackr:::pv_binding_target_content_from_brr_target(
      target = retained,
      manifest = fixture$current$binding_manifest,
      conf_level = 0.95,
      data = fixture$data,
      formula = fixture$formula
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(helper_mismatch, "pvstackr_binding_error")
  expect_identical(helper_mismatch$code, "PV_BIND_E080")
})

test_that("legacy revalidation rejects stale raw inputs and tampered primitives", {
  fixture <- binding_migration_fixture()
  changed <- fixture$data
  changed$PV1[[1L]] <- changed$PV1[[1L]] + 1
  stale_error <- tryCatch(
    pv_revalidate_brr_target(
      fixture$legacy,
      data = changed,
      formula = fixture$formula,
      conf_level = 0.9
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(stale_error, "pvstackr_binding_error")
  expect_identical(stale_error$code, "PV_BIND_E080")

  tampered <- fixture$legacy
  tampered$beta[[1L]] <- tampered$beta[[1L]] + 1
  tamper_error <- tryCatch(
    pv_revalidate_brr_target(
      tampered,
      data = fixture$data,
      formula = fixture$formula,
      conf_level = 0.9
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(tamper_error, "pvstackr_binding_error")
  expect_identical(tamper_error$code, "PV_BIND_E080")

  nonfinite <- fixture$data
  nonfinite$x[[1L]] <- Inf
  nonfinite_error <- tryCatch(
    pv_revalidate_brr_target(
      fixture$legacy,
      data = nonfinite,
      formula = fixture$formula,
      conf_level = 0.9
    ),
    pvstackr_binding_error = identity
  )
  expect_s3_class(nonfinite_error, "pvstackr_binding_error")
  expect_identical(nonfinite_error$code, "PV_BIND_E081")
})

test_that("revalidation is idempotent for current bound targets", {
  fixture <- binding_migration_fixture()
  before <- serialize(fixture$current, NULL, version = 3L)
  expect_identical(
    pv_revalidate_brr_target(fixture$current),
    fixture$current
  )
  expect_identical(serialize(fixture$current, NULL, version = 3L), before)
})

test_that("legacy migration supports Barnard-Rubin and explicitly allowed M equals one", {
  fixture <- binding_migration_fixture()
  barnard_current <- pv_brr_target(
    data = fixture$data,
    formula = fixture$formula,
    pv_cols = c("PV1", "PV2"),
    weight_col = "W",
    rep_weight_cols = c("RW1", "RW2"),
    fay_k = 0.5,
    id_cols = "id",
    conf_level = 0.9,
    df_method = "barnard_rubin",
    df_complete = 120,
    verbose = FALSE
  )
  barnard <- pv_revalidate_brr_target(
    binding_legacy_from_current(barnard_current),
    data = fixture$data,
    formula = fixture$formula,
    conf_level = 0.9
  )
  expect_identical(barnard$df_method, "barnard_rubin")
  expect_identical(barnard$df_complete, barnard_current$df_complete)
  expect_invisible(pvstackr:::validate_pvstackr_brr_target(barnard))

  m1_current <- pv_brr_target(
    data = fixture$data,
    formula = fixture$formula,
    pv_cols = "PV1",
    weight_col = "W",
    rep_weight_cols = c("RW1", "RW2"),
    fay_k = 0.5,
    id_cols = "id",
    conf_level = 0.9,
    allow_m1 = TRUE,
    verbose = FALSE
  )
  m1 <- pv_revalidate_brr_target(
    binding_legacy_from_current(m1_current),
    data = fixture$data,
    formula = fixture$formula,
    conf_level = 0.9
  )
  expect_identical(m1$M, 1L)
  expect_invisible(pvstackr:::validate_pvstackr_brr_target(m1))
})

test_that("bundled golden target supports a serialized schema-0.1 migration projection", {
  path <- system.file(
    "extdata", "examples", "pisa_tiny_stack_direct.rds",
    package = "pvstackr"
  )
  expect_true(nzchar(path))
  fixture <- readRDS(path)
  legacy_target <- if (identical(fixture$target$schema_version, "0.1.0")) {
    fixture$target
  } else {
    binding_legacy_from_current(fixture$target)
  }
  expect_identical(legacy_target$schema_version, "0.1.0")
  source_before <- serialize(legacy_target, NULL, version = 3L)
  migrated <- pv_revalidate_brr_target(
    target = legacy_target,
    data = fixture$design$data,
    formula = fixture$design$formula,
    conf_level = fixture$fit$control$conf_level
  )
  expect_identical(migrated$schema_version, "0.2.0")
  expect_invisible(pvstackr:::validate_pvstackr_brr_target(migrated))
  expect_identical(
    migrated$binding_manifest$legacy_hashes$target_hash,
    legacy_target$target_hash
  )
  expect_identical(serialize(legacy_target, NULL, version = 3L), source_before)
})

test_that("legacy algorithm identifier matches the frozen contract", {
  expect_identical(
    pvstackr:::pv_binding_legacy_hash_algorithm_id(),
    "adler32_r_serialize_v2"
  )
})
