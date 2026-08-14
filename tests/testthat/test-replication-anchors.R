replication_anchor_pkg_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
}

replication_anchor_dev_files <- function() {
  root <- replication_anchor_pkg_root()
  list(
    root = root,
    script = file.path(root, "dev", "verify-replication-anchors.R"),
    manifest = file.path(root, "dev", "replication-anchors.csv")
  )
}

replication_anchor_env <- function() {
  files <- replication_anchor_dev_files()
  skip_if_not(
    file.exists(files$script) && file.exists(files$manifest),
    "replication anchor audit runs only from the package source tree"
  )
  env <- new.env(parent = baseenv())
  sys.source(files$script, envir = env)
  list(files = files, env = env)
}

test_that("replication manifest pins only aggregate public tables", {
  audit <- replication_anchor_env()
  manifest <- audit$env$pv_replication_read_manifest(audit$files$manifest)

  expect_silent(audit$env$pv_replication_validate_manifest(manifest))
  expect_equal(nrow(manifest), 12L)
  expect_setequal(manifest$repo_id, c("companion", "lsae"))
  expect_setequal(
    unique(manifest$artifact_path),
    c(
      "output/tables/table3_simulation_headline.csv",
      "output/tables/tab3_rows.tex",
      "data/precomputed/tables/master_results.csv",
      "data/precomputed/tables/osm_survey_crosscheck.csv",
      "data/precomputed/tables/compare_agreement.csv",
      "data/precomputed/tables/psis_pareto_k.csv"
    )
  )
  expect_true(all(manifest$data_class == "aggregate_derived_table"))
  expect_true(all(!manifest$contains_unit_records))
})

test_that("replication manifest rejects privacy and path-boundary mutations", {
  audit <- replication_anchor_env()
  manifest <- audit$env$pv_replication_read_manifest(audit$files$manifest)

  bad <- manifest
  bad$contains_unit_records[[1L]] <- TRUE
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "exclude unit records"
  )

  bad <- manifest
  bad$artifact_path[[1L]] <- "data/pisa/student_rows.csv"
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "aggregate-table allowlist"
  )

  bad <- manifest
  bad$artifact_path[[1L]] <- "output/tables/../student_rows.csv"
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "aggregate-table allowlist"
  )

  bad <- manifest
  bad$artifact_path[[1L]] <- "output/tables/fit.rds"
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "aggregate-table allowlist"
  )

  bad <- manifest
  bad$artifact_path[[1L]] <- "output/tables/small_rows.csv"
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "aggregate-table allowlist"
  )

  bad <- manifest
  bad$acceptance_role[bad$repo_id == "companion"] <- "current_v02_frozen_evidence"
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "outside the allowlist"
  )

  bad <- manifest
  bad$tolerance[[1L]] <- Inf
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "finite and nonnegative"
  )

  bad <- manifest
  last_lsae <- tail(which(bad$repo_id == "lsae"), 1L)
  bad$repo_commit[[last_lsae]] <- paste(rep("a", 40L), collapse = "")
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "Conflicting repository metadata"
  )

  bad <- manifest
  bad$article_id[bad$repo_id == "lsae"] <- "companion-methods-v5"
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "outside the allowlist"
  )

  bad <- manifest
  bad$repo_id[] <- "lsae"
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "both required public repositories"
  )

  bad <- manifest
  bad$source_claim_id[bad$repo_id == "companion"] <- ""
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "requires a verified source claim"
  )

  bad <- manifest[manifest$anchor_id != "lsae_direct_agreement_z", , drop = FALSE]
  expect_error(
    audit$env$pv_replication_validate_manifest(bad),
    "anchor set is incomplete"
  )
})

test_that("replication table parser selects exact CSV and TeX anchors", {
  audit <- replication_anchor_env()
  tmp <- tempfile("replication-anchor-parser-")
  dir.create(tmp)

  csv <- file.path(tmp, "aggregate.csv")
  utils::write.csv(
    data.frame(group = c("A", "B"), method = c("cal", "cal"), value = c(0.95, 0.96)),
    csv,
    row.names = FALSE,
    quote = TRUE
  )
  csv_payload <- audit$env$pv_replication_read_raw(csv)
  csv_data <- audit$env$pv_replication_read_artifact(
    csv_payload, "csv", "group|method|value", "group|method", 2L, 3L
  )
  expect_equal(
    audit$env$pv_replication_select_anchor(csv_data, "group=B;method=cal", "value"),
    0.96
  )

  tex <- file.path(tmp, "aggregate.tex")
  writeLines(
    c(
      "% aggregate rows only",
      "Within & 25.35 & 2.31 \\\\",
      "Between & 80.73 & 7.36 \\\\"
    ),
    tex
  )
  tex_payload <- audit$env$pv_replication_read_raw(tex)
  tex_data <- audit$env$pv_replication_read_artifact(
    tex_payload, "tex_ampersand", "term|estimate|std_error", "term", 2L, 3L
  )
  expect_equal(
    audit$env$pv_replication_select_anchor(tex_data, "term=Between", "std_error"),
    7.36
  )
})

test_that("replication public-table scan fails closed on unit identifiers", {
  audit <- replication_anchor_env()
  tmp <- tempfile(fileext = ".csv")
  writeLines("CNTSTUID,estimate\n123,1.5", tmp)
  expect_error(
    audit$env$pv_replication_scan_public_table(
      audit$env$pv_replication_read_raw(tmp)
    ),
    "unit-record identifiers"
  )
})

test_that("replication parser rejects disguised large tables and duplicate keys", {
  audit <- replication_anchor_env()
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(key = seq_len(101L), estimate = rep(1, 101L)),
    tmp,
    row.names = FALSE
  )
  payload <- audit$env$pv_replication_read_raw(tmp)
  expect_error(
    audit$env$pv_replication_read_artifact(
      payload, "csv", "key|estimate", "key", 9L, 2L
    ),
    "dimensions have drifted"
  )

  writeLines("key,estimate\nA,1\nA,2", tmp)
  payload <- audit$env$pv_replication_read_raw(tmp)
  expect_error(
    audit$env$pv_replication_read_artifact(
      payload, "csv", "key|estimate", "key", 2L, 2L
    ),
    "duplicate keys"
  )
})

test_that("replication selector rejects zero, multiple, and nonnumeric matches", {
  audit <- replication_anchor_env()
  data <- data.frame(key = c("A", "A"), value = c(1, 2))
  attr(data, "pv_replication_format") <- "csv"
  expect_error(
    audit$env$pv_replication_select_anchor(data, "key=B", "value"),
    "matched 0 rows"
  )
  expect_error(
    audit$env$pv_replication_select_anchor(data, "key=A", "value"),
    "matched 2 rows"
  )
  data <- data.frame(key = "A", value = "1.0", stringsAsFactors = FALSE)
  attr(data, "pv_replication_format") <- "csv"
  expect_error(
    audit$env$pv_replication_select_anchor(data, "key=A", "value"),
    "natively numeric"
  )
})

test_that("replication path resolver rejects a symlink escape", {
  audit <- replication_anchor_env()
  root <- tempfile("replication-root-")
  outside <- tempfile("replication-outside-", fileext = ".csv")
  dir.create(file.path(root, "output", "tables"), recursive = TRUE)
  writeLines("key,value\nA,1", outside)
  link <- file.path(root, "output", "tables", "aggregate.csv")
  linked <- file.symlink(outside, link)
  skip_if_not(linked, "symbolic links are unavailable on this platform")
  expect_error(
    audit$env$pv_replication_resolve_file(
      root, "output/tables/aggregate.csv"
    ),
    "escapes its repository root"
  )
})

test_that("replication package root follows the real script through a symlink", {
  audit <- replication_anchor_env()
  fake <- tempfile("replication-script-link-")
  dir.create(file.path(fake, "dev"), recursive = TRUE)
  link <- file.path(fake, "dev", "verify-replication-anchors.R")
  linked <- file.symlink(audit$files$script, link)
  skip_if_not(linked, "symbolic links are unavailable on this platform")
  expect_identical(
    audit$env$pv_replication_root_from_script(link),
    audit$files$root
  )
})

test_that("replication Git child metadata cannot escape through symlinks", {
  audit <- replication_anchor_env()
  make_root <- function() {
    root <- tempfile("replication-git-root-")
    dir.create(file.path(root, ".git", "refs", "heads"), recursive = TRUE)
    root
  }

  root <- make_root()
  outside_head <- tempfile("outside-head-")
  writeLines(paste(rep("a", 40L), collapse = ""), outside_head)
  expect_true(file.symlink(outside_head, file.path(root, ".git", "HEAD")))
  expect_error(audit$env$pv_replication_git_head(root), "escapes .git")

  root <- make_root()
  writeLines("ref: refs/heads/main", file.path(root, ".git", "HEAD"))
  outside_ref <- tempfile("outside-ref-")
  writeLines(paste(rep("b", 40L), collapse = ""), outside_ref)
  expect_true(file.symlink(outside_ref, file.path(root, ".git", "refs", "heads", "main")))
  expect_error(audit$env$pv_replication_git_head(root), "escapes .git")

  root <- make_root()
  outside_config <- tempfile("outside-config-")
  writeLines(
    c("[remote \"origin\"]", "url = https://github.com/example/repo.git"),
    outside_config
  )
  expect_true(file.symlink(outside_config, file.path(root, ".git", "config")))
  expect_error(audit$env$pv_replication_git_origin(root), "escapes .git")

  root <- make_root()
  writeLines("ref: refs/heads/main", file.path(root, ".git", "HEAD"))
  outside_packed <- tempfile("outside-packed-")
  writeLines(
    paste(paste(rep("c", 40L), collapse = ""), "refs/heads/main"),
    outside_packed
  )
  expect_true(file.symlink(outside_packed, file.path(root, ".git", "packed-refs")))
  expect_error(audit$env$pv_replication_git_head(root), "escapes .git")
})

test_that("replication pinned reader brackets one payload with SHA-256", {
  audit <- replication_anchor_env()
  file <- tempfile(fileext = ".csv")
  writeLines("key,value\nA,1", file)
  hash <- unname(tools::sha256sum(file)[[1L]])
  payload <- audit$env$pv_replication_read_pinned(file, hash, "fixture")
  expect_type(payload, "raw")
  expect_identical(audit$env$pv_replication_payload_sha256(payload), hash)
  expect_error(
    audit$env$pv_replication_read_pinned(
      file, paste(rep("0", 64L), collapse = ""), "fixture"
    ),
    "mismatch before read"
  )
})

test_that("replication opt-in accepts only explicit roots and package-owned reads", {
  audit <- replication_anchor_env()
  expect_false(audit$env$pv_replication_parse_args(character())$run)
  parsed <- audit$env$pv_replication_parse_args(c(
    "--run", "--companion-root=/tmp/companion", "--lsae-root=/tmp/lsae"
  ))
  expect_true(parsed$run)
  expect_setequal(names(parsed$roots), c("companion", "lsae"))
  expect_error(
    audit$env$pv_replication_parse_args(c("--run", "--lsae-root=/tmp/lsae")),
    "requires one"
  )
  expect_silent(audit$env$pv_replication_require_vanilla("--vanilla"))
  expect_error(
    audit$env$pv_replication_require_vanilla("R"),
    "Rscript --vanilla"
  )

  risky <- c(
    "R_PROFILE_USER", "R_DEFAULT_PACKAGES", "R_SCRIPT_DEFAULT_PACKAGES", "R_LIBS"
  )
  old_risky <- Sys.getenv(risky, unset = NA_character_)
  on.exit({
    for (variable in risky) {
      if (is.na(old_risky[[variable]])) {
        Sys.unsetenv(variable)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old_risky[[variable]]), variable))
      }
    }
  }, add = TRUE)
  for (variable in risky) {
    Sys.unsetenv(risky)
    do.call(Sys.setenv, stats::setNames(list("/tmp/untrusted-startup"), variable))
    expect_error(
      audit$env$pv_replication_require_vanilla("--vanilla"),
      "startup injection variables",
      info = variable
    )
  }

  script <- paste(readLines(audit$files$script, warn = FALSE), collapse = "\n")
  expect_false(grepl(
    "\\b(source|download[.]file|url|writeLines|saveRDS|readRDS)\\s*\\(",
    script,
    perl = TRUE
  ))
  expect_false(grepl("system2", script, fixed = TRUE))
  expect_false(grepl("digest::", script, fixed = TRUE))
  expect_match(script, "upstream_replication_code_executed=FALSE", fixed = TRUE)
  expect_match(script, "algorithm_parity=NOT_CLAIMED", fixed = TRUE)
})

test_that("replication verifier is opt-in by default", {
  audit <- replication_anchor_env()
  out <- system2(
    "Rscript",
    shQuote(audit$files$script),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(out, "status"))
  expect_true(any(grepl("verification SKIP", out, fixed = TRUE)))
})
