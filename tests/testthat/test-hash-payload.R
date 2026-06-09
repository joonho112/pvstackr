test_that("pv_hash_payload is content-sensitive while preserving 8-hex shape", {
  payload1 <- list(beta = c(1.234, -2), x = seq_len(60))
  payload2 <- list(beta = c(1.235, -2), x = seq_len(60))
  tail1 <- list(x = as.numeric(seq_len(10000)))
  tail2 <- tail1
  tail2$x[[9000]] <- 99999

  hash1 <- pvstackr:::pv_hash_payload(payload1)

  expect_match(hash1, "^[0-9a-f]{8}$")
  expect_false(identical(hash1, pvstackr:::pv_hash_payload(payload2)))
  expect_false(identical(
    pvstackr:::pv_hash_payload(tail1),
    pvstackr:::pv_hash_payload(tail2)
  ))
})

test_that("pv_hash_payload is independent of display options", {
  payload <- list(
    beta = c(pi, sqrt(2), exp(1)),
    matrix = matrix(seq(0.1, 1.2, length.out = 12), nrow = 3),
    labels = c("alpha", "beta")
  )

  old_options <- options(digits = 7, width = 80)
  on.exit(options(old_options), add = TRUE)
  hash_default <- pvstackr:::pv_hash_payload(payload)

  options(digits = 3, width = 20)
  hash_compact <- pvstackr:::pv_hash_payload(payload)

  expect_identical(hash_compact, hash_default)
})

test_that("pv_hash_payload handles empty and NULL payloads deterministically", {
  null_hash <- pvstackr:::pv_hash_payload(NULL)
  empty_character_hash <- pvstackr:::pv_hash_payload(character())

  expect_match(null_hash, "^[0-9a-f]{8}$")
  expect_match(empty_character_hash, "^[0-9a-f]{8}$")
  expect_identical(null_hash, pvstackr:::pv_hash_payload(NULL))
  expect_false(identical(null_hash, pvstackr:::pv_hash_payload(list())))
})

test_that("pv_hash_payload is stable across ALTREP and materialized integer vectors", {
  altrep_payload <- data.frame(x = 1:10000)
  materialized_payload <- data.frame(x = as.integer(c(1:10000)))

  expect_identical(
    pvstackr:::pv_hash_payload(altrep_payload),
    pvstackr:::pv_hash_payload(materialized_payload)
  )
})
