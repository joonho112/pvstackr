test_sampler_payload <- function(chains = 2L, post_warmup = 5L, ...) {
  overrides <- list(...)
  payload <- list(
    rhat_max = 1.005,
    ess_bulk_min = 110 * as.integer(chains),
    ess_tail_min = 105 * as.integer(chains),
    divergences = 0L,
    chains = as.integer(chains),
    post_warmup_draws_per_chain = as.integer(post_warmup)
  )
  payload[names(overrides)] <- overrides
  payload
}

test_sampler_diagnose_function <- function(chains = 2L, post_warmup = 5L,
                                           extra = list(), ...) {
  payload <- test_sampler_payload(chains, post_warmup, ...)
  force(extra)
  function(fit) {
    c(extra, list(sampler = payload))
  }
}
