# pvstackr 0.2.0

* Revalidating a legacy target no longer requires its estimates to match
  bit-for-bit. `pv_revalidate_brr_target()` rebuilds the target from the raw
  inputs and compared every primitive with `identical()`, including the
  weighted least-squares estimates and a hash of them. Those last bits differ
  between BLAS implementations, so a target built on one machine could not be
  revalidated on another. Structure, declared inputs, and the input-derived
  design hash are still compared exactly; the estimates now go through the
  package's existing numeric policy (absolute 1e-12, relative 1e-10). A
  difference large enough to matter still fails closed.

* Content hashes and validation stamps no longer depend on the R version that
  wrote them. Both hashed the whole serialization stream, whose fourteen-byte
  header records the writing R version. The eight-hex content hashes
  (`design_hash`, `target_hash`, and the fixture digests) are compared
  fail-closed when a legacy target is revalidated, so that comparison could
  never succeed across two R versions, and the validation stamps made a saved
  fit refuse its own accessors after an R upgrade. The header is now excluded
  from both. **Hash values change**: content hashes computed by earlier
  releases will not match the ones this release derives, so a stored
  `design_hash` from 0.1.x or 0.2.0-dev is expected to differ.

* Validation stamps no longer depend on the R version that wrote them. The
  stamp hashed the whole serialization stream, whose fourteen-byte header
  records the writing R version, so a saved fit failed its own stamp when it
  was read back under a different R version and its accessors refused to
  return anything. The header is now excluded, and the four stamps -- fit,
  summary, method comparison, and comparison source projection -- are portable
  across R versions and machines. The bundled example fit was regenerated
  under the corrected scheme.

* The bundled brms adapter is now public: `pv_backend_brms_fit_function()`,
  `pv_backend_brms_draws_function()`, and
  `pv_backend_brms_sampler_diagnostics()`. Injecting all three into `pv_fit()`
  runs the same fitting, extraction, and diagnostic code the bundled backend
  runs, so the two routes are one implementation rather than two that have to be
  kept in step, and any one of the three can be replaced to attach a different
  engine. The routes are not interchangeable in what they record: an injected
  fit carries injected engine and diagnostic provenance, and the bundled route
  additionally checks backend availability before sampling. A reportable
  injected fit needs all three; an adapter that supplies only a fit and a draws
  function leaves the sampler evidence incomplete and the fit is blocked.

* `pv_backend_brms_fit_function()` now resolves the Stan backend itself when it
  is reached through an injected adapter, instead of rejecting every
  `pv_control(backend = )` value the bundled route would have resolved for it.
  Cache directories remain the injected adapter's own responsibility, as the
  recorded provenance has always said; create the directory first or pass
  `cache_dir = NULL`. Under `backend = "brms"` a supplied `diagnose_function` is
  still ignored in favour of the bundled sampler record, but the fit now always
  records that the override was ignored, not only when the override happened to
  return a `sampler` element.

* A population-level prior (`class = "b"` with no `coef`) is now expanded onto
  the individual slope coefficients instead of being refused. `stack_direct`
  materializes the model matrix and fits `0 + ...`, which turns the intercept
  into an ordinary column, so passing such a prior through unchanged would have
  silently widened it to cover the intercept. Expanding it restores the scope
  the prior had against the original formula. Priors whose scope cannot be
  reproduced exactly -- intercept-class, coefficient-specific, grouped, or
  opaque -- are still refused before the backend runs, and supplying any
  explicit prior still raises the `explicit_prior_warning` reason code.

* `stack_psis` now separates the weight input route, Pareto-k source, and
  smoothing provenance. Supplied or injected weights can enter the reportable
  caller-declared external-PSIS route only when the caller records both an
  external producer and version; this records a declaration, not package
  verification;
  otherwise they fail closed, as do raw self-normalized log ratios. Every path
  records per-PV Kish-style iid weight ESS, its draw-count fraction, and the
  maximum normalized weight. These are weight-concentration diagnostics, not
  MCMC ESS, `loo::relative_eff()`, or a replacement for the immutable Pareto-k
  gate. Mixed weight routes are rejected rather than silently prioritized.
  The legacy `fallback = "warn"` request now emits its scheduled deprecation
  warning and still behaves as immutable fail-closed `"block"`.

* New bundled brms backend for the stacked fit: `pv_fit(method =
  "stack_direct")` with `pv_control(backend = "brms")` now runs without an
  injected `fit_function`/`draws_function` pair. The bundled adapter fits the
  prepared stacked formula with `brms::brm()` (through cmdstanr when it is
  installed and CmdStan is configured, rstan when cmdstanr is absent; an
  installed cmdstanr without a working CmdStan is an error, not a fallback) and returns fixed-effect and residual-scale
  draws as a plain base matrix. brms is in Suggests; injected adapters remain
  the mechanism for other engines, and `per_pv`/`stack_psis` inputs are
  unchanged.

* Final `pvstackr_fit` objects now carry an exact validation record. A deep
  tier preserves the full semantic checks, while a cheap tier re-hashes the
  current package-owned payload with domain-separated SHA-256 and rejects
  stale stamps. This catches same-moment draw replacements, coordinated
  proposal/weight row changes, hidden leaf attributes, and self-rehashed
  data-free design snapshots that summary-only validation cannot distinguish.
  Opaque retained backend objects are explicit fast-tier exceptions and force
  deep validation. Public fit accessors and fit print/summary methods use the
  rehash tier, avoiding repeated method-specific recomputation while still
  reading every mutable retained payload byte.

* Blocked fits now use a generic fail-closed retention firewall. Effective
  `return_draws`, `keep_data`, `keep_backend_fit`, and `keep_log_lik` are always
  `FALSE`; reportable estimates and heavy design, stack-fit, CCC, draw, weight,
  pooling, backend, and log-likelihood payloads are absent. Blocked
  `stack_direct` fits may retain only a canonical independently valid external
  target plus exact slim sampler or scalar CCC gate evidence, while blocked
  `stack_psis` fits retain only the canonical Pareto-k/provenance decision
  record and bounded weight-concentration diagnostics. The
  recursive validator also rejects hidden payloads in nested fields or leaf
  attributes. Because `per_pv` has no typed blocked schema, relabeling a
  reportable reference fit as blocked is rejected.

* Historical PSIS results now have an explicit inspection-only migration path.
  `pv_migrate_legacy_psis_fit()` returns current validated fits unchanged but
  projects any non-current `stack_psis` fit to bounded Pareto-k evidence and a
  redaction record; saved estimates, draws, pooling, weights, backend objects,
  and data are not migrated, and estimate/draw accessors refuse the inspection
  object. New method comparisons and fit/comparison summaries carry deep-valid
  compact source reportability objects, canonical source projections, and
  owned-payload SHA-256 records. Pre-marker serialized comparisons or summaries
  containing PSIS are refused and must be rebuilt, while an independent semantic
  gate also rejects warning PSIS rows and any blocked row retaining numeric or
  pooling metadata.

# pvstackr 0.1.1

Patch release fixing two bugs found while preparing the LSAE software article.
Point estimates are unchanged for every method; the fixes affect `stack_psis`
interval width and the acceptance of `posterior::draws_matrix` inputs.

## Bug fixes

- **`stack_psis` PSIS-weighted covariance.** `pv_weighted_mean_cov()` applied
  `sqrt(weights)` to only one factor of the cross-product
  (`crossprod(centered * sqrt(weights), centered)`), weighting each draw by
  `sqrt(w)` instead of `w`. This inflated `stack_psis` within-imputation
  covariances by roughly the square root of the draw count. Weighted-mean point
  estimates were unaffected, but `stack_psis` standard errors and intervals
  reported by 0.1.0 are inflated and should be recomputed. Added a direct
  regression test for the weighted-covariance formula.

- **CCC validation of `posterior::draws_matrix` inputs.** When a
  `draws_function` returned a `posterior::draws_matrix` (the natural output of
  `posterior::as_draws_matrix()` on a `brmsfit`), `validate_pvstackr_ccc()`
  compared `draws_fe_cal` against the fixed-effect block of `draws_calibrated`
  with a strict `identical()`. The `draws_matrix` S3 class and draw-id row
  names made that check fail even when the values matched, aborting with
  "CCC `draws_fe_cal` must equal the fixed-effect block of `draws_calibrated`".
  Draw matrices are now normalized to plain base matrices at the ingestion
  boundary (`ccc_as_draw_matrix()`), so the identity check compares values and
  column names rather than S3 provenance. Added regression tests for a
  `draws_matrix`-classed input through calibration and validation, plus (under
  the optional backend-smoke suite) a genuine `posterior::as_draws_matrix()`
  through the full `stack_direct` CCC path.

# pvstackr 0.1.0

First public release. `pvstackr` calibrates a stacked Bayesian-backend fit of
plausible-value survey data to an external Rubin/BRR-Fay **fixed-effect** target,
anchoring the reportable fixed-effect estimates and intervals to a design-based
reference.

## New features

- **Declare a PISA-style design.** `pv_design()` describes a plausible-value
  survey design and auto-detects its structure with
  `detect_pisa_pv_columns()` (plausible-value sets) and
  `detect_pisa_brr_replicate_weights()` (BRR replicate weights), recording
  row-support and Fay/design fingerprints and validating duplicate, missing, or
  malformed columns.
- **Assemble an external target.** `pv_brr_target()` builds the external
  Rubin/BRR-Fay fixed-effect target (Rubin combining, the BRR-Fay sandwich
  variance, and Barnard-Rubin degrees of freedom) that `stack_direct`
  calibrates to.
- **Fit a method.** `pv_fit()` dispatches across three methods:
  `stack_direct` (the default), `per_pv`, and `stack_psis`. The engines are
  also exported directly as `pv_fit_direct()`, `pv_fit_reference()`, and
  `pv_fit_stack_psis()`. `pv_control()` collects fitting controls.
- **Compare methods.** `pv_compare_methods()` aligns the fixed-effect estimates
  from the three methods into a single comparison object with agreement
  diagnostics, retains blocked methods explicitly, and validates the result.
- **Read results.** Stable accessors return the pieces you report:
  `get_estimates()` (the reportable estimate table with interval metadata),
  `get_target()`, `get_draws()` (retained top-level fixed-effect draws), and
  `get_diagnostics()`.
- **Inspect at a glance.** Compact `print()` and `summary()` methods cover
  designs, targets, controls, fits, and method comparisons.
- **Bundled synthetic example.** A tiny **synthetic** `pisa_tiny` fixture and a
  cached `stack_direct` fit (under `inst/extdata/`) let every example and
  vignette run instantly with no live sampler. No real PISA records or
  OECD/PISA-distributed files are bundled.
- **Documented object contract.** `?pvstackr_object_contracts` specifies the
  structure of each object and the meaning of the interval-metadata columns.

## Documentation

- **Ten vignettes in two tracks.** An **Applied track** (workflow-first):
  A1 getting started, A2 the end-to-end workflow, A3 reading results and what
  to report, A4 comparing methods, and A5 real PISA data guidance. A
  **Method track** (rigorous, equation-backed): M1 foundations and notation,
  M2 the BRR-Fay fixed-effect target, M3 the stacked fractional bridge,
  M4 CCC (Cholesky Calibration Correction), and M5 methods, PSIS, and coverage.
- **pkgdown website** at <https://joonho112.github.io/pvstackr/>, with the
  reference index grouped by workflow stage (design -> target -> fit ->
  compare -> read).
- **Rewritten roxygen** documentation with `@family` groups and `@seealso`
  cross-links connecting each stage of the workflow.
- **Rewritten README** with a problem-first introduction, the method in three
  steps, a five-minute fixture example, and both vignette tables.
- **Hex logo** and favicons for the package and site.

## Scope and honesty

- **Fixed-effect only.** Reportable estimates and retained draws are
  fixed-effect quantities. Group/random-effect terms are out of scope for the
  v0.1 BRR-Fay target engine.
- **Coverage is claimable for `stack_direct` only.** `stack_direct` intervals
  are coverage-claimable only when backed by an external design-based BRR-Fay
  target with Barnard-Rubin degrees of freedom. The claim rests on target
  provenance and validation evidence, not on CCC arithmetic alone. `per_pv` and
  `stack_psis` results are descriptive, and the estimate table flags this per
  interval.
- **No efficiency or speed claims.** "One fit" refers to the stacked-fit
  architecture (topology), not to benchmarked runtime. The package makes no
  "faster" or "N-times" claims.
- **No real PISA microdata.** Only the synthetic `pisa_tiny` fixture ships.
  Real PISA data must be obtained under its own license; vignette A5 explains
  how to connect to it without bundling licensed files.
- **`stack_psis` expects supplied diagnostics.** It uses PSIS weights and
  Pareto-k diagnostics that you supply or inject; the package does not run live
  `loo::psis()` in v0.1.
- **Not affiliated.** This is an independent research package and is not
  affiliated with or endorsed by the OECD or the PISA programme.
- **Method paper forthcoming.** The methodology manuscript is in preparation;
  for now, please cite the package (see `citation("pvstackr")`).


# pvstackr 0.0.0.9000

## Development

- Initial development scaffold: established the package boundary as a clean
  applied-method package (not a paper-replication repository), fixed the
  `pvstackr` name and the `per_pv` / `stack_psis` / `stack_direct` method IDs,
  and built out the public surface, fixtures, and test/parity coverage that the
  0.1.0 release documents.
