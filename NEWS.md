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
