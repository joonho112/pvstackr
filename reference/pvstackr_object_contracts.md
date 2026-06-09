# pvstackr Object Contracts

This page records the package-level object contracts used by the current
public API. These contracts are intentionally conservative while
`pvstackr` is in its v0.1 method stage.

## Current Method Boundary

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
recognizes three public method IDs: `"stack_direct"`, `"stack_psis"`,
and `"per_pv"`. All three are implemented in this package stage through
injected or precomputed light paths. Only `stack_direct` carries a
design-based external Rubin/BRR-Fay target. `per_pv` carries a
`pvstackr_reference_pool` built from per-PV backend fixed-effect draws;
its within-PV variance is the model-based covariance of the selected
posterior draws. `stack_psis` carries no formal target; its estimates
are diagnostic/reference Rubin pooling of PSIS-weighted stacked draw
summaries using model-based weighted covariance.

## `stack_direct` v0.1 Boundary

`stack_direct` requires an external `pvstackr_brr_target` from
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
The wrapper checks formula RHS equality, fixed-effect name alignment,
target provenance, and target policy before fitting. Group terms such as
`(1 | school)` and `(1 || school)` are not accepted for the v0.1 BRR-Fay
target path because the current target engine is fixed-effect-only.

## `pvstackr_ccc` Contract

A `pvstackr_ccc` object stores deterministic fixed-effect CCC
calibration from raw stacked draws to the external BRR-Fay/Rubin
fixed-effect target. Its reportable fields include calibrated
fixed-effect draws, raw and target centers, raw and target covariance
matrices, the calibration matrix, target hash, and structured
diagnostics. The diagnostic `delta_c_rel` is the RMS target-SE-scaled
raw-center to target-center separation across fixed-effect terms;
`delta_c_max` is the maximum absolute target-SE-scaled separation and is
the NO-SEND gate metric. `a_matrix_fro_rel` is the relative Frobenius
distance of the calibration matrix from identity and is descriptive
rather than a NO-SEND gate. Calibration-matrix conditioning is gated by
`kappa_A`: green for `kappa_A < 1e6`, yellow for `1e6 <= kappa_A < 1e8`,
and red for `kappa_A >= 1e8`. Center separation bands are green for
`delta_c_max < 1e-2`, yellow for `1e-2 <= delta_c_max < 5e-2`, and red
for `delta_c_max >= 5e-2`. The diagnostics always compare the raw
stacked fixed-effect center to the external target center, even when the
lower-level `center = "posterior"` convention leaves the output draws
centered at the raw posterior mean. That posterior-centered CCC
convention is diagnostic/exploratory only; reportable `stack_direct`
`pvstackr_fit` output requires `control$center = "target"`. Yellow
center separation or conditioning becomes a top-level `pvstackr_fit`
warning; red center separation or conditioning blocks reportable
estimates.

## `pvstackr_fit` Contract

A reportable `stack_direct` `pvstackr_fit` with status `"ok"` or
`"warning"` must contain non-null `design`, `target`, `stack_fit`, and
`ccc` components; a non-empty fixed-effect estimate table; and
calibrated fixed-effect draws when `control$return_draws = TRUE`.
Reportable `stack_direct` fits must use `control$center = "target"` so
fixed-effect estimates are the external Rubin/BRR-Fay target centers.
Top-level `warnings` are reserved for interpretation-level diagnostics
such as explicit priors or yellow center separation. Routine backend
draw-column drops, such as sampler columns named `lp__`, are retained in
nested stack-fit diagnostics instead of becoming top-level fit warnings.

A reportable `per_pv` `pvstackr_fit` with status `"ok"` or `"warning"`
must contain a non-null `pvstackr_reference_pool` target and a non-empty
fixed-effect estimate table. `per_pv` fits do not carry `stack_fit` or
`ccc` components. The within-PV variance component is the model-based
posterior covariance of the selected per-PV fixed-effect draws, not a
BRR/Fay replicate variance. `weight_col`, `rep_weight_cols`, `fay_k`,
and `id_cols` are provenance/design metadata for this path; they are not
automatically supplied to `fit_function` and do not affect Rubin
pooling. Weighted backend fits must be handled by the backend adapter,
data supplied to the backend, or backend-specific `additional_args`.
Per-PV draws, when retained, remain nested in diagnostics rather than
top-level reportable draws because Rubin pooling does not synthesize a
single calibrated top-level draw matrix.

A reportable `stack_psis` `pvstackr_fit` with status `"ok"` or
`"warning"` must contain PSIS diagnostics, pooling diagnostics, weighted
per-PV summaries, and a non-empty fixed-effect estimate table. The
pooling is diagnostic/reference Rubin pooling of PSIS-weighted
fixed-effect summaries using model-based weighted covariance; this path
does not construct an external BRR-Fay target. Group terms such as
`(1 | school)` and `(1 || school)` are rejected for this v0.1
fixed-effect-only path. Failed Pareto-k diagnostics cannot be reported
as `status = "ok"`. With the default block fallback, failed PSIS
diagnostics produce a blocked fit with diagnostics but no reportable
estimates.

## `pvstackr_method_comparison` Contract

A `pvstackr_method_comparison` object compares two or more
already-created `pvstackr_fit` objects. It contains aligned fixed-effect
rows in `estimate_table`, method-level fields in `diagnostic_table`,
agreement diagnostics in `agreement`, and elapsed-time metadata in
`timing`. Method IDs must remain one of `"stack_direct"`,
`"stack_psis"`, or `"per_pv"`, while method labels are stable unique
labels used to identify compared fits. Blocked methods are retained with
reason codes and `NA` comparison statistics instead of being dropped.
The aligned `estimate_table` also preserves interval metadata from each
fit: `df_method`, `df_complete`, `conf_level`, `interval_role`, and
`coverage_claim_allowed`. Method-level diagnostics summarize the
interval role and count how many reportable rows are descriptive rather
than coverage-claimable. They also preserve available target and pooling
provenance (`target_source`, `target_hash`, `pooling_source`, and
`pooling_hash`) and shared-provenance flags (`shared_target_hash`,
`shared_pooling_hash`, and `shared_external_target`). The
`shared_external_target` flag is narrow: it is true only when two or
more methods share the same non-missing `external_brr_fay_rubin` target
hash. Broader overlap, including shared target-source families or shared
reference hashes, is summarized in `diagnostics$target_overlap`.
Agreement diagnostics are descriptive and should not be read as
automatic independent corroboration when compared methods share a target
hash, pooling hash, target source, or estimand construction.

## Accessor Contract

Public accessors expose stable object fields without changing their
semantics.
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
returns reportable fixed-effect estimate tables, including
interval/provenance columns such as `df_method`, `df_complete`,
`interval_role`, `coverage_claim_allowed`, `target_source`,
`target_hash`, `pooling_source`, and `pooling_hash` when those columns
are part of the fit contract.
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
returns the formal fit target object when a method has one, and `NULL`
for methods such as `stack_psis` that do not carry a formal target.
Estimate-row `target_source` labels are provenance metadata; they are
not a guarantee that
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
returns a non-null object.
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
returns retained top-level reportable draws only; per-PV reference draws
and PSIS weights remain diagnostic artifacts accessible through
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).

## Reporting Scope

In v0.1, calibrated reporting is fixed-effect-only. Variance-component
or sampler-diagnostic columns may be retained inside backend or
stack-fit components, but they are not calibrated to the Rubin/BRR-Fay
target and are not included in the reportable estimate table.

## Interval Metadata

The interval-role vocabulary is method-specific:

- `descriptive_classic_rubin`: `stack_direct` rows backed by an external
  `external_brr_fay_rubin` target using classic Rubin imputation df.
  These intervals are descriptive, and `coverage_claim_allowed = FALSE`.

- `coverage_barnard_rubin`: `stack_direct` rows backed by an external
  `external_brr_fay_rubin` target using Barnard-Rubin df with explicit
  `df_complete`. These rows set `coverage_claim_allowed = TRUE`.

- `reference_classic_rubin` and `reference_barnard_rubin`: `per_pv`
  reference-pooling rows based on model-based posterior-draw covariance.
  These rows always set `coverage_claim_allowed = FALSE`, even when
  Barnard-Rubin df and `df_complete` are supplied.

- `psis_classic_rubin` and `psis_barnard_rubin`: `stack_psis` diagnostic
  pooling rows based on PSIS-weighted model-based covariance. These rows
  always set `coverage_claim_allowed = FALSE`, even when Barnard-Rubin
  df and `df_complete` are supplied.

Thus `df_method = "barnard_rubin"` does not by itself create a
coverage-claimable interval. In this package stage, coverage claims are
enabled only for `stack_direct` rows backed by the external BRR-Fay
target with `interval_role = "coverage_barnard_rubin"`. Default fit and
method-comparison print methods emit a one-line interval note whenever
reportable rows are descriptive rather than coverage-claimable.

## Source Vocabulary

`target_source = "external_brr_fay_rubin"` denotes a formal external
design-based target object. `target_source = "per_pv_rubin_draws"`
denotes a `per_pv` reference pool built from per-PV backend draw
summaries. `stack_psis` fits do not carry a formal target object; their
estimate rows use `target_source = "stack_psis_rubin_pooling"` and
`pooling_source = "stack_psis_rubin_pooling"` to record diagnostic
pooling provenance. PSIS input-source diagnostics use a separate
vocabulary: `supplied_weights`, `psis_function`, and
`log_ratios_self_normalized`. The `log_ratios_self_normalized` path
self-normalizes the supplied log ratios and does not by itself run
Pareto smoothing.
