# pvstackr Object Contracts

This page records the package-level object contracts used by the current
public API. These contracts are intentionally conservative while
`pvstackr` is in its fixed-effect method stage.

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

## `stack_direct` Boundary

`stack_direct` requires an external `pvstackr_brr_target` from
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
The wrapper checks formula RHS equality, fixed-effect name alignment,
target provenance, and target policy before fitting. Group terms such as
`(1 | school)` and `(1 || school)` are not accepted for the BRR-Fay
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
Top-level calibrated fixed-effect draws are the only retained
individual-draw matrix: nested raw stacked draws, nuisance draws, full
calibrated draws, and duplicate calibrated fixed-effect matrices are
removed for either retention setting. Blocked fits retain no individual
draws and record effective `return_draws = FALSE`. Top-level `warnings`
are reserved for interpretation-level diagnostics such as explicit
priors or yellow center separation. Routine backend draw-column drops,
such as sampler columns named `lp__`, are retained in nested stack-fit
diagnostics instead of becoming top-level fit warnings.

Every current `pvstackr_fit` also carries an exact `validation` record
with a validation-schema ID, semantic-policy ID, canonicalizer ID,
fast-path eligibility flag, and a domain-separated SHA-256 stamp. Deep
validation evaluates the full method contract and then verifies this
stamp. Cheap validation still re-hashes the current package-owned
payload; it never trusts the stored stamp alone. Fits that retain opaque
backend objects are not fast-path eligible and fall back to deep
validation because backend implementation state is outside the portable
stamp boundary. The stamp is tamper-evident for stale or uncoordinated
mutations, not an external digital signature or proof of authorship.

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
Per-PV draws, when retained, remain fixed-effect-only matrices nested in
diagnostics rather than top-level reportable draws because Rubin pooling
does not synthesize a single calibrated top-level draw matrix.

A reportable `stack_psis` `pvstackr_fit` with status `"ok"` or
`"warning"` must contain PSIS diagnostics, pooling diagnostics, weighted
per-PV summaries, and a non-empty fixed-effect estimate table. The
pooling is diagnostic/reference Rubin pooling of PSIS-weighted
fixed-effect summaries using model-based weighted covariance; this path
does not construct an external BRR-Fay target. Group terms such as
`(1 | school)` and `(1 || school)` are rejected for this
fixed-effect-only path. Failed Pareto-k diagnostics cannot be reported
as `status = "ok"`. With the default block fallback, failed PSIS
diagnostics produce a blocked fit containing only the canonical
scalar/vector Pareto-k decision record and its redaction manifest, with
no reportable estimates, target, design, stack fit, CCC, or draws. All
heavy-retention controls are recorded as effective `FALSE` for a blocked
fit, irrespective of the original request. When `return_draws = TRUE` on
a reportable fit, diagnostics retain the fixed-effect-only proposal
matrix together with its normalized PV weight matrix; these two payloads
are retained or removed as a pair. Package-owned nuisance draws and
nested full stack draws are never retained in a final PSIS fit.
Explicitly authorized opaque backend objects and log-likelihood matrices
remain governed separately by `keep_backend_fit` and `keep_log_lik`.

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
and the PSIS fixed-effect proposal/weight pair remain diagnostic
artifacts accessible through
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
Current fit accessors and fit print/summary methods use the rehash
validation tier: they rescan the current package-owned payload and
compare its SHA-256 stamp without repeating method-specific semantic
recomputation. This is linear in retained payload bytes rather than an
O(1) cache lookup. Opaque-backend fits remain ineligible and fall back
to deep validation.

Historical `stack_psis` fits that do not satisfy the current fit
envelope must be passed to
[`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md).
The migrator never promotes saved estimates, pooling summaries, weights,
or draws. It returns an explicit inspection-only object containing
bounded Pareto-k decision evidence and a redaction manifest;
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
and
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
refuse that object, while
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
remains available. Current comparisons and fit or comparison summaries
carry source and owned-payload validation stamps. Pre-marker serialized
comparisons/summaries containing `stack_psis` are refused and must be
rebuilt, closing derived-table paths that no longer carry their source
fit.

Current live `stack_direct` fits expose normalized sampler diagnostics
and a frozen reportability gate through
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
A sampler- or CCC-blocked fit retains only its independently valid
external BRR-Fay target plus the exact slim gate evidence required to
reproduce the blocked status; CCC evidence is scalar-only. Design data,
stack/backend fit, full CCC matrices, estimates, and draws are absent,
and every heavy-retention control is recorded as effective `FALSE`. The
target is rebuilt from an exact recursive allowlist with empty warnings
and a safe formula environment; the slim preflight snapshot does not
retain a formula object. A recursive fail-closed validator rejects
result-like fields, matrices/arrays/data frames, raw or complex
payloads, executable/external objects, hidden leaf attributes, and
noncanonical diagnostic variants anywhere outside that independent
target. Current `per_pv` fits have no typed blocked-object path, so a
reportable reference fit cannot be relabeled as blocked. Legacy cached
schema-0.1 stack fits remain readable without the newer sampler gate.

## Reporting Scope

Calibrated reporting is fixed-effect-only. Variance-component or
sampler-diagnostic columns may be retained inside backend or stack-fit
components, but they are not calibrated to the Rubin/BRR-Fay target and
are not included in the reportable estimate table.

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
vocabulary: `supplied_psis_weights`, `injected_psis_function`, and
`self_normalized_log_ratios`. `pareto_k_source` distinguishes supplied
diagnostics from injected-function output, while `weight_method`
distinguishes caller-declared external PSIS from unspecified external
weights and self-normalized raw importance weights. Only the
caller-declared path can expose estimates; the declaration records
producer/version but is not package verification. Each path records
bounded per-PV Kish-style iid weight ESS (`weight_ess_iid`), its
draw-count fraction (`weight_ess_fraction`), and
`max_normalized_weight`. These concentration diagnostics are not MCMC
ESS and never replace the immutable Pareto-k gate. If normalized weights
are retained, deep validation recomputes these diagnostics. A compact or
blocked object instead records
`weight_diagnostic_authority = "owned_stamp_bounded_projection"`:
feasibility and the owned payload stamp remain checkable, but the
redacted original weights cannot be independently reconstructed.
