# Objects returned by pvstackr

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
and the method functions
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)
and
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)
return a fit, a list of class `pvstackr_fit`;
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
returns a list of class `pvstackr_method_comparison`. A fit holds the
estimate table, the target, the draws and the diagnostics of one method,
and a `status` that says whether estimates are reported. Read these
objects with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
rather than with `$`.

## Methods

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
has three methods, and all of them report only the fixed effects
([`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
gives the scope and the interval rule). Only `stack_direct` has a
bundled engine and a design-based target.

- `"stack_direct"` fits one model to the stacked data, with the bundled
  brms engine (`pv_control(backend = "brms")`) or with fitting functions
  that you supply, and calibrates its fixed-effect draws to the target
  from
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  (BRR-Fay replicate weights and Rubin's rules). All `stack_direct` fits
  must use `control$center = "target"` (the default), so the estimates,
  standard errors and degrees of freedom are those of the target. The
  other value, `"posterior"`, is accepted by
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  but used by no fitting function:
  [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  stops with an error when it is set.

- `"per_pv"` takes the draws of one fit per plausible value, made by
  your fitting functions or computed elsewhere, and combines them with
  Rubin's rules. The variance within each plausible value comes from the
  posterior draws, not from the replicate weights. The survey weights
  are not automatically passed to your `fit_function`; a weighted fit
  must take them from the data inside that function or through
  `additional_args`. See
  [`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md).

- `"stack_psis"` reweights one set of stacked draws, which you supply or
  which your fitting functions produce, toward each plausible value with
  importance weights that you supply, then combines the weighted results
  with Rubin's rules. See
  [`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md).

## The fit object (`pvstackr_fit`)

A fit is a list with these components:

- `method`: `"stack_direct"`, `"per_pv"` or `"stack_psis"`.

- `status`, `reason_codes`, `warnings`: the status of the fit, the codes
  that explain a `"warning"` or `"blocked"` status, and the matching
  messages (see "Status and checks"). An `"ok"` fit has no reason codes
  and no warnings.

- `estimates`: the estimate table (see "Estimate table"); empty when the
  fit is blocked.

- `target`: for `stack_direct`, the `pvstackr_brr_target` object from
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md);
  for `per_pv`, a `pvstackr_reference_pool` object that holds the
  Rubin's-rules combination of the per-plausible-value draws; `NULL` for
  `stack_psis`.

- `draws`: for `stack_direct` with `return_draws = TRUE` (the default of
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)),
  the matrix of calibrated fixed-effect draws; otherwise `NULL`. With
  `return_draws = TRUE`, `per_pv` keeps its draws in
  `get_diagnostics(fit)$reference$per_pv_draws`, and `stack_psis` keeps
  the stacked fixed-effect draws and the normalized weights in
  `get_diagnostics(fit)$weighted`.

- `design`: the design record (the plausible-value and weight columns,
  the numbers of plausible values and replicate weights, and checksums).
  It holds the data only when `keep_data = TRUE`, and it is `NULL` for a
  blocked fit and when `data` and `formula` were not both supplied.

- `stack_fit`: a record of the stacked fit (engine, settings and sampler
  diagnostics) for `stack_direct`, and for `stack_psis` when it fitted
  or received a stacked fit. It holds the engine's fitted model only
  when `keep_backend_fit = TRUE` and log-likelihood draws only when
  `keep_log_lik = TRUE`; the stacked draws are not kept.

- `ccc`: the calibration record of a `stack_direct` fit (see
  "Calibration record").

- `diagnostics`: a list whose names depend on the method. For
  `stack_direct`: `preflight` (the record of the check that the formula,
  the data and the target match), `sampler` and `sampler_gate` (the
  sampler diagnostics and their check), `stack_fit`,
  `stack_fit_warnings` and `ccc`. Messages from the stacked fit, such as
  a note that the engine's `lp__` column was dropped, go to
  `stack_fit_warnings`, not to `warnings`. For `per_pv`: `reference` and
  `pooling`. For `stack_psis`: `psis`, `pooling` and `weighted`.

- `control`: the
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  settings.

- `schema_version`, `provenance`, `validation`: see "Technical details".

## Status and checks

`status` is `"ok"` when no check gave a warning or a block; it is not a
statement that the sampler converged. With `"warning"`, estimates are
returned and `reason_codes` and `warnings` say which checks gave the
warning. With `"blocked"`, the estimate table is empty and
`reason_codes` says why. The most severe result of all checks sets the
status. The checks and their thresholds are pvstackr's rules:

- Sampler diagnostics (`stack_direct`, checked before calibration; a fit
  blocked here is not calibrated). The six values `rhat_max`,
  `ess_bulk_min`, `ess_tail_min`, `divergences`, `chains` and
  `post_warmup_draws_per_chain`, from the bundled engine or from your
  `diagnose_function`, must all be available, and the last two must
  equal `chains` and `iter - warmup` in `control`; otherwise the fit is
  blocked. The largest R-hat gives a warning above 1.01 and a block
  above 1.05. For the bulk and the tail effective sample size (ESS),
  pvstackr takes the smallest ESS over the parameters (a total over all
  chains) and divides it by the number of chains: a total below 100 or a
  per-chain value below 25 blocks, and a per-chain value below 100 gives
  a warning. Any divergent transition blocks. These reason codes start
  with `sampler_`.

- Center separation (`stack_direct`). `delta_c_max` is the largest
  absolute difference, over the fixed effects, between the target
  estimate and the mean of the stacked draws before calibration, divided
  by the target standard error. It gives a warning from 0.01
  (`center_separation_yellow`) and a block from 0.05
  (`center_separation_red`). `delta_c_rel`, the root mean square (RMS)
  of the same ratios, is recorded but does not change the status.

- Conditioning (`stack_direct`). `kappa_A`, the condition number of the
  calibration matrix (see "Calibration record"), gives a warning from
  1e6 (`ccc_conditioning_yellow`) and a block from 1e8
  (`ccc_conditioning_red`).

- Priors (`stack_direct`). Any explicit prior, even a flat one, gives a
  warning (`explicit_prior_warning`).

- Pareto k-hat (`stack_psis`). The fit is blocked (`psis_k_too_high` or
  `psis_k_not_evaluated`) unless every plausible value has a finite
  Pareto k-hat below `control$psis_k_threshold` (default 0.7, also the
  maximum), and blocked (`psis_weight_provenance_incomplete`) unless the
  weights come with `psis_producer` and `psis_producer_version`. Weights
  that pvstackr makes from `log_ratios` are not smoothed and are always
  blocked (`psis_smoothing_not_applied`). A `stack_psis` fit never has
  status `"warning"`.

- `per_pv` has no checks; its status is always `"ok"`.

A blocked fit keeps only what explains the block, and its `control`
records `return_draws`, `keep_data`, `keep_backend_fit` and
`keep_log_lik` as `FALSE`. A blocked `stack_direct` fit keeps the
target, the diagnostics `preflight`, `sampler` and `sampler_gate` (and
`ccc`, reduced to single values grouped as `center`, `conditioning`,
`residual` and `prior`, when the calibration check blocked it) and a
`redaction` list that names what was removed; it has no data, stacked
fit, calibration matrices, estimates or draws. The diagnostics of a
blocked `stack_psis` fit hold only `psis` and `redaction`; it has no
target, design, stacked fit, estimates, weights or draws.

## Calibration record (`pvstackr_ccc`)

The `ccc` component of a `stack_direct` fit, of class `pvstackr_ccc`,
records the Cholesky calibration correction (CCC). `psi_raw` and
`Sigma_raw` are the mean and covariance of the stacked fixed-effect
draws before calibration; `psi_target` and `Sigma_target` are the target
estimates and the target covariance (`T_MI`). With \\L\_{raw}\\ and
\\L\_{tgt}\\ the lower Cholesky factors of `Sigma_raw` and
`Sigma_target`, the calibration matrix is \\A = L\_{tgt} L\_{raw}^{-1}\\
(`A`). Each stacked draw is centered at `psi_raw`, multiplied by `A` and
shifted to `psi_target`, so the calibrated draws have mean `psi_target`
and covariance `Sigma_target`; they are kept only in the fit's `draws`.

`ccc$diagnostics`, which
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
also returns as `ccc`, holds `delta_c_max`, `delta_c_rel` and `kappa_A`
(see "Status and checks") and values that are recorded for inspection
and do not change the status: `a_matrix_fro_rel` (how far `A` is from
the identity matrix) and `rho1`, `rho2` and `empirical_fro_rel` (how far
the covariance of the calibrated draws is from `Sigma_target`).

## Estimate table

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
returns one row per fixed effect. All methods give these columns:

- `term`: the coefficient name with the prefix `b_`, such as
  `b_Intercept`.

- `estimate`, `se` and `df`: the estimate, its standard error and its
  degrees of freedom. `std.error` repeats `se`.

- `df_method`: `"classic"` or `"barnard_rubin"`. `df_complete`: the
  complete-data degrees of freedom that you gave for the Barnard-Rubin
  rule; `NA` with classic degrees of freedom (a `stack_psis` fit copies
  a value given with `"classic"` into this column without using it).

- `conf_level`, `conf_low` and `conf_high`: the interval level (from
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
  default 0.95) and the interval, `estimate` plus or minus a t quantile
  with `df` degrees of freedom times `se`. `conf.low` and `conf.high`
  repeat the interval.

- `interval_role` and `coverage_claim_allowed`: how the interval can be
  read (see below).

- `parameter_scope`: always `"fixed_effect"`.

- `target_source` and `target_hash`: see "Source labels".

`per_pv` and `stack_psis` tables add `pooling_source` and
`pooling_hash`; `stack_psis` tables also add `psis_status`,
`pareto_k_max`, `psis_k_threshold`, `psis_source`, `pareto_k_source`,
`weight_method`, `psis_producer` and `psis_producer_version`. There is
no column for the fraction of missing information; for `stack_direct`
and `per_pv` fits it is `get_target(fit)$fmi`. For `stack_direct`,
`estimate`, `se` and `df` are the target's (see "Methods"); for `per_pv`
and `stack_psis` they come from the Rubin's-rules combination of the
per-plausible-value results.

`coverage_claim_allowed` is `TRUE` when pvstackr's reporting rule lets
you read the interval as a confidence interval with nominal coverage and
`FALSE` when the interval is descriptive; the label records how the
interval was built and does not certify its coverage. pvstackr sets
`interval_role` and `coverage_claim_allowed` by the rule in
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md):

- `coverage_barnard_rubin`: `stack_direct` with a target built with
  `df_method = "barnard_rubin"` and `df_complete`;
  `coverage_claim_allowed = TRUE`.

- `descriptive_classic_rubin`: `stack_direct` with the classic degrees
  of freedom; `coverage_claim_allowed = FALSE`.

- `reference_classic_rubin` and `reference_barnard_rubin`: `per_pv`;
  always `coverage_claim_allowed = FALSE`.

- `psis_classic_rubin` and `psis_barnard_rubin`: `stack_psis`; always
  `coverage_claim_allowed = FALSE`.

The labels do not depend on the status: a fit with status `"warning"`
keeps them.

When some or all intervals are descriptive,
[`print()`](https://rdrr.io/r/base/print.html) of a fit or a comparison
adds a line beginning "interval note:".

## Source labels

`target_source` names where the estimates come from:

- `"external_brr_fay_rubin"` (`stack_direct`): the target that
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  builds from the BRR-Fay replicate weights and Rubin's rules.
  [`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
  returns it, and `target_hash` is its checksum.

- `"per_pv_rubin_draws"` (`per_pv`): Rubin's rules applied to the draws
  of the per-plausible-value fits; the result is the fit's `target`.

- `"stack_psis_rubin_pooling"` (`stack_psis`): Rubin's rules applied to
  the weighted results. This is only a label: a `stack_psis` fit has no
  target object, and
  [`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
  returns `NULL`.

For `per_pv` and `stack_psis`, `pooling_source` repeats this label, and
both `pooling_hash` and `target_hash` hold a checksum of the combined
result.

A `stack_psis` fit records how its weights were obtained, in the
estimate table and in `get_diagnostics(fit)$psis`:

- `psis_source`: `"supplied_psis_weights"` (you passed `psis_weights`
  and `pareto_k`), `"injected_psis_function"` (your `psis_function`
  computed the weights and the Pareto k-hat values from `log_ratios`) or
  `"self_normalized_log_ratios"` (you passed `log_ratios` and
  `pareto_k`, and pvstackr turned the ratios into normalized weights
  without smoothing).

- `pareto_k_source`: `"supplied"` or `"injected_function_output"`.

- `weight_method`: `"caller_declared_external_psis"` when you named the
  program that produced the weights (`psis_producer` and
  `psis_producer_version`), `"unspecified_external"` when you did not,
  and `"self_normalized_raw_importance"` for weights made from
  `log_ratios`. Only the first can give estimates; pvstackr records the
  program that you name but does not check it.

`get_diagnostics(fit)$psis` also records how concentrated the normalized
weights `w` of each plausible value are: `weight_ess_iid`
(`1 / sum(w^2)`, a Kish-type effective sample size),
`weight_ess_fraction` (that value divided by the number of draws) and
`max_normalized_weight`. These are not MCMC effective sample sizes and
do not change the Pareto k-hat check. `weight_diagnostic_authority` is
`"retained_weights_recomputed"` when the weights are kept
(`return_draws = TRUE` on a fit with estimates), so that pvstackr can
recompute the three values, and `"owned_stamp_bounded_projection"` when
the weights were removed: the values are then covered by the fit's
checksum but cannot be recomputed.

## Method comparisons

[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
returns a list of class `pvstackr_method_comparison` and describes its
tables; read it with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
Agreement between methods is descriptive.

## Reading results

Read a fit with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md);
they stop with an error if the fit was changed after it was created (see
"Technical details"). A `stack_psis` fit saved by an earlier version of
pvstackr that does not pass the current checks must go through
[`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md),
which returns an object for inspection only.

## Technical details

`validation` holds a SHA-256 checksum of the fit, which the four reading
functions, [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) recompute; it detects
changes but is not a digital signature. A model kept with
`keep_backend_fit = TRUE` is not covered by the checksum, so changes
inside that model are not detected. `schema_version` is the format
version of the object and `provenance` a record of how it was built.
