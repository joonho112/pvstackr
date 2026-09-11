# Package index

## Package overview

What pvstackr computes, and the objects that its functions return.

- [`pvstackr`](https://joonho112.github.io/pvstackr/reference/pvstackr-package.md)
  [`pvstackr-package`](https://joonho112.github.io/pvstackr/reference/pvstackr-package.md)
  : pvstackr: Bayesian Plausible-Value Analysis with One Calibrated
  Stacked Fit
- [`pvstackr_object_contracts`](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
  : Objects returned by pvstackr

## Check the columns

Find and check the plausible-value, final-weight and replicate-weight
columns of a PISA-style data set.

- [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
  [`print(`*`<pvstackr_design>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
  : Declare the plausible-value and weight columns
- [`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)
  : Find PISA-style plausible-value columns
- [`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md)
  : Find PISA-style replicate-weight columns

## Compute the target

Compute the design-based target for the fixed effects from the BRR-Fay
replicate weights and Rubin’s rules, and check or rebuild a saved
target.

- [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  [`print(`*`<pvstackr_brr_target>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  : Compute the design-based target for the fixed effects
- [`pv_revalidate_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_revalidate_brr_target.md)
  : Check a BRR-Fay target or rebuild it in the current format

## Fit a model

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
the functions for the three methods, and their settings.

- [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
  : Fit a model to plausible-value data
- [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  : Fit a stacked model and calibrate it to a BRR-Fay target
- [`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)
  : Combine one fit per plausible value with Rubin's rules
- [`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)
  : Reweight stacked draws toward each plausible value
- [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  [`print(`*`<pvstackr_control>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  : Settings for the fitting functions

## Fitting engine

The fitting functions of the bundled brms engine, which you can reuse
when you supply fitting functions of your own.

- [`pv_backend_brms_fit_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
  [`pv_backend_brms_draws_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
  [`pv_backend_brms_sampler_diagnostics()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
  : Bundled brms engine for the stacked fit

## Compare methods

Put fits made with different methods side by side, with agreement
diagnostics.

- [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
  : Compare fits made with different methods

## Read results

Read the estimate table, the target, the draws and the diagnostics of a
fit.

- [`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
  : Get the estimate table of a fit
- [`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
  : Get the target of a fit
- [`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
  : Get the calibrated draws of a fit
- [`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
  : Get the diagnostics of a fit or comparison

## Print and summary

Print methods for all pvstackr objects, and summary methods for fits and
comparisons.

- [`print(`*`<pvstackr_fit>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_fit_summary.md)
  [`summary(`*`<pvstackr_fit>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_fit_summary.md)
  [`print(`*`<summary.pvstackr_fit>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_fit_summary.md)
  : Print and summarize a fit
- [`print(`*`<pvstackr_method_comparison>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_method_comparison_summary.md)
  [`summary(`*`<pvstackr_method_comparison>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_method_comparison_summary.md)
  [`print(`*`<summary.pvstackr_method_comparison>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_method_comparison_summary.md)
  : Print and summarize a method comparison

## Earlier versions

Read `stack_psis` fits saved by an earlier version of pvstackr.

- [`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  [`print(`*`<pvstackr_legacy_psis_inspection>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  [`summary(`*`<pvstackr_legacy_psis_inspection>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  [`print(`*`<summary.pvstackr_legacy_psis_inspection>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  : Inspect a stack_psis fit saved by an earlier version of pvstackr
