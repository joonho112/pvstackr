# Package index

## Package overview

Orientation plus the formal object and interval contracts.

- [`pvstackr`](https://joonho112.github.io/pvstackr/reference/pvstackr-package.md)
  [`pvstackr-package`](https://joonho112.github.io/pvstackr/reference/pvstackr-package.md)
  : pvstackr: Stacked-Fit Calibration to Rubin/BRR-Fay Fixed-Effect
  Targets for Plausible Values
- [`pvstackr_object_contracts`](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
  : pvstackr Object Contracts

## Declare the design

Describe a PISA-style PV design and detect PV / BRR replicate columns.

- [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
  [`print(`*`<pvstackr_design>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
  : Declare a PISA-Style Plausible-Value Design
- [`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)
  : Detect PISA-Style Plausible-Value Columns
- [`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md)
  : Detect PISA-Style BRR Replicate-Weight Columns

## Assemble the target

Build the external Rubin / BRR-Fay fixed-effect target.

- [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  [`print(`*`<pvstackr_brr_target>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  : Assemble a Rubin/BRR-Fay Fixed-Effect Target
- [`pv_revalidate_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_revalidate_brr_target.md)
  : Revalidate a Legacy BRR-Fay Target Against Its Original Inputs

## Fit a method

The dispatcher, the three method engines, and fitting controls.

- [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
  : Fit a pvstackr Method
- [`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
  : Fit the Direct Stacked Plausible-Value Model
- [`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md)
  : Fit the Per-PV Bayesian/Backend Reference Method
- [`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)
  : Fit the PSIS-Reweighted Stacked Method
- [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  [`print(`*`<pvstackr_control>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  : Construct pvstackr Fitting Controls

## Attach a backend

The bundled brms adapter, exported so an injected adapter can reuse it
or replace one piece of it.

- [`pv_backend_brms_fit_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
  [`pv_backend_brms_draws_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
  [`pv_backend_brms_sampler_diagnostics()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
  : Bundled brms adapter for the stacked fit

## Compare methods

Aligned fixed-effect comparison across methods with agreement
diagnostics.

- [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
  : Compare pvstackr Method Fits

## Read results

Stable accessors for the reportable estimate table, target, draws, and
diagnostics.

- [`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
  : Access pvstackr Estimates
- [`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
  : Access pvstackr Targets
- [`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
  : Access pvstackr Draws
- [`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
  : Access pvstackr Diagnostics

## Display methods

Compact console displays for all object types; summary methods for fit
and comparison objects.

- [`print(`*`<pvstackr_fit>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_fit_summary.md)
  [`summary(`*`<pvstackr_fit>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_fit_summary.md)
  [`print(`*`<summary.pvstackr_fit>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_fit_summary.md)
  : Display Methods for pvstackr Fits
- [`print(`*`<pvstackr_method_comparison>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_method_comparison_summary.md)
  [`summary(`*`<pvstackr_method_comparison>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_method_comparison_summary.md)
  [`print(`*`<summary.pvstackr_method_comparison>`*`)`](https://joonho112.github.io/pvstackr/reference/pvstackr_method_comparison_summary.md)
  : Display Methods for pvstackr Method Comparisons

## Migration

Reading objects written by an earlier release.

- [`pv_migrate_legacy_psis_fit()`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  [`print(`*`<pvstackr_legacy_psis_inspection>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  [`summary(`*`<pvstackr_legacy_psis_inspection>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  [`print(`*`<summary.pvstackr_legacy_psis_inspection>`*`)`](https://joonho112.github.io/pvstackr/reference/pv_migrate_legacy_psis_fit.md)
  : Convert a Legacy PSIS Fit to a Safe Inspection Object
