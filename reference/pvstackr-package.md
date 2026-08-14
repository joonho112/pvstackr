# pvstackr: Stacked-Fit Calibration to Rubin/BRR-Fay Fixed-Effect Targets for Plausible Values

`pvstackr` currently implements the `stack_direct` path: one stacked
plausible-value fit plus fixed-effect CCC calibration to a design-based
external Rubin/BRR-Fay target. It also implements the `per_pv` reference
path: one backend fit or posterior-draw source per plausible value plus
model-based Rubin pooling, and the `stack_psis` diagnostic/reference
path: one stacked draw source plus caller-declared external PSIS
weights, model-based Rubin pooling, and Pareto-k gating.

## Details

### Typical workflow

A reportable analysis moves through the package API in five stages:
declare the plausible-value design with
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md);
assemble the external Rubin/BRR-Fay fixed-effect target with
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md);
fit a method with
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
(default `method = "stack_direct"`); optionally line the methods up side
by side with
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md);
and read the results with the accessor family
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).

### Where to learn more

The documentation is organized as two vignette tracks. The Applied track
is workflow-first and starts at
[`vignette("a1-getting-started", package = "pvstackr")`](https://joonho112.github.io/pvstackr/articles/a1-getting-started.md).
The Method track is the rigorous treatment of the notation, target,
stacked bridge, and calibration, and starts at
[`vignette("m1-foundations-and-notation", package = "pvstackr")`](https://joonho112.github.io/pvstackr/articles/m1-foundations-and-notation.md).

### Scope and honesty

Reportable output is fixed-effect-only (variance components are fit but
not calibrated); among the three methods only `stack_direct`, whose
intervals rest on the external Rubin/BRR-Fay target, is
coverage-claimable, while `per_pv` and `stack_psis` are
descriptive/reference paths.

## See also

Design:
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md),
[`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md),
[`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md).
Target:
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
Fit:
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md),
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).
Compare:
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md).
Read:
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md),
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
Contracts:
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]

Authors:

- JoonHo Lee <jlee296@ua.edu>
  ([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]
