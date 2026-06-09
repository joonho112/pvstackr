
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pvstackr <img src="man/figures/logo.png" align="right" height="139" alt="pvstackr logo" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/joonho112/pvstackr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/joonho112/pvstackr/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/joonho112/pvstackr/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/joonho112/pvstackr/actions/workflows/pkgdown.yaml)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

`pvstackr` is an R package for Bayesian-backend fixed-effect calibration
and reporting with plausible-value survey data.

## The problem

Large educational surveys such as PISA report each latent score not as a
single number but as **M plausible values** — M imputations that
together carry the measurement uncertainty in the construct. The honest
analysis fits the model M times and combines the results with Rubin’s
rules so that *between-imputation* variance is propagated into every
standard error.

The tempting shortcut is to **fit `PV1` only** and report that single
fit. It is fast and it looks fine, but it silently discards the
between-imputation variance. The point estimate drifts and, more
dangerously, the standard errors are too small, so confidence intervals
**under-cover**: a nominal 95% interval traps the truth far less than
95% of the time.

`pvstackr` produces an *honest* fixed-effect fit from **one stacked
fit** calibrated to an external **Rubin/BRR-Fay fixed-effect target**,
and records, **per row**, whether an interval is coverage-claimable. You
get a single set of reportable fixed-effect estimates and standard
errors that carry the design-based combining variance, together with
diagnostics that say plainly when an interval may be read as a coverage
statement and when it may not.

`pvstackr` is a focused, self-contained method package — not a paper-replication
repository. A companion methods paper is in preparation.

## Installation

You can install the development version of `pvstackr` from
[GitHub](https://github.com/joonho112/pvstackr) with:

``` r
# install.packages("pak")
pak::pak("joonho112/pvstackr")

# CRAN release coming soon:
# install.packages("pvstackr")
```

Optional live Bayesian backends (e.g. `cmdstanr`/`brms`) and their
toolchains are **not** required for installation, loading, examples, or
the default tests.

## The method in three steps

A reportable `stack_direct` analysis is three calls:

``` r
# 1. Declare the plausible-value design (PV columns, BRR replicate weights, Fay k)
design <- pv_design(data, formula, pv_suffix = "READ", ...)

# 2. Assemble the external Rubin/BRR-Fay fixed-effect target
target <- pv_brr_target(data, formula, pv_cols = design$pv_cols, ...)

# 3. Fit one stacked model and calibrate its draws to that target
fit <- pv_fit(data, formula, target = target, method = "stack_direct", ...)
```

`pv_design()` resolves the survey design, `pv_brr_target()` builds the
design-based combining target, and `pv_fit(method = "stack_direct")`
produces the calibrated, reportable fixed-effect fit. The Five-Minute
Synthetic Example below runs the whole shape end to end on a bundled
synthetic fixture.

## Method Names

The public method IDs are:

| Method | Role |
|----|----|
| `stack_direct` | Default calibrated-reporting path: a design-based external Rubin/BRR-Fay fixed-effect target for reported estimates and SEs, plus one stacked backend fit used to calibrate draws and diagnose agreement. |
| `stack_psis` | Diagnostic/reference method: one stacked draw source plus supplied/precomputed or injected PSIS weights, Pareto-k diagnostics, and model-based Rubin pooling of PSIS-weighted fixed-effect summaries. |
| `per_pv` | Per-PV Bayesian/backend reference: one fit or posterior-draw source per plausible value, with fixed-effect centers and model-based within-PV covariance combined by Rubin pooling. |

In the current package stage, all three method IDs are implemented
through injected or precomputed light paths. `stack_psis` currently
expects supplied or injected PSIS weights and Pareto-k diagnostics
rather than depending on a live PSIS backend.

## Scope

The `stack_direct` path uses a one-stacked-fit architecture for
Bayesian/backend fitting while preserving a design-based external
Rubin/BRR-Fay fixed-effect target. `per_pv` and `stack_psis` are
reference/diagnostic paths: they use model-based within-PV covariance in
Rubin pooling and do not construct or calibrate to the external BRR-Fay
target.

Interval metadata follows that method boundary. Coverage-claimable
intervals are reserved for `stack_direct` rows backed by an external
BRR-Fay target with Barnard-Rubin df and explicit `df_complete`;
`per_pv` and `stack_psis` intervals remain descriptive/reference
intervals even when they use Barnard-Rubin df.

The package should not be read as claiming unqualified “same inference
with one MCMC.” In v0.1, calibrated reporting is scoped to fixed effects
and is not a full posterior for all model parameters.

## Current `stack_direct` Boundary

`stack_direct` v0.1 requires an external `pvstackr_brr_target` from
`pv_brr_target()`. The public wrapper checks formula RHS equality,
fixed-effect name alignment, target provenance, and target policy before
fitting. Group terms such as `(1 | school)` and `(1 || school)` are
rejected in this version because the current BRR-Fay target engine is
fixed-effect-only.

Reportable `pvstackr_fit` objects contain the design, target, stacked
fit, CCC calibration object, fixed-effect estimate table, and calibrated
fixed-effect draws when requested by `control$return_draws`. Variance
components and sampler diagnostic columns may be retained inside nested
backend or stack-fit components, but they are not calibrated reportable
estimates.

Under the default `control$center = "target"`, the reported fixed-effect
estimate and standard error are the external Rubin/BRR-Fay target
values: `estimate` is the target coefficient after CCC centering, and
`se` is `sqrt(diag(target$T_MI))` with the target degrees of freedom and
interval metadata. The stacked fit supplies the raw draw cloud for
calibration, retained calibrated fixed-effect draws when requested, and
center-separation agreement diagnostics. It is not a replacement source
for the headline Rubin/BRR-Fay fixed-effect numbers.

Reportable `stack_direct` fits require `control$center = "target"`. The
posterior-centered CCC convention is reserved for diagnostic/exploratory
checks because it retains the raw stacked posterior center while using
the external target covariance.

Yellow CCC center-separation diagnostics and explicit priors become
top-level warning-status fits. Red center separation blocks reportable
estimates. Routine dropped draw columns such as `lp__` remain nested in
`stack_fit$warnings`.

For the formal object contract, see `?pvstackr_object_contracts`.

## Five-Minute Synthetic Example

The package ships a tiny, package-owned synthetic PISA-shaped fixture
for API smoke tests and examples. It contains no real PISA records or
OECD/PISA distributed files, and is not suitable for real inference,
coverage claims, or performance benchmarking.

``` r
library(pvstackr)

pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

list(
  pv_cols = detect_pisa_pv_columns(
    pisa_tiny,
    suffix = "READ",
    expected_M = 2L
  ),
  rep_weight_cols = detect_pisa_brr_replicate_weights(
    pisa_tiny,
    expected_R = 4L
  )
)
#> $pv_cols
#> [1] "PV1READ" "PV2READ"
#> 
#> $rep_weight_cols
#> [1] "W_FSTURWT1" "W_FSTURWT2" "W_FSTURWT3" "W_FSTURWT4"
```

``` r
design <- pv_design(
  data = pisa_tiny,
  formula = OUTCOME ~ x + female,
  pv_suffix = "READ",
  expected_M = 2L,
  expected_R = 4L,
  id_cols = "CNTSTUID"
)

design
#> pvstackr design
#>   rows: 12
#>   formula: OUTCOME ~ x + female
#>   plausible values: 2
#>   final weight: W_FSTUWT
#>   replicate weights: 4
#>   fay_k: 0.5
#>   design hash: 4776c0f1

target <- pv_brr_target(
  data = pisa_tiny,
  formula = OUTCOME ~ x + female,
  pv_cols = design$pv_cols,
  weight_col = design$weight_col,
  rep_weight_cols = design$rep_weight_cols,
  fay_k = design$fay_k,
  id_cols = design$id_cols
)

target
#> pvstackr BRR-Fay target
#>   fixed effects: 3
#>   plausible values: 2
#>   replicate weights: 4
#>   fay_k: 0.5
#>   df method: classic
#>   interval role: descriptive_classic_rubin
#>   source: external_brr_fay_rubin
```

The cached example fit was generated by the package's fixture-build script through
the lightweight injected backend. It lets the README demonstrate the
public object surface without running heavy MCMC.

``` r
cached <- readRDS(
  system.file(
    "extdata",
    "examples",
    "pisa_tiny_stack_direct.rds",
    package = "pvstackr"
  )
)

fit <- cached$fit
fit
#> pvstackr fit
#>   method: stack_direct
#>   status: ok
#>   fixed effects: 3
#>   target: external_brr_fay_rubin
#>   draws: not retained
#>   diagnostics: preflight, stack_fit, stack_fit_warnings, ccc
#>   interval note: intervals are descriptive rather than coverage-claimable.

get_estimates(fit)[
  ,
  c("term", "estimate", "se", "df", "conf_low", "conf_high")
]
#>          term   estimate        se       df  conf_low conf_high
#> 1 b_Intercept 457.894088 1.2873118 1.021194 442.31804 473.47013
#> 2         b_x  46.883361 0.3717929 1.402308  44.41457  49.35215
#> 3    b_female   2.143702 3.5550309 1.013730 -41.60687  45.89428

get_target(fit)$target_hash
#> [1] "4a4d40f8"
names(get_diagnostics(fit))
#> [1] "preflight"          "stack_fit"          "stack_fit_warnings"
#> [4] "ccc"
get_draws(fit)
#> NULL
```

## Conceptual Live Workflow

For a real analysis, regenerate the design, target, and fit from the
real data. The current `stack_direct` boundary is explicit: build the
external Rubin/BRR-Fay fixed-effect target first, then fit and calibrate
the stacked model to that target. For reportable output,
`center = "target"` is required and the fixed-effect estimate table is
target-based; the stacked model contributes calibrated draws and
agreement diagnostics around that external target. Backend wiring is
user- or adapter-specific in this package stage, so the
real-data/live-backend sketch below is not run by the README.

See `vignette("a5-real-pisa-guidance", package = "pvstackr")` for
real-data loading, licensing, memory, and runtime guidance. Real PISA
work can be large: stacked preparation creates about `N * M` rows, and
BRR-Fay target construction runs about `M * (R + 1)` fixed-effect fits
per model before any live Bayesian backend runs.

Modern PISA plausible values are often subject-suffixed. For columns
such as `PV1MATH`, set `pv_suffix = "MATH"`; leave `pv_suffix = ""` only
when the analytic extract intentionally uses bare plausible-value names
such as `PV1`, `PV2`, and so on.

``` r
design <- pv_design(
  data = pisa_country,
  formula = OUTCOME ~ escs + female,
  pv_suffix = "READ",
  expected_M = 10L,
  expected_R = 80L,
  id_cols = "CNTSTUID"
)

target <- pv_brr_target(
  data = design$data,
  formula = design$formula,
  pv_cols = design$pv_cols,
  weight_col = design$weight_col,
  rep_weight_cols = design$rep_weight_cols,
  fay_k = design$fay_k,
  id_cols = design$id_cols
)

fit <- pv_fit(
  data = design$data,
  formula = design$formula,
  target = target,
  method = "stack_direct",
  control = pv_control(
    method = "stack_direct",
    backend = "injected",
    iter = 2000L,
    warmup = 1000L,
    chains = 4L,
    seed = 20260607
  ),
  fit_function = fit_stacked_model,
  draws_function = extract_fixed_effect_draws
)

summary(fit)
```

## Vignettes

The documentation is organized into two tracks. The **Applied track** is
workflow-first with light math; the **Method track** derives the
estimators and states the coverage-claimability contract. Every vignette
runs on the bundled synthetic fixture — real PISA numbers are
*described*, never recomputed in-package.

### Applied track

| Article | What it covers |
|----|----|
| [A1 · Getting started](https://joonho112.github.io/pvstackr/articles/a1-getting-started.html) | The fastest path: install, the single `pv_fit(method = "stack_direct")` shape, load the cached fixture fit, and read it three ways. |
| [A2 · The end-to-end workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html) | The full pipeline on synthetic data — `pv_design()` → `pv_brr_target()` → `pv_fit()` — plus the conceptual live-backend sketch. |
| [A3 · Reading results & what to report](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html) | The accessor family, the estimate table, the interval-metadata columns, fraction of missing information, and a reporting checklist. |
| [A4 · Comparing methods](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.html) | `pv_compare_methods()` aligns the three paths; how to read agreement diagnostics and the two cautions; which method to use. |
| [A5 · Real PISA data guidance](https://joonho112.github.io/pvstackr/articles/a5-real-pisa-guidance.html) | Connecting to real PISA without bundling licensed files: licensing & non-affiliation, the design declaration, memory/runtime, reproducibility. |

### Method track

| Article | What it covers |
|----|----|
| [M1 · Foundations and notation](https://joonho112.github.io/pvstackr/articles/m1-foundations-and-notation.html) | Fixes the notation the track reuses: the PV survey setting, M/R/Fay k, the within/between model, and the fixed-effect estimand. |
| [M2 · The BRR–Fay target](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.html) | Derives the external target — Rubin combining, the BRR–Fay sandwich, Barnard–Rubin df, FMI — and the design-variance coverage result. |
| [M3 · The stacked fractional bridge](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.html) | What “stacked fractional” means, the stacked-MLE identity, why one stacked fit recovers the Rubin average, and where the identity stops. |
| [M4 · CCC calibration](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.html) | The Cholesky Calibration Correction map, the centering conventions, the deterministic algorithm, and the center-separation diagnostics. |
| [M5 · Methods, PSIS & coverage](https://joonho112.github.io/pvstackr/articles/m5-methods-and-coverage.html) | The precise contrast of the three methods, PSIS weights and Pareto-$\hat k$ thresholds, and why only `stack_direct` is coverage-claimable. |

## Citation

To cite `pvstackr` in publications, please cite the package:

``` r
citation("pvstackr")
```

JoonHo Lee · <jlee296@ua.edu> · [ORCID
0009-0006-4019-8703](https://orcid.org/0009-0006-4019-8703) · Assistant
Professor, The University of Alabama.

The method paper describing the BRR-Fay target, the stacked fractional
bridge, and CCC calibration is forthcoming; a method citation will be
added here once it is finalized.

## Related work

`pvstackr` is built to interoperate with the established survey and
multiple-imputation ecosystem in R. Its runtime BRR-Fay target engine is
dependency-free base-R weighted least squares; optional development
tests compare that engine to
[survey](https://CRAN.R-project.org/package=survey)’s replicate-weight
machinery. The package follows Rubin-combining conventions familiar from
[mitools](https://CRAN.R-project.org/package=mitools) and
[mice](https://CRAN.R-project.org/package=mice). For general PISA-style
plausible-value analysis, see
[intsvy](https://CRAN.R-project.org/package=intsvy) and
[BIFIEsurvey](https://CRAN.R-project.org/package=BIFIEsurvey);
`pvstackr` differs by calibrating a single stacked fit to an external
design-based fixed-effect target rather than pooling per-PV model fits.

## Getting help

Found a bug, or have a question about the workflow? Please open an issue
at <https://github.com/joonho112/pvstackr/issues>.

## Status

This package is in early development. Public APIs are being added step
by step with focused tests. Real PISA data are not bundled. Optional
live Bayesian backends and their toolchains are not required for package
installation, loading, README rendering, or default tests. `pvstackr` is
an independent research package and is not affiliated with or endorsed
by the OECD or the PISA programme.

## License

MIT © JoonHo Lee. See [LICENSE.md](LICENSE.md) for details.
