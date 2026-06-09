# Assemble a Rubin/BRR-Fay Fixed-Effect Target

`pv_brr_target()` assembles the external, design-based fixed-effect
covariance target that `stack_direct` calibrates against. For each
plausible value it forms a BRR-Fay sandwich covariance from the supplied
replicate weights, then Rubin-combines the per-plausible-value estimates
and covariances across plausible values into a single fixed-effect mean,
total covariance, standard errors, and degrees of freedom. The target
engine is a dependency-free weighted least-squares fit; in this package
stage the target is fixed-effect-only.

## Usage

``` r
pv_brr_target(
  data,
  formula,
  pv_cols = NULL,
  weight_col,
  rep_weight_cols = NULL,
  fay_k = 0.5,
  pv_prefix = "PV",
  pv_suffix = "",
  rep_weight_prefix = "W_FSTURWT",
  expected_M = NULL,
  expected_R = NULL,
  id_cols = NULL,
  conf_level = 0.95,
  allow_m1 = FALSE,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL,
  engine = "lm",
  verbose = FALSE
)

# S3 method for class 'pvstackr_brr_target'
print(x, ...)
```

## Arguments

- data:

  Data frame containing the outcome, predictors, plausible-value
  columns, weight column, and replicate-weight columns.

- formula:

  Two-sided formula with `OUTCOME` as the left-hand-side placeholder,
  for example `OUTCOME ~ x + female`. The right-hand side defines the
  fixed-effect design; group terms such as `(1 | school)` are not
  accepted in this fixed-effect-only stage.

- pv_cols:

  Character vector of plausible-value column names. If `NULL`, columns
  are detected with
  [`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)
  using `pv_prefix`, `pv_suffix`, and `expected_M`.

- weight_col:

  Character scalar naming the full-sample weight column. Required; there
  is no default.

- rep_weight_cols:

  Character vector of BRR replicate-weight column names (at least two,
  distinct from `weight_col`). If `NULL`, columns are detected with
  [`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md)
  using `rep_weight_prefix` and `expected_R`.

- fay_k:

  Numeric scalar Fay coefficient \\k\\ entering \\a_d = 1/(R(1-k)^2)\\.
  Must satisfy `0 <= fay_k < 1`. Default `0.5` (the PISA convention).

- pv_prefix, pv_suffix:

  Character scalars giving the prefix and suffix used when detecting
  plausible-value columns. Modern PISA files use subject-suffixed
  plausible values such as `PV1MATH`, so set, for example,
  `pv_suffix = "MATH"` rather than relying on the bare `pv_suffix = ""`
  default with `pv_prefix = "PV"`. Ignored when `pv_cols` is supplied.

- rep_weight_prefix:

  Character scalar prefix used when detecting replicate weights. Default
  `"W_FSTURWT"`. Ignored when `rep_weight_cols` is supplied.

- expected_M:

  Optional integer expected plausible-value count. If supplied, it is
  enforced against detected or supplied `pv_cols`. Default `NULL`.

- expected_R:

  Optional integer expected replicate-weight count. If supplied, it is
  enforced against detected or supplied `rep_weight_cols`. Default
  `NULL`.

- id_cols:

  Optional character vector of row-identifier columns. If supplied, they
  must jointly identify unique rows. Default `NULL`.

- conf_level:

  Numeric scalar interval level in `(0, 1)`, passed through to Rubin
  pooling. Default `0.95`.

- allow_m1:

  Logical scalar; whether to allow a single plausible value (`M = 1`).
  Default `FALSE`.

- df_method:

  Degrees-of-freedom rule for Rubin pooling, one of `"classic"` or
  `"barnard_rubin"`. The default `"classic"` is descriptive only;
  `"barnard_rubin"` requires `df_complete`.

- df_complete:

  Complete-data degrees of freedom used when
  `df_method = "barnard_rubin"`; a positive numeric scalar or a per-term
  named numeric vector. Ignored under `"classic"`. Default `NULL`.

- engine:

  Character scalar target engine. Only `"lm"` (dependency-free weighted
  least squares) is implemented in this package stage. Default `"lm"`.

- verbose:

  Logical scalar; whether to emit per-plausible-value progress messages.
  Default `FALSE`.

- x:

  A `pvstackr_brr_target` object from `pv_brr_target()`.

- ...:

  Ignored.

## Value

A `pvstackr_brr_target` object: a list carrying the external
fixed-effect target and its provenance, with class
`pvstackr_brr_target`. Reportable fields include:

- `beta_bar`, `U_bar`, `B`, `T_MI`:

  The Rubin mean \\\bar\beta\\, within-imputation covariance \\\bar U\\,
  between-imputation covariance \\B\\, and total (target) covariance
  \\T\_{\mathrm{MI}}\\ over the fixed-effect block.

- `se`:

  Target standard errors, \\\sqrt{\mathrm{diag}(T\_{\mathrm{MI}})}\\.

- `df`, `df_classic`, `df_method`, `df_complete`:

  Active Rubin degrees of freedom, the classic-df reference, the df rule
  (`"classic"` or `"barnard_rubin"`), and the complete-data df when
  Barnard-Rubin is used.

- `interval_role`, `coverage_claim_allowed`:

  Interval-policy fields for downstream fits (see the Interval metadata
  section and
  [pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)).

- `lambda`, `fmi`, `riv`:

  Per-term Rubin missing-information quantities: `lambda` (proportion of
  total variance attributable to the between-imputation component),
  `fmi` (fraction of missing information), and `riv` (relative increase
  in variance).

- `fe_names`, `M`, `R`, `fay_k`:

  Fixed-effect names, number of plausible values \\M\\, replicate count
  \\R\\, and the Fay coefficient \\k\\.

- `fay_variance_multiplier`:

  The BRR-Fay multiplier \\a_d = 1/(R(1-k)^2)\\.

- `pv_cols`, `weight_col`, `rep_weight_cols`, `id_cols`:

  Resolved column names defining the plausible values, full-sample
  weight, replicate weights, and row identifiers.

- `target_source`, `target_hash`, `design_hash`:

  Provenance: `target_source = "external_brr_fay_rubin"`, and stable
  content hashes of the target and its design manifest.

## Details

The construction follows the canonical equation set (per-plausible-value
BRR-Fay within-covariance, then Rubin combining).

### Per-plausible-value BRR-Fay sandwich

For plausible value \\m\\, the design-based within-covariance is the
replicate sandwich \$\$\hat U_m = a_d \sum\_{r=1}^{R}
(\hat\beta_m^{(r)} - \hat\beta_m)(\hat\beta_m^{(r)} - \hat\beta_m)^\top,
\qquad a_d = \frac{1}{R\\(1 - k)^2},\$\$ where \\\hat\beta_m\\ is the
full-sample weighted fixed-effect estimate, \\\hat\beta_m^{(r)}\\ is the
estimate under replicate weight \\r\\, \\R\\ is the replicate count, and
\\k\\ is the Fay coefficient (`fay_k`). The multiplier \\a_d\\ is
returned as `fay_variance_multiplier`.

### Rubin combining across plausible values

The per-plausible-value estimates and covariances are combined into the
Rubin mean \\\bar\beta\\, the within-imputation covariance \\\bar U\\,
and the between-imputation covariance \\B\\, \$\$\bar\beta =
\frac{1}{M}\sum\_{m=1}^{M}\hat\beta_m, \qquad \bar U =
\frac{1}{M}\sum\_{m=1}^{M}\hat U_m, \qquad B =
\frac{1}{M-1}\sum\_{m=1}^{M} (\hat\beta_m - \bar\beta)(\hat\beta_m -
\bar\beta)^\top,\$\$ giving the Rubin total covariance — the external
target — \$\$T\_{\mathrm{MI}} = \bar U + \left(1 + \tfrac{1}{M}\right)
B.\$\$ Standard errors are \\\sqrt{\mathrm{diag}(T\_{\mathrm{MI}})}\\.
The fraction of missing information \\\gamma_k = (1 +
1/M)\\B\_{kk}/T\_{\mathrm{MI},kk}\\ and the related Rubin quantities are
returned as `fmi`, `riv`, and `lambda`.

### Degrees of freedom

The default `df_method = "classic"` uses the classic Rubin imputation
degrees of freedom and is descriptive only.
`df_method = "barnard_rubin"` uses the Barnard-Rubin small-sample
degrees of freedom and requires `df_complete` (the complete-data degrees
of freedom); `df_classic` is retained alongside the active `df` for
reference.

### Computation

The target requires on the order of \\M\\(R+1)\\ weighted fixed-effect
fits: one full-sample fit plus \\R\\ replicate fits for each of the
\\M\\ plausible values. This is a description of the computational work,
not a benchmarked performance claim.

## Interval metadata

The target fixes the interval-policy fields that downstream fits carry
on reportable rows (consistent with
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)):

- `interval_role`:

  `"descriptive_classic_rubin"` under the default
  `df_method = "classic"`, or `"coverage_barnard_rubin"` under
  `df_method = "barnard_rubin"` with `df_complete`.

- `coverage_claim_allowed`:

  `FALSE` for the classic-df target; `TRUE` only under
  `"coverage_barnard_rubin"`. A coverage claim is actually realized
  downstream by a `stack_direct` fit backed by this external target;
  this object only declares the policy.

- `df_method`:

  `"classic"` or `"barnard_rubin"`. Selecting `"barnard_rubin"` does not
  by itself make a row coverage-claimable absent the `stack_direct`
  target provenance.

- `df_complete`:

  Complete-data degrees of freedom recorded when
  `df_method = "barnard_rubin"`; otherwise unset.

This is a design-based external target. Coverage is reserved for
`stack_direct` rows under Barnard-Rubin df (see
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md));
no speed or efficiency claim attaches to assembling the target.

## References

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.*
Wiley.

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948-955.

Judkins, D. R. (1990). Fay's method for variance estimation. *Journal of
Official Statistics*, 6(3), 223-239.

## See also

Build the design with
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md);
detect columns with
[`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)
and
[`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md).
Calibrate against this target with
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
or via the dispatcher
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md).

## Examples

``` r
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)
design <- pv_design(
  pisa_tiny, formula = OUTCOME ~ x + female,
  pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
)
target <- pv_brr_target(
  pisa_tiny, OUTCOME ~ x + female,
  pv_cols = design$pv_cols, weight_col = design$weight_col,
  rep_weight_cols = design$rep_weight_cols, fay_k = design$fay_k,
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
