# Compute the design-based target for the fixed effects

`pv_brr_target()` estimates the fixed effects of a survey-weighted
linear regression from plausible-value data, with their covariance
matrix, standard errors and degrees of freedom. It fits the regression
to each plausible value with the final weight, estimates the sampling
covariance of the coefficients from the BRR-Fay replicate weights
(balanced repeated replication with Fay's method), and combines the
plausible values with Rubin's rules. The result is the target of a
`stack_direct` fit
([`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)),
which reports the target's estimates, standard errors and degrees of
freedom and calibrates its fixed-effect draws to the target's mean and
covariance.

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

  A data frame with the plausible-value columns, the final and replicate
  weights, and the variables in `formula`. A `stack_direct` fit with
  this target must use the same data
  ([`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)).

- formula:

  A two-sided formula with the placeholder `OUTCOME` on the left-hand
  side, for example `OUTCOME ~ x + female`; `OUTCOME` stands for the
  plausible values. The right-hand side gives the fixed effects, as in
  [`lm()`](https://rdrr.io/r/stats/lm.html), and its variables must be
  columns of `data` without missing values. A random-effect term such as
  `(1 | school)` stops the function with an error, and so does a
  [`weights()`](https://rdrr.io/r/stats/weights.html) term: name the
  weights in `weight_col` and `rep_weight_cols`.

- pv_cols:

  A character vector with the names of the plausible-value columns, or
  `NULL` (default) to find them with
  [`detect_pisa_pv_columns()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_pv_columns.md)
  from `pv_prefix` and `pv_suffix`. At least two are needed unless
  `allow_m1 = TRUE`, and their values must be finite numbers.

- weight_col:

  The name of the final weight column, whose values must be positive.
  There is no default: give, for example, `design$weight_col` from
  [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md),
  or `"W_FSTUWT"` for PISA data.

- rep_weight_cols:

  A character vector with the names of at least two replicate-weight
  columns other than `weight_col`, or `NULL` (default) to find them with
  [`detect_pisa_brr_replicate_weights()`](https://joonho112.github.io/pvstackr/reference/detect_pisa_brr_replicate_weights.md)
  from `rep_weight_prefix`. The weights must be positive, so plain BRR
  weights, which are zero for half of the sample, are not accepted.

- fay_k:

  The Fay coefficient \\k\\ used to make the replicate weights, a number
  with `0 <= fay_k < 1`. Default `0.5`, the PISA value.

- pv_prefix, pv_suffix:

  The text before and after the number in the names of the
  plausible-value columns, used only when `pv_cols` is `NULL`. For PISA
  data give the subject as the suffix, for example `pv_suffix = "READ"`
  for `PV1READ`, `PV2READ`, and so on. The defaults, `"PV"` and `""`,
  match only bare names such as `PV1` and `PV2`.

- rep_weight_prefix:

  The text before the number in the names of the replicate-weight
  columns, used only when `rep_weight_cols` is `NULL`. Default
  `"W_FSTURWT"`, which matches the PISA names `W_FSTURWT1`,
  `W_FSTURWT2`, and so on.

- expected_M, expected_R:

  The numbers of plausible-value and replicate-weight columns you
  expect, as whole numbers (10 and 80 for PISA 2022). Each is compared
  with the columns found or given in `pv_cols` and `rep_weight_cols`,
  and a different number stops the function with an error. `NULL`
  (default) skips the check.

- id_cols:

  A character vector of columns that together identify each row, such as
  a student ID; their combined values must be unique and not missing.
  `NULL` (default) identifies the rows by their position.

- conf_level:

  A confidence level stored in the target, a number strictly between 0
  and 1. Default `0.95`. The estimates, standard errors and degrees of
  freedom do not depend on it, and it does not set the level of the
  intervals that a fit reports: that is the `conf_level` of
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).

- allow_m1:

  Whether one plausible-value column is accepted (`TRUE`); by default
  (`FALSE`) at least two are required. With one plausible value, `B` is
  zero, `T_MI` equals `U_bar`, `fmi` is 0 and the classic degrees of
  freedom are infinite.
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
  cannot use such a target, because a `stack_direct` fit needs at least
  two plausible values.

- df_method:

  The rule for the degrees of freedom: `"classic"` (default) or
  `"barnard_rubin"`, which needs `df_complete`. It also sets the
  interval labels `interval_role` and `coverage_claim_allowed`. See
  "Degrees of freedom and interval labels" in Details.

- df_complete:

  For `df_method = "barnard_rubin"`, the complete-data degrees of
  freedom of the replicate-weight variance, which you state: the degrees
  of freedom that the analysis would have if the outcome were observed
  directly instead of through plausible values. Give one positive number
  for all coefficients, or a vector with one value per coefficient,
  named by the coefficient names (`fe_names`, such as `b_Intercept`).
  The Barnard-Rubin rule matters when this number is small. `Inf` gives
  the classic degrees of freedom, but the target still has
  `coverage_claim_allowed = TRUE`. Default `NULL`. Giving a value with
  `df_method = "classic"` is an error.

- engine:

  How the regressions are fitted. Only `"lm"` (default), weighted least
  squares with [`lm.wfit()`](https://rdrr.io/r/stats/lmfit.html) from
  the stats package, is available; other values are an error.

- verbose:

  If `TRUE`, a message is shown as each plausible value is processed.
  Default `FALSE`.

- x:

  A `pvstackr_brr_target` object from `pv_brr_target()`.

- ...:

  Ignored.

## Value

A `pvstackr_brr_target` object, a list with these fields:

- `beta`, `beta_bar`:

  The estimates \\\bar\beta\\, one per coefficient; the two fields are
  identical.

- `U_bar`, `B`, `T_MI`, `total_var`:

  The matrices \\\bar U\\, \\B\\ and \\T\_{\mathrm{MI}}\\; `total_var`
  repeats `T_MI`.

- `se`:

  The standard errors, `sqrt(diag(T_MI))`.

- `df`, `df_classic`:

  The degrees of freedom by `df_method` and by the classic rule.

- `df_method`, `df_complete`, `conf_level`:

  The settings used. With `"barnard_rubin"`, `df_complete` has one value
  per coefficient; with `"classic"`, it is `NULL`.

- `interval_role`, `coverage_claim_allowed`:

  The interval labels that `df_method` sets (see Details).

- `lambda`, `fmi`, `riv`:

  The ratios defined in Details, one per coefficient; `lambda` and `fmi`
  are identical. The estimate table of a fit has no column for them:
  read them from the target, for example `get_target(fit)$fmi`.

- `fe_names`:

  The coefficient names, which also name the vectors and matrices above:
  the column names of the model matrix with the prefix `b_`, such as
  `b_Intercept` and `b_x`.

- `M`, `R`, `fay_k`, `fay_variance_multiplier`:

  The numbers of plausible values and replicate weights, the Fay
  coefficient, and the factor \\1/(R(1 - k)^2)\\.

- `per_pv`:

  A list with one element per plausible value: `beta` (\\\hat\beta_m\\),
  `U` (\\U_m\\), `replicate_beta` (the replicate estimates, one column
  per replicate weight), `replicate_diff` (their differences from
  `beta`), `pv_col` (the plausible-value column), and copies of
  `fe_names`, `R`, `fay_k` and `fay_variance_multiplier`.

- `formula`, `formula_string`, `rhs_string`:

  The formula, its text and the text of its right-hand side.

- `pv_cols`, `weight_col`, `rep_weight_cols`, `id_cols`:

  The columns used; `id_cols` is `character(0)` when none were given.

- `target_source`, `engine`:

  Always `"external_brr_fay_rubin"` and `"lm"`.

- `design_hash`, `target_hash`:

  SHA-256 checksums of the inputs (rows, plausible values, covariates,
  weights, `fay_k` and formula) and of the contents of the target.
  `design_hash` is not the same checksum as the `design_hash` of a
  [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
  object. The estimate table of a fit repeats `target_hash`.

- `binding_manifest`, `target_content`, `policy`:

  The records that `design_hash` and `target_hash` cover. With them,
  [`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
  and
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
  stop with an error when the target was changed after it was created,
  or when a `stack_direct` fit uses other data than the target. `policy`
  repeats `df_method` and the interval labels, with fixed settings such
  as `fixed_effects_only = TRUE`.

- `schema_version`, `provenance`, `warnings`:

  The format version of the object (`"0.2.0"`), a record of how it was
  made (the function, the time in UTC and the package), and warnings,
  which are empty.

[`print()`](https://rdrr.io/r/base/print.html) shows the numbers of
fixed effects, plausible values and replicate weights, `fay_k`,
`df_method`, `interval_role` and `target_source`, and returns the object
invisibly.

## Details

### Calculation

For each plausible value \\m = 1, \ldots, M\\, the function fits the
regression in `formula`, with the plausible value as the outcome, by
weighted least squares
([`lm.wfit()`](https://rdrr.io/r/stats/lmfit.html) from the stats
package): once with the final weight, giving the coefficients
\\\hat\beta_m\\, and once with each replicate weight \\r = 1, \ldots,
R\\, giving \\\hat\beta_m^{(r)}\\. That is \\M(R + 1)\\ fits in all. The
BRR-Fay replicate covariance of \\\hat\beta_m\\ is \$\$U_m =
\frac{1}{R(1 - k)^2} \sum\_{r=1}^{R} (\hat\beta_m^{(r)} -
\hat\beta_m)(\hat\beta_m^{(r)} - \hat\beta_m)^\top,\$\$ where \\k\\ is
the Fay coefficient `fay_k`. In each replicate, Fay's method multiplies
the weights of one half-sample by \\2 - k\\ and those of the other by
\\k\\, where plain BRR uses 2 and 0; the factor \\1/(1 - k)^2\\ corrects
the variance for this smaller change of the weights (Judkins 1990).

Rubin's rules then combine the \\M\\ results: \$\$\bar\beta =
\frac{1}{M} \sum\_{m=1}^{M} \hat\beta_m, \qquad \bar U = \frac{1}{M}
\sum\_{m=1}^{M} U_m, \qquad B = \frac{1}{M - 1} \sum\_{m=1}^{M}
(\hat\beta_m - \bar\beta)(\hat\beta_m - \bar\beta)^\top,\$\$
\$\$T\_{\mathrm{MI}} = \bar U + (1 + 1/M) B.\$\$ \\\bar U\\ averages the
replicate covariances over the plausible values, and \\B\\ is the
covariance matrix of the \\M\\ estimates \\\hat\beta_m\\. The term
\\B/M\\ is the extra variance of \\\bar\beta\\ that comes from averaging
a finite number of plausible values. The standard errors are the square
roots of the diagonal of \\T\_{\mathrm{MI}}\\.

The fields `lambda` and `fmi` both hold \$\$\lambda = (1 + 1/M) \\
\mathrm{diag}(B) / \mathrm{diag}(T\_{\mathrm{MI}}),\$\$ with the
division taken coefficient by coefficient: the share of the total
variance of each coefficient that is due to the variation between
plausible values. Barnard and Rubin (1999) use this ratio as an
approximation to the fraction of missing information. `riv` holds the
relative increase in variance, \\(1 + 1/M) \\ \mathrm{diag}(B) /
\mathrm{diag}(\bar U)\\.

### Degrees of freedom and interval labels

`df` holds the degrees of freedom of each coefficient by the rule in
`df_method`:

- `"classic"` (default): the rule of Rubin (1987), \\\nu = (M - 1) /
  \lambda^2\\, which is never smaller than \\M - 1\\.

- `"barnard_rubin"`: the rule of Barnard and Rubin (1999), \$\$\nu =
  \left(\frac{1}{\nu\_{\mathrm{classic}}} +
  \frac{1}{\nu\_{\mathrm{obs}}}\right)^{-1}, \qquad \nu\_{\mathrm{obs}}
  = \frac{\nu\_{\mathrm{com}} + 1}{\nu\_{\mathrm{com}} + 3} \\
  \nu\_{\mathrm{com}} \\ (1 - \lambda),\$\$ where
  \\\nu\_{\mathrm{com}}\\ is `df_complete`. This rule corrects the
  classic one when the complete-data degrees of freedom are small. Its
  result is smaller than both the classic value and
  \\\nu\_{\mathrm{com}}\\; with `df_complete = Inf` it equals the
  classic value.

`df_classic` holds the classic value under both rules.

`df_method` also sets the interval labels that the target stores and
that a `stack_direct` fit copies to its estimate table (see
[pvstackr_object_contracts](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)).
`"classic"` gives `interval_role = "descriptive_classic_rubin"` and
`coverage_claim_allowed = FALSE`. `"barnard_rubin"` gives
`interval_role = "coverage_barnard_rubin"` and
`coverage_claim_allowed = TRUE` for any positive `df_complete`, `Inf`
included. By pvstackr's reporting rule, which
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
states in full, only intervals with `coverage_claim_allowed = TRUE` can
be read as confidence intervals with nominal coverage.

## References

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.*
Wiley.

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948-955.

Judkins, D. R. (1990). Fay's method for variance estimation. *Journal of
Official Statistics*, 6(3), 223-239.

## See also

[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
declares the columns,
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
fits a model with this target, and
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md)
returns the target of a fit.

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
