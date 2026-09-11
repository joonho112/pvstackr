# Reading and reporting the results

Abstract

This article reads a fitted object with
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md),
[`get_target()`](https://joonho112.github.io/pvstackr/reference/get_target.md),
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
and
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
It explains the columns of the estimate table, including the columns
that say whether pvstackr’s reporting rule lets you read an interval as
a confidence interval with nominal coverage, and shows where to find the
fraction of missing information. It ends with a list of what to report.

``` r

library(pvstackr)
```

## Loading a fit

The examples read the `stack_direct` fit that ships with pvstackr:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
```

The example fit and its synthetic data (\\M = 2\\ plausible values, \\R
= 4\\ replicate weights) are described in [Getting started with
pvstackr](https://joonho112.github.io/pvstackr/articles/a1-getting-started.html#example-fit).

## Functions that read a fit

Four functions read a fit. Use them rather than `$`: they stop with an
error if the fit was changed after it was created.

- `get_estimates(fit)` returns the estimate table: one row per fixed
  effect, with the estimate, its standard error, degrees of freedom and
  interval, and the columns that say how the interval can be read (see
  [the estimate table](#estimate-table)).
- `get_target(fit)` returns the target that the fit was calibrated to,
  the `pvstackr_brr_target` object made by
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
  It holds the estimates combined over the plausible values, their
  covariance matrix \\T\_{\text{MI}}\\, the degrees of freedom and the
  fraction of missing information.
- `get_draws(fit)` returns the calibrated fixed-effect draws of the
  stacked fit, or `NULL` when the fit was saved without them.
- `get_diagnostics(fit)` returns the results of pvstackr’s checks and
  other records of the fit. For a `stack_direct` fit it is a list with
  six elements: `preflight`, `sampler`, `sampler_gate`, `stack_fit`,
  `stack_fit_warnings` and `ccc`.

[The target, the draws and the diagnostics](#target-and-diagnostics)
shows the last three in more detail. The example fit was saved without
its draws (`return_draws = FALSE` in
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)),
so
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
returns `NULL`:

``` r

get_draws(fit)
#> NULL
```

The estimate table does not need the draws, because for a `stack_direct`
fit its values come from the target.

## The estimate table

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
returns a data frame with one row per fixed effect. For a `stack_direct`
fit it has 17 columns:

``` r

est <- get_estimates(fit)

names(est)
#>  [1] "term"                   "estimate"               "se"                    
#>  [4] "std.error"              "df"                     "df_method"             
#>  [7] "df_complete"            "conf_level"             "conf_low"              
#> [10] "conf_high"              "conf.low"               "conf.high"             
#> [13] "interval_role"          "coverage_claim_allowed" "parameter_scope"       
#> [16] "target_source"          "target_hash"
```

- `term` names the coefficient: the column name of the model matrix with
  the prefix `b_`, here `b_Intercept`, `b_x` and `b_female`. Only the
  fixed effects (the intercept and the slopes) have rows; other
  parameters of the model, such as the residual standard deviation
  `sigma`, are not reported.
- `estimate`, `se` and `df` are the estimate, its standard error and its
  degrees of freedom.
- `conf_level`, `conf_low` and `conf_high` give the interval, and
  `df_method`, `df_complete`, `interval_role` and
  `coverage_claim_allowed` say how it was made and how it can be read
  (see [the interval columns](#interval-columns)).
- `parameter_scope`, `target_source` and `target_hash` say what each row
  is and where its numbers come from (see below).
- `std.error`, `conf.low` and `conf.high` repeat `se`, `conf_low` and
  `conf_high` under the names that the broom package uses for model
  summaries. pvstackr has no `tidy()` method; the repeated columns only
  let code written for these names read the table.

The first columns hold the numbers that most analyses start from:

``` r

est[, c("term", "estimate", "se", "df", "df_method")]
#>          term   estimate        se       df df_method
#> 1 b_Intercept 457.894088 1.2873118 1.021194   classic
#> 2         b_x  46.883361 0.3717929 1.402308   classic
#> 3    b_female   2.143702 3.5550309 1.013730   classic
```

In the example the intercept is 457.9, the slope of `x` is 46.9 and the
coefficient of `female` is 2.1, with standard errors 1.29, 0.37 and
3.56. The degrees of freedom are close to 1, which makes the intervals
very wide. With \\M = 2\\ plausible values the classic degrees of
freedom cannot be smaller than \\M - 1 = 1\\, and they are close to that
value here because most of the variance comes from the differences
between the two plausible values (see [the fraction of missing
information](#fmi)).

The last three columns record what each row is and where it comes from:

``` r

est[, c("term", "parameter_scope", "target_source", "target_hash")]
#>          term parameter_scope          target_source
#> 1 b_Intercept    fixed_effect external_brr_fay_rubin
#> 2         b_x    fixed_effect external_brr_fay_rubin
#> 3    b_female    fixed_effect external_brr_fay_rubin
#>                                                               target_hash
#> 1 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
#> 2 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
#> 3 sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa
```

`parameter_scope` is `"fixed_effect"` in every row. `target_source` is
`"external_brr_fay_rubin"` for a `stack_direct` fit: the numbers come
from the target that
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
computed with the BRR–Fay replicate weights and Rubin’s rules.
`target_hash` is a SHA-256 checksum of the contents of that target. The
target stores the same value, and pvstackr uses it to detect a target
that was changed after it was made. Keep it with your results, since it
identifies the target they come from; [Real PISA data: access, memory
and
reproducibility](https://joonho112.github.io/pvstackr/articles/a5-real-pisa-guidance.html#reproducibility-checklist)
lists the other records to keep.

For a `stack_direct` fit, `estimate`, `se` and `df` are copied from the
target ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)):
`estimate` is the target’s \\\bar\beta\\, the average of the \\M\\
survey-weighted estimates, one per plausible value, and `se` is
\\\sqrt{\operatorname{diag}(T\_{\text{MI}})}\\, the square root of the
diagonal of the target covariance matrix.

## Interval columns and pvstackr’s reporting rule

Each row has an interval from `conf_low` to `conf_high`: the estimate
plus or minus a \\t\\ quantile with `df` degrees of freedom times `se`,
at the level `conf_level` set in
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
(0.95 by default). Four columns say how the interval was made and how it
can be read:

``` r

est[, c("term", "interval_role", "df_method",
        "df_complete", "coverage_claim_allowed")]
#>          term             interval_role df_method df_complete
#> 1 b_Intercept descriptive_classic_rubin   classic          NA
#> 2         b_x descriptive_classic_rubin   classic          NA
#> 3    b_female descriptive_classic_rubin   classic          NA
#>   coverage_claim_allowed
#> 1                  FALSE
#> 2                  FALSE
#> 3                  FALSE
```

- `df_method` is the rule for the degrees of freedom, chosen when the
  target is made with
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
  `"classic"`, the default, is the rule of Rubin (1987), \\(M -
  1)/\lambda^2\\, where \\\lambda\\ is the fraction of missing
  information of the coefficient (see [the fraction of missing
  information](#fmi)). `"barnard_rubin"` is the correction of Barnard
  and Rubin (1999) for small complete-data degrees of freedom ([The full
  analysis
  workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#target)
  describes both rules).
- `df_complete` is used only with `"barnard_rubin"`: the complete-data
  degrees of freedom of the replicate-weight variance, which you state
  when you make the target. They are the degrees of freedom the analysis
  would have if the outcome were observed directly instead of through
  plausible values, and the Barnard–Rubin adjustment matters when this
  number is small. pvstackr does not choose this number and does not
  check it against the replicate weights:
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  accepts any positive value, `Inf` included. The value depends on how
  the survey formed its replicate weights, which the survey’s technical
  documentation describes; report it with the intervals. With
  `"classic"` the column is `NA`, as here.
- `coverage_claim_allowed` is `TRUE` when pvstackr’s reporting rule lets
  you read the interval as a confidence interval with nominal coverage
  (a coverage-claimable interval), and `FALSE` when pvstackr makes no
  statement about its coverage (a descriptive interval).
- `interval_role` names the case, as in the table below.

pvstackr sets `interval_role` and `coverage_claim_allowed` by its
reporting rule, from the method and `df_method`. For `stack_direct`,
`df_method` and `df_complete` are those of the target; `per_pv` and
`stack_psis` take them as their own arguments.

| Method | `df_method` | `interval_role` | `coverage_claim_allowed` |
|:---|:---|:---|:---|
| `stack_direct` | `"barnard_rubin"` | `coverage_barnard_rubin` | `TRUE` |
| `stack_direct` | `"classic"` | `descriptive_classic_rubin` | `FALSE` |
| `per_pv` | `"barnard_rubin"` | `reference_barnard_rubin` | `FALSE` |
| `per_pv` | `"classic"` | `reference_classic_rubin` | `FALSE` |
| `stack_psis` | `"barnard_rubin"` | `psis_barnard_rubin` | `FALSE` |
| `stack_psis` | `"classic"` | `psis_classic_rubin` | `FALSE` |

So `coverage_claim_allowed` is `TRUE` only for a `stack_direct` fit
whose target was made with `df_method = "barnard_rubin"` and a positive
`df_complete`
([`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
stops with an error without one). The labels do not depend on the status
of the fit: a fit with status `"warning"` keeps them, and a fit with
status `"blocked"` has an empty estimate table. `per_pv` and
`stack_psis` intervals are always descriptive, because their variance
within each plausible value comes from the model’s posterior draws, not
from the replicate weights. The label follows the rule, not the value of
the degrees of freedom: the classic rule assumes that the complete-data
degrees of freedom are infinite (Barnard and Rubin 1999, 948), and a
target made with the Barnard–Rubin rule keeps the label
`coverage_barnard_rubin` even with `df_complete = Inf`, which gives the
classic degrees of freedom ([The design-based target: BRR–Fay replicate
weights and Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.html#degrees-of-freedom)
explains both rules). The label records how the interval was built; it
does not certify coverage ([Choosing a method and reading its
intervals](https://joonho112.github.io/pvstackr/articles/m5-methods-and-coverage.html#interval-rule)
gives the reason for the rule and its limits).

The example fit is a `stack_direct` fit whose target uses the classic
degrees of freedom, so every row is `descriptive_classic_rubin` with
`coverage_claim_allowed = FALSE`, and `print(fit)` adds the note
“intervals are descriptive rather than coverage-claimable”. The method
alone does not decide the label; read it in the table. Report an
interval as a confidence interval with nominal coverage only when
`coverage_claim_allowed` is `TRUE`, and as a descriptive interval
otherwise.

The calibration diagnostics are not evidence about coverage. The
Cholesky calibration correction (CCC) transforms the draws of the
stacked fit so that their mean and covariance equal the target’s
\\\bar\beta\\ and \\T\_{\text{MI}}\\, whatever the target is. Its
diagnostics (see [the diagnostics](#target-and-diagnostics)) show how
far the stacked fit was from the target before the calibration and
whether the calibration was numerically stable; they cannot show that
the target’s intervals have nominal coverage. [Calibrating the
fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md)
describes the calibration, and [Choosing a method and reading its
intervals](https://joonho112.github.io/pvstackr/articles/m5-methods-and-coverage.md)
explains the intervals of each method.

## Fraction of missing information

For each coefficient \\\ell\\,
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
computes

\\ \lambda\_\ell = \frac{(1 +
1/M)\\B\_{\ell\ell}}{T\_{\text{MI},\ell\ell}}, \qquad T\_{\text{MI}} =
\bar U + (1 + 1/M)\\B, \\

where \\\bar U\\ is the BRR–Fay replicate covariance of the estimates
averaged over the plausible values, and \\B\\ is the covariance of the
estimates between plausible values. \\\lambda\_\ell\\ is the share of
the total variance of coefficient \\\ell\\ that comes from the variation
between plausible values; Barnard and Rubin (1999) use this ratio as an
approximation to the fraction of missing information (FMI). pvstackr
stores it in the target twice, as `fmi` and as `lambda`, which are
identical; the estimate table has no column for it. Since the classic
degrees of freedom are \\(M - 1)/\lambda\_\ell^2\\, an FMI close to 1
gives degrees of freedom close to their smallest value, \\M - 1\\.

The target also stores `riv`, the relative increase in variance \\(1 +
1/M)B\_{\ell\ell}/\bar U\_{\ell\ell}\\, and the degrees of freedom:

``` r

tg <- get_target(fit)

data.frame(
  term = tg$fe_names,
  fmi  = round(as.numeric(tg$fmi), 3),
  riv  = round(as.numeric(tg$riv), 2),
  df   = round(as.numeric(tg$df),  2)
)
#>          term   fmi    riv   df
#> 1 b_Intercept 0.990  94.86 1.02
#> 2         b_x 0.844   5.43 1.40
#> 3    b_female 0.993 146.17 1.01
```

In the example the FMI is between 0.84 and 0.99 because the variance
between the two plausible values is much larger than the
replicate-weight variance: `riv` is about 5 for `b_x` and about 95 and
146 for the intercept and `b_female`. Because `fmi` equals
`riv / (1 + riv)`, a large `riv` gives an FMI close to 1, and the
degrees of freedom are then close to \\M - 1 = 1\\. The small number of
plausible values is not the reason: with the same \\B\\ and \\\bar U\\
and ten plausible values, the FMI would still be between 0.80 and 0.99.

PISA 2022 provides ten plausible values per domain (OECD 2024). The FMI
of an analysis of real data depends on its data and model, and this
synthetic example says nothing about it.

## The target, the draws and the diagnostics

### The target

`tg`, made in the previous chunk, holds the target of the fit. The
estimate table takes its standard errors, degrees of freedom and
checksum from these fields:

``` r

diag(tg$T_MI)        # variance of each coefficient; se is its square root
#> b_Intercept         b_x    b_female 
#>   1.6571716   0.1382299  12.6382448
round(tg$df, 3)      # degrees of freedom of each coefficient
#> b_Intercept         b_x    b_female 
#>       1.021       1.402       1.014
tg$target_hash       # checksum of the target, also in the estimate table
#> [1] "sha256:f173650e9120742a1a6fc6406bfe3ab130e454b17f28e4822cb99e25c108bfaa"
```

The two comparisons below return `TRUE` because pvstackr copies
`estimate` from `beta_bar` and `se` from the square roots of the
diagonal of `T_MI` when `center = "target"`. The equality holds for
every `stack_direct` fit that has estimates; it shows how the table is
made, not that two methods agree.

``` r

all.equal(est$estimate, unname(tg$beta_bar))            # estimate = beta_bar
#> [1] TRUE
all.equal(est$se,       unname(sqrt(diag(tg$T_MI))))    # se = sqrt(diag(T_MI))
#> [1] TRUE
```

[`?pv_brr_target`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
lists all fields of the target. Besides `beta_bar` (also stored as
`beta`), `T_MI`, `df` and `target_hash`, they include `U_bar` and `B`,
the two parts of \\T\_{\text{MI}} = \bar U + (1 + 1/M)B\\; `fmi` and
`lambda`, which hold the same values, and `riv` (see [the fraction of
missing information](#fmi)); `M`, `R` and `fay_k`; and
`fay_variance_multiplier`, the factor \\1/(R(1 - k)^2)\\ of the BRR–Fay
replicate covariance, where \\k\\ is the Fay coefficient `fay_k`. In the
example \\R = 4\\ and \\k = 0.5\\, so the factor is 1. [The design-based
target: BRR–Fay replicate weights and Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.md)
explains how the target is computed: the BRR–Fay replicate covariance
for each plausible value, Rubin’s rules, and the two rules for the
degrees of freedom.

### The draws

A fit made with `return_draws = TRUE`, the default of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
keeps the calibrated fixed-effect draws: a matrix with one row per draw
of the stacked fit and one column per fixed effect, named as in the
estimate table. They can be used for quantities derived from the
coefficients: their mean and covariance are the target’s, so a linear
combination computed from them has the target’s estimate and variance.
They do not carry the target’s degrees of freedom, so their quantiles
give roughly normal-reference intervals, narrower than the \\t\\
intervals of the estimate table when the degrees of freedom are small
([Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.html#map)).
For the example fit, which was saved without them,
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
returns `NULL`, as shown above; it does the same for a blocked fit.

### The diagnostics

``` r

dg <- get_diagnostics(fit)

names(dg)
#> [1] "preflight"          "sampler"            "sampler_gate"      
#> [4] "stack_fit"          "stack_fit_warnings" "ccc"
```

For a `stack_direct` fit the list has six elements:

- `preflight`: the record of the check, made before fitting, that the
  formula, the data and the target match;
- `sampler`: the sampler diagnostics of the stacked fit (the largest
  R-hat, the smallest bulk and tail effective sample sizes, the number
  of divergent transitions, and the numbers of chains and of draws per
  chain);
- `sampler_gate`: the result of pvstackr’s check of those values;
- `stack_fit` and `stack_fit_warnings`: a record of the stacked fit and
  the notes it produced;
- `ccc`: the calibration diagnostics.

The thresholds of the checks are listed under “Status and checks” in
[`?pvstackr_object_contracts`](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md);
a fit has status `"ok"` when none of them gave a warning or a block. The
sampler record of the example fit was typed in by its fitting functions:
2 chains of 10 draws each, an R-hat of 1 and no divergent transitions.
The record reports 20 draws although the functions return 12: pvstackr
compares the reported numbers of chains and draws with the settings in
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
not with the draws, so the record passed pvstackr’s checks. The
center-separation values in `ccc` are:

``` r

dg$ccc[c("center_status", "delta_c_rel", "delta_c_max")]
#> $center_status
#> [1] "ok"
#> 
#> $delta_c_rel
#> [1] 2.549387e-14
#> 
#> $delta_c_max
#> [1] 4.415668e-14
```

`delta_c_max` is the largest, over the fixed effects, of the absolute
difference between the target estimate and the mean of the stacked draws
before calibration, divided by the target standard error:

\\ \max\_\ell \frac{\lvert \bar\beta\_\ell -
\bar\beta^{\text{raw}}\_\ell \rvert} {\sqrt{T\_{\text{MI},\ell\ell}}},
\\

where \\\bar\beta^{\text{raw}}\_\ell\\ is the mean of the stacked draws
of coefficient \\\ell\\ before calibration. By pvstackr’s thresholds,
`center_status` is `"ok"` when `delta_c_max` is below 0.01, `"warning"`
from 0.01 and `"blocked"` from 0.05. `delta_c_rel` is the root mean
square (RMS) of the same ratios; it is recorded but does not change the
status. Both are about \\10^{-14}\\ in the example because its draws
were built around the target. `ccc` also holds `kappa_A`, the condition
number of the calibration matrix, which has thresholds of its own.
[Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md)
explains these values.

## What to report

A checklist for reporting the fixed effects of a pvstackr fit:

- The method (`stack_direct`, `per_pv` or `stack_psis`); [Comparing the
  three fitting
  methods](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.html#choosing)
  gives pvstackr’s recommendation.
- The status shown by `print(fit)` and, for `"warning"`, the reason
  codes it lists.
- `estimate` and `se` from
  [`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md);
  for `stack_direct` they are the target’s ([The full analysis
  workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)
  explains why).
- `conf_low`, `conf_high` and `conf_level`, with the `interval_role` of
  each row.
- The name of each interval: “confidence interval” only where
  `coverage_claim_allowed` is `TRUE`, otherwise “descriptive interval”.
- `df`, `df_method` and, with the Barnard–Rubin rule, `df_complete`.
- The fraction of missing information, `get_target(fit)$fmi`
  (`stack_direct` and `per_pv`).
- `target_hash` and the pvstackr version, `packageVersion("pvstackr")`;
  [Real PISA data: access, memory and
  reproducibility](https://joonho112.github.io/pvstackr/articles/a5-real-pisa-guidance.html#reproducibility-checklist)
  lists the other records to keep.

For the example fit a report would read: “The slope of `x` was 46.9
(standard error 0.37; 95% descriptive interval 44.4 to 49.4, classic
Rubin degrees of freedom 1.40; `stack_direct`, status ok; fraction of
missing information 0.84).” The data are synthetic.

## References

Barnard, John, and Donald B. Rubin. 1999. “Small-Sample Degrees of
Freedom with Multiple Imputation.” *Biometrika* 86 (4): 948–55.
<https://doi.org/10.1093/biomet/86.4.948>.

OECD. 2024. *PISA 2022 Technical Report*. PISA. OECD Publishing.
<https://doi.org/10.1787/01820d6d-en>.

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.
