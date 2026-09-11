# The full analysis workflow

Abstract

This article goes through an analysis with pvstackr in order: checking
the plausible-value and weight columns with
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md),
computing the design-based target with
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
fitting the model with
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
and reading the estimate table. The first two steps run on the bundled
synthetic data; the fit is shown as code and read from the example fit
that ships with the package. The article also explains why
`stack_direct` requires `center = "target"` and what your own fitting
functions receive and must return.

``` r

library(pvstackr)
```

To analyze your own data with pvstackr, go through these steps in order:

| Step | Function | Result |
|:---|:---|:---|
| Check the plausible-value and weight columns (optional) | [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md) | the column names, the numbers of plausible values (\\M\\) and replicate weights (\\R\\), and the Fay coefficient |
| Compute the design-based target | [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md) | the estimates of the fixed effects, their covariance matrix, standard errors and degrees of freedom |
| Fit the model | [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md) with settings from [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md) | a fit whose estimate table holds the target’s numbers, with calibrated draws and diagnostics |
| Compare with fits made by other methods (optional) | [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md) | the estimates of two or more fits side by side |
| Read the result | [`summary()`](https://rdrr.io/r/base/summary.html), [`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md) | the estimate table |

The comparison is the subject of [Comparing the three fitting
methods](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.md).
The target and the fit use the same model, a linear regression of the
plausible values on the covariates, weighted by the survey weights, and
only its fixed effects are reported. pvstackr fits no school-level
variance, and the only other parameter of the stacked fit, the residual
standard deviation `sigma`, is neither calibrated nor reported.

The examples use `pisa_tiny`, a small synthetic data set included in the
package. The first two steps run below. Fitting needs a Bayesian
sampler, so the article shows the fitting calls without running them and
reads the example fit that comes with pvstackr. The example fit and its
synthetic data are described in [Getting started with
pvstackr](https://joonho112.github.io/pvstackr/articles/a1-getting-started.html#example-fit).

## Checking the columns with `pv_design()`

[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
finds and checks the columns that hold the plausible values, the final
survey weight and the replicate weights, and records the Fay coefficient
\\k\\ that was used to make the replicate weights. It fits no model. The
step is optional:
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
takes the same column arguments and finds PISA-style plausible-value and
replicate-weight columns itself; there too, give `pv_suffix` for names
such as `PV1READ`, and name the final weight, which has no default. A
design object is a quick check that a new file has the columns you
expect before anything is computed. Read the bundled data and check its
columns:

``` r

pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

design <- pv_design(
  pisa_tiny,
  formula     = OUTCOME ~ x + female,   # OUTCOME stands for each plausible value
  pv_suffix   = "READ",                 # finds PV1READ and PV2READ
  expected_M  = 2L,                     # stop unless 2 plausible values are found
  expected_R  = 4L,                     # stop unless 4 replicate weights are found
  id_cols     = "CNTSTUID"             # the column that identifies each row
)

design
#> pvstackr design
#>   rows: 12
#>   formula: OUTCOME ~ x + female
#>   plausible values: 2
#>   final weight: W_FSTUWT
#>   replicate weights: 4
#>   fay_k: 0.5
#>   design hash: fa78c04b
```

The print shows what was found: 12 rows, the formula
`OUTCOME ~ x + female`, 2 plausible values (\\M = 2\\), the final weight
`W_FSTUWT`, 4 replicate weights (\\R = 4\\), the Fay coefficient
`fay_k = 0.5`, and `design hash`, a checksum of the formula, the
declared columns and their values.

The formula has the placeholder `OUTCOME` on the left-hand side. You
write the model once, and pvstackr replaces `OUTCOME` by each
plausible-value column in turn (`PV1READ`, then `PV2READ`). The survey
weights do not go into the formula: a
[`weights()`](https://rdrr.io/r/stats/weights.html) term is an error.
Name the weight columns in `weight_col` and `rep_weight_cols`; the
defaults of
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
match the PISA names `W_FSTUWT` and `W_FSTURWT1`, `W_FSTURWT2`, and so
on.

PISA 2022 files name the plausible values by subject, such as `PV1READ`,
`PV2READ`, … for reading, so `pv_suffix = "READ"` selects the reading
values; the default suffix `""` matches only bare names such as `PV1`.
`expected_M` and `expected_R` are optional: they stop the function when
a different number of columns is found, for example because the suffix
is wrong. `id_cols` names the column that identifies each row (here the
student ID); its values must be unique. [Real PISA data: access, memory
and
reproducibility](https://joonho112.github.io/pvstackr/articles/a5-real-pisa-guidance.html#design)
shows the call for a PISA 2022 file.

[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
checks only the form of the formula; the covariates are checked by
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
which stops with an error if one of them has missing values, so remove
incomplete rows first. No other pvstackr function takes the design
object itself. Its fields hold the columns that were found, and you pass
them on to
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md):

``` r

design$pv_cols
#> [1] "PV1READ" "PV2READ"
design$weight_col
#> [1] "W_FSTUWT"
design$rep_weight_cols
#> [1] "W_FSTURWT1" "W_FSTURWT2" "W_FSTURWT3" "W_FSTURWT4"
```

[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
finds the columns with two functions that you can also call on their
own, for example to look at a new file:

``` r

detect_pisa_pv_columns(pisa_tiny, suffix = "READ")
#> [1] "PV1READ" "PV2READ"
detect_pisa_brr_replicate_weights(pisa_tiny)
#> [1] "W_FSTURWT1" "W_FSTURWT2" "W_FSTURWT3" "W_FSTURWT4"
```

## Computing the target with `pv_brr_target()`

The target is the design-based estimate that the fit is calibrated to
and whose numbers it reports.
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
computes it without a Bayesian model. For each plausible value it fits
the regression by weighted least squares, once with the final weight and
once with each of the \\R\\ replicate weights, and estimates the
covariance matrix of the coefficients from the spread of the replicate
estimates. It then combines the \\M\\ results with Rubin’s rules. That
makes \\M(R + 1)\\ weighted least-squares fits, 10 here. Pass the data,
the formula and the columns from the design:

``` r

target <- pv_brr_target(
  pisa_tiny,
  formula         = OUTCOME ~ x + female,
  pv_cols         = design$pv_cols,
  weight_col      = design$weight_col,
  rep_weight_cols = design$rep_weight_cols,
  fay_k           = design$fay_k,
  id_cols         = design$id_cols
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

The print shows 3 fixed effects, 2 plausible values, 4 replicate
weights, `fay_k = 0.5`, the degrees-of-freedom rule
(`df method: classic`), the interval label `descriptive_classic_rubin`
and the source label `external_brr_fay_rubin`. A `stack_direct` fit
copies the interval label into its estimate table and repeats the source
label in the column `target_source`.

For plausible value \\m\\, let \\\hat\beta_m\\ be the estimates with the
final weight and \\\hat\beta_m^{(r)}\\ those with replicate weight
\\r\\. The BRR–Fay replicate covariance of \\\hat\beta_m\\ is

\\ \hat U_m = \frac{1}{R(1-k)^2} \sum\_{r=1}^{R}
\big(\hat\beta_m^{(r)} - \hat\beta_m\big)\big(\hat\beta_m^{(r)} -
\hat\beta_m\big)^\top , \\

where \\k\\ is the Fay coefficient `fay_k`. In each replicate, Fay’s
method multiplies the weights of one half-sample by \\2 - k\\ and those
of the other by \\k\\, where plain BRR uses 2 and 0; the factor
\\1/(1-k)^2\\ corrects the variance for this smaller change of the
weights (Judkins 1990). Rubin’s rules then combine the \\M\\ results
(Rubin 1987): the estimate \\\bar\beta\\ is the average of the
\\\hat\beta_m\\, and its covariance matrix is

\\ T\_{\text{MI}} = \bar U + \Big(1 + \frac{1}{M}\Big) B , \\

where \\\bar U\\, the average of the \\\hat U_m\\, is the sampling
covariance within a plausible value, and \\B\\, the covariance matrix of
the \\M\\ estimates \\\hat\beta_m\\, is the variation between plausible
values. The standard errors are the square roots of the diagonal of
\\T\_{\text{MI}}\\.

The degrees of freedom of each coefficient \\\ell\\ follow the rule
chosen with `df_method`:

- `"classic"` (the default) is the rule of Rubin (1987), \\\nu\_\ell =
  (M - 1)/\lambda\_\ell^2\\, where \\\lambda\_\ell = (1 +
  1/M)\\B\_{\ell\ell}/T\_{\text{MI},\ell\ell}\\ is the share of the
  coefficient’s total variance that comes from the variation between
  plausible values (stored in the target as `fmi` and as `lambda`).
- `"barnard_rubin"` is the rule of Barnard and Rubin (1999), which
  corrects the classic value when the complete-data degrees of freedom
  are small. It needs `df_complete`, the complete-data degrees of
  freedom of the replicate-weight variance, which you state.

The rule also sets the interval labels that the fit copies. A target
made with `"classic"`, like this one, gives descriptive intervals; only
a target made with `"barnard_rubin"` gives
`interval_role = "coverage_barnard_rubin"` and
`coverage_claim_allowed = TRUE`. [Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
explains these columns.

The target keeps each piece as a field: `beta` (\\\bar\beta\\), `U_bar`
(\\\bar U\\), `B`, `T_MI`, `se`, `df`, `fmi` (\\\lambda\_\ell\\),
`fay_variance_multiplier` (the factor \\1/\\R(1-k)^2\\\\, which is 1
here) and `target_hash`, a checksum of its contents:

``` r

target$beta
#> b_Intercept         b_x    b_female 
#>  457.894088   46.883361    2.143702
target$T_MI
#>             b_Intercept        b_x  b_female
#> b_Intercept   1.6571716  0.4378132 -4.553496
#> b_x           0.4378132  0.1382299 -1.219646
#> b_female     -4.5534961 -1.2196459 12.638245
sqrt(diag(target$T_MI))   # these become the reported standard errors
#> b_Intercept         b_x    b_female 
#>   1.2873118   0.3717929   3.5550309
```

These are the standard errors that a `stack_direct` fit reports. [The
design-based target: BRR–Fay replicate weights and Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.md)
derives the formulas above.

The target has fixed effects only.
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
stops with an error when the formula has a random-effect term such as
`(1 | CNTSCHID)`, and so do `stack_direct` and `stack_psis` fits:

``` r

pv_brr_target(
  pisa_tiny,
  formula         = OUTCOME ~ x + (1 | CNTSCHID),
  pv_cols         = design$pv_cols,
  weight_col      = design$weight_col,
  rep_weight_cols = design$rep_weight_cols,
  fay_k           = design$fay_k
)
#> Error:
#> ! pvstackr: Random-effect formula terms are not supported by the base WLS BRR-Fay target engine yet.
```

The message means that in pvstackr 0.2.x a target, and therefore a
`stack_direct` fit, can have fixed effects only.

## Fitting the model with `pv_fit()`

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
with `method = "stack_direct"`, the default, fits one Bayesian model to
the stacked data: \\M\\ copies of the data, one for each plausible
value, with that plausible value as the outcome, so \\N \cdot M\\ rows
for \\N\\ students. Row \\i\\ of each copy has the weight \\\tilde
w_i/M\\, where \\\tilde w_i = W_i/\overline W\\ is the student’s final
weight \\W_i\\ divided by the mean final weight \\\overline W\\, so the
stacked log-likelihood is the average of the \\M\\ survey-weighted
log-likelihoods, one for each plausible value. With these weights, the
weighted least-squares estimate from the stacked data equals the
target’s \\\bar\beta\\, the average of the estimates for the separate
plausible values. The stacked fixed-effect point identity, Theorem 4.1
of the companion preprint (Lee et al. 2026), states the conditions for
this equality; [One stacked fit and the Rubin point
estimate](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.md)
explains it.

The identity concerns the point estimate only. The covariance of the
stacked draws does not contain the variation between plausible values,
\\B\\ (Lee et al. 2026), and its part within a plausible value comes
from the model, not from the replicate weights. pvstackr therefore
calibrates the stacked fit to the target. The Cholesky calibration
correction (CCC) transforms the fixed-effect draws so that their mean
and covariance equal the target’s \\\bar\beta\\ and \\T\_{\text{MI}}\\.

### Loading the example fit

The example fit that comes with pvstackr was made from the same data and
target:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit

fit
#> pvstackr fit
#>   method: stack_direct
#>   status: ok
#>   fixed effects: 3
#>   target: external_brr_fay_rubin
#>   draws: not retained
#>   diagnostics: preflight, sampler, sampler_gate, stack_fit, stack_fit_warnings, ccc
#>   interval note: intervals are descriptive rather than coverage-claimable.
```

The print shows the method, the status, the number of fixed effects, the
source label of the target, that the calibrated draws were not kept, the
names of the diagnostics, and a note that the intervals are descriptive;
the note appears because the target uses the classic rule for the
degrees of freedom (see the degrees-of-freedom rules in [Computing the
target](#target)). `status: ok` means that none of pvstackr’s checks
gave a warning or a block, not that a sampler converged; “Status and
checks” in
[`?pvstackr_object_contracts`](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
lists the checks and their thresholds.

### The call for your own data

On your own data, call
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
with the same data and formula that built the target;
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
stops with an error when they differ. The fit also needs an engine. The
bundled one, selected with `pv_control(backend = "brms")`, fits the
stacked model with the brms function `brm()`:

``` r

fit <- pv_fit(
  data    = pisa_tiny,
  formula = OUTCOME ~ x + female,
  target  = target,              # from pv_brr_target() above
  method  = "stack_direct",
  control = pv_control(method = "stack_direct", backend = "brms")
)
```

The bundled engine needs the brms and posterior packages, which pvstackr
suggests but does not require. If cmdstanr is installed, CmdStan must be
configured, otherwise the fit stops with an error instead of switching
to rstan; if cmdstanr is not installed, brms uses rstan.
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
also sets the sampler: by default 4 chains of 2000 iterations each, of
which 1000 are warm-up. `cache_dir`, an argument of
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
that is passed on to
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
defaults to `"cache"`: the bundled engine saves the brms fit in a folder
`cache` in the working directory, which it creates. brms reuses that
file as long as the model and the data are unchanged, even when
`chains`, `iter`, `warmup` or `seed` change: a fit with another number
of chains or of draws per chain is then blocked, and a new seed returns
the old draws without a warning. Set `cache_dir = NULL`, or give a new
`cache_stem`, when you change these settings (see
[`?pv_fit_direct`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)).

To use another engine, pass three functions of your own instead; [Using
your own fitting functions](#own-engine) describes what each one
receives and returns:

``` r

fit <- pv_fit(
  data              = pisa_tiny,
  formula           = OUTCOME ~ x + female,
  target            = target,
  method            = "stack_direct",
  control           = pv_control(method = "stack_direct"),
  fit_function      = your_fit_function,      # fits the stacked model
  draws_function    = your_draws_function,    # returns the draws of that fit
  diagnose_function = your_diagnose_function, # returns its sampler diagnostics
  cache_dir         = NULL                    # no cache file for your function
)
```

## Reading the result

[`summary()`](https://rdrr.io/r/base/summary.html) repeats the overview
of [`print()`](https://rdrr.io/r/base/print.html) and adds the main
columns of the estimate table:

``` r

summary(fit)
#> pvstackr fit summary
#>   method: stack_direct
#>   status: ok
#>   fixed effects: 3
#>   target: external_brr_fay_rubin
#>   draws: not retained
#>   diagnostics: preflight, sampler, sampler_gate, stack_fit, stack_fit_warnings, ccc
#>   interval note: intervals are descriptive rather than coverage-claimable.
#>         term   estimate        se       df  conf_low conf_high
#>  b_Intercept 457.894088 1.2873118 1.021194 442.31804 473.47013
#>          b_x  46.883361 0.3717929 1.402308  44.41457  49.35215
#>     b_female   2.143702 3.5550309 1.013730 -41.60687  45.89428
```

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
returns the estimate table as a data frame with one row per fixed effect
and 17 columns, of which the chunk below shows eight:

``` r

est <- get_estimates(fit)

est[, c("term", "estimate", "se", "df",
        "conf_low", "conf_high",
        "interval_role", "coverage_claim_allowed")]
#>          term   estimate        se       df  conf_low conf_high
#> 1 b_Intercept 457.894088 1.2873118 1.021194 442.31804 473.47013
#> 2         b_x  46.883361 0.3717929 1.402308  44.41457  49.35215
#> 3    b_female   2.143702 3.5550309 1.013730 -41.60687  45.89428
#>               interval_role coverage_claim_allowed
#> 1 descriptive_classic_rubin                  FALSE
#> 2 descriptive_classic_rubin                  FALSE
#> 3 descriptive_classic_rubin                  FALSE
```

These are the rows that [Getting started with
pvstackr](https://joonho112.github.io/pvstackr/articles/a1-getting-started.md)
reads. The estimates and standard errors are the target’s \\\bar\beta\\
and \\\sqrt{\operatorname{diag}(T\_{\text{MI}})}\\ computed above, for
the reason given in [the next section](#center-target), and `df` holds
the target’s classic degrees of freedom, between 1.0 and 1.4 here, which
is why the intervals are so wide. [Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
explains the last two columns.

## Why `stack_direct` requires `center = "target"`

The setting `center` of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
says where the calibrated draws are centered. `stack_direct` requires
the default, `center = "target"`: CCC gives the fixed-effect draws the
target’s mean \\\bar\beta\\ and covariance \\T\_{\text{MI}}\\, and the
estimate table copies the target’s numbers. `estimate` is \\\bar\beta\\,
`se` is \\\sqrt{\operatorname{diag}(T\_{\text{MI}})}\\, and `df`,
`df_method`, `df_complete`, `interval_role` and `coverage_claim_allowed`
come from the target; none of them is computed from the stacked fit. The
stacked fit supplies the calibrated draws, which you can use for
quantities derived from the coefficients, and the diagnostics that
compare the fit with the target.

This is pvstackr’s rule. It makes every number in a row of the estimate
table come from one design-based calculation: the estimate, its standard
error, the degrees of freedom and the interval label are all the
target’s. The mean of the stacked draws differs from \\\bar\beta\\ by
Monte Carlo error, and further when the prior is not flat, so pvstackr
does not report it. It uses it in a check instead: before the
calibration, `delta_c_max` measures the largest difference between a
target estimate and the mean of the stacked draws, in target standard
errors, and a value of 0.01 or more gives the fit the status
`"warning"`, 0.05 or more the status `"blocked"`. [Calibrating the
fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md)
describes the calibration and its diagnostics.

[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
also accepts `center = "posterior"`, which would keep the mean of the
stacked draws and calibrate only their covariance. No fitting function
uses it:
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
with `method = "stack_direct"` stops with an error when it is set, and
`per_pv` and `stack_psis` ignore `center`.

The interval labels come from the target as well, and the calibration
does not change them. The example target uses the classic rule, so every
row of the estimate table above is descriptive. [Reading and reporting
the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
gives the rule that sets these labels for every method.

## Using your own fitting functions

Instead of the bundled engine, `stack_direct` can use three functions
that you pass to
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md):
`fit_function` fits the stacked model, `draws_function` returns its
draws, and `diagnose_function` returns its sampler diagnostics. When you
pass `fit_function`, pvstackr uses it whatever the `backend` value
(leave the default; the value is passed to your function and recorded in
the fit), `draws_function` is required as well, and without
`diagnose_function` the fit is blocked. Without `fit_function`,
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
stops with an error unless `backend = "brms"` selects the bundled
engine. The skeleton below is not run; its comments say what each
function receives and must return:

``` r

# fit_function: fit the stacked model and return the fitted object. pvstackr
# calls it with the named arguments formula, data, family, prior, chains, iter,
# warmup, cores, seed, backend, file and file_refit, plus the elements of
# additional_args. There is no weights argument: the weights are in the
# formula term weights(.pvstackr_weight) and in the data column
# .pvstackr_weight.
my_fit_function <- function(formula, data, ...) {
  # ... call your Bayesian engine here ...
}

# draws_function: return the draws of that fit as a numeric matrix with one
# row per draw. Name the fixed-effect columns as in target$fe_names
# (b_Intercept, b_x, ...) or after the stacked formula (b_pvstackrMM001, ...),
# or pass param_map.
my_draws_function <- function(backend_fit, ...) {
  # ... extract a draws matrix from backend_fit ...
}

# diagnose_function: return the sampler diagnostics as a named list with six
# values. If it is missing, fails or leaves out a value, the fit is blocked.
my_diagnose_function <- function(backend_fit, ...) {
  # ... return list(rhat_max = , ess_bulk_min = , ess_tail_min = ,
  #                 divergences = , chains = , post_warmup_draws_per_chain = )
}

fit <- pv_fit(
  data              = your_real_pisa_data,
  formula           = OUTCOME ~ x + female,
  target            = your_target,
  method            = "stack_direct",
  control           = pv_control(method = "stack_direct"),
  fit_function      = my_fit_function,
  draws_function    = my_draws_function,
  diagnose_function = my_diagnose_function,
  cache_dir         = NULL
)
```

`fit_function` receives the stacked data: one copy of your data per
plausible value, with the outcome in `.pvstackr_y`, the row weight
\\\tilde w_i/M\\ in `.pvstackr_weight`, and the columns of the design
matrix, the intercept included, as `pvstackrMM001`, `pvstackrMM002`, and
so on. For the model of this article its formula, in brms syntax, is
`.pvstackr_y | weights(.pvstackr_weight) ~ 0 + pvstackrMM001 + pvstackrMM002 + pvstackrMM003`.
A function that ignores the weights fits an unweighted model. It may
return any object that your other two functions accept.

`draws_function` takes that object and returns the draws as a numeric
matrix or data frame with one row per draw (at least two). The columns
whose names start with `b_` are the fixed effects; a `sigma` column is
taken as the residual standard deviation, which is not reported, and
other columns such as `lp__` are dropped. When your fixed effects have
other names, say which columns they are with `param_map`.

`diagnose_function` takes the same object and returns a named list with
six values: `rhat_max` (the largest R-hat), `ess_bulk_min` and
`ess_tail_min` (the smallest bulk and tail effective sample sizes,
totals over all chains), `divergences` (the number of divergent
transitions), `chains` and `post_warmup_draws_per_chain`. If the
function is missing or fails, or any of the six values is missing, the
last two included, the fit is blocked (`sampler_diagnostics_incomplete`)
and has no estimates. The last two must also equal `chains` and
`iter - warmup` in
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
or the fit is blocked as well. pvstackr then checks R-hat, the effective
sample sizes and the divergences by the thresholds under “Status and
checks” in
[`?pvstackr_object_contracts`](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md).

The other arguments of `fit_function` carry settings: `chains`, `iter`,
`warmup`, `cores`, `seed` and `backend` come from
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
`family` is always
[`stats::gaussian()`](https://rdrr.io/r/stats/family.html), `prior` is
`NULL` unless you give one (see
[`?pv_fit_direct`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)),
and `additional_args`, a named list in the
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
call, adds arguments of your own. With the default
`cache_dir = "cache"`, your function receives
`file = "cache/pvstackr-stack-direct"` and `file_refit = "on_change"`,
but pvstackr does not create the folder; with `cache_dir = NULL` it
receives `file = NULL` and `file_refit = "never"`.

You do not have to write all three functions. The functions of the
bundled engine are exported as
[`pv_backend_brms_fit_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md),
[`pv_backend_brms_draws_function()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
and
[`pv_backend_brms_sampler_diagnostics()`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md),
so you can pass them and replace only the one that differs for your
engine.
[`?pv_fit_direct`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
gives the full lists of arguments and return values (see `fit_function`,
`draws_function`, `param_map`, `diagnose_function`, `cache_dir` and
`additional_args`), and
[`?pv_backend_brms_fit_function`](https://joonho112.github.io/pvstackr/reference/pv_backend_brms_fit_function.md)
describes the bundled functions.

Once you have a fit, [Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.md)
describes everything it contains and what to report. For real PISA
files, [Real PISA data: access, memory and
reproducibility](https://joonho112.github.io/pvstackr/articles/a5-real-pisa-guidance.md)
covers data access, memory and running time: the stacked data have \\N
\cdot M\\ rows, and the target takes \\M(R + 1)\\ weighted least-squares
fits.

## References

Barnard, John, and Donald B. Rubin. 1999. “Small-Sample Degrees of
Freedom with Multiple Imputation.” *Biometrika* 86 (4): 948–55.
<https://doi.org/10.1093/biomet/86.4.948>.

Judkins, David R. 1990. “Fay’s Method for Variance Estimation.” *Journal
of Official Statistics* 6 (3): 223–39.
<https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/fay39s-method-for-variance-estimation.pdf>.

Lee, JoonHo, Matthew R. Williams, and Terrance D. Savitsky. 2026. *One
Markov Chain Monte Carlo Fit for Many Plausible Values: A Calibrated
Stacked Posterior Workflow for Bayesian Multilevel Models of Large-Scale
Assessment Data*. Zenodo preprint, version 1.
<https://doi.org/10.5281/zenodo.22407935>.

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.
