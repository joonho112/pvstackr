# Choosing a method and reading its intervals

Abstract

pvstackr has three fitting methods, `stack_direct`, `per_pv` and
`stack_psis`. All three combine an estimate and a covariance matrix for
each plausible value with Rubin’s rules; they differ in how they compute
these. This article gives the formulas of the three methods, the
importance weights and Pareto k-hat values that `stack_psis` takes from
another program together with pvstackr’s check of them, and the reason
why pvstackr’s reporting rule labels intervals coverage-claimable only
for `stack_direct` fits whose target uses Barnard–Rubin degrees of
freedom.

``` r

library(pvstackr)
```

pvstackr estimates the fixed effects of a regression on plausible-value
data with one of three methods, set by the argument `method` of
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md):
`stack_direct` (the default), `per_pv` and `stack_psis`. [Comparing the
three fitting
methods](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.md)
compares example fits of the three methods and gives pvstackr’s
recommendation for choosing among them.

The code in this article reads the example fit that ships with pvstackr,
a `stack_direct` fit, and fits no model:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # the target of the fit, made by pv_brr_target()
```

The example fit and its synthetic data are described in [Getting started
with
pvstackr](https://joonho112.github.io/pvstackr/articles/a1-getting-started.html#example-fit).

## The three methods

Write \\\beta\\ for the \\p\\ fixed-effect coefficients and \\m = 1,
\dots, M\\ for the plausible values. Each method computes, for each
plausible value, an estimate \\\hat\beta_m\\ of \\\beta\\ and a
covariance matrix \\\hat U_m\\, and combines them with Rubin’s rules
(Rubin 1987):

\\ \bar\beta = \frac{1}{M}\sum\_{m=1}^{M} \hat\beta_m, \qquad \bar U =
\frac{1}{M}\sum\_{m=1}^{M} \hat U_m, \qquad B =
\frac{1}{M-1}\sum\_{m=1}^{M} (\hat\beta_m - \bar\beta)(\hat\beta_m -
\bar\beta)^\top, \\

\\ T\_{\text{MI}} = \bar U + \Big(1 + \frac{1}{M}\Big) B . \\

For coefficient \\\ell\\, the estimate table reports the estimate
\\\bar\beta\_\ell\\, the standard error
\\\sqrt{T\_{\text{MI},\ell\ell}}\\, and an interval: the estimate plus
or minus a \\t\\ quantile times the standard error, with the classic or
the Barnard–Rubin degrees of freedom and the level `conf_level` of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
(0.95 by default). The standard errors of every method therefore include
the variation between plausible values, \\B\\. [The design-based target:
BRR–Fay replicate weights and Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.html#rubin)
derives these formulas and the degrees of freedom. The methods differ in
\\\hat\beta_m\\ and \\\hat U_m\\:

| Method | \\\hat\beta_m\\ | \\\hat U_m\\ | Bayesian fits | Intervals |
|:---|:---|:---|:---|:---|
| `stack_direct` (default) | weighted least-squares estimate with the final weights | BRR–Fay replicate covariance | 1 | coverage-claimable only with an external Barnard–Rubin target, otherwise descriptive |
| `per_pv` | posterior mean of the fit to plausible value \\m\\ | posterior covariance of that fit | \\M\\ | always descriptive |
| `stack_psis` | weighted mean of the draws of one stacked fit | weighted covariance of the same draws | 1 | always descriptive |

A coverage-claimable interval is one that pvstackr’s reporting rule lets
you read as a confidence interval with nominal coverage; a descriptive
interval only describes the uncertainty of the estimate. The source of
\\\hat U_m\\ is the reason why only `stack_direct` intervals can be
coverage-claimable ([the reason for the rule](#interval-rule)).

The companion preprint (Lee et al. 2026) calls the corresponding
workflows the calibrated stack, the per-PV workflow and the
reweighted-stack variant, and its Appendix G.1 matches them to
`stack_direct`, `per_pv` and `stack_psis`. The correspondence is one of
role. The preprint’s calibrated stack uses by default a target from
model-based fits of each plausible value (Lee et al. 2026, sec. 4.4),
whereas
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
computes only the replicate-weight target, and the preprint’s analysis
model has a school random intercept, which pvstackr 0.2.x does not fit.

### `stack_direct`

`stack_direct`
([`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md))
takes \\\hat\beta_m\\ and \\\hat U_m\\ from the target that
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
computes before the fit. \\\hat\beta_m\\ is the weighted least-squares
estimate for plausible value \\m\\ with the final survey weights, and

\\ \hat U_m = \frac{1}{R(1-k)^2} \sum\_{r=1}^{R}
\big(\hat\beta_m^{(r)} - \hat\beta_m\big)\big(\hat\beta_m^{(r)} -
\hat\beta_m\big)^\top \\

is its BRR–Fay replicate covariance, where \\\hat\beta_m^{(r)}\\ is the
same estimate with replicate weight \\r = 1, \dots, R\\ and \\k\\ is the
Fay coefficient ([The design-based target: BRR–Fay replicate weights and
Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.html#replicate-covariance)).
The target holds \\\bar\beta\\, \\\bar U\\, \\B\\ and
\\T\_{\text{MI}}\\, and the estimate table of the fit copies the
estimates, the standard errors and the degrees of freedom from it.

The one Bayesian fit of `stack_direct` is a fit to the stacked data, one
copy of the data for each plausible value. Its fixed-effect draws
\\\beta_1, \dots, \beta_S\\ do not give the reported numbers: the
covariance of the draws does not contain \\B\\, and it comes from the
model, not from the replicate weights ([Calibrating the fixed-effect
draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.html#why-calibrate)).
The Cholesky calibration correction (CCC) maps each draw to
\\\beta^{\text{cal}}\_s = \bar\beta + A(\beta_s -
\bar\beta^{\text{raw}})\\, where \\\bar\beta^{\text{raw}}\\ is the mean
of the draws and \\A\\ a calibration matrix, so that the calibrated
draws have the mean \\\bar\beta\\ and the covariance \\T\_{\text{MI}}\\
of the target. The stacked fixed-effect point identity, Theorem 4.1 of
the preprint (Lee et al. 2026), gives conditions under which the
fixed-effect estimate of the stacked fit equals \\\bar\beta\\ ([One
stacked fit and the Rubin point
estimate](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.html#point-identity)).

### `per_pv`

`per_pv`
([`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md))
combines \\M\\ Bayesian fits, one for each plausible value. You supply
their draws, either computed elsewhere (`per_pv_draws`) or returned by
fitting functions of your own, which pvstackr calls once for each
plausible value. With \\\beta_1^{(m)}, \dots, \beta\_{S_m}^{(m)}\\ the
fixed-effect draws of the fit to plausible value \\m\\,

\\ \hat\beta_m = \frac{1}{S_m}\sum\_{s=1}^{S_m} \beta_s^{(m)}, \qquad
\hat U_m = \frac{1}{S_m - 1}\sum\_{s=1}^{S_m} \big(\beta_s^{(m)} -
\hat\beta_m\big)\big(\beta_s^{(m)} - \hat\beta_m\big)^\top , \\

the posterior mean and the posterior covariance of the draws. The model
is the one your fits use: pvstackr does not add the survey weights to
its formula and does not check the formula for random-effect terms
([Plausible values, survey weights and
notation](https://joonho112.github.io/pvstackr/articles/m1-foundations-and-notation.html#model)).

### `stack_psis`

`stack_psis`
([`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md))
uses the fixed-effect draws \\\beta_1, \dots, \beta_S\\ of one Bayesian
fit to the stacked data and, for each plausible value, importance
weights \\\omega_1^{(m)}, \dots, \omega_S^{(m)}\\ that reweight these
draws toward the posterior for plausible value \\m\\. pvstackr divides
the weights by their sum, \\\tilde\omega_s^{(m)} = \omega_s^{(m)} \big/
\sum\_{s'=1}^{S} \omega\_{s'}^{(m)}\\, and computes the weighted mean
and the weighted covariance of the draws:

\\ \hat\beta_m = \sum\_{s=1}^{S} \tilde\omega_s^{(m)} \beta_s, \qquad
\hat U_m = \frac{\sum\_{s=1}^{S} \tilde\omega_s^{(m)} \big(\beta_s -
\hat\beta_m\big)\big(\beta_s - \hat\beta_m\big)^\top} {1 -
\sum\_{s=1}^{S} \big(\tilde\omega_s^{(m)}\big)^2} . \\

With equal weights \\1/S\\, these are the mean and the usual sample
covariance of the draws. The weights and a Pareto \\\hat k\\ for each
plausible value come from another program ([importance weights and the
Pareto \\\hat k\\](#psis)).

## Importance weights and the Pareto \\\hat k\\

The importance weights of `stack_psis` reweight the draws of the stacked
fit, which serve as the proposal, toward the posterior for each
plausible value. The importance ratio of a draw is the ratio of that
posterior density to the proposal density at the draw. Suppose the
stacked fit samples the stacked fractional posterior

\\ q(\beta, \sigma) \propto p(\beta, \sigma) \prod\_{m=1}^{M} \tilde
L_m(\beta, \sigma)^{1/M} , \\

where \\\tilde L_m\\ is the survey-weighted likelihood of plausible
value \\m\\, \\\sigma\\ the residual standard deviation and \\p(\beta,
\sigma)\\ the prior ([One stacked fit and the Rubin point
estimate](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.html#stacked-objective)).
If the posterior for plausible value \\m\\ is \\p_m(\beta, \sigma)
\propto p(\beta, \sigma)\\ \tilde L_m(\beta, \sigma)\\, with the same
prior, the log ratio at draw \\(\beta_s, \sigma_s)\\ is

\\ \log \frac{p_m(\beta_s, \sigma_s)}{q(\beta_s, \sigma_s)} = \log
\tilde L_m(\beta_s, \sigma_s) - \frac{1}{M} \sum\_{m'=1}^{M} \log \tilde
L\_{m'}(\beta_s, \sigma_s) + c_m , \\

where \\c_m\\ does not depend on the draw and cancels when the weights
are divided by their sum.

A few draws with very large ratios can dominate raw importance weights.
Pareto smoothed importance sampling (PSIS) fits a generalized Pareto
distribution to the largest ratios, replaces them by values from the
fitted distribution, and reports the estimated shape parameter of that
distribution, the Pareto \\\hat k\\, as a diagnostic of how reliable the
importance-sampling estimates are (Vehtari et al. 2024).

pvstackr computes no part of PSIS: it fits no Pareto distribution,
smooths no weights and estimates no \\\hat k\\. You supply, for each
plausible value, the weights \\\omega_s^{(m)}\\ and the Pareto \\\hat
k_m\\ in one of two ways:

- `psis_weights`, a matrix with one row per stacked draw and one column
  per plausible value, with `pareto_k`, one value per plausible value;
- `log_ratios`, a matrix of log importance ratios of the same shape,
  with `psis_function`, a function of your own that pvstackr calls on
  `log_ratios` and that returns the weights and the \\\hat k_m\\, for
  example by applying the function `psis()` of the loo package to each
  column.

With `log_ratios` and `pareto_k` alone, pvstackr turns the ratios into
normalized weights without smoothing and blocks the fit, with the reason
code `psis_smoothing_not_applied` when the \\\hat k_m\\ pass the check
below (otherwise `psis_k_too_high` or `psis_k_not_evaluated`).
[`?pv_fit_stack_psis`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)
lists all inputs.

By pvstackr’s rule, a `stack_psis` fit has estimates (status `"ok"`)
only when the Pareto \\\hat k_m\\ of every plausible value is finite and
strictly below `psis_k_threshold` in
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
which is 0.7 by default and cannot be set higher, and when
`psis_producer` and `psis_producer_version` name the program that
produced the weights. Otherwise the fit is blocked and has no estimates.
Its reason code is `psis_k_not_evaluated` when a \\\hat k_m\\ is missing
or not finite, `psis_k_too_high` when all are finite but one is at or
above the cut-off, and `psis_weight_provenance_incomplete` when the
\\\hat k_m\\ pass but no program is named. pvstackr records the name of
the program but cannot check it. It also records in
`get_diagnostics(fit)$psis`, for each plausible value, the effective
sample size \\1/\sum_s (\tilde\omega_s^{(m)})^2\\ of the normalized
weights (`weight_ess_iid`) and their largest value
(`max_normalized_weight`); these values do not change the status.

The literature states the cut-off for \\\hat k\\ in two ways. Vehtari et
al. (2017, 1416) recommend warning the user when \\\hat k\\ exceeds 0.5,
although they observed good performance in practice for values up to
0.7, and above 0.7 advise considering alternatives to the
importance-sampling estimate. Vehtari et al. (2024, 5–6) consider
requiring \\\hat k \< 0.5\\ stricter than necessary and recommend a
warning when \\\hat k\\ exceeds \\\min(1 - 1/\log\_{10} S, 0.7)\\ for
\\S\\ draws, which is 0.67 for 1,000 draws and which they give as 0.7
for more than 2,000 draws. pvstackr uses a single cut-off and no level
at 0.5. To follow Vehtari et al. (2024) with fewer than 2,000 draws, set
`psis_k_threshold` to \\1 - 1/\log\_{10} S\\.

## Why only `stack_direct` intervals can be coverage-claimable

By pvstackr’s reporting rule, `coverage_claim_allowed` is `TRUE`, and
`interval_role` is `"coverage_barnard_rubin"`, only for the rows of a
`stack_direct` fit whose target was made with
`df_method = "barnard_rubin"` and a positive `df_complete`. The
intervals of every `per_pv` and `stack_psis` fit are descriptive.
[Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
lists the labels of every case and explains how to read and report them.
The rule is pvstackr’s own; it is not a result of the companion preprint
or of the sources of the formulas.

The rule follows the source of the variance within a plausible value.
All three methods add the variation between plausible values through
\\B\\, and the rule looks at \\\bar U\\, the average of the \\\hat
U_m\\, whose source differs between the methods; their \\\hat\beta_m\\,
and with them \\\bar\beta\\ and \\B\\, can differ as well. In a
`stack_direct` fit, \\\hat U_m\\ comes from the replicate weights: the
estimate is recomputed with each replicate weight, and the spread of the
replicate estimates around the full-sample estimate estimates its
sampling variance (Judkins 1990). PISA computes sampling variances in
this way, and the method accounts for the sampling of schools and of
students within them; for an estimate from plausible values, the
variation between them is added separately (OECD 2024, chs. 10–11).

In `per_pv` and `stack_psis`, \\\hat U_m\\ is the covariance of
posterior draws, so it depends on the model. It reflects the clustering
of the sample and the unequal selection probabilities only as far as the
model includes them; a single-level regression, such as the one pvstackr
fits, treats the students as independent. Weighting the likelihood by
the survey weights does not by itself give a design-based variance:
Williams and Savitsky (2021, 72–73) call the posterior of a
survey-weighted likelihood a pseudo-posterior and show that under a
complex sampling design its asymptotic distribution differs in scale and
shape from that of the weighted estimate, and that its credible
intervals can fall short of nominal coverage. pvstackr therefore labels
`per_pv` and `stack_psis` intervals descriptive, whichever rule for the
degrees of freedom they use. This is pvstackr’s reason, not a finding
about the preprint’s per-PV workflow: in the preprint’s simulation,
whose model has a school random intercept, that workflow’s model-based
intervals covered at 0.951 to 0.986 across strata and gradients (Lee et
al. 2026, sec. 5.2).

Among `stack_direct` fits, pvstackr also requires the Barnard–Rubin
degrees of freedom. Barnard and Rubin (1999, 948) note that the classic
degrees of freedom of Rubin (1987) are derived under the assumption that
the complete-data degrees of freedom \\\nu\_{\text{com}}\\ are infinite,
and that when \\\nu\_{\text{com}}\\ is small and little information is
missing they can be much larger than \\\nu\_{\text{com}}\\; their
adjusted degrees of freedom never exceed \\\nu\_{\text{com}}\\. With
`df_method = "barnard_rubin"` you state \\\nu\_{\text{com}}\\ of the
replicate-weight variance as `df_complete`; pvstackr does not choose or
check it ([The design-based target: BRR–Fay replicate weights and
Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.html#degrees-of-freedom)
gives the formulas). The labels follow `df_method` alone:
`df_complete = Inf` gives the classic degrees of freedom but keeps
`coverage_claim_allowed = TRUE`, and a fit with status `"warning"` keeps
its labels.

The rule rests on the source of the variance, not on a study of
coverage. The companion preprint’s simulation uses targets from
model-based fits, on simulated data without replicate weights, so it
does not evaluate a BRR–Fay target (Lee et al. 2026, sec. 5). Two parts
of a `stack_direct` fit say nothing about coverage either. The stacked
fixed-effect point identity concerns the point estimate only ([One
stacked fit and the Rubin point
estimate](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.html#scope)).
CCC gives the draws the mean and covariance of whatever target it
receives, and its diagnostics compare the stacked fit with that target.
The preprint calls CCC a moment calibration and treats the coverage of
intervals built in this way as an empirical question (Lee et al. 2026,
app. B.3); [Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md)
describes the calibration.

The label records how an interval was built; it does not certify
coverage. A coverage-claimable interval is built as in the standard
replicate-weight analysis, whose coverage is approximately nominal only
under that analysis’s assumptions: the replicate covariance estimates
the sampling covariance well, the regression is congenial to the model
that generated the plausible values, and `df_complete` suits the design.
pvstackr checks none of these. In the simulations of Judkins (1990, 233,
236), intervals for a regression coefficient based on replicate
variances, Fay’s method included, tended to cover less often than their
nominal level.

## Reading the Pareto \\\hat k\\ of a `stack_psis` fit

A Pareto \\\hat k_m\\ below the cut-off concerns the importance-sampling
step. Vehtari et al. (2024, 7) find PSIS estimates of an expectation
reliable, with low bias and low variance, when the Pareto \\\hat k\\ for
that expectation is below 0.7 at the numbers of draws typical in
practice; it can be larger than the \\\hat k\\ of the importance ratios
that `pareto_k` holds. When both are small, the weighted mean and
covariance of the stacked draws approximate the posterior mean and
covariance for plausible value \\m\\. A small \\\hat k_m\\ does not
change where that covariance comes from: when the weights reweight
toward the posteriors of the same model that a `per_pv` fit would use, a
`stack_psis` fit with small \\\hat k_m\\ approximates that `per_pv` fit,
with its model-based standard errors, and not the replicate-weight
target. Its intervals stay descriptive.

Whether the \\\hat k_m\\ are small depends on the data and the model.
With a common prior, the stacked fractional posterior \\q\\ is
proportional to the geometric mean of the \\M\\ posteriors for the
separate plausible values, \\\prod\_{m} \\p(\beta, \sigma)\\ \tilde
L_m(\beta, \sigma)\\^{1/M}\\, and each reweighting moves it toward one
of them. The companion preprint notes that when these posteriors are far
apart relative to their spread, the importance ratios have heavy tails
and \\\hat k\\ rises, and it treats the reliability of the reweighting
as something to check in each analysis (Lee et al. 2026, sec. 4.5 and
app. D.1). The preprint reports its reweighted stack only when \\\hat
k_m \< 0.7\\ for every plausible value, the cut-off that pvstackr uses
by default.

In the preprint’s application to PISA 2022 reading, the reweighting was
built in the same way, for the same model, in the United States and in
Korea. In the United States the \\\hat k\\ of all ten plausible values
were below 0.7, between 0.35 and 0.67; in Korea all ten were above it,
between 0.82 and 1.34 and above 1 for three of them, and the Korean
estimates were withheld (Lee et al. 2026, sec. 6.4 and table D.1). The
preprint explains the contrast by the larger and more strongly clustered
Korean sample, which makes each posterior narrow relative to the
distances between the posteriors, and calls this an interpretation of
the diagnostics rather than an established mechanism (Lee et al. 2026,
22–23). Its reweighting targets also differ from the per-PV posteriors
in the formula above: they place the survey weights in two stages, a
school factor on each school’s integrated likelihood and a student
factor inside the integral over the school intercept (its Eq. (6)),
whereas the stacked proposal carries the final student weight at the
student level. Each of its ratios therefore spans a change of weight
placement as well as a change of plausible value, a gap that the
preprint names as an additional difficulty for the proposal (Lee et al.
2026, sec. 4.5 and 6.4, app. D.1). Its \\\hat k\\ values are thus not
direct evidence about reweighting toward the per-PV posteriors of a
single-level model with the same weights. It also notes that with the
cut-off of Vehtari et al. (2024) for its 1,000 draws, 0.67, the largest
United States value, 0.674, would have failed as well (Lee et al. 2026,
app. D.1). The model of the application has a school random intercept,
which pvstackr 0.2.x does not fit.

[Comparing the three fitting
methods](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.html#cautions)
gives the practical cautions for reading `stack_psis` in a comparison,
and its section on [choosing a
method](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.html#choosing)
gives pvstackr’s recommendation, which uses `stack_psis` only as a check
on a stacked fit.

## One fit instead of \\M\\

`per_pv` needs \\M\\ Bayesian fits, one for each plausible value.
`stack_direct` and `stack_psis` need one, to the stacked data: \\M\\
copies of the \\N\\ rows of the data, in which row \\i\\ of copy \\m\\
has the outcome \\y_i^{(m)}\\. In a `stack_direct` fit that row has the
weight \\\tilde w_i/M\\, where \\\tilde w_i = W_i/\overline W\\ is the
final weight \\W_i\\ divided by the mean final weight \\\overline W\\
([One stacked fit and the Rubin point
estimate](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.html#stacked-objective)
derives the stacked objective from these weights). When `stack_psis`
stacks the data for your fitting function, the rows get the same weights
if you give `weight_col`, and the weight \\1/M\\ otherwise. With the ten
plausible values per domain of PISA 2022 (OECD 2024), the stacked
methods need one fit where `per_pv` needs ten.

The number of fits says how many samplers run and how many sets of
convergence diagnostics you check, not how long the analysis takes: the
stacked fit has \\M\\ times as many rows as a fit to one plausible
value, and `stack_direct` also needs its target, \\M(R + 1)\\ weighted
least-squares fits, which need no sampler. The companion preprint
likewise counts the gain as one fit, one set of convergence diagnostics
and one set of calibrated draws per model, and makes no claim about
processor time (Lee et al. 2026, sec. 7).
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
records the number of model fits of each method and, when you pass times
that you measured with its argument `timings`, those times; it does not
time the fits itself. [Real PISA data: access, memory and
reproducibility](https://joonho112.github.io/pvstackr/articles/a5-real-pisa-guidance.html#memory-runtime)
gives the size of these jobs for PISA 2022.

## Checks

The first chunk writes pvstackr’s check of the Pareto \\\hat k\\ values
in base R and applies it to made-up values: it illustrates the rule and
is not the function that pvstackr runs. The second reads the interval
labels that pvstackr gave the example fit.

### The Pareto \\\hat k\\ check

``` r

## pvstackr's check of the Pareto k-hat values of a stack_psis fit, written in
## base R to illustrate it; pvstackr applies it inside pv_fit_stack_psis().
## A value passes when it is finite and strictly below the cut-off.
k_passes <- function(k, threshold = 0.7) is.finite(k) & k < threshold

k <- c(0.20, 0.50, 0.69, 0.70, 0.95, NA)   # made-up k-hat values
data.frame(
  k            = k,
  passes_0.7   = k_passes(k),                               # the default cut-off
  passes_0.667 = k_passes(k, threshold = 1 - 1 / log10(1000))
)
#>      k passes_0.7 passes_0.667
#> 1 0.20       TRUE         TRUE
#> 2 0.50       TRUE         TRUE
#> 3 0.69       TRUE        FALSE
#> 4 0.70      FALSE        FALSE
#> 5 0.95      FALSE        FALSE
#> 6   NA      FALSE        FALSE
```

Each row applies the check to one made-up value. With the default
cut-off of 0.7, the values 0.20, 0.50 and 0.69 pass: pvstackr has no
separate level at 0.5. The value 0.70 fails because the check is strict,
and 0.95 and the missing value fail as well. A `stack_psis` fit has
estimates only when the values of all its plausible values pass and the
program that made the weights is named. A fit with a missing or
non-finite value is blocked with the reason code `psis_k_not_evaluated`,
and a fit whose values are finite but not all below the cut-off with
`psis_k_too_high`. The last column uses the cut-off that Vehtari et al.
(2024) recommend for 1,000 draws, \\1 - 1/\log\_{10} 1000 \approx
0.667\\, which you would set with
`pv_control(method = "stack_psis", psis_k_threshold = 1 - 1/log10(1000))`;
0.69 then fails too.

### The interval labels of the example fit

``` r

est <- get_estimates(fit)

unique(est$interval_role)
#> [1] "descriptive_classic_rubin"
unique(est$coverage_claim_allowed)
#> [1] FALSE
est[, c("term", "interval_role", "df_method", "coverage_claim_allowed")]
#>          term             interval_role df_method coverage_claim_allowed
#> 1 b_Intercept descriptive_classic_rubin   classic                  FALSE
#> 2         b_x descriptive_classic_rubin   classic                  FALSE
#> 3    b_female descriptive_classic_rubin   classic                  FALSE
```

All three rows have `interval_role = "descriptive_classic_rubin"`,
`df_method = "classic"` and `coverage_claim_allowed = FALSE`. The fit is
a `stack_direct` fit, but its target was made with the classic degrees
of freedom, so by pvstackr’s rule its intervals are descriptive; a
target made with `df_method = "barnard_rubin"` and a positive
`df_complete` would give `coverage_barnard_rubin` and `TRUE`. The labels
do not depend on the draws: the fit copies them from its target.

[Comparing the three fitting
methods](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.md)
lines up an example fit of each method with
[`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
and reads their labels side by side, and [Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
explains how to report each kind of interval.

The numbers above were computed in this session:

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] pvstackr_0.2.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.60         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  
#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.3.0       fs_2.1.0
```

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

OECD. 2024. *PISA 2022 Technical Report*. PISA. OECD Publishing.
<https://doi.org/10.1787/01820d6d-en>.

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC.”
*Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

Vehtari, Aki, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah
Gabry. 2024. “Pareto Smoothed Importance Sampling.” *Journal of Machine
Learning Research* 25 (72): 1–58.
<https://www.jmlr.org/papers/v25/19-556.html>.

Williams, Matthew R., and Terrance D. Savitsky. 2021. “Uncertainty
Estimation for Pseudo-Bayesian Inference Under Complex Sampling.”
*International Statistical Review* 89 (1): 72–107.
<https://doi.org/10.1111/insr.12376>.
