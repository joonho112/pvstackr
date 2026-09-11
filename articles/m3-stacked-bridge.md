# One stacked fit and the Rubin point estimate

Abstract

A `stack_direct` fit is one Bayesian fit to the stacked data: the
plausible values in one data set, with each row weighted by the
student’s normalized survey weight divided by the number of plausible
values. This article writes the log likelihood of that fit as the
average of the survey-weighted log-likelihoods of the plausible values
and its posterior as the stacked fractional posterior of the companion
preprint. It states the stacked fixed-effect point identity (Theorem 4.1
of the preprint) and its conditions, under which the fixed-effect
estimate of the stacked fit equals the Rubin point estimate, shows why
the identity holds exactly in the model that pvstackr fits, and lists
what it does not cover. A check on the example fit shows that the
reported estimates are those of the target.

``` r

library(pvstackr)
```

With \\M\\ plausible values, the standard analysis fits its model once
for each plausible value and combines the \\M\\ results with Rubin’s
rules (Rubin 1987). Its point estimate is the average of the \\M\\
estimates, \\\bar\beta = \frac{1}{M}\sum\_{m=1}^{M} \hat\beta_m\\. The
`stack_direct` method of pvstackr makes one Bayesian fit instead, to all
\\M\\ plausible values at once. Under conditions stated in the companion
preprint (Lee et al. 2026), the fixed-effect estimate of this weighted
fit to the stacked data equals \\\bar\beta\\.

In pvstackr, \\\hat\beta_m\\ is the weighted least-squares estimate for
plausible value \\m\\, and \\\bar\beta\\ is part of the design-based
target that
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
computes together with its covariance matrix \\T\_{\text{MI}}\\ ([The
design-based target: BRR–Fay replicate weights and Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.md)).
The equality concerns the point estimate only. The standard errors of a
`stack_direct` fit come from \\T\_{\text{MI}}\\, and the Cholesky
calibration correction (CCC) gives the draws of the stacked fit the mean
and covariance of the target ([Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md)).

The code in this article reads the example fit that ships with pvstackr
and its target; it fits no model:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # the BRR–Fay target of the fit, made by pv_brr_target()
```

The example fit and its synthetic data are described in [Getting started
with
pvstackr](https://joonho112.github.io/pvstackr/articles/a1-getting-started.html#example-fit).

## The stacked data

Write \\y_i^{(m)}\\ for plausible value \\m\\ of student \\i\\, \\i = 1,
\dots, N\\, and \\x_i\\ for the row of the \\N \times p\\ design matrix
\\X\\ that holds the intercept and the covariates of student \\i\\. The
stacked data put \\M\\ copies of the data one below the other, \\N \cdot
M\\ rows in all: row \\i\\ of copy \\m\\ has the outcome \\y_i^{(m)}\\
and the covariates of student \\i\\. The covariates are the same in
every copy; only the outcome changes.

Each row has a weight. Let \\W_i\\ be the final survey weight of student
\\i\\, \\\overline W\\ the mean of the final weights, and \\\tilde w_i =
W_i/\overline W\\ the normalized weight, which has mean 1. Row \\i\\ of
every copy has the weight

\\ \frac{\tilde w_i}{M} . \\

The \\M\\ copies of a student together carry the weight \\\tilde w_i\\,
so the weights of each copy sum to \\N/M\\ and those of all rows to
\\N\\. The factor \\1/M\\ keeps the fit from counting each student \\M\\
times, and \\\tilde w_i\\ brings the survey weights into the fit. With
the weight \\1/M\\ alone, the stacked fit would estimate the
coefficients of an unweighted regression, which in general differ from
the survey-weighted \\\bar\beta\\.

pvstackr builds the stacked data itself. `stack_direct` takes \\W_i\\
from the column `weight_col` of the target, stores \\\tilde w_i/M\\ in
the column `.pvstackr_weight`, and passes it to the fitting engine in
the formula term `weights(.pvstackr_weight)`. The brms engine bundled
with pvstackr multiplies the log-likelihood contribution of each row by
its weight; a fitting function of your own has to use the weight in the
same way ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#own-engine)).
The row weight is \\1/M\\ alone in one case only: when `stack_psis` fits
the stacked model with your `fit_function` and no `weight_col` is given.

## The stacked objective and the stacked fractional posterior

pvstackr fits a linear regression with normal errors, \\y_i^{(m)} =
x_i^\top \beta + \varepsilon_i^{(m)}\\ with \\\varepsilon_i^{(m)} \sim
N(0, \sigma^2)\\, where \\\beta\\ holds the \\p\\ fixed-effect
coefficients and \\\sigma\\ is the residual standard deviation
([Plausible values, survey weights and
notation](https://joonho112.github.io/pvstackr/articles/m1-foundations-and-notation.html#model)
describes the model). For plausible value \\m\\, its survey-weighted
likelihood is

\\ \tilde L_m(\beta, \sigma) = \prod\_{i=1}^{N} \phi\big(y_i^{(m)};\\
x_i^\top \beta,\\ \sigma^2\big)^{\tilde w_i}, \\

where \\\phi(\cdot\\; \mu, \sigma^2)\\ is the normal density with mean
\\\mu\\ and variance \\\sigma^2\\: each student’s log-likelihood
contribution is multiplied by \\\tilde w_i\\. This is the weighted
likelihood of the preprint (its Eq. (3), whose weights are also
normalized to mean 1) for a model without the school random intercept of
the preprint’s analysis model, a model that pvstackr 0.2.x does not fit.

Because the engine multiplies the log-likelihood of each stacked row by
\\\tilde w_i/M\\, the log likelihood of the stacked fit is

\\ \sum\_{m=1}^{M} \sum\_{i=1}^{N} \frac{\tilde w_i}{M} \log
\phi\big(y_i^{(m)};\\ x_i^\top \beta,\\ \sigma^2\big) = \frac{1}{M}
\sum\_{m=1}^{M} \log \tilde L_m(\beta, \sigma) , \\

the average of the \\M\\ survey-weighted log-likelihoods. The preprint
calls this average the stacked objective. With a prior \\p(\beta,
\sigma)\\, the posterior that the stacked fit samples is the stacked
fractional posterior of the preprint (its Eq. (4)),

\\ q(\beta, \sigma) \propto p(\beta, \sigma) \prod\_{m=1}^{M} \tilde
L_m(\beta, \sigma)^{1/M} . \\

The exponent \\1/M\\, the fractional weight, makes the \\M\\ likelihoods
enter as a geometric average rather than as \\M\\ reuses of the same
students, so that the posterior carries one sample’s worth of
information (Lee et al. 2026, sec. 4.3). The preprint places this
construction among the fractional and composite likelihoods (Lee et al.
2026, sec. 3.2).

### The prior on the fixed effects

With the brms engine bundled with pvstackr
(`pv_control(backend = "brms")`) and the default `prior = NULL`, every
fixed effect, the intercept included, has the flat prior that brms gives
by default to coefficients of its class `"b"`. The reason is that
`stack_direct` adds the columns of the design matrix, the intercept
column included, to the stacked data as ordinary covariates and fits the
formula
`.pvstackr_y | weights(.pvstackr_weight) ~ 0 + pvstackrMM001 + ...`, so
all coefficients belong to class `"b"`. The Student-\\t\\ prior that
brms gives to an intercept it builds from the formula does not arise.
The residual standard deviation \\\sigma\\ keeps brms’s default prior.

An explicit `prior` gives the fit the status `"warning"`, with the
reason code `explicit_prior_warning`, even when the prior is flat,
because pvstackr does not check whether a prior is flat;
[`?pv_fit_direct`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
lists the priors that `stack_direct` accepts. A fitting function of your
own receives `prior`, and the priors are then those of your engine. The
flat prior on \\\beta\\ is part of condition (R5) of the point identity
([the conditions](#conditions)).

## The stacked fixed-effect point identity

The stacked fixed-effect point identity, Theorem 4.1 of the preprint
(Lee et al. 2026, sec. 4.3), concerns the preprint’s analysis model, a
linear model with a school random intercept ([Plausible values, survey
weights and
notation](https://joonho112.github.io/pvstackr/articles/m1-foundations-and-notation.html#model)
shows it), at a fixed working covariance matrix \\V\\ of the outcome
vector. Write \\\log \tilde L_m(\beta; V)\\ for the weighted
log-likelihood of plausible value \\m\\ at this \\V\\, and
\\\hat\beta_m\\ for the generalized least-squares estimate that
maximizes it. Under conditions (R1)–(R5) below, the maximizer of the
stacked objective in \\\beta\\ is the average of these estimates,

\\ \hat\beta^{q} = \arg\max\_{\beta}\\ \frac{1}{M} \sum\_{m=1}^{M} \log
\tilde L_m(\beta; V) = \frac{1}{M} \sum\_{m=1}^{M} \hat\beta_m =
\bar\beta , \\

exactly. Because the stacked objective is then quadratic in \\\beta\\,
the marginal distribution of \\\beta\\ under the stacked fractional
posterior \\q\\ is normal, so its mode and its mean coincide.

### The conditions

The preprint states the conditions as follows (Lee et al. 2026, 11):

- (R1) “all \\M\\ PV analyses share the same analytic rows and cluster
  membership”;
- (R2) “the design matrix \\X\\ is common across PVs”;
- (R3) “the marginal working covariance \\V\\ of the outcome vector is
  common across PVs and held fixed”;
- (R4) “the PV weights are equal at \\1/M\\”;
- (R5) “the prior on \\\beta\\ is flat and \\X^\top V^{-1} X\\ is
  invertible”.

PV stands for plausible value. Appendix A.3 of the preprint restates the
theorem as Theorem A.1. There (R3) reads “the working covariance \\V\\
is common across PVs and held fixed”, and (R5) adds that any common
symmetric positive-definite working precision \\P\\ may take the place
of \\V^{-1}\\. The setting is the preprint’s two-level Gaussian model
with its variance components held fixed. The sampling weights enter by
scaling each student’s log-likelihood contribution, and the design
matrix is the same for all plausible values “because the covariates,
analytic rows, and cluster membership do not vary with \\m\\” (Appendix
A.2 of the preprint). A school mean of a covariate is a column of \\X\\,
so under (R2) it must be the same for every plausible value.

### Why the identity holds

The preprint calls the proof “one observation” (Lee et al. 2026, 11).
Let \\y^{(m)} = (y_1^{(m)}, \dots, y_N^{(m)})^\top\\ be the vector of
plausible value \\m\\ and \\P = V^{-1}\\. At a fixed \\V\\, the
log-likelihood of plausible value \\m\\ is, apart from terms without
\\\beta\\, \\-\frac12 (y^{(m)} - X\beta)^\top P\\ (y^{(m)} - X\beta)\\,
a quadratic in \\\beta\\ whose score

\\ X^\top P \big(y^{(m)} - X\beta\big) \\

is linear in \\y^{(m)}\\. Setting it to zero gives \\\hat\beta_m =
(X^\top P X)^{-1} X^\top P\\ y^{(m)}\\. By (R1)–(R4), the score of the
stacked objective is the average of the \\M\\ scores, \\X^\top P (\bar
y - X\beta)\\ with \\\bar y = \frac{1}{M}\sum\_{m=1}^{M} y^{(m)}\\, and
the flat prior of (R5) adds nothing to it. The stacked estimate is
therefore \\(X^\top P X)^{-1} X^\top P\\ \bar y\\, the average of the
\\\hat\beta_m\\, because the map from \\y^{(m)}\\ to \\\hat\beta_m\\ is
linear; the invertibility in (R5) makes it unique. If \\X\\ or \\P\\
differed between plausible values, the stacked estimate would be a
matrix-weighted average of the \\\hat\beta_m\\, which in general is not
\\\bar\beta\\. The full proof is in Appendix A.3 of the preprint, which
also reports a numerical check: over a grid of 24 test configurations,
the largest absolute difference between the stacked solution and the
average of the per-PV solutions is \\4.4 \times 10^{-16}\\.

### The identity in the model that pvstackr fits

In the survey-weighted linear model of pvstackr the conditions hold as
follows.

- (R1): every copy of the stacked data holds the same \\N\\ rows, and
  the model has no random effects, so there is no cluster membership to
  share.
- (R2): pvstackr copies every covariate unchanged into each copy, so
  \\X\\ is the same for all plausible values. This includes a school
  mean that you compute and store as a column; pvstackr computes no
  school means itself.
- (R3): at a given \\\sigma\\, the log-likelihood \\\log \tilde
  L_m(\beta, \sigma)\\ is, apart from terms without \\\beta\\,
  \\-\frac12 (y^{(m)} - X\beta)^\top P\\ (y^{(m)} - X\beta)\\ with \\P =
  \tilde W/\sigma^2\\ and \\\tilde W = \operatorname{diag}(\tilde w_1,
  \dots, \tilde w_N)\\. Every copy has the same weights, so \\P\\ is the
  same for all plausible values.
- (R4): every row of every copy has the factor \\1/M\\ in its weight.
- (R5): with the bundled engine and `prior = NULL`, the prior on
  \\\beta\\ is flat. \\X^\top \tilde W X\\ must be invertible; the
  target needs this as well, and
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  stops with an error when a coefficient cannot be estimated.

The theorem holds the working covariance fixed, while the stacked fit
estimates \\\sigma\\ together with \\\beta\\. In this model that does
not change the result: \\\sigma\\ only scales \\P\\, so the \\\beta\\
that maximizes the stacked objective is the same for every \\\sigma\\.
It is the weighted least-squares estimate from the stacked data,

\\ \hat\beta^{q} = \big(X^\top \tilde W X\big)^{-1} X^\top \tilde W\\
\bar y = \frac{1}{M} \sum\_{m=1}^{M} \hat\beta_m = \bar\beta , \\

where \\\hat\beta_m = (X^\top \tilde W X)^{-1} X^\top \tilde W\\
y^{(m)}\\ is the weighted least-squares estimate of the target for
plausible value \\m\\ (the scale of the weights does not change it). For
the example data the two sides differ by rounding error only, as the
first check below shows. Under the flat prior, the stacked fractional
posterior of \\\beta\\ given \\\sigma\\ is normal with mean
\\\bar\beta\\ for every \\\sigma\\, so the posterior mean of \\\beta\\
is \\\bar\beta\\. The mean of the draws of a stacked fit made with the
flat prior therefore differs from \\\bar\beta\\ only by Monte Carlo
error.

## What the identity does not cover

Appendix A.4 of the preprint lists four limits of the theorem (Lee et
al. 2026, 36):

1.  It is a point identity only. It equates maximizers and, under the
    flat prior, posterior means; it “makes no statement about posterior
    variances, interval widths, or frequentist coverage”.
2.  It covers the fixed effects only. The per-PV score for the variance
    components is not linear in the plausible values, so the stacked
    estimates of the variance components can differ from pooled per-PV
    estimates.
3.  It is “not a production theorem”. In the preprint’s production fits
    the working covariance is estimated together with \\\beta\\, which
    violates (R3); the priors are weakly informative, which violates the
    flat-prior part of (R5); and Monte Carlo error enters. There, the
    agreement between the stacked center and the pooled center “is an
    empirical quantity to be measured and reported, never asserted as a
    corollary of the theorem”.
4.  It does not carry over to reweighting: the importance-reweighted
    per-PV centers of the preprint’s reweighted-stack variant are
    approximations, and the identity says nothing about them.

For pvstackr these limits mean the following.

The identity gives no variance. The covariance of the stacked draws does
not contain the variation between plausible values, \\B\\, which exists
only across separate analyses of the plausible values (Lee et al. 2026,
sec. 4.4). The standard errors, degrees of freedom and interval labels
of a `stack_direct` fit therefore come from the target, and CCC gives
the fixed-effect draws the target’s mean \\\bar\beta\\ and covariance
\\T\_{\text{MI}}\\ ([Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.html#map)).
The identity alone does not make one stacked fit reproduce the inference
of \\M\\ fits: the preprint’s calibrated-stack workflow, like
`stack_direct`, combines the stacked fit with a target computed from the
\\M\\ plausible values and a calibration to that target. Which intervals
are labeled coverage-claimable is set by pvstackr’s reporting rule, not
by the identity, and the label does not certify coverage ([Reading and
reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)).

The stacked fit also estimates \\\sigma\\, and its estimate can differ
from those of the separate fits for each plausible value; pvstackr
neither calibrates nor reports it.

In pvstackr’s model, estimating \\\sigma\\ does not move the maximizer
in \\\beta\\, as shown above. Of the preprint’s reasons why the identity
is only approximate in production, two remain: a prior that is not flat,
and Monte Carlo error. With the bundled engine and `prior = NULL`, only
Monte Carlo error remains; an explicit prior, or a fitting function of
your own that uses other priors or ignores the weights, can move the
mean of the draws further. pvstackr measures the difference in every
`stack_direct` fit. `delta_c_max` is the largest difference, in target
standard errors, between a target estimate and the mean of the stacked
draws before calibration; by pvstackr’s thresholds a value of 0.01 or
more gives the fit the status `"warning"`, and 0.05 or more the status
`"blocked"`. Monte Carlo error alone can reach the first of these
thresholds ([Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.html#diagnostics)
gives its size). The estimate that pvstackr reports is \\\bar\beta\\
itself, copied from the target ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)).
The preprint’s workflow does the same: “the reported center is set to
the target directly” (Lee et al. 2026, 11).

`stack_psis` reweights the stacked draws with importance weights that
you supply, to approximate the posterior for each plausible value; the
identity says nothing about these reweighted centers ([Choosing a method
and reading its
intervals](https://joonho112.github.io/pvstackr/articles/m5-methods-and-coverage.html#psis)).

## Checks on the example fit

The first chunk below checks the point identity on the example data. The
second compares the estimate table of the example fit with its target;
it shows pvstackr’s reporting rule, not the point identity. With
`center = "target"`, which `stack_direct` requires, the estimate table
of every `stack_direct` fit holds the target’s \\\bar\beta\\ (field
`beta_bar`) as `estimate` and the square roots of the diagonal of
\\T\_{\text{MI}}\\ (field `T_MI`) as `se`, whatever the draws of the
stacked fit are ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)).

The identity itself can be checked without a fit. The first chunk stacks
the example data, weights each row by \\\tilde w_i/M\\ and computes the
weighted least-squares estimate from the stacked data:

``` r

pisa_tiny <- read.csv(system.file("extdata", "pisa_tiny.csv", package = "pvstackr"))
X <- model.matrix(~ x + female, data = pisa_tiny)     # the design matrix
w <- pisa_tiny$W_FSTUWT / mean(pisa_tiny$W_FSTUWT)    # normalized weights
M <- 2
stacked_y <- c(pisa_tiny$PV1READ, pisa_tiny$PV2READ)  # one copy per plausible value
stacked_X <- rbind(X, X)
b_stacked <- lm.wfit(stacked_X, stacked_y, rep(w / M, M))$coefficients
max(abs(unname(b_stacked) - unname(tg$beta_bar)))     # rounding error only
#> [1] 5.684342e-14

b_equal <- lm.wfit(stacked_X, stacked_y, rep(1 / M, M * nrow(X)))$coefficients
max(abs(unname(b_equal) - unname(tg$beta_bar)))       # rows weighted 1/M only
#> [1] 0.1071195
```

With the rows weighted by \\\tilde w_i/M\\ the stacked estimate equals
\\\bar\beta\\ up to rounding error, about \\10^{-13}\\. With every row
weighted \\1/M\\, that is, without the survey weights, it differs by
about 0.1 score points, because it then estimates the unweighted
regression.

The next chunk compares the estimate table of the fit with the target:

``` r

est <- get_estimates(fit)

all.equal(est$estimate, unname(tg$beta_bar))   # TRUE: the estimate is copied from beta_bar
#> [1] TRUE
data.frame(
  term             = est$term,
  stacked_estimate = est$estimate,
  rubin_beta_bar   = unname(tg$beta_bar)
)
#>          term stacked_estimate rubin_beta_bar
#> 1 b_Intercept       457.894088     457.894088
#> 2         b_x        46.883361      46.883361
#> 3    b_female         2.143702       2.143702
```

The comparison returns `TRUE`, and the two columns hold the same
numbers: about 457.89 for the intercept, 46.88 for `b_x` and 2.14 for
`b_female`. The column `stacked_estimate` is the `estimate` column of
the fit, which pvstackr copied from the target; it is not the mean of
the stacked draws. The difference between that mean and \\\bar\beta\\ is
what the point identity is about, and pvstackr records it as
`delta_c_max` (see [what the identity does not cover](#scope)). For the
example fit it is of the order of \\10^{-14}\\, because the draws were
built around the target ([Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#target-and-diagnostics)
shows the value).

[Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md)
describes how the stacked draws are given the mean and covariance of the
target, and [Choosing a method and reading its
intervals](https://joonho112.github.io/pvstackr/articles/m5-methods-and-coverage.md)
compares `stack_direct` with `per_pv` and `stack_psis`.

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

Lee, JoonHo, Matthew R. Williams, and Terrance D. Savitsky. 2026. *One
Markov Chain Monte Carlo Fit for Many Plausible Values: A Calibrated
Stacked Posterior Workflow for Bayesian Multilevel Models of Large-Scale
Assessment Data*. Zenodo preprint, version 1.
<https://doi.org/10.5281/zenodo.22407935>.

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.
