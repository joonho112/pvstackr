# Calibrating the fixed-effect draws (CCC)

Abstract

In a `stack_direct` fit, the Cholesky calibration correction (CCC)
transforms the fixed-effect draws of one Bayesian fit to the stacked
plausible-value data so that their mean and covariance equal the
estimates and the covariance matrix of the design-based target. This
article gives the map, its exact moment properties and their source, the
algorithm that pvstackr uses and when it stops with an error, and the
diagnostics `delta_c_max`, `delta_c_rel` and `kappa_A` with pvstackr’s
thresholds. A check on simulated draws confirms the moment properties;
neither the check nor the calibration says anything about the coverage
of intervals.

``` r

library(pvstackr)
```

A `stack_direct` fit rests on two calculations. The target, computed by
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
fits a survey-weighted linear regression by weighted least squares to
each of the \\M\\ plausible values, estimates the covariance of each
estimate from the BRR–Fay replicate weights, and combines the results
with Rubin’s rules into the estimate \\\bar\beta\\ of the \\p\\ fixed
effects and its total covariance matrix \\T\_{\text{MI}} = \bar U + (1 +
1/M)B\\. Here \\\bar U\\ is the average of the \\M\\ replicate
covariances and \\B\\ the covariance between the \\M\\ estimates ([The
design-based target: BRR–Fay replicate weights and Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.html#rubin)).
The stacked fit is one Bayesian fit of the same regression to the
stacked data: \\M\\ copies of the \\N\\ rows, one for each plausible
value, in which row \\i\\ has the weight \\\tilde w_i/M\\, where
\\\tilde w_i = W_i/\overline W\\ is the final weight \\W_i\\ divided by
the mean final weight \\\overline W\\.

The Cholesky calibration correction (CCC) connects the two. It
transforms the fixed-effect draws of the stacked fit so that their mean
is \\\bar\beta\\ and their covariance is \\T\_{\text{MI}}\\, and before
it does so it records how far the stacked fit was from the target. The
code below reads the example fit that ships with pvstackr and its
target:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # the target of the fit, made by pv_brr_target()

is.null(get_draws(fit))
#> [1] TRUE
```

The example fit and its synthetic data are described in [Getting started
with
pvstackr](https://joonho112.github.io/pvstackr/articles/a1-getting-started.html#example-fit).
The fit was saved without its calibrated draws, so [the check of the
map](#checks) uses simulated draws, and the diagnostics that the
calibration recorded are read from the fit.

## Why the stacked draws are calibrated

The stacked fit is built to agree with the target in its point estimate.
Its log likelihood is the average of the \\M\\ survey-weighted
log-likelihoods, one for each plausible value, and the stacked
fixed-effect point identity, Theorem 4.1 of the companion preprint (Lee
et al. 2026), gives conditions under which the value of \\\beta\\ that
maximizes this average equals \\\bar\beta\\. In pvstackr’s model the
weighted least-squares estimate from the stacked data equals
\\\bar\beta\\ exactly, so with a flat prior on the fixed effects the
mean of the draws differs from \\\bar\beta\\ only by Monte Carlo error
([One stacked fit and the Rubin point
estimate](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.html#point-identity)
explains the identity and its conditions).

The identity concerns the point estimate only and gives no variance. The
between-imputation covariance \\B\\ is the covariance of estimates from
separate analyses of the \\M\\ plausible values and exists only across
such analyses, so the draws of one stacked fit do not contain it. The
preprint notes that the covariance of the stacked fixed-effect draws
approximates a one-sample quantity that corresponds to the
within-imputation part \\\bar U\\ only, and so lacks the term \\(1 +
1/M)B\\ of \\T\_{\text{MI}}\\ (Lee et al. 2026, 11–12, 40). In pvstackr
the target’s \\\bar U\\ comes from the replicate weights, whereas the
spread of the stacked draws comes from the model: the likelihood of the
stacked fit carries the survey weights, but not the clustering of the
sample or the replicate weights.

pvstackr therefore takes the numbers it reports from the target ([The
full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)),
and CCC transforms the draws to agree with these numbers. The calibrated
draws, returned by `get_draws(fit)`, have the mean \\\bar\beta\\ and the
covariance \\T\_{\text{MI}}\\, so a linear combination of the
coefficients computed from them has the target’s estimate and variance.
They do not carry the target’s degrees of freedom: quantiles of the
draws give roughly a normal-reference interval, which is narrower than
the \\t\\ interval of the estimate table when the degrees of freedom are
small (see [the calibration map](#map)). Before the draws are moved, the
calibration compares the stacked fit with the target (see [the
diagnostics](#diagnostics)).

## The calibration map

Let \\\beta_1, \dots, \beta_S\\ be the \\S\\ draws of the \\p\\ fixed
effects from the stacked fit, with mean and covariance

\\ \bar\beta^{\text{raw}} = \frac{1}{S}\sum\_{s=1}^{S} \beta_s, \qquad
\Sigma\_{\text{raw}} = \frac{1}{S-1}\sum\_{s=1}^{S} \big(\beta_s -
\bar\beta^{\text{raw}}\big)\big(\beta_s -
\bar\beta^{\text{raw}}\big)^\top , \\

and let \\L\_{\text{raw}}\\ and \\L\_{\text{tgt}}\\ be the
lower-triangular Cholesky factors of \\\Sigma\_{\text{raw}}\\ and
\\T\_{\text{MI}}\\, so that \\L\_{\text{raw}} L\_{\text{raw}}^\top =
\Sigma\_{\text{raw}}\\ and \\L\_{\text{tgt}} L\_{\text{tgt}}^\top =
T\_{\text{MI}}\\. CCC maps each draw to

\\ \beta^{\text{cal}}\_s = \bar\beta + A\big(\beta_s -
\bar\beta^{\text{raw}}\big), \qquad A = L\_{\text{tgt}}
L\_{\text{raw}}^{-1} . \\

The preprint gives this map as its Eq. (5) (Lee et al. 2026, 12); in
pvstackr, \\\bar\beta\\ and \\T\_{\text{MI}}\\ are those of the target
from
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
The map centers each draw at the mean of the stacked draws, multiplies
it by the calibration matrix \\A\\, and adds the target estimate
\\\bar\beta\\. Multiplying the centered draws by
\\L\_{\text{raw}}^{-1}\\ gives them the identity matrix as covariance,
and multiplying them by \\L\_{\text{tgt}}\\ then gives them the
covariance \\T\_{\text{MI}}\\.

The calibrated draws have three properties, which make up exact moment
calibration, Proposition B.1 of the preprint (Lee et al. 2026, 37):

- their mean is \\\bar\beta\\;
- their covariance, computed like \\\Sigma\_{\text{raw}}\\, is
  \\T\_{\text{MI}}\\;
- the draws of parameters other than the fixed effects are unchanged.

The first holds because the centered draws have mean zero and \\A\\ maps
zero to zero. The second follows from the way a covariance matrix
changes under a linear map:

\\ A\\\Sigma\_{\text{raw}}\\A^\top = L\_{\text{tgt}}
L\_{\text{raw}}^{-1} \big(L\_{\text{raw}} L\_{\text{raw}}^\top\big)
L\_{\text{raw}}^{-\top} L\_{\text{tgt}}^\top = L\_{\text{tgt}}
L\_{\text{tgt}}^\top = T\_{\text{MI}} . \\

The third holds because the map acts on the fixed-effect columns only.
In pvstackr’s model the only other parameter of the stacked fit is the
residual standard deviation \\\sigma\\, whose draws pass through
unchanged; pvstackr neither calibrates nor reports it, and
`get_draws(fit)` returns the calibrated fixed-effect draws only.

These properties are algebraic. They hold for any draws whose covariance
matrix is positive definite, apart from the rounding error of
floating-point arithmetic, and they do not depend on how well the
stacked fit was sampled: because the map uses the mean and covariance of
the draws at hand, Monte Carlo error changes \\A\\ but not the mean and
covariance of the calibrated draws (Lee et al. 2026, 37).

The preprint describes CCC as the scale-and-shape correction that
Williams and Savitsky (2021) developed for pseudo-posteriors under
complex sampling, specialized to the fixed-effect block of the model and
to a covariance target from the combining rules (Lee et al. 2026, 7,
13).

CCC sets the mean and the covariance of the draws and nothing else. The
calibrated draws are an affine transformation of the stacked draws, so
they are normal only when the stacked draws are, and their other
features follow from those of the stacked draws through the same
transformation. The interval of a `stack_direct` fit does not use the
draws: it is the estimate plus or minus a \\t\\ quantile times the
standard error, all taken from the target ([Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)).
Quantiles of the calibrated draws therefore do not reproduce this
interval. In the example, whose classic degrees of freedom are between
1.0 and 1.4, the 95% \\t\\ intervals are three to six times as wide as
normal-reference intervals with the same standard errors; with \\M =
10\\ the classic degrees of freedom are at least 9, so the \\t\\
quantile is at most 2.26 against 1.96; Barnard–Rubin degrees of freedom
can be smaller, and the gap then larger.

## The center of the calibrated draws

The map above centers the calibrated draws at the target estimate
\\\bar\beta\\. This is the setting `center = "target"` of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
the default, and a `stack_direct` fit requires it:

``` r

fit$control$center   # "target": the calibrated draws are centered at beta_bar
#> [1] "target"
```

With this setting the calibrated draws have the mean and covariance of
the target, whose numbers the estimate table copies ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)
explains why pvstackr requires this setting), so the table and the draws
agree.
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
also accepts `center = "posterior"`, which would keep the mean of the
stacked draws, \\\beta^{\text{cal}}\_s = \bar\beta^{\text{raw}} +
A(\beta_s - \bar\beta^{\text{raw}})\\, and calibrate only the
covariance; no fitting function uses it ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)).

## The algorithm

CCC is deterministic: pvstackr computes it with plain Cholesky
factorizations, a triangular solve and matrix products, without random
numbers, iterations or tuning, so the same draws and target always give
the same calibrated draws. A fit that pvstackr’s check of the sampler
diagnostics blocks is not calibrated. For the other fits, pvstackr
proceeds as follows:

1.  It takes the fixed-effect columns of the draws of the stacked fit
    and computes their mean \\\bar\beta^{\text{raw}}\\ and covariance
    \\\Sigma\_{\text{raw}}\\ with
    [`colMeans()`](https://rdrr.io/r/base/colSums.html) and
    [`cov()`](https://rdrr.io/r/stats/cor.html).
2.  It checks the target: \\T\_{\text{MI}}\\ must be symmetric and
    positive definite, its rows and columns must name the same
    coefficients as the draws, in the same order, and the target must
    come from
    [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
    It also stops with an error when \\T\_{\text{MI}}\\ equals
    \\\Sigma\_{\text{raw}}\\ or is proportional to it.
3.  It factors \\\Sigma\_{\text{raw}}\\ and \\T\_{\text{MI}}\\ with
    [`chol()`](https://rdrr.io/r/base/chol.html) from base R, without
    pivoting, after making each matrix exactly symmetric as \\(\Sigma +
    \Sigma^\top)/2\\. [`chol()`](https://rdrr.io/r/base/chol.html)
    returns the upper-triangular factor \\U\\ with \\U^\top U =
    \Sigma\\, and pvstackr uses its transpose, \\L = U^\top\\.
4.  It computes \\L\_{\text{raw}}^{-1}\\ by forward substitution, with
    [`forwardsolve()`](https://rdrr.io/r/base/backsolve.html) applied to
    the identity matrix, and forms \\A = L\_{\text{tgt}}
    L\_{\text{raw}}^{-1}\\, a lower-triangular \\p \times p\\ matrix.
5.  It centers every draw at \\\bar\beta^{\text{raw}}\\, multiplies it
    by \\A\\ and adds \\\bar\beta\\, and records [the
    diagnostics](#diagnostics).

A positive definite matrix has exactly one lower-triangular Cholesky
factor with a positive diagonal, so \\A\\ is well defined. Other
matrices \\G\\ with \\G\\\Sigma\_{\text{raw}}\\G^\top = T\_{\text{MI}}\\
exist, for example matrices built from symmetric square roots; used in
place of \\A\\, they give the calibrated draws the same mean and
covariance but move each draw differently. The preprint’s algorithm,
like pvstackr, uses the lower-triangular factors of both matrices (Lee
et al. 2026, 37).

Both factorizations need positive definite matrices, and pvstackr
applies no automatic Tikhonov or near-PD repair: it adds no ridge (a
constant on the diagonal) and does not replace a matrix by a nearby
positive definite one. When \\\Sigma\_{\text{raw}}\\ or
\\T\_{\text{MI}}\\ is not positive definite, pvstackr stops with an
error and returns no fit. For \\\Sigma\_{\text{raw}}\\ the message is
“pvstackr: `Sigma_raw` must be positive definite for CCC calibration.”;
`Sigma_raw` and `Sigma_target` are the names of \\\Sigma\_{\text{raw}}\\
and \\T\_{\text{MI}}\\ in the calibration record of a fit (“Calibration
record” in
[`?pvstackr_object_contracts`](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)).
The calibration in the preprint also halts with an error when its
preconditions, positive definite matrices among them, fail; its
implementation keeps a nearest positive-semidefinite repair as an
explicitly logged fallback, which pvstackr does not have (Lee et al.
2026, 38).

\\\Sigma\_{\text{raw}}\\ is not positive definite when there are fewer
than \\p + 1\\ draws, or when the draws of one fixed effect are an exact
linear combination of the draws of others, as for a fixed effect that
the model does not identify. The remedy is to change the stacked fit,
for example with more draws or with fewer or less collinear fixed
effects. A target covariance that is not positive definite already stops
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
with an error when the target is made, and
`pv_control(allow_target_nearpd = TRUE)`, which would allow a repair of
the target, stops with an error as well.

A covariance matrix that is positive definite but close to singular
passes the factorization, and so, because of rounding, can a matrix that
is singular in exact arithmetic. When one of the two matrices is much
closer to singular than the other, the condition number \\\kappa_A\\ of
\\A\\, stored as `kappa_A`, is large and pvstackr’s thresholds apply
(see [the diagnostics](#diagnostics)); when both are close to singular
in the same direction, \\\kappa_A\\ can stay small.

## Diagnostics

The calibration records its diagnostics in `get_diagnostics(fit)$ccc`.
Two of them can give the fit the status `"warning"` or `"blocked"`: the
center separation and the condition number of \\A\\. They are combined
with pvstackr’s other checks, and the most severe result sets the status
(“Status and checks” in
[`?pvstackr_object_contracts`](https://joonho112.github.io/pvstackr/reference/pvstackr_object_contracts.md)
lists all checks). Their thresholds are pvstackr’s rules; the preprint
sets no thresholds for these quantities. The calibration also records
`a_matrix_fro_rel`, which measures how far \\A\\ is from the identity
matrix, and `rho1`, `rho2` and `empirical_fro_rel`, which measure how
far the covariance of the calibrated draws is from \\T\_{\text{MI}}\\
and are zero apart from rounding; these values do not change the status.

### Center separation

The center separation compares the mean of the stacked draws before
calibration with the target estimate, coefficient by coefficient, in
units of the target standard error:

\\ \Delta_c = \max\_{\ell = 1, \dots, p} \frac{\lvert \bar\beta\_\ell -
\bar\beta^{\text{raw}}\_\ell \rvert} {\sqrt{T\_{\text{MI},\ell\ell}}} .
\\

pvstackr stores \\\Delta_c\\ as `delta_c_max` and the ratio for each
coefficient as `delta_c_by_term`; `delta_c_rel` is the root mean square
(RMS) of the same ratios. Measured in standard errors, the ratios can be
compared across coefficients on different scales, such as an intercept
in score points and the coefficient of a 0/1 indicator.

With a flat prior on the fixed effects, which the bundled engine uses
when `prior = NULL`, the mean of the stacked draws differs from
\\\bar\beta\\ only by Monte Carlo error, as noted
[above](#why-calibrate). It differs further when the prior is not flat,
or when a fitting function of your own does not use the stacked weights.
In Appendix A.4, the preprint treats this agreement in its own setting
as an empirical quantity to be measured and reported, not as a
consequence of the point identity (Lee et al. 2026, 36).

Monte Carlo error alone is not small on this scale. The Monte Carlo
standard error of the mean of the draws of coefficient \\\ell\\ is their
standard deviation divided by \\\sqrt{\mathrm{ESS}\_\ell}\\, where
\\\mathrm{ESS}\_\ell\\ is their effective sample size. In target
standard errors it is therefore the ratio of that standard deviation to
\\\sqrt{T\_{\text{MI},\ell\ell}}\\, divided by
\\\sqrt{\mathrm{ESS}\_\ell}\\. With a ratio of 0.5 and 4,000 effective
draws it is about 0.008, close to the warning threshold of 0.01 below. A
fit that differs from the target only by Monte Carlo error can therefore
get the status `"warning"`, and with few effective draws even
`"blocked"`; more draws make this less likely.

pvstackr sets `center_status` from `delta_c_max`. Below 0.01 it is
`"ok"`. From 0.01 the fit gets the status `"warning"`, with the reason
code `center_separation_yellow`; its estimate table is unchanged, since
the table copies the target. From 0.05 the fit is `"blocked"`, with
`center_separation_red`, and its estimate table is empty; its
`get_diagnostics(fit)$ccc` then holds single values grouped as `center`,
`conditioning`, `residual` and `prior`, for example
`get_diagnostics(fit)$ccc$center$delta_c_max`. `delta_c_rel`, the RMS of
the ratios, is recorded and printed in the warning message but does not
change the status.

### Condition number of the calibration matrix

\\\kappa_A\\, stored as `kappa_A`, is the condition number of \\A\\: the
ratio of its largest to its smallest singular value. It is at least 1
and does not change when \\A\\ is multiplied by a nonzero constant. A
large \\\kappa_A\\ means that \\A\\ stretches the centered draws far
more in some directions than in others, as when \\\Sigma\_{\text{raw}}\\
is close to singular, so that small errors in the stacked draws can be
greatly magnified. From \\\kappa_A = 10^6\\ the fit gets the status
`"warning"`, with the reason code `ccc_conditioning_yellow`, and from
\\10^8\\ it is `"blocked"`, with `ccc_conditioning_red`.

## What the calibration does not show

CCC gives the draws the mean and covariance of the target it receives,
whatever that target is. It cannot show that the target is right, and it
does not correct a target that is wrong. The preprint calls the
correction a moment calibration that guarantees agreement with the
supplied target, and in its Appendix B.3 it treats the coverage of
intervals built from the calibrated draws as an empirical question,
which its simulation addresses for targets from model-based fits only
(Lee et al. 2026, 13, 38). The diagnostics compare the stacked fit with
the target and describe the matrix \\A\\; they say nothing about
coverage either.

The interval of a `stack_direct` fit is computed from the estimate, the
standard error and the degrees of freedom of the target. pvstackr’s
reporting rule labels it coverage-claimable only when the target uses
the Barnard–Rubin degrees of freedom; the label records how the interval
was built and does not certify its coverage. [Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
gives the rule, [Choosing a method and reading its
intervals](https://joonho112.github.io/pvstackr/articles/m5-methods-and-coverage.html#interval-rule)
gives its reason, and [The design-based target: BRR–Fay replicate
weights and Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.md)
shows how \\\bar\beta\\ and \\T\_{\text{MI}}\\ are computed.

## Checks

The first check applies the map to simulated draws; it is an arithmetic
check of the two moment properties. The second reads the diagnostics
recorded in the example fit, which hold by construction.

### The moment properties on simulated draws

The example fit has no draws, so the chunk below simulates them: \\S =
4000\\ draws of \\p = 3\\ coefficients from a normal distribution with a
mean and a covariance chosen for the check. It sets a target mean and
covariance in place of \\\bar\beta\\ and \\T\_{\text{MI}}\\ and applies
the map as pvstackr computes it, except that it inverts
\\L\_{\text{raw}}\\ with [`solve()`](https://rdrr.io/r/base/solve.html)
where pvstackr uses
[`forwardsolve()`](https://rdrr.io/r/base/backsolve.html), which gives
the same matrix up to rounding.

``` r

set.seed(5104)                      # fixed seed, so the simulated draws are the same in every run
p <- 3L; S <- 4000L

## simulated draws of p = 3 coefficients with a chosen mean and covariance
Sigma_raw_true <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
raw <- sweep(matrix(rnorm(S * p), S, p) %*% chol(Sigma_raw_true),
             2L, c(10, -3, 0.5), `+`)

## target mean and covariance for the check, in place of beta_bar and T_MI
c_tgt     <- c(458, 47, 2)
Sigma_tgt <- matrix(c(4.0, 0.5, 0.0,
                      0.5, 2.0, 0.3,
                      0.0, 0.3, 1.5), 3, 3, byrow = TRUE)

## the map: beta_cal = c_tgt + A (beta - colMeans(raw)), with A = L_tgt L_raw^{-1}
L_tgt <- t(chol(Sigma_tgt))         # lower Cholesky factor of the target covariance
L_raw <- t(chol(cov(raw)))          # lower Cholesky factor of the covariance of the draws
A     <- L_tgt %*% solve(L_raw)     # the calibration matrix A

cal <- sweep(raw, 2L, colMeans(raw), `-`) %*% t(A)   # center the draws, then multiply by A
cal <- sweep(cal, 2L, c_tgt, `+`)                    # add the target mean

all.equal(unname(colMeans(cal)), c_tgt)      # TRUE (largest difference about 1e-13)
#> [1] TRUE
all.equal(unname(cov(cal)),      Sigma_tgt)  # TRUE (largest difference about 1e-14)
#> [1] TRUE
```

Both comparisons return `TRUE`.
[`all.equal()`](https://rdrr.io/r/base/all.equal.html) accepts a
relative difference below about \\1.5 \times 10^{-8}\\; the largest
absolute differences in this run are about \\10^{-13}\\ for the mean,
whose entries are as large as 458, and \\10^{-14}\\ for the covariance,
the size of rounding errors. The agreement is exact apart from rounding
because the map is built from the mean and covariance of the same draws
that the check summarizes, as in pvstackr. The simulated covariance is
positive definite and well conditioned by construction (the term
`+ diag(p)` and 4000 draws for three coefficients), so the check reaches
neither the error for a matrix that is not positive definite nor the
thresholds of `kappa_A`.

### The diagnostics of the example fit

The example fit recorded its calibration diagnostics when it was made,
and they can be read without the draws:

``` r

dg <- get_diagnostics(fit)

names(dg)
#> [1] "preflight"          "sampler"            "sampler_gate"      
#> [4] "stack_fit"          "stack_fit_warnings" "ccc"
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

`center_status` is `"ok"`: `delta_c_max` is about \\4.4 \times
10^{-14}\\ and `delta_c_rel` about \\2.5 \times 10^{-14}\\, far below
0.01. These values hold by construction. The functions that made the
example fit returned draws built around the target, whose means equal
\\\bar\beta\\ apart from rounding ([Getting started with
pvstackr](https://joonho112.github.io/pvstackr/articles/a1-getting-started.html#example-fit)),
so the values show what the diagnostics look like, not how close the
draws of a sampler come to the target. The list also holds
`delta_c_by_term` and `kappa_A` (see [the diagnostics](#diagnostics)).

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

Williams, Matthew R., and Terrance D. Savitsky. 2021. “Uncertainty
Estimation for Pseudo-Bayesian Inference Under Complex Sampling.”
*International Statistical Review* 89 (1): 72–107.
<https://doi.org/10.1111/insr.12376>.
