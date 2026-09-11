# The design-based target: BRR–Fay replicate weights and Rubin's rules

Abstract

The function
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
computes the design-based target that a `stack_direct` fit reports: for
each plausible value, weighted least-squares estimates of the fixed
effects and their BRR–Fay replicate covariance, combined over the
plausible values with Rubin’s rules. This article gives these formulas,
the classic and Barnard–Rubin degrees of freedom and the fraction of
missing information, with the source of each. It then recomputes the
combined quantities of the example target from the estimates and
covariances stored in it for each plausible value.

``` r

library(pvstackr)
```

## The target and the fit that reports it

[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
estimates the fixed effects of a survey-weighted linear regression from
plausible-value data, together with their covariance matrix. It fits no
Bayesian model. For each plausible value it fits the regression by
weighted least squares, once with the final weight and once with each
replicate weight, and estimates the sampling covariance of the
coefficients from the spread of the replicate estimates. Rubin’s rules
then combine the results for the \\M\\ plausible values. The replicate
weights are assumed to come from balanced repeated replication (BRR)
with Fay’s method, as in PISA. Because the sampling covariance comes
from these replicate weights of the survey design, not from a model, the
target is called design-based.

The result, an object of class `pvstackr_brr_target`, holds the combined
estimates \\\bar\beta\\, their covariance matrix \\T\_{\text{MI}}\\, the
standard errors, the degrees of freedom and the fraction of missing
information;
[`?pv_brr_target`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
lists all its fields. The regression has fixed effects only: a
random-effect term such as `(1 | school)` in the formula stops
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
with an error ([Plausible values, survey weights and
notation](https://joonho112.github.io/pvstackr/articles/m1-foundations-and-notation.html#model)
describes the model and the symbols).

A `stack_direct` fit reports this target: its estimate table copies the
target’s numbers ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)),
and `get_target(fit)` returns the target.

pvstackr takes the covariance from the target because a single stacked
fit cannot supply it. The draws of the stacked fit do not contain the
variation between plausible values, \\B\\, which exists only across
separate analyses of the plausible values, as the companion preprint
notes (Lee et al. 2026, sec. 4.4). Their variance within a plausible
value comes from the model, not from the replicate weights. The Cholesky
calibration correction (CCC) gives the stacked draws the mean
\\\bar\beta\\ and the covariance \\T\_{\text{MI}}\\ of the target
([Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md)).

pvstackr’s reporting rule labels the intervals of a `stack_direct` fit
coverage-claimable, that is, readable as confidence intervals with
nominal coverage, only when its target uses the Barnard–Rubin degrees of
freedom described below ([Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
states the rule in full).

The code below reads the target of the example fit that ships with
pvstackr and recomputes parts of it with base R; it fits no model. The
target was computed by
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
from the bundled synthetic data, which [Getting started with
pvstackr](https://joonho112.github.io/pvstackr/articles/a1-getting-started.html#example-fit)
describes together with the example fit:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # the target of the fit, a pvstackr_brr_target object

class(tg)
#> [1] "pvstackr_brr_target" "list"
c(M = tg$M, R = tg$R, fay_k = tg$fay_k)
#>     M     R fay_k 
#>   2.0   4.0   0.5
```

The synthetic data have \\M = 2\\ plausible values, \\R = 4\\ replicate
weights and the Fay coefficient \\k = 0.5\\.

## The BRR–Fay replicate covariance

For plausible value \\m\\, let \\\hat\beta_m\\ be the weighted
least-squares estimate of the \\p\\ fixed-effect coefficients with the
final weights \\W_i\\, and \\\hat\beta_m^{(r)}\\ the estimate with
replicate weight \\W_i^{(r)}\\, for \\r = 1, \dots, R\\. pvstackr
computes them with [`lm.wfit()`](https://rdrr.io/r/stats/lmfit.html)
from the stats package, \\M(R + 1)\\ fits in all. The BRR–Fay replicate
covariance of \\\hat\beta_m\\ is

\\ \hat U_m = \frac{1}{R(1-k)^2} \sum\_{r=1}^{R}
\big(\hat\beta_m^{(r)} - \hat\beta_m\big)\big(\hat\beta_m^{(r)} -
\hat\beta_m\big)^\top , \\

where \\k\\ is the Fay coefficient `fay_k`. The sum adds the outer
products of the deviations of the replicate estimates from the estimate
with the final weights. \\\hat U_m\\ is a replication estimate of
variance: it needs only the replicate estimates, not a linearization of
the estimator. In the terms of multiple imputation, \\\hat\beta_m\\ and
\\\hat U_m\\ are the complete-data estimate and its variance for the
data set completed with plausible value \\m\\.

Fay’s method changes the weights less than plain BRR: in each replicate
it multiplies the weights of one half-sample by \\2 - k\\ and those of
the other by \\k\\, where plain BRR uses 2 and 0 ([Plausible values,
survey weights and
notation](https://joonho112.github.io/pvstackr/articles/m1-foundations-and-notation.html#pv-design)
describes the replicate weights). For a linear statistic the replicate
deviations are then \\1 - k\\ times those of plain BRR, so their mean
square is too small by the factor \\(1-k)^2\\, and multiplying by
\\1/(1-k)^2\\ corrects it (Judkins 1990, 224–25). The factor \\1/R\\
averages over the replicates. pvstackr requires positive replicate
weights, so it does not accept plain BRR weights, which are zero for
half of the sample.

The target stores the factor \\1/\\R(1-k)^2\\\\ as
`fay_variance_multiplier`. For PISA 2022, with \\R = 80\\ and \\k =
0.5\\, it is \\1/(80 \cdot 0.25) = 0.05\\, the constant in the variance
formula of the PISA 2022 Technical Report (OECD 2024, ch. 10). In the
example, \\R = 4\\ and \\k = 0.5\\ give \\1/(4 \cdot 0.25) = 1\\, a
value that follows from the small example design.

For each plausible value, `per_pv[[m]]` holds `beta` (\\\hat\beta_m\\),
`U` (\\\hat U_m\\), `replicate_beta` (the \\\hat\beta_m^{(r)}\\, one
column per replicate weight) and `replicate_diff` (the deviations
\\\hat\beta_m^{(r)} - \hat\beta_m\\). For the first plausible value:

``` r

tg$fay_variance_multiplier                 # 1/(R (1 - k)^2), 1 in the example
#> [1] 1
dim(tg$per_pv[[1]]$replicate_diff)         # 3 coefficients x 4 replicate weights
#> [1] 3 4
str(tg$per_pv[[1]]$replicate_diff)
#>  num [1:3, 1:4] 0.0191 0.0562 -0.1768 0.0274 -0.0516 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:3] "b_Intercept" "b_x" "b_female"
#>   ..$ : chr [1:4] "W_FSTURWT1" "W_FSTURWT2" "W_FSTURWT3" "W_FSTURWT4"
```

Each column holds the deviations of the three coefficients for one
replicate weight. \\\hat U_1\\, stored as `per_pv[[1]]$U`, is
`fay_variance_multiplier` times the sum of the outer products of these
columns.

## Rubin’s rules

Rubin’s rules combine the \\M\\ estimates and covariances (Rubin 1987,
76–77):

\\ \bar\beta = \frac{1}{M}\sum\_{m=1}^{M} \hat\beta_m, \qquad \bar U =
\frac{1}{M}\sum\_{m=1}^{M} \hat U_m, \qquad B =
\frac{1}{M-1}\sum\_{m=1}^{M} (\hat\beta_m - \bar\beta)(\hat\beta_m -
\bar\beta)^\top, \\

\\ T\_{\text{MI}} = \bar U + \Big(1 + \frac{1}{M}\Big) B . \\

- \\\bar\beta\\ (field `beta_bar`, also stored as `beta`) is the average
  of the \\M\\ estimates, the target’s estimate of the coefficients.
- \\\bar U\\ (`U_bar`) is the average of the replicate covariances: the
  complete-data covariance, that is, the sampling covariance an estimate
  would have if the plausible value were the observed outcome.
- \\B\\ (`B`) is the covariance matrix of the \\M\\ estimates, with
  divisor \\M - 1\\: the variation between plausible values, which comes
  from the uncertainty about each student’s proficiency.
- \\T\_{\text{MI}}\\ (`T_MI`, also stored as `total_var`) is the total
  covariance matrix of \\\bar\beta\\. The standard errors (`se`) are the
  square roots of its diagonal.

The factor \\1 + 1/M\\ adds \\B/M\\ to \\\bar U + B\\. It accounts for
the simulation variance of \\\bar\beta\\, which averages a finite number
of plausible values: \\B/M\\ estimates the variance of \\\bar\beta\\
around the average over infinitely many plausible values (Rubin 1987,
89–90; Barnard and Rubin 1999, 950).

For a single statistic this is the combination that the PISA 2022
Technical Report gives for estimates from plausible values, with the
sampling variance of each plausible value from the replicate weights
(OECD 2024, ch. 11). pvstackr applies it to the whole coefficient vector
and keeps the full matrices: CCC uses all of \\T\_{\text{MI}}\\, not
only its diagonal, as the covariance of the calibrated draws
([Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.html#map)).

## Degrees of freedom

A `stack_direct` fit forms each interval as the estimate plus or minus a
\\t\\ quantile times the standard error. The quantile uses the degrees
of freedom of the target and the level `conf_level` of
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
(0.95 by default); the `conf_level` stored in the target does not set
it. Both rules for the degrees of freedom work coefficient by
coefficient with the ratio

\\ \lambda\_\ell = \frac{(1 +
1/M)\\B\_{\ell\ell}}{T\_{\text{MI},\ell\ell}} , \\

the share of the total variance of coefficient \\\ell\\ that comes from
the variation between plausible values ([the fraction of missing
information](#fmi) returns to it). The argument `df_method` of
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
chooses the rule.

`df_method = "classic"`, the default, is the rule of Rubin (1987, 77):

\\ \nu\_\ell=(M-1)/\lambda\_\ell^{2} . \\

Rubin writes it as \\(M - 1)(1 + 1/\mathrm{riv}\_\ell)^2\\, with the
relative increase in variance \\\mathrm{riv}\_\ell = (1 +
1/M)B\_{\ell\ell}/\bar U\_{\ell\ell}\\; the two forms are equal because
\\\lambda\_\ell = \mathrm{riv}\_\ell/(1 + \mathrm{riv}\_\ell)\\. The
rule is derived under the assumption that the complete-data degrees of
freedom, those the analysis would have if no information were missing,
are infinite (Barnard and Rubin 1999, 948). Because \\\lambda\_\ell \le
1\\, \\\nu\_\ell\\ is never smaller than \\M - 1\\.

`df_method = "barnard_rubin"` is the rule of Barnard and Rubin (1999,
948–49). It needs the complete-data degrees of freedom
\\\nu\_{\text{com}}\\, which you give as `df_complete`, and combines the
classic value with an estimate of the observed-data degrees of freedom:

\\
\nu\_{\text{obs},\ell}=\frac{\nu\_{\text{com}}+1}{\nu\_{\text{com}}+3}\\
\nu\_{\text{com}}\\(1-\lambda\_\ell), \qquad
\nu\_{\text{BR},\ell}=\left(\frac{1}{\nu\_\ell}+\frac{1}{\nu\_{\text{obs},\ell}}\right)^{-1}
. \\

Barnard and Rubin derived this adjustment because the classic value can
be many times larger than \\\nu\_{\text{com}}\\ when
\\\nu\_{\text{com}}\\ is small and little information is missing. It
corrects for small complete-data degrees of freedom, not for a small
number of plausible values. \\\nu\_{\text{BR},\ell}\\ is never larger
than \\\nu\_\ell\\ or \\\nu\_{\text{com}}\\, increases with
\\\nu\_{\text{com}}\\, and equals \\\nu\_\ell\\ when
\\\nu\_{\text{com}}\\ is infinite. Unlike \\\nu\_\ell\\, it can be
smaller than \\M - 1\\: when \\\lambda\_\ell\\ is close to 1,
\\\nu\_{\text{obs},\ell}\\ and therefore \\\nu\_{\text{BR},\ell}\\ are
small.

`df_complete` is the complete-data degrees of freedom of the
replicate-weight variance: the degrees of freedom the analysis would
have if the outcome were observed directly instead of through plausible
values. You state it, as one positive number or as one number per
coefficient; pvstackr does not choose it and does not check it against
the replicate weights. `df_complete = Inf` is accepted and gives the
classic degrees of freedom. The adjustment matters when the number is
small.

`df` holds the degrees of freedom of the chosen rule, and `df_classic`
holds the classic values under both rules. With `"classic"`, the field
`df_complete` of the target is `NULL`, and the estimate table of a fit
shows `NA` in its column `df_complete`. `df_method` also sets the
interval labels that a `stack_direct` fit copies: `"classic"` gives
`interval_role = "descriptive_classic_rubin"` and
`coverage_claim_allowed = FALSE`, and `"barnard_rubin"` gives
`"coverage_barnard_rubin"` and `TRUE` for any positive `df_complete`,
`Inf` included. These labels follow pvstackr’s reporting rule (see
[where the formulas and the interval labels come from](#sources)).

The example target uses the classic rule and has no `df_complete`, so
this article prints the formula but does not evaluate it.

## Fraction of missing information

pvstackr stores \\\lambda\_\ell\\ in the target as `fmi` and again as
`lambda`; the two fields are identical, and the estimate table has no
column for them. Barnard and Rubin (1999, 948) use this ratio as an
approximation to the fraction of missing information (FMI) of
coefficient \\\ell\\. Rubin (1987, 77) estimates the fraction of missing
information with an additional term that depends on the degrees of
freedom:

\\ \gamma\_\ell = \frac{\mathrm{riv}\_\ell + 2/(\nu\_\ell +
3)}{\mathrm{riv}\_\ell + 1} , \\

where \\\mathrm{riv}\_\ell\\ is the relative increase in variance,
stored as `riv`, and \\\nu\_\ell\\ are the classic degrees of freedom.
Since \\\lambda\_\ell = \mathrm{riv}\_\ell/(1 + \mathrm{riv}\_\ell)\\,
\\\gamma\_\ell\\ exceeds \\\lambda\_\ell\\ by \\2/\\(\nu\_\ell +
3)(\mathrm{riv}\_\ell + 1)\\\\, and the two are close when \\\nu\_\ell\\
is large. pvstackr computes only \\\lambda\_\ell\\, so an FMI reported
from pvstackr is \\\lambda\_\ell\\.

\\\lambda\_\ell\\ is close to 1 when the variation between plausible
values is much larger than the replicate-weight variance, that is, when
\\\mathrm{riv}\_\ell\\ is large, whatever the number of plausible
values. The classic degrees of freedom are then close to their smallest
value, \\M - 1\\. This is the case in the example, where \\M - 1 = 1\\
and the relative increase in variance is large ([the checks on the
example target](#checks) show the values). The FMI of an analysis of
real data depends on its data and model; the synthetic example says
nothing about it.

## Where the formulas and the interval labels come from

The replicate covariance is the variance estimator of Fay’s method for
balanced repeated replication (Judkins 1990), and the combination is
Rubin’s rules (Rubin 1987) with the degrees of freedom of Rubin (1987)
or Barnard and Rubin (1999). The PISA 2022 Technical Report computes
standard errors from plausible values and replicate weights in the same
way, one statistic at a time (OECD 2024, chs. 10–11). pvstackr applies
these formulas to the coefficients of its survey-weighted linear
regression.

The companion preprint (Lee et al. 2026) builds its default target from
model-based fits of each plausible value and describes a target from
replicate weights as an alternative (Appendix C of the preprint). Its
simulation study compares estimators whose variances and targets come
from model-based fits, on simulated data without replicate weights (Lee
et al. 2026, sec. 5). It therefore provides no evidence about a BRR–Fay
target, and it does not compare model-based with design-based standard
errors. [Real PISA data: access, memory and
reproducibility](https://joonho112.github.io/pvstackr/articles/a5-real-pisa-guidance.html#companion-preprint)
summarizes what the preprint shows.

The interval labels `interval_role` and `coverage_claim_allowed` are
pvstackr’s reporting rule, not a result of any of these sources. Under
this rule only a `stack_direct` fit whose target was made with
`df_method = "barnard_rubin"` and a positive `df_complete` has intervals
labeled coverage-claimable. [Choosing a method and reading its
intervals](https://joonho112.github.io/pvstackr/articles/m5-methods-and-coverage.html#interval-rule)
gives the reason for the rule, and [Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
explains how to read the labels.

## Checks on the example target

The chunks below recompute the combined quantities of the example target
from the estimates and covariances of its two plausible values, stored
in `tg$per_pv`. These are arithmetic checks: they show that the stored
`beta_bar`, `U_bar`, `B`, `T_MI` and `fmi` are Rubin’s rules applied to
the stored pieces for each plausible value. They take those pieces as
given and do not repeat the weighted least-squares fits.

``` r

M     <- tg$M                                    # 2
betas <- sapply(tg$per_pv, function(p) p$beta)   # 3 x M: one column per plausible value
Us    <- lapply(tg$per_pv, function(p) p$U)      # the replicate covariances U_m

beta_bar <- rowMeans(betas)                      # average of the estimates
U_bar    <- Reduce(`+`, Us) / M                  # average replicate covariance
dev      <- betas - beta_bar
B        <- (dev %*% t(dev)) / (M - 1)           # covariance between plausible values
T_MI     <- U_bar + (1 + 1/M) * B                # total covariance

all.equal(unname(beta_bar), unname(tg$beta_bar))
#> [1] TRUE
all.equal(unname(U_bar),    unname(tg$U_bar))
#> [1] TRUE
all.equal(unname(B),        unname(tg$B))
#> [1] TRUE
all.equal(unname(T_MI),     unname(tg$T_MI))
#> [1] TRUE
```

The four comparisons return `TRUE`: the recomputed and the stored values
agree within the default tolerance of
[`all.equal()`](https://rdrr.io/r/base/all.equal.html), a relative
difference of about \\1.5 \times 10^{-8}\\. The ratio \\\lambda\_\ell\\
follows from \\B\\ and \\T\_{\text{MI}}\\ in the same way:

``` r

fmi <- (1 + 1/M) * diag(B) / diag(T_MI)          # lambda for each coefficient
all.equal(unname(fmi), unname(tg$fmi))
#> [1] TRUE
```

The target also stores the relative increase in variance and the degrees
of freedom, with the rule that gave them:

``` r

data.frame(
  term      = tg$fe_names,
  fmi       = round(unname(tg$fmi), 3),
  riv       = round(unname(tg$riv), 2),
  df        = round(unname(tg$df),  2),
  df_method = unname(tg$df_method)
)
#>          term   fmi    riv   df df_method
#> 1 b_Intercept 0.990  94.86 1.02   classic
#> 2         b_x 0.844   5.43 1.40   classic
#> 3    b_female 0.993 146.17 1.01   classic
```

`df_method` is `"classic"` for every coefficient, so `df` holds
\\(M-1)/\lambda\_\ell^2\\: with the rounded values of `fmi`, \\1/0.990^2
\approx 1.02\\, \\1/0.844^2 \approx 1.40\\ and \\1/0.993^2 \approx
1.01\\. The FMI is large because `riv` is large (about 95, 5 and 146),
so the degrees of freedom are near their smallest possible value, \\M -
1 = 1\\, most closely for the intercept and `b_female`.

The example target uses classic degrees of freedom, so it has no
complete-data degrees of freedom:

``` r

tg$df_complete   # NULL for a classic target
#> NULL
```

The estimate table of the fit shows `NA` in its column `df_complete`,
and by pvstackr’s reporting rule the intervals of the fit are labeled
`descriptive_classic_rubin`, with `coverage_claim_allowed = FALSE`.

For one weighted fit to the stacked data, in which every row has the
weight \\\tilde w_i/M\\ (with \\\tilde w_i = W_i/\overline W\\, the
final weight divided by its mean), the stacked fixed-effect point
identity, Theorem 4.1 of the preprint (Lee et al. 2026), gives the
conditions under which the fixed-effect estimate equals \\\bar\beta\\;
[One stacked fit and the Rubin point
estimate](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.html#point-identity)
explains it. [Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md)
shows how the stacked draws are given the covariance \\T\_{\text{MI}}\\.

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
