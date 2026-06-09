# A1: Getting started — your first calibrated plausible-value fit

Abstract

The fastest path into pvstackr. In a few lines you install the package,
see the single `pv_fit(method = "stack_direct")` call shape, load the
bundled cached fixture fit, and read it three ways: the console
[`print()`](https://rdrr.io/r/base/print.html) box, the tidy
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
table, and one figure. You leave knowing what a calibrated fixed-effect
fit needs, how to read its intervals honestly, and where the full
end-to-end workflow lives next.

``` r

library(pvstackr)
```

## 1. Why you are here

Large-scale assessments such as PISA do not report a single test score
per student. Because proficiency is a latent trait that is never
observed directly, the score is published as a set of *plausible values*
(PVs): \\M\\ draws from each student’s posterior proficiency
distribution (Mislevy 1991; Davier et al. 2009). To honour the
uncertainty those draws encode, an analysis must fit the model once per
PV and combine the \\M\\ results with multiple-imputation (MI) rules
(Rubin 1987). The common shortcut — fit one model to `PV1` and stop —
silently discards the between-imputation variance and reports intervals
that are too narrow.

**pvstackr** gives you an honest fixed-effect fit without that shortcut
and without the bookkeeping. Its default method, `stack_direct`,
produces a fixed-effect report from **one** stacked model fit that is
calibrated to an external Rubin / BRR–Fay target. The stacked fit
supplies the model-side draw cloud and agreement diagnostics; the
reportable fixed-effect estimate and standard error come from the
external target. Intervals are coverage-claimable only on the external
design-based BRR-Fay / Barnard-Rubin path; otherwise the metadata marks
them descriptive.

Two scope notes, stated up front so nothing later surprises you:

- **Fixed effects only.** v0.1 calibrates and reports the fixed-effect
  block \\\beta\_{\text{FE}}\\ (intercept and slopes). Variance
  components are estimated but **not** calibrated, and are not part of
  the reportable output.
- **“One fit” is an architecture statement, not a speed claim.** It
  describes the *topology* of the method (one stacked fit instead of
  \\M\\ separate ones); this vignette makes no benchmarked timing
  claims.

The road map for this article: install → the one call → load a cached
fit → read it three ways → read the intervals honestly → where to go
next.

## 2. What you will leave with

After this vignette you will be able to:

- recognise the single canonical `pv_fit(method = "stack_direct")` call
  and name the pieces it needs (an external target and a model backend);
- load the bundled cached fixture fit and read it three ways — the
  console [`print()`](https://rdrr.io/r/base/print.html) box, the tidy
  [`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
  table, and one coefficient figure;
- tell at a glance, from the estimate table itself, whether a row’s
  interval is **coverage-claimable** or merely **descriptive**, and say
  why;
- know exactly which companion vignette to open for the full workflow,
  the accessor tour, method comparison, real PISA data, and the
  underlying method.

This article runs in **seconds** because it never fits a model: it loads
a precomputed fit object and reads off it. There is no live MCMC
anywhere below.

## 3. Install and load

Once the package is on CRAN:

``` r

# install.packages("pvstackr")
```

Until then, install the development version from GitHub:

``` r

# install.packages("pak")
pak::pak("joonho112/pvstackr")
```

The one line to remember — already run in the setup above — is:

``` r

library(pvstackr)
```

Reading a cached fit (everything in this vignette) needs **no**
modelling backend and no extra dependencies. A *live* fit on your own
data additionally needs a fitting engine; that is the subject of
vignette A2.

## 4. The one call (conceptual)

Every pvstackr analysis funnels through a single dispatcher,
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md).
Its canonical `stack_direct` form looks like this. The chunk is **not
run** here — it only shows the *shape* of the call, since a real fit
needs an external target and a backend that this article does not
assemble:

``` r

fit <- pv_fit(
  data          = your_data,
  formula       = OUTCOME ~ x + female,           # OUTCOME is a placeholder
  target        = your_target,                    # from pv_brr_target()
  method        = "stack_direct",
  control       = pv_control(method = "stack_direct"),
  fit_function  = your_fit_function,              # the modelling backend
  draws_function = your_draws_function            # optional posterior draws
)
```

A few words on the pieces:

- **`OUTCOME` is a placeholder.** Internally pvstackr substitutes it
  with each plausible-value column in turn (`PV1READ`, `PV2READ`, …) so
  you write the model once, not \\M\\ times.
- **`target =`** supplies the external Rubin / BRR–Fay variance target
  built by
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md).
  `stack_direct` intervals are coverage-claimable only when backed by an
  external design-based BRR-Fay target with Barnard-Rubin degrees of
  freedom; without that coverage path the fit falls back to descriptive
  intervals.
- **`fit_function =`** is the modelling backend that actually estimates
  the stacked model.

Assembling the target and a backend end to end is the job of **vignette
A2 (the end-to-end workflow)**. Here we skip straight to a fit that has
already been run for you.

## 5. Load the cached fixture fit

The package ships a small synthetic dataset and a cached `stack_direct`
fit so that every example runs offline and deterministically. Load both:

``` r

pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
```

The fixture is a deliberately tiny, **synthetic** PISA-shaped table —
one make-believe country (`CNT == "SYN"`), 12 students nested in a
handful of schools. It carries one continuous covariate `x`, a binary
`female` indicator, **two** plausible reading values `PV1READ` /
`PV2READ` (so \\M = 2\\), **four** BRR replicate weights
`W_FSTURWT1`–`W_FSTURWT4` (so \\R = 4\\), a final weight `W_FSTUWT`, and
student / school identifiers:

``` r

str(pisa_tiny)
#> 'data.frame':    12 obs. of  12 variables:
#>  $ CNT       : chr  "SYN" "SYN" "SYN" "SYN" ...
#>  $ CNTSCHID  : chr  "sch01" "sch01" "sch01" "sch01" ...
#>  $ CNTSTUID  : chr  "stu01" "stu02" "stu03" "stu04" ...
#>  $ x         : num  -1.4 -1.1 -0.7 -0.4 -0.1 0.2 0.5 0.8 1 1.3 ...
#>  $ female    : int  1 0 1 0 1 0 0 1 1 0 ...
#>  $ PV1READ   : int  398 410 423 437 451 466 480 494 507 522 ...
#>  $ PV2READ   : int  402 408 427 435 454 463 483 492 511 519 ...
#>  $ W_FSTUWT  : num  1 1.08 0.94 1.12 1.03 1.19 0.91 1.25 1.1 0.98 ...
#>  $ W_FSTURWT1: num  0.82 1.28 1.01 1.09 0.96 1.35 0.84 1.42 1.21 0.89 ...
#>  $ W_FSTURWT2: num  1.18 0.91 0.85 1.31 1.12 1.02 0.76 1.5 0.99 1.08 ...
#>  $ W_FSTURWT3: num  1.05 1.02 1.12 0.98 1.25 1.14 0.87 1.3 1.34 1.03 ...
#>  $ W_FSTURWT4: num  0.95 1.17 0.93 1.24 1.01 1.31 1.07 1.08 0.88 1.19 ...
```

These are **synthetic, illustrative** numbers — included only so
examples and package checks run without licensed data. They are **not**
real PISA results and must never be cited as such. Real PISA reading
uses \\M = 10\\ plausible values and \\R = 80\\ replicate weights; here
\\M = 2\\ and \\R = 4\\ keep the object small. For working with genuine
PISA files, see vignette A5.

Everything below simply *reads* the `fit` object. Nothing is
re-estimated.

## 6. Read it three ways

### 6.1 The console — `print()`

The quickest look is the result-first
[`print()`](https://rdrr.io/r/base/print.html) box, which leads with the
answer and the things you most need to know before trusting it:

``` r

print(fit)
#> pvstackr fit
#>   method: stack_direct
#>   status: ok
#>   fixed effects: 3
#>   target: external_brr_fay_rubin
#>   draws: not retained
#>   diagnostics: preflight, stack_fit, stack_fit_warnings, ccc
#>   interval note: intervals are descriptive rather than coverage-claimable.
```

It tells you the method (`stack_direct`), that the fit converged
(`status: ok`), that there are 3 fixed effects, the target it was
calibrated against (`external_brr_fay_rubin`), that posterior draws were
not retained, which diagnostics are attached, and — crucially — an
interval note. For this fixture the note reads that *intervals are
descriptive rather than coverage-claimable*; Section 7 explains why, and
how to change it.

### 6.2 The tidy table — `get_estimates()`

[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
returns the reportable data frame: one row per fixed-effect term, with
the estimate, its standard error, the confidence interval, and the
metadata that says how far you may trust that interval. The full table
is wide; here is the compact view most analyses start from:

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

The three terms are `b_Intercept`, `b_x`, and `b_female`. On this
synthetic fixture the intercept sits near 458 score points, the slope on
`x` near 47 points per unit, and the `female` coefficient near 2 points
(with a wide interval — see below). The metadata columns
(`interval_role`, `coverage_claim_allowed`, and the degrees-of-freedom
machinery) are the heart of honest reporting; the full column-by-column
tour, including the *fraction of missing information*, lives in
**vignette A3**.

### 6.3 One figure

A dot-and-interval plot makes the slope estimates easy to read at a
glance. We plot **only the slope coefficients** (`b_x` and `b_female`)
and drop the intercept, whose scale (~458) would otherwise flatten
everything else. A reference line at zero marks the no-effect value:

``` r

slopes <- est[est$term != "b_Intercept", ]
slopes <- slopes[order(slopes$term), ]

y    <- seq_len(nrow(slopes))
xlim <- range(c(slopes$conf_low, slopes$conf_high, 0))

op <- par(mar = c(4.5, 7, 1, 1))
plot(
  slopes$estimate, y,
  xlim = xlim, ylim = c(0.5, nrow(slopes) + 0.5),
  yaxt = "n", ylab = "",
  xlab = "Coefficient (synthetic reading-score points)",
  pch = 19, cex = 1.4, col = "#1f6f9c"
)
abline(v = 0, lty = 2, col = "grey50")
segments(slopes$conf_low, y, slopes$conf_high, y, lwd = 2, col = "#1f6f9c")
points(slopes$estimate, y, pch = 19, cex = 1.4, col = "#1f6f9c")
axis(2, at = y, labels = slopes$term, las = 1)
```

![A horizontal dot-and-interval plot of two slope coefficients from the
synthetic fixture fit. The coefficient on x is about 47 score points
with a narrow interval well to the right of zero. The coefficient on
female is about 2 score points with a wide interval spanning zero from
roughly minus 42 to plus 46. A dashed vertical reference line is drawn
at zero.](a1-getting-started_files/figure-html/coef-figure-1.png)

Slope coefficients from the cached synthetic stack_direct fit, with 95%
descriptive intervals. The intercept is omitted because its scale (~458
score points) would dominate the axis. The dashed line marks zero (no
effect). These are illustrative synthetic values, not real PISA
estimates.

``` r

par(op)
```

The slope on `x` is clearly separated from zero; the `female`
coefficient is small and its interval comfortably spans zero —
unsurprising in a 12-row synthetic toy.

For completeness: this fit retained no posterior draws, so the draw
accessor returns `NULL`. That is expected for a cached `stack_direct`
object and does not affect the fixed-effect report:

``` r

get_draws(fit)
#> NULL
```

## 7. How to read the intervals honestly

pvstackr never hides whether an interval can bear a coverage claim — it
records the answer in the table so you do not have to remember a rule.
Read it directly:

``` r

est[, c("term", "interval_role", "coverage_claim_allowed")]
#>          term             interval_role coverage_claim_allowed
#> 1 b_Intercept descriptive_classic_rubin                  FALSE
#> 2         b_x descriptive_classic_rubin                  FALSE
#> 3    b_female descriptive_classic_rubin                  FALSE
```

For this fixture every row shows
`interval_role == "descriptive_classic_rubin"` and
`coverage_claim_allowed == FALSE`. The intervals are still useful for
*describing* the spread of the estimate — but they are **not** to be
reported as calibrated confidence intervals.

**The rule** is simple and is enforced by the metadata, not by trust:

> A row is `coverage_claim_allowed = TRUE` **only** when the fit is
> `stack_direct` *and* it was calibrated against an external
> **Barnard–Rubin** BRR–Fay target (Barnard and Rubin 1999; Rubin 1987).
> The `per_pv` and `stack_psis` methods are always descriptive. The
> bundled fixture, although `stack_direct`, uses **classic** Rubin
> degrees of freedom (`df_method == "classic"`) rather than the
> Barnard–Rubin small-sample df, so it stays descriptive by design — a
> faithful demonstration of the object surface, not a coverage-claimable
> result.

Two practical consequences are visible in the table above:

- **Read the column, do not assume.** `coverage_claim_allowed` is the
  single field that licenses (or withholds) a confidence claim for each
  row. Let it decide.
- **Tiny \\M\\ means tiny df means wide intervals.** With only \\M = 2\\
  plausible values, the between-imputation degrees of freedom are near
  1, so the \\t\\-multiplier is enormous and the intervals — most
  visibly for `b_female` in the figure — are very wide. This is the MI
  machinery being honest about how little information two imputations
  carry, not a bug.

One related quantity worth a sentence: the **fraction of missing
information** (FMI) measures how much of a coefficient’s uncertainty
comes from imputation rather than sampling. It drives the degrees of
freedom and is reported alongside the estimates; its full treatment —
and how it interacts with `coverage_claim_allowed` — is in **vignette
A3**.

## 8. Where to next

You have loaded a fit and read it honestly. From here:

- **A2 · The end-to-end workflow** — declare the design with
  [`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md),
  assemble the target with
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
  run `pv_fit(method = "stack_direct")`, and report; plus the conceptual
  live-backend sketch and why reportable fits need `center = "target"`.
- **A3 · Reading results & what to report** — the full accessor family,
  every estimate-table column, the interval metadata in depth, the
  fraction of missing information, and a copy-paste reporting checklist.
- **A4 · Comparing methods** —
  [`pv_compare_methods()`](https://joonho112.github.io/pvstackr/reference/pv_compare_methods.md)
  across `stack_direct`, `per_pv`, and `stack_psis`; agreement
  diagnostics; and the two cautions (design variance is not optional; a
  small Pareto-\\\hat k\\ does not guarantee correct variance).
- **A5 · Real PISA data guidance** — connecting to genuine PISA files
  without bundling licensed data: licensing and non-affiliation, the
  design declaration, memory and runtime, and a reproducibility
  checklist.

For the mathematics skipped here — the PV survey setting, the notation,
and the fixed-effect estimand — start the **Method track** with **M1 ·
Foundations and notation**.

Bug reports and feature requests:
<https://github.com/joonho112/pvstackr/issues>.

### Session info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
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
#> [1] pvstackr_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.58         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.0    tools_4.6.0       ragg_1.5.2        bslib_0.11.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.2.0       fs_2.1.0
```

## References

Barnard, John, and Donald B. Rubin. 1999. “Small-Sample Degrees of
Freedom with Multiple Imputation.” *Biometrika* 86 (4): 948–55.
<https://doi.org/10.1093/biomet/86.4.948>.

Davier, Matthias von, Eugenio Gonzalez, and Robert J. Mislevy. 2009.
“What Are Plausible Values and Why Are They Useful?” *IERI Monograph
Series* 2: 9–36.

Mislevy, Robert J. 1991. “Randomization-Based Inference about Latent
Variables from Complex Samples.” *Psychometrika* 56 (2): 177–96.
<https://doi.org/10.1007/BF02294457>.

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.
