# M2: The BRR–Fay fixed-effect target — Rubin combining, design variance, and degrees of freedom

Abstract

The second page of the Method track derives the external variance answer
the rest of the machinery is built to report. It assembles the
per-plausible-value design covariance from the
balanced-repeated-replication Fay sandwich, combines the M per-value
fits with the Rubin rules into the external total covariance, and
explains the two degrees-of-freedom regimes — the classic Rubin formula
the bundled fixture carries, and the small-sample Barnard–Rubin
correction that a coverage-claimable target supplies. It relates the
fraction of missing information to the collapse of degrees of freedom at
small M, states the design-variance coverage result qualitatively, and
re-derives the whole Rubin pooling chain on the cached target so every
quantity can be checked against the stored object. Where M1 named these
equations, M2 derives them.

``` r

library(pvstackr)
```

## 1. The target as an external object

M1 fixed the vocabulary; this page derives the first piece of machinery
that vocabulary describes — the **external fixed-effect target**. The
single most important fact about it is structural: the variance answer
for \\\beta\_{\text{FE}}\\ does **not** come out of the model fit. It is
assembled separately, from the plausible values and the replicate-weight
design, and the stacked fit is then *calibrated to report it* (the
calibration is M4’s subject). The target is an object in its own right —
class `pvstackr_brr_target`, reached with `get_target(fit)` — carrying
the pooled coefficients, the total covariance, the degrees of freedom,
and the fraction of missing information.

Why insist on this separation? Because **provenance is what licenses a
coverage claim**, not the arithmetic of the fit. A correct point
estimate paired with a model-based standard error can under-cover badly
(Section 6); what restores nominal coverage is a variance built from the
*design* — the BRR–Fay sandwich (Section 2) combined by the Rubin rules
(Section 3). A3 §6.1 met this object from the reading side: the standard
errors in
[`get_estimates()`](https://joonho112.github.io/pvstackr/reference/get_estimates.md)
are \\\sqrt{\operatorname{diag}(T\_{\text{MI}})}\\ and the point
estimates are the pooled `beta_bar`, both read straight off the target
under the reportable convention `center = "target"`. This vignette opens
that object up and derives every field in it.

**Scope of this page.** M2 derives the target. It estimates nothing —
there is no MCMC below, only base-R linear algebra and read-only
accessor calls against the bundled synthetic fixture. The calibration
that maps this target onto a stacked fit is **M4**; why only
`stack_direct` fits may carry a coverage claim is **M5**.

The shared light-path load — the same one used across the Method track —
gives us the fit and its target to point at:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # external Rubin / BRR–Fay target (class pvstackr_brr_target)

class(tg)
#> [1] "pvstackr_brr_target" "list"
c(M = tg$M, R = tg$R, fay_k = tg$fay_k)   # 2, 4, 0.5
#>     M     R fay_k 
#>   2.0   4.0   0.5
```

The design is the tiny synthetic one M1 introduced: \\M = 2\\ plausible
values, \\R = 4\\ replicate weights, Fay coefficient \\k = 0.5\\.

## 2. Per-PV design variance: the BRR–Fay sandwich

The first ingredient is the per-plausible-value covariance \\\hat U_m\\.
In this package it is **design-based**: rather than trust the model’s
nominal covariance, it is built from PISA’s replicate-weight machinery.
As M1 §2.1 set out, alongside each student’s final weight \\W\_{ij}\\
come \\R\\ replicate weights \\W\_{ij}^{(r)}\\, each a balanced
perturbation of the sample. Re-estimating the fixed effects on every
replicate and measuring how much they move is the design-based sampling
variance.

For plausible value \\m\\, write \\\hat\beta_m\\ for the full-weight
estimate and \\\hat\beta_m^{(r)}\\ for the estimate on the \\r\\-th
replicate weight. The **BRR–Fay sandwich** — **tag EQ-BRRFAY** — is

\\ \hat
U_m^{\text{BRR-Fay}}=a_d\sum\_{r=1}^{R}\big(\hat\beta_m^{(r)}-\hat\beta_m\big)\big(\hat\beta_m^{(r)}-\hat\beta_m\big)^\top,\quad
a_d=\tfrac1{R(1-k)^2}. \\

Read it as a sum of outer products of replicate deviations. Each term
\\\big(\hat\beta_m^{(r)}-\hat\beta_m\big)\\ records how far the \\r\\-th
replicate pulls the estimate from its full-weight value; squaring and
summing those deviations across the \\R\\ replicates measures the
sampling spread, and the **Fay multiplier** \\a_d\\ rescales it to the
correct variance (Judkins 1990).

The multiplier carries the whole Fay correction. PISA uses Fay’s variant
of BRR, which shrinks each perturbation by the coefficient \\k\\ rather
than zeroing out half-samples; the factor \\1/(1-k)^2\\ undoes that
shrinkage and \\1/R\\ averages over the replicates. For PISA reading
(\\R = 80\\, \\k = 0.5\\) this gives \\a_d = 1/\[80\cdot 0.25\] =
0.05\\.

**The fixture’s \\a_d = 1\\ is an arithmetic coincidence — do not
generalise.** With the fixture’s \\R = 4\\ and \\k = 0.5\\, \\a_d =
1/\[4\\(1-0.5)^2\] = 1/(4\cdot 0.25) = 1\\. That is a quirk of this one
tiny design — it matches the `fay_variance_multiplier` field on the
target — **not** a typical value. Real PISA reading has \\a_d = 0.05\\
(M1 §2.1).

The replicate ingredients are stored on the target, one entry per
plausible value. We can peek at the deviation matrix
\\\big(\hat\beta_m^{(r)}-\hat\beta_m\big)\\ that EQ-BRRFAY sums over —
for the first plausible value it is a \\3 \times R\\ matrix (three
fixed-effect terms, \\R = 4\\ replicates):

``` r

tg$fay_variance_multiplier                 # a_d = 1 on this fixture (coincidence)
#> [1] 1
dim(tg$per_pv[[1]]$replicate_diff)         # 3 x 4 : (beta_hat_m^(r) - beta_hat_m)
#> [1] 3 4
str(tg$per_pv[[1]]$replicate_diff)
#>  num [1:3, 1:4] 0.0191 0.0562 -0.1768 0.0274 -0.0516 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:3] "b_Intercept" "b_x" "b_female"
#>   ..$ : chr [1:4] "W_FSTURWT1" "W_FSTURWT2" "W_FSTURWT3" "W_FSTURWT4"
```

Each column is one replicate’s deviation from the full-weight estimate;
the sandwich in EQ-BRRFAY forms \\a_d\sum_r\\ of the outer products of
these columns to give \\\hat U_m\\. (Section 7 reads the assembled
\\\hat U_m\\ back as `per_pv[[m]]$U` and pools it.)

## 3. Rubin combining

With a per-PV point estimate \\\hat\beta_m\\ and its design covariance
\\\hat U_m\\ in hand for each \\m = 1,\dots,M\\, the **Rubin rules**
combine the \\M\\ fits into a single external answer (Rubin 1987). **Tag
EQ-TMI** is the three combining quantities followed by the boxed total:

\\ \bar\beta=\tfrac1M\textstyle\sum_m\hat\beta_m,\quad \bar
U=\tfrac1M\textstyle\sum_m\hat U_m,\quad
B=\tfrac1{M-1}\textstyle\sum_m(\hat\beta_m-\bar\beta)(\hat\beta_m-\bar\beta)^\top,
\\

\\ \boxed{\\T\_{\text{MI}}=\bar U+\big(1+\tfrac1M\big)B\\}. \\

The three pieces decompose the uncertainty:

- \\\bar\beta\\ is the **Rubin mean** — the simple average of the per-PV
  estimates, and the reportable point. (M3 shows one stacked fit
  recovers exactly this value under stated conditions; here it is the
  combining rule that defines it.)
- \\\bar U\\ is the **within-imputation covariance** — the average of
  the per-PV design covariances \\\hat U_m\\. It is the sampling
  uncertainty you would have if the plausible values agreed perfectly.
- \\B\\ is the **between-imputation covariance** — the spread of the
  per-PV estimates *around* their mean. It is the extra uncertainty the
  measurement leaves, the part that exists precisely because the
  proficiency is imputed rather than observed.

The total \\T\_{\text{MI}}\\ adds them, inflating the between-component
by the finite- \\M\\ factor \\(1+1/M)\\ that corrects for estimating
\\B\\ from only \\M\\ imputations. This \\T\_{\text{MI}}\\ is the
external target covariance — the variance answer the whole package is
engineered to report, and the \\\Sigma\_{\text{tgt}}\\ that M4’s
calibration imposes on the stacked draws. Section 7 re-derives every
quantity here from the stored per-PV pieces and checks it against the
target.

## 4. Small-sample degrees of freedom

A variance is only half of an interval; the other half is the degrees of
freedom that set the reference distribution’s tails. Two regimes matter,
and the distinction governs whether an interval may carry a coverage
claim.

The **classic Rubin** degrees of freedom derive from the relative
increase in variance alone and assume the complete-data analysis has
effectively unlimited df (Rubin 1987). They are what the bundled fixture
carries (Section 7).

The small-sample correction is the **Barnard–Rubin** degrees of freedom
\\\nu\_{\text{BR}}\\ — **tag EQ-BARNARD** (Barnard and Rubin 1999). It
computes the combined-data df from two inputs: the **fraction of missing
information** \\\lambda_k\\ (the FMI \\\gamma_k\\ of Section 5) and the
**complete-data df** \\\nu\_{\text{com}}\\, the degrees of freedom the
analysis would have had with no missing information (package field
`df_complete`). When the complete-data analysis is itself small-sample,
\\\nu\_{\text{BR}}\\ corrects the classic formula downward; it is the df
a *coverage-claimable* interval requires. Concretely, with
\\\lambda_k=\gamma_k=(1+1/M)\\B\_{kk}/T\_{\text{MI},kk}\\ and
\\\nu\_{\text{com},k}\\ the complete-data df (`df_complete`), EQ-BARNARD
combines the classic df \\\nu\_{\text{old},k}=(M-1)/\lambda_k^{2}\\ with
an **observed-data df** \\\nu\_{\text{obs},k}\\ harmonically:

\\
\nu\_{\text{obs},k}=\dfrac{\nu\_{\text{com},k}+1}{\nu\_{\text{com},k}+3}\\
\nu\_{\text{com},k}\\(1-\lambda_k),\qquad
\nu\_{\text{BR},k}=\left(\nu\_{\text{old},k}^{-1}+\nu\_{\text{obs},k}^{-1}\right)^{-1}.
\\

The observed-data term pulls \\\nu\_{\text{BR}}\\ below the classic
value when \\\nu\_{\text{com}}\\ is small; as
\\\nu\_{\text{com}}\to\infty\\ it sends \\\nu\_{\text{obs}}\to\infty\\
and \\\nu\_{\text{BR},k}\to\nu\_{\text{old},k}\\, recovering classic
Rubin. The package evaluates this **only** on the Barnard–Rubin path;
the bundled fixture carries classic df, so this vignette prints the
formula but does **not** numerically evaluate it or fabricate a
\\\nu\_{\text{BR}}\\ number for it.

**The cached fixture carries classic df, not \\\nu\_{\text{BR}}\\.**
Every term on this target has `df_method == "classic"` and
`df_complete == NA`, so the `df` it reports is the **classic Rubin**
value, *not* the Barnard–Rubin \\\nu\_{\text{BR}}\\. The complete-data
df that EQ-BARNARD needs is simply not present (`df_complete == NA`),
because the classic path does not use it. The Barnard–Rubin correction
is the path a **reportable, coverage-claimable** target would take — it
is *not* exercised on `pisa_tiny`. Accordingly, this vignette reads the
**classic** `df` / `df_method` / `df_complete` off the fixture
(mirroring A3 §4) and does **not** compute or display a
\\\nu\_{\text{BR}}\\ number here.

Section 7 reads `df`, `df_method`, and `df_complete` directly off the
target so the regime is visible in the object, exactly as A3 §4 reads
them off the estimate table.

## 5. Fraction of missing information

The quantity that drives the degrees of freedom is the **fraction of
missing information** (FMI). For coefficient \\k\\ it is the share of
the total uncertainty contributed by the between-imputation component —
**tag EQ-FMI**:

\\ \gamma_k=\dfrac{(1+1/M)\\B\_{kk}}{T\_{\text{MI},kk}}, \\

where \\B\_{kk}\\ and \\T\_{\text{MI},kk}\\ are the \\k\\-th diagonal
entries of the between-imputation covariance and the total (Rubin 1987).
The numerator is exactly the term EQ-TMI adds to \\\bar U\_{kk}\\, so
\\\gamma_k\\ is the proportion of \\T\_{\text{MI},kk}\\ that exists
*because* the proficiency is imputed.

FMI is what collapses the degrees of freedom at small \\M\\. When
\\\gamma_k\\ is near one almost all of a coefficient’s uncertainty is
between-imputation, there is little information per coefficient, and the
df fall toward their floor — which is why a small-\\M\\ design produces
wide intervals and df near 1. Section 7 re-derives \\\gamma_k\\ from
\\B\\ and \\T\_{\text{MI}}\\ and confirms it against the stored `fmi`.

**High FMI here is \\M = 2\\ mechanics, not a finding.** The fixture’s
FMI values run roughly 0.84 to 0.99 because, with only \\M = 2\\
plausible values, the \\(1 + 1/M) = 1.5\\ factor multiplies a
between-imputation covariance estimated from a single pair of fits, so
the between-share dominates and the df collapse toward 1. This is the MI
machinery being honest about how little information two imputations
carry. Real PISA reading uses \\M = 10\\, where FMI is far smaller and
the Barnard–Rubin correction operates in its intended regime. Do not
read this fixture’s FMI as a substantive result.

## 6. The design-variance coverage result

Everything above is arithmetic on a target. The reason that arithmetic
is built the way it is — design-based \\\hat U_m\\, design-based
\\T\_{\text{MI}}\\ — is a **coverage** result, and that result is
empirical rather than something this vignette can compute on a
12-student synthetic fixture.

**Design variance is not optional (companion paper, in preparation).** A
correct weighted point estimate paired with **model-based** standard
errors can **under-cover** — its nominal-95% intervals capture the truth
far less than 95% of the time — because the model’s covariance does not
see the survey’s clustering and weighting. Substituting the
**design-based** BRR–Fay \\T\_{\text{MI}}\\ derived above **restores
nominal coverage**. This is established by simulation in the companion
methods paper (in preparation) and is reported there *qualitatively*;
the package neither bundles that evidence nor reproduces its numbers.
The lesson is about the *provenance* of the target — external and
design-based — which is precisely what Section 1 said licenses a
coverage claim, and what M5 ties to the `coverage_claim_allowed` flag.

So the target’s value is not that its identities check to high precision
(they do, in Section 7 — but that is algebra), but that its variance is
built from the design. The precise statement of which fits may therefore
carry a coverage claim is M5; the calibration that delivers this
target’s variance to a stacked fit is M4.

## 7. Verify Rubin pooling on the cached target

We now re-derive the entire Rubin pooling chain — EQ-TMI and EQ-FMI —
from the stored per-PV pieces, and check each quantity against the value
the package computed. This is a read-and-recompute exercise: `tg$per_pv`
holds, for each plausible value \\m\\, the estimate \\\hat\beta_m\\
(`beta`) and the assembled BRR–Fay covariance \\\hat U_m\\ (`U`); we
rebuild \\\bar\beta\\, \\\bar U\\, \\B\\, \\T\_{\text{MI}}\\, and the
FMI from them.

``` r

M     <- tg$M                                    # 2
betas <- sapply(tg$per_pv, function(p) p$beta)   # 3 x M : per-PV beta_hat_m
Us    <- lapply(tg$per_pv, function(p) p$U)      # list of per-PV U_hat_m (BRR–Fay)

beta_bar <- rowMeans(betas)                      # Rubin mean
U_bar    <- Reduce(`+`, Us) / M                  # within-imputation covariance
dev      <- betas - beta_bar
B        <- (dev %*% t(dev)) / (M - 1)           # between-imputation covariance
T_MI     <- U_bar + (1 + 1/M) * B                # EQ-TMI

all.equal(unname(beta_bar), unname(tg$beta_bar)) # TRUE
#> [1] TRUE
all.equal(unname(U_bar),    unname(tg$U_bar))    # TRUE
#> [1] TRUE
all.equal(unname(B),        unname(tg$B))        # TRUE
#> [1] TRUE
all.equal(unname(T_MI),     unname(tg$T_MI))     # TRUE
#> [1] TRUE
```

All four reconstructions reproduce the stored target — the package’s
`beta_bar`, `U_bar`, `B`, and `T_MI` are exactly the Rubin combining
rules applied to the per-PV pieces. The fraction of missing information
follows the same way, from \\B\\ and \\T\_{\text{MI}}\\:

``` r

fmi <- (1 + 1/M) * diag(B) / diag(T_MI)          # EQ-FMI
all.equal(unname(fmi), unname(tg$fmi))           # TRUE
#> [1] TRUE
```

Finally, read the target’s FMI, relative increase in variance, and
degrees of freedom together with the df method — exactly the regime
block of Section 4:

``` r

data.frame(
  term      = tg$fe_names,
  fmi       = round(tg$fmi, 3),
  riv       = round(tg$riv, 2),
  df        = round(tg$df,  2),
  df_method = tg$df_method
)
#>                    term   fmi    riv   df df_method
#> b_Intercept b_Intercept 0.990  94.86 1.02   classic
#> b_x                 b_x 0.844   5.43 1.40   classic
#> b_female       b_female 0.993 146.17 1.01   classic
```

Two things to read off this table. First, the FMI values (~0.84–0.99)
and the df (~1) are the small-\\M\\ mechanics of Section 5: little
information per coefficient, df at the floor. Second — and this is the
point of D-1 above — `df_method` is `"classic"` for every term. The `df`
shown is the classic Rubin value; the Barnard–Rubin \\\nu\_{\text{BR}}\\
of Section 4 is **not** what appears here, and the complete-data df it
would need is absent:

``` r

tg$df_complete   # NA for every term -> the classic path, not Barnard–Rubin
#> NULL
```

**This target is descriptive, not coverage-claimable.** The agreements
above are *algebraic* — the Rubin rules reproduced from their inputs —
and confirm only that the pooling is internally consistent, not that the
intervals cover. The fixture is **descriptive-classic**
(`coverage_claim_allowed == FALSE`; A3 §4): classic df, no complete-data
df, so no Barnard–Rubin small-sample correction. We make no “verified to
machine precision” coverage claim and report no \\\nu\_{\text{BR}}\\
number — both would over-read a synthetic fixture built to demonstrate
the object surface, not to establish coverage.

## 8. Where to next

You now hold the external target in full — the BRR–Fay sandwich
(EQ-BRRFAY) that builds each \\\hat U_m\\, the Rubin rules (EQ-TMI) that
combine them into \\T\_{\text{MI}}\\, the two degrees-of-freedom regimes
(EQ-BARNARD vs classic), and the fraction of missing information
(EQ-FMI) that drives the df. The target is the variance answer; the rest
of the Method track is about *delivering* it from a single fit.

- **M3 · The stacked fractional bridge** — the natural next step: it
  shows that one model fit on the stacked \\N\cdot M\\ rows, each
  weighted \\1/M\\, recovers the Rubin mean \\\bar\beta\\ derived here
  (Theorem 2.2), and is explicit that this is a fixed-effect *point*
  identity, not the variance answer — that answer is this page’s
  \\T\_{\text{MI}}\\.
- Then **M4** (the CCC calibration that makes the stacked draw cloud’s
  covariance equal this target’s \\T\_{\text{MI}}\\) and **M5** (why
  coverage-claimability rests on an external design-based BRR-Fay target
  with Barnard-Rubin degrees of freedom, plus validation evidence rather
  than CCC arithmetic alone) complete the chain.

To see this target from the reporting side — how its fields surface in
the estimate table, what the interval-metadata columns mean, and how FMI
and the classic-vs-Barnard–Rubin df read in practice — return to **A3 ·
Reading results** (§5–§6).

Bug reports and feature requests:
<https://github.com/joonho112/pvstackr/issues>.

### Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
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
#> [1] pvstackr_0.2.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.60         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
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

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.
