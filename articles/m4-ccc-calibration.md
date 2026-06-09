# M4: CCC — Cholesky Calibration Correction and the center-separation diagnostics

Abstract

The fourth page of the Method track derives the calibration map that
makes a single stacked fit report the external variance answer. M3
matched the point — the stacked fixed-effect estimate equals the Rubin
mean; M4 matches the variance, forcing the first two moments of the raw
stacked draw cloud to equal the external target’s center and covariance
through an affine Cholesky map. It states the map in full, distinguishes
the reportable target-centering convention from the diagnostic posterior
one, describes the deterministic algorithm and its conditioning
diagnostic, and defines the center-separation diagnostic with both
summaries the package tracks. A self-contained synthetic panel verifies
the moment-match invariant to machine precision. Throughout it is
explicit that CCC is an affine moment match, not a coverage theorem: the
invariants are algebraic by construction, and coverage comes from the
provenance of the target, never from the calibration itself.

``` r

library(pvstackr)
```

M3 settled the **point**. Its stacked fractional bridge showed that one
model fit on the \\N\cdot M\\ stacked rows, each weighted \\1/M\\,
recovers the Rubin mean \\\bar\beta\\ under the regularity conditions of
Theorem 2.2 (EQ-THM22) — but it was emphatic that this is a fixed-effect
*point* identity and says **nothing** about the variance. The variance
answer is a separate object: M2’s external Rubin / BRR–Fay total
\\T\_{\text{MI}}\\ (EQ-TMI), built from the *design*, not from the fit
(Rubin 1987). This page closes the gap. CCC — the **Cholesky Calibration
Correction** — is the deterministic affine map that takes the raw
stacked draw cloud and forces its first two moments to equal that
external target’s. Where Theorem 2.2 matched the point, CCC matches the
**variance**.

The shared light-path load — the same one used across the Method track —
gives us the cached fit and its external target to point at. As in A3
§6.2, this fit was produced with `return_draws = FALSE`, so it retains
**no** posterior draws and `get_draws(fit)` is `NULL`:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # external Rubin / BRR–Fay target (M2's object); Sigma_tgt = T_MI

c(M = tg$M, R = tg$R, fay_k = tg$fay_k)   # 2, 4, 0.5
#>     M     R fay_k 
#>   2.0   4.0   0.5
is.null(get_draws(fit))                   # TRUE -- no retained draws to calibrate
#> [1] TRUE
```

The design is the tiny synthetic one M1 introduced: \\M = 2\\ plausible
values, \\R = 4\\ replicate weights, Fay coefficient \\k = 0.5\\. The
empty draw matrix is the reason Section 6’s invariant must be
**self-contained synthetic** — there are no real draws on this fixture
to run the map on. We can still read the *realized* calibration
diagnostics off the fit, and we do so in Section 6.

**Scope of this page.** M4 derives the calibration map (EQ-CCC) and the
center-separation diagnostic (EQ-DELTAC), and verifies the moment-match
invariant on a self-contained synthetic draw cloud. It estimates nothing
— there is no MCMC below, only base-R linear algebra and read-only
accessor calls against the bundled synthetic fixture. M4 is the only
Method vignette with any randomness, and that is confined to the
synthetic cloud of Section 6 (with a fixed `set.seed`). The *point* side
is **M3** (EQ-THM22); the provenance of
\\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\ is **M2**; which fits may carry
a coverage claim is **M5**.

## 1. What CCC is for

The raw stacked fit hands you a posterior draw cloud for the fixed
effects — call its draws \\\beta_1,\dots,\beta_S\\. That cloud has
*some* mean \\\bar\beta^{\text{raw}}\\ and *some* covariance
\\\Sigma\_{\text{raw}}\\, but neither is the answer you want to report.
Its mean is the raw stacked posterior center — close to \\\bar\beta\\
under the bridge (M3), but not the external object — and its covariance
is the model’s nominal posterior spread, which M2 §6 warned can
**under-cover** because it never sees the survey’s clustering and
weighting. The variance you are entitled to report lives outside the
fit: the design-based \\T\_{\text{MI}}\\ from M2 (Rubin 1987).

CCC reconciles the two. It is the map that takes the raw cloud and
**moves its first two moments onto the external target’s**: it shifts
the mean to a target center \\c\\ and stretches the covariance to
\\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\. After CCC, the calibrated draw
cloud has — by construction — exactly the center and covariance of the
external target, so a fit that reports the calibrated cloud’s moments
*reports the external answer*. That is the whole purpose: M3’s bridge
makes one stacked fit recover the right **point**; CCC makes that same
fit carry the right **variance**. Together they are the “one fit instead
of \\M\\” machinery the package’s honest headline rests on (facts §5) —
an algebraic equivalence, never a speed claim.

Two pieces of the division of labour are worth stating now, because the
rest of the page keeps returning to them. First, CCC supplies an
*algebraic* match of moments; the *coverage* that match is useful for
comes entirely from the **provenance** of \\\Sigma\_{\text{tgt}}\\
(Section 5, and the central caveat). Second, the match is of the **first
two moments only** — center and covariance — not of the full
distributional shape; CCC is an affine map, and affine maps preserve,
they do not create, higher-order structure.

## 2. The map

The map itself is one line. Writing \\\beta_s\\ for the \\s\\-th raw
posterior draw and \\\beta_s^{\text{cal}}\\ for its calibrated image,
**CCC — tag EQ-CCC — is**

\\
\beta_s^{\text{cal}}=c+L\_{\text{tgt}}L\_{\text{raw}}^{-1}\big(\beta_s-\bar\beta^{\text{raw}}\big),\qquad
L L^\top=\Sigma\\ \text{(lower-triangular)}, \\

so by construction \\\operatorname{mean}\\\beta_s^{\text{cal}}\\=c\\ and
\\\operatorname{Cov}\\\beta_s^{\text{cal}}\\=\Sigma\_{\text{tgt}}\\ (to
machine precision, when \\\Sigma\_{\text{raw}}\\ is positive definite).

Read it as three moves applied to every draw, right to left:

- **Recenter.** \\\big(\beta_s-\bar\beta^{\text{raw}}\big)\\ subtracts
  the raw cloud’s own mean, so the recentered draws are centered at the
  origin. This strips off whatever the raw center happened to be,
  leaving only the cloud’s *shape* about its mean.
- **Linear map.** \\L\_{\text{tgt}}L\_{\text{raw}}^{-1}\\ is the affine
  transform that converts the raw covariance into the target covariance.
  Here \\L\_{\text{raw}}\\ and \\L\_{\text{tgt}}\\ are the
  **lower-triangular Cholesky factors** of, respectively, the raw
  covariance \\\Sigma\_{\text{raw}}\\ and the target
  \\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\ — each satisfies \\L
  L^\top=\Sigma\\. Applying \\L\_{\text{raw}}^{-1}\\ “whitens” the
  recentered draws (turns their covariance into the identity), and
  applying \\L\_{\text{tgt}}\\ then “colours” them to have exactly the
  target covariance.
- **Shift to \\c\\.** Adding \\c\\ relocates the now-correctly-shaped
  cloud so its mean is the target center \\c\\ (the reportable choice
  \\c=\bar\beta\\ is Section 3).

The covariance claim is a one-line consequence of how covariance
transforms under a linear map. For a linear transform \\\beta\mapsto
A\beta\\, the covariance becomes \\A\\\Sigma\_{\text{raw}}\\A^\top\\.
With \\A=L\_{\text{tgt}}L\_{\text{raw}}^{-1}\\ and
\\\Sigma\_{\text{raw}}=L\_{\text{raw}}L\_{\text{raw}}^\top\\,

\\ A\\\Sigma\_{\text{raw}}\\A^\top =
L\_{\text{tgt}}L\_{\text{raw}}^{-1}\\\big(L\_{\text{raw}}L\_{\text{raw}}^\top\big)\\L\_{\text{raw}}^{-\top}L\_{\text{tgt}}^\top
= L\_{\text{tgt}}L\_{\text{tgt}}^\top = \Sigma\_{\text{tgt}}, \\

the internal \\L\_{\text{raw}}^{-1}L\_{\text{raw}}\\ and
\\L\_{\text{raw}}^\top L\_{\text{raw}}^{-\top}\\ collapsing to
identities. The mean claim is even simpler: the recenter makes the
pre-shift mean zero, the linear map sends zero to zero, and the shift
sets the mean to \\c\\. Both invariants are therefore *exact by
construction* — which is precisely the fact Section 6 checks numerically
and the central caveat insists must **not** be over-read as coverage.

## 3. Centering conventions

EQ-CCC leaves one quantity to choose: the center \\c\\ the calibrated
draws are shifted to. The choice is what separates a **reportable** fit
from a merely **diagnostic** one, and it is the single most important
control in this whole chapter — exactly the point A2 §6 makes from the
workflow side.

**`center = "target"` is the reportable convention.** It sets
\\c=\bar\beta\\ — the external target’s pooled point from EQ-TMI, the
same \\\bar\beta\\ M3’s Theorem 2.2 recovers. With \\c=\bar\beta\\ the
calibrated cloud’s mean is the external point and its covariance is the
external \\T\_{\text{MI}}\\, so the fit reports the external answer in
full. This is the convention the bundled fixture uses, recorded in
`fit$control$center`, and the one A2 §6 and A3 §3 both insist on for any
reportable result.

``` r

fit$control$center   # "target" -- the reportable convention (c = beta_bar)
#> [1] "target"
```

**`center = "posterior"` is diagnostic only.** It leaves the
fixed-effect center at the raw stacked *posterior* mean
\\\bar\beta^{\text{raw}}\\ while still applying the covariance
calibration and computing the center-separation diagnostics. That hybrid
is a useful CCC sanity check — it lets you *see* how far the raw center
sat from \\\bar\beta\\ — but it is **not** a reportable `stack_direct`
result;
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
refuses to emit a reportable estimate table from it (A2 §6). Treat
`"posterior"` as a diagnostic, never a deliverable.

The reportable choice \\c=\bar\beta\\ is exactly why M3’s verification
panel found the stacked estimate term-by-term equal to \\\bar\beta\\ and
A3 §6.1 found the estimate column equal to the target’s pooled `beta`:
under target centering the reported point *is* the external point by
construction. M4 is the variance counterpart of that same construction —
the covariance leg of EQ-CCC is what makes the SE column equal
\\\sqrt{\operatorname{diag}(T\_{\text{MI}})}\\.

## 4. The deterministic algorithm

CCC is an **algorithm**, not a sampler: given the raw cloud and the
target, every quantity in EQ-CCC is computed by deterministic linear
algebra, and re-running it on the same inputs returns bit-identical
output. There is no Monte Carlo step inside the map. Three
implementation facts complete the picture beyond the bare equation.

**Lower-triangular Cholesky factors.** Both \\L\_{\text{raw}}\\ and
\\L\_{\text{tgt}}\\ are the *lower-triangular* Cholesky factors of their
covariances (\\LL^\top=\Sigma\\). The triangular structure is what makes
the map cheap and stable to apply: solving with a triangular factor is a
back-substitution rather than a general inverse, and the factorization
is unique for a positive-definite \\\Sigma\\. EQ-CCC’s parenthetical
“(lower-triangular)” is fixing this convention so that “\\L\\” is
unambiguous wherever the track cites it.

**Positive-definite inputs, fail-fast — no automatic repair.** The map
needs \\L\_{\text{raw}}^{-1}\\, which exists only if
\\\Sigma\_{\text{raw}}\\ is positive definite, and \\L\_{\text{tgt}}\\
likewise requires a positive-definite
\\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\. In v0.1.0, CCC factors both
covariances **as given** and applies **no** automatic Tikhonov or
near-PD repair to either. When the raw draw cloud is **ill-conditioned**
— near-collinear fixed effects, or too few effective draws relative to
the dimension — that condition surfaces in one of two ways rather than
being silently patched: a matrix that is *not* positive definite is
rejected outright on the Cholesky path (the factorization aborts with a
“must be positive definite for CCC calibration” error), while a matrix
that is positive definite but poorly conditioned clears the Cholesky and
is instead flagged by the \\\kappa_A\\ conditioning diagnostic below.
The remedy is upstream — increase effective draws, reduce collinearity,
or simplify the fixed-effect block before reporting — never an internal
ridge.

**The conditioning diagnostic \\\kappa_A\\.** Because the map’s
stability hinges on \\L\_{\text{raw}}^{-1}\\, the package records a
conditioning diagnostic, **\\\kappa_A\\** — a numerical-stability
summary of the transform \\A=L\_{\text{tgt}}L\_{\text{raw}}^{-1}\\ (it
lives in the `ccc` diagnostic block as `kappa_A`). A well-conditioned
raw cloud gives a small \\\kappa_A\\ and a trustworthy map; a large
\\\kappa_A\\ flags that the inverse was delicate and the calibration
should be read with care. It is a *trust* diagnostic for the algorithm,
read alongside the center-separation numbers of the next section — never
a coverage statistic.

## 5. Center-separation diagnostics

The covariance leg of EQ-CCC is exact by construction, but the
**center** leg carries information worth auditing: *how far apart were
the raw stacked center and the external point \\\bar\beta\\ before CCC
shifted things?* If the bridge (M3) is behaving, the raw stacked
posterior mean should already sit close to \\\bar\beta\\; a large gap is
a warning that the stacked fit and the external target disagree on the
point, which no amount of covariance calibration would fix. The
center-separation diagnostic measures exactly that gap.

**EQ-DELTAC — the center-separation diagnostic.** \\\Delta_c\\ is the
**relative gap between the raw stacked center and \\\bar\beta\\,
measured on the target-SE scale** — the per-coefficient difference
\\\big(\text{raw center}-\bar\beta\big)\\ divided by that coefficient’s
target standard error \\\sqrt{\operatorname{diag}(T\_{\text{MI}})}\\,
then summarised across terms. Putting it on the SE scale makes it
dimensionless and comparable across coefficients on wildly different
scales (an intercept near 458 versus a slope near 2). The summary is
then read against **green / yellow / red tiers** that gate the fit:
green is a negligible separation (the bridge and target agree on the
point), yellow warns, and red blocks a reportable result.

The package computes **two summaries** of the per-term gaps:

- **`delta_c_rel`** is the root-mean-square (RMS) target-SE-scaled
  separation across terms. It is a descriptive diagnostic that
  summarizes the overall size of the center shift.
- **`delta_c_max`** is the maximum absolute target-SE-scaled separation
  — the single worst-separated coefficient. This is the reportable
  **NO-SEND gate** used for the green / yellow / red center-separation
  bands.

So when this vignette says “the reportable gate,” it means
**`delta_c_max`**. A3 §6.3 previews the same pair: an RMS diagnostic
(`delta_c_rel`) and a max-over-terms gate (`delta_c_max`).

**\\\Delta_c\\ has two summaries — state both.** `delta_c_rel` is the
RMS diagnostic; `delta_c_max` is the reportable **NO-SEND gate**.
Center-separation bands are keyed to `delta_c_max`: green below `1e-2`,
yellow from `1e-2` to below `5e-2`, and red at or above `5e-2`. Neither
summary is a coverage statistic; they are algorithmic trust diagnostics
for deciding whether the calibrated fit may be reported.

## 6. Verify the CCC invariant in code

The natural way to demonstrate EQ-CCC would be to apply it to the fit’s
own raw draws and confirm the calibrated cloud has mean \\c\\ and
covariance \\\Sigma\_{\text{tgt}}\\. We cannot do that here for a
concrete reason stated at the top of the page: **this fixture retains no
draws** — it was produced with `return_draws = FALSE`, so
`get_draws(fit)` is `NULL` (A3 §6.2), and there is no raw cloud to
calibrate. The invariant is therefore demonstrated on a **self-contained
synthetic** draw cloud, built only to exercise the algebra of the map.
This is the only randomness in the Method track, and it is seeded.

We manufacture an arbitrary raw cloud with a non-trivial mean and
covariance, choose target moments to impose (standing in for
\\c=\bar\beta\\ and \\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\), and apply
EQ-CCC exactly as written:

``` r

set.seed(5104)                      # M4 is the only Method vignette with randomness
p <- 3L; S <- 4000L

## an arbitrary raw draw cloud with a non-trivial mean and covariance
Sigma_raw_true <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
raw <- sweep(matrix(rnorm(S * p), S, p) %*% chol(Sigma_raw_true),
             2L, c(10, -3, 0.5), `+`)

## target moments to impose (stand in for c = beta_bar and Sigma_tgt = T_MI)
c_tgt     <- c(458, 47, 2)
Sigma_tgt <- matrix(c(4.0, 0.5, 0.0,
                      0.5, 2.0, 0.3,
                      0.0, 0.3, 1.5), 3, 3, byrow = TRUE)

## EQ-CCC: beta_cal = c + L_tgt L_raw^{-1} (beta - colMeans(raw))
L_tgt <- t(chol(Sigma_tgt))         # lower-triangular Cholesky factor of the target
L_raw <- t(chol(cov(raw)))          # lower-triangular Cholesky factor of the (empirical) raw cov
A     <- L_tgt %*% solve(L_raw)     # the affine transform L_tgt L_raw^{-1}

cal <- sweep(raw, 2L, colMeans(raw), `-`) %*% t(A)   # recenter, then linear map
cal <- sweep(cal, 2L, c_tgt, `+`)                    # shift to the target center

all.equal(unname(colMeans(cal)), c_tgt)      # TRUE  (max|diff| ~ 1e-13)
#> [1] TRUE
all.equal(unname(cov(cal)),      Sigma_tgt)  # TRUE  (max|diff| ~ 1e-15)
#> [1] TRUE
```

Both invariants hold: the calibrated cloud’s mean equals the imposed
center \\c\\ and its covariance equals the imposed target
\\\Sigma\_{\text{tgt}}\\, each to machine precision. That is EQ-CCC
realised — the affine map does exactly what Section 2’s algebra said it
must.

**Construction note — why these hit machine precision, and what
production requires.** In this panel \\L\_{\text{raw}}\\ is the Cholesky
factor of the **empirical** `cov(raw)` and the draws are recentered on
the **empirical** `colMeans(raw)` — the very same cloud’s own moments.
Because the map is built from and applied to the identical empirical
moments, the cancellation in Section 2’s covariance derivation is exact
and both invariants land at ~\\10^{-14}\\. Production CCC adds no ridge
to this: it factors \\\Sigma\_{\text{raw}}\\ exactly as this panel does,
by plain Cholesky, and **requires** only that \\\Sigma\_{\text{raw}}\\
and the external target both be positive definite (Section 4). The
difference is one of guardrails, not arithmetic — in production a
non-positive-definite \\\Sigma\_{\text{raw}}\\ is rejected on the
Cholesky path and a poorly conditioned one is flagged through
\\\kappa_A\\, whereas this synthetic cloud is well-conditioned by
construction (`Sigma_raw_true` is positive-definite with a `+ diag(p)`
term), so those guardrails never fire. The plain inverse used here is
exactly the production map, not a simplified stand-in for it.

We can still read the fit’s **realized** center-separation numbers
directly off the fixture, even though there are no draws — they were
recorded in the `ccc` diagnostic block when the fixture was calibrated.
This is read-only:

``` r

dg <- get_diagnostics(fit)

names(dg)                                                    # preflight, stack_fit, stack_fit_warnings, ccc
#> [1] "preflight"          "stack_fit"          "stack_fit_warnings"
#> [4] "ccc"
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

`center_status == "ok"` says the center separation fell in the green
tier — the raw stacked center and \\\bar\beta\\ agree. The two summaries
are both effectively zero: `delta_c_rel` (the **RMS** diagnostic of
Section 5) is about \\2.5\times10^{-14}\\ and `delta_c_max` (the
**max-over-terms reportable gate**) about \\4.4\times10^{-14}\\. Both
are reported, and the green / yellow / red gate is keyed to
`delta_c_max`. The same `ccc` block carries the conditioning diagnostic
`kappa_A` of Section 4 and a per-term breakdown `delta_c_by_term`,
should you wish to inspect the algorithm’s behaviour further.

**CCC is an affine moment match, NOT a coverage theorem.** This is the
central honesty point of the page, and it mirrors A3 §4’s warning
exactly. Everything verified above is **algebraic by construction**:
EQ-CCC *forces* \\\operatorname{mean}=c\\ and
\\\operatorname{Cov}=\Sigma\_{\text{tgt}}\\, so the ~\\10^{-14}\\
agreement merely confirms the linear algebra **ran cleanly** — it is
not, and must never be read as, evidence that the resulting intervals
*cover*. Coverage comes from the **provenance** of
\\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\ — external, design-based,
BRR–Fay (Judkins 1990), pooled by the Rubin rules (Rubin 1987) (M2) —
together with the companion paper’s simulation oracle, **not** from CCC.
Concretely, this fixture’s near-zero \\\Delta_c\\ sits alongside
`coverage_claim_allowed == FALSE` (it carries classic Rubin df, not
Barnard–Rubin; A3 §4, M5): the calibration ran perfectly and the fit is
still, honestly, **descriptive**. A clean CCC is necessary plumbing, not
a coverage certificate.

## 7. Where to next

You now have the variance counterpart to M3’s point identity in full:
the CCC affine map (EQ-CCC) that forces the raw stacked draw cloud’s
first two moments onto the external target’s \\c=\bar\beta\\ and
\\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\, the target-versus-posterior
centering convention that decides whether a fit is reportable, the
deterministic fail-fast Cholesky algorithm and its \\\kappa_A\\
conditioning check, and the center-separation diagnostic (EQ-DELTAC)
with both the RMS diagnostic (`delta_c_rel`) and the max-over-terms gate
(`delta_c_max`). And its hard boundary: an algebraic moment match, exact
by construction, that confers **no** coverage on its own.

- **M5 · Methods, PSIS, and coverage** — the natural next step, and
  where the caveat of this page is cashed out. CCC’s clean moment match
  is *not* the source of coverage; M5 states precisely why only
  `stack_direct` — calibrated by CCC to an external design-based target
  — is coverage-claimable, while `per_pv` and `stack_psis` carry
  model-based covariance and are descriptive always.
- **M2 · The BRR–Fay fixed-effect target** — the provenance of
  \\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\ that this page imposes: the
  BRR–Fay sandwich (EQ-BRRFAY) and Rubin pooling (EQ-TMI) that build the
  external object, and the design-variance coverage result that — unlike
  CCC — actually licenses coverage.
- **M3 · The stacked fractional bridge** — the *point* leg this page’s
  variance leg completes: Theorem 2.2 (EQ-THM22) matched \\\bar\beta\\,
  CCC matches the covariance, and the two together are the “one fit
  instead of \\M\\” machinery.

To see these diagnostics from the reporting side — the
`get_diagnostics(fit)$ccc` peek, the `center_status` / `delta_c_rel` /
`delta_c_max` fields, and the reminder that their tiny values do not
confer coverage — return to **A3 · Reading results** (§6.3).

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

Judkins, David R. 1990. “Fay’s Method for Variance Estimation.” *Journal
of Official Statistics* 6 (3): 223–39.

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.
