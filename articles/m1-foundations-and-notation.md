# M1: Foundations and notation — the plausible-value setting, model, and estimand

Abstract

The opening page of the Method track, and the dictionary the rest of it
cites by tag. This vignette fixes the vocabulary every later method page
reuses without redefinition: what a plausible value is and why a survey
publishes M of them per student, the balanced-repeated-replication
weight design and its Fay variance multiplier, the two-level
within-between model, the fixed-effect estimand, and the
per-plausible-value likelihood that the stacked bridge later averages.
It displays one equation in full and lists the other nine by tag, then
orients you on the bundled synthetic fixture so the symbols have
something concrete to point at. It is a lookup, not a derivation; the
workflow lives in the Applied track.

``` r

library(pvstackr)
```

## 1. Why a notation track

The Applied track (A1–A5) shows you *how to run* a pvstackr analysis:
load a fit, read its intervals, compare methods, point the workflow at
real PISA files. This Method track (M1–M5) is the other half — it ties
each piece of machinery back to the mathematics and, wherever the light
path allows, checks the identity in code.

This page is the **dictionary** for that track. It establishes a
*vocabulary contract*: every symbol and every displayed equation
introduced here is reused **verbatim** by M2–M5. No later page redefines
a symbol, re-letters an equation, or invents a synonym for a quantity
already named here. When M2 needs the BRR–Fay sandwich it writes “see
EQ-BRRFAY”; when M3 invokes the stacked-MLE identity it writes “Theorem
2.2 / EQ-THM22”. So a reader who has internalised M1 should never meet
an undefined symbol downstream.

That division of labour means M1 is deliberately **a lookup, not a
derivation**. It displays exactly one equation in full — the model
itself — and catalogues the other nine by tag, with a one-line gloss and
a pointer to the page that derives each. The derivations live where they
belong:

- **M2** — the external BRR–Fay fixed-effect *target* (Rubin combining,
  the Fay sandwich, small-sample degrees of freedom, the fraction of
  missing information).
- **M3** — the *stacked fractional bridge* and the point identity
  (Theorem 2.2).
- **M4** — *CCC*, the Cholesky Calibration Correction, and its
  center-separation diagnostics.
- **M5** — the precise contrast of the three methods, PSIS tail-shape
  thresholds, and what makes a fit coverage-claimable.

**Scope of this page.** M1 fixes notation and orients you on the bundled
fixture. It estimates nothing — there is no MCMC anywhere below, only
the synthetic data table and read-only accessor calls. If you want to
*run* a fit first, start with **A1 (getting started)** and **A2 (the
end-to-end workflow)**, then return here for the mathematics those pages
defer to.

## 2. The plausible-value survey setting

Large-scale assessments such as PISA never observe a student’s
proficiency directly: it is a latent trait. Rather than publish a single
point score (which would hide the measurement uncertainty), the survey
publishes a set of **plausible values** — \\M\\ draws from each
student’s posterior proficiency distribution, conditioning on their item
responses and background variables (Mislevy 1991; Davier et al. 2009). A
single plausible value is therefore not “the score”; it is one
imputation of an unobserved quantity, and the spread *across* the \\M\\
values encodes how much the assessment leaves uncertain.

The practical consequence is the rule that organises this entire track:
an analysis must fit its model **once per plausible value** and then
combine the \\M\\ results with multiple-imputation rules (Rubin 1987).
Using only `PV1` and stopping silently throws away the
between-imputation variance and reports intervals that are too narrow.

### 2.1 The replicate-weight design and the Fay coefficient

PISA also supplies its sampling-variance machinery as **balanced
repeated replication (BRR)** weights: alongside each student’s final
survey weight \\W\_{ij}\\ come \\R\\ replicate weights
\\W\_{ij}^{(r)}\\, each obtained by perturbing the sample in a balanced
way. Re-estimating a quantity on every replicate and measuring how much
it moves gives a design-based sampling variance.

PISA uses **Fay’s** variant of BRR, which shrinks the perturbation by a
coefficient \\k\\ rather than zeroing out half-samples. The replicate
spread is then rescaled by the **BRR–Fay multiplier**

\\ a_d=\frac{1}{R\\(1-k)^2}, \\

so that the sandwich estimator returns the correct sampling variance (M2
assembles the full estimator as EQ-BRRFAY). Two sets of constants matter
for this track:

| Setting | \\M\\ (PVs) | \\R\\ (replicates) | \\k\\ (Fay) | \\a_d=\tfrac{1}{R(1-k)^2}\\ |
|----|:--:|:--:|:--:|:--:|
| **PISA reading** (reference) | \\10\\ | \\80\\ | \\0.5\\ | \\0.05\\ |
| **Bundled fixture** `pisa_tiny` | \\2\\ | \\4\\ | \\0.5\\ | \\1\\ |

The PISA-reading row is the design the package is built to serve (OECD
2024); the fixture row is the deliberately tiny synthetic design that
makes every example here run offline.

**The fixture’s \\a_d=1\\ is an arithmetic coincidence — do not
generalise.** With \\R=4\\ and \\k=0.5\\ the multiplier works out to
\\1/\[4\\(1-0.5)^2\]=1/(4\cdot 0.25)=1\\. That is a quirk of this one
tiny design (it equals the `fay_variance_multiplier` field you will read
in Section 6), **not** a typical value. Real PISA reading has
\\a_d=0.05\\. Never carry the fixture’s \\a_d=1\\ into a real-data
setting.

## 3. The two-level hybrid model

Students (\\i\\) are nested in schools or clusters (\\j\\), so the model
is two-level. For each plausible value \\m\\ — **tag EQ-MODEL**, the one
equation this page displays in full:

\\ y\_{ij}^{(m)} = \beta_0 + \beta^{W}\big(x\_{ij}-\bar x\_{\cdot
j}\big) + \beta^{B}\\\bar x\_{\cdot j} + \zeta_j +
\varepsilon\_{ij}^{(m)},\quad \zeta_j\sim N(0,\sigma^2\_{\text{sch}}),\\
\varepsilon\sim N(0,\sigma^2). \\

Here \\y\_{ij}^{(m)}\\ is the \\m\\-th plausible value for student
\\ij\\, \\\zeta_j\\ is a school random effect, and
\\\varepsilon\_{ij}^{(m)}\\ is the student-level residual.

The defining feature is the **within–between (Mundlak / hybrid)** split
of the covariate \\x\\ (Mundlak 1978). Rather than entering \\x\_{ij}\\
with a single slope, the model separates it into:

- a **within-cluster** deviation \\\big(x\_{ij}-\bar x\_{\cdot j}\big)\\
  with slope \\\beta^{W}\\ — how an outcome moves with a student’s
  standing *relative to their own school’s mean*; and
- a **between-cluster** mean \\\bar x\_{\cdot j}\\ with slope
  \\\beta^{B}\\ — how the outcome moves with the *school’s average*
  level of \\x\\.

This is what lets a single specification distinguish a student-level
association from a school-level one. The cluster means \\\bar x\_{\cdot
j}\\ are computed **once** from the full sample and reused across every
plausible value; they are **PV-invariant**, so they are not recomputed
for each \\m\\. That invariance is one of the regularity conditions
behind the stacking identity in M3.

Two scope facts, stated here and inherited by the whole track:

- **Reportable scope is \\\beta\_{\text{FE}}\\ only.** The fixed-effect
  block — the intercept and the slopes — is the calibrated, reportable
  output. The variance components \\(\sigma^2,\sigma^2\_{\text{sch}})\\
  are *fit* but **not** calibrated, and are deferred to later work.

**The fixture’s `b_x` is a single-covariate demo, not a within/between
pair.** EQ-MODEL above shows the general hybrid parameterization with
both \\\beta^{W}\\ and \\\beta^{B}\\. The bundled `pisa_tiny` fit,
however, uses the simpler formula `OUTCOME ~ x + female` and reports
fixed-effect terms `b_Intercept`, `b_x`, `b_female`. That single `b_x`
slope is the *simplest* instance of EQ-MODEL — it does **not** exercise
the explicit within/between split. Read `b_x` as a plain
single-covariate stand-in, not as a \\(\beta^{W},\beta^{B})\\ pair.

## 4. The fixed-effect estimand and the per-PV likelihood

The estimand this package calibrates and reports is the fixed-effect
block \\\beta\_{\text{FE}}\\ from EQ-MODEL — and only that block.
Everything in M2–M4 is in service of getting \\\beta\_{\text{FE}}\\ and
its uncertainty right.

Working one plausible value at a time gives the per-PV building blocks.
Fitting EQ-MODEL to plausible value \\m\\ yields a fixed-effect estimate
\\\hat\beta_m\\ and an estimated covariance \\\hat U_m\\. Combining
these across \\m=1,\dots,M\\ with the Rubin rules produces the
reportable quantities — the Rubin mean \\\bar\beta\\, the within- and
between-imputation covariances \\\bar U\\ and \\B\\, and the total
\\T\_{\text{MI}}\\ — which **M2** assembles as EQ-TMI, with the
design-based \\\hat U_m\\ coming from the BRR–Fay sandwich (EQ-BRRFAY).

There is one more object to name now, because the stacked machinery
downstream is built on it. For each plausible value \\m\\, write
\\L_m^{0}(\psi)\\ for the **unweighted per-PV likelihood** — the
likelihood of the parameters \\\psi\\ under EQ-MODEL fit to plausible
value \\m\\ alone, with no survey weighting applied. M1 only *names* it.
Its role is what matters for the rest of the track: **M3** shows that
stacking all \\M\\ plausible values into one model with a \\1/M\\ weight
per row makes the stacked objective the *average* of these per-PV
likelihoods (the bridge EQ-BRIDGE), and that the resulting fixed-effect
estimate coincides with the Rubin mean \\\bar\beta\\ under stated
conditions (Theorem 2.2 / EQ-THM22). In other words, \\L_m^{0}\\ is the
term the stacked bridge averages — keep it in view as you read M3.

Two pointers, no derivations here: the *variance* answer (not just the
point) comes from the external target of **M2**, mapped onto the stacked
fit by the calibration of **M4** (EQ-CCC); and which fits may carry a
coverage claim is the subject of **M5**.

## 5. Symbol table & equation inventory at a glance

This section is the reference card. The first table is the symbol lock;
the second is the catalogue of tagged equations. **M1 prints only
EQ-MODEL (Section 3) in full** — every other equation is listed here by
tag and derived on the page named in the last column. This is the
“lookup, not a derivation” promise made concrete. The math typeset in
the cells renders via MathJax.

### 5.1 Symbol table

| Symbol | Meaning |
|----|----|
| \\i,j\\ | student \\i\\ in school/cluster \\j\\ |
| \\M\\ | number of plausible values (PVs); PISA reading \\M=10\\; **fixture \\M=2\\** |
| \\R\\ | number of BRR replicate weights; PISA \\R=80\\; **fixture \\R=4\\** |
| \\k\\ | Fay coefficient; PISA \\k=0.5\\; **fixture \\k=0.5\\** |
| \\a_d\\ | BRR–Fay multiplier \\=\dfrac{1}{R(1-k)^2}\\; PISA \\\Rightarrow a_d=0.05\\; **fixture \\\Rightarrow a_d=1\\** (arithmetic coincidence of \\R{=}4,k{=}0.5\\; matches `fay_variance_multiplier == 1`; **do not generalize**) |
| \\y\_{ij}^{(m)}\\ | the \\m\\-th plausible value for student \\ij\\ |
| \\\beta\_{\text{FE}}\\ | fixed-effect parameter block (the **only** calibrated/reportable block) |
| \\\beta^{W},\beta^{B}\\ | within- and between-cluster slopes (hybrid/Mundlak parameterization) |
| \\W\_{ij}\\, \\W\_{ij}^{(r)}\\ | final survey weight and \\r\\-th replicate weight |
| \\\hat\beta_m,\hat U_m\\ | per-PV fixed-effect estimate and its covariance |
| \\\bar\beta,\bar U,B\\ | Rubin mean, within-, and between-imputation covariance |
| \\T\_{\text{MI}}\\ | Rubin total (the external target covariance) |
| \\\nu\_{\text{BR}}\\ | Barnard–Rubin degrees of freedom |
| \\\gamma_k\\ | fraction of missing information for coefficient \\k\\ |
| \\\hat k\\ (or \\\hat\kappa\\) | PSIS Pareto tail-shape diagnostic (`stack_psis`) |

A few **auxiliary symbols** appear inside individual equations
downstream; they are listed here once so M2–M5 may cite them without
reintroduction:

| Symbol | Meaning | First used |
|----|----|----|
| \\\hat\beta_m^{(r)}\\ | per-PV estimate on the \\r\\-th replicate weight | EQ-BRRFAY (M2) |
| \\L,\\ L L^\top=\Sigma\\ | lower-triangular Cholesky factor of a covariance \\\Sigma\\ | EQ-CCC (M4) |
| \\L\_{\text{raw}},L\_{\text{tgt}}\\ | Cholesky factors of \\\Sigma\_{\text{raw}}\\ (raw draw covariance) and \\\Sigma\_{\text{tgt}}=T\_{\text{MI}}\\ (target) | EQ-CCC (M4) |
| \\\beta_s,\\ \beta_s^{\text{cal}}\\ | the \\s\\-th raw posterior draw and its CCC-calibrated image | EQ-CCC (M4) |
| \\c,\\ \bar\beta^{\text{raw}}\\ | target center the calibrated draws are shifted to (\\c=\bar\beta\\ for reporting), and the mean of the raw draw cloud | EQ-CCC (M4) |
| \\\Delta_c\\ | center-separation diagnostic (raw center vs \\\bar\beta\\, on the target-SE scale) | EQ-DELTAC (M4) |
| \\L_m^{0}(\psi),\\ q\_{\text{SWL}}(\psi)\\ | unweighted per-PV likelihood for PV \\m\\, and the stacked fractional objective | EQ-BRIDGE (M3) |
| \\\hat\beta\_{\text{FE}}^{q}\\ | fixed-effect estimate from maximizing the stacked objective | EQ-THM22 (M3) |
| \\B\_{kk},\\T\_{\text{MI},kk}\\ | the \\k\\-th diagonal entries of \\B\\ and \\T\_{\text{MI}}\\ | EQ-FMI (M2) |
| \\N\\ | number of students (the stacked design has \\N\cdot M\\ rows) | EQ-BRIDGE (M3) |

### 5.2 Equation inventory

Each row names a tagged equation, glosses it in one plain-English line,
and points to the vignette that derives it. Only EQ-MODEL is displayed
in full on this page (Section 3).

| Tag | What it says (one line) | Derived in |
|----|----|:--:|
| **EQ-MODEL** | the two-level within–between (hybrid) model for each plausible value | **M1** (Section 3) |
| **EQ-TMI** | Rubin combining: \\\bar\beta\\, \\\bar U\\, \\B\\, and the external total \\T\_{\text{MI}}=\bar U+(1+1/M)B\\ | **M2** |
| **EQ-BRRFAY** | the BRR–Fay replicate-weight sandwich for each per-PV covariance \\\hat U_m\\ | **M2** |
| **EQ-BARNARD** | Barnard–Rubin small-sample degrees of freedom \\\nu\_{\text{BR}}\\ | **M2** |
| **EQ-FMI** | fraction of missing information \\\gamma_k\\ for each coefficient | **M2** |
| **EQ-BRIDGE** | the stacked fractional objective: one fit on \\N\cdot M\\ rows, weight \\1/M\\ per row, averaging the \\L_m^{0}\\ | **M3** |
| **EQ-THM22** | Theorem 2.2: under regularity R0–R5 the stacked fixed-effect estimate equals \\\bar\beta\\ | **M3** |
| **EQ-CCC** | the Cholesky Calibration Correction — an affine map matching the draw cloud’s first two moments to the target | **M4** |
| **EQ-DELTAC** | the center-separation diagnostic \\\Delta_c\\ that gates the calibrated fit | **M4** |
| **EQ-PSIS** | Pareto-smoothed importance-sampling tail-shape \\\hat k\\ and its good/borderline/unreliable thresholds | **M5** |

**How to read a tag downstream.** When M3 writes “by EQ-THM22”, it means
the identity in the row above, derived in M3. M1 does not prove any of
these; it only fixes their names and symbols so the later pages need
not.

## 6. Orientation panel

Now we attach the symbols to something concrete. The package ships a
tiny **synthetic** PISA-shaped table, `pisa_tiny`, and a cached
`stack_direct` fit on it, so that every example runs offline and
deterministically. The panel below reads only — it detects the design
and reads the fitted target’s metadata. **Nothing is estimated.**

First, load the raw table and see its shape:

``` r

pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

dim(pisa_tiny)        # 12 rows x 12 columns
#> [1] 12 12
names(pisa_tiny)
#>  [1] "CNT"        "CNTSCHID"   "CNTSTUID"   "x"          "female"    
#>  [6] "PV1READ"    "PV2READ"    "W_FSTUWT"   "W_FSTURWT1" "W_FSTURWT2"
#> [11] "W_FSTURWT3" "W_FSTURWT4"
```

It is one make-believe country with 12 students: a continuous covariate
`x`, a binary `female` indicator, **two** plausible reading values
`PV1READ` / `PV2READ`, **four** BRR replicate weights
`W_FSTURWT1`–`W_FSTURWT4`, a final weight `W_FSTUWT`, and student /
school identifiers.

The detectors recover the design constants from the column names. The
plausible values carry a `READ` suffix, so you must pass
`suffix = "READ"`:

``` r

detect_pisa_pv_columns(pisa_tiny, suffix = "READ")   # -> c("PV1READ","PV2READ")  => M = 2
#> [1] "PV1READ" "PV2READ"
detect_pisa_brr_replicate_weights(pisa_tiny)         # -> 4 W_FSTURWT* cols        => R = 4
#> [1] "W_FSTURWT1" "W_FSTURWT2" "W_FSTURWT3" "W_FSTURWT4"
```

**The bare detector call errors on this fixture.** Calling
`detect_pisa_pv_columns(pisa_tiny)` *without* the suffix raises “No
plausible-value columns detected …”, because the default looks for
unsuffixed `PV` columns and this fixture’s columns are subject-suffixed
(`PV1READ`, `PV2READ`). You must pass `suffix = "READ"`, as above.
Modern PISA files are suffixed this way, so this is the common case, not
an edge case.

Finally, read the design straight off the cached fit’s external target.
The same two constants (\\M=2\\, \\R=4\\) reappear, together with the
Fay coefficient and the fixed-effect term names:

``` r

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)

c(M = tg$M, R = tg$R, fay_k = tg$fay_k)   # 2, 4, 0.5
#>     M     R fay_k 
#>   2.0   4.0   0.5
tg$fe_names                               # "b_Intercept" "b_x" "b_female"
#> [1] "b_Intercept" "b_x"         "b_female"
```

So the fixture realises \\M=2\\, \\R=4\\, \\k=0.5\\ — and therefore the
BRR–Fay multiplier \\a_d=1/\[4\\(1-0.5)^2\]=1\\, exactly the
`fay_variance_multiplier` field. The three fixed-effect terms are
`b_Intercept`, `b_x`, and `b_female`, the reportable
\\\beta\_{\text{FE}}\\ block from Section 4.

**Honesty notes for this fixture (carried by the whole track).**

- **It is synthetic and illustrative.** One make-believe country
  (`CNT == "SYN"`), 12 students, \\M=2\\, \\R=4\\. The numbers are
  included only so examples and package checks run without licensed data
  — they are **not** real PISA results and must never be cited as such.
  The package is **not** OECD-affiliated and bundles no real PISA.
- **`b_x` is a single-covariate demo**, not a within/between
  \\(\beta^{W}, \beta^{B})\\ pair (Section 3).
- **\\a_d=1\\ is the tiny-design coincidence** of \\R=4,k=0.5\\ —
  distinct from PISA’s \\a_d=0.05\\. Do not generalise it (Section 2.1).
- **Fixed-effect-only scope.** Variance components are fit but not
  calibrated and are not part of the reportable output (Section 3).

## 7. Where to next

You now hold the vocabulary the rest of the Method track cites by tag.
From here:

- **M2 · The BRR–Fay fixed-effect target** — the natural next step: it
  derives EQ-TMI, EQ-BRRFAY, EQ-BARNARD, and EQ-FMI, building the
  external target that \\\beta\_{\text{FE}}\\ is calibrated against, and
  re-derives Rubin pooling on the cached target in code.
- Then **M3** (the stacked fractional bridge and Theorem 2.2), **M4**
  (CCC and the center-separation diagnostics), and **M5** (the three
  methods, PSIS, and what makes a fit coverage-claimable) complete the
  provenance chain.

To see the workflow these equations sit underneath, the Applied track is
the place: **A1 · Getting started** (load a fit and read it three ways)
and **A2 · The end-to-end workflow** (declare the design, assemble the
target, run the fit, report).

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

Davier, Matthias von, Eugenio Gonzalez, and Robert J. Mislevy. 2009.
“What Are Plausible Values and Why Are They Useful?” *IERI Monograph
Series* 2: 9–36.

Mislevy, Robert J. 1991. “Randomization-Based Inference about Latent
Variables from Complex Samples.” *Psychometrika* 56 (2): 177–96.
<https://doi.org/10.1007/BF02294457>.

Mundlak, Yair. 1978. “On the Pooling of Time Series and Cross Section
Data.” *Econometrica* 46 (1): 69–85. <https://doi.org/10.2307/1913646>.

OECD. 2024. *PISA 2022 Technical Report*. OECD Publishing.

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.
