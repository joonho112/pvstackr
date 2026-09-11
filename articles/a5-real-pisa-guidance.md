# Real PISA data: access, memory and reproducibility

Abstract

This article describes what to do when you analyze real PISA files with
pvstackr: obtaining the files under their terms, declaring the PISA 2022
design, writing the target and fit calls, planning for the memory and
computing time they need, and recording what is needed to repeat the
analysis. pvstackr ships no real PISA records, so the calls are shown
without being run. The last section says what the companion preprint
shows about PISA and what it does not.

``` r

library(pvstackr)
```

pvstackr does not ship real PISA records. The only data in the package
are the synthetic `pisa_tiny` data (one made-up country, `CNT` is
`"SYN"`, with \\M = 2\\ plausible values and \\R = 4\\ replicate
weights) and an example fit made from them. They exist so that the
examples and the package tests run without licensed files. pvstackr is
not affiliated with or endorsed by OECD/PISA. Its documentation uses the
PISA column names, such as `PV1MATH`, `W_FSTUWT`, `W_FSTURWT1` and
`CNTSTUID`, because they are the names in the files that users analyze
(OECD 2024).

Each call below is written for a PISA 2022 file and is not run, since
the package has no real data; [The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.md)
runs the same steps on the synthetic data.

## Data access and licensing

### Obtaining the files

Obtain real PISA files only from an authorized source: the OECD’s
official PISA website or a copy that your institution has approved.
Before you use a file, read the license, citation, redistribution, and
privacy terms of the cycle and of the files you downloaded. The terms
differ between cycles and releases, and the terms that come with the
files are the ones that apply.

### The pvstackr license and your data

The MIT license for `pvstackr` does not license users’ real PISA data,
the extracts derived from them or any other files processed with the
package. It covers the software only. Your data are governed by the
terms of their source, including its rules on attribution,
redistribution, third-party rights and privacy, and following them is
your responsibility.

### Keeping the data out of your repository

Do not commit real PISA microdata, large extracts or the cache files of
a modeling engine to a package or project repository under version
control. Keep them in a folder outside the source tree and refer to the
folder through an environment variable, so that its path does not appear
in your code:

``` r

pisa_dir <- Sys.getenv("PISA_DATA_DIR")
stopifnot(nzchar(pisa_dir))

student_file <- file.path(pisa_dir, "student_file_for_your_cycle")
```

[`Sys.getenv()`](https://rdrr.io/r/base/Sys.getenv.html) returns an
empty string when `PISA_DATA_DIR` is not set, and
[`stopifnot()`](https://rdrr.io/r/base/stopifnot.html) then stops the
script with an error. You can set the variable in your `.Renviron` file.

### Reading the file

File names and formats differ between cycles and releases. Read your
authorized file with the reader for its format, whether an SPSS, SAS,
transport, or converted local extract, and keep only the columns that
the model needs: the country, the identifiers, the covariates, the
plausible values, the final weight and the replicate weights. Keeping
only these columns reduces the memory that the later steps need (see
[memory and runtime](#memory-runtime)).

``` r

# An example: use the reader for the format of your file, such as
# haven::read_sav() for SPSS, a reader for SAS or SAS transport files, or
# read.csv() for a converted CSV file. pvstackr has no PISA reader.
pisa_raw <- haven::read_sav(student_file)

pisa_country <- subset(pisa_raw, CNT == "USA")
pisa_country <- pisa_country[
  ,
  c(
    "CNT",
    "CNTSCHID",
    "CNTSTUID",
    "ESCS",
    "ST004D01T",
    paste0("PV", 1:10, "MATH"),   # M = 10 plausible values in mathematics
    "W_FSTUWT",                    # final student weight
    paste0("W_FSTURWT", 1:80)     # R = 80 replicate weights
  )
]
```

The extract keeps the country code `CNT`, the school and student IDs
`CNTSCHID` and `CNTSTUID`, two covariates (`ESCS`, the PISA index of
economic, social and cultural status, and `ST004D01T`, the student’s
gender (OECD 2024, ch. 19)), the ten plausible values in mathematics,
the final weight `W_FSTUWT` and the 80 replicate weights.

## Declaring the design

[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
finds and checks the plausible-value and weight columns before anything
is computed, as described in [The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#check-columns).
On a real file, check each argument that describes the design against
the technical documentation of the cycle: the plausible-value suffix
`pv_suffix`, the number of plausible values `expected_M`, the number of
replicate weights `expected_R`, the Fay coefficient `fay_k` and the
identifier columns `id_cols`. In PISA 2022, \\M = 10\\, \\R = 80\\ and
\\k = 0.5\\ (OECD 2024).
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
stops with an error when the number of plausible-value or
replicate-weight columns it finds differs from `expected_M` or
`expected_R`, for example because the suffix is wrong. It does not check
`fay_k` against the replicate weights, so take that value from the
documentation.

``` r

design <- pv_design(
  data       = pisa_country,
  formula    = OUTCOME ~ ESCS + ST004D01T,   # OUTCOME stands for each plausible value
  pv_suffix  = "MATH",                        # finds PV1MATH to PV10MATH
  expected_M = 10L,                           # stop unless 10 plausible values are found
  expected_R = 80L,                           # stop unless 80 replicate weights are found
  fay_k      = 0.5,                           # the Fay coefficient of PISA 2022
  id_cols    = "CNTSTUID"                     # the column that identifies each row
)

design
```

The formula has the placeholder `OUTCOME` on the left-hand side, which
pvstackr replaces by each plausible value in turn, `PV1MATH` to
`PV10MATH`. `pv_suffix = "MATH"` selects the mathematics values, since
PISA 2022 names the plausible values by subject. The default,
`pv_suffix = ""`, matches only bare names such as `PV1`; on a PISA 2022
file it finds no plausible-value columns and stops with an error. The
survey weights stay out of the formula; the defaults of
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
find the final weight `W_FSTUWT` and the replicate weights `W_FSTURWT1`
to `W_FSTURWT80` by their PISA names.

[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md)
checks only the form of the formula. The covariates are checked by
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md),
which stops with an error when one of them has a missing value, is a
character column (make it a factor) or has a class that pvstackr does
not accept. `haven::read_sav()` returns the columns that carry SPSS
value labels with the class `haven_labelled`, and such a covariate stops
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
with the message “A formula predictor has an unsupported type or class”.
Since the target is computed from `design$data`, prepare the extract
before you call
[`pv_design()`](https://joonho112.github.io/pvstackr/reference/pv_design.md):
remove the incomplete rows, for example with
[`complete.cases()`](https://rdrr.io/r/stats/complete.cases.html), and
convert the labeled columns to plain numbers, for example with
`haven::zap_labels()`.

## Computing the target

[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
computes the design-based target from the columns in the design. For
each plausible value it fits the regression by weighted least squares,
once with the final weight and once with each replicate weight, and
computes the BRR–Fay replicate covariance of the coefficients from the
spread of the replicate estimates. Rubin’s rules then combine the \\M\\
results into the estimates \\\bar\beta\\ and their covariance matrix
\\T\_{\text{MI}} = \bar U + (1 + 1/M)\\B\\, where \\\bar U\\ is the
average of the \\M\\ replicate covariances and \\B\\ is the covariance
matrix of the \\M\\ estimates (Rubin 1987). [The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#target)
gives the formulas. The target is fixed-effect-only:
[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
stops with an error when the formula contains a random-effect term such
as `(1 | CNTSCHID)`, and so does a `stack_direct` fit.

``` r

target <- pv_brr_target(
  data            = design$data,
  formula         = design$formula,
  pv_cols         = design$pv_cols,
  weight_col      = design$weight_col,
  rep_weight_cols = design$rep_weight_cols,
  fay_k           = design$fay_k,
  id_cols         = design$id_cols
)

target
```

The call uses the default rule for the degrees of freedom,
`df_method = "classic"`, so the intervals of a fit calibrated to this
target are descriptive (`coverage_claim_allowed = FALSE`). By pvstackr’s
reporting rule, only a target made with `df_method = "barnard_rubin"`
and a positive `df_complete` gives `stack_direct` intervals that can be
read as confidence intervals with nominal coverage. [Reading and
reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#interval-columns)
explains the rule and `df_complete`.

## Fitting the model

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
with `method = "stack_direct"` fits one Bayesian model to the stacked
data, which hold one copy of the data for each plausible value, and
calibrates its fixed-effect draws to the target with the Cholesky
calibration correction (CCC). The fit needs an engine: the bundled one,
selected with `pv_control(backend = "brms")`, or three functions of your
own, `fit_function`, `draws_function` and `diagnose_function`. The chunk
below uses functions of your own. With them, leave `backend` at its
default: pvstackr passes the value to your `fit_function` and records it
in the fit, but the value selects nothing. [The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#own-engine)
describes what the three functions receive and must return; without
`diagnose_function` the fit is blocked.

``` r

fit <- pv_fit(
  data           = design$data,
  formula        = design$formula,
  target         = target,
  method         = "stack_direct",
  control        = pv_control(
                     method           = "stack_direct",
                     center           = "target",     # the default; stack_direct requires it
                     iter             = 2000L,
                     warmup           = 1000L,
                     chains           = 4L,
                     cores            = 4L,
                     seed             = 20260607,
                     return_draws     = FALSE,
                     keep_data        = FALSE,
                     keep_backend_fit = FALSE,
                     keep_log_lik     = FALSE
                   ),
  fit_function      = your_fit_function,      # fits the stacked model
  draws_function    = your_draws_function,    # returns the draws of that fit
  diagnose_function = your_diagnose_function, # returns its sampler diagnostics
  cache_dir         = NULL                    # no cache file for your function
)

summary(fit)
get_estimates(fit)
```

To use the bundled engine instead, leave out the three functions and set
`backend = "brms"` in
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md).
It fits the stacked model with the brms function `brm()` and needs the
brms and posterior packages. It saves the brms fit, which contains the
stacked outcome and weights, as a cache file in the folder `cache_dir`
(by default `cache` in the working directory), creating the folder if
needed. On real data, set `cache_dir` to a folder outside your
repository, for example under `PISA_DATA_DIR`. brms reuses a cached fit
when only `chains`, `iter`, `warmup` or `seed` change, so give the cache
file a new name with `cache_stem`, or set `cache_dir = NULL`, when you
change them (see
[`?pv_fit_direct`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)).

### Reading the result

`get_estimates(fit)` returns the estimate table, which you read as in
[Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#estimate-table).
With `control$center = "target"`, the default, which `stack_direct`
requires, its numbers are copied from the target ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#center-target)).
The other value, `center = "posterior"`, is accepted by
[`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md),
but no fitting function uses it:
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
stops with an error when a `stack_direct` fit sets it. The diagnostics
include the center separation, which shows how far the mean of the
stacked draws was from the target before the calibration ([Reading and
reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#target-and-diagnostics)).

## Memory and runtime

Start with one country, one domain and a small formula, and keep only
the columns that the analysis needs (see [data access and
licensing](#data-access)). Two counts then tell you the size of a job
before you start it: the amount of data and the number of fits.

The stacked data hold one copy of the data for each plausible value, so
a design with `N` rows and `M` plausible values gives `N * M` stacked
rows. Each copy keeps every column of the data, not only those in the
formula, so a column that you drop from the extract is not repeated `M`
times. Row \\i\\ of each copy has the weight \\\tilde w_i/M\\, where
\\\tilde w_i\\ is the student’s final weight divided by the mean final
weight ([The full analysis
workflow](https://joonho112.github.io/pvstackr/articles/a2-the-workflow.html#fitting)
explains why). With \\M = 10\\, as in PISA 2022, the engine fits the
model to ten times as many rows as the extract has.

[`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
fits the regression by weighted least squares once with the final weight
and once with each of the `R` replicate weights, for each of the `M`
plausible values: `M * (R + 1)` fits. For PISA 2022 (\\M = 10\\, \\R =
80\\) that is `10 * 81 = 810` weighted least-squares fits, before any
Bayesian sampling. Each is an ordinary weighted regression, but on a
large national extract it is worth trying the call on a smaller part of
the data first.

The Bayesian fit usually needs the most memory and time. These settings
help:

- Set `seed` and record every setting of
  [`pv_control()`](https://joonho112.github.io/pvstackr/reference/pv_control.md)
  (see the [reproducibility checklist](#reproducibility-checklist)), so
  that the fit can be repeated.
- Keep the engine’s cache files outside the repository, with the data
  (see [data access and licensing](#data-access)).
- Keep the defaults `keep_data = FALSE`, `keep_backend_fit = FALSE` and
  `keep_log_lik = FALSE` unless you need the data, the engine’s fit
  object or the log-likelihood draws in the saved fit; each makes the
  saved fit larger.
- Set `return_draws = FALSE` when you do not need the calibrated draws:
  the estimate table comes from the target and does not use them
  ([Reading and reporting the
  results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#target-and-diagnostics)).
- Increase `chains`, `iter` and `cores` only once the design and the
  target are final, that is, when their checksums `design_hash` and
  `target_hash` no longer change from one run to the next; a long run on
  inputs that still change has to be repeated.

## Reproducibility checklist

Record these items with each analysis of real data, so that you and your
readers can repeat it and see which inputs produced which numbers:

- the PISA cycle, the release date and the authorized source of the
  files;
- the filters applied to the extract, such as the country and the
  domain;
- the model variables and the plausible-value suffix `pv_suffix`;
- the final weight, the prefix and number of the replicate weights, and
  the Fay coefficient `fay_k`;
- the pvstackr version, `utils::packageVersion("pvstackr")`, and the R
  session, printed by
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html);
- the checksums of the design, `design_hash`, `row_support_hash`,
  `pv_value_hash` and `weight_design_hash`, and of the target,
  `target_hash`;
- the engine (`backend`), the seed, and `chains`, `iter`, `warmup` and
  `cores`;
- what the fit keeps: `return_draws`, `keep_backend_fit`, `keep_data`
  and `keep_log_lik`.

The design checksums do not cover the values of the covariates.
`target_hash` is a checksum of the contents of the target, including its
estimates, so a change in the covariates changes it as well; it is also
a column of the estimate table ([Reading and reporting the
results](https://joonho112.github.io/pvstackr/articles/a3-reading-results.html#estimate-table)).
The chunk below collects the records in one list from the objects made
above:

``` r

# A sketch: one list with the records above, taken from the objects made in
# this article. Rename the fields to suit your project.
analysis_record <- list(
  pisa_cycle         = "2022",                  # record the release date too
  data_source        = "OECD PISA portal",      # the authorized source
  filters            = "CNT == 'USA'; domain = MATH",
  variables          = all.vars(design$formula),
  pv_suffix          = "MATH",
  weight_col         = design$weight_col,
  rep_weight_prefix  = "W_FSTURWT",
  R                  = length(design$rep_weight_cols),
  fay_k              = design$fay_k,
  pkg_version        = as.character(utils::packageVersion("pvstackr")),
  design_hash        = design$design_hash,
  row_support_hash   = design$row_support_hash,
  pv_value_hash      = design$pv_value_hash,
  weight_design_hash = design$weight_design_hash,
  target_hash        = target$target_hash,
  backend            = "none",                  # "brms" with the bundled engine
  seed               = 20260607,
  chains             = 4L,
  iter               = 2000L,
  warmup             = 1000L,
  cores              = 4L,
  return_draws       = FALSE,
  keep_backend_fit   = FALSE,
  keep_data          = FALSE,
  keep_log_lik       = FALSE
)

str(analysis_record)
```

[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) prints the
versions of R, of the operating system and of the loaded packages. Here
it runs and shows the session that built this article, which does not
attach pvstackr; in your analysis, call it after
[`library(pvstackr)`](https://joonho112.github.io/pvstackr/) so that the
pvstackr version is recorded:

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
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.60         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  
#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.3.0       fs_2.1.0
```

## The companion preprint

The stacked fit and the calibration that pvstackr uses are described in
a companion preprint (Lee et al. 2026). The preprint gives the
conditions under which the fixed-effect estimate of one stacked fit
equals the average of the estimates for the separate plausible values:
the stacked fixed-effect point identity, Theorem 4.1 of the preprint,
which [One stacked fit and the Rubin point
estimate](https://joonho112.github.io/pvstackr/articles/m3-stacked-bridge.md)
explains. It also shows that CCC gives the calibrated draws exactly the
mean and covariance of the target (Proposition B.1 of the preprint), as
described in [Calibrating the fixed-effect draws
(CCC)](https://joonho112.github.io/pvstackr/articles/m4-ccc-calibration.md).
The preprint builds its default target from model-based fits of each
plausible value and describes a target from replicate weights, the only
kind that pvstackr computes, as an alternative ([The design-based
target: BRR–Fay replicate weights and Rubin’s
rules](https://joonho112.github.io/pvstackr/articles/m2-brr-fay-target.md)
describes pvstackr’s target).

The preprint’s application to PISA 2022 reading, in the United States
and Korea, checks whether one MCMC run, together with the target and one
calibration, reproduces the results of ten MCMC runs, one for each
plausible value (Lee et al. 2026, sec. 6). The preprint presents its
PISA values as context rather than as findings, and three points matter
when you read them:

- In the application, the stacked fit is calibrated to the combined
  output of the ten fits, so it has their estimates and standard errors
  by construction, as a `stack_direct` fit has those of its target. The
  preprint calls this agreement target-calibrated by design and does not
  offer it as independent evidence.
- The application fits a model with a random intercept for schools and
  takes its target from the Bayesian fits of that model to each
  plausible value. pvstackr 0.2.x cannot repeat this analysis: its
  target uses the BRR–Fay replicate weights, and
  [`pv_brr_target()`](https://joonho112.github.io/pvstackr/reference/pv_brr_target.md)
  and `stack_direct` stop with an error on a random-effect term.
- In the preprint, the per-PV workflow, one fit for each plausible value
  combined with Rubin’s rules, is the reference that defines the target.
  The package’s `per_pv` method combines posterior draws from fits that
  you supply, and its intervals are always descriptive, because the
  variance within each plausible value comes from the model’s posterior
  draws (model-based), not from the replicate weights.

The simulation of the preprint compares the per-PV workflows and the
calibrated stack with the generating values; their variances and targets
come from model-based fits, and the simulated data have no replicate
weights (Lee et al. 2026, sec. 5). It therefore does not test the
BRR–Fay target that pvstackr computes, and pvstackr’s reporting rule for
intervals is the package’s own rule, not a result of the preprint.

In a qualified secondary analysis (Appendix F.4 of the preprint), the
difference between Korea and the United States in the between-school
gradient of `ESCS` shrinks by more than half and loses nominal
significance when the survey weights are placed in two stages rather
than at the student level, but it keeps its sign. The preprint draws no
conclusion about the two countries from it.

[Comparing the three fitting
methods](https://joonho112.github.io/pvstackr/articles/a4-comparing-methods.html#cautions)
gives the cautions for comparing pvstackr’s methods on your own data,
including the preprint’s PISA result for the Pareto \\\hat k\\ check.

## References

Lee, JoonHo, Matthew R. Williams, and Terrance D. Savitsky. 2026. *One
Markov Chain Monte Carlo Fit for Many Plausible Values: A Calibrated
Stacked Posterior Workflow for Bayesian Multilevel Models of Large-Scale
Assessment Data*. Zenodo preprint, version 1.
<https://doi.org/10.5281/zenodo.22407935>.

OECD. 2024. *PISA 2022 Technical Report*. PISA. OECD Publishing.
<https://doi.org/10.1787/01820d6d-en>.

Rubin, Donald B. 1987. *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons. <https://doi.org/10.1002/9780470316696>.
