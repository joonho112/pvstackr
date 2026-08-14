# Construct pvstackr Fitting Controls

`pv_control()` validates package-level options used by the fitting,
calibration, diagnostics, and object-retention layers, returning a
frozen `pvstackr_control` object that every `pv_fit*()` entry point
consumes. It centralizes backend policy, interval level, the calibration
centering convention, and how heavy a fitted object is allowed to be.

## Usage

``` r
pv_control(
  method = "stack_direct",
  chains = 4L,
  iter = 2000L,
  warmup = NULL,
  cores = 1L,
  seed = NULL,
  backend = "none",
  conf_level = 0.95,
  psis_k_threshold = 0.7,
  center = "target",
  allow_target_nearpd = FALSE,
  return_draws = TRUE,
  keep_data = FALSE,
  keep_backend_fit = FALSE,
  keep_log_lik = FALSE,
  verbose = FALSE
)

# S3 method for class 'pvstackr_control'
print(x, ...)
```

## Arguments

- method:

  Public method identifier. Character scalar; must be one of
  `"stack_direct"`, `"stack_psis"`, or `"per_pv"`. Default
  `"stack_direct"`. A control's `method` must equal the `method` passed
  to
  [`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md)
  and to the dispatched fitter.

- chains:

  Number of MCMC chains requested by a live backend. Integer-valued
  scalar, `>= 1`.

- iter:

  Total iterations per chain. Integer-valued scalar, `>= 2`.

- warmup:

  Warmup iterations per chain. Integer-valued scalar, `>= 0` and
  strictly less than `iter`. If `NULL`, defaults to `floor(iter / 2)`.

- cores:

  Number of cores requested by a live backend. Integer-valued scalar,
  `>= 1`.

- seed:

  Optional random seed. Integer-valued scalar `>= 0`, or `NULL` for no
  fixed seed.

- backend:

  Backend policy. Character scalar; one of `"none"`, `"injected"`,
  `"brms"`, or `"cmdstanr"`. Default `"none"`. In this package stage
  `"brms"` selects the bundled brms adapter (with deterministic
  cmdstanr/rstan resolution), while other live engines use an injected
  `fit_function` adapter or precomputed draws.

- conf_level:

  Confidence or credible-interval level for report tables. Numeric
  scalar in `(0, 1)`. Default `0.95`.

- psis_k_threshold:

  Pareto-k reportability threshold for `stack_psis`. Numeric scalar in
  `(0, 0.7]`. Default `0.7`; values below `0.7` may impose a stricter
  gate, but the package-level ceiling cannot be relaxed. Every plausible
  value must have a finite Pareto-k strictly below this threshold.

- center:

  Calibration centering convention, `"target"` or `"posterior"`.
  Reportable `stack_direct` output requires `"target"`. `"posterior"` is
  reserved for CCC diagnostic/exploratory checks: it leaves fixed-effect
  draws at the raw stacked posterior mean while still computing
  target-covariance diagnostics against the external target. Default
  `"target"`.

- allow_target_nearpd:

  Reserved for future target-covariance repair. Logical scalar; must be
  `FALSE`. Automatic target repair is not currently supported. Default
  `FALSE`.

- return_draws:

  Logical scalar. Whether fitted objects retain their method-specific
  fixed-effect draw representation. For `stack_direct`, read the
  calibrated matrix via
  [`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md).
  Per-PV matrices and the PSIS proposal/weight pair remain in
  [`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md).
  Default `TRUE`.

- keep_data:

  Logical scalar. Whether fitted objects may retain the user data frame.
  Default `FALSE` (fits stay light).

- keep_backend_fit:

  Logical scalar. Whether fitted objects may retain the heavy, opaque
  backend fit object. This requires `keep_data = TRUE` on public
  composite fits because a backend object may contain analysis data.
  Default `FALSE` (fits stay light).

- keep_log_lik:

  Logical scalar. Whether fitted objects may retain log-likelihood
  draws. Default `FALSE` (fits stay light).

- verbose:

  Logical scalar. Whether functions emit progress messages. Default
  `FALSE`.

- x:

  A `pvstackr_control` object.

- ...:

  Ignored.

## Value

A validated `pvstackr_control` object: a named list of the resolved
options above, with class `c("pvstackr_control", "list")`. Pass it to a
`pv_fit*()` function via the `control` argument; print it for a compact
summary.

## Details

### Centering and reportable output

The `center` convention decides whether a fit is *reportable* or merely
*diagnostic*. Reportable `stack_direct` output requires
`center = "target"`: fixed-effect draws are CCC-calibrated so their mean
and covariance match the external Rubin/BRR-Fay target, and the estimate
table inherits the target's standard errors, degrees of freedom, and
interval metadata. `center = "posterior"` is diagnostic and exploratory
only: it leaves fixed-effect draws at the raw stacked posterior mean
while still computing target-covariance diagnostics against the external
target.
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md)
rejects a control with `center = "posterior"` rather than emit a
non-reportable `stack_direct` estimate table; treat `"posterior"` as a
CCC check, never a deliverable.

### Object retention

The retention flags govern how much a fitted `pvstackr_fit` carries, and
the defaults keep fits light. `return_draws` (default `TRUE`) retains
only the method-specific fixed-effect draw representation: calibrated
top-level draws for `stack_direct`, fixed-effect-only per-PV matrices in
diagnostics for `per_pv`, and a fixed-effect proposal matrix paired with
normalized PV weights in diagnostics for `stack_psis`.
[`get_draws()`](https://joonho112.github.io/pvstackr/reference/get_draws.md)
exposes only the synthesized top-level `stack_direct` matrix; use
[`get_diagnostics()`](https://joonho112.github.io/pvstackr/reference/get_diagnostics.md)
for the other two representations. Package-owned raw full stacked draws,
nuisance draws, and duplicate calibrated matrices are never retained in
a final composite fit. This does not inspect or override explicitly
authorized opaque backend objects (`keep_backend_fit`) or log-likelihood
matrices (`keep_log_lik`), which have separate retention authorities.
Blocked fits fail closed and record effective `return_draws`,
`keep_data`, `keep_backend_fit`, and `keep_log_lik` values of `FALSE`,
irrespective of the original retention request. `keep_data` retains the
user data frame in the top-level design; otherwise the fit stores a
metadata/hash-only design snapshot with a base-environment formula.
Package-owned stacked-data copies are not retained inside composite
fits. `keep_backend_fit` retains an opaque backend object and therefore
requires `keep_data = TRUE`, because the package cannot prove that an
arbitrary backend object is data-free. `keep_log_lik` retains
log-likelihood draws. Enable the `FALSE`-by-default flags only when you
need the extra payload (for example, re-extraction or model checking),
as each materially increases the size of the saved object.

## See also

[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)

Other pvstackr-fitting:
[`pv_fit()`](https://joonho112.github.io/pvstackr/reference/pv_fit.md),
[`pv_fit_direct()`](https://joonho112.github.io/pvstackr/reference/pv_fit_direct.md),
[`pv_fit_reference()`](https://joonho112.github.io/pvstackr/reference/pv_fit_reference.md),
[`pv_fit_stack_psis()`](https://joonho112.github.io/pvstackr/reference/pv_fit_stack_psis.md)

## Examples

``` r
# Default controls target reportable stack_direct output.
ctrl <- pv_control()
ctrl
#> pvstackr control
#>   method: stack_direct
#>   backend: none
#>   iter/warmup/chains: 2000/1000/4
#>   target repair: unsupported (disabled)

# A control's method must match the method you fit with.
ctrl_psis <- pv_control(method = "stack_psis", psis_k_threshold = 0.7)
ctrl_psis$method
#> [1] "stack_psis"

# Retain the heavy backend fit and log-likelihood draws when you need them.
ctrl_heavy <- pv_control(keep_backend_fit = TRUE, keep_log_lik = TRUE)
c(ctrl_heavy$keep_backend_fit, ctrl_heavy$keep_log_lik)
#> [1] TRUE TRUE
```
