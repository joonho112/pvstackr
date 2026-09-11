## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----build-fits---------------------------------------------------------------
fe <- c("b_Intercept", "b_x", "b_female")

fit_direct <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit

# Made-up draws for the per_pv fit: for each plausible value, 300 random values
# per fixed effect, with standard deviation 2, around 458, 47 and 2, near the
# example estimates. They stand in for posterior draws; no model was fitted.
set.seed(1)
mk <- function(mu) {
  d <- matrix(rnorm(300 * 3, sd = 2), ncol = 3, dimnames = list(NULL, fe))
  sweep(d, 2L, mu, "+")
}

fit_per_pv <- pv_fit_reference(
  per_pv_draws = list(PV1READ = mk(c(458, 47, 2)),
                      PV2READ = mk(c(458, 47, 2))),
  control      = pv_control(method = "per_pv")
)

# Made-up stacked draws with equal placeholder weights and Pareto k-hat values
# of 0.2. No PSIS was run, so no program is named, and pvstackr blocks the fit
# (its PSIS status is provenance_incomplete).
sd_mat <- mk(c(458, 47, 2))
fit_psis <- pv_fit_stack_psis(
  stacked_draws = sd_mat,
  pv_cols       = c("PV1READ", "PV2READ"),
  psis_weights  = matrix(1 / nrow(sd_mat), nrow = nrow(sd_mat), ncol = 2),
  pareto_k      = c(0.2, 0.2),
  control       = pv_control(method = "stack_psis")
)

## ----compare------------------------------------------------------------------
cmp <- pv_compare_methods(
  stack_direct = fit_direct,
  per_pv       = fit_per_pv,
  stack_psis   = fit_psis
)

cmp

## ----aligned-estimates--------------------------------------------------------
est <- get_estimates(cmp)

est[, c("method_label", "term", "estimate", "se",
        "conf_low", "conf_high",
        "interval_role", "coverage_claim_allowed")]

## ----interval-widths----------------------------------------------------------
bx <- est[est$term == "b_x",
          c("method_label", "se", "conf_low", "conf_high")]
bx$width <- bx$conf_high - bx$conf_low
bx

## ----overlay-figure, fig.width = 7, fig.height = 3.4, fig.cap = "Estimates and 95% intervals of the two slopes, b_x and b_female, from the stack_direct and per_pv fits. The blocked stack_psis fit has no estimates to draw. The dashed line marks zero. The per_pv values come from made-up draws.", fig.alt = "Dot-and-interval plot with one row for b_female, at the bottom, and one for b_x, at the top. In each row the stack_direct estimate, in blue, is drawn slightly below the per_pv estimate, in yellow, each with its 95 percent interval. The stack_direct interval for b_female runs from about −42 to 46 and covers most of the horizontal axis; the per_pv interval for b_female, from about −2 to 6, and the two intervals for b_x, between about 43 and 51, are short. A dashed vertical line marks zero. There is no point or interval for stack_psis."----
slopes  <- est[est$term != "b_Intercept",
               c("method_label", "term", "estimate", "conf_low", "conf_high")]
slopes  <- slopes[is.finite(slopes$estimate), , drop = FALSE]
methods <- unique(slopes$method_label)
terms   <- sort(unique(slopes$term))

offsets <- setNames(seq(-0.18, 0.18, length.out = length(methods)), methods)
cols    <- setNames(c("#1f6f9c", "#f2c14e", "grey45")[seq_along(methods)], methods)

xlim <- range(c(slopes$conf_low, slopes$conf_high, 0), na.rm = TRUE)

op <- par(mar = c(4.5, 7, 1, 1))
plot(
  NA, xlim = xlim, ylim = c(0.5, length(terms) + 0.5),
  yaxt = "n", ylab = "",
  xlab = "Coefficient (synthetic reading-score points)"
)
abline(v = 0, lty = 2, col = "grey50")
for (mlab in methods) {
  sub <- slopes[slopes$method_label == mlab, ]
  yy  <- match(sub$term, terms) + offsets[[mlab]]
  segments(sub$conf_low, yy, sub$conf_high, yy, lwd = 2, col = cols[[mlab]])
  points(sub$estimate, yy, pch = 19, cex = 1.3, col = cols[[mlab]])
}
axis(2, at = seq_along(terms), labels = terms, las = 1)
legend("topleft", legend = methods, col = cols, pch = 19, lwd = 2,
       bty = "n", cex = 0.9)
par(op)

## ----diag-keys----------------------------------------------------------------
dg <- get_diagnostics(cmp)
names(dg)

## ----agreement----------------------------------------------------------------
dg$agreement

## ----target-overlap-----------------------------------------------------------
str(dg$target_overlap)

## ----method-diagnostics-------------------------------------------------------
dg$method_diagnostics[, c("method_label", "interval_role",
                          "coverage_claim_allowed", "target_source",
                          "psis_status", "pareto_k_max", "n_fits")]

