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

# Synthetic per-PV draw clouds, centred near the fixture's coefficients
# (~458 / 47 / 2). These are random noise for illustration, NOT real posteriors.
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

# Synthetic stacked draw cloud + equal (placeholder) PSIS weights and small
# Pareto-k. Again injected for illustration; nothing is importance-sampled live.
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

## ----overlay-figure, fig.width = 7, fig.height = 3.4, fig.cap = "The two slope coefficients (b_x, b_female) under all three methods, with their 95% intervals overlaid. The dashed line marks zero. Because the per_pv and stack_psis rows are built from synthetic injected draws, their points and widths are illustrative only — the figure shows the comparison structure and the per-method interval semantics, not a substantive agreement result.", fig.alt = "A horizontal dot-and-interval plot with two rows, one for the coefficient on x near 47 and one for the coefficient on female near 2. Within each row the three methods stack_direct, per_pv, and stack_psis are drawn at slightly offset heights in different colours, each with its own 95 percent interval. A dashed vertical reference line is drawn at zero."----
slopes  <- est[est$term != "b_Intercept",
               c("method_label", "term", "estimate", "conf_low", "conf_high")]
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
legend("topright", legend = methods, col = cols, pch = 19, lwd = 2,
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

## ----session-info-------------------------------------------------------------
sessionInfo()

