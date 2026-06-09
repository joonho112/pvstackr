## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----install-cran, eval = FALSE-----------------------------------------------
# # install.packages("pvstackr")

## ----install-dev, eval = FALSE------------------------------------------------
# # install.packages("pak")
# pak::pak("joonho112/pvstackr")

## ----load-again, eval = FALSE-------------------------------------------------
# library(pvstackr)

## ----the-one-call, eval = FALSE-----------------------------------------------
# fit <- pv_fit(
#   data          = your_data,
#   formula       = OUTCOME ~ x + female,           # OUTCOME is a placeholder
#   target        = your_target,                    # from pv_brr_target()
#   method        = "stack_direct",
#   control       = pv_control(method = "stack_direct"),
#   fit_function  = your_fit_function,              # the modelling backend
#   draws_function = your_draws_function            # optional posterior draws
# )

## ----load-fixture-------------------------------------------------------------
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit

## ----fixture-peek-------------------------------------------------------------
str(pisa_tiny)

## ----print-fit----------------------------------------------------------------
print(fit)

## ----get-estimates------------------------------------------------------------
est <- get_estimates(fit)

est[, c("term", "estimate", "se", "df",
        "conf_low", "conf_high",
        "interval_role", "coverage_claim_allowed")]

## ----coef-figure, fig.width = 7, fig.height = 3, fig.cap = "Slope coefficients from the cached synthetic stack_direct fit, with 95% descriptive intervals. The intercept is omitted because its scale (~458 score points) would dominate the axis. The dashed line marks zero (no effect). These are illustrative synthetic values, not real PISA estimates.", fig.alt = "A horizontal dot-and-interval plot of two slope coefficients from the synthetic fixture fit. The coefficient on x is about 47 score points with a narrow interval well to the right of zero. The coefficient on female is about 2 score points with a wide interval spanning zero from roughly minus 42 to plus 46. A dashed vertical reference line is drawn at zero."----
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
par(op)

## ----get-draws----------------------------------------------------------------
get_draws(fit)

## ----honest-columns-----------------------------------------------------------
est[, c("term", "interval_role", "coverage_claim_allowed")]

## ----session-info-------------------------------------------------------------
sessionInfo()

