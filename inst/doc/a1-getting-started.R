## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----install-dev, eval = FALSE------------------------------------------------
# # install.packages("pak")
# pak::pak("joonho112/pvstackr")

## ----load-again, eval = FALSE-------------------------------------------------
# library(pvstackr)

## ----pv-fit-call, eval = FALSE------------------------------------------------
# fit <- pv_fit(
#   data    = your_data,
#   formula = OUTCOME ~ x + female,   # OUTCOME stands for each plausible value
#   target  = your_target,            # made by pv_brr_target()
#   method  = "stack_direct",
#   control = pv_control(method = "stack_direct", backend = "brms")
# )

## ----load-example-------------------------------------------------------------
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit

## ----example-data-------------------------------------------------------------
str(pisa_tiny)

## ----print-fit----------------------------------------------------------------
print(fit)

## ----get-estimates------------------------------------------------------------
est <- get_estimates(fit)

est[, c("term", "estimate", "se", "df",
        "conf_low", "conf_high",
        "interval_role", "coverage_claim_allowed")]

## ----coef-figure, fig.width = 7, fig.height = 3, fig.cap = "Slopes of the example fit, from synthetic data, with their 95% descriptive intervals. The intercept, about 458 points, is not shown. The dashed line marks zero.", fig.alt = "Dot-and-interval plot of the two slopes of the example fit, on a horizontal axis in synthetic reading-score points. b_x, at the top, is about 47 with a short interval from about 44 to 49, well to the right of zero. b_female, at the bottom, is about 2 with a long interval from about minus 42 to 46, which contains zero. A dashed vertical line marks zero."----
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

