## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----example-data-shape-------------------------------------------------------
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

dim(pisa_tiny)        # 12 students, 12 columns
names(pisa_tiny)

## ----detect-design------------------------------------------------------------
detect_pisa_pv_columns(pisa_tiny, suffix = "READ")   # the M = 2 plausible values
detect_pisa_brr_replicate_weights(pisa_tiny)         # the R = 4 replicate weights

## ----read-target--------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)

c(M = tg$M, R = tg$R, fay_k = tg$fay_k)   # M, R and the Fay coefficient k
tg$fe_names                               # the p = 3 fixed-effect coefficients

