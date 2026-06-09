## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----fixture-shape------------------------------------------------------------
pisa_tiny <- read.csv(
  system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
)

dim(pisa_tiny)        # 12 rows x 12 columns
names(pisa_tiny)

## ----detect-design------------------------------------------------------------
detect_pisa_pv_columns(pisa_tiny, suffix = "READ")   # -> c("PV1READ","PV2READ")  => M = 2
detect_pisa_brr_replicate_weights(pisa_tiny)         # -> 4 W_FSTURWT* cols        => R = 4

## ----read-target--------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)

c(M = tg$M, R = tg$R, fay_k = tg$fay_k)   # 2, 4, 0.5
tg$fe_names                               # "b_Intercept" "b_x" "b_female"

## ----session-info-------------------------------------------------------------
sessionInfo()

