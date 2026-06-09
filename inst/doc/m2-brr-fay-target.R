## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center", out.width = "85%")

## ----library------------------------------------------------------------------
library(pvstackr)

## ----common-load--------------------------------------------------------------
fit <- readRDS(
  system.file("extdata", "examples", "pisa_tiny_stack_direct.rds",
              package = "pvstackr")
)$fit
tg <- get_target(fit)   # external Rubin / BRR–Fay target (class pvstackr_brr_target)

class(tg)
c(M = tg$M, R = tg$R, fay_k = tg$fay_k)   # 2, 4, 0.5

## ----brrfay-peek--------------------------------------------------------------
tg$fay_variance_multiplier                 # a_d = 1 on this fixture (coincidence)
dim(tg$per_pv[[1]]$replicate_diff)         # 3 x 4 : (beta_hat_m^(r) - beta_hat_m)
str(tg$per_pv[[1]]$replicate_diff)

## ----rubin-pooling------------------------------------------------------------
M     <- tg$M                                    # 2
betas <- sapply(tg$per_pv, function(p) p$beta)   # 3 x M : per-PV beta_hat_m
Us    <- lapply(tg$per_pv, function(p) p$U)      # list of per-PV U_hat_m (BRR–Fay)

beta_bar <- rowMeans(betas)                      # Rubin mean
U_bar    <- Reduce(`+`, Us) / M                  # within-imputation covariance
dev      <- betas - beta_bar
B        <- (dev %*% t(dev)) / (M - 1)           # between-imputation covariance
T_MI     <- U_bar + (1 + 1/M) * B                # EQ-TMI

all.equal(unname(beta_bar), unname(tg$beta_bar)) # TRUE
all.equal(unname(U_bar),    unname(tg$U_bar))    # TRUE
all.equal(unname(B),        unname(tg$B))        # TRUE
all.equal(unname(T_MI),     unname(tg$T_MI))     # TRUE

## ----fmi-recompute------------------------------------------------------------
fmi <- (1 + 1/M) * diag(B) / diag(T_MI)          # EQ-FMI
all.equal(unname(fmi), unname(tg$fmi))           # TRUE

## ----df-table-----------------------------------------------------------------
data.frame(
  term      = tg$fe_names,
  fmi       = round(tg$fmi, 3),
  riv       = round(tg$riv, 2),
  df        = round(tg$df,  2),
  df_method = tg$df_method
)

## ----df-complete--------------------------------------------------------------
tg$df_complete   # NA for every term -> the classic path, not Barnard–Rubin

## ----session-info-------------------------------------------------------------
sessionInfo()

