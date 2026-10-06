# ============================================================
# Beyond Shortfall — 03: Single-Index expectile model (Jiang et al., 2020)
#
# The other ten models are produced by step 02 in MATLAB
# (see 02_run_matlab_forecasts.R). The Single-Index model has no MATLAB
# implementation — the authors' code is in R — so it is run here and its
# forecasts are written in the same format as the MATLAB output, ready for
# step 04.
#
# Design, matching the other models in step 02:
#   estimate ONCE on the in-sample period (first R = 1095 obs), then predict
#   forward over all P = 2099 out-of-sample days. The previous paper re-fitted
#   the model for every single day, which is not feasible over 2099 days and
#   is also inconsistent with how GARCH/GAS forecasts are produced there.
#
#   For every tau in {0.01, 0.025, 0.05, 0.1} and both tails:
#     1. calibrate the expectile level tau_e so the fitted in-sample
#        exceedance rate equals tau                     (SP500_find_alpha.R)
#     2. fit the single-index expectile model            (Jiang et al., 2020)
#     3. VaR = predicted expectile, ES from Taylor (2008)
#
#   Upper tail is handled exactly as in step 02: run everything on NEGATED
#   returns, so outputs live in the same negated space.
#
# Output (into the MATLAB folder, next to the other forecasts):
#   VaR_SI_<tail>_<suffix>.txt , ES_SI_<tail>_<suffix>.txt   (P x 1)
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
source("R/single_index.R")

R_IN   <- 1095                    # in-sample length, as in step 02
TAUS   <- c(0.01, 0.025, 0.05, 0.1)
SUFFIX <- c("001","0025","005","01")
TAILS  <- c(low = 1, up = -1)     # multiplier applied to the return series
set.seed(123)                     # MAVE starting values

# ---- data: the same filtered series step 02 uses --------------------------
raw <- read.table(file.path(MATLAB_DIR, "BTC.txt"), col.names = c("date","r"))
raw <- raw[raw$r != 0, ]
y_all <- raw$r
T <- length(y_all); P <- T - R_IN
stopifnot(T == 3194, P == 2099)

# alignment check against the MATLAB out-of-sample returns
oos_ref <- scan(file.path(MATLAB_DIR, "oos_ret.txt"), quiet = TRUE)
stopifnot(length(oos_ref) == P, max(abs(y_all[(R_IN+1):T] - oos_ref)) < 1e-8)
cat(sprintf("T=%d  in-sample R=%d  OOS P=%d  (aligned with step 02)\n\n", T, R_IN, P))

tr <- 3:R_IN                      # training rows (need 2 lags)
te <- (R_IN + 1):T                # out-of-sample rows

summary_tbl <- list()
for (tl in names(TAILS)) {
  s <- TAILS[[tl]]
  y <- s * y_all                  # upper tail = negated returns
  X_tr <- si_design(y[tr - 1], y[tr - 2]); y_tr <- y[tr]
  X_te <- si_design(y[te - 1], y[te - 2]); y_te <- y[te]
  mu_in <- mean(y[1:R_IN])        # in-sample mean, avoids look-ahead in ES

  for (k in seq_along(TAUS)) {
    tau <- TAUS[k]; suf <- SUFFIX[k]
    cat(sprintf("[%s] tau=%.3f: calibrating expectile level ...\n", tl, tau))
    t0 <- Sys.time()

    cal   <- si_calibrate_tau(y_tr, X_tr, alpha = tau, verbose = TRUE)
    fit   <- cal$fit
    tau_e <- cal$tau

    v <- si_predict(fit, X_te)                              # VaR forecasts
    e <- si_es_from_expectile(v, tau_e, tau, mu_in)         # ES forecasts

    # ES must sit beyond VaR; both must be negative for the FZ0/AL losses
    n_bad <- sum(e > v) + sum(v >= 0)
    if (n_bad > 0)
      cat(sprintf("   !! %d observations violate ES <= VaR < 0\n", n_bad))

    write.table(v, file.path(MATLAB_DIR, paste0("VaR_SI_", tl, "_", suf, ".txt")),
                row.names = FALSE, col.names = FALSE)
    write.table(e, file.path(MATLAB_DIR, paste0("ES_SI_",  tl, "_", suf, ".txt")),
                row.names = FALSE, col.names = FALSE)

    summary_tbl[[paste(tl, suf)]] <- data.frame(
      tail = tl, tau = tau, tau_e = round(tau_e, 5),
      hit_in = round(fit$hit_rate, 4), hit_oos = round(mean(y_te < v), 4),
      VaR_mean = round(mean(v), 3), ES_mean = round(mean(e), 3),
      secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
    cat(sprintf("   tau_e=%.5f  hit: in=%.4f oos=%.4f  VaR=%.3f ES=%.3f  (%.0fs)\n\n",
        tau_e, fit$hit_rate, mean(y_te < v), mean(v), mean(e),
        as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

res <- do.call(rbind, summary_tbl); rownames(res) <- NULL
cat("\n=== Single-Index summary ===\n"); print(res)
saveRDS(res, "data/single_index_summary.rds")
cat("\nForecasts written to", MATLAB_DIR, "as VaR_SI_<tail>_<suffix>.txt / ES_SI_...\n")
