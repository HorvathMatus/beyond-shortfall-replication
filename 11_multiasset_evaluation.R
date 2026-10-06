# ============================================================
# Beyond Shortfall — 11: multi-asset forecast evaluation
#
# Extends the Bitcoin comparison (steps 03, 04, 08) to all seven risky assets:
#   1. Single-Index expectile forecasts for each asset (as in step 03)
#   2. FZ0 and AL losses, MCS p-values (full elimination path, as in step 04)
#      and backtests (Kupiec UC, Christoffersen CC, Bayer-Dimitriadis ESR, as in
#      step 08) for every asset, tail, quantile level and model
#   3. Cross-asset summary tables for the paper
#
# Inputs : step 10 forecasts in MATLAB/Beyond_Shortfall/assets/<ASSET>/
# Outputs: assets/<ASSET>/VaR_SI_*.txt, ES_SI_*.txt
#          data/multiasset_eval.rds          (long table of all results)
#          LaTeX/tables/MultiAsset_fixed.tex       (by model; appendix)
#          LaTeX/tables/MultiAssetByAsset_fixed.tex (by asset; appendix)
#
# Assets are processed in parallel (one process per asset).
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages({
  library(MCS); library(rugarch); library(esback); library(parallel); library(data.table)
})
source("R/single_index.R")
source("R/forecast_eval.R")       # losses, MCS p-values, read_fc, eval_asset, si_block
source("R/multiasset_tables.R")   # write_multiasset_tables

ROOT    <- file.path(MATLAB_DIR, "assets")
ASSETS  <- c("BNB", "BTC", "DOGE", "ETH", "SOL", "TRX", "XRP")
R_IN    <- 1095
NCORES  <- min(length(ASSETS), max(1, detectCores() - 2))

# ---- 1. Single-Index forecasts (as in step 03) -------------------------------
si_asset <- function(a) {
  set.seed(123)                                   # same RNG path as step 03
  dir <- file.path(ROOT, a)
  raw <- read.table(file.path(dir, "ret.txt"), col.names = c("date", "r"))
  raw <- raw[raw$r != 0, ]
  y_all <- raw$r; T <- length(y_all)
  oos <- scan(file.path(dir, "oos_ret.txt"), quiet = TRUE)
  stopifnot(length(oos) == T - R_IN, max(abs(y_all[(R_IN + 1):T] - oos)) < 1e-8)
  for (tl in c("low", "up")) {
    y <- if (tl == "low") y_all else -y_all
    for (k in seq_along(TAUS)) {
      f <- si_block(y, R_IN + 1, T, R_IN, TAUS[k])  # one block = estimated once
      write.table(f$v, file.path(dir, paste0("VaR_SI_", tl, "_", SUFFIX[k], ".txt")),
                  row.names = FALSE, col.names = FALSE)
      write.table(f$e, file.path(dir, paste0("ES_SI_", tl, "_", SUFFIX[k], ".txt")),
                  row.names = FALSE, col.names = FALSE)
    }
  }
  a
}

need_si <- ASSETS[!sapply(ASSETS, function(a)
  file.exists(file.path(ROOT, a, "ES_SI_up_01.txt")))]
if (length(need_si)) {
  cat("Single-Index model for:", paste(need_si, collapse = ", "), "\n")
  invisible(mclapply(need_si, si_asset, mc.cores = min(NCORES, length(need_si))))
}

# BTC check: the multi-asset SI forecasts must reproduce step 03
btc_old <- scan(file.path(MATLAB_DIR, "VaR_SI_low_005.txt"), quiet = TRUE)
btc_new <- scan(file.path(ROOT, "BTC", "VaR_SI_low_005.txt"), quiet = TRUE)
cat(sprintf("BTC SI reproduction: max abs difference %.2e\n", max(abs(btc_old - btc_new))))

# ---- 2. losses, MCS and backtests per asset (R/forecast_eval.R) --------------

cat("Evaluating", length(ASSETS), "assets on", NCORES, "cores ...\n")
t0 <- Sys.time()
ev <- rbindlist(mclapply(ASSETS, function(a) eval_asset(ROOT, a), mc.cores = NCORES))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
ev[, rank_fz := frank(fz, ties.method = "first", na.last = "keep"), by = .(asset, tail, tau)]
ev[, rank_fzh := frank(fzh, ties.method = "first", na.last = "keep"), by = .(asset, tail, tau)]
saveRDS(ev, "data/multiasset_eval.rds")

n_bad <- ev[is.na(fz), .N]
if (n_bad) { cat("cases with missing forecasts:\n"); print(ev[is.na(fz), .(asset, tail, tau, model)]) }

# ---- 3. summary tables (R/multiasset_tables.R); the main tables of Section 5
#         use the re-estimated models (step 15), these go to the appendix -------
write_multiasset_tables(ev, TEX_DIR, suffix = "_fixed", scheme_note = ", models estimated once")
