# ============================================================
# Beyond Shortfall — 15: robustness — models re-estimated on a rolling window
#
# Steps 10-11 estimate every model once and filter it forward. Step 14
# (MATLAB) re-estimates the ten MATLAB models every K = 91 days on a rolling
# 1095-day window. This step
#   1. produces the Single-Index forecasts with the same rolling scheme,
#   2. evaluates the rolling forecasts exactly as step 11 (FZ0/AL losses, MCS
#      p-values, UC/CC/ESR backtests) for all seven assets,
#   3. compares the two schemes model by model: average rank, MCS membership,
#      backtest rejections, exceedance ratio, and Diebold-Mariano-type tests of
#      the FZ loss of the rolling against the fixed-parameter forecasts.
#
# Inputs : assets_roll/<A>/ (step 14), assets/<A>/ret.txt, data/multiasset_eval.rds
# Outputs: assets_roll/<A>/VaR_SI_*.txt, ES_SI_*.txt
#          data/rolling_eval.rds, LaTeX/tables/Rolling.tex
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages({
  library(MCS); library(rugarch); library(esback); library(parallel)
  library(data.table); library(sandwich)
})
source("R/single_index.R")
source("R/forecast_eval.R")
source("R/multiasset_tables.R")   # write_multiasset_tables

ROOT_FIX   <- file.path(MATLAB_DIR, "assets")
ROOT_ROLL  <- file.path(MATLAB_DIR, "assets_roll")
ASSETS     <- c("BNB", "BTC", "DOGE", "ETH", "SOL", "TRX", "XRP")
R_IN <- 1095
K    <- 91                                   # must match K_ of step 14
NCORES <- max(1, detectCores() - 2)

# ---- 1. Single-Index forecasts, re-estimated every K days --------------------
series <- lapply(setNames(ASSETS, ASSETS), function(a) {
  raw <- read.table(file.path(ROOT_FIX, a, "ret.txt"), col.names = c("date", "r"))
  y <- raw$r[raw$r != 0]
  oos <- scan(file.path(ROOT_ROLL, a, "oos_ret.txt"), quiet = TRUE)
  stopifnot(length(oos) == length(y) - R_IN, max(abs(y[(R_IN + 1):length(y)] - oos)) < 1e-8)
  y
})
need_si <- ASSETS[!sapply(ASSETS, function(a) file.exists(file.path(ROOT_ROLL, a, "ES_SI_up_01.txt")))]
if (length(need_si)) {
  jobs <- rbindlist(lapply(need_si, function(a)
    data.table(asset = a, s = seq(R_IN + 1, length(series[[a]]), by = K))))
  cat(sprintf("Single-Index model: %d blocks for %s on %d cores ...\n",
              nrow(jobs), paste(need_si, collapse = ", "), NCORES))
  t0 <- Sys.time()
  out <- mclapply(seq_len(nrow(jobs)), function(j) {
    set.seed(123)
    y_all <- series[[jobs$asset[j]]]; s <- jobs$s[j]; e <- min(s + K - 1, length(y_all))
    res <- list()
    for (tl in c("low", "up")) {
      y <- if (tl == "low") y_all else -y_all
      for (k in seq_along(TAUS)) res[[paste(tl, k)]] <- si_block(y, s, e, R_IN, TAUS[k])
    }
    res
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  stopifnot(!any(sapply(out, inherits, "try-error")))
  cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  for (a in need_si) {
    idx <- which(jobs$asset == a)                      # blocks in date order
    for (tl in c("low", "up")) for (k in seq_along(TAUS)) {
      key <- paste(tl, k)
      v <- unlist(lapply(out[idx], function(o) o[[key]]$v))
      e <- unlist(lapply(out[idx], function(o) o[[key]]$e))
      stopifnot(length(v) == length(series[[a]]) - R_IN)
      write.table(v, file.path(ROOT_ROLL, a, paste0("VaR_SI_", tl, "_", SUFFIX[k], ".txt")),
                  row.names = FALSE, col.names = FALSE)
      write.table(e, file.path(ROOT_ROLL, a, paste0("ES_SI_", tl, "_", SUFFIX[k], ".txt")),
                  row.names = FALSE, col.names = FALSE)
    }
  }
}

# ---- 2. evaluation of the rolling forecasts (as in step 11) -----------------
cat("Evaluating rolling forecasts ...\n")
ev <- rbindlist(mclapply(ASSETS, function(a) eval_asset(ROOT_ROLL, a),
                         mc.cores = min(length(ASSETS), NCORES)))
ev[, rank_fz := frank(fz, ties.method = "first", na.last = "keep"), by = .(asset, tail, tau)]
ev[, rank_fzh := frank(fzh, ties.method = "first", na.last = "keep"), by = .(asset, tail, tau)]
saveRDS(ev, "data/rolling_eval.rds")
# main cross-asset tables of Section 5
write_multiasset_tables(ev, TEX_DIR, suffix = "", scheme_note = "")
if (ev[is.na(fz), .N]) { cat("cases with missing forecasts:\n"); print(ev[is.na(fz), .(asset, tail, tau, model)]) }

# ---- 3. rolling vs fixed: loss-differential tests ----------------------------
# d_t = L_t(rolling) - L_t(fixed), FZ loss; t = mean(d) / sqrt(HAC variance of
# the mean), Andrews (1991) bandwidth with prewhitening, as in R/sharpe_test.R
dm <- rbindlist(lapply(ASSETS, function(a) {
  r_oos <- scan(file.path(ROOT_ROLL, a, "oos_ret.txt"), quiet = TRUE)
  rbindlist(lapply(c("low", "up"), function(tl) {
    y <- if (tl == "low") r_oos else -r_oos
    rbindlist(lapply(seq_along(TAUS), function(k) {
      f0 <- read_fc(file.path(ROOT_FIX, a), tl, SUFFIX[k]); f1 <- read_fc(file.path(ROOT_ROLL, a), tl, SUFFIX[k])
      rbindlist(lapply(setdiff(MODEL_ORDER, c("RW-125", "RW-250", "RW-500")), function(m) {
        d <- fz_loss(y, f1$v[, m], f1$e[, m], TAUS[k]) - fz_loss(y, f0$v[, m], f0$e[, m], TAUS[k])
        tstat <- if (all(is.finite(d))) mean(d) / sqrt(lrvar(d, type = "Andrews", prewhite = TRUE)) else NA_real_
        data.table(asset = a, tail = tl, tau = TAUS[k], model = m, dloss = mean(d), t = tstat)
      }))
    }))
  }))
}))

# ---- 4. comparison table -----------------------------------------------------
ev0 <- readRDS("data/multiasset_eval.rds")
SIG <- 0.05
summ <- function(x) x[, .(avgrank = mean(rank_fz, na.rm = TRUE),
                           mcs = sum(p_fz >= 0.10, na.rm = TRUE),
                           rej = sum(uc < SIG, na.rm = TRUE) + sum(cc < SIG, na.rm = TRUE) + sum(esr < SIG, na.rm = TRUE),
                           hitr = mean(hit / tau, na.rm = TRUE)), by = .(tail, model)]
cmp <- merge(summ(ev0), summ(ev), by = c("tail", "model"), suffixes = c("_f", "_r"))
cmp <- merge(cmp, dm[, .(better = sum(t < -1.96, na.rm = TRUE), worse = sum(t > 1.96, na.rm = TRUE),
                         lower = mean(dloss < 0, na.rm = TRUE)), by = .(tail, model)],
             by = c("tail", "model"), all.x = TRUE)
cmp[, model := factor(model, levels = MODEL_ORDER)]; setorder(cmp, tail, model)
cat("\n=== Fixed (f) vs rolling (r) estimation ===\n"); print(cmp, digits = 3)

best_model <- function(dt) dt[, .(ar = mean(rank_fz, na.rm = TRUE)), by = model][which.min(ar), model]
rho_tau <- function(x, a, tt) {
  lo <- x[asset == a & tail == "low" & abs(tau - tt) < 1e-9 & model != "Single Index", .(model, rank_fz)]
  up <- x[asset == a & tail == "up"  & abs(tau - tt) < 1e-9 & model != "Single Index", .(model, rank_fz)]
  z <- merge(lo, up, by = "model")[complete.cases(rank_fz.x, rank_fz.y)]
  cor(z$rank_fz.x, z$rank_fz.y, method = "spearman")
}
ba <- rbindlist(lapply(ASSETS, function(a) data.table(asset = a,
  best_es = best_model(ev[asset == a & tail == "low"]), best_eu = best_model(ev[asset == a & tail == "up"]),
  rho05 = rho_tau(ev, a, 0.05), rho_mean = mean(sapply(TAUS, function(t) rho_tau(ev, a, t))))))
cat("\n=== Rolling scheme, by asset ===\n"); print(ba, digits = 3)
cat(sprintf("mean rho: %.2f (tau = 0.05: %.2f)\n", mean(ba$rho_mean), mean(ba$rho05)))
saveRDS(list(eval = ev, dm = dm, compare = cmp, by_asset = ba), "data/rolling_eval.rds")

fmt_row <- function(x) sprintf("%-12s & %.2f & %.2f & %d & %d & %d & %d & %.2f & %.2f & %s \\\\",
  as.character(x$model), x$avgrank_r, x$avgrank_f, x$mcs_r, x$mcs_f, x$rej_r, x$rej_f,
  x$hitr_r, x$hitr_f, if (is.na(x$better)) "--" else sprintf("%d / %d", x$better, x$worse))
panel <- function(tl, title) {
  x <- cmp[tail == tl]
  c(sprintf("\\multicolumn{10}{l}{\\textit{%s}} \\\\", title),
    unlist(lapply(seq_len(nrow(x)), function(i) {
      r <- fmt_row(x[i]); if (x$model[i] %in% c("RW-500", "GARCH-FZ", "Hybrid")) c(r, "\\hline") else r })))
}
ncase <- length(ASSETS) * length(TAUS)
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Forecasts with models re-estimated every 91 days (Rolling) and estimated once (Fixed).}",
  "\\label{tab:Rolling}", "\\begin{center}", "\\begin{tabular}{lrrrrrrrrc}", "\\toprule",
  "& \\multicolumn{2}{c}{Avg rank} & \\multicolumn{2}{c}{MCS} & \\multicolumn{2}{c}{Rej} & \\multicolumn{2}{c}{Hit$/\\tau$} & DM \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "Model & Rolling & Fixed & Rolling & Fixed & Rolling & Fixed & Rolling & Fixed & $+$ / $-$ \\\\",
  "\\midrule", panel("low", "Panel A: lower tail (ES)"), "\\midrule",
  panel("up", "Panel B: upper tail (EU)"), "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: Results pooled over the seven assets and four quantile levels (", ncase,
    " cases per tail). Rolling: models re-estimated every ", K, " days on the most recent 1{,}095 ",
    "observations (Table \\ref{tab:MultiAsset}). Fixed: models estimated once on the first 1{,}095 observations and ",
    "filtered forward (Table \\ref{tab:MultiAsset_fixed}). Avg rank is the average rank by mean FZ loss among the eleven models of the same scheme; ",
    "MCS counts the cases in the 90\\% model confidence set; Rej counts the rejections of the UC, CC and ESR ",
    "backtests at the 5\\% level (maximum ", 3 * ncase, "); Hit$/\\tau$ is the average ratio of the exceedance rate ",
    "to the nominal level. DM reports the number of cases in which the FZ loss of the rolling forecasts is ",
    "significantly lower ($+$) or higher ($-$) than that of the fixed-parameter forecasts, by the two-sided ",
    "test of \\citet{Diebold1995} at the 5\\% level with the HAC variance of \\citet{Andrews1992}. The rolling-window estimators do not ",
    "depend on the estimation scheme.}"),
  "\\end{table}")
writeLines(L, file.path(TEX_DIR, "Rolling.tex"))
cat("\nTable written:", file.path(TEX_DIR, "Rolling.tex"), "\n")
