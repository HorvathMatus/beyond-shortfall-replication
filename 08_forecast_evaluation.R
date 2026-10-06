# ============================================================
# Beyond Shortfall — 08: absolute calibration of the tail forecasts and
#                       the asymmetry of model rankings across tails
#
# The MCS of step 04 ranks the models RELATIVE to one another. A model can be
# the best of the set and still be miscalibrated, so this step adds absolute
# backtests for every model, tail and quantile level:
#   UC   Kupiec (1995) unconditional coverage of the VaR
#   CC   Christoffersen (1998) conditional coverage (coverage + independence)
#   ESR  Bayer & Dimitriadis (2022) strict expected shortfall regression test
# It also summarises how the model ranking changes between the two tails.
#
# Inputs : step 02/03 forecasts in the MATLAB folder, step 04 MCS results
# Outputs: LaTeX/tables/Backtest.tex       (tau = 0.05, p-values)
#          LaTeX/tables/BacktestCount.tex  (rejections across all tau)
#          LaTeX/tables/RankShift.tex      (average rank, ES vs EU)
# ============================================================
# SCHEME: "fixed" (default; models estimated once, step 02; tables with the
# suffix _fixed, reported in the appendix) or "rolling" (models re-estimated
# every 91 days, steps 14-15; the main tables of Section 5, produced by step 21).
SCHEME <- if (exists("SCHEME")) SCHEME else "fixed"
rm(list = setdiff(ls(), "SCHEME"))
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages({ library(rugarch); library(esback) })
source("R/forecast_eval.R")   # read_fc with the real-time guard, MODEL_ORDER, TAUS

FC_DIR  <- if (SCHEME == "rolling") file.path(MATLAB_DIR, "assets_roll", "BTC") else MATLAB_DIR
SUF_OUT <- if (SCHEME == "rolling") "" else "_fixed"
CAP_NOTE <- if (SCHEME == "rolling") "" else ", models estimated once"
SCHEME_NOTE <- if (SCHEME == "rolling") "Models re-estimated every 91 days on a rolling window of 1{,}095 days. " else "Models estimated once on the first 1{,}095 days and filtered forward. "
TAILS  <- c(low = "Shortfall", up = "Upside")
SIG    <- 0.05

r_oos <- scan(file.path(FC_DIR, "oos_ret.txt"), quiet = TRUE)

# ---- run all backtests -------------------------------------------------------
bt <- list()
for (tl in names(TAILS)) {
  y <- if (tl == "low") r_oos else -r_oos       # upper tail in negated space
  for (k in seq_along(TAUS)) {
    tau <- TAUS[k]; fc <- read_fc(FC_DIR, tl, SUFFIX[k])
    for (m in MODEL_ORDER) {
      v <- fc$v[, m]; e <- fc$e[, m]
      if (anyNA(v)) {                       # no usable forecasts for this model
        bt[[length(bt) + 1]] <- data.frame(tail = tl, tau = tau, model = m,
          hit = NA, uc = NA, cc = NA, esr = NA); next }
      vt  <- VaRTest(alpha = tau, actual = y, VaR = v, conf.level = 1 - SIG)
      esr <- tryCatch(esr_backtest(y, v, e, alpha = tau, version = 1)$pvalue_twosided_asymptotic,
                      error = function(err) NA_real_)
      bt[[length(bt) + 1]] <- data.frame(tail = tl, tau = tau, model = m,
        hit = mean(y <= v), uc = vt$uc.LRp, cc = vt$cc.LRp, esr = esr)
    }
  }
}
bt <- do.call(rbind, bt)
saveRDS(bt, paste0("data/backtests", SUF_OUT, ".rds"))

cat("\n=== Backtest p-values at tau = 0.05 ===\n")
print(transform(bt[bt$tau == 0.05, ], hit = round(hit, 4), uc = round(uc, 3),
                cc = round(cc, 3), esr = round(esr, 3)), row.names = FALSE)

# ---- LaTeX helpers ------------------------------------------------------------
pfmt <- function(p) {
  if (is.na(p)) return("--")
  s <- if (p < 0.001) "0.000" else sprintf("%.3f", p)
  if (p < SIG) paste0("\\textbf{", s, "}") else s
}
wrap <- function(rows, label, caption, notes, hdr, align) {
  c("\\begin{table}[ht]", "\\centering", paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"), "\\begin{center}",
    paste0("\\begin{tabular}{", align, "}"), "\\toprule", hdr, "\\midrule", rows,
    "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
    paste0("\\textit{Notes: ", notes, "}"), "\\end{table}")
}
sep_after <- c("RW-500", "GARCH-FZ", "Hybrid")
dir.create(TEX_DIR, showWarnings = FALSE)

# ---- Table: detailed p-values at tau = 0.05 ----------------------------------
rows <- character(0)
for (m in MODEL_ORDER) {
  a <- bt[bt$tail == "low" & bt$tau == 0.05 & bt$model == m, ]
  b <- bt[bt$tail == "up"  & bt$tau == 0.05 & bt$model == m, ]
  hf <- function(h) if (is.na(h)) "--" else sprintf("%.3f", h)
  rows <- c(rows, sprintf("%-12s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", m,
    hf(a$hit), pfmt(a$uc), pfmt(a$cc), pfmt(a$esr), hf(b$hit), pfmt(b$uc), pfmt(b$cc), pfmt(b$esr)))
  if (m %in% sep_after) rows <- c(rows, "\\hline")
}
writeLines(wrap(rows, paste0("tab:Backtest", SUF_OUT),
  paste0("Backtests of the VaR and tail-expectation forecasts at $\\tau = 0.05$", CAP_NOTE, "."),
  paste0(SCHEME_NOTE, "Hit is the empirical exceedance rate (nominal 0.05). UC and CC are p-values of the ",
    "unconditional coverage test of \\citet{Kupiec1995} and the conditional coverage test of ",
    "\\citet{Christoffersen1998}; ESR is the p-value of the strict expected shortfall regression ",
    "test of \\citet{Bayer2020}. The upper-tail tests are applied to the negated return series. ",
    "Bold entries indicate rejection of correct calibration at the 5\\% level."),
  paste0("& \\multicolumn{4}{c}{Lower tail (ES)} & \\multicolumn{4}{c}{Upper tail (EU)} \\\\\n",
         "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\n",
         "Model & Hit & UC & CC & ESR & Hit & UC & CC & ESR \\\\"),
  "lrrrrrrrr"), file.path(TEX_DIR, paste0("Backtest", SUF_OUT, ".tex")))

# ---- Table: number of rejections across the four quantile levels ------------
rows <- character(0)
for (m in MODEL_ORDER) {
  cnt <- function(tl, col) sum(bt[bt$tail == tl & bt$model == m, col] < SIG, na.rm = TRUE)
  rows <- c(rows, sprintf("%-12s & %d & %d & %d & %d & %d & %d & %d & %d \\\\", m,
    cnt("low","uc"), cnt("low","cc"), cnt("low","esr"),
    cnt("low","uc") + cnt("low","cc") + cnt("low","esr"),
    cnt("up","uc"),  cnt("up","cc"),  cnt("up","esr"),
    cnt("up","uc") + cnt("up","cc") + cnt("up","esr")))
  if (m %in% sep_after) rows <- c(rows, "\\hline")
}
writeLines(wrap(rows, paste0("tab:BacktestCount", SUF_OUT),
  paste0("Number of backtest rejections across the four quantile levels", CAP_NOTE, "."),
  paste0(SCHEME_NOTE, "Each entry counts the quantile levels $\\tau \\in \\{0.01, 0.025, 0.05, 0.1\\}$ at which ",
    "the corresponding test rejects correct calibration at the 5\\% level (maximum 4 per test, ",
    "12 per tail). UC, CC and ESR as in Table \\ref{tab:Backtest", SUF_OUT, "}."),
  paste0("& \\multicolumn{4}{c}{Lower tail (ES)} & \\multicolumn{4}{c}{Upper tail (EU)} \\\\\n",
         "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\n",
         "Model & UC & CC & ESR & Total & UC & CC & ESR & Total \\\\"),
  "lrrrrrrrr"), file.path(TEX_DIR, paste0("BacktestCount", SUF_OUT, ".tex")))

# ---- Table: average rank in the two tails ------------------------------------
load(paste0("data/mcs_results_beyond", SUF_OUT, ".RData"))   # res, tables (step 04)
ar <- function(key) tables[[key]][MODEL_ORDER, "AvgRank"]
fzE <- ar("Shortfall_fz"); fzU <- ar("Upside_fz")
alE <- ar("Shortfall_fzh"); alU <- ar("Upside_fzh")
rho <- sapply(seq_along(TAUS), function(k)
  cor(tables$Shortfall_fz[MODEL_ORDER, paste0("Rank_", TAUS[k])],
      tables$Upside_fz[MODEL_ORDER, paste0("Rank_", TAUS[k])], method = "spearman"))
# excluding the Single Index, which is last in every column and inflates rho
rho10 <- sapply(seq_along(TAUS), function(k) {
  mm <- setdiff(MODEL_ORDER, "Single Index")
  cor(tables$Shortfall_fz[mm, paste0("Rank_", TAUS[k])],
      tables$Upside_fz[mm, paste0("Rank_", TAUS[k])], method = "spearman") })
rows <- character(0)
for (i in seq_along(MODEL_ORDER)) {
  rows <- c(rows, sprintf("%-12s & %.2f & %.2f & %+.2f & %.2f & %.2f & %+.2f \\\\", MODEL_ORDER[i],
    fzE[i], fzU[i], fzU[i] - fzE[i], alE[i], alU[i], alU[i] - alE[i]))
  if (MODEL_ORDER[i] %in% sep_after) rows <- c(rows, "\\hline")
}
writeLines(wrap(rows, paste0("tab:RankShift", SUF_OUT),
  paste0("Average rank of each model in the lower and the upper tail", CAP_NOTE, "."),
  paste0("Average rank over $\\tau \\in \\{0.01, 0.025, 0.05, 0.1\\}$ from Tables ",
    "\\ref{tab:FZLoss_ES", SUF_OUT, "}--\\ref{tab:FZHLoss_EU", SUF_OUT, "}; $\\Delta$ = EU rank minus ES rank, so a positive ",
    "value indicates a model that forecasts the upper tail worse than the lower tail. ",
    sprintf("Spearman rank correlations between the ES and EU rankings (FZ loss) at $\\tau = 0.01, 0.025, 0.05, 0.1$ are %s for all eleven models and %s when the Single Index model is excluded.",
      paste(sprintf("%.2f", rho), collapse = ", "), paste(sprintf("%.2f", rho10), collapse = ", "))),
  paste0("& \\multicolumn{3}{c}{FZ$^{(0)}$ loss} & \\multicolumn{3}{c}{FZ$^{(1/2)}$ score} \\\\\n",
         "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n",
         "Model & ES & EU & $\\Delta$ & ES & EU & $\\Delta$ \\\\"),
  "lrrrrrr"), file.path(TEX_DIR, paste0("RankShift", SUF_OUT, ".tex")))

cat("\nSpearman rho ES vs EU ranks (FZ), all models:  ", round(rho, 2), "\n")
cat("Spearman rho ES vs EU ranks (FZ), excl. SI:     ", round(rho10, 2), "\n")
cat(sprintf("\nTables written: Backtest%s.tex, BacktestCount%s.tex, RankShift%s.tex\n", SUF_OUT, SUF_OUT, SUF_OUT))
