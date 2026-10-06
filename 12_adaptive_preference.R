# ============================================================
# Beyond Shortfall — 12: ex-ante choice of the preference coefficient c
#
# Steps 06-07 report each strategy for a fixed c and identify the best c ex
# post. Here c is chosen ex ante: at every rebalancing date, the investor
# selects the value of c whose portfolio would have achieved the highest
# Sharpe ratio over the preceding LOOKBACK days, using only returns realised
# before that date, and holds the corresponding weights until the next date.
# During the initial burn-in, when no such history exists, the conservative
# minimum-ES choice c = 0 is used. No information from the evaluation period
# beyond the decision date enters the choice.
#
# The weights are those of step 06 (no new optimisation is needed): the
# adaptive strategy at date l holds the step-06 weights of the selected c at l.
#
# Inputs : weights_<model>[_cap<kk>]_c<cc>.txt (weekly), data/portfolio_input.rds
# Outputs: data/adaptive_c.rds, LaTeX/tables/AdaptiveC.tex
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
source("R/sharpe_test.R")

W_EST <- 1095; A_WIN <- 7; PPY <- 365
CS    <- seq(0, 1, 0.1)
MODS  <- c("RW125", "GARCH-Skt", "Hybrid", "GARCH-FZ")
VARIANTS <- c(none = "", "0.5" = "_cap05", "0.3" = "_cap03")
LOOKBACKS <- c(91, 182, 365)           # trailing windows (days) for the choice
OOS_N <- 1080                          # common out-of-sample span of step 07

inp <- readRDS("data/portfolio_input.rds")
R <- inp$ret / 100; colnames(R) <- inp$coins
N <- ncol(R); NPER <- floor((nrow(R) - W_EST) / A_WIN)
holds <- lapply(1:NPER, function(l) (W_EST + (l - 1) * A_WIN + 1):(W_EST + l * A_WIN))

# daily returns, turnover, and the drifted weights, for a target-weight matrix
# (same accounting as step 07)
run_strategy <- function(W) {
  r_out <- numeric(0); to <- numeric(0); w_drift <- NULL
  for (l in 1:NPER) {
    w <- W[l, ]; Rl <- R[holds[[l]], , drop = FALSE]
    if (!is.null(w_drift)) to <- c(to, sum(abs(w - w_drift)))
    wc <- w
    for (d in 1:nrow(Rl)) {
      rp <- sum(wc * Rl[d, ]); r_out <- c(r_out, rp)
      wc <- wc * (1 + Rl[d, ]) / (1 + rp)
    }
    w_drift <- wc
  }
  list(r = r_out, TO = mean(to))
}
perf <- function(r, TO) {
  r <- r[1:OOS_N]; cw <- cumprod(1 + r); sd_a <- sd(r) * sqrt(PPY); mu_a <- mean(r) * PPY
  data.frame(TO = TO, CW = tail(cw, 1), SD = sd_a, MD = max(1 - cw / cummax(cw)),
             SR = mu_a / sd_a, CR = mu_a / max(1 - cw / cummax(cw)))
}
sr_window <- function(x) { s <- sd(x); if (!is.finite(s) || s == 0) -Inf else mean(x) / s }

naive <- run_strategy(matrix(1 / N, NPER, N))

res <- list()
for (m in MODS) for (v in names(VARIANTS)) {
  tag <- paste0(m, VARIANTS[[v]])
  files <- file.path(MATLAB_DIR, sprintf("weights_%s_c%02d.txt", tag, round(CS * 10)))
  if (!all(file.exists(files))) next
  Ws <- lapply(files, function(f) as.matrix(read.table(f, sep = ",")))
  rets <- lapply(Ws, function(W) run_strategy(W)$r)     # daily returns of each fixed-c strategy

  for (lb in LOOKBACKS) {
    chosen <- integer(NPER)
    for (l in 1:NPER) {
      past_end <- (l - 1) * A_WIN                          # days realised before date l
      if (past_end < lb) { chosen[l] <- 1; next }          # burn-in: c = 0
      idx <- (past_end - lb + 1):past_end
      chosen[l] <- which.max(sapply(rets, function(r) sr_window(r[idx])))
    }
    W_ad <- t(sapply(1:NPER, function(l) Ws[[chosen[l]]][l, ]))
    s <- run_strategy(W_ad)
    p <- perf(s$r, s$TO)
    t_sr <- sharpe_diff_test(s$r[1:OOS_N], naive$r[1:OOS_N])
    res[[length(res) + 1]] <- cbind(data.frame(model = m, cap = v, lookback = lb,
      mean_c = mean(CS[chosen]), share_c0 = mean(chosen == 1)), p,
      pSR = t_sr["p"], dSR = t_sr["dSR"])
  }
}
res <- do.call(rbind, res); rownames(res) <- NULL
pn <- perf(naive$r, naive$TO)
cat(sprintf("Naive: CW %.3f SD %.3f SR %.3f\n\n", pn$CW, pn$SD, pn$SR))
print(res, digits = 3)
saveRDS(list(adaptive = res, naive = pn), "data/adaptive_c.rds")

# ---- LaTeX: main lookback (182 days) for all models and caps ---------------
LB_MAIN <- 182
x <- res[res$lookback == LB_MAIN, ]
pf <- function(p) if (is.na(p)) "--" else if (p < 0.001) "0.000" else sprintf("%.3f", p)
rows <- character(0)
for (m in MODS) {
  xm <- x[x$model == m, ]
  for (i in seq_len(nrow(xm))) {
    rows <- c(rows, sprintf("%-10s & %s & %.2f & %.3f & %.3f & %.3f & %.3f & %.3f & %s \\\\",
      if (i == 1) sub("RW125", "RW-125", m) else "", xm$cap[i], xm$mean_c[i],
      xm$TO[i], xm$CW[i], xm$SD[i], xm$MD[i], xm$SR[i], pf(xm$pSR[i])))
  }
  rows <- c(rows, "\\hline")
}
rows <- c(rows, sprintf("%-10s & -- & -- & %.3f & %.3f & %.3f & %.3f & %.3f & -- \\\\",
  "Naive", pn$TO, pn$CW, pn$SD, pn$MD, pn$SR))
rb <- sapply(LOOKBACKS, function(lb) { y <- res[res$lookback == lb, ]; sprintf("%.2f--%.2f", min(y$SR), max(y$SR)) })
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{MQR portfolios with the preference coefficient chosen ex ante.}",
  "\\label{tab:AdaptiveC}", "\\begin{center}", "\\begin{tabular}{llrrrrrrr}", "\\toprule",
  "Model & $\\kappa$ & $\\bar{c}$ & TO & CW & SD & MD & SR & $p_{SR}$ \\\\", "\\midrule", rows,
  "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: At each weekly rebalancing date the value of $c \\in \\{0, 0.1, \\dots, 1\\}$ whose ",
    "portfolio attained the highest Sharpe ratio over the preceding ", LB_MAIN, " days is selected, ",
    "using only information available at that date; $c = 0$ is used until ", LB_MAIN, " days of history ",
    "are available. $\\bar{c}$ is the average selected value. $p_{SR}$ is the p-value of the ",
    "\\citet{Ledoit2008} test of equal Sharpe ratios against the naive portfolio. With look-back windows ",
    "of ", paste(LOOKBACKS, collapse = ", "), " days the Sharpe ratios range over ",
    paste(rb, collapse = ", "), " respectively.}"),
  "\\end{table}")
writeLines(L, file.path(TEX_DIR, "AdaptiveC.tex"))
cat("\nTable written:", file.path(TEX_DIR, "AdaptiveC.tex"), "\n")
