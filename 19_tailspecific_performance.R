# ============================================================
# Beyond Shortfall — 19: performance of the tail-specific MQR portfolios
#
# Evaluates the portfolios of step 18 with the accounting of step 07 (weekly
# rebalancing, drifting weights, common out-of-sample span of 1080 days) and
# compares, for every value of c and both weight caps,
#   TailSpec   ES and EU from the models chosen per tail (step 17)
#   SameModel  one model for both tails, chosen by the same rule
#   GARCH-EDF  static ex-ante choice (best in both tails before July 2023)
# with the single-model portfolios of GARCH-Skt and GARCH-FZ (step 06) and
# the naive portfolio. The value of choosing the model tail by tail is the
# difference between TailSpec and SameModel, tested with the Ledoit-Wolf (2008)
# test of equal Sharpe ratios.
#
# Inputs : weights_*.txt (steps 06, 18), data/portfolio_input.rds,
#          data/tail_selection.rds
# Outputs: data/tailspecific_performance.rds, LaTeX/tables/TailSpecific.tex,
#          LaTeX/figures/tailspecific_sr.png
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages(library(data.table))
source("R/sharpe_test.R")

W_EST <- 1095; A_WIN <- 7; PPY <- 365; OOS_N <- 1080
CS    <- seq(0, 1, 0.1)
DESIGNS <- c("TailSpec" = "Tail-specific", "SameModel" = "Same model", "GARCH-EDF" = "GARCH-EDF (static)",
             "GARCH-Skt" = "GARCH-Skt", "GARCH-FZ" = "GARCH-FZ")
CAPS <- c(none = "", "0.3" = "_cap03")

inp <- readRDS("data/portfolio_input.rds")
R <- inp$ret / 100; colnames(R) <- inp$coins
N <- ncol(R); NPER <- floor((nrow(R) - W_EST) / A_WIN)
holds <- lapply(1:NPER, function(l) (W_EST + (l - 1) * A_WIN + 1):(W_EST + l * A_WIN))

run_strategy <- function(W) {                       # as in steps 07, 12, 13
  r_out <- numeric(0); to <- numeric(0); w_drift <- NULL
  for (l in 1:NPER) {
    w <- W[l, ]; Rl <- R[holds[[l]], , drop = FALSE]
    if (!is.null(w_drift)) to <- c(to, sum(abs(w - w_drift)))
    wc <- w
    for (d in 1:nrow(Rl)) { rp <- sum(wc * Rl[d, ]); r_out <- c(r_out, rp); wc <- wc * (1 + Rl[d, ]) / (1 + rp) }
    w_drift <- wc
  }
  list(r = r_out[1:OOS_N], TO = mean(to))
}
perf <- function(r, TO) {
  cw <- cumprod(1 + r); sd_a <- sd(r) * sqrt(PPY); mu_a <- mean(r) * PPY; md <- max(1 - cw / cummax(cw))
  data.table(TO = TO, CW = tail(cw, 1), SD = sd_a, MD = md, SR = mu_a / sd_a)
}

naive <- run_strategy(matrix(1 / N, NPER, N)); pn <- perf(naive$r, naive$TO)
res <- list(); rets <- list()
for (k in names(CAPS)) for (d in names(DESIGNS)) for (c in CS) {
  f <- file.path(MATLAB_DIR, sprintf("weights_%s%s_c%02d.txt", d, CAPS[[k]], round(c * 10)))
  if (!file.exists(f)) next
  W <- as.matrix(read.table(f, sep = ",")); s <- run_strategy(W)
  key <- paste(d, k, c); rets[[key]] <- s$r
  res[[key]] <- cbind(data.table(design = d, cap = k, c = c, stable = mean(rowSums(W[, c(7, 8)])),
                                 effN = mean(1 / rowSums(W^2))), perf(s$r, s$TO),
                      p_naive = sharpe_diff_test(s$r, naive$r)["p"])
}
res <- rbindlist(res)
# tail-specific vs same model, c by c
cmp <- rbindlist(lapply(names(CAPS), function(k) rbindlist(lapply(CS, function(c) {
  a <- rets[[paste("TailSpec", k, c)]]; b <- rets[[paste("SameModel", k, c)]]
  if (is.null(a) || is.null(b)) return(NULL)
  t <- sharpe_diff_test(a, b)
  data.table(cap = k, c = c, dSR = t["dSR"], p = t["p"])
}))))
if (!nrow(cmp)) cmp <- data.table(cap = character(0), c = numeric(0), dSR = numeric(0), p = numeric(0))
cat(sprintf("Naive: SR %.3f\n\n", pn$SR))
print(dcast(res, design + cap ~ c, value.var = "SR"), digits = 3)
cat("\nTail-specific minus same-model Sharpe ratio, by c:\n"); print(cmp, digits = 3)
saveRDS(list(results = res, compare = cmp, naive = pn), "data/tailspecific_performance.rds")

# ---- figure: Sharpe ratio against c ---------------------------------------------
cols <- c(TailSpec = "#d62728", SameModel = "#1f77b4", "GARCH-EDF" = "#2ca02c",
          "GARCH-Skt" = "grey40", "GARCH-FZ" = "grey70")
ltys <- c(TailSpec = 1, SameModel = 1, "GARCH-EDF" = 1, "GARCH-Skt" = 2, "GARCH-FZ" = 3)
png(file.path(FIG_DIR, "tailspecific_sr.png"), width = 1600, height = 720, res = 140)
layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(1, 0.13))
par(mar = c(4, 4, 2.2, 0.8))
for (k in c("0.3", "none")) {
  x <- res[cap == k]
  plot(NA, xlim = c(0, 1), ylim = range(c(x$SR, pn$SR), na.rm = TRUE), xlab = "c  (weight on the upper tail)",
       ylab = "Sharpe ratio", main = if (k == "none") "No weight cap" else "Weight cap 0.3", cex.main = 0.95)
  abline(h = pn$SR, lty = 3, col = "black")
  for (d in names(DESIGNS)) { y <- x[design == d][order(c)]
    if (nrow(y)) lines(y$c, y$SR, col = cols[d], lty = ltys[d], lwd = 2, type = "b", pch = 16, cex = 0.6) }
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = c(DESIGNS, "Naive"), col = c(cols, "black"), lty = c(ltys, 3), lwd = 2,
       bty = "n", cex = 0.8, horiz = TRUE)
dev.off()

# ---- LaTeX table: Sharpe ratio over the c grid -------------------------------------
CSHOW <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
pf <- function(p) if (is.na(p)) "--" else if (p < 0.001) "0.000" else sprintf("%.2f", p)
rows <- character(0)
for (k in c("0.3", "none")) {
  rows <- c(rows, sprintf("\\multicolumn{9}{l}{\\textit{%s}} \\\\",
                          if (k == "none") "Panel B: no weight cap" else "Panel A: weight cap $\\kappa = 0.3$"))
  for (d in names(DESIGNS)) {
    x <- res[design == d & cap == k]
    if (!nrow(x)) next
    sr <- sapply(CSHOW, function(cv) { v <- x[abs(c - cv) < 1e-9, SR]; if (length(v)) sprintf("%.2f", if (abs(v) < 0.005) 0 else v) else "--" })
    rows <- c(rows, sprintf("%-20s & %s & %.2f & %.1f \\\\", DESIGNS[[d]], paste(sr, collapse = " & "),
                            mean(x$SR), 100 * mean(x$stable)))
  }
  y <- cmp[cap == k]
  rows <- c(rows, sprintf("%-20s & %s & & \\\\", "$p$: tail-specific vs same model",
                          paste(sapply(CSHOW, function(cv) { v <- y[abs(c - cv) < 1e-9, p]; pf(if (length(v)) v else NA) }), collapse = " & ")))
  rows <- c(rows, "\\midrule")
}
rows <- rows[-length(rows)]
sel <- readRDS("data/tail_selection.rds")$selection
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{MQR portfolios with the forecasting model chosen per tail: Sharpe ratios.}",
  "\\label{tab:TailSpecific}", "\\begin{center}", "\\small", "\\begin{tabular}{lrrrrrrrr}", "\\toprule",
  "& \\multicolumn{6}{c}{$c$} & & \\\\", "\\cmidrule(lr){2-7}",
  "Design & 0 & 0.2 & 0.4 & 0.6 & 0.8 & 1.0 & Mean & Stable (\\%) \\\\", "\\midrule", rows,
  "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: Annualised Sharpe ratios of the MQR portfolios, weekly rebalancing, July 1, 2023 to June 14, ",
    "2026. Tail-specific: at each rebalancing date the ES is forecast by the model with the lowest average rank by ",
    "FZ$^{(0)}$ loss in the lower tail over the preceding 365 days, and the EU by the corresponding model for the ",
    "upper tail; ranks are averaged over the seven assets, at $\\tau = 0.05$, using the forecasts of the models ",
    "re-estimated every 91 days, and the candidates are RW-125, GARCH-N, GARCH-Skt, GARCH-EDF, GARCH-FZ and Hybrid. ",
    "The two tails receive different models at ", sum(sel$es != sel$eu), " of the ", nrow(sel), " weekly rebalancing dates ",
    "(", format(min(sel$date), "%B %d, %Y"), " to ", format(max(sel$date), "%B %d, %Y"), "; the holding period of the last date lies after ",
    "the evaluation period, which is covered by the first ", ceiling(OOS_N / A_WIN), "). Same model: the ",
    "candidate with the best average of the two ranks forecasts both tails. GARCH-EDF (static): the best model in ",
    "both tails over October 2020 to June 2023. Mean: average over $c \\in \\{0, 0.1, \\dots, 1\\}$. Stable: average ",
    "stablecoin weight over all $c$. $p$: p-value of the \\citet{Ledoit2008} test of equal Sharpe ratios of the ",
    "tail-specific and same-model portfolios. The naive portfolio has a Sharpe ratio of ", sprintf("%.2f", pn$SR), ". ",
    "Without a cap, the portfolios with $c \\leq 0.5$ hold almost only stablecoins, whose Sharpe ratios reflect very small ",
    "returns and volatilities and carry little economic meaning.}"),
  "\\end{table}")
writeLines(L, file.path(TEX_DIR, "TailSpecific.tex"))
cat("\nWritten:", file.path(TEX_DIR, "TailSpecific.tex"), "and", file.path(FIG_DIR, "tailspecific_sr.png"), "\n")
