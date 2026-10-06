# ============================================================
# Beyond Shortfall — 07: out-of-sample portfolio performance
#
# Reads the weights produced by step 06 (MATLAB) and evaluates them
# out-of-sample, together with the benchmarks:
#   Naive    1/N, rebalanced on the same dates
#   MinVar   minimum-variance on the same rolling window
#   Bitcoin  buy and hold
#
# Measures follow the previous paper (DeMiguel et al. 2009; Petukhina et al.
# 2021): Turnover, Target Turnover, Cumulative Wealth, annualised SD,
# Maximum Drawdown, Sharpe Ratio, Calmar Ratio.
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
library(data.table)
source("R/sharpe_test.R")

W_EST <- 1095; PPY <- 365                  # crypto trades every day
FREQS <- c(7, 14, 30, 91)                  # rebalancing frequencies (days)

inp    <- readRDS("data/portfolio_input.rds")
R      <- inp$ret                          # T x N, percent
dates  <- inp$dates
coins  <- inp$coins
T_all  <- nrow(R); N <- ncol(R)

nper_of  <- function(f) floor((T_all - W_EST) / f)
# holding block for rebalancing l at frequency f
holds_of <- function(f) lapply(1:nper_of(f), function(l)
                          (W_EST + (l-1)*f + 1):(W_EST + l*f))

# Frequencies imply slightly different out-of-sample spans (e.g. 30-day gives
# 1080 days against 1092 for the others). All strategies are evaluated over
# the common span so that the comparison across frequencies is like-for-like.
OOS_N <- min(sapply(FREQS, function(f) nper_of(f) * f))
A_WIN <- 7                                 # baseline used for benchmarks
NPER  <- nper_of(A_WIN)
hold_idx <- holds_of(A_WIN)
oos_idx  <- unlist(hold_idx)[1:OOS_N]
cat(sprintf("OOS: %s .. %s  (%d days, common span across frequencies)\n",
            as.character(dates[min(oos_idx)]), as.character(dates[max(oos_idx)]), OOS_N))
for (f in FREQS) cat(sprintf("   %2d-day: %3d rebalancings\n", f, nper_of(f)))

# ---- portfolio returns and turnover for a weight matrix -------------------
# Wmat: NPER x N target weights. Returns the daily OOS return series,
# turnover (TO, includes drift) and target turnover (TTO).
run_strategy <- function(Wmat, freq = A_WIN) {
  holds <- holds_of(freq); np <- nper_of(freq)
  r_out <- numeric(0); to <- numeric(0); tto <- numeric(0)
  to_day <- integer(0)                          # day index of each rebalancing
  w_prev_target <- NULL; w_prev_drift <- NULL
  for (l in 1:np) {
    w  <- Wmat[l, ]
    idx <- holds[[l]]
    Rl <- R[idx, , drop = FALSE] / 100          # decimal returns
    if (!is.null(w_prev_target)) {
      tto <- c(tto, sum(abs(w - w_prev_target)))
      to  <- c(to,  sum(abs(w - w_prev_drift)))
      to_day <- c(to_day, length(r_out) + 1L)   # cost charged on the first day
    }
    # daily portfolio returns with weights drifting inside the block
    wc <- w
    for (d in 1:nrow(Rl)) {
      rp <- sum(wc * Rl[d, ])
      r_out <- c(r_out, rp)
      wc <- wc * (1 + Rl[d, ]) / (1 + rp)       # drift
    }
    w_prev_target <- w; w_prev_drift <- wc
  }
  keep <- to_day <= OOS_N                       # truncate to the common span
  list(r = r_out[1:OOS_N], TO = mean(to), TTO = mean(tto),
       to_series = to[keep], to_day = to_day[keep])
}

# net-of-cost return series: charge `bp` basis points per unit of turnover
net_returns <- function(s, bp) {
  r <- s$r
  r[s$to_day] <- r[s$to_day] - (bp/1e4) * s$to_series
  r
}

# ---- performance measures --------------------------------------------------
perf <- function(r, TO = NA, TTO = NA) {
  cw   <- cumprod(1 + r)
  sd_a <- sd(r) * sqrt(PPY)
  mu_a <- mean(r) * PPY
  dd   <- 1 - cw / cummax(cw)
  md   <- max(dd)
  data.frame(TO = TO, TTO = TTO, CW = tail(cw, 1), SD = sd_a, MD = md,
             SR = mu_a / sd_a, CR = mu_a / md)
}

# ---- benchmarks ------------------------------------------------------------
W_naive <- matrix(1/N, nrow = NPER, ncol = N)

W_mv <- t(sapply(1:NPER, function(l) {
  win <- R[((l-1)*A_WIN + 1):(W_EST + (l-1)*A_WIN), , drop = FALSE]
  S   <- cov(win); S <- S + diag(1e-8, N)
  # long-only minimum variance by projected gradient on the simplex
  w <- rep(1/N, N)
  for (it in 1:2000) {
    g  <- as.vector(S %*% w)
    w2 <- w - g / max(abs(g)) * 0.01
    w2 <- pmax(w2, 0); w2 <- w2 / sum(w2)
    if (max(abs(w2 - w)) < 1e-10) { w <- w2; break }
    w <- w2
  }
  w
}))

res <- list(); rets <- list(); strat <- list()
s <- run_strategy(W_naive); res[["Naive"]]  <- perf(s$r, s$TO, s$TTO); rets[["Naive"]]  <- s$r; strat[["Naive"]] <- s
s <- run_strategy(W_mv);    res[["MinVar"]] <- perf(s$r, s$TO, s$TTO); rets[["MinVar"]] <- s$r; strat[["MinVar"]] <- s

btc_col <- which(coins == "BTC")
r_btc <- R[oos_idx, btc_col] / 100
res[["Bitcoin"]] <- perf(r_btc, 0, 0); rets[["Bitcoin"]] <- r_btc
strat[["Bitcoin"]] <- list(r = r_btc, to_series = 0, to_day = 1L)

# ---- MQR portfolios from step 06 ------------------------------------------
files <- list.files(MATLAB_DIR, pattern = "^weights_.*\\.txt$")
if (length(files) == 0) stop("no weight files found — run step 06 first")
info <- data.table(file = files)
info[, tag   := sub("^weights_(.*)\\.txt$", "\\1", file)]
info[, model := sub("_c\\d+$", "", tag)]
info[, c     := as.numeric(sub("^.*_c(\\d+)$", "\\1", tag)) / 10]
# frequency is encoded as _r<days>; its absence means the 7-day baseline
info[, freq := 7]
info[grepl("_r\\d+", model), freq := as.numeric(sub("^.*_r(\\d+).*$", "\\1", model))]
info[, model := sub("_r\\d+", "", model)]
setorder(info, model, freq, c)
cat(sprintf("found %d weight files: %s | frequencies: %s\n", nrow(info),
            paste(unique(info$model), collapse = ", "),
            paste(sort(unique(info$freq)), collapse = ", ")))

# benchmarks at every frequency present (rebalancing dates differ by frequency)
for (f in setdiff(sort(unique(info$freq)), A_WIN)) {
  np <- nper_of(f)
  s <- run_strategy(matrix(1/N, nrow = np, ncol = N), f)
  nm <- sprintf("Naive r%d", f)
  res[[nm]] <- perf(s$r, s$TO, s$TTO); rets[[nm]] <- s$r; strat[[nm]] <- s
}

wsum <- list(); conc <- list(); wmat <- list()
for (i in 1:nrow(info)) {
  Wm <- as.matrix(read.table(file.path(MATLAB_DIR, info$file[i]), sep = ","))
  f  <- info$freq[i]
  if (nrow(Wm) != nper_of(f)) { cat("skip (wrong rows):", info$file[i], "\n"); next }
  colnames(Wm) <- coins
  nm <- sprintf("%s%s c=%.1f", info$model[i],
                if (f == 7) "" else sprintf(" r%d", f), info$c[i])
  s  <- run_strategy(Wm, f)
  res[[nm]]  <- perf(s$r, s$TO, s$TTO)
  rets[[nm]] <- s$r
  strat[[nm]] <- s
  wsum[[nm]] <- colMeans(Wm)
  mx <- apply(Wm, 1, max)
  conc[[nm]] <- c(maxw = mean(mx), hhi = mean(rowSums(Wm^2)), single = mean(mx > 0.9),
                  stable = mean(rowSums(Wm[, intersect(c("USDC","USDT"), coins), drop = FALSE])),
                  effN = mean(1 / rowSums(Wm^2)))
  wmat[[nm]] <- Wm
}

tbl <- rbindlist(lapply(names(res), function(n) cbind(Portfolio = n, res[[n]])))
tbl <- as.data.frame(tbl)
num <- sapply(tbl, is.numeric); tbl[num] <- lapply(tbl[num], function(x) round(x, 4))
cat("\n=== Portfolio performance (out-of-sample) ===\n"); print(tbl, row.names = FALSE)

saveRDS(list(table = tbl, returns = rets, wmean = wsum, conc = conc, dates = dates[oos_idx]),
        "data/portfolio_results.rds")

# ---- LaTeX table -----------------------------------------------------------
dir.create("tables", showWarnings = FALSE)
tex <- c("\\begin{table}[htbp]", "\\centering",
  "\\caption{Out-of-sample portfolio performance. $c$ is the investor's",
  "risk-return preference in $\\min_w ES(w)-c\\,EU(w)$.}",
  "\\label{tab:portfolio}", "\\begin{tabular}{lrrrrrrr}", "\\toprule",
  "Portfolio & TO & TTO & CW & SD & MD & SR & CR \\\\", "\\midrule")
for (i in 1:nrow(tbl)) {
  x <- tbl[i, ]
  tex <- c(tex, sprintf("%s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
    gsub("_", " ", x$Portfolio), x$TO, x$TTO, x$CW, x$SD, x$MD, x$SR, x$CR))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, "tables/portfolio_performance.tex")

# ---- figures (written to figures/ and copied into FIG_DIR) ----------------
dir.create("figures", showWarnings = FALSE)
dts  <- dates[oos_idx]
MODS <- c("RW125","GARCH-Skt","Hybrid","GARCH-FZ")
mcol <- c("RW125"="#1f77b4","GARCH-Skt"="#d62728","Hybrid"="#2ca02c","GARCH-FZ"="#9467bd")
cw   <- function(n) cumprod(1 + rets[[n]])
disp <- function(x) gsub("RW125", "RW-125", x, fixed = TRUE)
fig_files <- character(0)

## (a) c-frontier: CW, SR and MD against c for each model, benchmarks as lines
png("figures/c_frontier.png", width = 1600, height = 620, res = 130)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.6, 0.8), mgp = c(2.3, 0.7, 0))
for (metric in c("CW", "SR", "MD")) {
  vals <- lapply(MODS, function(m) sapply(seq(0, 1, 0.1), function(cc) {
    p <- res[[sprintf("%s c=%.1f", m, cc)]]; if (is.null(p)) NA else p[[metric]] }))
  yl <- range(c(unlist(vals), res[["Naive"]][[metric]], res[["Bitcoin"]][[metric]]), na.rm = TRUE)
  plot(NA, xlim = c(0, 1), ylim = yl, xlab = "c  (weight on the upper tail)", ylab = metric,
       main = c(CW = "Cumulative wealth", SR = "Sharpe ratio", MD = "Maximum drawdown")[metric])
  abline(h = res[["Naive"]][[metric]],   col = "black",   lty = 2, lwd = 1.5)
  abline(h = res[["Bitcoin"]][[metric]], col = "orange3", lty = 3, lwd = 1.5)
  for (j in seq_along(MODS)) lines(seq(0, 1, 0.1), vals[[j]], col = mcol[MODS[j]],
                                   lwd = 2, type = "b", pch = 16, cex = 0.65)
  if (metric == "CW") legend("topleft", c(disp(MODS), "Naive", "Bitcoin"),
      col = c(mcol, "black", "orange3"), lty = c(1, 1, 1, 1, 2, 3), lwd = 2, bty = "n", cex = 0.78)
  grid(col = "grey92")
}
dev.off(); fig_files <- c(fig_files, "figures/c_frontier.png")

## (b) cumulative wealth by model for selected c, unconstrained and kappa = 0.3
plot_cw_panel <- function(nms, labs, cols, ttl) {
  ok <- nms %in% names(rets); nms <- nms[ok]; labs <- labs[ok]; cols <- cols[ok]
  yl <- range(unlist(lapply(c(nms, "Naive", "Bitcoin"), cw)))
  plot(dts, cw("Naive"), type = "n", log = "y", ylim = yl, xlab = "", ylab = "Cumulative wealth",
       main = ttl)
  abline(h = 1, col = "grey80", lty = 2)
  lines(dts, cw("Naive"), col = "black", lwd = 1.8, lty = 2)
  lines(dts, cw("Bitcoin"), col = "orange3", lwd = 1.8, lty = 3)
  for (i in seq_along(nms)) lines(dts, cw(nms[i]), col = cols[i], lwd = 1.7)
  legend("topleft", c(labs, "Naive", "Bitcoin"), col = c(cols, "black", "orange3"),
         lty = c(rep(1, length(nms)), 2, 3), lwd = 1.8, bty = "n", cex = 0.68)
  grid(col = "grey92")
}
cs <- c(0.0, 0.5, 0.8, 1.0); ccol <- c("grey55", "#1f77b4", "#d62728", "#9467bd")
png("figures/cumulative_wealth_models.png", width = 1500, height = 1100, res = 125)
par(mfrow = c(2, 2), mar = c(3.4, 4, 2.4, 0.8), mgp = c(2.3, 0.7, 0))
for (m in MODS) plot_cw_panel(sprintf("%s c=%.1f", m, cs), sprintf("c = %.1f", cs), ccol, disp(m))
dev.off(); fig_files <- c(fig_files, "figures/cumulative_wealth_models.png")

png("figures/cumulative_wealth_cap.png", width = 1500, height = 1100, res = 125)
par(mfrow = c(2, 2), mar = c(3.4, 4, 2.4, 0.8), mgp = c(2.3, 0.7, 0))
for (m in MODS) {
  best <- function(tag) { nn <- grep(paste0("^", tag, " c="), names(res), value = TRUE)
                          if (!length(nn)) NA else nn[which.max(sapply(nn, function(n) res[[n]]$SR))] }
  nms <- c(best(m), best(paste0(m, "_cap05")), best(paste0(m, "_cap03")))
  labs <- c(sprintf("none, %s", sub(".* ", "", nms[1])), sprintf("kappa = 0.5, %s", sub(".* ", "", nms[2])),
            sprintf("kappa = 0.3, %s", sub(".* ", "", nms[3])))
  plot_cw_panel(nms[!is.na(nms)], labs[!is.na(nms)], c("#d62728", "#1f77b4", "#2ca02c")[!is.na(nms)], disp(m))
}
dev.off(); fig_files <- c(fig_files, "figures/cumulative_wealth_cap.png")

## (c) portfolio weights over time
pal <- c(BNB = "grey60", BTC = "red", DOGE = "green3", ETH = "orange", SOL = "purple",
         TRX = "black", USDC = "gold", USDT = "darkred", XRP = "pink")
rdates <- dates[W_EST + (0:(NPER - 1)) * A_WIN + 1]
stack_plot <- function(W, ttl) {
  plot(NA, xlim = range(rdates), ylim = c(0, 1), xlab = "", ylab = "weight", main = ttl,
       cex.main = 0.95, xaxs = "i", yaxs = "i", xaxt = "n")
  axis.Date(1, at = seq(as.Date(format(min(rdates), "%Y-%m-01")), max(rdates), by = "6 months"),
            format = "%b %Y", cex.axis = 0.8)
  cum <- cbind(0, t(apply(W, 1, cumsum)))
  for (j in seq_len(ncol(W))) polygon(c(rdates, rev(rdates)), c(cum[, j], rev(cum[, j + 1])),
                                      col = pal[coins[j]], border = NA)
}
weights_fig <- function(file, keys, titles, nrow, ncol) {
  png(file, width = 1500, height = 380 * nrow + 120, res = 125)
  layout(rbind(matrix(seq_len(nrow * ncol), nrow, ncol, byrow = TRUE), nrow * ncol + 1),
         heights = c(rep(1, nrow), 0.25))
  par(mar = c(2.2, 3.4, 2, 0.6), mgp = c(2.1, 0.6, 0))
  for (i in seq_along(keys)) {
    if (!is.null(wmat[[keys[i]]])) stack_plot(wmat[[keys[i]]], titles[i]) else plot.new()
  }
  for (i in seq_len(nrow * ncol - length(keys))) plot.new()
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = coins, fill = pal[coins], horiz = TRUE, bty = "n", cex = 1.0)
  dev.off()
}
weights_fig("figures/weights_models.png",
            as.vector(rbind(sprintf("%s c=0.8", MODS), sprintf("%s_cap03 c=0.8", MODS))),
            as.vector(rbind(sprintf("%s, unconstrained, c = 0.8", disp(MODS)),
                            sprintf("%s, kappa = 0.3, c = 0.8", disp(MODS)))), 4, 2)
weights_fig("figures/weights_rw125.png", sprintf("RW125 c=%.1f", seq(0, 1, 0.1)),
            sprintf("RW-125, c = %.1f", seq(0, 1, 0.1)), 4, 3)
fig_files <- c(fig_files, "figures/weights_models.png", "figures/weights_rw125.png")

file.copy(fig_files, FIG_DIR, overwrite = TRUE)
cat("\nFigures written and copied to", FIG_DIR, ":", paste(basename(fig_files), collapse = ", "), "\n")

# ---- paper-style tables written into the LaTeX project ---------------------
MODS    <- c("RW125","GARCH-Skt","Hybrid","GARCH-FZ")
BMK     <- c("Naive","MinVar","Bitcoin")
# p-value of the Ledoit-Wolf (2008) test of equal Sharpe ratios against Naive
psr <- function(nm, ref = "Naive") {
  if (is.null(rets[[nm]]) || nm == ref) return(NA)
  as.numeric(sharpe_diff_test(rets[[nm]], rets[[ref]])["p"])
}
pfmt <- function(p) if (is.na(p)) "--" else if (p < 0.001) "0.000" else sprintf("%.3f", p)
fmtrow  <- function(lab, p, nm = NA) sprintf("%-22s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %s \\\\",
                                    lab, p$TO, p$TTO, p$CW, p$SD, p$MD, p$SR, p$CR,
                                    pfmt(if (is.na(nm)) NA else psr(nm)))
tex_tab <- function(rows, label, caption, notes,
                    hdr = "Portfolio& TO & TTO & CW & SD & MD & SR & CR \\\\",
                    align = "lrrrrrrr") {
  c("\\begin{table}[ht]", "", "\\centering", "",
    paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"), "",
    "\\begin{center}", paste0("\\begin{tabular}{", align, "}"), "\\toprule", hdr, "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
    paste0("\\textit{Notes: ", notes, "}"), "", "\\end{table}")
}

if (dir.exists(TEX_DIR)) {
  avail <- MODS[sapply(MODS, function(m) any(grepl(paste0("^", m, " c="), names(res))))]

  ## main table: benchmarks + c in {0, 0.5, 1} + best-Sharpe c per model
  rows <- unlist(lapply(BMK, function(b) fmtrow(b, res[[b]], b)))
  rows <- c(rows, "\\hline")
  for (m in avail) {
    nms <- grep(paste0("^", m, " c="), names(res), value = TRUE)
    cv  <- as.numeric(sub(".*c=", "", nms))
    for (cc in c(0, 0.5, 1)) {
      k <- nms[which(abs(cv - cc) < 1e-9)]
      if (length(k)) rows <- c(rows, fmtrow(sprintf("%s $c=%.1f$", m, cc), res[[k]], k))
    }
    kb <- nms[which.max(sapply(nms, function(n) res[[n]]$SR))]
    rows <- c(rows, fmtrow(sprintf("%s $c^*=%.1f$", m, cv[which.max(sapply(nms, function(n) res[[n]]$SR))]),
                           res[[kb]], kb))
    if (m != tail(avail,1)) rows <- c(rows, "\\hline")
  }
  writeLines(tex_tab(rows, "tab:Performance", hdr = "Portfolio & TO & TTO & CW & SD & MD & SR & CR & $p_{SR}$ \\\\",
    align = "lrrrrrrrr",
    "Out-of-sample portfolio performance measures.",
    paste0("TO and TTO are turnover and target turnover, CW cumulative wealth, SD the annualised ",
      "standard deviation, MD the maximum drawdown, SR the Sharpe ratio and CR the Calmar ratio. ",
      "$c$ is the quantile proportion coefficient of Equation (\\ref{eq:MQRmin}); $c=0$ is the pure ",
      "minimum-ES portfolio. $c^*$ denotes the value of $c$ maximising the Sharpe ratio for that ",
      "model and is selected ex post. $p_{SR}$ is the p-value of the test of \\citet{Ledoit2008} ",
      "of equal Sharpe ratios against the naive portfolio. Weekly rebalancing, no constraint on ",
      "individual weights.")),
    file.path(TEX_DIR, "Performance.tex"))

  ## capped vs uncapped, one block per model, weekly frequency
  caps_present <- sort(unique(as.numeric(sub("^.*_cap(\\d+)$", "\\1",
                     grep("_cap\\d+$", sub(" c=.*", "", names(res)), value = TRUE))) / 10))
  if (length(caps_present)) {
    getc <- function(m, kap, cc) {
      nm <- if (is.na(kap)) sprintf("%s c=%.1f", m, cc)
            else sprintf("%s_cap%02d c=%.1f", m, round(kap*10), cc)
      res[[nm]]
    }
    rows <- character(0)
    for (m in avail) {
      has <- sapply(caps_present, function(k) !is.null(getc(m, k, 0.7)))
      if (!any(has)) next
      rows <- c(rows, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\",
                              1 + 2*(1+sum(has)), m))
      for (cc in seq(0, 1, 0.1)) {
        p0 <- getc(m, NA, cc)
        pk <- lapply(caps_present[has], function(k) getc(m, k, cc))
        if (is.null(p0) || any(sapply(pk, is.null))) next
        rows <- c(rows, sprintf("$c=%.1f$ & %s & %s \\\\", cc,
          paste(c(sprintf("%.4f", p0$CW), sapply(pk, function(p) sprintf("%.4f", p$CW))), collapse = " & "),
          paste(c(sprintf("%.4f", p0$SR), sapply(pk, function(p) sprintf("%.4f", p$SR))), collapse = " & ")))
      }
      rows <- c(rows, "\\hline")
    }
    nk  <- length(caps_present) + 1
    lab <- paste(c("none", sprintf("$\\kappa=%.1f$", caps_present)), collapse = " & ")
    writeLines(tex_tab(rows, "tab:PerformanceCap",
      "Effect of the concentration constraint on the MQR portfolios.",
      paste0("The table compares the unconstrained portfolios with those constrained by Equation ",
        "(\\ref{eq:cap}). CW is cumulative wealth and SR the Sharpe ratio. Weekly rebalancing."),
      hdr = paste0("& \\multicolumn{", nk, "}{c}{CW} & \\multicolumn{", nk, "}{c}{SR} \\\\\n",
                   sprintf("\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}\n", nk+1, nk+2, 2*nk+1),
                   " & ", lab, " & ", lab, " \\\\"),
      align = paste0("l", strrep("r", 2*nk))),
      file.path(TEX_DIR, "PerformanceCap.tex"))
  }

  ## summary of the concentration constraint: best c per model and cap
  best_of <- function(tag) {
    nn <- grep(paste0("^", tag, " c="), names(res), value = TRUE)
    if (!length(nn)) return(NULL)
    nn[which.max(sapply(nn, function(n) res[[n]]$SR))]
  }
  rows <- character(0)
  for (m in avail) {
    for (k in c("", "_cap05", "_cap03")) {
      nm <- best_of(paste0(m, k)); if (is.null(nm)) next
      p  <- res[[nm]]; cc <- sub(".*c=", "", nm)
      lab <- c("none", "0.5", "0.3")[match(k, c("", "_cap05", "_cap03"))]
      pv <- as.numeric(var_ratio_test(rets[[nm]], rets[["Naive"]])["p"])
      rows <- c(rows, sprintf("%-10s & %s & %s & %.2f & %.3f & %.3f & %.3f & %.3f & %.3f & %s & %s \\\\",
        if (k == "") m else "", lab, cc, conc[[nm]]["effN"], p$TO, p$CW, p$SD, p$MD, p$SR,
        pfmt(psr(nm)), pfmt(pv)))
    }
    if (m != tail(avail, 1)) rows <- c(rows, "\\hline")
  }
  rows <- c(rows, "\\hline", sprintf("%-10s & -- & -- & %.2f & %.3f & %.3f & %.3f & %.3f & %.3f & -- & -- \\\\",
    "Naive", length(coins), res$Naive$TO, res$Naive$CW, res$Naive$SD, res$Naive$MD, res$Naive$SR))
  writeLines(tex_tab(rows, "tab:CapSummary",
    "Portfolio performance with and without the concentration constraint.",
    paste0("For each forecasting model and weight cap $\\kappa$ of Equation (\\ref{eq:cap}) the ",
      "row reports the preference coefficient $c^*$ with the highest Sharpe ratio and the ",
      "corresponding performance. $N_{\\mathrm{eff}}$ is the effective number of assets, the inverse ",
      "Herfindahl index of the weights averaged over rebalancing dates. $p_{SR}$ is the p-value of ",
      "the \\citet{Ledoit2008} test of equal Sharpe ratios and $p_{SD}$ that of the corresponding ",
      "HAC test of equal variances \\citep{Ledoit2011}, both against the naive portfolio. Weekly rebalancing."),
    hdr = "Model & $\\kappa$ & $c^*$ & $N_{\\mathrm{eff}}$ & TO & CW & SD & MD & SR & $p_{SR}$ & $p_{SD}$ \\\\",
    align = "llrrrrrrrrr"), file.path(TEX_DIR, "CapSummary.tex"))

  ## concentration of the unconstrained portfolios
  rows <- character(0)
  for (m in avail) {
    f1 <- function(cc) conc[[sprintf("%s c=%.1f", m, cc)]]
    a <- f1(0.5); b <- f1(1.0)
    if (is.null(a) || is.null(b)) next
    rows <- c(rows, sprintf("%-10s & %.3f & %.1f & %.1f & %.2f & %.3f & %.1f & %.1f & %.2f \\\\", m,
      a["maxw"], 100 * a["single"], 100 * a["stable"], a["effN"],
      b["maxw"], 100 * b["single"], 100 * b["stable"], b["effN"]))
  }
  writeLines(tex_tab(rows, "tab:Concentration",
    "Concentration of the unconstrained MQR portfolios.",
    paste0("Averages over the 156 weekly rebalancing dates. Max $w$ is the largest single weight, ",
      "Single the percentage of dates on which one asset receives more than 90\\% of the portfolio, ",
      "Stable the percentage invested in the two stablecoins and $N_{\\mathrm{eff}}$ the effective ",
      "number of assets (inverse Herfindahl index)."),
    hdr = paste0("& \\multicolumn{4}{c}{$c = 0.5$} & \\multicolumn{4}{c}{$c = 1.0$} \\\\\n",
                 "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\n",
                 "Model & Max $w$ & Single (\\%) & Stable (\\%) & $N_{\\mathrm{eff}}$ & ",
                 "Max $w$ & Single (\\%) & Stable (\\%) & $N_{\\mathrm{eff}}$ \\\\"),
    align = "lrrrrrrrr"), file.path(TEX_DIR, "Concentration.tex"))

  ## net of transaction costs
  rows <- character(0)
  sel  <- c(BMK, unlist(lapply(avail, function(m) {
    nms <- grep(paste0("^", m, " c="), names(res), value = TRUE)
    nms[which.max(sapply(nms, function(n) res[[n]]$SR))] })))
  for (n in sel) {
    s <- strat[[n]]
    g <- res[[n]]
    p10 <- perf(net_returns(s, 10)); p25 <- perf(net_returns(s, 25))
    rows <- c(rows, sprintf("%-24s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
      gsub("_", " ", n), g$TO, g$CW, p10$CW, p25$CW, g$SR, p10$SR, p25$SR))
    if (n == tail(BMK,1)) rows <- c(rows, "\\hline")
  }
  writeLines(tex_tab(rows, "tab:PerformanceNet",
    "Portfolio performance net of proportional transaction costs.",
    paste0("Costs are charged at each rebalancing date in proportion to turnover. CW and SR are ",
      "reported gross and net of 10 and 25 basis points per unit of turnover. For each forecasting ",
      "model the value of $c$ maximising the gross Sharpe ratio is shown. Weekly rebalancing."),
    hdr = paste0("& & \\multicolumn{3}{c}{CW} & \\multicolumn{3}{c}{SR} \\\\\n",
                 "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\n",
                 "Portfolio & TO & gross & 10 bp & 25 bp & gross & 10 bp & 25 bp \\\\")),
    file.path(TEX_DIR, "PerformanceNet.tex"))

  ## ---- rebalancing frequency: summary table + heat map -------------------
  fr_all <- sort(unique(info$freq))
  if (length(fr_all) > 1) {
    getp <- function(m, f, cc) {
      nm <- sprintf("%s%s c=%.1f", m, if (f == 7) "" else sprintf(" r%d", f), cc)
      res[[nm]]
    }
    ## table: best-Sharpe c for each model x frequency
    rows <- character(0)
    for (m in MODS) {
      for (f in fr_all) {
        cand <- lapply(seq(0,1,0.1), function(cc) getp(m, f, cc))
        ok   <- !sapply(cand, is.null)
        if (!any(ok)) next
        srs  <- sapply(cand[ok], function(p) p$SR)
        bi   <- which.max(srs); bc <- seq(0,1,0.1)[ok][bi]; p <- cand[ok][[bi]]
        rows <- c(rows, sprintf("%-14s & %d & %.1f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
                  if (f == fr_all[1]) m else "", f, bc, p$TO, p$CW, p$SD, p$MD, p$SR))
      }
      rows <- c(rows, "\\hline")
    }
    for (f in fr_all) {
      nm <- if (f == 7) "Naive" else sprintf("Naive r%d", f)
      p  <- res[[nm]]
      if (!is.null(p)) rows <- c(rows, sprintf("%-14s & %d & --- & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
                                 if (f == fr_all[1]) "Naive" else "", f, p$TO, p$CW, p$SD, p$MD, p$SR))
    }
    writeLines(tex_tab(rows, "tab:PerformanceRebal",
      "Portfolio performance across rebalancing frequencies.",
      paste0("For each forecasting model and rebalancing frequency the table reports the value of ",
        "$c$ maximising the Sharpe ratio and the corresponding performance. All strategies are ",
        "evaluated over the same out-of-sample span. No constraint on individual weights."),
      hdr = "Model & Freq. & $c^*$ & TO & CW & SD & MD & SR \\\\",
      align = "llrrrrrr"),
      file.path(TEX_DIR, "PerformanceRebal.tex"))

    ## heat map of Sharpe ratios: rows = model x frequency, cols = c
    cs   <- seq(0, 1, 0.1)
    grid <- expand.grid(model = MODS, freq = fr_all, stringsAsFactors = FALSE)
    M    <- t(sapply(1:nrow(grid), function(i)
              sapply(cs, function(cc) { p <- getp(grid$model[i], grid$freq[i], cc)
                                        if (is.null(p)) NA else p$SR })))
    keep <- rowSums(!is.na(M)) > 0
    M <- M[keep, , drop = FALSE]; grid <- grid[keep, , drop = FALSE]
    rn <- sprintf("%s  %dd", gsub("RW125", "RW-125", grid$model, fixed = TRUE), grid$freq)

    png("figures/rebalancing_heatmap.png", width = 1250, height = 950, res = 130)
    par(mar = c(4.2, 8.5, 3.2, 5.5))
    pal <- colorRampPalette(c("#b2182b","#f7f7f7","#2166ac"))(101)
    rng <- range(M, na.rm = TRUE)
    image(x = seq_along(cs), y = seq_len(nrow(M)), z = t(M)[, nrow(M):1],
          col = pal, zlim = rng, axes = FALSE, xlab = "c  (upside preference)", ylab = "")
    axis(1, at = seq_along(cs), labels = sprintf("%.1f", cs), cex.axis = 0.8)
    axis(2, at = seq_len(nrow(M)), labels = rev(rn), las = 2, cex.axis = 0.75)
    title("Sharpe ratio by model, rebalancing frequency and preference c")
    for (i in seq_len(nrow(M))) for (j in seq_along(cs)) {
      v <- M[nrow(M) - i + 1, j]
      if (!is.na(v)) text(j, i, sprintf("%.2f", v), cex = 0.55,
                          col = if (v > mean(rng)) "white" else "black")
    }
    box()
    lv <- pretty(rng, 4)
    par(new = TRUE, mar = c(4.2, 0, 3.2, 1.5), fig = c(0.9, 0.97, 0.12, 0.88))
    image(1, seq(rng[1], rng[2], length.out = 101), matrix(1:101, 1), col = pal, axes = FALSE,
          xlab = "", ylab = "")
    axis(4, at = lv, labels = lv, las = 2, cex.axis = 0.7); box()
    dev.off()
    file.copy("figures/rebalancing_heatmap.png", FIG_DIR, overwrite = TRUE)
    cat("figures/rebalancing_heatmap.png written (and copied to", FIG_DIR, ")\n")
  }

  for (f in file.path(TEX_DIR, c("Performance.tex", "CapSummary.tex", "Concentration.tex",
                                  "PerformanceCap.tex", "PerformanceRebal.tex", "PerformanceNet.tex")))
    if (file.exists(f)) writeLines(gsub("RW125", "RW-125", readLines(f), fixed = TRUE), f)
  cat("paper-style portfolio tables written to", TEX_DIR, "\n")
}
