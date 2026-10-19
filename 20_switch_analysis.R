# ============================================================
# Beyond Shortfall — 20: what happens after the unconstrained portfolios switch
#
# Section 6 conjectures that the responsive conditional models concentrate in
# an asset that has just become calm, shortly before its volatility rises.
# This step tests the conjecture on the weights of step 06 (unconstrained,
# weekly rebalancing). A switch is a rebalancing date at which the asset with
# the largest weight (at least 50%) is a risky asset that differs from the
# previous date's. For every switch we compare the volatility of the newly
# entered asset over the 30 days before the decision with its volatility over
# the 30 days after it, and record its return over the following week.
# The benchmark is the same ratio for every risky asset at every rebalancing
# date, which captures the mean reversion of volatility that any rule picking
# calm assets would face.
# Inference: the 30-day windows of nearby switches overlap and switches cluster
# in time, so the intervals come from a moving-block bootstrap over calendar
# time. Blocks of BLOCK consecutive rebalancing periods (63 days, the length of
# the pre- plus post-switch windows) are resampled, and the switches and the
# benchmark observations of the drawn periods enter the replication together.
#
# Inputs : weights_<model>_c<cc>.txt (step 06), data/portfolio_input.rds
# Outputs: data/switch_analysis.rds, LaTeX/tables/SwitchAnalysis.tex
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages(library(data.table))

W_EST <- 1095; A_WIN <- 7; WIN <- 30
BLOCK <- 9; B_BOOT <- 2000                          # block length in rebalancing periods (63 days)
MODS <- c("RW125" = "RW-125", "GARCH-Skt" = "GARCH-Skt", "GARCH-FZ" = "GARCH-FZ", "Hybrid" = "Hybrid")
CS <- c(0.8, 1.0)
set.seed(7)

inp <- readRDS("data/portfolio_input.rds")
R <- inp$ret / 100; colnames(R) <- inp$coins
RISKY <- setdiff(colnames(R), c("USDC", "USDT"))
NPER <- floor((nrow(R) - W_EST) / A_WIN)
dec <- W_EST + (0:(NPER - 1)) * A_WIN               # last day of each estimation window
vol <- function(j, rows) sd(R[rows, j]) * sqrt(365)
pre_post <- function(j, t) {
  pre <- (t - WIN + 1):t; post <- (t + 1):min(t + WIN, nrow(R))
  c(pre = vol(j, pre), post = vol(j, post), ret7 = prod(1 + R[(t + 1):min(t + 7, nrow(R)), j]) - 1)
}

# benchmark: every risky asset at every rebalancing date with 30 days ahead
bench <- rbindlist(lapply(seq_len(NPER), function(l) {
  t <- dec[l]; if (t + WIN > nrow(R)) return(NULL)
  rbindlist(lapply(RISKY, function(j) as.data.table(as.list(c(period = l, pre_post(j, t))))))
}))
bench[, lr := log(post / pre)]

sw <- rbindlist(lapply(names(MODS), function(m) rbindlist(lapply(CS, function(c) {
  W <- as.matrix(read.table(file.path(MATLAB_DIR, sprintf("weights_%s_c%02d.txt", m, round(c * 10))), sep = ","))
  colnames(W) <- colnames(R)
  dom <- apply(W, 1, function(w) if (max(w) >= 0.5) names(w)[which.max(w)] else NA_character_)
  rbindlist(lapply(2:NPER, function(l) {
    t <- dec[l]
    if (is.na(dom[l]) || !(dom[l] %in% RISKY) || identical(dom[l], dom[l - 1]) || t + WIN > nrow(R)) return(NULL)
    x <- pre_post(dom[l], t)
    # rank of the entered asset's pre-switch volatility among the risky assets (1 = calmest)
    rk <- rank(sapply(RISKY, function(j) vol(j, (t - WIN + 1):t)))[dom[l]]
    data.table(model = MODS[[m]], c = c, period = l, asset = dom[l], pre = x["pre"], post = x["post"],
               ret7 = x["ret7"], calm_rank = rk)
  }))
}))))
sw[, lr := log(post / pre)]

# ---- moving-block bootstrap over rebalancing periods --------------------------
# One set of resampled periods per replication, shared by all models and by the
# benchmark. cnt[b, l] = how often period l enters replication b.
NP_EFF <- max(bench$period)                         # last period with 30 days ahead
nb <- ceiling(NP_EFF / BLOCK)
cnt <- t(replicate(B_BOOT, {
  st <- sample.int(NP_EFF - BLOCK + 1, nb, replace = TRUE)
  tabulate(unlist(lapply(st, function(s) s:(s + BLOCK - 1)))[seq_len(NP_EFF)], NP_EFF)
}))
b_sum <- tapply(bench$lr, bench$period, sum)[as.character(seq_len(NP_EFF))]; b_sum[is.na(b_sum)] <- 0
b_n   <- tabulate(bench$period, NP_EFF)
bench_bs <- as.vector(cnt %*% b_sum) / as.vector(cnt %*% b_n)

# mean log volatility ratio relative to the benchmark and mean 7-day return
summ <- sw[, {
  w  <- cnt[, period, drop = FALSE]                 # weights of the switches in each replication
  nw <- rowSums(w)
  ex_bs  <- ifelse(nw > 0, as.vector(w %*% lr) / nw - bench_bs, NA)
  r7_bs  <- ifelse(nw > 0, as.vector(w %*% ret7) / nw, NA)
  .(n = .N, calm = mean(calm_rank), pre = mean(pre), post = mean(post),
    excess = mean(lr) - mean(bench$lr),
    lo = quantile(ex_bs, 0.025, na.rm = TRUE), hi = quantile(ex_bs, 0.975, na.rm = TRUE),
    ret7 = mean(ret7),
    ret7_lo = quantile(r7_bs, 0.025, na.rm = TRUE), ret7_hi = quantile(r7_bs, 0.975, na.rm = TRUE),
    empty = mean(nw == 0))
}, by = .(model, c)]
cat(sprintf("Benchmark (all risky assets, all dates): mean log(post/pre) = %.3f, n = %d\n\n",
            mean(bench$lr), nrow(bench)))
print(summ, digits = 3)
saveRDS(list(switches = sw, summary = summ, bench = bench), "data/switch_analysis.rds")

# ---- LaTeX ------------------------------------------------------------------
rows <- character(0)
for (cv in CS) {
  for (m in MODS) {
    x <- summ[model == m & abs(c - cv) < 1e-9]
    if (!nrow(x)) next
    e <- sprintf("%+.2f", x$excess); if (x$lo > 0 | x$hi < 0) e <- paste0("\\textbf{", e, "}")
    r <- sprintf("%+.1f", 100 * x$ret7); if (x$ret7_lo > 0 | x$ret7_hi < 0) r <- paste0("\\textbf{", r, "}")
    rows <- c(rows, sprintf("%-10s & %.1f & %d & %.1f & %.2f & %.2f & %s & [%+.2f, %+.2f] & %s & [%+.1f, %+.1f] \\\\",
      m, cv, x$n, x$calm, x$pre, x$post, e, x$lo, x$hi, r, 100 * x$ret7_lo, 100 * x$ret7_hi))
  }
  rows <- c(rows, "\\hline")
}
rows <- rows[-length(rows)]
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Switches of the unconstrained MQR portfolios: volatility before and after entry.}",
  "\\label{tab:SwitchAnalysis}", "\\begin{center}", "\\small", "\\begin{tabular}{lrrrrrrcrc}", "\\toprule",
  "Model & $c$ & Switches & Calm rank & SD pre & SD post & Excess & 95\\% CI & Ret$_7$ (\\%) & 95\\% CI \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: A switch is a weekly rebalancing date at which the asset with the largest weight, at least ",
    "50\\%, is a risky asset different from that of the previous date (unconstrained portfolios). Calm rank is the ",
    "average rank of the entered asset's volatility over the preceding ", WIN, " days among the seven risky assets ",
    "(1 = calmest). SD pre and SD post are its annualised volatilities over the ", WIN, " days before and after the ",
    "decision. Excess is the average of $\\log(\\mathrm{SD\\ post}/\\mathrm{SD\\ pre})$ minus its average over all risky ",
    "assets and all rebalancing dates (", sprintf("%.2f", mean(bench$lr)), "), which removes the mean reversion of ",
    "volatility common to all assets. Ret$_7$ is the average return of the entered asset over the following seven ",
    "days. The intervals are percentile intervals of a moving-block bootstrap over calendar time (", B_BOOT,
    " replications, blocks of ", BLOCK, " consecutive rebalancing periods, i.e.\\ ", BLOCK * A_WIN, " days, resampled ",
    "jointly for the switches and the benchmark), which allows for the overlap of the windows and the clustering ",
    "of switches; bold entries exclude zero. With few switches the intervals are imprecise and the table is ",
    "descriptive.}"),
  "\\end{table}")
writeLines(L, file.path(TEX_DIR, "SwitchAnalysis.tex"))
cat("\nTable written:", file.path(TEX_DIR, "SwitchAnalysis.tex"), "\n")
