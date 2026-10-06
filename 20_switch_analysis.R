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
#
# Inputs : weights_<model>_c<cc>.txt (step 06), data/portfolio_input.rds
# Outputs: data/switch_analysis.rds, LaTeX/tables/SwitchAnalysis.tex
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages(library(data.table))

W_EST <- 1095; A_WIN <- 7; WIN <- 30
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
  rbindlist(lapply(RISKY, function(j) as.data.table(as.list(c(asset = 0, pre_post(j, t))))))
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

# mean log volatility ratio relative to the benchmark, with a bootstrap interval
summ <- sw[, {
  d <- lr - mean(bench$lr)
  bs <- replicate(2000, mean(sample(d, replace = TRUE)))
  .(n = .N, calm = mean(calm_rank), pre = mean(pre), post = mean(post),
    excess = mean(d), lo = quantile(bs, 0.025), hi = quantile(bs, 0.975), ret7 = mean(ret7),
    ret7_lo = quantile(replicate(2000, mean(sample(ret7, replace = TRUE))), 0.025),
    ret7_hi = quantile(replicate(2000, mean(sample(ret7, replace = TRUE))), 0.975))
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
    sig <- x$lo > 0 | x$hi < 0
    e <- sprintf("%+.2f", x$excess); if (sig) e <- paste0("\\textbf{", e, "}")
    rows <- c(rows, sprintf("%-10s & %.1f & %d & %.1f & %.2f & %.2f & %s & [%+.2f, %+.2f] & %+.1f \\\\",
      m, cv, x$n, x$calm, x$pre, x$post, e, x$lo, x$hi, 100 * x$ret7))
  }
  rows <- c(rows, "\\hline")
}
rows <- rows[-length(rows)]
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Switches of the unconstrained MQR portfolios: volatility before and after entry.}",
  "\\label{tab:SwitchAnalysis}", "\\begin{center}", "\\small", "\\begin{tabular}{lrrrrrrcr}", "\\toprule",
  "Model & $c$ & Switches & Calm rank & SD pre & SD post & Excess & 95\\% CI & Ret$_7$ (\\%) \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: A switch is a weekly rebalancing date at which the asset with the largest weight, at least ",
    "50\\%, is a risky asset different from that of the previous date (unconstrained portfolios). Calm rank is the ",
    "average rank of the entered asset's volatility over the preceding ", WIN, " days among the seven risky assets ",
    "(1 = calmest). SD pre and SD post are its annualised volatilities over the ", WIN, " days before and after the ",
    "decision. Excess is the average of $\\log(\\mathrm{SD\\ post}/\\mathrm{SD\\ pre})$ minus its average over all risky ",
    "assets and all rebalancing dates (", sprintf("%.2f", mean(bench$lr)), "), which removes the mean reversion of ",
    "volatility common to all assets; the interval is a bootstrap percentile interval over switches, and bold ",
    "entries exclude zero. Ret$_7$ is the average return of the entered asset over the following seven days.}"),
  "\\end{table}")
writeLines(L, file.path(TEX_DIR, "SwitchAnalysis.tex"))
cat("\nTable written:", file.path(TEX_DIR, "SwitchAnalysis.tex"), "\n")
