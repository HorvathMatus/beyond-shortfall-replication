# ============================================================
# Beyond Shortfall — 13: historical two-sided benchmarks from the literature
#
# The MQR portfolios of steps 06-07 use FORECASTS of both tails. The
# established alternatives in the portfolio literature use HISTORICAL returns.
# This step builds three of them on exactly the same rolling scheme (1095-day
# window, weekly rebalancing, tau = 0.05, long-only, optional weight cap):
#
#   Min-CVaR  minimum conditional value-at-risk, Rockafellar and Uryasev
#             (2000); exact linear programme. This is the historical
#             counterpart of the MQR portfolio at c = 0.
#   Omega     maximum Omega ratio at threshold 0, Keating and Shadwick (2002);
#             exact linear programme after the Charnes-Cooper transformation.
#   Rachev    maximum Rachev ratio EU/ES, Biglova et al. (2004); the problem
#             is not convex, so it is solved by Nelder-Mead from many starting
#             points (equal weights, all feasible vertices and random draws).
#
# Inputs : data/portfolio_input.rds
# Outputs: data/reward_risk_benchmarks.rds, LaTeX/tables/RewardRisk.tex
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages(library(lpSolve))
source("R/sharpe_test.R")

W_EST <- 1095; A_WIN <- 7; PPY <- 365; TAU <- 0.05
OOS_N <- 1080                                  # common span of step 07
CAPS  <- c(none = 1, "0.3" = 0.3)
set.seed(2024)

inp <- readRDS("data/portfolio_input.rds")
R <- inp$ret / 100; colnames(R) <- inp$coins
N <- ncol(R); NPER <- floor((nrow(R) - W_EST) / A_WIN)
window_of <- function(l) R[((l - 1) * A_WIN + 1):(W_EST + (l - 1) * A_WIN), , drop = FALSE]

# ---- Min-CVaR: Rockafellar-Uryasev LP --------------------------------------
# variables: w (N), zeta+ , zeta- , u (T);  min zeta + 1/(tau T) sum u
# s.t. u_t + R_t'w + zeta >= 0,  sum w = 1,  w <= kappa
min_cvar <- function(X, kappa) {
  T <- nrow(X)
  obj <- c(rep(0, N), 1, -1, rep(1 / (TAU * T), T))
  A1  <- cbind(X, 1, -1, diag(T))                         # >= 0
  A2  <- c(rep(1, N), 0, 0, rep(0, T))                     # == 1
  A3  <- cbind(diag(N), 0, 0, matrix(0, N, T))             # <= kappa
  sol <- lp("min", obj, rbind(A1, A2, A3),
            c(rep(">=", T), "=", rep("<=", N)), c(rep(0, T), 1, rep(kappa, N)))
  w <- pmax(sol$solution[1:N], 0); w / sum(w)
}

# ---- Omega ratio (threshold 0): Charnes-Cooper LP ---------------------------
# max mu'z  s.t. (1/T) sum d = 1, d_t >= -R_t'z, d >= 0, sum z = s, z <= kappa s
# variables: z (N), s, d (T);  w = z / s
max_omega <- function(X, kappa) {
  T <- nrow(X); mu <- colMeans(X)
  obj <- c(mu, 0, rep(0, T))
  A1 <- c(rep(0, N), 0, rep(1 / T, T))                     # == 1
  A2 <- cbind(X, 0, diag(T))                               # d_t + R_t'z >= 0
  A3 <- c(rep(1, N), -1, rep(0, T))                        # sum z - s == 0
  A4 <- cbind(diag(N), -kappa, matrix(0, N, T))            # z - kappa s <= 0
  sol <- lp("max", obj, rbind(A1, A2, A3, A4),
            c("=", rep(">=", T), "=", rep("<=", N)), c(1, rep(0, T), 0, rep(0, N)))
  z <- sol$solution[1:N]; s <- sol$solution[N + 1]
  if (sol$status != 0 || s <= 0) return(rep(1 / N, N))
  w <- pmax(z / s, 0); w / sum(w)
}

# ---- Rachev ratio: multi-start Nelder-Mead ----------------------------------
tail_mags <- function(y) {
  k <- max(1, floor(TAU * length(y))); s <- sort(y)
  c(es = -mean(s[1:k]), eu = mean(s[(length(y) - k + 1):length(y)]))
}
max_rachev <- function(X, kappa) {
  softmax <- function(th) { e <- exp(th - max(th)); e / sum(e) }
  f <- function(th) {
    w <- softmax(th); pen <- 1e3 * sum(pmax(w - kappa, 0)^2)
    tm <- tail_mags(X %*% w)
    if (tm["es"] <= 0) return(1e6)
    -(tm["eu"] / tm["es"]) + pen
  }
  starts <- list(rep(0, N))
  kv <- ceiling(1 / kappa)
  if (choose(N, kv) <= 40) for (cm in combn(N, kv, simplify = FALSE)) {
    th <- rep(-4, N); th[cm] <- 0; starts[[length(starts) + 1]] <- th }
  for (i in 1:15) starts[[length(starts) + 1]] <- rnorm(N, sd = 1.5)
  best <- NULL; bv <- Inf
  for (th0 in starts) {
    o <- optim(th0, f, method = "Nelder-Mead", control = list(maxit = 1500, reltol = 1e-9))
    if (o$value < bv) { bv <- o$value; best <- o$par }
  }
  w <- softmax(best)
  if (kappa < 1) {                                  # project tiny cap violations
    for (it in 1:50) { ex <- pmax(w - kappa, 0); if (sum(ex) < 1e-10) break
      w <- pmin(w, kappa); free <- w < kappa; w[free] <- w[free] + sum(ex) * w[free] / sum(w[free]) }
  }
  w
}

# ---- weights for every rebalancing date -------------------------------------
crit <- list("Min-CVaR" = min_cvar, "Omega" = max_omega, "Rachev" = max_rachev)
W <- list()
t0 <- Sys.time()
for (k in names(CAPS)) for (cn in names(crit)) {
  nm <- sprintf("%s|%s", cn, k)
  W[[nm]] <- t(sapply(1:NPER, function(l) crit[[cn]](window_of(l), CAPS[[k]])))
  colnames(W[[nm]]) <- colnames(R)
  cat(sprintf("%-16s done (%.1f min)\n", nm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ---- out-of-sample evaluation (same accounting as step 07) ------------------
holds <- lapply(1:NPER, function(l) (W_EST + (l - 1) * A_WIN + 1):(W_EST + l * A_WIN))
run_strategy <- function(Wm) {
  r_out <- numeric(0); to <- numeric(0); w_drift <- NULL
  for (l in 1:NPER) {
    w <- Wm[l, ]; Rl <- R[holds[[l]], , drop = FALSE]
    if (!is.null(w_drift)) to <- c(to, sum(abs(w - w_drift)))
    wc <- w
    for (d in 1:nrow(Rl)) { rp <- sum(wc * Rl[d, ]); r_out <- c(r_out, rp); wc <- wc * (1 + Rl[d, ]) / (1 + rp) }
    w_drift <- wc
  }
  list(r = r_out[1:OOS_N], TO = mean(to))
}
perf <- function(r, TO) {
  cw <- cumprod(1 + r); sd_a <- sd(r) * sqrt(PPY); mu_a <- mean(r) * PPY; md <- max(1 - cw / cummax(cw))
  data.frame(TO = TO, CW = tail(cw, 1), SD = sd_a, MD = md, SR = mu_a / sd_a, CR = mu_a / md)
}
naive <- run_strategy(matrix(1 / N, NPER, N))
out <- do.call(rbind, lapply(names(W), function(nm) {
  s <- run_strategy(W[[nm]]); Wm <- W[[nm]]
  cbind(data.frame(criterion = sub("\\|.*", "", nm), cap = sub(".*\\|", "", nm),
    effN = mean(1 / rowSums(Wm^2)), stable = mean(rowSums(Wm[, c("USDC", "USDT")]))),
    perf(s$r, s$TO), pSR = sharpe_diff_test(s$r, naive$r)["p"],
    pSD = var_ratio_test(s$r, naive$r)["p"])
}))
rownames(out) <- NULL
pn <- perf(naive$r, naive$TO)
cat(sprintf("\nNaive: CW %.3f SD %.3f SR %.3f\n", pn$CW, pn$SD, pn$SR))
print(out, digits = 3)
saveRDS(list(results = out, weights = W, naive = pn), "data/reward_risk_benchmarks.rds")

# ---- LaTeX ------------------------------------------------------------------
pf <- function(p) if (is.na(p)) "--" else if (p < 0.001) "0.000" else sprintf("%.3f", p)
rows <- character(0)
for (k in names(CAPS)) {
  for (i in which(out$cap == k)) {
    x <- out[i, ]
    rows <- c(rows, sprintf("%-9s & %s & %.2f & %.1f & %.3f & %.3f & %.3f & %.3f & %.3f & %s & %s \\\\",
      x$criterion, k, x$effN, 100 * x$stable, x$TO, x$CW, x$SD, x$MD, x$SR, pf(x$pSR), pf(x$pSD)))
  }
  rows <- c(rows, "\\hline")
}
rows <- c(rows, sprintf("%-9s & -- & %.2f & %.1f & %.3f & %.3f & %.3f & %.3f & %.3f & -- & -- \\\\",
  "Naive", N, 100 * 2 / N, pn$TO, pn$CW, pn$SD, pn$MD, pn$SR))
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Historical reward--risk benchmarks.}", "\\label{tab:RewardRisk}", "\\begin{center}",
  "\\begin{tabular}{llrrrrrrrrr}", "\\toprule",
  "Criterion & $\\kappa$ & $N_{\\mathrm{eff}}$ & Stable (\\%) & TO & CW & SD & MD & SR & $p_{SR}$ & $p_{SD}$ \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: Portfolios estimated on historical returns over the same 1{,}095-day rolling ",
    "window, weekly rebalancing and $\\tau = 0.05$ as the MQR portfolios. Min-CVaR minimises the ",
    "historical conditional value-at-risk \\citep{Rockafellar2000}, Omega maximises the Omega ratio at ",
    "a zero threshold \\citep{Keating2002}, and Rachev maximises the ratio of the historical expected ",
    "tail gain to the expected tail loss \\citep{Biglova2004}. $N_{\\mathrm{eff}}$ is the effective number ",
    "of assets and Stable the average stablecoin weight. $p_{SR}$ and $p_{SD}$ are the p-values of the ",
    "\\citet{Ledoit2008, Ledoit2011} tests against the naive portfolio.}"), "\\end{table}")
writeLines(L, file.path(TEX_DIR, "RewardRisk.tex"))
cat("\nTable written:", file.path(TEX_DIR, "RewardRisk.tex"), "\n")
