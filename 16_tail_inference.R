# ============================================================
# Beyond Shortfall — 16: inference on the differences between the tails
#
# Steps 11 and 15 show that the tails rank the models differently. Differences
# in ranks can, however, arise from small and noisy loss differences. This step
# tests directly whether the RELATIVE accuracy of the models changes between
# the tails, for both estimation schemes (models estimated once, step 10, and
# re-estimated every 91 days, step 14):
#
#   1. Rank shift: for each model, the average over assets and quantile levels
#      of (rank in the upper tail - rank in the lower tail), by mean FZ0 loss.
#   2. Tail interaction for model pairs (A, B): on every day,
#        d_t = [L_t^low(A) - L_t^low(B)] - [L_t^up(A) - L_t^up(B)],
#      the change between the tails in A's loss advantage over B. FZ0 loss
#      differences are unit-free, so d_t can be pooled across assets and tails.
#      Its mean is tested asset by asset and level by level with a HAC t-test,
#      and pooled across assets and levels.
#      Alongside the interaction, the loss differences of A and B in each tail,
#      which show whether the better of the two models changes between the tails.
#   3. Spearman correlation between the ES and EU rankings (models other than
#      the Single Index), averaged over assets and levels (printed and saved; not
#      tabulated, because its percentile interval is biased downwards).
#
# Multiple testing: within each scheme the 11 rank shifts and the 4
# interactions form one family of 15 hypotheses. Bootstrap p-values are
# adjusted with the method of Holm (1979), and simultaneous 95% intervals are
# computed from the bootstrap distribution of the maximum studentised deviation.
#
# Confidence intervals come from a moving-block bootstrap that resamples blocks
# of calendar days jointly for all assets, tails, levels and models, so that
# the dependence across them is preserved (block length 20 days; 50 as a check).
#
# Inputs : forecasts of steps 10/11 (assets/) and 14/15 (assets_roll/)
# Outputs: data/tail_inference.rds, LaTeX/tables/TailInference.tex
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages({ library(data.table); library(sandwich) })
source("R/forecast_eval.R")

SCHEMES <- c(Rolling = file.path(MATLAB_DIR, "assets_roll"), Fixed = file.path(MATLAB_DIR, "assets"))
ASSETS  <- c("BNB", "BTC", "DOGE", "ETH", "SOL", "TRX", "XRP")
PAIRS   <- list(c("GARCH-Skt", "GARCH-FZ"), c("GARCH-Skt", "GARCH-EDF"),
                c("GARCH-Skt", "Hybrid"), c("GARCH-Skt", "RW-250"))
B_BOOT  <- 999
BLOCK   <- 20
set.seed(2026)

# ---- FZ0 losses on a common calendar ------------------------------------------
# L[[scheme]][[asset]][[k]] = list(low = P x 11, up = P x 11), with the calendar
# position of every row in pos[[scheme]][[asset]]
cal <- sort(unique(unlist(lapply(ASSETS, function(a)
  scan(file.path(SCHEMES[1], a, "oos_dates.txt"), quiet = TRUE)))))
load_scheme <- function(root) {
  out <- list(); pos <- list()
  for (a in ASSETS) {
    dates <- scan(file.path(root, a, "oos_dates.txt"), quiet = TRUE)
    r <- scan(file.path(root, a, "oos_ret.txt"), quiet = TRUE)
    pos[[a]] <- match(dates, cal)
    out[[a]] <- lapply(seq_along(TAUS), function(k) {
      sapply(c("low", "up"), function(tl) {
        y <- if (tl == "low") r else -r
        fc <- read_fc(file.path(root, a), tl, SUFFIX[k])
        list(sapply(MODEL_ORDER, function(m) fz_loss(y, fc$v[, m], fc$e[, m], TAUS[k])))
      })
    })
  }
  list(L = out, pos = pos)
}
DAT <- lapply(SCHEMES, load_scheme)

# ---- statistics for one set of calendar positions -----------------------------
# row of each asset for every calendar position (NA if the asset has no data)
make_lookup <- function(D) lapply(D$pos, function(p) { m <- rep(NA_integer_, length(cal)); m[p] <- seq_along(p); m })

stats_fast <- function(D, lk, cal_idx) {
  ranks <- list(); inter <- list(); rhos <- c()
  for (a in ASSETS) {
    rows <- lk[[a]][cal_idx]; rows <- rows[!is.na(rows)]
    if (length(rows) < 50) next
    for (k in seq_along(TAUS)) {
      mlo <- colMeans(D$L[[a]][[k]]$low[rows, , drop = FALSE])
      mup <- colMeans(D$L[[a]][[k]]$up[rows, , drop = FALSE])
      rlo <- rank(ifelse(is.finite(mlo), mlo, NA), na.last = "keep", ties.method = "first")
      rup <- rank(ifelse(is.finite(mup), mup, NA), na.last = "keep", ties.method = "first")
      ranks[[length(ranks) + 1]] <- rup - rlo
      keep <- setdiff(MODEL_ORDER, "Single Index"); ok <- keep[is.finite(rlo[keep]) & is.finite(rup[keep])]
      rhos <- c(rhos, cor(rank(mlo[ok]), rank(mup[ok]), method = "spearman"))
      inter[[length(inter) + 1]] <- c(
        sapply(PAIRS, function(p) (mlo[p[1]] - mlo[p[2]]) - (mup[p[1]] - mup[p[2]])),
        sapply(PAIRS, function(p) mlo[p[1]] - mlo[p[2]]),
        sapply(PAIRS, function(p) mup[p[1]] - mup[p[2]]))
    }
  }
  pn <- sapply(PAIRS, paste, collapse = "|")
  c(setNames(rowMeans(do.call(cbind, ranks), na.rm = TRUE), paste0("shift:", MODEL_ORDER)),
    setNames(rowMeans(do.call(cbind, inter), na.rm = TRUE), c(paste0("inter:", pn), paste0("dlow:", pn), paste0("dup:", pn))),
    rho = mean(rhos))
}

mbb_index <- function(n, l) {
  starts <- sample.int(n - l + 1, ceiling(n / l), replace = TRUE)
  as.vector(sapply(starts, function(s) s:(s + l - 1)))[1:n]
}

res <- list()
for (sc in names(SCHEMES)) {
  D <- DAT[[sc]]; lk <- make_lookup(D)
  est <- stats_fast(D, lk, seq_along(cal))
  for (l in c(BLOCK, 50)) {
    t0 <- Sys.time()
    boot <- t(replicate(B_BOOT, stats_fast(D, lk, mbb_index(length(cal), l))))
    cat(sprintf("%s, block %d: %d replications in %.1f s\n", sc, l, B_BOOT,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    # bootstrap p-values (percentile) and the family of 15 tests: rank shifts and interactions
    pval <- apply(boot, 2, function(b) min(1, 2 * min(mean(b <= 0, na.rm = TRUE), mean(b >= 0, na.rm = TRUE))))
    se   <- apply(boot, 2, sd, na.rm = TRUE)
    fam  <- grepl("^shift:|^inter:", names(est))    # 11 rank shifts + 4 interactions
    p_holm <- rep(NA_real_, length(est)); p_holm[fam] <- p.adjust(pval[fam], method = "holm")
    tmax <- apply(abs(sweep(boot[, fam, drop = FALSE], 2, est[fam])) / matrix(se[fam], nrow(boot), sum(fam), byrow = TRUE), 1, max, na.rm = TRUE)
    cval <- quantile(tmax, 0.95)
    res[[paste(sc, l)]] <- data.table(scheme = sc, block = l, stat = names(est), est = est,
      lo95 = apply(boot, 2, quantile, 0.025, na.rm = TRUE), hi95 = apply(boot, 2, quantile, 0.975, na.rm = TRUE),
      p = pval, p_holm = p_holm, slo = ifelse(fam, est - cval * se, NA), shi = ifelse(fam, est + cval * se, NA))
  }
}
res <- rbindlist(res)

# ---- case-by-case HAC tests of the tail interaction ---------------------------
hac <- rbindlist(lapply(names(SCHEMES), function(sc) {
  D <- DAT[[sc]]
  rbindlist(lapply(ASSETS, function(a) rbindlist(lapply(seq_along(TAUS), function(k) {
    lo <- D$L[[a]][[k]]$low; up <- D$L[[a]][[k]]$up
    rbindlist(lapply(PAIRS, function(p) {
      d <- (lo[, p[1]] - lo[, p[2]]) - (up[, p[1]] - up[, p[2]])
      t <- if (all(is.finite(d))) mean(d) / sqrt(lrvar(d, type = "Andrews", prewhite = TRUE)) else NA_real_
      data.table(scheme = sc, asset = a, tau = TAUS[k], pair = paste(p, collapse = "|"), mean = mean(d), t = t)
    }))
  }))))
}))
hac_sum <- hac[, .(neg = sum(t < -1.96, na.rm = TRUE), pos = sum(t > 1.96, na.rm = TRUE), n = sum(!is.na(t))),
               by = .(scheme, pair)]

cat("\n=== Bootstrap results (block 20) ===\n"); print(res[block == BLOCK], digits = 3)
cat("\n=== Block 50 for comparison ===\n"); print(res[block == 50], digits = 3)
cat("\n=== HAC interaction tests (t < -1.96 : A relatively better in the lower tail) ===\n"); print(hac_sum)
saveRDS(list(boot = res, hac = hac, hac_sum = hac_sum), "data/tail_inference.rds")

# ---- LaTeX ------------------------------------------------------------------
num <- function(x, d) sprintf(paste0("%+.", d, "f"), x)
pf  <- function(p) if (is.na(p)) "" else if (p < 0.001) "0.000" else sprintf("%.3f", p)
cell <- function(sc, st, d = 2, holm = TRUE) {
  x <- res[scheme == sc & block == BLOCK & stat == st]
  sig <- if (holm) isTRUE(x$p_holm < 0.05) else (x$lo95 > 0 | x$hi95 < 0)
  e <- num(x$est, d); if (sig) e <- paste0("\\textbf{", e, "}")
  sprintf("%s & [%s, %s] & %s", e, num(x$lo95, d), num(x$hi95, d), if (holm) pf(x$p_holm) else "")
}
rowsA <- unlist(lapply(MODEL_ORDER, function(m) {
  r <- sprintf("%-12s & %s & %s \\\\", m, cell("Rolling", paste0("shift:", m)), cell("Fixed", paste0("shift:", m)))
  if (m %in% c("RW-500", "GARCH-FZ", "Hybrid")) c(r, "\\hline") else r }))
rowsB <- unlist(lapply(seq_along(PAIRS), function(i) {
  p <- PAIRS[[i]]; pn <- paste(p, collapse = "|")
  c(sprintf("\\multicolumn{7}{l}{%s vs %s} \\\\", p[1], p[2]),
    sprintf("\\quad lower tail & %s & %s \\\\", cell("Rolling", paste0("dlow:", pn), 3, FALSE), cell("Fixed", paste0("dlow:", pn), 3, FALSE)),
    sprintf("\\quad upper tail & %s & %s \\\\", cell("Rolling", paste0("dup:", pn), 3, FALSE), cell("Fixed", paste0("dup:", pn), 3, FALSE)),
    sprintf("\\quad interaction & %s & %s \\\\", cell("Rolling", paste0("inter:", pn), 3), cell("Fixed", paste0("inter:", pn), 3)),
    if (i < length(PAIRS)) "\\hline")
}))
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Inference on the differences between the lower- and upper-tail model comparisons.}",
  "\\label{tab:TailInference}", "\\begin{center}", "\\small", "\\begin{tabular}{lrcrrcr}", "\\toprule",
  "& \\multicolumn{3}{c}{Rolling} & \\multicolumn{3}{c}{Fixed} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "& Estimate & 95\\% CI & $p_{\\mathrm{Holm}}$ & Estimate & 95\\% CI & $p_{\\mathrm{Holm}}$ \\\\", "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Panel A: average rank shift (EU rank $-$ ES rank)}} \\\\",
  rowsA, "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Panel B: FZ loss difference of A minus B, by tail, and its change between the tails}} \\\\",
  rowsB,
  "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: Results pooled over the seven assets and four quantile levels, for the models re-estimated ",
    "every 91 days (Rolling) and estimated once (Fixed). Panel A: average over assets and levels of a model's ",
    "rank by mean FZ$^{(0)}$ loss in the upper tail minus its rank in the lower tail; a positive value indicates ",
    "a model that is relatively worse in the upper tail. Panel B: average difference between the mean FZ$^{(0)}$ ",
    "losses of A and B in each tail (negative: A more accurate) and the interaction ",
    "$[L^{low}_t(A) - L^{low}_t(B)] - [L^{up}_t(A) - L^{up}_t(B)]$ (negative: A's advantage over B is larger in ",
    "the lower tail). Confidence intervals are percentile intervals from ", B_BOOT, " moving-block bootstrap ",
    "replications with blocks of ", BLOCK, " calendar days resampled jointly across assets, tails, levels and models; ",
    "blocks of 50 days give similar intervals. $p_{\\mathrm{Holm}}$: bootstrap p-values adjusted with the method of ",
    "\\citet{Holm1979} within the family of the 15 rank shifts and interactions of each scheme. Bold: ",
    "$p_{\\mathrm{Holm}} < 0.05$ for the rank shifts and interactions; an unadjusted 95\\% interval that excludes ",
    "zero for the losses in each tail. Simultaneous 95\\% intervals for the family (max-$t$ bootstrap) lead to the ",
    "same conclusions, except that they also exclude zero for the rank shift of GARCH-Skt under the fixed scheme. ",
    "Because resampling adds noise to the rankings and so pulls rank shifts towards zero, the inference in Panel A ",
    "is conservative.}"),
  "\\end{table}")
writeLines(L, file.path(TEX_DIR, "TailInference.tex"))
cat("\nTable written:", file.path(TEX_DIR, "TailInference.tex"), "\n")
