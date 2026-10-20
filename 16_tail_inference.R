# ============================================================
# Beyond Shortfall — 16: inference on the differences between the tails
#
# Steps 11 and 15 show that the tails rank the models differently. Differences
# in ranks can, however, arise from small and noisy loss differences. This step
# tests directly whether the RELATIVE accuracy of the models changes between
# the tails, for both estimation schemes (models estimated once, step 10, and
# re-estimated every 91 days, step 14), and whether the change differs between
# the schemes:
#
#   1. Loss differences of model pairs (A, B) in each tail and their tail
#      interaction: on every day
#        d_t = [L_t^low(A) - L_t^low(B)] - [L_t^up(A) - L_t^up(B)],
#      the change between the tails in A's loss advantage over B. FZ0 loss
#      differences are unit-free, so d_t can be pooled across assets and tails.
#   2. Reversal: the better of A and B changes between the tails. A reversal is
#      an intersection of two one-tail hypotheses, so its p-value is the larger
#      of the two tail p-values when the two estimates have opposite signs and
#      1 otherwise (intersection-union test, Berger 1982).
#   3. Scheme contrast: q_t = d_t(Fixed) - d_t(Rolling), and the same contrast
#      for the loss differences in each tail. The two schemes are evaluated on
#      the same days, and every bootstrap replication resamples the same
#      calendar blocks for both.
#   4. Rank shift (appendix table): for each model, the average over assets and
#      quantile levels of (rank in the upper tail - rank in the lower tail).
#   5. Spearman correlation between the ES and EU rankings (models other than
#      the Single Index), averaged over assets and levels (printed and saved).
#   6. Check without SOL, whose evaluation period starts in April 2023: the
#      statistics of GARCH-Skt vs GARCH-FZ for the six assets with a common
#      evaluation period, and the yearly sums behind the figure (printed).
#
# Pooling: every statistic is the equally weighted average over the 28
# asset-level combinations of a per-combination mean over days. SOL enters
# with its shorter evaluation period (from April 2023) and ETH without the two
# days removed as data errors; in a bootstrap replication an asset contributes
# the resampled days on which it has data.
#
# Confidence intervals come from a moving-block bootstrap that resamples blocks
# of calendar days jointly for all assets, tails, levels, models and both
# schemes (block length 20 days; 50 as a check). Bootstrap p-values invert the
# percentile intervals: for an estimate theta with replications theta*_1..B,
#   p = min(1, 2 * min(1 + #{theta*_b <= 0}, 1 + #{theta*_b >= 0}) / (B + 1)),
# with the usual finite-B correction (a p-value is never below 2/(B+1)).
#
# Multiple testing (Holm, 1979) within families of related hypotheses: for
# each scheme, the 8 loss differences in the tails, the 4 interactions, the 4
# reversals and the 11 rank shifts; across the schemes, the 8 contrasts of the
# tail differences and the 4 contrasts of the interactions.
#
# Inputs : forecasts of steps 10/11 (assets/) and 14/15 (assets_roll/)
# Outputs: data/tail_inference.rds, LaTeX/tables/TailInference.tex,
#          LaTeX/tables/TailInferenceRanks.tex, LaTeX/figures/tail_cumdiff.png
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages({ library(data.table); library(sandwich) })
source("R/forecast_eval.R")

SCHEMES <- c(Rolling = file.path(MATLAB_DIR, "assets_roll"), Fixed = file.path(MATLAB_DIR, "assets"))
ASSETS  <- c("BNB", "BTC", "DOGE", "ETH", "SOL", "TRX", "XRP")
PAIRS   <- list(c("GARCH-Skt", "GARCH-FZ"), c("GARCH-Skt", "GARCH-EDF"),
                c("GARCH-Skt", "Hybrid"), c("GARCH-Skt", "RW-250"))
PN      <- sapply(PAIRS, paste, collapse = "|")
B_BOOT  <- 9999
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
# the scheme contrast requires the same evaluation days in both schemes
stopifnot(identical(DAT$Rolling$pos, DAT$Fixed$pos))

# ---- statistics for one set of calendar positions -----------------------------
# row of each asset for every calendar position (NA if the asset has no data)
lk <- lapply(DAT$Rolling$pos, function(p) { m <- rep(NA_integer_, length(cal)); m[p] <- seq_along(p); m })

stats_fast <- function(D, cal_idx) {
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
  c(setNames(rowMeans(do.call(cbind, ranks), na.rm = TRUE), paste0("shift:", MODEL_ORDER)),
    setNames(rowMeans(do.call(cbind, inter), na.rm = TRUE), c(paste0("inter:", PN), paste0("dlow:", PN), paste0("dup:", PN))),
    rho = mean(rhos))
}
# both schemes and their contrasts (Fixed - Rolling) for one set of positions
stats_joint <- function(cal_idx) {
  r <- stats_fast(DAT$Rolling, cal_idx); f <- stats_fast(DAT$Fixed, cal_idx)
  ct <- grepl("^(inter|dlow|dup):", names(r))
  c(setNames(r, paste0("Rolling/", names(r))), setNames(f, paste0("Fixed/", names(f))),
    setNames(f[ct] - r[ct], paste0("Contrast/", names(r)[ct])))
}

mbb_index <- function(n, l) {
  starts <- sample.int(n - l + 1, ceiling(n / l), replace = TRUE)
  as.vector(sapply(starts, function(s) s:(s + l - 1)))[1:n]
}
boot_p <- function(b) { b <- b[is.finite(b)]; min(1, 2 * min(1 + sum(b <= 0), 1 + sum(b >= 0)) / (length(b) + 1)) }

est <- stats_joint(seq_along(cal))
res <- list()
for (l in c(BLOCK, 50)) {
  t0 <- Sys.time()
  boot <- t(replicate(B_BOOT, stats_joint(mbb_index(length(cal), l))))
  cat(sprintf("block %d: %d joint replications in %.1f s\n", l, B_BOOT,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  pval <- apply(boot, 2, boot_p)
  x <- data.table(scheme = sub("/.*", "", names(est)), stat = sub("^[^/]+/", "", names(est)), est = est,
                  lo95 = apply(boot, 2, quantile, 0.025, na.rm = TRUE), hi95 = apply(boot, 2, quantile, 0.975, na.rm = TRUE),
                  p = pval)
  # reversals: intersection-union p-values
  rev <- rbindlist(lapply(c("Rolling", "Fixed"), function(sc) rbindlist(lapply(PN, function(q) {
    lo <- x[scheme == sc & stat == paste0("dlow:", q)]; up <- x[scheme == sc & stat == paste0("dup:", q)]
    data.table(scheme = sc, stat = paste0("rev:", q), est = NA_real_, lo95 = NA_real_, hi95 = NA_real_,
               p = if (sign(lo$est) != sign(up$est)) max(lo$p, up$p) else 1) }))))
  x <- rbind(x, rev)
  x[, family := fcase(grepl("^(dlow|dup):", stat), "tail", grepl("^inter:", stat), "inter",
                      grepl("^rev:", stat), "rev", grepl("^shift:", stat), "shift", default = NA_character_)]
  x[!is.na(family), p_holm := p.adjust(p, method = "holm"), by = .(scheme, family)]
  x[, block := l]
  res[[as.character(l)]] <- x
}
res <- rbindlist(res)

# ---- check: GARCH-Skt vs GARCH-FZ without SOL (six assets, common period) -----
ASSETS_ALL <- ASSETS; ASSETS <- setdiff(ASSETS_ALL, "SOL")       # stats_fast loops over ASSETS
est6 <- stats_joint(seq_along(cal))
boot6 <- t(replicate(1999, stats_joint(mbb_index(length(cal), BLOCK))))
ASSETS <- ASSETS_ALL
nosol <- rbindlist(lapply(c("Rolling", "Fixed", "Contrast"), function(sc) rbindlist(lapply(c("dlow:", "dup:", "inter:"), function(ty) {
  st <- paste0(sc, "/", ty, PN[1]); b <- boot6[, st]
  data.table(scheme = sc, stat = paste0(ty, PN[1]), est = est6[st], lo95 = quantile(b, 0.025), hi95 = quantile(b, 0.975), p = boot_p(b)) }))))
cat("\n=== GARCH-Skt vs GARCH-FZ without SOL (1,999 replications) ===\n"); print(nosol, digits = 3)

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

cat("\n=== Bootstrap results (block 20) ===\n"); print(res[block == BLOCK & stat != "rho"], digits = 3)
cat("\n=== Spearman correlation, Rolling / Fixed: ", sprintf("%.2f / %.2f", est["Rolling/rho"], est["Fixed/rho"]), "\n")
cat("\n=== HAC interaction tests (t < -1.96 : A relatively better in the lower tail) ===\n"); print(hac_sum)
saveRDS(list(boot = res, hac = hac, hac_sum = hac_sum, nosol = nosol), "data/tail_inference.rds")

# ---- figure: cumulative loss differences of GARCH-Skt and GARCH-FZ ------------
# daily average over the asset-level combinations with data on that day
daily_diff <- function(D, p, tl) {
  s <- numeric(length(cal)); n <- numeric(length(cal))
  for (a in ASSETS) for (k in seq_along(TAUS)) {
    L <- D$L[[a]][[k]][[tl]]; pos <- D$pos[[a]]
    s[pos] <- s[pos] + (L[, p[1]] - L[, p[2]]); n[pos] <- n[pos] + 1
  }
  s / n
}
cd <- as.Date(as.character(cal), "%Y%m%d")
CUM <- lapply(names(SCHEMES), function(sc) {
  lo <- daily_diff(DAT[[sc]], PAIRS[[1]], "low"); up <- daily_diff(DAT[[sc]], PAIRS[[1]], "up")
  list(low = cumsum(lo), up = cumsum(up), inter = cumsum(lo - up))
}); names(CUM) <- names(SCHEMES)
png(file.path(FIG_DIR, "tail_cumdiff.png"), width = 1600, height = 660, res = 140)
layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(1, 0.12))
par(mar = c(3.2, 4.2, 2.2, 0.8))
cols <- c(Rolling = "#1f77b4", Fixed = "#d62728")
yl <- range(unlist(lapply(CUM, function(z) c(z$low, z$up))))
plot(NA, xlim = range(cd), ylim = yl, xaxt = "n", xlab = "", ylab = "Cumulative loss difference",
     main = "GARCH-Skt minus GARCH-FZ, by tail", cex.main = 0.95)
axis.Date(1, cd); abline(h = 0, col = "grey60")
for (sc in names(SCHEMES)) { lines(cd, CUM[[sc]]$low, col = cols[sc], lwd = 2); lines(cd, CUM[[sc]]$up, col = cols[sc], lwd = 2, lty = 2) }
plot(NA, xlim = range(cd), ylim = range(unlist(lapply(CUM, `[[`, "inter")), 0), xaxt = "n", xlab = "",
     ylab = "Cumulative interaction", main = "Change between the tails (lower minus upper)", cex.main = 0.95)
axis.Date(1, cd); abline(h = 0, col = "grey60")
for (sc in names(SCHEMES)) lines(cd, CUM[[sc]]$inter, col = cols[sc], lwd = 2)
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = c("Re-estimated every 91 days", "Estimated once", "lower tail", "upper tail"),
       col = c(cols, "black", "black"), lty = c(1, 1, 1, 2), lwd = 2, bty = "n", cex = 0.8, horiz = TRUE)
dev.off()

# check: the same yearly sums without SOL (six assets with a common evaluation period)
yearly <- function(D, assets) {
  s <- list(low = numeric(length(cal)), up = numeric(length(cal))); n <- numeric(length(cal))
  for (a in assets) for (k in seq_along(TAUS)) {
    pos <- D$pos[[a]]; n[pos] <- n[pos] + 1
    for (tl in c("low", "up")) { L <- D$L[[a]][[k]][[tl]]; s[[tl]][pos] <- s[[tl]][pos] + (L[, PAIRS[[1]][1]] - L[, PAIRS[[1]][2]]) }
  }
  ok <- n > 0; yr <- format(cd, "%Y")[ok]
  round(cbind(low = tapply((s$low / n)[ok], yr, sum), up = tapply((s$up / n)[ok], yr, sum)), 1)
}
for (sc in names(SCHEMES)) {
  cat("\n", sc, ": yearly sums of the daily GARCH-Skt minus GARCH-FZ loss differences, 7 assets | 6 assets without SOL\n")
  print(cbind(yearly(DAT[[sc]], ASSETS), yearly(DAT[[sc]], setdiff(ASSETS, "SOL"))))
}

# ---- LaTeX ------------------------------------------------------------------
num <- function(x, d) { z <- sprintf(paste0("%+.", d, "f"), x); sub("^[+-](0\\.0+)$", "\\1", z) }  # no signed zero
pf  <- function(p) if (is.na(p)) "" else if (p < 0.001) "$<$0.001" else sprintf("%.3f", p)
getr <- function(sc, st) res[scheme == sc & block == BLOCK & stat == st]
cell3 <- function(sc, st, d = 3) {              # estimate, CI, p_Holm
  x <- getr(sc, st)
  if (!nrow(x)) return(" & & ")
  if (is.na(x$est)) return(sprintf(" & & %s", pf(x$p_holm)))
  e <- num(x$est, d); if (isTRUE(x$p_holm < 0.05)) e <- paste0("\\textbf{", e, "}")
  sprintf("%s & [%s, %s] & %s", e, num(x$lo95, d), num(x$hi95, d), pf(x$p_holm))
}
cell4 <- function(sc, st, d = 2) {              # estimate, CI, p, p_Holm
  x <- getr(sc, st)
  e <- num(x$est, d); if (isTRUE(x$p_holm < 0.05)) e <- paste0("\\textbf{", e, "}")
  sprintf("%s & [%s, %s] & %s & %s", e, num(x$lo95, d), num(x$hi95, d), pf(x$p), pf(x$p_holm))
}
lab <- function(st) { ty <- sub(":.*", "", st); q <- sub("\\|", " and ", sub("^[a-z]+:", "", st))
  switch(ty, shift = paste("rank shift of", q), inter = paste("interaction of", q),
         dlow = paste("lower-tail difference of", q), dup = paste("upper-tail difference of", q),
         rev = paste("reversal of", q)) }
block_note <- function(fams) {
  bl <- dcast(res[family %in% fams], scheme + stat ~ block, value.var = "p_holm")
  bl <- bl[(get(as.character(BLOCK)) < 0.05) != (`50` < 0.05)]
  sc_lab <- function(sc) if (sc == "Contrast") "scheme contrast" else paste(tolower(sc), "scheme")
  if (!nrow(bl)) return("blocks of 50 days lead to the same conclusions")
  paste0("blocks of 50 days lead to the same conclusions, except for the ",
    paste(sprintf("%s under the %s ($p_{\\mathrm{Holm}} = %.3f$)", sapply(bl$stat, lab), sapply(bl$scheme, sc_lab), bl$`50`),
          collapse = " and the "))
}
PFORM <- paste0("$p = \\min\\{1,\\, 2\\min(1 + \\#\\{\\hat\\theta^*_b \\leq 0\\},\\, 1 + \\#\\{\\hat\\theta^*_b \\geq 0\\})/(B+1)\\}$, ",
  "which inverts the percentile interval and cannot fall below $2/(B+1) = ", sprintf("%.4f", 2 / (B_BOOT + 1)), "$")
PNOTE <- paste0("$p_{\\mathrm{Holm}}$ is the bootstrap p-value ", PFORM, ", adjusted with the method of \\citet{Holm1979}")
POOLNOTE <- paste0("Every statistic is the equally weighted average over the seven assets and four quantile levels of ",
  "a mean over days; SOL enters with its shorter evaluation period and ETH without the two days removed as data errors. ",
  "Confidence intervals are percentile intervals from ", format(B_BOOT, big.mark = ","), " moving-block bootstrap ",
  "replications with blocks of ", BLOCK, " calendar days, resampled jointly across assets, tails, levels, models ",
  "and both schemes")

# main table: loss differences, interactions, reversals and scheme contrasts
rowsB <- unlist(lapply(seq_along(PAIRS), function(i) {
  p <- PAIRS[[i]]; q <- PN[i]
  c(sprintf("\\multicolumn{10}{l}{%s vs %s} \\\\", p[1], p[2]),
    sprintf("\\quad lower tail & %s & %s & %s \\\\", cell3("Rolling", paste0("dlow:", q)), cell3("Fixed", paste0("dlow:", q)), cell3("Contrast", paste0("dlow:", q))),
    sprintf("\\quad upper tail & %s & %s & %s \\\\", cell3("Rolling", paste0("dup:", q)), cell3("Fixed", paste0("dup:", q)), cell3("Contrast", paste0("dup:", q))),
    sprintf("\\quad interaction & %s & %s & %s \\\\", cell3("Rolling", paste0("inter:", q)), cell3("Fixed", paste0("inter:", q)), cell3("Contrast", paste0("inter:", q))),
    sprintf("\\quad reversal & %s & %s & & & \\\\", cell3("Rolling", paste0("rev:", q)), cell3("Fixed", paste0("rev:", q))),
    if (i < length(PAIRS)) "\\hline")
}))
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Inference on the differences between the lower- and upper-tail model comparisons.}",
  "\\label{tab:TailInference}", "\\begin{center}", "\\footnotesize", "\\setlength{\\tabcolsep}{3.5pt}",
  "\\begin{tabular}{lrcrrcrrcr}", "\\toprule",
  "& \\multicolumn{3}{c}{Rolling} & \\multicolumn{3}{c}{Fixed} & \\multicolumn{3}{c}{Fixed $-$ Rolling} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
  "& Estimate & 95\\% CI & $p_{\\mathrm{Holm}}$ & Estimate & 95\\% CI & $p_{\\mathrm{Holm}}$ & Estimate & 95\\% CI & $p_{\\mathrm{Holm}}$ \\\\",
  "\\midrule", rowsB,
  "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: Models re-estimated every 91 days (Rolling), estimated once (Fixed), and the difference ",
    "between the two schemes. Lower and upper tail: average difference between the mean FZ$^{(0)}$ losses of A and B ",
    "(negative: A more accurate). Interaction: the mean of $d_t(A,B) = [L^{low}_t(A) - L^{low}_t(B)] - [L^{up}_t(A) - L^{up}_t(B)]$ ",
    "(negative: A's advantage over B is larger in the lower tail). Reversal: test that the more accurate of A and B ",
    "differs between the tails, which rejects only if both tail differences are significant with opposite signs ",
    "(intersection--union test); its p-value is the larger of the two tail p-values, or 1 if the two estimates have ",
    "the same sign. Fixed $-$ Rolling: contrast between the schemes, e.g.\\ the mean of ",
    "$q_t = d^{\\mathrm{Fixed}}_t - d^{\\mathrm{Rolling}}_t$ for the interaction. ", POOLNOTE, "; ",
    block_note(c("tail", "inter", "rev")), ". ", PNOTE, " within families of related tests: for each scheme, ",
    "the 8 tail differences, the 4 interactions and the 4 reversals; for the contrasts, the 8 tail differences and ",
    "the 4 interactions. Bold: $p_{\\mathrm{Holm}} < 0.05$.}"),
  "\\end{table}")
writeLines(L, file.path(TEX_DIR, "TailInference.tex"))

# appendix table: rank shifts
rowsA <- unlist(lapply(MODEL_ORDER, function(m) {
  r <- sprintf("%-12s & %s & %s \\\\", m, cell4("Rolling", paste0("shift:", m)), cell4("Fixed", paste0("shift:", m)))
  if (m %in% c("RW-500", "GARCH-FZ", "Hybrid")) c(r, "\\hline") else r }))
L <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Average rank shift between the tails (EU rank minus ES rank), seven cryptocurrencies.}",
  "\\label{tab:TailInferenceRanks}", "\\begin{center}", "\\small", "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{lrcrrrcrr}", "\\toprule",
  "& \\multicolumn{4}{c}{Rolling} & \\multicolumn{4}{c}{Fixed} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
  "& Estimate & 95\\% CI & $p$ & $p_{\\mathrm{Holm}}$ & Estimate & 95\\% CI & $p$ & $p_{\\mathrm{Holm}}$ \\\\", "\\midrule",
  rowsA, "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
  paste0("\\textit{Notes: Average over assets and levels of a model's rank by mean FZ$^{(0)}$ loss in the upper tail ",
    "minus its rank in the lower tail, for the models re-estimated every 91 days (Rolling) and estimated once (Fixed); ",
    "a positive value indicates a model that is relatively worse in the upper tail. ", POOLNOTE, "; ",
    block_note("shift"), ". $p$ is the bootstrap p-value ", PFORM, ", and $p_{\\mathrm{Holm}}$ adjusts it with the method of ",
    "\\citet{Holm1979} within the 11 rank shifts of each scheme. Bold: $p_{\\mathrm{Holm}} < 0.05$.}"),
  "\\end{table}")
writeLines(L, file.path(TEX_DIR, "TailInferenceRanks.tex"))
cat("\nTables written:", file.path(TEX_DIR, "TailInference.tex"), file.path(TEX_DIR, "TailInferenceRanks.tex"), "\n")
