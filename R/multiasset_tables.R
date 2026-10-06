# ============================================================
# Beyond Shortfall — cross-asset summary tables (by model and by asset)
#
# Used by 11_multiasset_evaluation.R (models estimated once; tables with the
# suffix "_fixed", reported in the appendix) and 15_rolling_evaluation.R
# (models re-estimated every 91 days; the main tables of Section 5).
# ev: output of eval_asset() for all assets, with rank_fz and rank_fzh.
# ============================================================
write_multiasset_tables <- function(ev, tex_dir, suffix = "", scheme_note = "") {
  ASSETS <- unique(ev$asset)
  # ---- 3. summary by model ----------------------------------------------------
  SIG <- 0.05
  sm <- ev[, .(mcs90 = sum(p_fz >= 0.10, na.rm = TRUE),
               mcs90_fzh = sum(p_fzh >= 0.10, na.rm = TRUE),
               avgrank = mean(rank_fz, na.rm = TRUE),
               avgrank_fzh = mean(rank_fzh, na.rm = TRUE),
               rej = sum(uc < SIG, na.rm = TRUE) + sum(cc < SIG, na.rm = TRUE) + sum(esr < SIG, na.rm = TRUE),
               ncases = sum(!is.na(fz))),
           by = .(tail, model)]
  cat("\n=== Summary by model (FZ loss) ===\n")
  print(dcast(sm, model ~ tail, value.var = c("mcs90", "avgrank", "rej")))
  cat("\nFZ1/2 robustness (cases in 90% MCS, average rank):\n")
  print(dcast(sm, model ~ tail, value.var = c("mcs90_fzh", "avgrank_fzh")))

  ncase <- length(ASSETS) * length(TAUS)
  rows <- character(0)
  for (m in MODEL_ORDER) {
    a <- sm[tail == "low" & model == m]; b <- sm[tail == "up" & model == m]
    rows <- c(rows, sprintf("%-12s & %d & %.2f & %d & %d & %.2f & %d & %+.2f \\\\", m,
      a$mcs90, a$avgrank, a$rej, b$mcs90, b$avgrank, b$rej, b$avgrank - a$avgrank))
    if (m %in% c("RW-500", "GARCH-FZ", "Hybrid")) rows <- c(rows, "\\hline")
  }
  wrap <- function(rows, label, caption, notes, hdr, align) {
    c("\\begin{table}[ht]", "\\centering", paste0("\\caption{", caption, "}"),
      paste0("\\label{", label, "}"), "\\begin{center}",
      paste0("\\begin{tabular}{", align, "}"), "\\toprule", hdr, "\\midrule", rows,
      "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
      paste0("\\textit{Notes: ", notes, "}"), "\\end{table}")
  }
  writeLines(wrap(rows, paste0("tab:MultiAsset", suffix),
    paste0("Forecast accuracy and calibration across the seven risky cryptocurrencies", scheme_note, "."),
    paste0("Results pooled over BNB, BTC, DOGE, ETH, SOL, TRX and XRP and the four quantile levels ",
      "$\\tau \\in \\{0.01, 0.025, 0.05, 0.1\\}$, that is ", ncase, " cases per tail. ",
      "MCS counts the cases in which the model belongs to the 90\\% model confidence set (FZ loss, ",
      "$p_{\\mathrm{MCS}} \\geq 0.10$); Avg rank is the average rank by mean FZ loss; Rej counts the ",
      "rejections at the 5\\% level of the UC, CC and ESR backtests (maximum ", 3 * ncase, "). ",
      "$\\Delta$ is the EU average rank minus the ES average rank. The upper tail is evaluated on the ",
      "negated return series."),
    paste0("& \\multicolumn{3}{c}{Lower tail (ES)} & \\multicolumn{3}{c}{Upper tail (EU)} & \\\\\n",
           "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n",
           "Model & MCS & Avg rank & Rej & MCS & Avg rank & Rej & $\\Delta$ \\\\"),
    "lrrrrrrr"), file.path(tex_dir, paste0("MultiAsset", suffix, ".tex")))

  # ---- 4. summary by asset -----------------------------------------------------
  best_model <- function(dt) {
    r <- dt[, .(ar = mean(rank_fz, na.rm = TRUE)), by = model]
    r[which.min(ar), model]
  }
  rho_tau <- function(a, tt, drop_si = TRUE) {      # tt, not tau: tau is a column of ev
    lo <- ev[asset == a & tail == "low" & abs(tau - tt) < 1e-9, .(model, rank_fz)]
    up <- ev[asset == a & tail == "up"  & abs(tau - tt) < 1e-9, .(model, rank_fz)]
    x  <- merge(lo, up, by = "model", suffixes = c("_lo", "_up"))
    if (drop_si) x <- x[model != "Single Index"]
    x <- x[complete.cases(x)]
    cor(x$rank_fz_lo, x$rank_fz_up, method = "spearman")
  }
  ba <- rbindlist(lapply(ASSETS, function(a) {
    rh <- sapply(TAUS, function(t) rho_tau(a, t))
    data.table(asset = a, P = ev[asset == a, P[1]],
      best_es = best_model(ev[asset == a & tail == "low"]),
      best_eu = best_model(ev[asset == a & tail == "up"]),
      rho01 = rh[1], rho05 = rh[3], rho_mean = mean(rh),
      gas_es = ev[asset == a & tail == "low" & model %in% GAS, sum(p_fz < 0.10, na.rm = TRUE)],
      gas_eu = ev[asset == a & tail == "up"  & model %in% GAS, sum(p_fz < 0.10, na.rm = TRUE)],
      mcs_es = ev[asset == a & tail == "low", sum(p_fz >= 0.10, na.rm = TRUE)] / length(TAUS),
      mcs_eu = ev[asset == a & tail == "up",  sum(p_fz >= 0.10, na.rm = TRUE)] / length(TAUS))
  }))
  cat("\n=== Summary by asset ===\n"); print(ba, digits = 3)

  rows <- sprintf("%-5s & %d & %s & %s & %.2f & %.2f & %.2f & %.1f & %.1f & %d & %d \\\\",
    ba$asset, ba$P, ba$best_es, ba$best_eu, ba$rho01, ba$rho05, ba$rho_mean,
    ba$mcs_es, ba$mcs_eu, ba$gas_es, ba$gas_eu)
  tot <- sprintf("All   & & & & %.2f & %.2f & %.2f & %.1f & %.1f & %d & %d \\\\",
    mean(ba$rho01), mean(ba$rho05), mean(ba$rho_mean), mean(ba$mcs_es), mean(ba$mcs_eu),
    sum(ba$gas_es), sum(ba$gas_eu))
  writeLines(wrap(c(rows, "\\hline", tot), paste0("tab:MultiAssetByAsset", suffix),
    paste0("Agreement between the lower- and upper-tail model comparisons, by asset", scheme_note, "."),
    paste0("$P$ is the number of out-of-sample days. Best ES and Best EU are the models with the lowest ",
      "average rank by FZ loss over the four quantile levels. $\\rho$ is the Spearman correlation between ",
      "the ES and EU rankings of the ten models other than the Single Index model, at $\\tau = 0.01$, at ",
      "$\\tau = 0.05$ and averaged over the four levels. $|M^*|$ is the average size of the 90\\% model ",
      "confidence set. GAS out counts the cases (of 12) in which one of the three score-driven models ",
      "(FZ-1F, FZ-2F, Hybrid) is excluded from the 90\\% set. The last row reports averages of $\\rho$ ",
      "and $|M^*|$ and totals of the exclusions."),
    paste0("& & & & \\multicolumn{3}{c}{$\\rho$ (ES vs EU)} & \\multicolumn{2}{c}{$|M^*|$} & ",
           "\\multicolumn{2}{c}{GAS out} \\\\\n",
           "\\cmidrule(lr){5-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}\n",
           "Asset & $P$ & Best ES & Best EU & 0.01 & 0.05 & Mean & ES & EU & ES & EU \\\\"),
    "lrllrrrrrrr"), file.path(tex_dir, paste0("MultiAssetByAsset", suffix, ".tex")))
  cat(sprintf("\nTables written: MultiAsset%s.tex, MultiAssetByAsset%s.tex\n", suffix, suffix))

  invisible(list(by_model = sm, by_asset = ba))
}
