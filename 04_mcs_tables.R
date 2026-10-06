# ============================================================
# Beyond Shortfall — 04: Model Confidence Set for ES (Shortfall) and EU (Upside)
#
# Pipeline: 01_download_data.R  ->  02 (MATLAB: 02_BTC_forecast_tails.m,
# see 02_run_matlab_forecasts.R)  ->  03_single_index.R  ->  this script.
#
# Reads the step-02 MATLAB outputs (VaR_<tail>_<suf>.txt, ES_<tail>_<suf>.txt,
# P x 10, columns: RW-125, RW-250, RW-500, GARCH-N, GARCH-Skt, GARCH-EDF,
# FZ2F, FZ1F, GARCH-FZ, Hybrid) plus the step-03 Single-Index forecasts
# (VaR_SI_<tail>_<suf>.txt, ES_SI_...), computes the FZ^(0) loss and the
# FZ^(1/2) score (R/forecast_eval.R) and runs the MCS per tau, separately for:
#   Shortfall (lower tail, original returns)
#   Upside    (upper tail; forecasts and returns in NEGATED space, so the
#              same lower-tail loss functions apply unchanged)
#
# MCS (Hansen, Lunde and Nason, 2011): T_max statistic, B = 5000 block
# bootstrap replications, fixed seed. We report the MCS p-value of every
# model, i.e. the running maximum of the test p-values along the full
# elimination path, and the standard sets
#     90% MCS = {models with p_MCS >= 0.10},  75% MCS = {p_MCS >= 0.25}.
#
# NB on the R package `MCS`: MCSprocedure(alpha = a) keeps a model whenever
# the test p-value at ITS OWN elimination step is >= a, i.e. it uses the raw
# step p-values rather than their running maximum. For a <= 0.25 this
# coincides with the Hansen et al. definition in our data, but for large a
# (e.g. 0.75, 0.90 as in the earlier letter) it does not, and those values
# also correspond to 25% and 10% sets rather than 75% and 90%. We therefore
# run the procedure to the end (alpha close to one) and compute the MCS
# p-values ourselves from the printed elimination path.
# ============================================================
# SCHEME: "fixed" (default; models estimated once, step 02; tables with the
# suffix _fixed, reported in the appendix) or "rolling" (models re-estimated
# every 91 days, steps 14-15; the main tables of Section 5, produced by step 21).
SCHEME <- if (exists("SCHEME")) SCHEME else "fixed"
rm(list = setdiff(ls(), "SCHEME"))
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
if (!requireNamespace("MCS", quietly = TRUE)) install.packages("MCS")
library(MCS)
source("R/forecast_eval.R")   # fz_loss, fzh_loss, mcs_pvalues, read_fc (real-time guard)

FC_DIR  <- if (SCHEME == "rolling") file.path(MATLAB_DIR, "assets_roll", "BTC") else MATLAB_DIR
SUF_OUT <- if (SCHEME == "rolling") "" else "_fixed"
SCHEME_NOTE <- if (SCHEME == "rolling") "Models re-estimated every 91 days on a rolling window of 1{,}095 days. " else "Models estimated once on the first 1{,}095 days and filtered forward. "
cat("Scheme:", SCHEME, "\n")
TAILS  <- c(Shortfall = "low", Upside = "up")
P90 <- 0.10; P75 <- 0.25                     # MCS p-value thresholds

# ---- data ------------------------------------------------------------------
r_oos <- scan(file.path(FC_DIR, "oos_ret.txt"), quiet = TRUE)  # percent, original sign
cat("OOS observations:", length(r_oos), "\n")

# MCS p-values along the full elimination path (mcs_pvalues, R/forecast_eval.R):
# the procedure is run until a single model remains; each eliminated model
# receives the MCS p-value printed at its elimination (the running maximum of
# the step p-values), and the last surviving model receives 1.

# ---- main loop -------------------------------------------------------------
res <- list()
for (tl in names(TAILS)) {
  tail_id <- TAILS[[tl]]
  y <- if (tail_id == "low") r_oos else -r_oos   # negated space for Upside
  res[[tl]] <- vector("list", length(TAUS))

  for (k in seq_along(TAUS)) {
    tau <- TAUS[k]; suf <- SUFFIX[k]
    # step 02 MATLAB forecasts + step 03 Single-Index forecasts; forecasts
    # violating ES <= VaR < 0 are replaced in real time (R/forecast_eval.R)
    fc <- read_fc(FC_DIR, tail_id, suf)
    v <- fc$v; e <- fc$e
    stopifnot(nrow(v) == length(y))
    nrep <- attr(fc, "replaced")
    if (any(nrep > 0)) cat(sprintf("  [%s tau=%.3f] forecasts replaced in real time: %s\n", tl, tau,
                                   paste(names(nrep[nrep > 0]), nrep[nrep > 0], collapse = ", ")))

    fz_mat  <- sapply(MODEL_ORDER, function(m) fz_loss(y, v[, m], e[, m], tau))
    fzh_mat <- sapply(MODEL_ORDER, function(m) fzh_loss(y, v[, m], e[, m], tau))

    bad <- colSums(!is.finite(fz_mat))
    if (any(bad > 0)) {
      cat(sprintf("  [%s tau=%.3f] non-finite losses: %s\n", tl, tau,
          paste(names(bad[bad>0]), bad[bad>0], collapse=", ")))
    }

    cat(sprintf("[%s] tau=%.3f: running MCS (FZ0, FZ1/2) ...\n", tl, tau))
    res[[tl]][[k]] <- list(tau = tau, fz = fz_mat, fzh = fzh_mat, replaced = nrep,
                           p_fz = mcs_pvalues(fz_mat), p_fzh = mcs_pvalues(fzh_mat))
  }
}

# ---- table builder ---------------------------------------------------------
build_table <- function(res_tail, loss = "fz") {
  pk <- paste0("p_", loss)
  means_list <- lapply(res_tail, function(x) colMeans(x[[loss]]))
  ranks_list <- lapply(means_list, function(m) rank(m, ties.method = "first"))
  avg_rank   <- rowMeans(do.call(cbind, ranks_list))
  rows <- lapply(seq_along(TAUS), function(k) {
    data.frame(Mean = round(means_list[[k]], 3),
               Pmcs = round(res_tail[[k]][[pk]][MODEL_ORDER], 3),
               Rank = as.integer(ranks_list[[k]]),
               stringsAsFactors = FALSE, row.names = MODEL_ORDER)
  })
  df <- do.call(cbind, rows)
  colnames(df) <- paste0(rep(c("Mean","Pmcs","Rank"), length(TAUS)), "_", rep(TAUS, each = 3))
  df$AvgRank <- round(avg_rank, 2)
  df
}

print_table <- function(df, caption) {
  cat(sprintf("\n%s\n%s\n", caption, strrep("-", 104)))
  cat(sprintf("%-12s", "Model"))
  for (tau in TAUS) cat(sprintf("  %8s  pMCS  R |", paste0("tau=", tau)))
  cat("  AvgRank\n")
  for (m in MODEL_ORDER) {
    cat(sprintf("%-12s", m))
    for (tau in TAUS)
      cat(sprintf("  %8.3f %5.2f %2d |", df[m, paste0("Mean_", tau)],
          df[m, paste0("Pmcs_", tau)], df[m, paste0("Rank_", tau)]))
    cat(sprintf("  %7.2f\n", df[m, "AvgRank"]))
  }
}

# ---- output ----------------------------------------------------------------
tables <- list()
for (tl in names(TAILS)) {
  for (loss in c("fz","fzh")) {
    key <- paste0(tl, "_", loss)
    tables[[key]] <- build_table(res[[tl]], loss)
    lname <- if (loss=="fz") "FZ0 loss" else "FZ1/2 score"
    print_table(tables[[key]], sprintf("Table [%s / %s]: mean loss, MCS p-value, ranks", tl, lname))
  }
}
save(res, tables, file = paste0("data/mcs_results_beyond", SUF_OUT, ".RData"))

# ---- paper-style tables written straight into the LaTeX project -----------
# Each tau block reports the mean loss, the MCS p-value and the rank. MCS
# p-values of at least 0.10 (members of the 90% MCS) are printed in bold.
pmcs_fmt <- function(p) {
  s <- if (p >= 0.995) "1.00" else sprintf("%.2f", p)
  if (p >= P90) paste0("\\textbf{", s, "}") else s
}
paper_latex <- function(df, label, caption, notes) {
  L <- c("\\begin{table}[ht]", "", "\\centering", "",
    paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"), "",
    "\\begin{center}", "\\begin{tabular}{lcccccccccccccc}", "\\toprule",
    paste0("\\multirow{2}{*}{Model} & \\multicolumn{3}{c}{$\\tau = 0.01$} & ",
           "\\multicolumn{3}{c}{$\\tau = 0.025$} & \\multicolumn{3}{c}{$\\tau = 0.05$} & ",
           "\\multicolumn{3}{c}{$\\tau = 0.1$} & \\multirow{2}{*}{Avg rank} \\\\"),
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10} \\cmidrule(lr){11-13}",
    paste0(" & Mean & $p_{\\mathrm{MCS}}$ & R & Mean & $p_{\\mathrm{MCS}}$ & R & Mean & ",
           "$p_{\\mathrm{MCS}}$ & R & Mean & $p_{\\mathrm{MCS}}$ & R \\\\"),
    "\\midrule")
  for (m in MODEL_ORDER) {
    L <- c(L, sprintf("%-12s & %.3f & %s & %-2d & %.3f & %s & %-2d & %.3f & %s & %-2d & %.3f & %s & %-2d & %.2f \\\\",
      m,
      df[m,"Mean_0.01"],  pmcs_fmt(df[m,"Pmcs_0.01"]),  df[m,"Rank_0.01"],
      df[m,"Mean_0.025"], pmcs_fmt(df[m,"Pmcs_0.025"]), df[m,"Rank_0.025"],
      df[m,"Mean_0.05"],  pmcs_fmt(df[m,"Pmcs_0.05"]),  df[m,"Rank_0.05"],
      df[m,"Mean_0.1"],   pmcs_fmt(df[m,"Pmcs_0.1"]),   df[m,"Rank_0.1"], df[m, "AvgRank"]))
    if (m %in% c("RW-500","GARCH-FZ","Hybrid")) L <- c(L, "\\hline")
  }
  c(L, "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
    paste0("\\textit{Notes: ", notes, "}"), "", "\\end{table}")
}

if (dir.exists(TEX_DIR)) {
  meas <- list(fz = "FZ$^{(0)}$ loss", fzh = "FZ$^{(1/2)}$ score")
  for (tl in names(TAILS)) for (loss in c("fz","fzh")) {
    q   <- if (tl == "Shortfall") "Expected Shortfall" else "Expected Upside"
    tag <- if (tl == "Shortfall") "ES" else "EU"
    writeLines(paper_latex(tables[[paste0(tl,"_",loss)]],
      label   = sprintf("tab:%sLoss_%s%s", toupper(loss), tag, SUF_OUT),
      caption = sprintf("Out-of-sample %s and model confidence set: %s of Bitcoin%s.", meas[[loss]], q,
                        if (SCHEME == "rolling") "" else ", models estimated once"),
      notes   = paste0(SCHEME_NOTE, "Mean is the average ", meas[[loss]], " of the ", q, " forecasts of each model over ",
        "the out-of-sample period (October 1, 2020 to June 30, 2026; 2,099 days)",
        if (tl == "Upside") ", computed on the negated return series" else "", ". ",
        "$p_{\\mathrm{MCS}}$ is the MCS p-value of \\citet{Hansen2011} ($T_{\\max}$ statistic, ",
        "5{,}000 bootstrap replications); models with $p_{\\mathrm{MCS}} \\geq 0.10$ belong to the ",
        "90\\% model confidence set and are shown in bold, and models with $p_{\\mathrm{MCS}} \\geq 0.25$ ",
        "also belong to the 75\\% set. R is the rank by mean loss and Avg rank its average over the four ",
        "quantile levels.")),
      file.path(TEX_DIR, sprintf("%sLoss_%s%s.tex", toupper(loss), tag, SUF_OUT)))
  }
  cat("\npaper-style MCS tables written to", TEX_DIR, "\n")
}
cat("results saved to", paste0("data/mcs_results_beyond", SUF_OUT, ".RData"), "\n")
