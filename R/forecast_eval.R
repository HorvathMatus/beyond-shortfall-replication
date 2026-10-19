# ============================================================
# Beyond Shortfall — shared forecast-evaluation functions
#
# Used by 04_mcs_tables.R, 08_forecast_evaluation.R (Bitcoin),
# 11_multiasset_evaluation.R (models estimated once),
# 15_rolling_evaluation.R (models re-estimated on a rolling window) and
# 16_tail_inference.R.
# Requires: MCS, rugarch, esback
# ============================================================

RAW_COLS <- c("RW-125","RW-250","RW-500","GARCH-N","GARCH-Skt",
              "GARCH-EDF","FZ-2F","FZ-1F","GARCH-FZ","Hybrid")      # MATLAB column order
MODEL_ORDER <- c("RW-125","RW-250","RW-500","GARCH-N","GARCH-Skt",
                 "GARCH-EDF","GARCH-FZ","FZ-1F","FZ-2F","Hybrid","Single Index")
GAS    <- c("FZ-1F", "FZ-2F", "Hybrid")
TAUS   <- c(0.01, 0.025, 0.05, 0.1)
SUFFIX <- c("001","0025","005","01")

# ---- scoring functions --------------------------------------------------------
# Both are members of the Fissler-Ziegel class, strictly consistent for
# (VaR, ES) when ES < 0:  S = G2(e) (e - v + 1{r<=v}(v - r)/tau) - calG2(e).
#   FZ0     G2(e) = -1/e,          calG2(e) = -log(-e)      homogeneous of degree 0
#           (Patton, Ziegel and Chen, 2019)
#   FZ(1/2) G2(e) = (-e)^(-1/2),   calG2(e) = -2 (-e)^(1/2) homogeneous of degree 1/2
#           (Nolde and Ziegel, 2017); puts more weight on high-risk periods.
# The AL log score of Taylor (2019) is not used: AL = FZ0 + 1 - log(1 - tau) - r/e,
# so it differs from FZ0 only by a term with zero conditional mean when returns
# have zero conditional mean, and is not strictly consistent otherwise.
fz_loss <- function(r, var, es, tau) {
  hit <- as.integer(r <= var)
  (hit * (r - var)) / (tau * es) + var / es + log(-es) - 1
}
fzh_loss <- function(r, var, es, tau, b = 0.5) {
  hit <- as.integer(r <= var)
  (-es)^b * (1 / b - 1) + (-es)^(b - 1) * (hit * (var - r) / tau - var)
}

# ---- MCS p-values along the full elimination path (as in step 04) -----------
mcs_pvalues <- function(L, B = 5000, seed = 123) {
  set.seed(seed)
  out <- capture.output(m <- MCS::MCSprocedure(Loss = L, alpha = 0.9999, B = B,
                                                statistic = "Tmax"))
  p <- setNames(rep(NA_real_, ncol(L)), colnames(L))
  el <- grep("is eliminated", out); pl <- grep("current MCS p-Value", out)
  for (i in el) {
    nm <- sub("^Model (.*) is eliminated.*$", "\\1", out[i])
    p[nm] <- as.numeric(sub(".*is ", "", out[pl[pl > i][1]]))
  }
  p[is.na(p)] <- 1
  ord <- order(p); p[ord] <- cummax(p[ord])
  p
}

# ---- real-time usability of the forecasts -----------------------------------
# The scoring functions require ES <= VaR < 0. A forecast for day t that is
# non-finite or violates this condition is replaced by the model's most recent
# usable forecast, a decision that uses only information available at t - 1.
# Leading unusable forecasts (no earlier usable one) are set to NA, which
# excludes the model from that comparison. Returns the cleaned matrices and the
# number of replaced days per model (attribute "replaced").
realtime_guard <- function(v, e) {
  rep_n <- setNames(integer(ncol(v)), colnames(v))
  for (j in seq_len(ncol(v))) {
    ok <- is.finite(v[, j]) & is.finite(e[, j]) & v[, j] < 0 & e[, j] <= v[, j]
    if (all(ok)) next
    last <- 0L
    for (t in seq_len(nrow(v))) {
      if (ok[t]) { last <- t
      } else if (last > 0L) { v[t, j] <- v[last, j]; e[t, j] <- e[last, j]; rep_n[j] <- rep_n[j] + 1L
      } else { v[t, j] <- NA; e[t, j] <- NA }
    }
  }
  structure(list(v = v, e = e), replaced = rep_n)
}

# ---- forecasts of one tail and level from one folder --------------------------
# dir: folder with the MATLAB files VaR_<tail>_<suf>.txt, ES_<tail>_<suf>.txt
# (P x 10) and the Single-Index files VaR_SI_*, ES_SI_*. Upper-tail forecasts
# are in the negated space, so ES <= VaR < 0 applies to both tails.
read_fc <- function(dir, tail, suf) {
  v <- as.matrix(read.table(file.path(dir, paste0("VaR_", tail, "_", suf, ".txt")), sep = ","))
  e <- as.matrix(read.table(file.path(dir, paste0("ES_",  tail, "_", suf, ".txt")), sep = ","))
  colnames(v) <- colnames(e) <- RAW_COLS
  v <- cbind(v, "Single Index" = scan(file.path(dir, paste0("VaR_SI_", tail, "_", suf, ".txt")), quiet = TRUE))
  e <- cbind(e, "Single Index" = scan(file.path(dir, paste0("ES_SI_",  tail, "_", suf, ".txt")), quiet = TRUE))
  realtime_guard(v[, MODEL_ORDER], e[, MODEL_ORDER])
}

# ---- losses, MCS and backtests for every tail, level and model of one asset --
# Models with missing forecasts are excluded from the MCS and backtests.
eval_asset <- function(root, a) {
  dir <- file.path(root, a)
  r_oos <- scan(file.path(dir, "oos_ret.txt"), quiet = TRUE)
  rows <- list()
  for (tl in c("low", "up")) {
    y <- if (tl == "low") r_oos else -r_oos
    for (k in seq_along(TAUS)) {
      tau <- TAUS[k]; fc <- read_fc(dir, tl, SUFFIX[k]); nrep <- attr(fc, "replaced")
      fz  <- sapply(MODEL_ORDER, function(m) fz_loss(y, fc$v[, m], fc$e[, m], tau))
      fzh <- sapply(MODEL_ORDER, function(m) fzh_loss(y, fc$v[, m], fc$e[, m], tau))
      ok <- colSums(!is.finite(fz)) == 0 & colSums(!is.finite(fzh)) == 0
      p_fz <- p_fzh <- setNames(rep(NA_real_, length(MODEL_ORDER)), MODEL_ORDER)
      p_fz[ok]  <- mcs_pvalues(fz[, ok, drop = FALSE])
      p_fzh[ok] <- mcs_pvalues(fzh[, ok, drop = FALSE])
      mfz <- ifelse(ok, colMeans(fz), NA); mfzh <- ifelse(ok, colMeans(fzh), NA)
      for (m in MODEL_ORDER) {
        v <- fc$v[, m]; e <- fc$e[, m]
        bt <- if (ok[m]) {
          vt  <- tryCatch(rugarch::VaRTest(alpha = tau, actual = y, VaR = v, conf.level = 0.95),
                          error = function(err) list(uc.LRp = NA, cc.LRp = NA))
          set.seed(123)   # the ESR test uses random starting values (esreg); fixed for reproducibility
          esr <- tryCatch(esback::esr_backtest(y, v, e, alpha = tau, version = 1)$pvalue_twosided_asymptotic,
                          error = function(err) NA_real_)
          c(hit = mean(y <= v), uc = vt$uc.LRp, cc = vt$cc.LRp, esr = esr)
        } else c(hit = NA, uc = NA, cc = NA, esr = NA)
        rows[[length(rows) + 1]] <- data.frame(asset = a, P = length(y), tail = tl, tau = tau,
          model = m, fz = mfz[m], fzh = mfzh[m], p_fz = p_fz[m], p_fzh = p_fzh[m],
          replaced = nrep[m], hit = bt["hit"], uc = bt["uc"], cc = bt["cc"], esr = bt["esr"])
      }
    }
  }
  do.call(rbind, rows)
}

# ---- Single-Index forecasts for one block (as in step 03) --------------------
# y: full return series of one tail (negated for the upper tail);
# s..e: forecast dates; the model is estimated on the R observations before s.
# With s = R + 1 and e = length(y) this is exactly step 03.
si_block <- function(y, s, e, R, alpha) {
  tr <- (s - R + 2):(s - 1); te <- s:e
  X_tr <- si_design(y[tr - 1], y[tr - 2]); X_te <- si_design(y[te - 1], y[te - 2])
  cal <- si_calibrate_tau(y[tr], X_tr, alpha = alpha)
  v <- si_predict(cal$fit, X_te)
  list(v = v, e = si_es_from_expectile(v, cal$tau, alpha, mean(y[(s - R):(s - 1)])))
}
