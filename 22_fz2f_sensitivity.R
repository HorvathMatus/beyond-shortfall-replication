# ============================================================
# Beyond Shortfall — 22: sensitivity of the cross-asset rankings to FZ-2F
#
# The two-factor GAS model (FZ-2F) produces forecasts that violate
# ES <= VaR < 0 on up to 11% of the days of an asset-tail-level combination;
# they are replaced in real time by the model's last usable forecast
# (realtime_guard in R/forecast_eval.R). This step checks whether the main
# cross-asset results (models re-estimated every 91 days, step 15) depend on
# these replacements:
#   Baseline   all models, all days (as in Tables MultiAsset, MultiAssetByAsset)
#   Without    FZ-2F removed from the comparison
#   Clean days all models, evaluated only on the days on which the FZ-2F
#              forecast did not have to be replaced
# For each variant: average rank by mean FZ0 loss in each tail (pooled over
# the seven assets and four levels), the number of assets whose best model
# (average rank over the levels, ties broken in the order of MODEL_ORDER, as in
# Table MultiAssetByAsset) differs between the tails, and the average Spearman
# correlation between the lower- and upper-tail rankings (without the Single
# Index model, as in Table MultiAssetByAsset).
#
# Inputs : MATLAB_DIR/assets_roll/<asset>/ (steps 14 and 15)
# Outputs: data/fz2f_sensitivity.rds (printed summary)
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
source("R/forecast_eval.R")

ROOT   <- file.path(MATLAB_DIR, "assets_roll")
ASSETS <- c("BNB", "BTC", "DOGE", "ETH", "SOL", "TRX", "XRP")

# usable raw FZ-2F forecasts (before the real-time replacement)
fz2f_ok <- function(dir, tl, suf) {
  v <- as.matrix(read.table(file.path(dir, paste0("VaR_", tl, "_", suf, ".txt")), sep = ","))
  e <- as.matrix(read.table(file.path(dir, paste0("ES_",  tl, "_", suf, ".txt")), sep = ","))
  j <- match("FZ-2F", RAW_COLS)
  is.finite(v[, j]) & is.finite(e[, j]) & v[, j] < 0 & e[, j] <= v[, j]
}

# mean FZ0 losses by asset, tail, level and variant
ML <- list()
for (a in ASSETS) {
  dir <- file.path(ROOT, a); r <- scan(file.path(dir, "oos_ret.txt"), quiet = TRUE)
  for (tl in c("low", "up")) for (k in seq_along(TAUS)) {
    y <- if (tl == "low") r else -r
    fc <- read_fc(dir, tl, SUFFIX[k])
    L  <- sapply(MODEL_ORDER, function(m) fz_loss(y, fc$v[, m], fc$e[, m], TAUS[k]))
    clean <- fz2f_ok(dir, tl, SUFFIX[k]) & rowSums(!is.finite(L)) == 0
    ML[[length(ML) + 1]] <- data.frame(asset = a, tail = tl, tau = TAUS[k], model = MODEL_ORDER,
      base = colMeans(L), clean = colMeans(L[clean, , drop = FALSE]), share_clean = mean(clean))
  }
}
ML <- do.call(rbind, ML)

rank_by <- function(x, col, drop = NULL) {
  x <- x[!x$model %in% drop, ]
  x$rank <- ave(x[[col]], x$asset, x$tail, x$tau, FUN = rank)
  x
}
summ_variant <- function(x) {
  avg  <- tapply(x$rank, list(x$model, x$tail), mean)
  byas <- aggregate(rank ~ asset + tail + model, x, mean)
  byas <- byas[order(match(byas$model, MODEL_ORDER)), ]
  best <- sapply(ASSETS, function(a) sapply(c("low", "up"), function(tl) {
    b <- byas[byas$asset == a & byas$tail == tl, ]; b$model[which.min(b$rank)] }))
  rho <- mapply(function(a, t) {
    lo <- x[x$asset == a & x$tau == t & x$tail == "low" & x$model != "Single Index", ]
    up <- x[x$asset == a & x$tau == t & x$tail == "up" & x$model != "Single Index", ]
    cor(lo$rank, up$rank[match(lo$model, up$model)], method = "spearman") },
    rep(ASSETS, each = length(TAUS)), rep(TAUS, length(ASSETS)))
  list(avg = avg, best = best, n_differ = sum(best["low", ] != best["up", ]), rho = mean(rho))
}
V <- list(Baseline  = summ_variant(rank_by(ML, "base")),
          Without   = summ_variant(rank_by(ML, "base", drop = "FZ-2F")),
          CleanDays = summ_variant(rank_by(ML, "clean")))

cat(sprintf("Share of days with a usable FZ-2F forecast: min %.3f, mean %.3f\n\n",
            min(ML$share_clean), mean(ML$share_clean)))
for (nm in names(V)) {
  cat("=== ", nm, " ===\n", sep = "")
  print(round(V[[nm]]$avg[intersect(MODEL_ORDER, rownames(V[[nm]]$avg)), c("low", "up")], 2))
  print(V[[nm]]$best)
  cat(sprintf("assets with a different best model: %d of 7; mean Spearman correlation: %.2f\n\n",
              V[[nm]]$n_differ, V[[nm]]$rho))
}
saveRDS(list(losses = ML, variants = V), "data/fz2f_sensitivity.rds")
