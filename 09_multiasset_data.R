# ============================================================
# Beyond Shortfall — 09: return files for the multi-asset forecasting comparison
#
# The model comparison of steps 02-04 and 08 uses Bitcoin. This step extends it
# to all seven risky assets of the portfolio universe (the two stablecoins have
# near-zero returns and are not meaningful targets for tail forecasting).
#
# Each asset uses its own available history up to the sample end, with the
# first 1095 observations (three years) as the estimation period, exactly as
# for Bitcoin. Six assets start in October 2017 and therefore share Bitcoin's
# out-of-sample period (October 2020 - June 2026); Solana starts in April 2020
# and has a correspondingly shorter out-of-sample period (from April 2023).
#
# Output: MATLAB_DIR/assets/<ASSET>/ret.txt  (MATLAB_DIR set in config.R)
#         (space separated: YYYYMMDD, daily return in PERCENT), the input of
#         10_multiasset_forecasts.m
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
library(data.table)

END    <- as.Date("2026-06-30")
ASSETS <- c("BNB", "BTC", "DOGE", "ETH", "SOL", "TRX", "XRP")
R_IN   <- 1095

ret <- fread("data/coinReturns.csv")
ret[, timestamp := as.Date(timestamp)]

info <- rbindlist(lapply(ASSETS, function(a) {
  x <- ret[!is.na(get(paste0("r", a))) & timestamp <= END,
           .(date = format(timestamp, "%Y%m%d"), r = 100 * get(paste0("r", a)))]
  dir.create(file.path(MATLAB_DIR, "assets", a), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, file.path(MATLAB_DIR, "assets", a, "ret.txt"), sep = " ", col.names = FALSE)
  nz <- sum(x$r != 0)                          # step 02 drops zero returns
  data.table(asset = a, first = x$date[1], last = x$date[nrow(x)], T = nrow(x),
             T_nonzero = nz, P = nz - R_IN, oos_start = x$date[x$r != 0][R_IN + 1])
}))
print(info)

# consistency check: the BTC file must reproduce the single-asset input of step 02
b_new <- fread(file.path(MATLAB_DIR, "assets", "BTC", "ret.txt"))
b_old <- fread(file.path(MATLAB_DIR, "BTC.txt"))
stopifnot(nrow(b_new) == nrow(b_old), max(abs(b_new$V2 - b_old$V2)) < 1e-10)
cat("\nBTC file identical to the single-asset input of step 02.\n")
saveRDS(info, "data/multiasset_info.rds")
