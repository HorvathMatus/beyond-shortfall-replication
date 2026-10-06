# ============================================================
# Beyond Shortfall — 05: build the asset-return file for the portfolio part
#
# Writes coin_ret.txt (space separated, returns in PERCENT, no header/date)
# into the MATLAB folder, where step 06 optimises the portfolio weights.
#
# Sample: the last SIX YEARS of the data set, 2020-07-01 .. 2026-06-30.
# The paper's asset filter is "market cap >= USD 10bn AND at least six years
# of history", so six years is exactly the window over which all nine
# selected coins exist (SOL, the shortest, starts 2020-04-10).
# With a three-year (1095-day) estimation window this leaves three years of
# out-of-sample portfolio performance, mirroring the previous paper.
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
library(data.table)

START <- as.Date("2020-07-01")
END   <- as.Date("2026-06-30")

ret <- fread("data/coinReturns.csv")
ret[, timestamp := as.Date(timestamp)]

# the nine coins selected in 01_download_data.R
COINS <- c("ADA","BNB","BTC","DOGE","ETH","SOL","TRX","USDC","USDT","XRP")
COINS <- intersect(COINS, sub("^r", "", grep("^r", names(ret), value = TRUE)))
cat("Coins:", paste(COINS, collapse = ", "), "\n")

rcols <- paste0("r", COINS)
dt <- ret[timestamp >= START & timestamp <= END, c("timestamp", rcols), with = FALSE]

# returns in percent
for (cc in rcols) dt[[cc]] <- 100 * dt[[cc]]

n_na <- colSums(is.na(dt[, ..rcols]))
if (any(n_na > 0)) {
  cat("Missing values per coin:\n"); print(n_na)
  stop("missing returns inside the portfolio sample")
}

cat(sprintf("Sample %s .. %s, T = %d days, N = %d assets\n",
            as.character(min(dt$timestamp)), as.character(max(dt$timestamp)),
            nrow(dt), length(rcols)))
cat(sprintf("Rebalancing every 7 days, 1095-day window -> %d periods\n",
            floor((nrow(dt) - 1095) / 7)))

fwrite(dt[, ..rcols], file.path(MATLAB_DIR, "coin_ret.txt"),
       sep = " ", col.names = FALSE)
fwrite(dt[, .(date = format(timestamp, "%Y%m%d"))],
       file.path(MATLAB_DIR, "coin_dates.txt"), col.names = FALSE)
saveRDS(list(dates = dt$timestamp, coins = COINS,
             ret = as.matrix(dt[, ..rcols])), "data/portfolio_input.rds")

cat("\nWritten:", file.path(MATLAB_DIR, "coin_ret.txt"), "and coin_dates.txt\n")
cat("Column order:", paste(COINS, collapse = " "), "\n")

# ---- descriptive statistics for the paper (both tails) ---------------------
# Portfolio sample for all nine assets, plus Bitcoin over the longer
# forecasting sample used in the model comparison. Tail measures are the
# empirical 5% VaR/ES of the lower tail and VaR/EU of the upper tail.
desc_row <- function(x) {
  x <- x[!is.na(x)]
  ql <- quantile(x, 0.05, type = 7); qu <- quantile(x, 0.95, type = 7)
  es <- mean(x[x <= ql]); eu <- mean(x[x >= qu])
  m <- mean(x); s <- sd(x)
  c(N = length(x), Mean = m, SD = s, Min = min(x), Max = max(x),
    Skew = mean((x - m)^3) / s^3, ExKurt = mean((x - m)^4) / s^4 - 3,
    ES = es, EU = eu, Ratio = abs(eu / es))
}
TEX_TAB <- TEX_DIR
if (dir.exists(TEX_TAB)) {
  M <- t(sapply(rcols, function(cc) desc_row(dt[[cc]])))
  rownames(M) <- COINS
  btc_full <- 100 * ret[!is.na(rBTC) & timestamp <= END, rBTC]
  M <- rbind(M, "BTC$^{\\dagger}$" = desc_row(btc_full))
  L <- c("\\begin{table}[ht]", "\\centering",
    "\\caption{Descriptive statistics of daily returns (\\%).}",
    "\\label{tab:descriptive}", "\\begin{center}",
    "\\begin{tabular}{lrrrrrrrrrr}", "\\toprule",
    "Asset & $N$ & Mean & SD & Min & Max & Skew & Ex.\\ kurt. & ES$_{0.05}$ & EU$_{0.05}$ & $|$EU/ES$|$ \\\\",
    "\\midrule")
  for (i in seq_len(nrow(M))) {
    if (i == nrow(M)) L <- c(L, "\\midrule")
    v <- M[i, ]
    L <- c(L, sprintf("%s & %d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\",
      rownames(M)[i], as.integer(v["N"]), v["Mean"], v["SD"], v["Min"], v["Max"],
      v["Skew"], v["ExKurt"], v["ES"], v["EU"], v["Ratio"]))
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\small",
    paste0("\\textit{Notes: Daily returns in percent over the portfolio sample from ",
      format(START, "%B %e, %Y"), " to ", format(END, "%B %e, %Y"),
      ". BTC$^{\\dagger}$ refers to Bitcoin over the longer forecasting sample from October 2, 2017. ",
      "ES$_{0.05}$ and EU$_{0.05}$ are the empirical Expected Shortfall and Expected Upside at the 5\\% level, ",
      "that is the mean return in the worst and in the best 5\\% of days; $|$EU/ES$|$ above one indicates ",
      "an upper tail larger in magnitude than the lower tail at the 5\\% level. Ex.\\ kurt.\\ is excess kurtosis.}"),
    "\\end{table}")
  writeLines(gsub("  +", " ", L), file.path(TEX_TAB, "Descriptive.tex"))
  cat("\nDescriptive table written to", file.path(TEX_TAB, "Descriptive.tex"), "\n")
  print(round(M, 2))
}
