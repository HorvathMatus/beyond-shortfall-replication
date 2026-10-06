# ============================================================
# Beyond Shortfall — 01: download crypto data
#
# Sample end:            2026-06-30
# Selection criteria (as in the No Shortfall paper, updated):
#   - market cap >= USD 10 billion (current, at download date)
#   - at least 6 years of daily history before the sample end,
#     i.e. first observation on or before 2020-06-30
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
library(data.table)
source('R/coinMarketCap.R')

END_DATE  <- as.Date('2026-06-30')   # sample end
MIN_START <- as.Date('2020-06-30')   # first obs must be <= this (6y of history)
CAP_MIN   <- 10e9

dir.create('data', showWarnings = FALSE)

# --== 1) current top-100 by market cap =====================================
if (!file.exists('data/coinList.csv')) fwrite(getCoinMarketCapList(),'data/coinList.csv')
clist <- fread('data/coinList.csv')
cap10 <- clist[marketCap >= CAP_MIN][order(-marketCap)]
cat(sprintf('Coins with market cap >= USD %.0f bn: %d\n', CAP_MIN/1e9, nrow(cap10)))

# --== 2) history check: first obs on/before MIN_START? =====================
# One narrow window query per coin around the 6-year threshold.
# firstObs = NA means no data until 2020-07-15 -> history too short.
if (!file.exists('data/historyCheck.csv')) {
  hist_check <- rbindlist(lapply(cap10$id, function(i) {
    tmp <- tryCatch(
      getCoinMarketCap(ids=i, timeStart='2020-06-01', timeEnd='2020-07-15', verbose=FALSE),
      error=function(e) NULL)
    data.table(id = i,
               firstObs = if (is.null(tmp) || nrow(tmp)==0) as.Date(NA) else min(tmp$timestamp))
  }))
  fwrite(hist_check,'data/historyCheck.csv')
}
hist_check <- fread('data/historyCheck.csv')
hist_check[, firstObs := as.Date(firstObs)]

cap10 <- merge(cap10, hist_check, by='id', all.x=TRUE)[order(-marketCap)]
cap10[, has6y := !is.na(firstObs) & firstObs <= MIN_START]

cat('\n--- Screening table (market cap >= 10 bn) ---\n')
print(cap10[,.(symbol, name, marketCap_bn=round(marketCap/1e9,1), firstObs, has6y)])
sel <- cap10[has6y==TRUE]
cat(sprintf('\nQualifying coins (cap >= 10bn & history >= 6y): %d\n', nrow(sel)))
cat(paste(sel$symbol, collapse=', '), '\n')

# --== 3) full daily history for qualifying coins ===========================
if (!file.exists('data/coinPrices.csv')) {
  fwrite(getCoinMarketCap(ids=sel$id, timeStart='2017-10-01', timeEnd=END_DATE),
         'data/coinPrices.csv')
}
cprc <- merge(fread('data/coinPrices.csv')[,c('id','timestamp','close')],
              clist[,c('id','symbol')], by='id', all.x=TRUE, all.y=FALSE)[,c('symbol','timestamp','close')]
cprc <- dcast(cprc, timestamp~symbol, value.var='close')

# coins listed after 2017-10-01 have leading NAs; returns start at first obs
for (i in names(cprc)[-1])
  cprc[[paste0('r',i)]] <- c(NA, cprc[[i]][2:nrow(cprc)]/cprc[[i]][1:(nrow(cprc)-1)]-1)
fwrite(cprc,'data/coinReturns.csv')

ret_cols <- paste0('r', sel$symbol)
data <- cprc[-1, c('timestamp', ret_cols), with=FALSE]
setnames(data, c('Date', sel$symbol))
saveRDS(data, 'data/coin_return.rds')

fwrite(cprc[-1, c('timestamp','rBTC')], 'data/BTC.txt', col.names=FALSE)
# input of step 02 (MATLAB): date YYYYMMDD and daily return in percent
btc <- cprc[-1, .(date = format(as.Date(timestamp), "%Y%m%d"), r = 100 * rBTC)]
fwrite(btc, file.path(MATLAB_DIR, "BTC.txt"), sep = " ", col.names = FALSE)

cat(sprintf('\nSaved: data/coinPrices.csv, data/coinReturns.csv, data/coin_return.rds (until %s)\n',
            as.character(END_DATE)))
