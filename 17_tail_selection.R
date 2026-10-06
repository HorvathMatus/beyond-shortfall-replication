# ============================================================
# Beyond Shortfall — 17: ex-ante choice of the forecasting model for each tail
#
# Experiment suggested in review: does it pay to forecast the ES with the
# model best suited to the lower tail and the EU with the model best suited to
# the upper tail? The choice must use only information available at each
# rebalancing date. At every weekly rebalancing date of the portfolio
# application we rank the candidate models by their average FZ0 loss at
# tau = 0.05 over the preceding LOOKBACK days, separately in each tail and for
# each of the seven risky assets, using the forecasts of step 14 (models
# re-estimated every 91 days, i.e. produced in real time), and average the
# ranks over the assets. Two schedules result:
#
#   tail-specific  ES from the model with the best lower-tail rank,
#                  EU from the model with the best upper-tail rank
#   same-model     the model with the best average of the two ranks, for both
#
# Comparing the two isolates the value of choosing the model tail by tail.
# A third, static benchmark uses the model with the best ranks over the whole
# pre-evaluation period (October 2020 to June 2023), which is GARCH-EDF in
# both tails (printed below).
#
# Candidates: the models that the portfolio code (portfolio_tails.m)
# implements: RW-125, GARCH-N, GARCH-Skt, GARCH-EDF, GARCH-FZ, Hybrid.
#
# Inputs : assets_roll/<A>/ (step 14), data/portfolio_input.rds
# Outputs: MATLAB/Beyond_Shortfall/schedule_tailspec.txt, schedule_samemodel.txt
#          (one model string per rebalancing date, read by 18_*.m),
#          data/tail_selection.rds
# ============================================================
rm(list=ls())
source("config.R")   # MATLAB_DIR, TEX_DIR, FIG_DIR
suppressPackageStartupMessages(library(data.table))
source("R/forecast_eval.R")

ROOT_ROLL  <- file.path(MATLAB_DIR, "assets_roll")
ASSETS     <- c("BNB", "BTC", "DOGE", "ETH", "SOL", "TRX", "XRP")
CAND       <- c("RW-125" = "RW125", "GARCH-N" = "GARCH-N", "GARCH-Skt" = "GARCH-Skt",
                "GARCH-EDF" = "GARCH-EDF", "GARCH-FZ" = "GARCH-FZ", "Hybrid" = "Hybrid")
TAU        <- 0.05
LOOKBACK   <- 365          # days of past forecast losses used for the choice
MIN_DAYS   <- 100          # an asset enters the ranking with at least this many days
W_EST <- 1095; A_WIN <- 7

# ---- daily FZ0 losses of the candidates, by asset and tail -------------------
k <- match(TAU, TAUS)
LOSS <- lapply(setNames(ASSETS, ASSETS), function(a) {
  d <- file.path(ROOT_ROLL, a)
  dates <- as.Date(as.character(scan(file.path(d, "oos_dates.txt"), quiet = TRUE)), "%Y%m%d")
  r <- scan(file.path(d, "oos_ret.txt"), quiet = TRUE)
  out <- lapply(c(low = "low", up = "up"), function(tl) {
    y <- if (tl == "low") r else -r
    fc <- read_fc(d, tl, SUFFIX[k])
    sapply(names(CAND), function(m) fz_loss(y, fc$v[, m], fc$e[, m], TAU))
  })
  list(dates = dates, low = out$low, up = out$up)
})

# ---- rebalancing dates (as in step 06) -----------------------------------------
inp <- readRDS("data/portfolio_input.rds")
NPER <- floor((nrow(inp$ret) - W_EST) / A_WIN)
dec_dates <- inp$dates[W_EST + (0:(NPER - 1)) * A_WIN]      # last day of each estimation window

avg_ranks <- function(from, to) {
  rk <- list()
  for (a in ASSETS) {
    i <- which(LOSS[[a]]$dates > from & LOSS[[a]]$dates <= to)
    if (length(i) < MIN_DAYS) next
    for (tl in c("low", "up")) rk[[tl]][[a]] <- rank(colMeans(LOSS[[a]][[tl]][i, , drop = FALSE]), ties.method = "first")
  }
  sapply(c(low = "low", up = "up"), function(tl) rowMeans(do.call(cbind, rk[[tl]])))
}

sel <- rbindlist(lapply(seq_len(NPER), function(l) {
  ar <- avg_ranks(dec_dates[l] - LOOKBACK, dec_dates[l])
  es_m <- names(which.min(ar[, "low"])); eu_m <- names(which.min(ar[, "up"]))
  same <- names(which.min(rowMeans(ar)))
  data.table(period = l, date = dec_dates[l], es = es_m, eu = eu_m, same = same)
}))

cat("Choices of the tail-specific schedule (ES model x EU model):\n")
print(sel[, .N, by = .(es, eu)][order(-N)])
cat("\nChoices of the same-model schedule:\n"); print(sel[, .N, by = same][order(-N)])
cat(sprintf("\nPeriods in which the two tails get different models: %d of %d\n", sel[es != eu, .N], NPER))

# static benchmark: whole pre-evaluation period
pre <- avg_ranks(as.Date("2020-09-30"), as.Date("2023-06-30"))
cat("\nAverage ranks over October 2020 - June 2023 (static benchmark):\n"); print(round(pre, 2))

# ---- schedules for MATLAB ---------------------------------------------------------
tailspec <- ifelse(sel$es == sel$eu, CAND[sel$es], paste0(CAND[sel$es], "|", CAND[sel$eu]))
writeLines(tailspec, file.path(MATLAB_DIR, "schedule_tailspec.txt"))
writeLines(unname(CAND[sel$same]), file.path(MATLAB_DIR, "schedule_samemodel.txt"))
saveRDS(list(selection = sel, static_ranks = pre, lookback = LOOKBACK), "data/tail_selection.rds")
cat("\nSchedules written to", MATLAB_DIR, "\n")
