# ============================================================
# Beyond Shortfall — 21: Bitcoin tables of Section 5 for the main scheme
#
# The main forecasting results use models re-estimated every 91 days on a
# rolling window of 1095 days (step 14) and the matching Single-Index
# forecasts (step 15). Steps 04 and 08, run in their place, produce the
# Bitcoin tables for the models estimated once (step 02), with the suffix
# _fixed, which are reported in the appendix. This step re-runs 04 and 08 on
# the re-estimated forecasts and writes the main tables:
#   LaTeX/tables/FZLoss_ES.tex, FZLoss_EU.tex, FZHLoss_ES.tex, FZHLoss_EU.tex
#   LaTeX/tables/Backtest.tex, BacktestCount.tex, RankShift.tex
#   data/mcs_results_beyond.RData, data/backtests.rds
# ============================================================
SCHEME <- "rolling"; source("04_mcs_tables.R")
SCHEME <- "rolling"; source("08_forecast_evaluation.R")
