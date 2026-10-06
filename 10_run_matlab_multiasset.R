# ============================================================
# Beyond Shortfall — 10: multi-asset OOS VaR/ES forecasts  >>> RUN IN MATLAB <<<
#
# THIS STEP IS NOT EXECUTED IN R. It documents step 10 so the numbering stays
# complete (same pattern as 02_run_matlab_forecasts.R).
#
# ------------------------------------------------------------
# WHAT TO RUN
# ------------------------------------------------------------
#   Folder: MATLAB_DIR, set in config.R (the folder with the MATLAB code)
#   File:   10_multiasset_forecasts.m      <-- open in MATLAB and run it
#
#   Head-less:
#     cd <MATLAB_DIR>
#     matlab -batch "eval(fileread('10_multiasset_forecasts.m'))"
#   (eval(fileread(...)), not run(...) — the file name starts with a digit)
#
#   To run a subset of assets, define ASSETS_ first, e.g.
#     ASSETS_ = {'ETH','SOL'}; eval(fileread('10_multiasset_forecasts.m'))
#
# ------------------------------------------------------------
# WHAT IT DOES
# ------------------------------------------------------------
#   Loops over the seven risky assets prepared by 09_multiasset_data.R
#   (BNB, BTC, DOGE, ETH, SOL, TRX, XRP) and runs 02_BTC_forecast_tails.m for
#   each of them, with IN_FILE and OUT_DIR pointing to assets/<ASSET>/.
#   Models, quantile levels, tails and the estimation scheme are identical to
#   the Bitcoin comparison: ten models (RW-125, RW-250, RW-500, GARCH-N,
#   GARCH-Skt, GARCH-EDF, FZ-2F, FZ-1F, GARCH-FZ, Hybrid), tau in
#   {0.01, 0.025, 0.05, 0.1}, lower tail ('low') and upper tail ('up', on the
#   negated returns), estimated once on the first 1095 observations and
#   filtered forward.
#
#   The BTC run reproduces the files of step 02 and serves as a check.
#
# ------------------------------------------------------------
# OUTPUT (per asset, in assets/<ASSET>/)
# ------------------------------------------------------------
#   VaR_<low|up>_<001|0025|005|01>.txt, ES_<...>.txt   P x 10 forecasts
#   oos_dates.txt, oos_ret.txt                         OOS dates and returns
#
#   Runtime: about one minute per asset.
#
# NEXT STEP
#   11_multiasset_evaluation.R — Single-Index model, MCS, backtests and the
#   cross-asset summary tables for the paper
# ============================================================

stop("Step 10 runs in MATLAB: open 10_multiasset_forecasts.m in MATLAB_DIR (config.R) — see the comments in this file.")
