# ============================================================
# Beyond Shortfall — 14: forecasts with rolling re-estimation  >>> RUN IN MATLAB <<<
#
# THIS STEP IS NOT EXECUTED IN R. It documents step 14 so the numbering stays
# complete (same pattern as 02_run_matlab_forecasts.R).
#
# ------------------------------------------------------------
# WHAT TO RUN
# ------------------------------------------------------------
#   Folder: MATLAB_DIR, set in config.R (the folder with the MATLAB code)
#   File:   14_rolling_forecasts.m      <-- open in MATLAB and run it
#           (uses the function file fit_forecast_block.m in the same folder)
#
#   Head-less:
#     cd <MATLAB_DIR>
#     matlab -batch "eval(fileread('14_rolling_forecasts.m'))"
#
#   Optional overrides, defined before the call:
#     ASSETS_   = {'BTC'};    % subset of assets (default: all seven)
#     K_        = 91;         % re-estimation interval in days (default 91)
#     NWORKERS_ = 8;          % parallel workers (default 8)
#
# ------------------------------------------------------------
# WHAT IT DOES
# ------------------------------------------------------------
#   Robustness check for the estimation scheme. Steps 02 and 10 estimate each
#   model once, on the first 1095 observations, and filter it forward with
#   fixed parameters. Here all ten MATLAB models are re-estimated every K_ = 91
#   days on a rolling window of the most recent 1095 observations and filtered
#   forward until the next re-estimation (blocks run in parallel with parfor).
#   Estimation is done by fit_forecast_block.m, the same function as in step
#   02, and every choice uses in-sample information only: a window whose
#   estimate has unusable in-sample fitted values (estimate_ok.m) would fall
#   back on the parameters of the most recent earlier window (this never
#   occurs in our data). Forecasts that violate ES <= VaR < 0 on a given day
#   are replaced in real time in R (R/forecast_eval.R).
#
#   Check: with K_ larger than the evaluation period (a single block) the
#   script reproduces step 02.
#
# ------------------------------------------------------------
# OUTPUT (per asset, in assets_roll/<ASSET>/)
# ------------------------------------------------------------
#   VaR_<low|up>_<001|0025|005|01>.txt, ES_<...>.txt   P x 10 forecasts
#   oos_dates.txt, oos_ret.txt                         as in step 10
#
#   Runtime: about 3-4 minutes per asset with 8 workers.
#
# NEXT STEP
#   15_rolling_evaluation.R — rolling Single-Index forecasts, MCS, backtests
#   and the comparison with the fixed-parameter scheme
# ============================================================

stop("Step 14 runs in MATLAB: open 14_rolling_forecasts.m in MATLAB_DIR (config.R) — see the comments in this file.")
