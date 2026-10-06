# ============================================================
# Beyond Shortfall — 06: portfolio weights  >>> RUN IN MATLAB <<<
#
# THIS STEP IS NOT EXECUTED IN R. It documents step 06 so the numbering
# 01 -> 07 stays complete. (Same pattern as 02_run_matlab_forecasts.R.)
#
# ------------------------------------------------------------
# WHAT TO RUN
# ------------------------------------------------------------
#   Folder: MATLAB_DIR, set in config.R (the folder with the MATLAB code)
#   File:   06_portfolio_weights.m       <-- open in MATLAB and run it
#
#   Head-less:
#     cd <MATLAB_DIR>
#     matlab -batch "eval(fileread('06_portfolio_weights.m'))"
#   (eval(fileread(...)), not run(...) — the file name starts with a digit)
#
#   To run a subset, define the two variables before the eval, e.g.
#     MODELS = {'GARCH-Skt'}; CVALS = 0:0.1:1;
#     eval(fileread('06_portfolio_weights.m'))
#   Existing weight files are skipped, so the sweep can be resumed.
#
#   Configurations used in the paper (four models, c = 0, 0.1, ..., 1):
#     MODELS = {'RW125','GARCH-Skt','GARCH-FZ','Hybrid'}; CVALS = 0:0.1:1;
#     for WCAP = [1 0.5 0.3]      % weekly rebalancing; no cap, caps 0.5 and 0.3
#         A_WIN = 7;  eval(fileread('06_portfolio_weights.m'))
#     end
#     for A_WIN = [14 30 91]      % other rebalancing frequencies, no cap
#         WCAP = 1;  eval(fileread('06_portfolio_weights.m'))
#     end
#   Approximate time per value of c with 8-10 workers: RW125 seconds,
#   GARCH-Skt 2-3 minutes, GARCH-FZ 10-14 minutes, Hybrid 20-26 minutes.
#
# ------------------------------------------------------------
# THE OPTIMISATION
# ------------------------------------------------------------
#   Thesis eq. (5.1) — maximising quantile returns (MQR):
#       max_w  E[y_p | y_p < q_lower] + c * E[y_p | y_p > q_upper]
#   With ES and EU written as positive magnitudes this is
#       min_w  ES(w) - c * EU(w)     s.t. sum(w) = 1, 0 <= w <= 1
#     c = 0 : pure minimum-ES portfolio (as in the previous paper)
#     c = 1 : downside and upside weighted equally
#     c is swept over 0, 0.1, ..., 1.0 (investor risk-return preference)
#
#   ES and EU are conditional FORECASTS of the portfolio return series
#   (not historical quantiles) from the models compared in step 04:
#   RW125, GARCH-Skt, Hybrid, GARCH-FZ.
#
#   Rolling scheme: 1095-day estimation window, rebalanced every 7 days,
#   alpha = 0.05, 156 rebalancing dates — as in the previous paper.
#   Each objective evaluation refits the model on the portfolio return
#   series implied by the candidate weights, so the optimisation is run in
#   parallel over rebalancing dates (parpool, 10 workers).
#
#   Starting values: 1/N plus every single-asset corner, best one kept and
#   then refined by fmincon. ES and EU both scale with volatility, so the
#   objective is frequently minimised at a vertex and a corner scan is a
#   more reliable start than 1/N alone.
#
# ------------------------------------------------------------
# INPUT / OUTPUT
# ------------------------------------------------------------
#   Input : coin_ret.txt  (written by 05_portfolio_data.R)
#           2191 x 9 daily returns in percent, 2020-07-01 .. 2026-06-30
#   Output: weights_<model>_c<cc>.txt   156 x 9 optimal weights
#           objval_<model>_c<cc>.txt    156 x 3  (objective, ES, EU)
#           cc is c*10, e.g. c00 = 0.0, c05 = 0.5, c10 = 1.0
#
#   Runtime (12-core laptop, 10 workers): about 6 s per (model, c) for
#   RW125 and 2-4 min for the GARCH/GAS models; roughly 1.7 h for the full
#   4 models x 11 c values sweep.
#
# NEXT STEP
#   07_portfolio_performance.R — out-of-sample performance, tables, figures
# ============================================================

stop("Step 06 runs in MATLAB: open 06_portfolio_weights.m in MATLAB_DIR (config.R) — see the comments in this file.")
