# ============================================================
# Beyond Shortfall — 18: MQR portfolios with the model chosen per tail  >>> RUN IN MATLAB <<<
#
# THIS STEP IS NOT EXECUTED IN R. It documents step 18 so the numbering stays
# complete (same pattern as 02_run_matlab_forecasts.R).
#
# ------------------------------------------------------------
# WHAT TO RUN
# ------------------------------------------------------------
#   Folder: MATLAB_DIR, set in config.R (the folder with the MATLAB code)
#   File:   18_tailspecific_portfolios.m      <-- open in MATLAB and run it
#
#   Head-less:
#     cd <MATLAB_DIR>
#     matlab -batch "eval(fileread('18_tailspecific_portfolios.m'))"
#
# ------------------------------------------------------------
# WHAT IT DOES
# ------------------------------------------------------------
#   Runs 06_portfolio_weights.m for three designs, each for c = 0, 0.1, ..., 1,
#   weekly rebalancing, with the weight cap 0.3 and without a cap:
#     TailSpec   ES and EU forecast by the models chosen per tail at every
#                rebalancing date (schedule_tailspec.txt, from step 17)
#     SameModel  one model for both tails, chosen by the same rule
#                (schedule_samemodel.txt, from step 17)
#     GARCH-EDF  static ex-ante choice (best in both tails before July 2023)
#   Existing weight files are skipped, so the step can be resumed.
#
# ------------------------------------------------------------
# OUTPUT (in MATLAB_DIR)
# ------------------------------------------------------------
#   weights_<TailSpec|SameModel|GARCH-EDF>[_cap03]_c<cc>.txt, objval_*.txt
#
#   Runtime: about 8-12 minutes per value of c for TailSpec and SameModel and
#   about 1 minute for GARCH-EDF, with 10 parallel workers (about 6 hours).
#
# NEXT STEP
#   19_tailspecific_performance.R
# ============================================================

stop("Step 18 runs in MATLAB: open 18_tailspecific_portfolios.m in MATLAB_DIR (config.R) — see the comments in this file.")
