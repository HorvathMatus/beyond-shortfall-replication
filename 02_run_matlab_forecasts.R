# ============================================================
# Beyond Shortfall — 02: OOS VaR/ES forecasts  >>> RUN IN MATLAB <<<
#
# THIS STEP IS NOT EXECUTED IN R.
# This file only documents step 02 of the pipeline so that the numbering
# 01 -> 02 -> 03 -> 04 stays complete and readable.
#
# ------------------------------------------------------------
# WHAT TO RUN
# ------------------------------------------------------------
#   Folder: MATLAB_DIR, set in config.R (the folder with the MATLAB code)
#   File:   02_BTC_forecast_tails.m        <-- open this in MATLAB and run it
#
#   From a terminal it can also be run head-less with:
#     cd <MATLAB_DIR>
#     matlab -batch "eval(fileread('02_BTC_forecast_tails.m'))"
#
#   NOTE: use eval(fileread(...)), NOT run('02_...m'). MATLAB cannot run()
#   a script whose file name starts with a digit (it is not a valid MATLAB
#   identifier); from the Editor the Run button works fine.
#
# ------------------------------------------------------------
# WHY MATLAB
# ------------------------------------------------------------
#   The VaR/ES models come from the original Patton, Ziegel and Chen (2019)
#   replication code, which is written in MATLAB. The R re-implementation
#   did not reproduce it correctly, so all
#   forecasting is done with Patton's original MATLAB routines.
#
# ------------------------------------------------------------
# WHAT 02_BTC_forecast_tails.m DOES
# ------------------------------------------------------------
#   Input : BTC.txt  (written by 01_download_data.R;
#                     columns = date YYYYMMDD, daily return in PERCENT)
#   Sample: 2017-10-02 .. 2026-06-30, T = 3194 obs
#           in-sample R = 1095 (3 years), out-of-sample P = 2099
#
#   Loops over tau in {0.01, 0.025, 0.05, 0.1} and over BOTH tails:
#     'low' = lower tail -> Expected Shortfall (ES), original returns
#     'up'  = upper tail -> Expected Upside   (EU), NEGATED returns
#   The upper tail simply applies the identical lower-tail machinery to -r,
#   so outputs are stored in the negated space. To read them as upside:
#     VaR_up = -VaR_saved ,  EU = -ES_saved
#
#   Models (10, in this column order in every output file):
#     RW-125, RW-250, RW-500, GARCH-N, GARCH-Skt, GARCH-EDF,
#     FZ-2F, FZ-1F, GARCH-FZ, Hybrid
#
#   Requires these toolboxes on the MATLAB path (added by the script itself):
#     mfe-toolbox-main, ucsd_garch (garchlikelihood), etoolbox (tdis_inv)
#
# ------------------------------------------------------------
# OUTPUT (written into the same MATLAB folder)
# ------------------------------------------------------------
#   VaR_low_001.txt  VaR_low_0025.txt  VaR_low_005.txt  VaR_low_01.txt
#   ES_low_001.txt   ES_low_0025.txt   ES_low_005.txt   ES_low_01.txt
#   VaR_up_001.txt   VaR_up_0025.txt   VaR_up_005.txt   VaR_up_01.txt
#   ES_up_001.txt    ES_up_0025.txt    ES_up_005.txt    ES_up_01.txt
#     -> each is a P x 10 csv matrix (rows = OOS days, cols = models above)
#   oos_dates.txt , oos_ret.txt
#     -> OOS dates and realised returns (ORIGINAL sign), used by steps 03-04
#
#   Runtime: about 1 minute for all 8 tail/tau combinations.
#
# ------------------------------------------------------------
# NEXT STEPS
# ------------------------------------------------------------
#   03_single_index.R  — Single Index expectile model (Jiang et al., 2020),
#                        run in R, same taus and both tails
#   04_mcs_tables.R    — FZ and AL losses + Model Confidence Set tables
# ============================================================

stop("Step 02 runs in MATLAB: open 02_BTC_forecast_tails.m in MATLAB_DIR (config.R) — see the comments in this file.")
