# Beyond Shortfall — replication package

Code for **"Beyond Shortfall: Expected Shortfall and Upside Estimators in
Cryptocurrency Portfolios"** by Matúš Horváth and Tomáš Výrost.

The package reproduces every table and figure of the paper, starting from the
download of daily cryptocurrency prices. The pipeline consists of numbered
steps: the R scripts in the root folder and the MATLAB scripts in `matlab/`.
Each MATLAB step has an R file with the same number that documents it and
stops with a message if it is sourced.

## Requirements

**R** (tested with R 4.4.2) and the packages

| package    | version | used for                                       |
|------------|---------|------------------------------------------------|
| data.table | 1.16.4  | data handling                                  |
| rjson      | 0.2.23  | download from the CoinMarketCap data API       |
| MCS        | 0.2.0   | model confidence set                           |
| rugarch    | 1.5.3   | VaR coverage tests (`VaRTest`)                 |
| esback     | 0.3.1   | expected shortfall regression test             |
| sandwich   | 3.1.1   | HAC variances (`lrvar`)                        |
| lpSolve    | 5.6.23  | historical minimum-CVaR and Omega portfolios   |
| MAVE, MASS | 1.3.12, 7.3.61 | single-index expectile model            |

```r
install.packages(c("data.table", "rjson", "MCS", "rugarch", "esback",
                   "sandwich", "lpSolve", "MAVE", "MASS"))
```

The base package `parallel` is used to run assets and blocks in parallel
(`mclapply`, which forks processes on macOS and Linux).

**MATLAB** (tested with R2023b) with the Optimization Toolbox, the Statistics
and Machine Learning Toolbox and the Parallel Computing Toolbox, plus
third-party code that is not redistributed here: the replication files of
Patton, Ziegel and Chen (2019), the UCSD GARCH toolbox, the Econometrics
Toolbox of J.P. LeSage and the MFE Toolbox. `matlab/thirdparty/README.md`
lists the files and where to obtain them.

## Data

Daily prices are downloaded by `01_download_data.R` from the public
CoinMarketCap data API (no key required) for October 1, 2017 to June 30, 2026,
and are not redistributed. CoinMarketCap occasionally revises historical
prices, so a download made at a different date can differ slightly from the
data used in the paper.

## Folder layout

```
config.R          paths (MATLAB_DIR, TEX_DIR, FIG_DIR); edit here only
01_...22_*.R      numbered steps (R)
R/                helper functions used by the R steps
matlab/           MATLAB steps 02, 06, 10, 14, 18 and their functions;
                  also holds the MATLAB input and output files
matlab/thirdparty third-party MATLAB code (to be added, see its README)
data/             intermediate R results (created by the scripts)
output/tables     LaTeX tables (created by the scripts)
output/figures    figures (created by the scripts)
```

Run the R scripts from the root folder, for example `Rscript 04_mcs_tables.R`,
and the MATLAB scripts from `matlab/`, for example

```
cd matlab
matlab -batch "eval(fileread('02_BTC_forecast_tails.m'))"
```

(`eval(fileread(...))` is needed because the file names start with a digit.)
`bs_paths.m` locates the code folder and adds the third-party code to the
path; machine-specific toolbox folders can be listed in `matlab/local_paths.m`
(see `matlab/local_paths_example.m`).

## Pipeline

Approximate run times refer to a 12-core Apple M-series machine with 8 to 10
parallel MATLAB workers.

| step | language | what it does | time |
|------|----------|--------------|------|
| 01 | R | download prices, select the nine cryptocurrencies, write the Bitcoin input of step 02 (`matlab/BTC.txt`) | minutes |
| 02 | MATLAB | Bitcoin VaR/ES forecasts of ten models, both tails, four levels, models estimated once | 1 min |
| 03 | R | Single-Index expectile forecasts for Bitcoin | 5 min |
| 04 | R | Bitcoin MCS tables, models estimated once (appendix tables `*_fixed`) | 3 min |
| 05 | R | portfolio data (`matlab/coin_ret.txt`), Table of descriptive statistics | seconds |
| 06 | MATLAB | MQR portfolio weights; configurations listed in `06_run_matlab_portfolio.R` | about 2 days in total |
| 07 | R | portfolio performance tables and figures | 2 min |
| 08 | R | Bitcoin backtests and rank shifts, models estimated once (`*_fixed`) | 1 min |
| 09 | R | return files of the seven risky assets (`matlab/assets/`) | seconds |
| 10 | MATLAB | forecasts for the seven assets, models estimated once | 6 min |
| 11 | R | Single-Index forecasts, MCS and backtests for the seven assets, models estimated once (`*_fixed`) | 40 min |
| 12 | R | ex-ante choice of the preference coefficient | 1 min |
| 13 | R | historical minimum-CVaR, Omega and Rachev portfolios | 20 min |
| 14 | MATLAB | forecasts for the seven assets, models re-estimated every 91 days (`matlab/assets_roll/`) | 20 min |
| 15 | R | Single-Index forecasts with re-estimation, evaluation, main cross-asset tables, comparison of the schemes | 55 min |
| 16 | R | bootstrap inference (9,999 replications, Holm adjustment) on the differences between the tails and between the two schemes | 6 min |
| 17 | R | ex-ante choice of the forecasting model for each tail (schedules for step 18) | 1 min |
| 18 | MATLAB | MQR portfolios with the model chosen per tail | 6 h |
| 19 | R | performance of the tail-specific portfolios | 1 min |
| 20 | R | volatility and returns around the switches of the unconstrained portfolios (calendar block bootstrap) | seconds |
| 21 | R | main Bitcoin tables (models re-estimated every 91 days): steps 04 and 08 on the step-14 forecasts | 4 min |
| 22 | R | sensitivity of the cross-asset rankings to the replaced FZ-2F forecasts (numbers in Section 5.5) | 1 min |

Run the steps in the order 01–22. The MATLAB weight files of steps 06 and 18 are
skipped if they already exist, so these steps can be interrupted and resumed.

## Tables and figures

| paper | file in `output/` | produced by |
|-------|-------------------|-------------|
| descriptive statistics | `tables/Descriptive.tex` | 05 |
| Bitcoin, lower and upper tail, FZ(0) and FZ(1/2) scores | `tables/FZLoss_ES.tex`, `FZLoss_EU.tex`, `FZHLoss_ES.tex`, `FZHLoss_EU.tex` | 21 |
| Bitcoin backtests and rank shift | `tables/Backtest.tex`, `BacktestCount.tex`, `RankShift.tex` | 21 |
| seven cryptocurrencies, by model and by asset | `tables/MultiAsset.tex`, `MultiAssetByAsset.tex` | 15 |
| models re-estimated vs estimated once | `tables/Rolling.tex` | 15 |
| inference on the differences between the tails and between the schemes, cumulative loss differences, rank shifts (appendix) | `tables/TailInference.tex`, `figures/tail_cumdiff.png`, `tables/TailInferenceRanks.tex` | 16 |
| portfolio performance, concentration, weight caps, rebalancing, transaction costs | `tables/Performance.tex`, `Concentration.tex`, `CapSummary.tex`, `PerformanceCap.tex`, `PerformanceRebal.tex`, `PerformanceNet.tex` | 07 |
| switches of the unconstrained portfolios | `tables/SwitchAnalysis.tex` | 20 |
| sensitivity to the replaced FZ-2F forecasts | printed by step 22 (`data/fz2f_sensitivity.rds`) | 22 |
| historical reward–risk benchmarks | `tables/RewardRisk.tex` | 13 |
| model chosen per tail | `tables/TailSpecific.tex`, `figures/tailspecific_sr.png` | 19 |
| ex-ante preference coefficient | `tables/AdaptiveC.tex` | 12 |
| portfolio figures | `figures/c_frontier.png`, `weights_models.png`, `rebalancing_heatmap.png`, `cumulative_wealth_models.png`, `cumulative_wealth_cap.png`, `weights_rw125.png` | 07 |
| Appendix B (models estimated once) | `tables/*_fixed.tex` | 04, 08, 11 |

## Reproducibility notes

* Random numbers: the MCS bootstrap uses `set.seed(123)` for every comparison,
  the Single-Index model `set.seed(123)`, the ESR backtest (`esback`, which
  draws random starting values) `set.seed(123)` before every test, the
  tail-inference bootstrap `set.seed(2026)`, the Rachev multistart
  `set.seed(2024)` and the block bootstrap of the switch analysis `set.seed(7)`.
* Bootstrap p-values in step 16 invert the percentile intervals with
  the finite-sample correction p = 2 min(1 + #{θ* ≤ 0}, 1 + #{θ* ≥ 0}) / (B + 1);
  they are adjusted with the Holm method within families of related tests.
* All estimation choices use in-sample information only. Forecasts that are
  non-finite or violate ES ≤ VaR < 0 are replaced in real time by the model's
  most recent usable forecast (`realtime_guard` in `R/forecast_eval.R`).
* With the forecasting window covering a single block, step 14 reproduces
  step 02 exactly; the multi-asset run of step 10 reproduces the Bitcoin
  forecasts of step 02 exactly.
* The MATLAB optimisers (`fminsearch`, `fminunc`, `fmincon`) are deterministic,
  but results can differ in the last digits across MATLAB versions and
  hardware.

## References

Patton, A.J., Ziegel, J.F. and Chen, R. (2019). Dynamic semiparametric models
for expected shortfall (and value-at-risk). *Journal of Econometrics* 211(2),
388–413.

Horváth, M. and Výrost, T. (2025). No shortfall of ES estimators: Insights from
cryptocurrency portfolios. *Finance Research Letters* 73, 106685.

## License

The code in this repository is released under the MIT License (see
`LICENSE`). The license does not cover the third-party MATLAB code listed in
`matlab/thirdparty/README.md` or the CoinMarketCap data, which are not part of
the repository.
