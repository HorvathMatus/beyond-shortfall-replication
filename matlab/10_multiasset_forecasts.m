% Beyond Shortfall - 10: OOS VaR/ES forecasts for all seven risky assets
%
% Runs 02_BTC_forecast_tails.m for each asset prepared by 09_multiasset_data.R,
% with the input file and output folder overridden. Every model, quantile
% level and tail is estimated exactly as for Bitcoin.
%
% Input : assets/<ASSET>/ret.txt
% Output: assets/<ASSET>/VaR_<tail>_<suf>.txt, ES_<tail>_<suf>.txt,
%         oos_dates.txt, oos_ret.txt   (same format as the BTC files)
%
% Run with:  matlab -batch "eval(fileread('10_multiasset_forecasts.m'))"
%
% NB: 02_BTC_forecast_tails.m defines its own variable `base` and clears
% nothing, so the driver variables carry a trailing underscore.

SCRIPT_PATH_ = [bs_paths(), '02_BTC_forecast_tails.m'];
ROOT_        = [bs_paths(), 'assets', filesep];
if ~exist('ASSETS_','var')
    ASSETS_ = {'BNB','BTC','DOGE','ETH','SOL','TRX','XRP'};
end

for ii_ = 1:numel(ASSETS_)
    IN_FILE = [ROOT_, ASSETS_{ii_}, '/ret.txt'];
    OUT_DIR = [ROOT_, ASSETS_{ii_}, '/'];
    fprintf('\n################ %s ################\n', ASSETS_{ii_});
    t0_ = tic;
    try
        eval(fileread(SCRIPT_PATH_));
    catch err_
        fprintf('!! %s failed: %s\n', ASSETS_{ii_}, err_.message);
    end
    fprintf('%s finished in %.1f s\n', ASSETS_{ii_}, toc(t0_));
    clearvars -except SCRIPT_PATH_ ROOT_ ASSETS_ ii_
end
fprintf('\nMULTI-ASSET FORECASTS DONE\n');
