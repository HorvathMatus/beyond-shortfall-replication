% Beyond Shortfall - 02: OOS VaR/ES forecasts for BTC, lower AND upper tail
%
% Based on Patton, Ziegel and Chen's replication code, extended to:
%   - loop over alpha in {0.01, 0.025, 0.05, 0.1}
%   - loop over both tails:
%       'low' = lower tail (Expected Shortfall), original returns
%       'up'  = upper tail (Expected Upside), NEGATED returns
%     Upper-tail outputs are saved in the NEGATED-return space, i.e. they are
%     the lower-tail VaR/ES of (-r). Flip signs to interpret as upside:
%     VaR_up = -VaR_saved, EU = -ES_saved.
%
% Every model is estimated ONCE on the first R = 1095 observations and filtered
% forward with fixed parameters (the scheme of Patton et al.). The estimation
% itself is in fit_forecast_block.m, which is also used, window by window, by
% the rolling re-estimation of step 14; here the whole evaluation period is a
% single block. All choices made there (starting values, usability of an
% estimate) use in-sample information only.
%
% Data: BTC.txt = date (YYYYMMDD), daily return in PERCENT,
%       2017-10-02 .. 2026-06-30 (3194 obs), in-sample = first 1095 obs.
%
% Outputs (per tail t and alpha suffix s in {001,0025,005,01}):
%   VaR_<t>_<s>.txt, ES_<t>_<s>.txt  - P x 10 matrices, models in columns:
%   RW-125, RW-250, RW-500, GARCH-N, GARCH-skt, GARCH-EDF, FZ2F, FZ1F,
%   GARCH-FZ, Hybrid   (same column order as the original Patton code)
%   oos_dates.txt, oos_ret.txt       - OOS dates and returns (original sign)

base = bs_paths();   % code folder; adds the third-party code to the path
% Optional overrides for the multi-asset run (10_multiasset_forecasts.m):
%   IN_FILE  - return file (YYYYMMDD, percent);  default BTC.txt
%   OUT_DIR  - folder for the forecast files (with trailing /); default base
if ~exist('IN_FILE','var'), IN_FILE = [base,'BTC.txt']; end
if ~exist('OUT_DIR','var'), OUT_DIR = base;            end

raw = load(IN_FILE,'ascii');
raw = raw(raw(:,2)~=0,:);    % dropping zero returns (as in the original code)
T = size(raw,1);
R = 1095;                    % 3-year in-sample period
P = T - R;
fprintf('T=%d  in-sample R=%d  OOS P=%d (%d .. %d)\n', T, R, P, raw(R+1,1), raw(end,1));

writematrix(raw(R+1:end,1), [OUT_DIR,'oos_dates.txt']);
writematrix(raw(R+1:end,2), [OUT_DIR,'oos_ret.txt']);

AA    = [0.01, 0.025, 0.05, 0.10];
SUF   = {'001','0025','005','01'};
TAILS = {'low','up'};
MODEL_NAMES = {'RW-125','RW-250','RW-500','GARCH-N','GARCH-Skt','GARCH-EDF','FZ-2F','FZ-1F','GARCH-FZ','Hybrid'};

for kk = 1:numel(TAILS)
    tail = TAILS{kk};
    if kk==1, y1 = raw(:,2); else, y1 = -raw(:,2); end
    t0 = tic;
    [V, E] = fit_forecast_block(y1, R, AA);       % one block = estimated once
    for aa = 1:numel(AA)
        writematrix(V(:,:,aa), [OUT_DIR,'VaR_',tail,'_',SUF{aa},'.txt']);
        writematrix(E(:,:,aa), [OUT_DIR,'ES_', tail,'_',SUF{aa},'.txt']);
        miss = find(all(isnan(V(:,:,aa)),1));
        if ~isempty(miss)
            fprintf('  !! [%s] alpha=%.3f: no usable estimate for %s\n', tail, AA(aa), strjoin(MODEL_NAMES(miss), ', '));
        end
    end
    fprintf('[%s] done in %.1f s\n', tail, toc(t0));
end

fprintf('ALL DONE\n');
