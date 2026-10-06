% Beyond Shortfall - 14: forecasts with ROLLING RE-ESTIMATION (robustness)
%
% Step 02/10 estimate every model once, on the first R = 1095 observations,
% and filter it forward with fixed parameters. Here the models are
% re-estimated every K_ days on a rolling window of the most recent R
% observations and filtered forward only until the next re-estimation:
%
%   block b covers the forecast dates s_b .. s_b+K_-1, s_b = R+1+(b-1)*K_,
%   and uses the estimation window s_b-R .. s_b-1.
%
% Estimation is done by fit_forecast_block.m, as in step 02. For the four
% models estimated with the FZ loss (FZ-2F, FZ-1F, Hybrid, GARCH-FZ), a window
% without a usable estimate - judged at estimation time from the in-sample
% fitted values only (estimate_ok.m) - falls back on the parameters of the most
% recent earlier window whose in-sample fit over the current window is usable.
% The rolling-window estimators (RW-125/250/500) are unaffected by construction.
%
% Inputs : assets/<A>/ret.txt   (from 09_multiasset_data.R)
% Outputs: assets_roll/<A>/VaR_<tail>_<suf>.txt, ES_<tail>_<suf>.txt
%          (same layout as step 10), oos_dates.txt, oos_ret.txt
%
% Optional overrides (set before running): ASSETS_, K_, OUTROOT_, NWORKERS_
% Run (MATLAB batch, digit-prefixed name):
%   matlab -batch "eval(fileread('14_rolling_forecasts.m'))"

warning off;

base_ = bs_paths();   % code folder; adds the third-party code to the path
if ~exist('ASSETS_','var'),   ASSETS_   = {'BTC','BNB','DOGE','ETH','SOL','TRX','XRP'}; end
if ~exist('K_','var'),        K_        = 91;                         end  % re-estimation interval (days)
if ~exist('OUTROOT_','var'),  OUTROOT_  = [base_,'assets_roll/'];     end
if ~exist('NWORKERS_','var'), NWORKERS_ = 8;                          end
R_    = 1095;
AA_   = [0.01, 0.025, 0.05, 0.10];
SUF_  = {'001','0025','005','01'};
TAILS_ = {'low','up'};

if isempty(gcp('nocreate')), parpool(NWORKERS_); end

for ia_ = 1:numel(ASSETS_)
    a_ = ASSETS_{ia_};
    t_asset_ = tic;
    raw_ = load([base_,'assets/',a_,'/ret.txt'],'ascii');
    raw_ = raw_(raw_(:,2)~=0,:);                  % as in step 02
    T_ = size(raw_,1); P_ = T_ - R_;
    starts_ = R_+1 : K_ : T_;                       % first forecast date of each block
    nb_ = numel(starts_);
    out_ = [OUTROOT_, a_, '/'];
    if ~exist(out_,'dir'), mkdir(out_); end
    writematrix(raw_(R_+1:end,1), [out_,'oos_dates.txt']);
    writematrix(raw_(R_+1:end,2), [out_,'oos_ret.txt']);
    fprintf('%s: T=%d P=%d, %d blocks of %d days\n', a_, T_, P_, nb_, K_);

    for kk_ = 1:numel(TAILS_)
        tail_ = TAILS_{kk_};
        if kk_==1, y_ = raw_(:,2); else, y_ = -raw_(:,2); end

        Vb_ = cell(nb_,1); Eb_ = cell(nb_,1); Pr_ = cell(nb_,1);
        parfor b = 1:nb_
            s = starts_(b); e = min(s+K_-1, T_);
            [Vb_{b}, Eb_{b}, Pr_{b}] = fit_forecast_block(y_(s-R_:e), R_, AA_);
        end

        % Fallback for the four models estimated with the FZ loss: if a
        % window has no usable estimate (estimate_ok.m, judged from in-sample
        % fitted values only), use the parameters of earlier windows, most
        % recent first, whose in-sample fitted values over the CURRENT window
        % are usable. The decision is taken at estimation time; forecasts that
        % turn out to be unusable on a given day are handled in real time in R.
        FB_ = {'FZ2F',7; 'FZ1F',8; 'GFZ',9; 'HYB',10};
        nfb_ = 0; nfail_ = 0;
        for im_ = 1:size(FB_,1)
            mdl_ = FB_{im_,1}; col_ = FB_{im_,2};
            for aa_ = 1:numel(AA_)
                prev_ = {};
                for b = 1:nb_
                    if ~isempty(Pr_{b}.(mdl_){aa_}), prev_{end+1} = Pr_{b}.(mdl_){aa_}; continue; end %#ok<SAGROW>
                    ok_ = false;
                    s = starts_(b); e = min(s+K_-1, T_); seg_ = y_(s-R_:e);
                    for jj_ = numel(prev_):-1:1
                        try
                            [fo_, fin_] = filter_forecast_block(mdl_, prev_{jj_}, seg_, R_, AA_(aa_));
                        catch
                            continue
                        end
                        if estimate_ok(fin_, strcmp(mdl_,'FZ2F'))
                            Vb_{b}(:,col_,aa_) = fo_(:,1); Eb_{b}(:,col_,aa_) = fo_(:,2); ok_ = true;
                            break
                        end
                    end
                    if ok_, nfb_ = nfb_+1;
                    else
                        nfail_ = nfail_+1;
                        fprintf('  !! %s unresolved: %s, alpha=%.3f, block %d (%d..%d)\n', ...
                            mdl_, tail_, AA_(aa_), b, raw_(s,1), raw_(e,1));
                    end
                end
            end
        end

        V_ = cat(1, Vb_{:}); E_ = cat(1, Eb_{:});
        assert(size(V_,1)==P_);
        for aa_ = 1:numel(AA_)
            writematrix(V_(:,:,aa_), [out_,'VaR_',tail_,'_',SUF_{aa_},'.txt']);
            writematrix(E_(:,:,aa_), [out_,'ES_', tail_,'_',SUF_{aa_},'.txt']);
        end
        nanc_ = squeeze(any(isnan(V_),1));          % models x levels with any NaN
        fprintf('  [%s] done; fallbacks %d, unresolved %d; NaN (model,level): %s\n', ...
            tail_, nfb_, nfail_, mat2str(find(nanc_)'));
    end
    fprintf('%s finished in %.1f min\n', a_, toc(t_asset_)/60);
end
fprintf('ROLLING FORECASTS DONE\n');
