% Beyond Shortfall - 18: MQR portfolios with the forecasting model chosen per tail
%
% Experiment: ES forecast by the model best suited to the lower tail and EU by
% the model best suited to the upper tail, both chosen at every rebalancing
% date from past forecast losses only (17_tail_selection.R). Three designs:
%   TailSpec   tail-specific schedule  (schedule_tailspec.txt, entries 'A|B')
%   SameModel  one model for both tails, chosen by the same rule
%              (schedule_samemodel.txt) - isolates the value of tail-specificity
%   GARCH-EDF  static ex-ante choice: best model in both tails over Oct 2020 -
%              Jun 2023, before the portfolio evaluation period
% Each for c = 0, 0.1, ..., 1, weekly rebalancing, with the weight cap
% kappa = 0.3 and without a cap. Uses 06_portfolio_weights.m (SCHEDULE option)
% and portfolio_tails.m (pairs 'A|B', GARCH-N, GARCH-EDF).
%
% Outputs: weights_<TailSpec|SameModel|GARCH-EDF>[_cap03]_c<cc>.txt, objval_*.txt
% Run:  matlab -batch "eval(fileread('18_tailspecific_portfolios.m'))"

base_ = bs_paths();
rd_ = @(f) strtrim(splitlines(strtrim(fileread([base_, f]))));
SCH_TS_ = rd_('schedule_tailspec.txt');
SCH_SM_ = rd_('schedule_samemodel.txt');
if ~exist('CAPS_','var'), CAPS_ = [0.3, 1.0]; end
if ~exist('NWORKERS','var'), NWORKERS = 10; end

for ic_ = 1:numel(CAPS_)
    for id_ = 1:3
        WCAP = CAPS_(ic_); CVALS = 0:0.1:1; A_WIN = 7;
        clear SCHEDULE SCHED_TAG MODELS
        switch id_
            case 1, SCHEDULE = SCH_TS_; SCHED_TAG = 'TailSpec';
            case 2, SCHEDULE = SCH_SM_; SCHED_TAG = 'SameModel';
            case 3, MODELS = {'GARCH-EDF'};
        end
        fprintf('\n======== design %d, cap %.1f ========\n', id_, WCAP);
        eval(fileread([base_, '06_portfolio_weights.m']));
    end
end
fprintf('\nTAIL-SPECIFIC PORTFOLIOS DONE\n');
