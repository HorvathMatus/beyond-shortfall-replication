% Beyond Shortfall - 06: MQR portfolio weights
%
% Objective (thesis eq. 5.1), for each rebalancing date:
%     min_w   ES(w) - c * EU(w)      s.t.  sum(w)=1,  0 <= w <= 1
% where ES and EU are conditional FORECASTS of the lower/upper tail of the
% portfolio return series, produced by one of the models compared in the MCS.
%   c = 0   pure minimum-ES portfolio (previous paper)
%   c = 1   downside and upside weighted equally
%
% Rolling scheme: 1095-day (3y) estimation window, rebalanced every 7 days,
% alpha = 0.05, matching the previous paper.
%
% Input  : coin_ret.txt (from 05_portfolio_data.R), T x N returns in percent
% Output : weights_<model>_c<cc>.txt   (n_periods x N)
%          objval_<model>_c<cc>.txt    (n_periods x 3: obj, ES, EU)
%
% Run with:  matlab -batch "eval(fileread('06_portfolio_weights.m'))"

base   = bs_paths();   % code folder; adds the third-party code to the path
data   = load([base,'coin_ret.txt'],'ascii');
[T, N] = size(data);

W_EST  = 1095;      % estimation window
if ~exist('A_WIN','var'), A_WIN = 7; end   % rebalancing frequency (days)
H_FC   = 7;         % forecast block at the end of the window
ALPHA  = 0.05;
NPER   = floor((T - W_EST) / A_WIN);

if ~exist('MODELS','var'), MODELS = {'RW125','GARCH-Skt'}; end
if ~exist('CVALS','var'),  CVALS  = 0:0.1:1;               end
if ~exist('WCAP','var'),   WCAP   = 1.0;                   end  % max weight per asset

% Optional model schedule (step 18): SCHEDULE is a cell array with one model
% string per rebalancing period, e.g. 'GARCH-Skt|GARCH-FZ' (ES from the first
% model, EU from the second), and SCHED_TAG names the output files. When it is
% set, MODELS is replaced by the single entry SCHED_TAG.
useSched = exist('SCHEDULE','var') && ~isempty(SCHEDULE);
if useSched
    assert(numel(SCHEDULE) == NPER, 'SCHEDULE must have one entry per rebalancing period');
    MODELS = {SCHED_TAG};
else
    SCHEDULE = cell(NPER, 1);   % sliced by the parfor loop, so it needs NPER entries; not used
end

fprintf('T=%d N=%d  window=%d  step=%d  periods=%d\n', T, N, W_EST, A_WIN, NPER);
fprintf('models: %s\n', strjoin(MODELS,', '));
fprintf('c values: %s\n\n', mat2str(CVALS));

if ~exist('NWORKERS','var'), NWORKERS = 8; end
pool = gcp('nocreate');
if isempty(pool), pool = parpool('local', min(NWORKERS, feature('numcores'))); end
fprintf('parallel workers: %d\n\n', pool.NumWorkers);

Aeq = ones(1,N); beq = 1; lb = zeros(N,1); ub = WCAP*ones(N,1);

% Vertices of {sum(w)=1, 0<=w<=WCAP}: kv = ceil(1/WCAP) assets, all but one
% at the cap. WCAP=1 gives the single-asset corners, WCAP=0.5 all pairs at
% (0.5,0.5). The count grows as nchoosek(N,kv), so for tighter caps full
% enumeration becomes the dominant cost; beyond VERT_MAX we fall back to a
% greedy vertex built from the ranking of the assets evaluated individually.
VERT_MAX = 40;
kv    = ceil(1/WCAP);
combs = nchoosek(1:N, kv);
if size(combs,1) <= VERT_MAX
    VERT = zeros(N, size(combs,1));
    for k = 1:size(combs,1)
        v = zeros(N,1);
        v(combs(k,1:end-1)) = WCAP;
        v(combs(k,end))     = 1 - WCAP*(kv-1);
        VERT(:,k) = v;
    end
    GREEDY = false;
    fprintf('weight cap = %.2f -> %d vertices scanned per period\n\n', WCAP, size(VERT,2));
else
    VERT = zeros(N,0);
    GREEDY = true;
    fprintf(['weight cap = %.2f -> %d vertices would be needed; using a greedy\n' ...
             'start from the single-asset ranking instead (%d evaluations)\n\n'], ...
             WCAP, size(combs,1), N+1);
end
opts = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
        'MaxFunctionEvaluations',300,'OptimalityTolerance',1e-6, ...
        'StepTolerance',1e-8,'ConstraintTolerance',1e-8, ...
        'FiniteDifferenceStepSize',1e-3);   % coarse step: objective is a
                                            % numerically-refitted model, so
                                            % very small steps are pure noise

for mi = 1:numel(MODELS)
  model = MODELS{mi};
  for ci = 1:numel(CVALS)
    c   = CVALS(ci);
    % file tag: <model>[_cap<kk>][_r<freq>]_c<cc>.  The 7-day baseline carries
    % no frequency suffix, so the files already computed remain valid.
    tag = model;
    if WCAP  <  1, tag = sprintf('%s_cap%02d', tag, round(WCAP*10)); end
    if A_WIN ~= 7, tag = sprintf('%s_r%02d',   tag, A_WIN);          end
    tag = sprintf('%s_c%02d', tag, round(c*10));
    fw  = [base,'weights_',tag,'.txt'];
    if exist(fw,'file'), fprintf('[%s] exists, skipping\n', tag); continue; end

    Wopt = zeros(NPER, N);  Ovals = zeros(NPER, 3);
    t0 = tic;
    parfor t = 1:NPER
      win = data((t-1)*A_WIN+1 : W_EST+(t-1)*A_WIN, :);
      if useSched, model_t = SCHEDULE{t}; else, model_t = model; end
      objf = @(w) portfolio_objective(w, win, ALPHA, model_t, H_FC, c);

      % starting points: 1/N and every vertex of the feasible set. The
      % objective is often minimised at a vertex (ES and EU both scale with
      % volatility), so scanning them first avoids poor local solutions.
      % Under a tight cap the vertex set is too large to enumerate, so the
      % assets are instead ranked by their standalone objective (evaluated
      % without the cap, purely as a ranking device) and the cap is filled
      % greedily from the best of them.
      if GREEDY
          fs = inf(1,N);
          for j = 1:N
              e = zeros(N,1); e(j) = 1;
              try, fs(j) = objf(e); catch, fs(j) = inf; end
          end
          [~, ord] = sort(fs, 'ascend');
          g = zeros(N,1); rem = 1; p = 1;
          while rem > 1e-12 && p <= N
              g(ord(p)) = min(WCAP, rem); rem = rem - g(ord(p)); p = p + 1;
          end
          cands = [ones(N,1)/N, g];
      else
          cands = [ones(N,1)/N, VERT];
      end
      fc = inf(1, size(cands,2));
      for k = 1:size(cands,2)
        try, fc(k) = objf(cands(:,k)); catch, fc(k) = inf; end
      end
      [~, kb] = min(fc);
      w0 = cands(:,kb);

      try
        [wo, fv] = fmincon(objf, w0, [],[], Aeq, beq, lb, ub, [], opts);
        if ~isfinite(fv) || fv > fc(kb), wo = w0; fv = fc(kb); end
      catch
        wo = w0; fv = fc(kb);
      end
      wo = max(wo,0); wo = wo / sum(wo);

      [es, eu] = portfolio_tails(win*wo, ALPHA, model_t, H_FC);
      Wopt(t,:)  = wo';
      Ovals(t,:) = [fv, es, eu];
    end
    writematrix(Wopt,  fw);
    writematrix(Ovals, [base,'objval_',tag,'.txt']);
    fprintf('[%s] %d periods in %.1f min  (mean ES=%.3f EU=%.3f)\n', ...
            tag, NPER, toc(t0)/60, mean(Ovals(:,2)), mean(Ovals(:,3)));
  end
end
fprintf('\nPORTFOLIO WEIGHTS DONE\n');
