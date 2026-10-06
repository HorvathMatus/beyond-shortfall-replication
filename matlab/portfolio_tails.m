function [esmag, eumag] = portfolio_tails(rp, alpha, model, h, which)
% Conditional lower- and upper-tail risk of a portfolio return series.
%
%  rp    : W x 1 portfolio returns (percent) over the estimation window
%  alpha : tail probability (e.g. 0.05)
%  model : 'RW125' | 'GARCH-N' | 'GARCH-Skt' | 'GARCH-EDF' | 'GARCH-FZ' | 'Hybrid'
%          or a pair 'A|B': ES forecast by model A, EU forecast by model B
%          (tail-specific choice of the forecasting model, step 18)
%  h     : length of the forecast block at the end of the window
%  which : 'both' (default), 'low' or 'up' - compute only the tails needed
%
%  esmag : POSITIVE magnitude of Expected Shortfall  (lower tail of rp)
%  eumag : POSITIVE magnitude of Expected Upside     (upper tail of rp,
%          obtained as the lower tail of -rp)
%
% The model is fitted on the first W-h observations and filtered over the
% last h, exactly as in the portfolio_*.m routines of the previous paper;
% the reported figure is the average forecast over that block.

if nargin < 4 || isempty(h), h = 7; end
if nargin < 5 || isempty(which), which = 'both'; end
doLo = ~strcmp(which, 'up');
doUp = ~strcmp(which, 'low');

% ---- tail-specific pair 'A|B' -------------------------------------------
bar = strfind(model, '|');
if ~isempty(bar)
    A = model(1:bar-1); B = model(bar+1:end);
    esmag = NaN; eumag = NaN;
    if doLo, [esmag, ~] = portfolio_tails(rp, alpha, A, h, 'low'); end
    if doUp, [~, eumag] = portfolio_tails(rp, alpha, B, h, 'up');  end
    return
end

rp = rp(:);
W  = numel(rp);
R  = W - h;
opts = optimset('Display','off','MaxFunEvals',2000);
esmag = NaN; eumag = NaN;

switch model

case 'RW125'
    M = 125;
    lo = nan(h,1); up = nan(h,1);
    for tt = R+1:W
        if doLo, a = sample_VE( rp(tt-M:tt-1), alpha);  lo(tt-R) = a(2); end
        if doUp, b = sample_VE(-rp(tt-M:tt-1), alpha);  up(tt-R) = b(2); end
    end
    esmag = -mean(lo);  eumag = -mean(up);

case {'GARCH-Skt', 'GARCH-N', 'GARCH-EDF'}
    % one GARCH(1,1) fit serves both tails (the variance process is
    % sign-symmetric); only the innovation tail differs between the tails
    din = rp(1:R);  dout = rp(R+1:W);
    mu  = mean(din);
    e   = din - mu;
    [par, ~, hh] = garchpq(e, 1, 1);
    z   = e ./ sqrt(hh);
    p3 = par(:)';                       % garchpq returns a column vector
    hO = nan(h,1);
    hO(1) = p3 * [1; (din(end)-mu)^2; hh(end)];
    for tt = 2:h
        hO(tt) = p3 * [1; (dout(tt-1)-mu)^2; hO(tt-1)];
    end
    switch model
        case 'GARCH-Skt'
            lb = [2.1;-0.99]; ub = [50;0.99]; th0 = [5;-0.1];
            if doLo
                sk_lo = fmincon('skewtdis_LL', th0, [],[],[],[], lb, ub, [], opts,  z);
                ve_lo =  mu + sqrt(hO) * skewt_VE(sk_lo, alpha);  esmag = -mean(ve_lo(:,2));
            end
            if doUp
                sk_up = fmincon('skewtdis_LL', th0, [],[],[],[], lb, ub, [], opts, -z);
                ve_up = -mu + sqrt(hO) * skewt_VE(sk_up, alpha);  eumag = -mean(ve_up(:,2));
            end
        case 'GARCH-N'
            qn = [norminv(alpha), -normpdf(norminv(alpha))/alpha];
            if doLo, ve_lo =  mu + sqrt(hO) * qn;  esmag = -mean(ve_lo(:,2)); end
            if doUp, ve_up = -mu + sqrt(hO) * qn;  eumag = -mean(ve_up(:,2)); end
        case 'GARCH-EDF'
            if doLo, ve_lo =  mu + sqrt(hO) * sample_VE( z, alpha);  esmag = -mean(ve_lo(:,2)); end
            if doUp, ve_up = -mu + sqrt(hO) * sample_VE(-z, alpha);  eumag = -mean(ve_up(:,2)); end
    end

case 'GARCH-FZ'
    if doLo, esmag = -mean(garchfz_block( rp, alpha, R, h, opts)); end
    if doUp, eumag = -mean(garchfz_block(-rp, alpha, R, h, opts)); end

case 'Hybrid'
    if doLo, esmag = -mean(hybrid_block( rp, alpha, R, h, opts)); end
    if doUp, eumag = -mean(hybrid_block(-rp, alpha, R, h, opts)); end

otherwise
    error('unknown model %s', model);
end

if ~isfinite(esmag), esmag = 1e6; end   % keep fmincon away from failures
if ~isfinite(eumag), eumag = 0;   end
end

% ---------------------------------------------------------------
function es = garchfz_block(y, alpha, R, h, opts)
y   = y(:);  din = y(1:R);  mu = mean(din);  e = din - mu;
[par, ~, hh] = garchpq(e, 1, 1);
VEbar = sample_VE(e ./ sqrt(hh), alpha)';
th0   = [norminv(par(3)); log(par(2)); log(-VEbar(2)); norminv(VEbar(1)/VEbar(2))];
[t1, L1] = fminsearch('garch_FZ_LL', th0, opts, e, alpha, [], par(1));
[t2, ~ ] = fminunc(  'garch_FZ_LL', th0, opts, e, alpha,  5, par(1));
[t3, ~ ] = fminunc(  'garch_FZ_LL', t2,  opts, e, alpha, 20, par(1));
[t4, L4] = fminsearch('garch_FZ_LL', t3,  opts, e, alpha, [], par(1));
if L1 < L4, th = t1; else, th = t4; end
[~, fhat] = garch_FZ_LL(th, y - mu, alpha, [], par(1), cov(e));
es = mu + fhat(R+1:end, 2);
end

% ---------------------------------------------------------------
function es = hybrid_block(y, alpha, R, h, opts)
y   = y(:);  din = y(1:R);
VEbar = sample_VE(din, alpha)';
th0 = [norminv(0.95); log(0.005); log(-VEbar(2)); norminv(VEbar(1)/VEbar(2))];
[a2, ~] = fminunc(  'GAS_onefactor_LL3', th0, opts, din, alpha,  5);
[a3, ~] = fminunc(  'GAS_onefactor_LL3', a2,  opts, din, alpha, 20);
[a4, ~] = fminsearch('GAS_onefactor_LL3', a3, opts, din, alpha, []);
th0h = [norminv(0.95); log(0.005); log(0.005); log(-VEbar(2)); norminv(VEbar(1)/VEbar(2))];
[b4, L4] = fminsearch('GAS_hybrid_LL3', th0h, opts, din, alpha, []);
[b5, L5] = fminsearch('GAS_hybrid_LL3', [a4(1:2); -10; a4(3:4)], opts, din, alpha, []);
if L4 < L5, th = b4; else, th = b5; end
[~, fhat] = GAS_hybrid_LL3(th, y, alpha, [], [], mean(log(abs(din(din~=0)))));
es = fhat(R+1:end, 2);
end
