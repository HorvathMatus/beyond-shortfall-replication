function ok = estimate_ok(fin, signrule)
% ESTIMATE_OK  Usability of a parameter estimate, judged at estimation time
% from its IN-SAMPLE fitted [VaR, ES] values only, so that the decision uses
% no out-of-sample information (used by fit_forecast_block.m and
% 14_rolling_forecasts.m).
%
%   fin      - R x 2 in-sample fitted [VaR, ES]
%   signrule - true to also require ES <= VaR < 0 throughout the sample
%              (two-factor GAS model, whose optimiser can diverge)
ok = all(isfinite(fin(:)));
if ok && signrule
    ok = all(fin(:,1)<0) && all(fin(:,2)<=fin(:,1));
end
end
