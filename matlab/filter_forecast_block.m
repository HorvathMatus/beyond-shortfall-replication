function [fo, fin] = filter_forecast_block(model, param, seg, R, alpha)
% FILTER_FORECAST_BLOCK  [VaR, ES] from GIVEN parameters of one of the four
% models estimated with the FZ loss, filtered over seg: the in-sample fitted
% values fin (first R observations) and the forecasts fo for the block that
% follows. Used by 14_rolling_forecasts.m to fall back on the parameters of
% an earlier window; the filter starts at the beginning of the current window
% and is initialised from in-sample information only.
%
%   model - 'FZ2F', 'FZ1F', 'HYB' or 'GFZ'
%   param - parameter vector as stored by fit_forecast_block.m
%           (for 'GFZ': [theta; omega])
seg = seg(:); data2 = seg(1:R);
switch model
    case 'FZ2F'
        [~,fhat] = GAS_twofactor_LL3(param, seg, alpha, VEhatTAU(data2, alpha, -1));
        f = fhat(:,1:2);
    case 'FZ1F'
        [~,fhat] = GAS_onefactor_LL3(param, seg, alpha);
        f = fhat(:,1:2);
    case 'HYB'
        [~,fhat] = GAS_hybrid_LL3(param, seg, alpha, [], [], mean(log(abs(data2))));
        f = fhat(:,1:2);
    case 'GFZ'
        resids2 = data2 - mean(data2);
        [~,fhat] = garch_FZ_LL(param(1:end-1), seg-mean(data2), alpha, [], param(end), cov(resids2));
        f = mean(data2) + fhat(:,1:2);
end
fin = f(1:R,:);
fo  = f(R+1:end,:);
end
