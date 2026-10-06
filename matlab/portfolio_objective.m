function v = portfolio_objective(w, data, alpha, model, h, c)
% Objective of the MQR (maximising quantile returns) portfolio.
%
% Thesis eq. (5.1):   max_w  E[y_p | y_p < q_lower] + c * E[y_p | y_p > q_upper]
% In signed return space the first term is negative (a loss) and the second
% positive (a gain). Writing both as POSITIVE magnitudes, ES and EU, that
% maximisation is the minimisation returned here:
%
%        v(w) = ES(w) - c * EU(w)
%
%   c = 0 : pure minimum-ES portfolio (the previous paper)
%   c = 1 : downside and upside weighted equally (thesis baseline)
%
% ES and EU are conditional FORECASTS from `model`, not historical averages.

[es, eu] = portfolio_tails(data * w(:), alpha, model, h);
v = es - c * eu;
end
