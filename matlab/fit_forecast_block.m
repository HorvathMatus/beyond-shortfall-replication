function [VaR, ES, prm] = fit_forecast_block(seg, R, AA)
% FIT_FORECAST_BLOCK  Estimate the ten models on one window and forecast the
% block that follows it (used by 02_BTC_forecast_tails.m with a single block
% and by 14_rolling_forecasts.m with blocks of 91 days).
%
%   seg - returns of ONE tail (already negated for the upper tail): the R
%         estimation observations followed by the Pb forecast-block observations
%   R   - length of the estimation window
%   AA  - quantile levels
%
%   VaR, ES - Pb x 10 x numel(AA) one-step-ahead forecasts for the block, models
%             in the column order of Patton's code: RW-125, RW-250, RW-500,
%             GARCH-N, GARCH-skt, GARCH-EDF, FZ2F, FZ1F, GARCH-FZ, Hybrid
%   prm     - struct with the parameters of the four models estimated with the
%             FZ loss (fields FZ2F, FZ1F, HYB, GFZ; cell per level, [] if no
%             usable estimate) so that the rolling driver can fall back on an
%             earlier window's estimate
%
% The estimation code is that of Patton, Ziegel and Chen's replication files.
% Every decision is taken at estimation time from in-sample information only
% (estimate_ok.m): an estimate is usable if its in-sample fitted VaR and ES are
% finite, and for the two-factor GAS model, whose optimiser can diverge, the
% starting value is chosen as the best in-sample fit among those whose fitted
% values satisfy ES <= VaR < 0. The two-factor filter is initialised from the
% estimation sample (Patton's code initialises it from the series passed in,
% which would include the forecast period). Forecasts that are unusable on a
% given day are handled in real time in R (R/forecast_eval.R).

warning off;
options = optimset('Display','off','MaxFunEvals',2000);

T  = numel(seg);
Pb = T - R;
y1 = seg(:);
data2 = y1(1:R);
data3 = y1(R+1:end);
nA = numel(AA);
VaR = nan(Pb,10,nA); ES = nan(Pb,10,nA);
prm.FZ2F = cell(nA,1); prm.FZ1F = cell(nA,1); prm.HYB = cell(nA,1); prm.GFZ = cell(nA,1);

% ---- mean and GARCH estimation (level-independent) ----------------------
resids2 = data2 - mean(data2);
[parameters, ~, hhat2] = garchpq(resids2, 1, 1);
GARCHparamALL = nan(5,1);
GARCHparamALL(1:3,1) = parameters';
lower = [2.1,-0.99]'; upper = [50,0.99]'; theta0 = [5,-0.1]';
skewtparams = fmincon('skewtdis_LL',theta0,[],[],[],[],lower,upper,[],options,resids2./sqrt(hhat2));
GARCHparamALL(4:5,1) = skewtparams;

muhatOOS = mean(data2)*ones(Pb,1);
hhatOOS  = nan(Pb,1);
hhatOOS(1) = GARCHparamALL(1:3,1)'*[1;(data2(end)-mean(data2))^2;hhat2(end)];
for tt = 2:Pb
    hhatOOS(tt) = GARCHparamALL(1:3,1)'*[1;(data3(tt-1)-muhatOOS(tt-1))^2;hhatOOS(tt-1)];
end

for aa = 1:nA
    alpha = AA(aa);
    VEhatALL = nan(Pb,10,2);

    % ---- rolling windows ------------------------------------------------
    MM = [125;250;500];
    for mm = 1:length(MM)
        M = MM(mm);
        for tt = R+1:T
            VEhatALL(tt-R,mm,:) = sample_VE(y1(tt-M:tt-1),alpha);
        end
    end

    % ---- GARCH-N / GARCH-skt / GARCH-EDF ----------------------------------
    VEhatALL(:,4,:) = muhatOOS + sqrt(hhatOOS)*[norminv(alpha),-normpdf(norminv(alpha))/alpha];
    VEhatALL(:,5,:) = muhatOOS + sqrt(hhatOOS)*skewt_VE(GARCHparamALL(4:5),alpha);
    VEhatALL(:,6,:) = muhatOOS + sqrt(hhatOOS)*sample_VE(resids2./sqrt(hhat2),alpha);

    % ---- GARCH-FZ ---------------------------------------------------------
    try
        VEbar  = sample_VE( resids2./sqrt(hhat2) , alpha)';
        theta0 = [norminv(GARCHparamALL(3)); log(GARCHparamALL(2)); log(-VEbar(2)); norminv(VEbar(1)/VEbar(2))];
        [thetahat1,LL1] = fminsearch('garch_FZ_LL',theta0,options,resids2,alpha,[],GARCHparamALL(1));
        [thetahat2,~  ] = fminunc(  'garch_FZ_LL',theta0,   options,resids2,alpha, 5,GARCHparamALL(1));
        [thetahat3,~  ] = fminunc(  'garch_FZ_LL',thetahat2,options,resids2,alpha,20,GARCHparamALL(1));
        [thetahat4,LL4] = fminsearch('garch_FZ_LL',thetahat3,options,resids2,alpha,[],GARCHparamALL(1));
        if LL1<LL4, garchFZparam = thetahat1; else, garchFZparam = thetahat4; end
        p = [garchFZparam; GARCHparamALL(1)];
        [fo, fin] = filter_forecast_block('GFZ', p, y1, R, alpha);
        if estimate_ok(fin, false)
            VEhatALL(:,9,:) = fo; prm.GFZ{aa} = p;
        end
    catch
    end

    % ---- one-factor GAS (FZ1F) ----------------------------------------------
    try
        VEbar  = sample_VE( data2 , alpha)';
        theta0 = [norminv(0.95); log(0.005); log(-VEbar(2)); norminv(VEbar(1)/VEbar(2))];
        [thetahat2,~  ] = fminunc(  'GAS_onefactor_LL3',theta0,   options,data2,alpha, 5);
        [thetahat3,~  ] = fminunc(  'GAS_onefactor_LL3',thetahat2,options,data2,alpha,20);
        [thetahat4,~  ] = fminsearch('GAS_onefactor_LL3',thetahat3,options,data2,alpha,[]);
        FZ1Fparam = thetahat4;
        [fo, fin] = filter_forecast_block('FZ1F', FZ1Fparam, y1, R, alpha);
        if estimate_ok(fin, false)
            VEhatALL(:,8,:) = fo; prm.FZ1F{aa} = FZ1Fparam;
        end
    catch
        FZ1Fparam = [];
    end

    % ---- Hybrid GAS/GARCH -------------------------------------------------
    try
        theta0 = [norminv(0.95); log(0.005); log(0.005); log(-VEbar(2)); norminv(VEbar(1)/VEbar(2))];
        [thetahat4,LL4] = fminsearch('GAS_hybrid_LL3', theta0, options, data2, alpha, []);
        if ~isempty(FZ1Fparam)
            [thetahat5,LL5] = fminsearch('GAS_hybrid_LL3', [FZ1Fparam(1:2);-10;FZ1Fparam(3:4)], options, data2, alpha, []);
        else
            LL5 = Inf;
        end
        if LL4<LL5, hybridFZparam = thetahat4; else, hybridFZparam = thetahat5; end
        [fo, fin] = filter_forecast_block('HYB', hybridFZparam, y1, R, alpha);
        if estimate_ok(fin, false)
            VEhatALL(:,10,:) = fo; prm.HYB{aa} = hybridFZparam;
        end
    catch
    end

    % ---- two-factor GAS (FZ2F) ----------------------------------------------
    try
        theta00 = [...
            [ (1-0.97)*VEbar(1) ; (1-0.97)*VEbar(2) ; norminv(0.97) ; norminv(0.97) ; 0 ; 0 ; 0 ; 0 ] , ...
            [ (1-0.98)*VEbar(1) ; (1-0.98)*VEbar(2) ; norminv(0.98) ; norminv(0.98) ; 0 ; 0 ; 0 ; 0 ] , ...
            [ (1-0.99)*VEbar(1) ; (1-0.99)*VEbar(2) ; norminv(0.99) ; norminv(0.99) ; 0 ; 0 ; 0 ; 0 ] ];
        outSTART = nan(size(theta00,1)+1,size(theta00,2));
        for ss = 1:size(theta00,2)
            theta0 = theta00(:,ss);
            temp = GAS_twofactor_LL3(theta0,data2,alpha,[],5);
            if temp<1e6 && ~isnan(temp) && ~isinf(temp)
                [thetahat2,~] = fminunc('GAS_twofactor_LL3',theta0,options,data2,alpha,[],5);
            else
                thetahat2 = theta0;
            end
            temp = GAS_twofactor_LL3(thetahat2,data2,alpha,[],20);
            if temp<1e6 && ~isnan(temp) && ~isinf(temp)
                thetahat2a = thetahat2;
            else
                thetahat2a = theta0;
            end
            [thetahat3,~] = fminunc('GAS_twofactor_LL3',thetahat2a,options,data2,alpha,[],20);
            temp = GAS_twofactor_LL3(thetahat3,data2,alpha,[],[]);
            if temp<1e6 && ~isnan(temp) && ~isinf(temp)
                thetahat3a = thetahat3;
            else
                thetahat3a = theta0;
            end
            [thetahat4,LL4] = fminsearch('GAS_twofactor_LL3',thetahat3a,options,data2,alpha,[],[]);
            outSTART(:,ss) = [thetahat4;LL4];
        end
        % best in-sample fit among the starting values whose IN-SAMPLE fitted
        % values satisfy ES <= VaR < 0 (decided at estimation time)
        [~,ordLL] = sort(outSTART(end,:),'ascend');
        for ss = ordLL
            if ~isfinite(outSTART(end,ss)), continue; end
            cand = outSTART(1:end-1,ss);
            [fo, fin] = filter_forecast_block('FZ2F', cand, y1, R, alpha);
            if estimate_ok(fin, true)
                prm.FZ2F{aa} = cand;
                VEhatALL(:,7,:) = fo;
                break
            end
        end
    catch
    end

    VaR(:,:,aa) = VEhatALL(:,:,1);
    ES(:,:,aa)  = VEhatALL(:,:,2);
end
end
