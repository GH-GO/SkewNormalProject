function [Xsim,garchlags,archlags,AR,MA,Fsim]=ArimaGarch(X)

% This function provides the ARIMA-GARCH 
% (in sample and out of sample) simulations for a given time series

% INPUT: 
% X: time-series 

% OUTPUT: 
% Xsim: in sample simulations from the selected ARIMA-GARCH model
% garchlags, archlags: optimal Garch and Arch lags
% AR, MA: optimal AR and MA lags
% Fsim: out of sample simulations from the selected ARIMA-GARCH model


% Find the optimal garchlags and archlags 

lag = 2;
gloglike = zeros(lag,lag);                                                   
PQ = zeros(lag,lag);
for p = 1:lag
    for q = 1:lag
        garchlagmdl = egarch(p,q);                                               
        [fit,~,logL] = estimate(garchlagmdl,X,'Display','off');         
        gloglike(p,q) = logL;
        PQ(p,q) = p+q;
    end
end
gloglike = reshape(gloglike,lag^2,1);
PQ = reshape(PQ,lag^2,1);
[~,bic] = aicbic(gloglike,PQ+1,100);
ind = find(bic == min(bic));
if ind == 1
	garchlags = 1;
    archlags = 1;
elseif ind == 2
    garchlags = 2;
    archlags = 1;
elseif ind == 3
    garchlags = 1;
    archlags = 2;
elseif ind == 4
    garchlags = 2;
    archlags = 2;
end

% Find the optimal AR and MA lags

lag = 2;
arimaloglike = zeros(lag,lag);
PQ = zeros(lag,lag);
for p = 1:lag
	for q = 1:lag
        mod = arima(p,0,q);
        [fit,~,logL] = estimate(mod,X,'Display','off');
        arimaloglike(p,q) = logL;
        PQ(p,q) = p+q;
     end
end
arimaloglike = reshape(arimaloglike,lag^2,1);           
PQ = reshape(PQ,lag^2,1);                               
[~,bic] = aicbic(arimaloglike,PQ+1,100);
ind = find(bic == min(bic));
if ind == 1
    AR = 1;
    MA = 1;
elseif ind == 2
    AR = 2;
    MA = 1;
elseif ind == 3
    AR = 1;
    MA = 2;
elseif ind == 4
    AR = 2;
    MA = 2;
end

% Simulate the selected ARIMA-GARCH model

model= egarch('GARCHLags',garchlags,'ARCHLags',archlags);
VarMdl = egarch(garchlags,archlags);
Mdl = arima('ARLags',AR,'MALags',MA,'Variance',VarMdl);
EstMdl = estimate(Mdl,X,'Display','off');
[res,v,logL] = infer(EstMdl,X);
% in sample simulations
Xsim=X-res;

L=12;
Fsim(1:L)=NaN;
for i=L:length(X)-1
	f=forecast(EstMdl,1,X(i-L+1:i)); 
	% out of sample simulations
	Fsim(i+1)=f(end);
end
