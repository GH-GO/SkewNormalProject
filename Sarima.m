function [Xsim,SAR,SMA,Fsim]=Sarima(X,ARV,MAV,L)

% This function provides the SARIMA 
% (in sample and out of sample) simulations for a given time series

% INPUT: 
% X: time-series 
% ARV, MAV: estimated AR and MA from the ARIMA model 
% L: seasonal lag

% OUTPUT: 
% Xsim: in sample simulations from the selected ARIMA-GARCH model
% SAR, SMA: sarima AR and MA lags
% Fsim: out of sample simulations from the selected SARIMA model

% Find the optimal SAR e SMA 

lag = 2;
arimaloglike = zeros(lag,lag);
PQ = zeros(lag,lag);
for p = 1:lag
    for q = 1:lag
        mod = arima('ARLags',ARV,'SARLags',p,'Seasonality',L,'MALags',MAV,'SMALags',q);
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
    SAR = 1;
    SMA = 1;
elseif ind == 2
    SAR = 2;
    SMA = 1;
elseif ind == 3
    SAR = 1;
    SMA = 2;
elseif ind == 4
    SAR = 2;
    SMA = 2;
end

% Simulate the selected SARIMA model


Mdl = arima('ARLags',ARV,'SARLags',p,'Seasonality',L,'MALags',MAV,'SMALags',q);
EstMdl = estimate(Mdl,X,'Display','off');
[res,v,logL] = infer(EstMdl,X);

% in sample simulations
Xsim=X-res;  
Plus=4;  % SARIMA needs of an extra lag 

% out of sample simulations
Fsim(1:L)=NaN;
for i=L+Plus:length(X)-1
	f=forecast(EstMdl,1,X(i-L-Plus+1:i)); 
	Fsim(i+1)=f(end);
end

