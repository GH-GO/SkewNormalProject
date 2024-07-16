function MyNRMSE=NRMSE(X,Xfitted)

% This function computes the normalized root mean squared error (NRMSE) statistics
% INPUT:
% X: time series
% Xfitted: simulations of X
% OUTPUT:
% MyNRMSE: NRMSE statistics

format long

R=RMSE(X,Xfitted);
MyNRMSE=R/(max(X)-min(X));