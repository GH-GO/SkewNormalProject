function MyRMSE=RMSE(X,Xfitted)

% This function computes the root mean squared error (RMSE) statistics
% INPUT:
% X: time series
% Xfitted: simulations of X
% OUTPUT:
% MyRMSE: RMSE statistics

format long

n=length(X);
s=0;
for i=1:n
    s = s + (X(i)-Xfitted(i))^2;
end
MyMop=s./n;

MyRMSE=sqrt(MyMop);

    