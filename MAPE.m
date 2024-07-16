function M=MAPE(X,F)

% This function computes the mean absolute percentage error (MAPE) statistics
% INPUT:
% X: time series
% F: simulations of X
% OUTPUT:
% M: MAPE statistics

n=length(X);
s=0;
for i=1:n
    s=s+abs((X(i)-F(i))/X(i));
end
M=s/n;