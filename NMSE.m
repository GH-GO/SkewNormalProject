function MyNMSE=NMSE(X,Xfitted)

% This function computes the normalized mean squared error (NMSE) statistics
% INPUT:
% X: time series
% Xfitted: simulations of X
% OUTPUT:
% MyNMSE: NMSE statistics


format long

n=length(X);
s1=0; s2=0;
for i=1:n
    s1 = s1 + (X(i)-Xfitted(i))^2;
	s2=s2+(X(i)-mean(X))^2;
end
MyNMSE=(s1/s2)/n;



    