function Parameters=fitdist_sn0(X,start,t,N)

% This function fits the parameters of a skew-Brownian motion through the MLE approach 

% INPUT
% X: time series (dependent variable)
% start, t: inital and final time over which calibrate the parameters  
% N: current time index
% OUTPUT:
% Parameters: estimated (real) parameters

X=X(start:t);

p = @(x,m,s) normpdf(x,m,s);
Phi=@(x,m,s,beta)  normcdf(beta*x,beta*m,s);
f = @(x,m,s,beta) (2/s)*p(x,m,s).*Phi(x,m,s,beta);

opt = statset('MaxIter',1e5,'MaxFunEvals',1e4,'FunValCheck','off');
parameters = mle(X,'pdf',f,'start',[mean(X),std(X),skewness(X)],'Options',opt,'LowerBound',[-Inf,0,-Inf],'UpperBound',[Inf,Inf,Inf]);


Parameters(1)=parameters(1)/N;
Parameters(2)=parameters(2)/sqrt(N);
Parameters(3)=parameters(3);