function Parameters=fitdist_sni(X,start,t,beta0,rho,N)

% This function fits the parameters of a skew-Brownian motion through the MLE approach 

% INPUT:
% X: time series (i-th regressor) 
% start, t: inital and final time over which calibrate the parameters  
% beta0: estimated shape parameter of the dependent variable
% rho: estimated correlation between the i-th regressor and the dependent variable
% N: current time index
% OUTPUT:
% Parameters: estimated (real) parameters

X=X(start:t);

delta0=beta0/sqrt(1+beta0^2);
betai=real(beta0/sqrt(1+(1+beta0^2)*(1-2*delta0^2/pi)*(1-rho^2)/rho^2));


p = @(x,m,s) normpdf(x,m,s);
Phi = @(x,m,s) normcdf(betai*x,betai*m,s);  
f = @(x,m,s) (2/s)*p(x,m,s).*Phi(x,m,s);

opt = statset('MaxIter',1e5,'MaxFunEvals',1e4,'FunValCheck','off','Display','off');
warning('off','all'); 
parameters = mle(X,'pdf',f,'start',[mean(X),std(X)],'Options',opt,'LowerBound',[-Inf,0],'UpperBound',[Inf,Inf]);

Parameters(1)=parameters(1)/N;
Parameters(2)=parameters(2)/sqrt(N);
Parameters(3)=betai;