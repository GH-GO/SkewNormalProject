function [mu0,sigma0,beta0,rho,mu,sigma,beta]=Calibration_SkewGBM(X,start,t,N)

% This function calibrates the parameters:
% - rho, 
% - mu0, sigma0, beta0 (for the dependent variable),
% - mu, sigma, beta (for the regressors) 
% see also cal_Liseo3.m

% INPUT: 
% X: matrix of data. 
% The first column is the dependet variable and the others represent its regressors.
% start, t: inital and final time over which calibrate the parameters  
% N: current time index
% OUTPUT:
% mu0,sigma0,beta0,rho,mu,sigma,beta: estimated (real) parameters 


[n,m]=size(X);
delta_t=t;
% mu0, sigma0, beta0 
parameters0=fitdist_sn0(X(:,1),start,delta_t,N);
mu0=parameters0(1); sigma0=parameters0(2); beta0=parameters0(3);

for i=1:m-1 
	% rho
	rho(i)=corr(X(start:t,1),X(start:t,i+1),'Type','Spearman');
	parametersi=fitdist_sni(X(:,i+1),start,delta_t,beta0,rho(i),N);
	% mu, sigma, beta
	mu(i)=parametersi(1); sigma(i)=parametersi(2); beta(i)=parametersi(3);
end
