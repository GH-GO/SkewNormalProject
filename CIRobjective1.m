function lnL = CIRobjective1(Params,Data)

% Function given from: 
% https://it.mathworks.com/matlabcentral/fileexchange/37297-maximum-likelihood-estimation-of-the-cox-ingersoll-ross-process-the-matlab-implementation?s_tid=prof_contriblnk

% =========================================================================
% PURPOSE : Log-likelihood objective function (multiplied by -1) for the
% CIR process using the MATLAB besseli function
% =========================================================================
% USAGE : Model.TimeStep = Delta t
% Model.Data = Time series of interest rates observations
% Params = Model parameters (alpha, mu, sigma)
% =========================================================================
% RETURNS : lnL = Objective function value
% =========================================================================

% CIR initial parameters estimation
% TimeStep=1/12;
% x = Data(1:end-1); % Time series of interest rates observations
% dx = diff(Data);
% dx = dx./x.^0.5;
% regressors = [TimeStep./x.^0.5, TimeStep*x.^0.5];
% drift = regressors\dx; % OLS regressors coefficients estimates
% res = regressors*drift - dx;
% alpha = -drift(2);
% mu = -drift(1)/drift(2);
% sigma = sqrt(var(res, 1)/TimeStep);
% Params = [alpha mu sigma]; % Vector of initial parameters


DataF = Data(2:end);
DataL = Data(1:end-1);
Nobs = length(Data);
alpha = Params(1);
mu = Params(2);
sigma = Params(3);
TimeStep=1/12;
c = 2*alpha/(sigma^2*(1-exp(-alpha*TimeStep)));
q = 2*alpha*mu/sigma^2-1;
u = c*exp(-alpha*TimeStep)*DataL;
v = c*DataF;
z = 2*sqrt(u.*v);
bf = besseli(q,z,1);
lnL= -(Nobs-1)*log(c) + sum(u + v - 0.5*q*log(v./u) - log(bf) - z);


end


