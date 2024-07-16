function Results = CIRestimation(Data,TimeStep,Method) 

% Function given from: 
% https://it.mathworks.com/matlabcentral/fileexchange/37297-maximum-likelihood-estimation-of-the-cox-ingersoll-ross-process-the-matlab-implementation?s_tid=prof_contriblnk

% =========================================================================
% PURPOSE : CIR process maximum likelihood estimation
%            1) OLS initial parameters estimation
%            2) Log-likelihood function optimization
%            3) Printing Results   
% =========================================================================
% USAGE   : Data       = Time series of interest rates observations
%         : TimeStep   = Delta t; recommended: 1/250 for daily data
%                                                    1/12 for monthly data                                                  
%         
%         : Method     = 'ncx2pdf' | 'besseli' (default: 'besseli')
%         
% =========================================================================
% RETURNS : Results.Params   = Estimated parameters (alpha, mu, sigma)
%         : Results.Fval     = Objective function value
% =========================================================================
% Kamil Kladivko; Technical Computing Prague 2007
% Date: October 2007 
% Questions: kladivko@gmail.com

% Initial parameters using OLS
Nobs = length(Data);

x = Data(1:end-1);
dx = diff(Data);           
dx = dx./x.^0.5;
regressors = [TimeStep./x.^0.5, TimeStep*x.^0.5];
drift = regressors\dx;
res = regressors*drift - dx;
alpha = -drift(2);
mu = -drift(1)/drift(2);
sigma = sqrt(var(res, 1)/TimeStep);
InitialParams = [alpha mu sigma];

	fprintf('\n initial alpha = %+3.6f\n initial mu    = %+3.6f\n initial sigma = %+3.6f\n', alpha, mu, sigma);


% Optimization using fminsearch

options = optimset('LargeScale', 'off', 'MaxIter', 300, 'MaxFunEvals', 300, 'TolFun', 1e-4, 'TolX', 1e-4, 'TolCon', 1e-4); 

if Method==1
	[Params, Fval, Exitflag, iterations] =  fminsearch(@(Params) CIRobjective1(Params, Data), InitialParams, options)   
else 
	if Method==2
		[Params, Fval, Exitflag,iterations] =  fminsearch(@(Params) CIRobjective2(Params, Data), InitialParams, options)   
    end
end

Results.Params = Params;
Results.Fval = -Fval/Nobs;
Results.Exitflag = Exitflag;


fprintf('\n alpha = %+3.6f\n mu    = %+3.6f\n sigma = %+3.6f\n', Params(1), Params(2), Params(3));
fprintf(' log-likelihood = %+3.6f\n', -Fval/Nobs);
    
end
