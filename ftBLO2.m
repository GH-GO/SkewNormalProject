function fhat = ftBLO2(alpha,z)

% Objective function to be used in calLiseo3.m 
% INPUT: 
% alpha: estimated parameters alpha_i
% z: matrix of regressors
% OUTPUT: objective function for the proposed model

a1 = alpha(1);
a2 = alpha(2);
a3 = alpha(3);

n=length(z);


z1 = z(:,1);
z2 = z(:,2);
z3= z(:,3);


fhat = exp(-a1*z1).*(a2*z2+a3)+z3;