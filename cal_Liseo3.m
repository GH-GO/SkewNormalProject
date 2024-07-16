function [param,MSE]=cal_Liseo3(X,V)

% This function calibrates the parametrs alpha_i (i=1,2,3,4) for our model 
% see also Calibration_SkewGBM.m 

% INPUT: 
% V: time series denoting the dependent variable 
% X: matrix whose columns represent the regeressors of X

% OUTPUT
% param: estimated parameters alpha_i
% MSE: mean squared error produced from the calibration 


n=length(X(:,1));
Z1=X(:,1);  Z2=X(:,2)+X(:,3); 
Z3=[0;V(1:end-1)]; 
Z=[Z1,Z2,Z3];

% initialization 
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
[alpha,Res,J,Cov,MSE]=nlinfit(Z(end,:),V(end),@ftBLO2,[0,0,0],opts);
param(1)=alpha(1);
param(2)=alpha(2);
param(3)=alpha(3);

% starting parameters
param0=alpha;

% find the optimal parameters
syms a b c 
eq1=0; eq2=0; eq3=0;
 
for i=n:n   
	eq1=eq1+Z1(i)*exp(-a*Z1(i))*(b*Z2(i)+c)*(V(i)-exp(-a*Z1(i))*(b*Z2(i)+c)-Z3(i));
	eq2=eq2+Z2(i)*exp(-a*Z1(i))*(V(i)-exp(-a*Z1(i))*(b*Z2(i)+c)-Z3(i));
	eq3=eq3+exp(-a*Z1(i))*(V(i)-exp(-a*Z1(i))*(b*Z2(i)+c)-Z3(i));
end

eqns=[eq1==0 eq2==0 eq3==0];
S=vpasolve(eqns,[a b c],param0);
% alpha_i (i=1,2,3)
param(1)=double(S.a);
param(2)=double(S.b);
param(3)=double(S.c);

tau=0;
for i=n-1:n
	tau=tau+(V(i)-exp(-param(1)*Z1(i))*(param(2)*Z2(i)+param(3))-Z3(i))^2;
end
% alpha_4
param(4)=tau/n;