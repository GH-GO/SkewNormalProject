function [F1,F2]=in_sample_for3(par1,par2,i)

% This function provides the in sample expectations from the proposed model
% INPUT:
% par1=[alpha1,alpha2,alpha3,mu,sigma,a1,b1,a2,b2,rho1,rho2]  estimated parameters array
%%        1         2    3     4   5    6  7  8  9  10   11
% par2=[alpha1,alpha2,alpha3,mu,sigma,a1,b1,a2,b2,rho1,rho2,delta] estimated parameters array
%%        1         2    3     4   5    6  7  8  9  10   11    12
% OUTPUT:
% F1, F2: in sample expectations in the normal and skew-normal cases

  
% normal case       
F1=exp(-par1(1)*par1(4)+0.5*par1(1)^2*par1(5)^2*i)*(par1(3)+ par1(2)*par1(6)*i-par1(1)*par1(2)*par1(10)*par1(5)*par1(7)*i...
   +par1(2)*par1(8)*i-par1(1)*par1(2)*par1(11)*par1(5)*par1(9)*i );

% skew-normal case
F2=exp(-par1(1)*par1(4)*i)...
   *(2*par2(3)*exp(0.5*par2(1)^2*par2(5)^2*i)*(1-normcdf(par2(1)*par2(5)*par2(12)*sqrt(i)))...
   +2*par2(2)*exp(0.5*par2(1)^2*par2(5)^2*i)*par2(6)*i*(1-normcdf(par2(1)*par2(5)*par2(12)*sqrt(i)))...
   +2*par2(2)*exp(0.5*par2(1)^2*par2(5)^2*i)*par2(8)*i*(1-normcdf(par2(1)*par2(5)*par2(12)*sqrt(i)))...
   -2*sqrt(1-par2(12)^2)*exp(0.5*par2(1)^2*par2(5)^2*i)*i*par2(1)*par2(2)*par2(5)*par2(7)*par2(10)*(1-normcdf(par2(1)*par2(5)*par2(12)*sqrt(i)))...
   -2*sqrt(1-par2(12)^2)*exp(0.5*par2(1)^2*par2(5)^2*i)*i*par2(1)*par2(2)*par2(5)*par2(9)*par2(11)*(1-normcdf(par2(1)*par2(5)*par2(12)))...
   -par2(2)*par2(7)*par2(10)*par2(12)*(2*par2(1)*par2(12)*par2(5)*exp(0.5*par2(1)^2*par2(5)^2)*(1-normcdf(par2(1)*par2(5)*par2(12)*sqrt(i)))-sqrt(2*i/pi) )...
   -par2(2)*par2(9)*par2(11)*par2(12)*(2*par2(1)*par2(12)*par2(5)*exp(0.5*par2(1)^2*par2(5)^2)*(1-normcdf(par2(1)*par2(5)*par2(12)*sqrt(i)))-sqrt(2*i/pi) ));

    
