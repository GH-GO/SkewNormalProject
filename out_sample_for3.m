function [F1,F2]=out_sample_for3(par1,par2,i)

% This function provides the out of sample expectations from the proposed model
% INPUT:
% par1=[alpha1,alpha2,alpha3,mu,sigma,a1,b1,a2,b2,rho1,rho2]  estimated parameters array
%%        1         2    3     4   5    6  7  8  9  10   11
% par2=[alpha1,alpha2,alpha3,mu,sigma,a1,b1,a2,b2,rho1,rho2,delta] estimated parameters array
%%        1         2    3     4   5    6  7  8  9  10   11    12
% OUTPUT:
% F1, F2: out of sample expectations in the normal and skew-normal cases
%i=1;

% randn realizations 
N=10^5;
s=rng; % random state 
% fix the random vectors r1,r2, for the sake of reproducibility
   r1=mean(randn(1,N)); r2=mean(randn(1,N)); r3=mean(randn(1,N)); r4=mean(randn(1,N)); r5=mean(randn(1,N));  
rng(s);

% normal case
F1=exp(-par1(1)*par1(4)*i-par1(1)*r1+0.5*par1(1)^2*par1(5)^2*i)*(par1(3)+ par1(2)*par1(6)*i-par1(1)*par1(2)*par1(10)*par1(5)*par1(7)...
   +par1(2)*par1(8)*i-par1(1)*par1(2)*par1(11)*par1(5)*par1(9)...
   +par1(2)*par1(7)*par1(10)*r1+par1(2)*par1(9)*par1(11)*r1+par1(2)*par1(7)*sqrt(1-par1(10)^2)*r2+par1(2)*par1(9)*sqrt(1-par1(11)^2)*r2);


% skew-normal case
A=exp(2*par2(1)*par2(5)*par2(12)*r4)...
  *(1-normcdf(r4+par2(1)*par2(5)*par2(12)))+normcdf(r4-par2(1)*par2(5)*par2(12));


E1=exp(-par2(1)*par2(5)*sqrt(1-par2(12)^2)*r3+par2(1)*par2(5)*par2(12)*r4+0.5*par2(1)^2*par2(5)^2);
E21=exp(0.5*(r4+par2(1)*par2(5)*par2(12))^2-0.5*r4^2); E22=exp(0.5*(r4-par2(1)*par2(5)*par2(12))^2-0.5*r4^2);
E3=exp(-par2(1)*par2(5)*sqrt(1-par2(12)^2)*r3+0.5*(1-par2(12)^2)*par2(1)^2*par2(5)^2);
S1=(par2(1)*par2(5)*par2(12)+r4); S2=(par2(1)*par2(5)*par2(12)-r4);  S3=(r3-par2(1)*par2(5)*sqrt(1-par2(12)^2));
N1=(normcdf(par2(1)*par2(5)*par2(12)+r4)-1); N2=(normcdf(par2(1)*par2(5)*par2(12)-r4)-1);
P1=sqrt(2/pi)*exp(-0.5*r4^2);

F2=exp(-par1(1)*par1(4)*i)*A...
    *( par2(3)*E1...
    +par2(2)*par2(7)*par2(10)*sqrt(1-par2(12)^2)*E1*S3...
    +par2(2)*par2(9)*par2(11)*sqrt(1-par2(12)^2)*E1*S3...
    +par2(2)*par2(6)*i*E1+par2(2)*par2(8)*i*E1...
    +par2(2)*par2(7)*par2(10)*(par2(12)/A)*E3...
    *(S1*N1*E21+S2*N2*E22+P1)...
    +par2(2)*par2(9)*par2(11)*(par2(12)/A)*E3...
    *(S1*N1*E21+S2*N2*E22+P1)...
    +par2(2)*par2(7)*sqrt(1-par2(10)^2)*E1*r4...
    +par2(2)*par2(9)*sqrt(1-par2(11)^2)*E1*r4);