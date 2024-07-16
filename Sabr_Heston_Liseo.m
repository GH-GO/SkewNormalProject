% Sabr_Heston_Liseo

% This procedure performs the SABR-HESTON benchmark for the time series Comm_volBrent 

clear all
close all
clc


% load time series
b1=Comm_volBrent;
m=12;
n=length(b1);
N=10^5;


% pointwise volatility of b1
sb1=point_vol(b1,m);
Delta=1/m;

% in sample simulations

% calibration of sigma
for i=m:n
	Results = CIRestimation(nonzeros(sb1(i-m+1:i)),Delta,1); 
	r=corr(b1(i-m+1:i),sb1(i-m+1:i),'Type','Spearman');
	k(i)=Results.Params(1); theta(i)=Results.Params(2);  nu(i)=Results.Params(3);
	rho(i)=r;
end

s=rng; % random state 
% fix the random vectors eps0,eps2, for the sake of reproducibility
   eps0=mean(randn(1,N)); eps2=mean(randn(1,N));   
rng(s);

% Eulero-Milnstein scheme

% initialization
sigm_sim(1)=sb1(1);
h=0;
for i=m:n-1
   h=h+1;
   eps1=sqrt(1-rho(i+1)^2)*eps0+rho(i+1)*eps2;

	sigm_sim(i+1)=sb1(i)+k(i+1)*(theta(i+1)-sb1(i))+nu(i+1)*sqrt(sb1(i))*eps2+0.25*(nu(i+1)^2)*(eps2^2-1); 
	sigm_sim(i+1)=abs(sigm_sim(i+1));

end

% objective function    
funS=@(par,data) data(:,1)+par(1)*data(:,1)+sqrt(data(:,2)).*((data(:,1)).^par(2))*eps1+0.5*par(2)*((data(:,1)).^(2*par(2)-1))*(eps1^2-1);
% parameters bounds
lb=[-Inf,1]; ub=[Inf,1];
% set lb=[-Inf,1] for the Heston model
% set lb=[-Inf,0] for the SABR-Heston model
for i=m:n
	xdata1=[0;b1(i-m+1:i-1)]; xdata2=[0;sb1(i-m+1:i-1)];
	xdata=[xdata1,xdata2];
	par=lsqcurvefit(funS,[0,1],xdata,b1(i-m+1:i),lb,ub);  
	% set par0=[0,1] for Heston model
	% set par0=[0,0.5] for SABR-Heston model
	zeta(i)=par(1); gamma(i)=par(2);
end

% Eulero-Milnstein scheme

% initialization
S_simj(1)=b1(1);   sigm_simj(1)=sb1(1); S_sim(1)=b1(1);   sigm_sim(1)=sb1(1);
h=0;
for i=m:n-1
   h=h+1;
	eps1=sqrt(1-rho(h)^2)*eps0+rho(h)*eps2;
	S_sim(i+1)=b1(i)+zeta(i+1)*b1(i)+sqrt(sb1(i))*((b1(i))^gamma(i+1))*eps1+0.5*gamma(i+1)*((b1(i))^(2*gamma(i+1)-1))*(eps1^2-1);
end

%% out of sample simulations
% initialization

F_simj(1)=b1(1);   sigm_simj(1)=sb1(1); F_sim(1)=b1(1);   
h=0;
for i=m:n-1
	h=h+1;
	F_simj(1)=b1(i); signm_simj(1)=sb1(i);
	for j=1:m  % Delta=1/m
		sigm_simj(j+1)=abs(sigm_simj(j)+k(i)*(theta(i)-sigm_simj(j))+nu(i)*sqrt(sigm_simj(j))*eps2+0.25*(nu(i)^2)*(eps2^2-1));
		F_simj(j+1)=F_simj(j)+zeta(i)*F_simj(j)*Delta+sqrt(sb1(i)*Delta)*((F_simj(j))^gamma(i))*eps1+0.5*gamma(i)*((F_simj(j))^(2*gamma(i)-1))*((sqrt(Delta)*eps1)^2-Delta);
	end
	F_sim(i+1)=F_simj(end);  
end


% statistics 
R(1)=NRMSE(b1(m+1:n),S_sim(m+1:n)); 
R(2)=NRMSE(b1(m+1:n),F_sim(m+1:n)); 
R
Ma(1)=MAPE(b1(m+1:n),S_sim(m+1:n)); 
Ma(2)=MAPE(b1(m+1:n),F_sim(m+1:n)); 
Ma
ma(1)=mae(b1(m+1:n),S_sim(m+1:n)');
ma(2)=mae(b1(m+1:n),F_sim(m+1:n)'); 
ma/(max(b1)-min(b1))
N(1)=NMSE(b1(m+1:n),S_sim(m+1:n)); 
N(2)=NMSE(b1(m+1:n),F_sim(m+1:n)); 
N

% plot

subplot(2,1,1)
plot(b1,'LineWidth',1.5);
hold on
plot([m+1:n],S_sim(m+1:n),'LineWidth',1.5);
legend('Brent ret.','SABR-Heston Sim.');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);

subplot(2,1,2)
plot(b1,'LineWidth',1.5);
hold on
plot([m+1:n],F_sim(m+1:n),'LineWidth',1.5);
legend('Brent ret.','SABR-Heston for.');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);
