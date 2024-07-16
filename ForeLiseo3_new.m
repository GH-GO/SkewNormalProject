% ForeLiseo3_new

% This procedure simulates the time series Comm_volBrent.mat through the proposed model
% see also ForeLiseo30_new.m

clear all
close all
clc

N=185;
N0=1;

% load time series
load Comm_volBrent
B1=Comm_volBrent(N0:N);
b1=Comm_volBrent(N0:N);
load Comm_fre
b2=Comm_fre(N0:N);
load Comm_press
b3=Comm_press(N0:N);
load Comm_volgold
b4=Comm_volgold(N0:N);


n=length(b1);
m=12+N0-1;  


for i=m:n-1
    i
	X(:,1)=b2(i-m+1:i); X(:,2)=b3(i-m+1:i); X(:,3)=b4(i-m+1:i); v=b1(i-m+1:i);
	% Calibration
	m0=mean(X(:,1))/i; sigm0=std(X(:,1))/sqrt(i);
	for j=1:2
		M(j)=mean(X(:,j))/i; sigm(j)=std(X(:,j))/sqrt(i);
	end
	[mu0,sigma0,beta0,rho,mu,sigma,beta]=Calibration_SkewGBM(X,1,m,i);
	[Alphaa]=cal_Liseo3(X,v);
	Alpha(i,:)=abs(Alphaa);
	alpha=Alphaa;
	par1=[alpha(1),alpha(2),alpha(3),m0,sigm0,M(1),sigm(1),M(2),sigm(2),rho(1),rho(2)];
	delta=beta0/sqrt(1+beta0^2);
	par2=[alpha(1),alpha(2),alpha(3),mu0,sigma0,mu(1),sigma(1),mu(2),sigma(2),rho(1),rho(2),delta];

	% in sample simulations
	[s1,s2]=in_sample_for3(par1,par2,1);
	S1o(i+1)=s1; S2o(i+1)=s2;
	% out of sample simulations
	[f1,f2]=out_sample_for3(par1,par2,1);
	F1o(i+1)=f1; F2o(i+1)=f2;

	S1(i+1)=B1(i+1)+S1o(i+1);
	S2(i+1)=B1(i+1)+S2o(i+1);
	F1(i+1)=B1(i)+F1o(i+1);
	F2(i+1)=B1(i)+F2o(i+1);

end

% Standard errors
for i=m:n-1
	C=Fish_matrix([b2,b3,b4],b1,Alpha(i,:));
	SE(:,i)=sqrt(diag(inv(abs(C))));
end

S1(1:m)=NaN; S2(1:m)=NaN; 
F1(1:m)=NaN; F2(1:m)=NaN; 

% statistics 
R(1)=NRMSE(b1(m+1:n),S1(m+1:n)); R(2)=NRMSE(b1(m+1:n),S2(m+1:n)); R(3)=NRMSE(b1(m+1:n),F1(m+1:n)); R(4)=NRMSE(b1(m+1:n),F2(m+1:n));
R
Ma(1)=MAPE(b1(m+1:n),S1(m+1:n)); Ma(2)=MAPE(b1(m+1:n),S2(m+1:n)); Ma(3)=MAPE(b1(m+1:n),F1(m+1:n)); Ma(4)=MAPE(b1(m+1:n),F2(m+1:n));
Ma
ma(1)=mae(b1(m+1:n),S1(m+1:n)'); ma(2)=mae(b1(m+1:n),S2(m+1:n)'); ma(3)=mae(b1(m+1:n),F1(m+1:n)'); ma(4)=mae(b1(m+1:n),F2(m+1:n)');
ma/(max(b1)-min(b1))
N(1)=NMSE(b1(m+1:n),S1(m+1:n)); N(2)=NMSE(b1(m+1:n),S2(m+1:n)); N(3)=NMSE(b1(m+1:n),F1(m+1:n)); N(4)=NMSE(b1(m+1:n),F2(m+1:n));
N

% plot

subplot(2,2,1)
plot(B1,'LineWidth',1.5);
hold on
plot(S1,'LineWidth',1.5);
title('Normal case');
legend('Brent vol.','Simulations');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1 185-m]);
set(gca,'FontSize',12);

subplot(2,2,2)
plot(B1,'LineWidth',1.5);
hold on
plot(S2,'LineWidth',1.5);
title('Skew-normal case');
legend('Brent vol.','Simulations');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1 185-m]);
set(gca,'FontSize',12);

subplot(2,2,3)
plot(B1,'LineWidth',1.5);
hold on
plot(F1,'LineWidth',1.5);
title('Normal case');
legend('Brent vol.','Forecasts');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1 185-m]);
set(gca,'FontSize',12);

subplot(2,2,4)
plot(B1,'LineWidth',1.5);
hold on
plot(F2,'LineWidth',1.5);
title('Skew-normal case');
legend('Brent vol.','Forecasts');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1 185-m]);
set(gca,'FontSize',12);

figure
% SE plot
E=real(SE);

subplot(2,2,1)
plot(E(1,m+1:184),'LineWidth',1.5);
ylabel({'$SE(\alpha_1)$'},'FontSize',16,'Interpreter','latex');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1 185-m]);
set(gca,'FontSize',12);

subplot(2,2,2)
plot(E(2,m+1:184),'LineWidth',1.5);
ylabel({'$SE(\alpha_2)$'},'FontSize',16,'Interpreter','latex');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1 185-m]);
set(gca,'FontSize',12);

subplot(2,2,3)
plot(E(3,m+1:184),'LineWidth',1.5);
ylabel({'$SE(\alpha_3)$'},'FontSize',16,'Interpreter','latex');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1 185-m]);
set(gca,'FontSize',12);

subplot(2,2,4)
plot(E(4,m+1:184),'LineWidth',1.5);
ylabel({'$SE(\alpha_4)$'},'FontSize',16,'Interpreter','latex');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1 185-m]);
set(gca,'FontSize',12);