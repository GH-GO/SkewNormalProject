% SARIMA_Liseo

% This procedure performs the SARIMA benchmark for the time series Comm_volBrent 
clear all
close all
clc

% load time series
load Comm_volBrent
b1=Comm_volBrent;


L=12;
n=length(b1);

% in sample simulations
[b1sim,garchlagsb1,archlagsb1,ARb1,MAb1]=ArimaGarch(b1);
[Ssim1,SARb1,SMAb1,Fsim1]=Sarima(b1,ARb1,MAb1,L);


R(1,1)=NRMSE(b1(L+1:end),Ssim1(L+1:end));  R(1,2)=NRMSE(b1(L+1:end),Fsim1(L+1:end)); 
R
MA(1,1)=MAPE(b1(L+1:end),Ssim1(L+1:end));  MA(1,2)=MAPE(b1(L+1:end),Fsim1(L+1:end)); 
MA
M(1,1)=mae(b1(L+1:end),Ssim1(L+1:end));  M(1,2)=mae(b1(L+1:end),Fsim1(L+1:end)); 
M/(max(b1)-min(b1))
N(1,1)=NMSE(b1(L+1:end),Ssim1(L+1:end));  N(1,2)=NMSE(b1(L+1:end),Fsim1(L+1:end));
N

% plot
subplot(2,1,1)
plot(b1,'LineWidth',1.5);
hold on
plot(Ssim1,'LineWidth',1.5);
legend('Brent ret.','SARIMA Sim.');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);

subplot(2,1,2)
plot(b1,'LineWidth',1.5);
hold on
plot(Fsim1,'LineWidth',1.5);
legend('Brent ret.','SARIMA for.');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);

