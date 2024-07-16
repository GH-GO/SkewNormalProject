% egarch_volbrent

% This procedure performs the i-garch model for the time series Comm_volBrent
% see, also igarch.m

clear all
close all
clc

% load time series

load Comm_volBrent
b1=Comm_volBrent;
m=12;
n=length(b1);

% initialization
[Xsim,garchlags,archlags,AR,MA,Fsim]=ArimaGarch(b1);

% in sample simulations
% simulation with normal distribution
[parameters, LL, qSb1, VCVrobust, VCV, scores, diagnostics] = igarch(b1, garchlags,archlags,'NORMAL',2,1);
Sb1=sqrt(qSb1);
% simulations with skew-t distribution
[parameters, LL, qSb12, VCVrobust, VCV, scores, diagnostics] = igarch(b1, garchlags,archlags,'SKEWT',2,1);
Sb12=sqrt(qSb12);

% statistics
RS(1)=NRMSE(b1,Sb1);      RS(2)=NRMSE(b1,Sb12);
MS(1)=MAPE(b1,Sb1);   MS(2)=MAPE(b1,Sb12);
maS(1)=mae(b1,Sb1);      maS(2)=mae(b1,Sb12);
maS/(max(b1)-min(b1))
NS(1)=NMSE(b1,Sb1);      NS(2)=NMSE(b1,Sb12);

% out of sample simulations
% simulation with normal distribution
[parameters, LL, qFb1, VCVrobust, VCV, scores, diagnostics] = igarch(Fsim(m+1:end)', garchlags,archlags,'NORMAL',2,1);
Fb1=sqrt(qFb1);
% simulation with skew-t distribution
[parameters, LL, qFb12, VCVrobust, VCV, scores, diagnostics] = igarch(Fsim(m+1:end)', garchlags,archlags,'SKEWT',2,1);
Fb12=sqrt(qFb12);


% statistics
RO(1)=NRMSE(b1(m+1:end),Fb1);       RO(2)=NRMSE(b1(m+1:end),Fb12);
MO(1)=MAPE(b1(m+1:end),Fb1);        MO(2)=MAPE(b1(m+1:end),Fb12);
maO(1)=mae(b1(m+1:end),Fb1);       maO(2)=mae(b1(m+1:end),Fb12);
maO/(max(b1)-min(b1))
NO(1)=NMSE(b1(m+1:end),Fb1);           NO(2)=NMSE(b1(m+1:end),Fb12);

% plot
subplot(2,1,1)
plot(b1,'LineWidth',1.5);
hold on
plot(Sb12,'LineWidth',1.5);
legend('Brent ret.','IGARCH (Skew-t) Sim.');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);

subplot(2,1,2)
plot(b1,'LineWidth',1.5);
hold on
plot([m+1:n],Fb12,'LineWidth',1.5);
legend('Brent ret.','IGARCH (Skew-t) for.');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);

