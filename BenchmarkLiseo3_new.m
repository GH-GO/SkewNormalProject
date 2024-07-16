% BenchmarkLiseo3_new

% This procedure performs the NRM benchmark for the time series Comm_volBrent (and its regressors) 

close all;
clear all;
clc


% load time series

load Comm_volBrent
b1=Comm_volBrent;
load Comm_fre
b2=Comm_fre;
load Comm_press
b3=Comm_press;
load Comm_volgold
b4=Comm_volgold;
n=length(b1);
m=12;  

% in sample simulations
V(:,1)=b1; V(:,2)=b2;  V(:,3)=b3+b4; V(:,4)=[0;V(1:end-1,1)];
modelfun1 = @(b,v) exp(-b(1).*v(:,2)).*(b(2).*v(:,3)+b(3))+v(:,4);
b0=[0,0,0];
mdl1 = fitnlm(V(m:n,:),b1(m:n),modelfun1,b0) 
Ssim=mdl1.Fitted;

% out of sample simulations
Xfor(1:m)=NaN;
for i=m:n-1
    i
	V(i-m+2:i,1)=b1(i-m+2:i); V(i-m+2:i,2)=b2(i-m+2:i); V(i-m+2:i,3)=b3(i-m+2:i)+b4(i-m+2:i); V(i-m+2:i,4)=V(i-m+1:i-1,1);
	mdl1 = fitnlm(V(i-m+1:i,:),b1(i-m+1:i),modelfun1,b0);
	xpred = predict(mdl1,V(i-m+1:i,:));
	Xfor(i+1)=xpred(end);
end

% statistics
R(1)=NRMSE(b1(m:n),Ssim); R(2)=NRMSE(b1(m+1:n),Xfor(m+1:n));
R
Ma(1)=MAPE(b1(m:n),Ssim); Ma(2)=MAPE(b1(m+1:n),Xfor(m+1:n));
Ma
ma(1)=mae(b1(m:n),Ssim); ma(2)=mae(b1(m+1:n),Xfor(m+1:n));
ma/(max(b1)-min(b1))
N(1)=NMSE(b1(m:n),Ssim); N(2)=NMSE(b1(m+1:n),Xfor(m+1:n));
N

% plot

subplot(2,1,1)
plot(b1(m:n),'LineWidth',1.5);
hold on
plot(Ssim,'LineWidth',1.5);
legend('Brent ret.','NRM Sim.');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);

subplot(2,1,2)
plot(b1,'LineWidth',1.5);
hold on
plot(Xfor,'LineWidth',1.5);
legend('Brent ret.','NRM for.');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);

