% ForeLiseo30_new

% This procedure simulates the regressor time series Comm_fre.mat, Comm_press.mat, Comm_volgold.mat through the proposed model
% see also ForeLiseo3_new.m

clear all
close all
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
s=rng; % random state 
% fix the random vectors r1,r2, for the sake of reproducibility
r1=mean(randn(1,10^5)); r2=abs(mean(randn(1,10^5)));
rng(s);

n=length(b1);
m=12;  

for i=m:n
    i
	X(:,1)=b2(i-m+1:i); X(:,2)=b3(i-m+1:i); X(:,3)=b4(i-m+1:i); v=b1(i-m+1:i);
	% Calibration
	m0(i)=mean(X(:,1))/i; sigm0(i)=std(X(:,1))/sqrt(i);
	for j=1:2
		M(i,j)=mean(X(:,j+1))/i; sigm(i,j)=std(X(:,j+1))/sqrt(i);
	end
	[smu0,ssigma0,sbeta0,srho,smu,ssigma,sbeta]=Calibration_SkewGBM(X,1,m,i);
	mu0(i)=smu0; sigma0(i)=ssigma0; mu(i,:)=smu; sigma(i,:)=ssigma; delta0(i)=sbeta0/sqrt(1+sbeta0^2);
	for j=1:2
		delta(i,:)=sbeta./sqrt(1+sbeta.^2);
	end
end

for i=m:n-1
	% (S,F) transit
	r2t=delta0(i+1)*r2;
	St1(i+1)=m0(i+1)*i;  % normal
	St2(i+1)=mu0(i+1)*i+sigma0(i+1)*delta0(i+1)*sqrt(2*i/pi);   % skew-normal
	Ft1(i+1)=b2(i)+m0(i);  % normal
	Ft2(i+1)=b2(i)+mu0(i)+sigma0(i)*(2*r2t*(normcdf(r2t/delta0(i))-1)+delta0(i)*sqrt(2/pi)*normpdf(r2t/delta0(i))); % skew-normal

	if delta0(i)==0
		Ft2(i+1)=Ft1(i+1);  
	end

	% (S,F) press
	r2g=delta(i+1,1)*r2;
	Sg1(i+1)=M(i+1,1)*i;  % normal
	Sg2(i+1)=mu(i+1,1)*i+sigma(i+1,1)*delta(i+1,1)*sqrt(2*i/pi);   % skew-normal
	Fg1(i+1)=b3(i)+M(i,1);  % normal
	Fg2(i+1)=b3(i)+mu(i,1)+sigma(i,1)*(2*r2g*(normcdf(r2g/delta(i,1))-1)+delta(i,1)*sqrt(2/pi)*normpdf(r2g/delta(i,1))); % skew-normal

	if delta(i,1)==0
		Fg2(i+1)=Fg1(i+1);  
	end

	% (S,F) vol gold
	r2v=delta(i+1,2)*r2;
	Sv1(i+1)=M(i+1,2)*i;  % normal
	Sv2(i+1)=mu(i+1,2)*i+sigma(i+1,2)*delta(i+1,2)*sqrt(2*i/pi);   % skew-normal
	Fv1(i+1)=b4(i)+M(i,2);  % normal
	Fv2(i+1)=b4(i)+mu(i,2)+sigma(i,2)*(2*r2v*(normcdf(r2v/delta(i,2))-1)+delta(i,2)*sqrt(2/pi)*normpdf(r2v/delta(i,2))); % skew-normal

	if delta(i,2)==0
		Fv2(i+1)=Fv1(i+1);  
	end
end

St1(1:m)=NaN;  St2(1:m)=NaN; Ft1(1:m)=NaN; Ft2(1:m)=NaN;
Sg1(1:m)=NaN;  Sg2(1:m)=NaN; Fg1(1:m)=NaN; Fg2(1:m)=NaN;
Sv1(1:m)=NaN;  Sv2(1:m)=NaN; Fv1(1:m)=NaN; Fv2(1:m)=NaN;

% statistics  transit
Rt(1)=NRMSE(b2(m+1:n),St1(m+1:n)); Rt(2)=NRMSE(b2(m+1:n),St2(m+1:n)); Rt(3)=NRMSE(b2(m+1:n),Ft1(m+1:n)); Rt(4)=NRMSE(b2(m+1:n),Ft2(m+1:n));
Rt
Mat(1)=MAPE(b2(m+1:n)+1,St1(m+1:n)+1); Mat(2)=MAPE(b2(m+1:n)+1,St2(m+1:n)+1); Mat(3)=MAPE(b2(m+1:n)+1,Ft1(m+1:n)+1); Mat(4)=MAPE(b2(m+1:n)+1,Ft2(m+1:n)+1);
Mat


% statistics  press
Rg(1)=NRMSE(b3(m+1:n),Sg1(m+1:n)); Rg(2)=NRMSE(b3(m+1:n),Sg2(m+1:n)); Rg(3)=NRMSE(b3(m+1:n),Fg1(m+1:n)); Rg(4)=NRMSE(b3(m+1:n),Fg2(m+1:n));
Rg
Mag(1)=MAPE(b3(m+1:n),Sg1(m+1:n)); Mag(2)=MAPE(b3(m+1:n),Sg2(m+1:n)); Mag(3)=MAPE(b3(m+1:n),Fg1(m+1:n)); Mag(4)=MAPE(b3(m+1:n),Fg2(m+1:n));
Mag

% statistics  vol gold
Rv(1)=NRMSE(b4(m+1:n),Sv1(m+1:n)); Rv(2)=NRMSE(b4(m+1:n),Sv2(m+1:n)); Rv(3)=NRMSE(b4(m+1:n),Fv1(m+1:n)); Rv(4)=NRMSE(b4(m+1:n),Fv2(m+1:n));
Rv
Mav(1)=MAPE(b4(m+1:n),Sv1(m+1:n)); Mav(2)=MAPE(b4(m+1:n),Sv2(m+1:n)); Mav(3)=MAPE(b4(m+1:n),Fv1(m+1:n)); Mav(4)=MAPE(b4(m+1:n),Fv2(m+1:n));
Mav


% Due to their similarity, we show the differences, i.e. RMSE(S-N)-RMSE(N), MAPE(S-N)-MAPE(N) 

dRt(1)=-Rt(2)+Rt(1); dRt(2)=-Rt(4)+Rt(3);
dRt
dMat(1)=-Mat(2)+Mat(1); dMat(2)=-Mat(4)+Mat(3);
dMat


dRg(1)=-Rg(2)+Rg(1); dRg(2)=-Rg(4)+Rg(3);
dRg
dMag(1)=-Mag(2)+Mag(1); dMag(2)=-Mag(4)+Mag(3);
dMag

dRv(1)=-Rv(2)+Rv(1); dRv(2)=-Rv(4)+Rv(3);
dRv
dMav(1)=-Mav(2)+Mav(1); dMav(2)=-Mav(4)+Mav(3);
dMav

% plot

subplot(3,2,1)
plot(b2,'LineWidth',1.5);
hold on
plot(St1,'LineWidth',1.5);
hold on
plot(Ft1,'--','LineWidth',1.5);
title('Normal case');
legend('Freight','Simulations','Forecasts');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1, length(b2)+1]);
set(gca,'FontSize',12);

subplot(3,2,2)
plot(b2,'LineWidth',1.5);
hold on
plot(St2,'LineWidth',1.5);
hold on
plot(Ft2,'--','LineWidth',1.5);
title('Skew-normal case');
legend('Freight','Simulations','Forecasts');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1, length(b2)+1]);
set(gca,'FontSize',12);

subplot(3,2,3)
plot(b3,'LineWidth',1.5);
hold on
plot(Sg1,'LineWidth',1.5);
hold on
plot(Fg1,'--','LineWidth',1.5);
title('Normal case');
legend('Price pressure','Simulations','Forecasts');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1, length(b2)+1]);
set(gca,'FontSize',12);

subplot(3,2,4)
plot(b3,'LineWidth',1.5);
hold on
plot(Sg2,'LineWidth',1.5);
hold on
plot(Fg2,'--','LineWidth',1.5);
title('Skew-normal case');
legend('Price pressure','Simulations','Forecasts');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1, length(b2)+1]);
set(gca,'FontSize',12);


subplot(3,2,5)
plot(b4(3*m:end),'LineWidth',2);
hold on
plot(Sv1(3*m:end),'LineWidth',1.5);
hold on
plot(Fv1(3*m:end),'--','LineWidth',1.5);
title('Normal case');
legend('Gold vol.','Simulations','Forecasts');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1, length(b2)-3*m+1]);
set(gca,'FontSize',12);



subplot(3,2,6)
plot(b4(3*m:end),'LineWidth',1.5);
hold on
plot(Sv2(3*m:end),'LineWidth',1.5);
hold on
plot(Fv2(3*m:end),'--','LineWidth',1.5);
title('Skew-normal case');
legend('Gold vol.','Simulations','Forecasts');
xlabel({'$t$\,(months)'},'FontSize',16,'Interpreter','latex');
xlim([1, length(b2)-3*m+1]);
set(gca,'FontSize',12);