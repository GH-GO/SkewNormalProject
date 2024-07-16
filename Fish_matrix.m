function C=Fish_matrix(X,V,a)

% This function calculates the Fisher information matrix 
% INPUT:
% V: time series denoting the dependent variable 
% X: matrix whose columns represent the regeressors of X
% a: model parameters (a1,a2,a3,a4) estimated by cal_Liseo3.m
% OUTPUT: 
% C: Fisher information matrix


n=length(V);

Z1=X(:,1); Z2=X(:,2)+X(:,3); Z3=[0;V(1:end-1,1)];

Start=1;

A11=0;
for i=Start:n
	A11=A11+(Z1(i)*exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)))^2/a(4)+(V(i)-Z3(i)-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)) )*(-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3))*Z1(i)^2 )/a(4);    
end

A12=0;
for i=Start:n
	A12=A12+(Z1(i)*exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)))*(-exp(-a(1)*Z1(i))*Z2(i))/a(4)+(V(i)-Z3(i)-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)) )*(-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3))*Z1(i)*Z2(i) )/a(4);    
end
A21=A12;

A13=0;
for i=Start:n
	A13=A13+(Z1(i)*exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)))*(-exp(-a(1)*Z1(i)))/a(4)+(V(i)-Z3(i)-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)) )*(-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3))*Z1(i)*Z2(i) )/a(4);    
end
A31=A13;

A14=0;
for i=Start:n
	A14=A14-(V(i)-Z3(i)-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)) )*(exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3))*Z1(i) )/a(4)^2;    
end
A41=A14;

A22=0;
for i=Start:n
	A22=A22+(exp(-a(1)*Z1(i))*Z2(i) )^2/a(4);    
end

A23=0;
for i=Start:n
A23=A23+(-exp(-a(1)*Z1(i)) )*(-exp(-a(1)*Z1(i))*Z2(i) )/a(4);    
end
A32=A23;

A24=0;
for i=Start:n
	A24=A24-(-exp(-a(1)*Z1(i))*Z2(i) )*(V(i)-Z3(i)-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)) )/a(4)^2;    
end
A42=A24;

A33=0;
for i=Start:n
A33=A33+(exp(-a(1)*Z1(i)) )^2/a(4);    
end

A34=0;
for i=Start:n
	A34=A34-(-exp(-a(1)*Z1(i)) )*(V(i)-Z3(i)-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)) )/a(4)^2;    
end
A43=A34;

A44=0;
for i=Start:n
	A44=A44-0.5*n/a(4)^2+(V(i)-Z3(i)-exp(-a(1)*Z1(i))*(a(2)*Z2(i)+a(3)) )^2/a(4)^3;    
end

C(1,1)=A11; C(1,2)=A12; C(1,3)=A13; C(1,4)=A14;
C(2,1)=A21; C(2,2)=A22; C(2,3)=A23; C(2,4)=A24;
C(3,1)=A31; C(3,2)=A32; C(3,3)=A33; C(3,4)=A34;
C(4,1)=A41; C(4,2)=A42; C(4,3)=A43; C(4,4)=A44;
