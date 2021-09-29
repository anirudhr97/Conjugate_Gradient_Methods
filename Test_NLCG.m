clear all
close all
clc


% Ackley Function, global minimum at (0,0)

syms x
a = 20; b = 0.2; c = 2*pi;
f     = @(x) -a*exp(-b*sqrt(0.5*(x(1)^2 + x(2)^2))) -exp(0.5*(cos(c*x(1))+cos(c*x(2)))) + a + exp(1);
gradf = @(x)[(c*exp(cos(c*x(1))/2 + cos(c*x(2))/2)*sin(c*x(1)))/2 + (a*b*x(1)*exp(-b*(x(1)^2/2 + x(2)^2/2)^(1/2)))/(2*(x(1)^2/2 + x(2)^2/2)^(1/2));
(c*exp(cos(c*x(1))/2 + cos(c*x(2))/2)*sin(c*x(2)))/2 + (a*b*x(2)*exp(-b*(x(1)^2/2 + x(2)^2/2)^(1/2)))/(2*(x(1)^2/2 + x(2)^2/2)^(1/2))];

x0 = [1;1];
c1 = 0.1;
c2 = 0.5;
rho = 2;

f1 = @(x) 4*x(1)^2+x(2)^2;
gradf1 =@(x) [8*x(1);2*x(2)];
pk = -gradf1(x0);
phi_xk = f1(x0);
gradphi_prev = gradf1(x0)'*pk;
[x_k,val,f_hist,x_hist] = NL_CG1(f1,gradf1,x0,pk,c1,c2,rho,'fletcher');

syms x1 y1
% f1     = @(x1,y1) -a*exp(-b*sqrt(0.5*(x1^2 + y1^2))) -exp(0.5*(cos(c*x1)+cos(c*y1))) + a + exp(1);
fsurf(f1,[-1,1,-1,1])
shading interp
hold on
for k = 1 : 10
   scatter3(x_hist(1,k),x_hist(2,k),f_hist(k),'o','fill'); hold on
end