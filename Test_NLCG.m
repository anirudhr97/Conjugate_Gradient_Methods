clear all
close all
clc


% Ackley Function, global minimum at (0,0)

syms x
a = 20; b = 0.2; c = 2*pi;
f     = @(x) -a*exp(-b*sqrt(0.5*(x(1)^2 + x(2)^2))) -exp(0.5*(cos(c*x(1))+cos(c*x(2)))) + a + exp(1);
gradf = @(x)[(c*exp(cos(c*x(1))/2 + cos(c*x(2))/2)*sin(c*x(1)))/2 + (a*b*x(1)*exp(-b*(x(1)^2/2 + x(2)^2/2)^(1/2)))/(2*(x(1)^2/2 + x(2)^2/2)^(1/2));
(c*exp(cos(c*x(1))/2 + cos(c*x(2))/2)*sin(c*x(2)))/2 + (a*b*x(2)*exp(-b*(x(1)^2/2 + x(2)^2/2)^(1/2)))/(2*(x(1)^2/2 + x(2)^2/2)^(1/2))];

x0 = [0.02;0.02];
pk = -gradf(x0)
c1 = 0.001;
c2 = 0.9;
rho = 2;
[x_k,val] = NL_CG(f,gradf,x0,pk,c1,c2,rho,'fletcher');