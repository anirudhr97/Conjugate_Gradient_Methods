clear all; close all; clc

method = 'fletcher';
% Ackley Function, global minimum at (0,0)
syms x
a = 20; b = 0.2; c = 2*pi;
f1     = @(x) -a*exp(-b*sqrt(0.5*(x(1)^2 + x(2)^2))) -exp(0.5*(cos(c*x(1))+cos(c*x(2)))) + a + exp(1);
gradf1 = @(x)[(c*exp(cos(c*x(1))/2 + cos(c*x(2))/2)*sin(c*x(1)))/2 + (a*b*x(1)*exp(-b*(x(1)^2/2 + x(2)^2/2)^(1/2)))/(2*(x(1)^2/2 + x(2)^2/2)^(1/2));
    (c*exp(cos(c*x(1))/2 + cos(c*x(2))/2)*sin(c*x(2)))/2 + (a*b*x(2)*exp(-b*(x(1)^2/2 + x(2)^2/2)^(1/2)))/(2*(x(1)^2/2 + x(2)^2/2)^(1/2))];

% Parameters
c1 = 1e-4; c2 = 0.4; rho = 2;

% starting point for the non-linear CG 
% Choose any one and comment the rest 
x0 = [1.2;0.27]; 
%x0 = [-0.55;-0.03]; 
%x0 = [0.57;0.5]; 
%x0 = [0.57;0.58];
%x0 = [1;0.2];
x0 = [1.2;0.24];
%x0 = [0.03; 0.31]
% Test for a simple convex function, uncomment following lines
% f1 = @(x) 4*x(1)^2+x(2)^2;
% gradf1 =@(x) [8*x(1);2*x(2)];

pk = -gradf1(x0); % inital search direction
phi_xk = f1(x0); gradphi_prev = gradf1(x0)'*pk;
[~,~,f_hist,x_hist,k,alpha_hist] = NL_CG1(f1,gradf1,x0,pk,c1,c2,rho,method);


% Plotting
plotting = true;

if plotting
    syms x1 y1
    f2  = @(x1,y1) -a*exp(-b.*sqrt(0.5*(x1.^2 + y1.^2))) -exp(0.5.*(cos(c.*x1)+cos(c.*y1))) + a + exp(1);
    sampling = 0.01;
    xc  = -1.5 : sampling : 1.5;  
    [Xc,Yc] = meshgrid(xc);
    f_out = f2(Xc,Yc);
    figure
    surf(Xc,Yc,f_out) 
    shading interp
    hold on
    plot3(x_hist(1,:),x_hist(2,:),f_hist(:),'*-r','MarkerSize',10,'LineWidth',4); hold on   

    figure
    [~, contourObj] = contourf(f_out);
    % This is the secret that 'keeps' the transparency.
    eventFcn = @(srcObj, e) updateTransparency(srcObj);
    addlistener(contourObj, 'MarkedClean', eventFcn);
    hold on
    plot(x_hist(1,:)*(1/sampling)+floor(length(xc)/2),x_hist(2,:)*(1/sampling)+floor(length(xc)/2),'*-','MarkerSize',5,'LineWidth',2); hold on   
    
end

display = true;
if display 

fprintf('------------------------------------------------------------------\n')

fprintf('Initial Point : [');
fprintf('%g, ', x0(1:end-1));
fprintf('%g]\n', x0(end));

fprintf('Minimum Point : [');
fprintf('%g, ', x_hist(1:end-1,end));
fprintf('%g]\n', x_hist(end,end));

fprintf('Minimum value : %g \n',f_hist(end));
fprintf('Method : %s \n',method);

fprintf('Steps taken %d \n',k);

fprintf('------------------------------------------------------------------\n')

end

comparison = true;
if comparison

method = 'SD';    
[~,~,f_hist_sd,x_hist_sd,k_sd,alpha_hist_sd] = NL_CG1(f1,gradf1,x0,pk,c1,c2,rho,method);
method = 'fletcher';    
[~,~,f_hist_fr,x_hist_fr,k_fr,alpha_hist_fr] = NL_CG1(f1,gradf1,x0,pk,c1,c2,rho,method);
method = 'polak';    
[~,~,f_hist_pr,x_hist_pr,k_pr,alpha_hist_pr] = NL_CG1(f1,gradf1,x0,pk,c1,c2,rho,method);

fprintf('------------------------------------------------------------------\n')

fprintf('Minimum value : %g \n',f_hist_sd(end));
fprintf('Method : %s \n','SD');
fprintf('Minimum value : %g \n',f_hist_fr(end));
fprintf('Method : %s \n','fletcher');
fprintf('Minimum value : %g \n',f_hist_pr(end));
fprintf('Method : %s \n','polak');



fprintf('Steps taken SD %d \n',k_sd);
fprintf('Steps taken FR %d \n',k_fr);
fprintf('Steps taken PR %d \n',k_pr);

fprintf('------------------------------------------------------------------\n')

end
    
plotting_cmp = true;    
if plotting_cmp
   figure
   plot(f_hist_sd); hold on
   plot(f_hist_fr); hold on
   plot(f_hist_pr); hold on
end



