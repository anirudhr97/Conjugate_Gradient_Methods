

function [x_k,val,f_hist,x_hist,k,alpha_hist] = NL_CG1(f,gradf,x0,pk,c1,c2,rho,method)

% This code performs conjugate gradient on Non - linear functions. 
% The update equation for beta differs by the method of choice 
% (a) Fletcher Reeves (b) Polak Riebere 

% Initialization
tol = 1e-4;
gradf_x0     = gradf(x0);
gradphi_prev = gradf_x0'*pk;
x_prev = x0; phi_xk = f(x0); phi_prev = phi_xk + 10;
x_k = x0; gradf_prev = 1;
gradf_xk = gradf_x0; k = 0;
x_hist = x_k; f_hist = f(x_k); alpha_hist = [];  

while norm(gradf_xk) > tol && norm(phi_xk - phi_prev) > tol 

    alpha_k    = lineSearch(f,gradf,gradphi_prev,c1,c2,x_k,phi_xk,pk,rho);
    x_k        = x_prev + alpha_k*pk;
    gradf_prev = gradf_xk;
    gradf_xk   = gradf(x_k);
    phi_xk     = f(x_k);
    phi_prev   = f(x_prev);
    if strcmp(method,'fletcher')
        beta_k     = (gradf_xk'*gradf_xk)/(gradf_prev'*gradf_prev);
    elseif strcmp(method,'polak')
        beta_k     = (gradf_xk'*(gradf_xk - gradf_prev))/(gradf_prev'*gradf_prev);
        if beta_k < 0
           beta_k = 0;
        end
    elseif strcmp(method,'SD')
        beta_k = 0;
    end
    alpha_hist = [alpha_hist,alpha_k];
    if alpha_k < 1e-5
        beta_k = 0;
    end
    pk         = -gradf_xk + beta_k*pk; 
    k          =  k + 1;
    x_prev     = x_k;
    gradphi_prev = gradf_xk'*pk;
    x_hist = [x_hist,x_k];
    f_hist = [f_hist,f(x_k)];
end

    val = f(x_k);
  
end
