


function alpha = lineSearch(f,gradf,gradphi_0,c1,c2,x,phi_x,pk,rho)

% This code performs Line search for finding optimal alpha (step length) 
% using strong Wolfe conditions. The algorithm for line search method is 
% used from " Numerical optimization " by Jorge Nocedal Stephen J. Wright
% Chapter 3, Algorithm 3.5.  

% The parameters passed to the search algorithm are 
% f       - objective function
% gradf   - gradient of the function 
% c1      - sufficient decrease constant
% c2      - curvature condition constant
% x       - solution at current iterate
% pk      - search direction  
% phi(alpha) = f(x + alpha p_k) 
% phi_x   - function value at current x 
% gradphi_0 - gradient of function at x


% Parameters used in the algorithm 
% xi        - the updated x (x_prev +  alpha_i*pk) 
% alpha_lo  - alpha whoich supplied to zoom function
% alpha_hi  - higher alpha 
 

% Initialize the alpha by specifying the range it can take. alpha > 0
alpha_0  = 0; alpha_max = 200; alpha_1  = 1; % alpha_1 \in (0,alpha_max)
alpha_i = alpha_1; alpha_prev = alpha_0;
phi_prev  = phi_x; phi_0 = phi_x;
x_prev = x; i = 0;
 
while true && (abs(alpha_prev - alpha_i) > 1e-4)
    
    xi          = x_prev + alpha_i*pk;
    phi_xi      = f(xi);
    gradphi_xi  = gradf(xi)'*pk;
      
    % alpha_1 violates the sufficient condition 
    if  (phi_xi  > phi_0 + c1*alpha_i*gradphi_0) || ((phi_xi >= phi_prev ) && (i > 0))
        alpha = zoom(f,gradf,gradphi_0,c1,c2,x,phi_0,phi_prev,pk,alpha_prev,alpha_i);
        break;
    end
    
    % strong Wolfe check
    if(abs(gradphi_xi) <= -c2*gradphi_0)
        alpha = alpha_i;
        break;
    end
    
    % Check if we moved ahead
    if (gradphi_xi >= 0)
        alpha = zoom(f,gradf,gradphi_0,c1,c2,x,phi_0,phi_prev,pk,alpha_prev,alpha_i);
        break;
    end

    alpha_prev = alpha_i;
    alpha_i    = min(rho*alpha_i,alpha_max);
    phi_prev   = phi_xi;
    i = i + 1;

end

end



