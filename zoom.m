

function alpha = zoom(f,gradf,gradphi_0,c1,c2,x,phi_0,phi_prev,pk,alpha_lo,alpha_hi)

% Zoom function - each iteration generates alpha between alpha_lo and
% alpha_hi, and then replaces one of these endpoints by alpha, in such a
% way that Wolfe condition is satisfied in the interval!

while true
    
    % Choose an alpha in (alpha_lo, alpha_hi). There are many ways to
    % select the alpha (a) interpolation or (b) just select a midpoint or (c)
    
    alpha  = (alpha_lo+alpha_hi)/2;
    xi     = x + alpha*pk;
    phi_xi = f(xi);
    gradphi_xi = gradf(xi)'*pk;

    
    if phi_xi > phi_0 + c1*alpha*gradphi_0 || phi_xi >= phi_prev
        alpha_hi = alpha;
    else
        gradphi_xi = gradf(xi)'*pk;
        % Check strong Wolfe condition
        if abs(gradphi_xi) <= -c2*gradphi_0
            return;
        end
        
        if gradphi_xi*(alpha_hi-alpha_lo) >= 0
            alpha_hi = alpha_lo;
            alpha_lo = alpha;
            phi_prev = phi_xi;
        end
    end
    
    if alpha < 1e-5
        return
    end
        
        
    end
    
end
