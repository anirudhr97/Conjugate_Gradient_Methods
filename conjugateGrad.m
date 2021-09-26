function [x_hist, gf_hist] = conjugateGrad(A, b, x0, max_iter, tolerance)
	%%%%%%%%%%%%%%
	% Function to run the Conjugate Gradient Method(CGM).
	% Argument: 
	% 	A			Symmetric PD matrix describing the linear system
	% 	b			Target of the linear system
	% 	x0			Initial value of 'x' to start the iterative algorithm
	% 	max_iter	Maximum number of iterations allowed
	% 	tolerance	Tolerance for norm of gradient of loss function at convergence
	% Return:
	% 	x_hist		History of values taken by x through the iterations
	%   gf_hist     History of norm of gradient of loss function through the iterations
	%%%%%%%%%%%%%%

	%The gradient of loss function
	gradf = @(x) A*x - b;

	%x = x + alpha p, where alpha = r'*r/(p'*A*p)
	x = x0;
	r = gradf(x);
	p = -r;
	%Storing x and norm of gradient in the history
	x_hist = [x];
	gf_hist = [norm(r)];
	k = 0;
	%Displaying values before any iteration
	disp(['CG: k=0, gf=' num2str(norm(r))]);

	while norm(r) > tolerance && k < max_iter
		%Step size calculated according to CGM		
		step_size = r'*r/(p'*A*p);
		%Update of x
		x_new = x + step_size * p;

		%x is updated. Updating parameters for next iteration
		r_new = r + step_size * A * p;
		beta_ = r_new'*r_new/(r'*r);
		p = -r_new + beta_ *p;
		r = r_new;
		x = x_new;

		%Storing the norm of grad and updated x
		gf_hist = [gf_hist, norm(r_new)];
		x_hist = [x_hist, x_new];
		%Updating variables for next iteration
		k = k+1;
		disp(['CG: k=' num2str(k) ', gf=' num2str(norm(r_new))]);
	end
end
