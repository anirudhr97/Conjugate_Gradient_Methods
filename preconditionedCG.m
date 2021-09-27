function [x_hist, gf_hist, time_taken, k] = preconditionedCG(A, b, x0, max_iter, tolerance)
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
	%	time_taken	Time taken for execution of each iteration
	%	k			Number of iterations that actually occurred
	%%%%%%%%%%%%%%

    %     Preconditioning
    C = ichol(sparse(A)); % CHECK THIS
    Minv = inv(C*C');
    
	%Starting measurement of time for initialization tasks
	tic

	%The gradient of loss function
	gradf = @(x) A*x - b;

	%x = x + alpha p, where alpha = r'*r/(p'*A*p)
	x = x0;
	r = gradf(x);
    y = Minv*r;
	p = -y;
	%Storing x and norm of gradient in the history
	x_hist = [x];
	gf_hist = [norm(r)];
	k = 0;

	%Logging the time taken before iterating
	time_taken = toc;

	%Displaying values before the iterations
	disp(['CG: k=0, gf=' num2str(norm(r)) '  , time elapsed: ' num2str(time_taken*10^6) ' micro seconds']);

	while norm(r) > tolerance && k < max_iter
		%Starting measurement of time for this iteration
		tic

		%Step size calculated according to CGM		
		step_size = r'*y/(p'*A*p);
		%Update of x
		x_new = x + step_size * p;
        
		%x is updated. Updating parameters for next iteration
		r_new = r + step_size * A * p;
        y_new = Minv*r_new;

		beta_ = r_new'*y_new/(r'*y);
		p = -y_new + beta_ *p;
		r = r_new;
		x = x_new;
        y = y_new;

		%Storing the norm of grad and updated x
		gf_hist = [gf_hist, norm(r_new)];
		x_hist = [x_hist, x_new];
		%Updating iterate for next iteration
		k = k+1;

		%Ending measurement of time and logging
		time_taken = [time_taken, toc];

		disp(['CG: k=' num2str(k) ', gf=' num2str(norm(r_new)) '  , time elapsed: ' num2str(time_taken(end)*10^6) ' micro seconds']);
	end
	%Displaying relevant details of the CGM execution
	fprintf('\n');
	disp(['Total number of iterations: ' num2str(k) ' , Total time taken: ' num2str(sum(time_taken)*10^6) ' micro seconds']);
end
