%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EE5121 - Convex Optimization
% Conjugate Gradient Method
%
% Authors: Anirudh R, Sreekar Sai R, Chandan Bhat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generating a symmetric positive definite matrix
mat = generatePDMatrix(5);
%Extracting eigen values of B
eig_values = eig(mat);
%Checking if eigen values are all >0
if eig_values > 0
    disp(['Matrix generated is positive definite.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%'size' determines the size of A and b
N = 50;
%Amplitude of noise added to original b
noise = 0.0;
%Generating matrix A and column vector b for calculations
A = generatePDMatrix(N);
b_orig = rand(N,1);
b = b_orig + noise * rand(N,1);

%Optimal solution to the linear equation Ax = b
x_opt = linsolve(A, b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters for running the iterative method
%Random initial point
x0 = rand(N, 1);
%Maximum number of iterations
max_iter = N;
%Tolerance for norm of gradient of loss at convergence
tolerance = 1e-6;

%Calling the conjugateGrad function
[x_hist, gf_hist] = conjugateGrad(A, b, x0, max_iter, tolerance);

%Plotting log of norm of gradient of loss function
plot(log(gf_hist),'r','LineWidth',1);
grid;
title('Norm of Gradient of Loss vs Iterations');
ylabel('$log( \| \nabla \phi (x) \| )$','interpreter','latex');
xlabel('Iterations');










