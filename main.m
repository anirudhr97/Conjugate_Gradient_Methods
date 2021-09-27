%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EE5121 - Convex Optimization
% Conjugate Gradient Method
%
% Authors: Anirudh R, Sreekar Sai R, Chandan Bhat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing and closing all figures
clear
close all
%Clear screen
clc
%Seeding the random number generator for reproducibility
rng(5, 'twister');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generating a symmetric positive definite matrix
mat = generatePDMatrix(5, 25);
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
noise = 0;
%Factor multiplied with b_orig
mult = 1;
%Maximum eigen value we are willing to allow
max_eig = N^2;
%Generating matrix A and column vector b for calculations
A = generatePDMatrix(N, max_eig);
b_orig = mult * rand(N,1);
b = b_orig + noise * rand(N,1);

%Optimal solution to the linear equation Ax = b
x_opt = linsolve(A, b_orig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters for running the iterative method
%Random initial point
x0 = rand(N, 1);
%Maximum number of iterations
max_iter = N;
%Tolerance for norm of gradient of loss at convergence
tolerance = 1e-6;

%Calling the conjugateGrad function
[x_hist, gf_hist, time_taken, num_iters] = conjugateGrad(A, b, x0, max_iter, tolerance);

%Plotting log of norm of gradient of loss function
figure;
subplot(2,2,1);
plot(0:num_iters, log(gf_hist),'r','LineWidth',1);
grid;
title('Norm of Gradient of Loss vs Iterations');
ylabel('$log( \| \nabla \phi (x) \| )$','interpreter','latex');
xlabel('Iterations');

%Plotting eigen values of A
subplot(2,2,2);
plot(sort(eig(A),'descend'),'x-','LineWidth',1);
grid;
title('Eigenvalues of A');

%Plotting log of norm of gradient of loss function
subplot(2,2,3);
plot(0:num_iters, log( vecnorm( x_hist- x_opt ) ),'b-x','LineWidth',1); grid;
ylabel('$log(\|x_k-x^*\| )$','interpreter','latex');
xlabel('Iterations');
title('Norm of error vs Iterations');

%Plotting time taken for the CGM
subplot(2,2,4);
plot(0:num_iters, 10^6*time_taken,'b','LineWidth',1);
grid;
title('Time Taken vs Iterations');
ylabel('Time Taken ($\mu s$)','interpreter','latex');
xlabel('Iterations');

%Plotting the values taken by the loss function for each x in x_hist
%Defining the loss function
f = @(x) 0.5*x'*A*x - b'*x;

%Calculating loss value at columns in x_hist
func_vals = zeros(num_iters+1, 1);
for i=1:(num_iters+1)
    func_vals(i, 1) = f(x_hist(:, i));
end

figure;
plot(0:num_iters, func_vals,'r','LineWidth',1);
grid;
title('Loss Function vs Iterations');
ylabel('$\phi (x)$', 'interpreter', 'latex');
xlabel('Iterations');

