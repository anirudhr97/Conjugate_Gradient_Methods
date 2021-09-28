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
close all;
%Clear screen
clc;
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

%Optimal solution to the linear equation Ax = b_orig
x_opt = linsolve(A, b_orig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters for running the iterative method
%Random initial point
x0 = rand(N, 1);
%Maximum number of iterations
max_iter = N;
%Tolerance for norm of gradient of loss at convergence
tolerance = 1e-6;

%Calling the conjugateGrad and preconditionedCG functions
[x_hist, gf_hist, time_taken, num_iters] = conjugateGrad(A, b, x0, max_iter, tolerance);
[x_hist_pre, gf_hist_pre, time_taken_pre, num_iters_pre] = preconditionedCG(A, b, x0, max_iter, tolerance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name1 = 'Direct';
name2 = 'Preconditioned';

%Plotting log of norm of gradient of loss function
figure;
subplot(2,3,1);
plot(0:num_iters, log(gf_hist),'r','LineWidth',1,'DisplayName',name1);
hold on;
plot(0:num_iters_pre, log(gf_hist_pre),'b','LineWidth',1,'DisplayName',name2);
grid;
title('Norm of Gradient of Loss vs Iterations', 'FontSize', 10);
ylabel('$log( \| \nabla \phi (x) \| )$','interpreter','latex', 'FontSize', 10);
xlabel('Iterations', 'FontSize', 8);
legend();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting eigen values of A
subplot(2,3,2);
plot(sort(eig(A),'descend'),'Marker', 'x','LineStyle', '-', 'Color', [0.4940 0.1840 0.5560],'LineWidth',1);
grid;
title('Eigenvalues of A', 'FontSize', 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting log of norm of gradient of loss function
subplot(2,3,3);
plot(0:num_iters, log( vecnorm( x_hist- x_opt ) ),'r-x','LineWidth',1,'DisplayName',name1);
hold on;
plot(0:num_iters_pre, log( vecnorm( x_hist_pre- x_opt ) ),'b-x','LineWidth',1,'DisplayName',name2);
grid;
title('Norm of error vs Iterations', 'FontSize', 10);
ylabel('$log(\|x_k-x^*\| )$','interpreter','latex', 'FontSize', 10);
xlabel('Iterations', 'FontSize', 8);
legend();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting time taken
subplot(2,3,4);
plot(0:num_iters, 10^6*time_taken,'r','LineWidth',1,'DisplayName',name1);
hold on;
plot(0:num_iters_pre, 10^6*time_taken_pre,'b','LineWidth',1,'DisplayName',name2);
grid;
title('Time Taken vs Iterations', 'FontSize', 10);
ylabel('Time Taken ($\mu s$)','interpreter','latex', 'FontSize', 10);
xlabel('Iterations', 'FontSize', 8);
legend();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting cumulative time taken vs iterations
%Calculating cumulative times taken
time_taken_cum = zeros(num_iters+1, 1);
time_taken_cum(1) = time_taken(1);
for i=2:num_iters+1
    time_taken_cum(i) = time_taken_cum(i-1) + time_taken(i);
end

time_taken_cum_pre = zeros(num_iters_pre+1, 1);
time_taken_cum_pre(1) = time_taken_pre(1);
for i=2:num_iters_pre+1
    time_taken_cum_pre(i) = time_taken_cum_pre(i-1) + time_taken_pre(i);
end

subplot(2,3,5);
plot(0:num_iters, 10^6*time_taken_cum,'r','LineWidth',1,'DisplayName',name1);
hold on;
plot(0:num_iters_pre, 10^6*time_taken_cum_pre,'b','LineWidth',1,'DisplayName',name2);
grid;
title('Cumulative Time Taken vs Iterations', 'FontSize', 10);
ylabel('Cumulative Time Taken ($\mu s$)','interpreter','latex', 'FontSize', 10);
xlabel('Iterations', 'FontSize', 8);
legend();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting the values taken by the loss function for each x in x_hist and x_hist_pre
%Defining the loss function
f = @(x) 0.5*x'*A*x - b'*x;

%Calculating loss value at columns in x_hist and x_hist_pre
func_vals = zeros(num_iters+1, 1);
for i=1:(num_iters+1)
    func_vals(i, 1) = f(x_hist(:, i));
end

func_vals_pre = zeros(num_iters_pre+1, 1);
for i=1:(num_iters_pre+1)
    func_vals_pre(i, 1) = f(x_hist_pre(:, i));
end

subplot(2,3,6);
plot(0:num_iters, func_vals,'r','LineWidth',1,'DisplayName',name1);
hold on;
plot(0:num_iters_pre, func_vals_pre,'b','LineWidth',1,'DisplayName',name2);
grid;
title('Loss Function vs Iterations', 'FontSize', 10);
ylabel('$\phi (x)$', 'interpreter', 'latex', 'FontSize', 10);
xlabel('Iterations', 'FontSize', 8);
legend();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Saving the plot generated
% print('plot','-dpng');

