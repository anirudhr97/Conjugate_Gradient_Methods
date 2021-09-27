[33mcommit 356dbdae2a946c4c3618edcc9438ef6bfffcdab2[m
Author: Anirudh R <r.anirudhsrinivas@gmail.com>
Date:   Mon Sep 27 00:34:11 2021 +0530

    Updated conjugateGrad function and wrote main to test

[1mdiff --git a/conjugateGrad.m b/conjugateGrad.m[m
[1mindex 00771fc..1e15d2a 100644[m
[1m--- a/conjugateGrad.m[m
[1m+++ b/conjugateGrad.m[m
[36m@@ -1,7 +1,49 @@[m
[31m-function [outputArg1,outputArg2] = conjugateGrad(inputArg1,inputArg2)[m
[31m-%UNTITLED2 Summary of this function goes here[m
[31m-%   Detailed explanation goes here[m
[31m-outputArg1 = inputArg1;[m
[31m-outputArg2 = inputArg2;[m
[31m-end[m
[32m+[m[32mfunction [x_hist, gf_hist] = conjugateGrad(A, b, x0, max_iter, tolerance)[m
[32m+[m	[32m%%%%%%%%%%%%%%[m
[32m+[m	[32m% Function to run the Conjugate Gradient Method(CGM).[m
[32m+[m	[32m% Argument:[m[41m [m
[32m+[m	[32m% 	A			Symmetric PD matrix describing the linear system[m
[32m+[m	[32m% 	b			Target of the linear system[m
[32m+[m	[32m% 	x0			Initial value of 'x' to start the iterative algorithm[m
[32m+[m	[32m% 	max_iter	Maximum number of iterations allowed[m
[32m+[m	[32m% 	tolerance	Tolerance for norm of gradient of loss function at convergence[m
[32m+[m	[32m% Return:[m
[32m+[m	[32m% 	x_hist		History of values taken by x through the iterations[m
[32m+[m	[32m%   gf_hist     History of norm of gradient of loss function through the iterations[m
[32m+[m	[32m%%%%%%%%%%%%%%[m
[32m+[m
[32m+[m	[32m%The gradient of loss function[m
[32m+[m	[32mgradf = @(x) A*x - b;[m
[32m+[m
[32m+[m	[32m%x = x + alpha p, where alpha = r'*r/(p'*A*p)[m
[32m+[m	[32mx = x0;[m
[32m+[m	[32mr = gradf(x);[m
[32m+[m	[32mp = -r;[m
[32m+[m	[32m%Storing x and norm of gradient in the history[m
[32m+[m	[32mx_hist = [x];[m
[32m+[m	[32mgf_hist = [norm(r)];[m
[32m+[m	[32mk = 0;[m
[32m+[m	[32m%Displaying values before any iteration[m
[32m+[m	[32mdisp(['CG: k=0, gf=' num2str(norm(r))]);[m
 [m
[32m+[m	[32mwhile norm(r) > tolerance && k < max_iter[m
[32m+[m		[32m%Step size calculated according to CGM[m[41m		[m
[32m+[m		[32mstep_size = r'*r/(p'*A*p);[m
[32m+[m		[32m%Update of x[m
[32m+[m		[32mx_new = x + step_size * p;[m
[32m+[m
[32m+[m		[32m%x is updated. Updating parameters for next iteration[m
[32m+[m		[32mr_new = r + step_size * A * p;[m
[32m+[m		[32mbeta_ = r_new'*r_new/(r'*r);[m
[32m+[m		[32mp = -r_new + beta_ *p;[m
[32m+[m		[32mr = r_new;[m
[32m+[m		[32mx = x_new;[m
[32m+[m
[32m+[m		[32m%Storing the norm of grad and updated x[m
[32m+[m		[32mgf_hist = [gf_hist, norm(r_new)];[m
[32m+[m		[32mx_hist = [x_hist, x_new];[m
[32m+[m		[32m%Updating variables for next iteration[m
[32m+[m		[32mk = k+1;[m
[32m+[m		[32mdisp(['CG: k=' num2str(k) ', gf=' num2str(norm(r_new))]);[m
[32m+[m	[32mend[m
[32m+[m[32mend[m
[1mdiff --git a/main.m b/main.m[m
[1mindex e69de29..f9f884c 100644[m
[1m--- a/main.m[m
[1m+++ b/main.m[m
[36m@@ -0,0 +1,61 @@[m
[32m+[m[32m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[m
[32m+[m[32m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[m
[32m+[m[32m% EE5121 - Convex Optimization[m
[32m+[m[32m% Conjugate Gradient Method[m
[32m+[m[32m%[m
[32m+[m[32m% Authors: Anirudh R, Sreekar Sai R, Chandan Bhat[m
[32m+[m[32m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[m
[32m+[m[32m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[m
[32m+[m
[32m+[m[32m%Generating a symmetric positive definite matrix[m
[32m+[m[32mmat = generatePDMatrix(5);[m
[32m+[m[32m%Extracting eigen values of B[m
[32m+[m[32meig_values = eig(mat);[m
[32m+[m[32m%Checking if eigen values are all >0[m
[32m+[m[32mif eig_values > 0[m
[32m+[m[32m    disp(['Matrix generated is positive definite.']);[m
[32m+[m[32mend[m
[32m+[m
[32m+[m[32m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[m
[32m+[m
[32m+[m[32m%'size' determines the size of A and b[m
[32m+[m[32mN = 50;[m
[32m+[m[32m%Amplitude of noise added to original b[m
[32m+[m[32mnoise = 0.0;[m
[32m+[m[32m%Generating matrix A and column vector b for calculations[m
[32m+[m[32mA = generatePDMatrix(N);[m
[32m+[m[32mb_orig = rand(N,1);[m
[32m+[m[32mb = b_orig + noise * rand(N,1);[m
[32m+[m
[32m+[m[32m%Optimal solution to the linear equation Ax = b[m
[32m+[m[32mx_opt = linsolve(A, b);[m
[32m+[m
[32m+[m[32m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[m
[32m+[m
[32m+[m[32m%Parameters for running the iterative method[m
[32m+[m[32m%Random initial point[m
[32m+[m[32mx0 = rand(N, 1);[m
[32m+[m[32m%Maximum number of iterations[m
[32m+[m[32mmax_iter = N;[m
[32m+[m[32m%Tolerance for norm of gradient of loss at convergence[m
[32m+[m[32mtolerance = 1e-6;[m
[32m+[m
[32m+[m[32m%Calling the conjugateGrad function[m
[32m+[m[32m[x_hist, gf_hist] = conjugateGrad(A, b, x0, max_iter, tolerance);[m
[32m+[m
[32m+[m[32m%Plotting log of norm of gradient of loss function[m
[32m+[m[32mplot(log(gf_hist),'r','LineWidth',1);[m
[32m+[m[32mgrid;[m
[32m+[m[32mtitle('Norm of Gradient of Loss vs Iterations');[m
[32m+[m[32mylabel('$log( \| \nabla \phi (x) \| )$','interpreter','latex');[m
[32m+[m[32mxlabel('Iterations');[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m
