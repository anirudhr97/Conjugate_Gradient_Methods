function [matrix] = generatePDMatrix(n, max_eig)
	%%%%%%%%%%%%%%
	% Function to generate a nxn positive definite matrix.
	% Argument: 
	% 	n			Size of the square matrix
	% Return:
	% 	matrix		nxn positive definite matrix
	%%%%%%%%%%%%%%
	
	%Generating a random nxn matrix
	temp = rand(n, n);
	%Making a nxn symmetric matrix from random nxn matrix
	temp = 0.5*(temp + temp');
	%Eigen decomposition of nxn symmetric matrix
	[V, D] = eig(temp);
	%Choosing positive eigen values to construct our positive definite matrix
	eig_vals = randi(max_eig, n, 1);
	%Constructing the positive definite matrix using eig_vals and V corresponding to matrix temp
	matrix = V*diag(eig_vals)*V';
end
 
