function [B, C, u, v] = mmd_matrix_decompose(A)
% MMD_MATRIX_DECOMPOSE Decomposes an MMD matrix into rank-one and residual components
%   [B, C, u, v] = mmd_matrix_decompose(A) decomposes the MMD matrix A into a
%   rank-one component B = u*v' and a residual component C = A - B using MAP
%   estimation with a flat prior on u,v and independent Laplace errors on C.
%
% Inputs:
%   A - T1 x T2 matrix of MMD values
%
% Outputs:
%   B - T1 x T2 rank-one matrix (u*v')
%   C - T1 x T2 residual matrix (A - B)
%   u - T1 x 1 vector
%   v - T2 x 1 vector

    % Get dimensions
    [T1, T2] = size(A);
    
    % Initialize u and v (e.g., using SVD for a good starting point)
    [U, S, V] = svd(A, 'econ');
    u = U(:,1) * sqrt(S(1,1));
    v = V(:,1) * sqrt(S(1,1));
    
    % Objective function: L1-norm of residuals
    objective = @(x) sum(abs(A(:) - kron(x(T1+1:end), x(1:T1))),'all');
    
    % Optimize u and v using fminunc
    options = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton');
    x0 = [u; v]; % Initial guess
    x_opt = fminunc(objective, x0, options);
    
    % Extract u and v
    u = x_opt(1:T1);
    v = x_opt(T1+1:end);
    
    % Compute rank-one and residual components
    B = u * v';
    C = A - B;
end