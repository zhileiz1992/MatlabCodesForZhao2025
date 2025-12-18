function [mmd] = ZZfunc_mmd_gaussian(X, Y, sigma)
% MMD_GAUSSIAN Computes Maximum Mean Discrepancy between two sets of samples
%   mmd = mmd_gaussian(X, Y, sigma) calculates the MMD between two sets of
%   samples X and Y using a Gaussian (RBF) kernel with bandwidth sigma.
%
% Inputs:
%   X - n x d matrix of n samples from the first distribution in d dimensions
%   Y - m x d matrix of m samples from the second distribution in d dimensions
%   sigma - scalar, bandwidth of the Gaussian kernel
%
% Output:
%   mmd - scalar, the Maximum Mean Discrepancy value

    % Get number of samples
    n = size(X, 1);
    m = size(Y, 1);
    
    % Ensure inputs are valid
    if size(X, 2) ~= size(Y, 2)
        error('X and Y must have the same number of dimensions');
    end
    if sigma <= 0
%         error('Kernel bandwidth sigma must be positive');
      % estimate the sigma using median pairwise distance
      all_data = [X; Y];
      D = pdist2(all_data, all_data);
      sigma = median(D(:));
    end
    
    % Compute pairwise squared Euclidean distances
    % XX term: sum of kernel evaluations between X samples
    XX = 0;
    for i = 1:n
        for j = 1:n
            XX = XX + exp(-sum((X(i,:) - X(j,:)).^2) / (2 * sigma^2));
        end
    end
    XX = XX / (n * n);
    
    % YY term: sum of kernel evaluations between Y samples
    YY = 0;
    for i = 1:m
        for j = 1:m
            YY = YY + exp(-sum((Y(i,:) - Y(j,:)).^2) / (2 * sigma^2));
        end
    end
    YY = YY / (m * m);
    
    % XY term: sum of kernel evaluations between X and Y samples
    XY = 0;
    for i = 1:n
        for j = 1:m
            XY = XY + exp(-sum((X(i,:) - Y(j,:)).^2) / (2 * sigma^2));
        end
    end
    XY = XY / (n * m);
    
    % Compute MMD (squared MMD for stability)
    mmd = sqrt(max(XX + YY - 2 * XY, 0));
%     mmd = XX + YY - 2 * XY;
end
