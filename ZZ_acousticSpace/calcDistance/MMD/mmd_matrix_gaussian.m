function mmd_matrix = mmd_matrix_gaussian(d1, d2, sigma, standardize)
% MMD_MATRIX_GAUSSIAN Computes MMD between all pairs of time points in two 3D arrays
%   mmd_matrix = mmd_matrix_gaussian(d1, d2, sigma, standardize) calculates the MMD
%   between samples at each pair of time points from d1 and d2 using a Gaussian (RBF)
%   kernel, returning a T1 x T2 matrix where mmd_matrix(i,j) is the MMD between
%   time point i of d1 and time point j of d2. Assumes no subsampling is needed.
%   If standardize is true, each feature is standardized across all time points
%   and samples in the combined dataset.
%
% Inputs:
%   d1 - T1 x F x N array, T1 time points, F features, N samples
%   d2 - T2 x F x M array, T2 time points, F features, M samples
%   sigma - scalar or T1 x T2 matrix, Gaussian kernel bandwidth
%   standardize - logical, whether to standardize each feature (default: false)
%
% Output:
%   mmd_matrix - T1 x T2 matrix of MMD values

    % Default value for standardize
    if nargin < 4
        standardize = false;
    end
    
    % Get dimensions
    [T1, F, N] = size(d1);
    [T2, F2, M] = size(d2);
    if F ~= F2
        error('Feature dimensions of d1 and d2 must match');
    end
    
    % Validate sigma
    if isscalar(sigma)
        sigma = sigma * ones(T1, T2); % Broadcast scalar sigma
    elseif ~isequal(size(sigma), [T1, T2])
        error('Sigma must be a scalar or a %d x %d matrix', T1, T2);
    end
    if any(sigma(:) <= 0)
        error('Sigma values must be positive');
    end
    
    % Reshape arrays for easier processing (time points as last dimension)
    d1 = permute(d1, [3, 2, 1]); % N x F x T1
    d2 = permute(d2, [3, 2, 1]); % M x F x T2
    
    % Standardize data if requested (across all time points and samples)
    if standardize
        % Combine all samples and time points
        Z1 = reshape(d1, [N * T1, F]); % (N*T1) x F
        Z2 = reshape(d2, [M * T2, F]); % (M*T2) x F
        Z = [Z1; Z2]; % (N*T1 + M*T2) x F
        mu = mean(Z, 1); % 1 x F, mean per feature
        std_dev = std(Z, 0, 1); % 1 x F, std per feature
        std_dev(std_dev == 0) = 1; % Avoid division by zero
        
        % Standardize d1 and d2
        d1 = (d1 - reshape(mu, [1, F, 1])) ./ reshape(std_dev, [1, F, 1]);
        d2 = (d2 - reshape(mu, [1, F, 1])) ./ reshape(std_dev, [1, F, 1]);
    end
    
    % Initialize MMD matrix
    mmd_matrix = zeros(T1, T2);
    
    % Compute XX terms for d1 (T1 x 1)
    XX = zeros(T1, 1);
    for i = 1:T1
        X = squeeze(d1(:, :, i)); % N x F
        dist_XX = squareform(pdist(X, 'squaredeuclidean')); % N x N
        XX(i) = sum(sum(exp(-dist_XX / (2 * sigma(i, 1)^2)))) / (N * N);
    end
    
    % Compute YY terms for d2 (T2 x 1)
    YY = zeros(T2, 1);
    for j = 1:T2
        Y = squeeze(d2(:, :, j)); % M x F
        dist_YY = squareform(pdist(Y, 'squaredeuclidean')); % M x M
        YY(j) = sum(sum(exp(-dist_YY / (2 * sigma(1, j)^2)))) / (M * M);
    end
    
    % Compute XY terms for all pairs (i,j)
    for i = 1:T1
        X = squeeze(d1(:, :, i)); % N x F
        for j = 1:T2
            Y = squeeze(d2(:, :, j)); % M x F
            dist_XY = pdist2(X, Y, 'squaredeuclidean'); % N x M
            XY = sum(sum(exp(-dist_XY / (2 * sigma(i, j)^2)))) / (N * M);
            mmd_matrix(i, j) = sqrt(max(XX(i) + YY(j) - 2 * XY, 0));
        end
    end
end