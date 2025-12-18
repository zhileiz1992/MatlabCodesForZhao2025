function interp_d = ZZfunc_interp_window(d, n_interp)
% interp_window  Linearly interpolates a 2D array along the first dimension.
%
%   interp_d = interp_window(d, n_interp)
%
%   INPUT:
%       d         - n × m array, where n is time or sequence length, m is features
%       n_interp  - number of interpolation points (e.g., 100)
%
%   OUTPUT:
%       interp_d  - n_interp × m array, linearly interpolated over [0, 1]

    % Validate input
    if isempty(d) || size(d, 1) < 2
        error('Input array must have at least two rows for interpolation.');
    end

    % Original and target time axes (normalized to [0, 1])
    n_orig = size(d, 1);
    orig_t = linspace(0, 1, n_orig);
    target_t = linspace(0, 1, n_interp);

    % Preallocate output
    interp_d = zeros(n_interp, size(d, 2));

    % Interpolate each column (feature) separately
    for col = 1:size(d, 2)
        interp_d(:, col) = interp1(orig_t, d(:, col), target_t, 'linear');
    end
end
