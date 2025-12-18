function [f, x_grid, y_grid] = ZZfunc_calcSpikeDensity_v1(x, y, grid_size, x_lim, y_lim)
%KDE2D_SIMPLE  Simple 2D kernel density estimate on a regular grid.
%
%   [f, x_grid, y_grid] = kde2d_simple(x, y, grid_size, x_lim, y_lim)
%
%   Inputs:
%     x, y      : coordinate vectors
%     grid_size : step size for grid
%     x_lim     : [xmin xmax]
%     y_lim     : [ymin ymax]
%
%   Outputs:
%     f         : estimated density on the grid (pdf values)
%     x_grid    : x coordinates of grid
%     y_grid    : y coordinates of grid

    % Make grid
    x_grid = x_lim(1):grid_size:x_lim(2);
    y_grid = y_lim(1):grid_size:y_lim(2);
    [Xq, Yq] = meshgrid(x_grid, y_grid);

    % KDE using MATLAB's ksdensity
    [f, ~] = ksdensity([x(:), y(:)], [Xq(:), Yq(:)]);

    % Reshape back into grid
    f = reshape(f, size(Xq));

end


