function [dens_s, x_edges, y_edges, x_centers, y_centers] = ZZfunc_calcSpikeDensity_count_V1(x, y, num_bin, x_lim, y_lim, sigma_bins)
%POINT_DENSITY_SMOOTH_BINS  Count-based 2D density with fixed number of bins,
%   then Gaussian-smooth.
%
%   [dens_s, x_edges, y_edges, x_centers, y_centers] =
%       point_density_smooth_bins(x, y, num_bin, x_lim, y_lim, sigma_bins)
%
%   Inputs:
%     x, y        : point coordinates (vectors, same length)
%     num_bin     : scalar or [Nx Ny], number of bins along x and y
%     x_lim       : [xmin xmax]
%     y_lim       : [ymin ymax]
%     sigma_bins  : Gaussian sigma in *bins* (default = 1)
%
%   Outputs:
%     dens_s      : smoothed counts per bin (Ny-by-Nx)
%     x_edges     : bin edges along x
%     y_edges     : bin edges along y
%     x_centers   : bin centers along x
%     y_centers   : bin centers along y

    if nargin < 6 || isempty(sigma_bins), sigma_bins = 1; end

    % --- number of bins along x and y ---
    if isscalar(num_bin)
        Nx = num_bin;
        Ny = num_bin;
    else
        Nx = num_bin(1);
        Ny = num_bin(2);
    end

    % --- bin edges ---
    x_edges = linspace(x_lim(1), x_lim(2), Nx+1);
    y_edges = linspace(y_lim(1), y_lim(2), Ny+1);

    % --- 2D counts ---
    counts = histcounts2(x, y, x_edges, y_edges);

    % --- Gaussian kernel in BIN units ---
    w = max(1, ceil(3*sigma_bins));   % half-width
    [Xk, Yk] = meshgrid(-w:w, -w:w);
    G = exp(-(Xk.^2 + Yk.^2) / (2*sigma_bins^2));
    G = G / sum(G(:));

    % --- smooth counts ---
    dens_s = conv2(counts, G, 'same');

    % --- bin centers ---
    x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
    y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
end
