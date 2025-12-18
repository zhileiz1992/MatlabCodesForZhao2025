function overlay_density_images(density_all, cmap_array, v_thre, clim_array)
%OVERLAY_DENSITY_IMAGES Overlay multiple 2-D arrays with distinct colormaps.
%   OVERLAY_DENSITY_IMAGES(density_all, cmap_array, v_thre, clim_array)
%   - density_all : 1xK or Kx1 cell; each cell is an MxN numeric matrix (same size)
%   - cmap_array  : 1xK or Kx1 cell; each cell is a colormap spec for that layer:
%                   * an Nx3 numeric colormap (e.g., parula(256))
%                   * a function handle (e.g., @parula)
%                   * a char/string name (e.g., 'parula')
%   - v_thre      : scalar threshold; values < v_thre become transparent
%   - clim_array  : (optional) 1xK or Kx1 cell; each cell is [vmin vmax] for that layer.
%                   If omitted/empty for a layer, it uses [min(Z(:)) max(Z(:))].
%
%   Each layer is drawn in order density_all{1}, density_all{2}, ... (topmost last).

    assert(iscell(density_all) && ~isempty(density_all), 'density_all must be a non-empty cell array.');
    K = numel(density_all);

    if nargin < 4, clim_array = []; end
    if isempty(cmap_array) || numel(cmap_array) ~= K
        error('cmap_array must be a cell array with the same length as density_all.');
    end
    if ~isempty(clim_array) && numel(clim_array) ~= K
        error('clim_array must match density_all in length, or be empty.');
    end

    % Prepare figure/axes
    fig = figure('Renderer','opengl'); %#ok<NASGU>  % OpenGL supports per-pixel alpha
    ax  = axes; hold(ax,'on'); axis(ax,'image'); set(ax,'YDir','normal');
    set(ax,'Color','w');  % background seen through transparent pixels

    % Draw each layer
    for k = 1:K
        Z = density_all{k};
        validateattributes(Z, {'numeric'}, {'2d','real','nonempty'});
        if k == 1
            [M,N] = size(Z);
        else
            assert(isequal(size(Z), [M,N]), 'All matrices must have identical size.');
        end

        % Determine colormap for this layer
        cmap_k = resolve_cmap(cmap_array{k});

        % Determine color limits
        if ~isempty(clim_array) && ~isempty(clim_array{k})
            clim_k = clim_array{k};
        else
            % default to the data range of this layer (ignoring NaNs)
            znz = Z(~isnan(Z));
            if isempty(znz)
                clim_k = [0 1];  % fallback
            else
                clim_k = [min(znz), max(znz)];
                if clim_k(1) == clim_k(2)
                    clim_k = clim_k + [-0.5 0.5]; % avoid zero span
                end
            end
        end

        % Map scalar Z to truecolor RGB using the chosen colormap/CLim
        RGB = scalar_to_rgb(Z, cmap_k, clim_k);

        % Alpha mask: transparent where Z < v_thre or isnan(Z)
        alpha = ~(Z < v_thre | isnan(Z));  % logical -> converted to 0/1 automatically

        % Draw this layer
        h = image(ax, RGB, 'AlphaData', alpha);
%         set(h, 'AlphaData', 0.75);
    end

    % Axes cosmetics
    box(ax,'on'); axis(ax,'tight');
end

% --- Helpers --------------------------------------------------------------

function cmap = resolve_cmap(spec)
% Accepts:
%   - numeric Nx3 colormap
%   - function handle (e.g., @parula)
%   - char/string name (e.g., 'parula')
    if isnumeric(spec)
        validateattributes(spec, {'numeric'}, {'ncols',3,'>=',0,'<=',1});
        cmap = spec;
    elseif isa(spec,'function_handle')
        cmap = spec(256);
    elseif ischar(spec) || isstring(spec)
        cmap = feval(char(spec), 256);
    else
        error('Unsupported colormap specification for one layer.');
    end
end

function RGB = scalar_to_rgb(Z, cmap, clim)
% Map Z (MxN) to RGB (MxNx3) via given colormap and [vmin vmax].
    Zs = (Z - clim(1)) ./ max(eps, (clim(2) - clim(1)));  % normalize
    Zs = max(0, min(1, Zs));                              % clip to [0,1]
    idx = 1 + floor(Zs * (size(cmap,1)-1));               % 1..size(cmap,1)
    % Handle NaNs: keep them, will be transparent via alpha
    idx(isnan(Z)) = 1;                                    % any valid row
    RGB = reshape(cmap(idx, :), [size(Z,1), size(Z,2), 3]);
end
