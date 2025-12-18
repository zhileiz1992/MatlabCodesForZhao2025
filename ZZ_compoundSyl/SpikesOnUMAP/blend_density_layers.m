function [RGB_out, A_out] = blend_density_layers(density_all, cmap_array, v_thre, clim_array, method, alpha_max)
    if nargin < 5 || isempty(method), method = 'weighted'; end
    if nargin < 6 || isempty(alpha_max), alpha_max = 0.9; end

    K = numel(density_all);
    if nargin < 4 || isempty(clim_array), clim_array = cell(1,K); end

    [M0,N0] = size(density_all{1});
    RGB_stack = zeros(M0, N0, 3, K);
    A_stack   = zeros(M0, N0, K);

    for k = 1:K
        Z = density_all{k};
        if ~isequal(size(Z), [M0,N0]), Z = imresize(Z, [M0 N0], 'bicubic'); end

        cmap_k = resolve_cmap(cmap_array{k});
        if ~isempty(clim_array{k}), clim_k = clim_array{k};
        else
            znz = Z(~isnan(Z));
            if isempty(znz), clim_k = [0 1];
            else
                vmin = min(znz); vmax = max(znz);
                if vmin==vmax, clim_k = vmin + [-0.5 0.5]; else, clim_k = [vmin vmax]; end
            end
        end

        RGB_stack(:,:,:,k) = scalar_to_rgb(Z, cmap_k, clim_k);

        t = (Z - v_thre) ./ max(eps, (clim_k(2) - v_thre));
        t = max(0, min(1, t));
        a = t.^2 .* (3 - 2*t);
        a(isnan(Z)) = 0;
        A_stack(:,:,k) = alpha_max * a;
    end

    % ---- Correct expansion (key fix) ----
    A3 = repmat(reshape(A_stack, [M0, N0, 1, K]), [1 1 3 1]);  % M×N×3×K

    switch lower(method)
        case 'weighted'
            num = sum(RGB_stack .* A3, 4);         % M×N×3
            W   = sum(A_stack, 3);                 % M×N
            Wsafe = W; Wsafe(Wsafe==0) = 1;
            RGB_out = bsxfun(@rdivide, num, Wsafe);
            A_out   = 1 - prod(1 - A_stack, 3);

        case 'screen'
            one_minus = 1 - A3 .* RGB_stack;       % M×N×3×K
            RGB_out   = 1 - prod(one_minus, 4);    % M×N×3
            A_out     = 1 - prod(1 - A_stack, 3);

        otherwise
            error('Unknown method. Use ''weighted'' or ''screen''.');
    end

    RGB_out = max(0, min(1, RGB_out));
    A_out   = max(0, min(1, A_out));
end

function cmap = resolve_cmap(spec)
    if isnumeric(spec)
        if size(spec,2)~=3, error('Numeric colormap must be N×3.'); end
        cmap = max(0, min(1, spec));
    elseif isa(spec,'function_handle')
        cmap = spec(256);
    elseif ischar(spec) || isstring(spec)
        cmap = feval(char(spec), 256);
    else
        error('Unsupported colormap specification.');
    end
end

function RGB = scalar_to_rgb(Z, cmap, clim)
    Zs = (Z - clim(1)) ./ max(eps, diff(clim));
    Zs = max(0, min(1, Zs));
    idx = 1 + floor(Zs * (size(cmap,1)-1));
    idx(isnan(Z)) = 1;
    RGB = reshape(cmap(idx,:), [size(Z,1), size(Z,2), 3]);
end
