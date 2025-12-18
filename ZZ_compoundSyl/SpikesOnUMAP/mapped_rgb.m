function RGB = mapped_rgb(Z, cmap, clim)
    if nargin < 3 || isempty(clim), clim = [min(Z(:)) max(Z(:))]; end
    Zs = (Z - clim(1)) / max(eps, diff(clim));          % normalize to [0,1]
    Zs = max(0, min(1, Zs));
    idx = 1 + floor(Zs * (size(cmap,1)-1));             % 1..N
    RGB = reshape(cmap(idx,:), [size(Z,1) size(Z,2) 3]);% HxWx3
end
