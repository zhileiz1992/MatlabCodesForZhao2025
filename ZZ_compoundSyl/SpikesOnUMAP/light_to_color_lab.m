function cmap = light_to_color_lab(baseRGB, n, white_frac, t_light)
%WHITE_LIGHT_TO_COLOR  Colormap from white -> light tint -> base color.
%   cmap = WHITE_LIGHT_TO_COLOR(baseRGB, n, white_frac, t_light)
%   - baseRGB    : 1x3 target RGB in [0,1], e.g., [0.89 0.10 0.11]
%   - n          : number of colors (default 256)
%   - white_frac : how "white" the light tint is (0..1), default 0.90
%                  (light = (1-white_frac)*base + white_frac*white)
%   - t_light    : position of the light tint in the ramp (0..1), default 0.35
%
%   Resulting ramp smoothly goes: white -> light tint -> base color.

    if nargin < 2 || isempty(n), n = 256; end
    if nargin < 3 || isempty(white_frac), white_frac = 0.90; end
    if nargin < 4 || isempty(t_light), t_light = 0.35; end

    % sanitize inputs
    base = max(0, min(1, baseRGB(:)'));         % 1x3 in [0,1]
    t_light = max(0, min(1, t_light));
    white = [1 1 1];
    light = (1 - white_frac) * base + white_frac * white;

    % anchor positions and colors
    T = [0, t_light, 1];
    C = [white; light; base];                   % 3x3

    % interpolate
    t = linspace(0, 1, n).';
    r = interp1(T, C(:,1), t, 'linear');
    g = interp1(T, C(:,2), t, 'linear');
    b = interp1(T, C(:,3), t, 'linear');
    cmap = [r g b];

    % numeric safety
    cmap = max(0, min(1, cmap));
end
