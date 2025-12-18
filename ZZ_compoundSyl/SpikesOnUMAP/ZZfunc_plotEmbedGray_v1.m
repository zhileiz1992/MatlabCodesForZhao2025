function [ax] = ZZfunc_plotEmbedGray_v1(ax, info, umap_res, i_rd, x_lim, y_lim, marker_size, marker_color, marker_alpha)
%ZZFUNC_PLOTEMBEDGRAY_V1 Summary of this function goes here
%   Detailed explanation goes here
hold(ax, 'on');
u_all = [];
for ii=1:length(i_rd)
    ustart = info.ustart(i_rd(ii));
    uend = info.uend(i_rd(ii));
    u = umap_res(ustart:uend,:);
    scatter(ax, u(:,1), u(:,2), marker_size, 'filled', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', marker_alpha);
    u_all = [u_all; u];
end
xlim(ax, x_lim);
ylim(ax, y_lim);
end

