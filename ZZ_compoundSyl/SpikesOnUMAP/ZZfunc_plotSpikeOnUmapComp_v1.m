function [ax] = ZZfunc_plotSpikeOnUmapComp_v1(ax, spike, umap_res, marker_size, marker_color, marker_alpha)
% plot spikes on top of UMAP
for ri=1:size(spike,1)
  mat_loc = spike.mat_loc{ri};
  if ~isempty(mat_loc)
    ustart = spike.ustart(ri);
    uend = spike.uend(ri);
    u = umap_res(ustart:uend,:);
    %       scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', '#a6761d', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.25);
    % plot where the spikes are
    scatter(ax, u(mat_loc,1), u(mat_loc,2), marker_size, 'filled', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', marker_alpha);
  end
end

end

