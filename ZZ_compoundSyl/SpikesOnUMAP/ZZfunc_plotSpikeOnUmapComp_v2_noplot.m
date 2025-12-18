function [spike_loc] = ZZfunc_plotSpikeOnUmapComp_v2_noplot(ax, spike, umap_res)
% plot spikes on top of UMAP
spike_loc = [];
for ri=1:size(spike,1)
  mat_loc = spike.mat_loc{ri};
  if ~isempty(mat_loc)
    ustart = spike.ustart(ri);
    uend = spike.uend(ri);
    u = umap_res(ustart:uend,:);
    %       scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', '#a6761d', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.25);
    % plot where the spikes are
    x = u(mat_loc,1);
    y = u(mat_loc,2);
%     scatter(ax, x, y, marker_size, 'filled', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', marker_alpha);
    temp = [x y];
    spike_loc = [spike_loc; temp];
  end
end

end

