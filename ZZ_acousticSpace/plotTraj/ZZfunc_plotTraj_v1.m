function [img, ax] = ZZfunc_plotTraj_v1(ax, embed, col_dict, marg, col_alpha, end_alpha)
% given a embedding data table, plot trajectories of each syllable
% rasterized the scatter plot for easy edit in Illustrator

% determine xy limits
x_lim = [min(embed.umap1)-marg, max(embed.umap1)+marg];
y_lim = [min(embed.umap2)-marg, max(embed.umap2)+marg];
% loop through syllable renditions then plot
syl_id = unique(embed.syl_id, 'sorted');  % get all syllable ids
v_all = {};
for si=1:size(syl_id, 1)
  % for si=1:50
  ss = syl_id{si};
  embed_s = embed(strcmp(embed.syl_id, ss), :);
  % sort by the window index
  embed_s = sortrows(embed_s, 'i_i', 'ascend');
  % get color
  v = embed_s.call_subtype{1};
  v_all = [v_all v];
  col_hex = col_dict.(v);
  col_rgb = sscanf(col_hex(2:end), '%2x%2x%2x', [1 3]) / 255;
  x = embed_s.umap1;
  y = embed_s.umap2;
  % plot as dots on lines
  %   plot(x, y, '-', 'Color', [col_rgb col_alpha], 'LineWidth', 1); hold on;
  scatter(ax, x, y, 15, 'filled', 'MarkerFaceColor', col_rgb, 'MarkerFaceAlpha', col_alpha, 'MarkerEdgeColor', 'none'); hold(ax, 'on');
  % mark the start and end
  %   scatter(x(1), y(1), 40, 'Marker', '^', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
  %   scatter(x(end), y(end), 40, 'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
  %   scatter(x(end), y(end), 40, 'Marker', 's', 'MarkerFaceColor', [0.25 0.25 0.25], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
  scatter(ax, x(1), y(1), 30, 'filled', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold(ax, 'on');
  scatter(ax, x(end), y(end), 30, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold(ax, 'on');
end
xlim(ax, x_lim); ylim(ax, y_lim);
% ax.XTick = [];  % remove x ticks 
%  ax.YTick = [];
ax.TickLength = [0, 0];
% box(ax, 'off');
frame = getframe(ax);
img = frame.cdata;
cla(ax);
xlim(ax, x_lim); ylim(ax, x_lim);
% add the scatter background
image(ax, ax.XLim, ax.YLim, flipud(img));
ax.TickLength = [0.01, 0.025];
ax.XTickMode = 'auto';
ax.YTickMode = 'auto';
ax.XTickLabelMode = 'auto';
ax.YTickLabelMode = 'auto';
% uistack(h_img, 'bottom');
% set(h_img, 'AlphaData', 1);
xlabel(ax, 'UMAP axis 1', 'FontSize', 16);
ylabel(ax, 'UMAP axis 2', 'FontSize', 16);
% Create dummy plot handles for legend
hold(ax, 'on');
syl = unique(v_all, 'sorted');
h_legend = gobjects(length(syl), 1);
for i = 1:length(syl)
  c = col_dict.(syl{i});
  h_legend(i) = plot(ax, nan, nan, 'o', 'MarkerSize', 30, 'Color', c, 'MarkerFaceColor', c, 'DisplayName', syl{i});
end
hold(ax, 'off');
legend(ax, h_legend, 'Location', 'southeast', 'FontSize', 16);
end

