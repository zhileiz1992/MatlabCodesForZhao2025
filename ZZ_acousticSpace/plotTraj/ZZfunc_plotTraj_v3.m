function [img, ax] = ZZfunc_plotTraj_v3(ax, embed, info, num_plot, batches, r_seed, col_dict, marg, col_alpha, end_alpha, x_lim, y_lim)
% given a embedding data table, plot trajectories of each syllable
% rasterized the scatter plot for easy edit in Illustrator
% differ from v2: use the UMAP trained on all syllable renditions, embed is now just an array
% meta info is in the info table; add option to sample syllables from certain batch of sorting

% determine xy limits if not set
if isempty(x_lim); x_lim = [min(embed(:,1))-marg, max(embed(:,1))+marg]; end
if isempty(y_lim); y_lim = [min(embed(:,2))-marg, max(embed(:,2))+marg]; end

% sample syllables
syl_type = unique(info.call_subtype, 'sorted');
syl_select = {};
rng(r_seed);
for ii=1:size(syl_type,1)
  i_this = find((strcmp(info.call_subtype, syl_type{ii})) & (ismember(info.batch, batches)));
  to_sample = min([num_plot, length(i_this)]);
  i_rd = randsample(i_this, to_sample);
  syl_select = [syl_select; info.syl_ID(i_rd)];
end
% shuffle the order, so one syllable types are mixed
syl_select = syl_select(randperm(numel(syl_select)));

% loop through syllable renditions then plot
v_all = {};
for si=1:size(syl_select, 1)
  % for si=1:50
  ss = syl_select{si};
  this_meta = info(strcmp(info.syl_ID, ss),:);
  ii_start = this_meta.count_start + 1;  % note python is 0-based
  ii_end = this_meta.count_end;
  this_umap = embed(ii_start:ii_end,:);
  % get color
  v = this_meta.call_subtype{1};
  v_all = [v_all v];
  col_hex = col_dict.(v);
  col_rgb = sscanf(col_hex(2:end), '%2x%2x%2x', [1 3]) / 255;
  x = this_umap(:,1);
  y = this_umap(:,2);
  % plot as dots on lines
  %   plot(x, y, '-', 'Color', [col_rgb col_alpha], 'LineWidth', 1); hold on;
  scatter(ax, x, y, 10, 'filled', 'MarkerFaceColor', col_rgb, 'MarkerFaceAlpha', col_alpha, 'MarkerEdgeColor', 'none'); hold(ax, 'on');
  % mark the start and end
  %   scatter(x(1), y(1), 40, 'Marker', '^', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
  %   scatter(x(end), y(end), 40, 'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
  %   scatter(x(end), y(end), 40, 'Marker', 's', 'MarkerFaceColor', [0.25 0.25 0.25], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
  scatter(ax, x(1), y(1), 20, 'filled', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold(ax, 'on');
  scatter(ax, x(end), y(end), 20, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold(ax, 'on');
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

