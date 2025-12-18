function [ax] = ZZfunc_showMMDmatrix_v3(ax, d, clim, tick_loc, tick_label, title_str, plot_boundary, line_loc)
% given a MMD matrix, plot as image
% differ from v2: explicitly set x/y direction
nan_mask = ~isnan(d);
if isempty(clim)
  img = imagesc(ax, d);
else
  img = imagesc(ax, d, clim);
end
img.AlphaData = nan_mask;
colormap(ax, 'gray'); colorbar(ax);
hold(ax, 'on');
set(ax, 'YDir', 'reverse');

if plot_boundary
  for si=1:length(line_loc)
    %     xline(ax, line_loc(si), 'LineStyle', '--', 'LineWidth', 0.25, 'Color', '#99d8c9');
    %     yline(ax, line_loc(si), 'LineStyle', '--', 'LineWidth', 0.25, 'Color', '#99d8c9');
    plot(ax, [1 1],  [0 line_loc(1)], 'LineStyle', '--', 'LineWidth', 1, 'Color', '#99d8c9');
    if si<length(line_loc)
    stop_x = find(isnan(d(:, line_loc(si)+1)));
    if isempty(stop_x)
      xline(ax, line_loc(si), 'LineStyle', '--', 'LineWidth', 1, 'Color', '#99d8c9');
    else
      plot(ax, [line_loc(si) line_loc(si)], [0 stop_x(1)], 'LineStyle', '--', 'LineWidth', 1, 'Color', '#99d8c9');
    end
    end
    stop_y = find(isnan(d(line_loc(si), :)));
    if isempty(stop_y)
      yline(ax, line_loc(si), 'LineStyle', '--', 'LineWidth', 1, 'Color', '#99d8c9');
    else
      plot(ax, [stop_y(end) size(d,2)], [line_loc(si) line_loc(si)], 'LineStyle', '--', 'LineWidth', 1, 'Color', '#99d8c9');
    end
  end
end


axis(ax, 'equal');
axis(ax, 'tight');
ax.XAxisLocation = 'top';
ax.YAxisLocation = 'right';
title(ax, title_str, 'FontSize', 12);
xticks(ax, tick_loc);
xticklabels(ax, tick_label);
yticks(ax, tick_loc);
yticklabels(ax, tick_label);

end

