function [ax] = ZZfunc_showMMDmatrix_v1(ax, d, clim, tick_loc, tick_label, title_str, plot_boundary, line_loc)
% given a MMD matrix, plot as image
if isempty(clim)
  imagesc(ax, d);
else
  imagesc(ax, d, clim);
end
colormap(ax, 'gray'); colorbar(ax);

if plot_boundary
  for si=1:length(line_loc)
    xline(ax, line_loc(si), 'LineStyle', '--', 'LineWidth', 0.25, 'Color', '#99d8c9');
    yline(ax, line_loc(si), 'LineStyle', '--', 'LineWidth', 0.25, 'Color', '#99d8c9');
  end
end
  
title(ax, title_str, 'FontSize', 12);
xticks(ax, tick_loc);
xticklabels(ax, tick_label);
yticks(ax, tick_loc);
yticklabels(ax, tick_label);

end

