function [ax] = ZZfunc_showMMDmatrix_v5(ax, d, clim, tick_loc, tick_label, title_str, plot_boundary, line_loc, line_col, line_style, line_width, c_map, updown)
% given a MMD matrix, plot as image
% differ from v4: option to set the colormap and plot in upside-down
cla(ax);
nan_mask = ~isnan(d);
if isempty(clim)
  img = imagesc(ax, d);
else
  img = imagesc(ax, d, clim);
end
img.AlphaData = nan_mask;
colormap(ax, c_map);
colorbar(ax);
hold(ax, 'on');
set(ax, 'YDir', 'reverse');

if plot_boundary
  if ~updown
    for sii=1:length(line_loc)
      for sjj=sii:length(line_loc)
        % plot horizontal line
        if sjj==1; ori_x=1; else; ori_x=line_loc(sjj-1);end
        if sii==1; ori_y=1; else; ori_y=line_loc(sii-1);end
        plot(ax, [ori_x line_loc(sjj)],  [line_loc(sii) line_loc(sii)], 'LineStyle', line_style, 'LineWidth', line_width, 'Color', line_col);
        % plot vertical line
        plot(ax, [ori_x ori_x], [ori_y line_loc(sii)], 'LineStyle', line_style, 'LineWidth', line_width, 'Color', line_col);
        
      end
    end
  else
    for sii=1:length(line_loc)
      for sjj=1:sii
        % plot horizontal line
        if sjj==1; ori_x=1; else; ori_x=line_loc(sjj-1);end
        if sii==1; ori_y=1; else; ori_y=line_loc(sii-1);end
        plot(ax, [ori_x line_loc(sjj)],  [ori_y ori_y], 'LineStyle', line_style, 'LineWidth', line_width, 'Color', line_col);
        % plot vertical line
        plot(ax, [line_loc(sjj) line_loc(sjj)],  [ori_y line_loc(sii)], 'LineStyle', line_style, 'LineWidth', line_width, 'Color', line_col);
      end    
    end
  end
end

axis(ax, 'equal');
axis(ax, 'tight');
if ~updown
  ax.XAxisLocation = 'top';
  ax.YAxisLocation = 'right';
else
  ax.XAxisLocation = 'bottom';
  ax.YAxisLocation = 'left';
end
title(ax, title_str, 'FontSize', 12);
xticks(ax, tick_loc);
xticklabels(ax, tick_label);
yticks(ax, tick_loc);
yticklabels(ax, tick_label);

end

