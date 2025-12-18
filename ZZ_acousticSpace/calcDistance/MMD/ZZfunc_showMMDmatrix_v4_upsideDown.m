function [ax] = ZZfunc_showMMDmatrix_v4_upsideDown(ax, d, clim, tick_loc, tick_label, title_str, plot_boundary, line_loc)
% given a MMD matrix, plot as image
% differ from v3: use the flipped 'gray' colormap
cla(ax);
nan_mask = ~isnan(d);
if isempty(clim)
  img = imagesc(ax, d);
else
  img = imagesc(ax, d, clim);
end
img.AlphaData = nan_mask;
colormap(ax, flipud(gray));
colorbar(ax);
hold(ax, 'on');
set(ax, 'YDir', 'reverse');

if plot_boundary
  for sii=1:length(line_loc)
    for sjj=1:sii
      % plot horizontal line
      if sjj==1; ori_x=1; else; ori_x=line_loc(sjj-1);end
      if sii==1; ori_y=1; else; ori_y=line_loc(sii-1);end
      plot(ax, [ori_x line_loc(sjj)],  [ori_y ori_y], 'LineStyle', '--', 'LineWidth', 1, 'Color', '#99d8c9');
      % plot vertical line
      plot(ax, [line_loc(sjj) line_loc(sjj)],  [ori_y line_loc(sii)], 'LineStyle', '--', 'LineWidth', 1, 'Color', '#99d8c9');
  
    end
    
  end
end


axis(ax, 'equal');
axis(ax, 'tight');
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
title(ax, title_str, 'FontSize', 12);
xticks(ax, tick_loc);
xticklabels(ax, tick_label);
yticks(ax, tick_loc);
yticklabels(ax, tick_label);

end

