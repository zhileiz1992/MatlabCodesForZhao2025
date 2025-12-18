function [ax, y_start] = ZZfunc_plotRasterNoWrap_v1(ax, d_plot, y_start, fs, marker, marker_size, marker_color, marker_alpha, x_lim, add_bar, bar_color, add_label, label_str)
% plot the raster of spikes, no time-wrapping
hold(ax, 'on');
y_start_old = y_start; 
for ri=1:size(d_plot, 1)
  if ~isempty(d_plot.spike_i{ri})
    % calculate relative position of syllable onset
    rel_ori = d_plot.seg_start_ori(ri) - d_plot.seg_start(ri);
    rel_pos = d_plot.spike_i{ri} - rel_ori;
    spike_t = rel_pos / fs;
    scatter(ax, spike_t, zeros(1, length(spike_t))+y_start, marker_size, 'Marker', marker, 'MarkerFaceColor', marker_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeColor', 'none');
  end
  y_start = y_start + 1;
end

% add a vertical bar if specified
if add_bar
  plot(ax, [x_lim(end) x_lim(end)], [y_start_old y_start], 'Color', bar_color, 'LineWidth', 4, 'LineStyle', '-');
end

if add_label
  text(ax, x_lim(end)+0.03, (y_start_old+y_start)/2, label_str, 'FontSize', 12, 'Color', bar_color);
end

xlim(ax, x_lim);



