function [ax, y_start] = ZZfunc_plotRasterRelTime_v1(ax, d_plot, y_start, marker, marker_size, marker_color, marker_alpha, x_lim, add_bar, bar_color, add_label, label_str, label_offset, offset_color, sort_dur)
% plot the raster of spikes, no time-wrapping
% differ from v2: add option to mark syllable offset and sort renditions by duration

hold(ax, 'on');

if sort_dur
  [~, sort_i] = sort(d_plot.dur, 'ascend');
  d_plot = d_plot(sort_i,:);
end

y_start_old = y_start; 
for ri=1:size(d_plot, 1)
  rel_ori = d_plot.seg_start_ori(ri) - d_plot.seg_start(ri);
  dur_pt = d_plot.seg_end_ori(ri) - d_plot.seg_start_ori(ri);
  if ~isempty(d_plot.spike_i{ri})
    % calculate relative position of syllable onset
    rel_pos = d_plot.spike_i{ri} - rel_ori;
    spike_t = rel_pos / dur_pt;
    scatter(ax, spike_t, zeros(1, length(spike_t))+y_start, marker_size, 'Marker', marker, 'MarkerFaceColor', marker_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeColor', 'none');
  end
  if ~isempty(offset_color)
    scatter(ax, 1, y_start, marker_size, 'Marker', marker, 'MarkerFaceColor', offset_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeColor', 'none');
  end
  y_start = y_start + 1;
end

% add a vertical bar if specified
if add_bar
  plot(ax, [x_lim(1) x_lim(1)], [y_start_old y_start], 'Color', bar_color, 'LineWidth', 4, 'LineStyle', '-');
end

if add_label
  text(ax, x_lim(1)-label_offset, (y_start_old+y_start)/2, label_str, 'FontSize', 12, 'Color', bar_color);
end



xlim(ax, x_lim);



