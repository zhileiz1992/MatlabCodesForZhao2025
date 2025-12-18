function ax = GPT_plot_mean_std_3d(ax, data3d, line_color, line_width, shade_color, shade_alpha)
if size(data3d, 3) ~= 2
  error('Input data must have size N × T × 2 (last dimension must be 2).');
end

[N, T, ~] = size(data3d);
hold(ax, 'on');

% Extract x and y trajectories
x = squeeze(data3d(:, :, 1));  % N × T
y = squeeze(data3d(:, :, 2));  % N × T

% Compute mean and std
mean_x = mean(x, 1);  % 1 × T
mean_y = mean(y, 1);  % 1 × T
std_x  = std(x, 0, 1); % 1 × T
std_y  = std(y, 0, 1); % 1 × T

% Create ribbon around mean trajectory
upper_x = mean_x + std_x;
lower_x = mean_x - std_x;
upper_y = mean_y + std_y;
lower_y = mean_y - std_y;

% Create closed polygon for shaded area
fill_x = [upper_x, fliplr(lower_x)];
fill_y = [upper_y, fliplr(lower_y)];
fill(ax, fill_x, fill_y, shade_color, ...
  'FaceAlpha', shade_alpha, 'EdgeColor', 'none');

% Plot the mean trajectory
plot(ax, mean_x, mean_y, 'Color', line_color, 'LineWidth', line_width);

% axis(ax, 'equal');
xlabel(ax, 'X');
ylabel(ax, 'Y');
% title(ax, 'Mean 2D Trajectory ± Std');
end


