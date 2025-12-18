function [fig, axs] = ZZfunc_UpSetPlot_v1(X, labels)
% ZZfunc_UpSetPlot - Generate an UpSet plot for binary matrix X
% Inputs:
%   X      : n × m binary matrix (rows = samples, columns = features)
%   labels : 1 × m cell array of strings (names of features)

if ~islogical(X)
    X = logical(X);
end

[n, m] = size(X);
if length(labels) ~= m
    error('Length of labels must match number of columns in X.');
end

% Unique binary combinations and counts
[unique_combos, ~, ic] = unique(X, 'rows');
counts = accumarray(ic, 1);
[sorted_counts, sort_idx] = sort(counts, 'descend');
sorted_combos = unique_combos(sort_idx, :);
num_patterns = size(sorted_combos, 1);

% --- Plotting ---
[fig, axs] = generatePanelGrid_v2(2, 1, [0.45;0.15], [0.02], [0.05;0.05], [0.15;0.05], 0.05, [0;0], [10 10 600 550]);

% Top barplot
% subplot(2,1,1);
ax1 = axs(1);
bar(ax1, sorted_counts, 'FaceColor', [0.3 0.3 0.3]);
x = 1:length(sorted_counts);
for i = 1:length(sorted_counts)
  text(ax1, x(i), sorted_counts(i) + 2, num2str(sorted_counts(i)), 'HorizontalAlignment', 'center', 'FontSize', 9);
end
ylabel(ax1, 'No. of neurons');
xlim(ax1, [0.5 num_patterns+0.5]);
set(ax1, 'XTick', []);
box(ax1, 'off');

% Matrix plot
% subplot(2,1,2);
ax2 = axs(2);
hold(ax2, 'on');

for i = 1:num_patterns
    active = find(sorted_combos(i,:));
    y_idx = m - active + 1;  % flip y-axis
    plot(ax2, i * ones(size(active)), y_idx, 'ko', ...
         'MarkerFaceColor', 'k', 'MarkerSize', 9);
    if numel(active) > 1
        plot(ax2, [i i], [min(y_idx) max(y_idx)], 'k-', 'LineWidth', 1);
    end
end

% Axes formatting
ylim(ax2, [0.5 m + 0.5]);
xlim(ax2, [0.5 num_patterns + 0.5]);
yticks(ax2, 1:m);
yticklabels(ax2, flip(labels));
xlabel(ax2, 'Property combinations');
set(ax2,'YDir','normal');
box(ax2, 'off');
ax2.FontSize=10;

end
