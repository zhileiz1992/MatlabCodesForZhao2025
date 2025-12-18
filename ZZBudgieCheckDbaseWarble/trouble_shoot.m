% subplots with different heights
plot_row = 2; 
plot_col = 10;
width = 0.075; 
heights = [0.1 0.1 0.02];

close all;
figure;
x = randn(100); 
y = randn(100);
ii = 1; 
ax1 = subplot(3, 2, 1); plot(ax1, x, y); 
% calculate the heights
xgap = 0.004; ygap = 0.004; 
left_margin = 0.1; 
bottom_margin = 0.1; 
x_pos = left_margin; 
y_pos = 1 - heights(1) - ygap; 
set(ax1, 'Position', [x_pos y_pos width heights(1)]);

% 2nd and 3rd plot
ax2 = subplot(3, 2, 3); plot(ax2, x, y); 
set(ax2, 'Position', [x_pos y_pos-heights(2)-ygap width heights(2)]);


ax3 = subplot(3, 2, 5); plot(ax3, x, y); 
set(ax3, 'Position', [x_pos y_pos-heights(2)-heights(3)-ygap*2 width heights(3)]);