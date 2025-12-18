function [fig, axes_all, distm] = ZZfunc_plotCrossCorreTwoSpec_v3_callLeft(d1, d2, spec1, sound1, clim2, cmap1, cmap2)
% differ from v3: show call on the left panel

% read wav then plot spectrogram
fs = 20000;

% close all;
fig_size = [10 10 900 900];
left_margin = 0.05;       % Space from left edge of figure
bottom_margin = 0.05;     % Space from bottom edge of figure
left_col_width = 0.2;    % Width of left column (for B)
% right_col_width = 0.7;   % Width of right column (for A and C; A and C will have same width)
top_row_height = 0.2;    % Height of top row (for A)
% bottom_row_height = 0.5; 
% bottom_row_height = right_col_width * (size(d2,1) / size(d1,1));
bottom_row_height = 0.2;
right_col_width = bottom_row_height * (size(d1,1) / size(d2,1));

[fig, axes_all] = ZZ_emptyCrossCorrePlot_v1(fig_size, left_margin, bottom_margin, left_col_width, right_col_width, top_row_height, bottom_row_height);
% imagesc(axes_all(1), spec1, [0.15 max(spec1(:))]);
[ax1, ~, ~,  ~, ~,  t2] = showAudioSpectrogramZZ_flexible_v1(sound1, fs, axes_all(1), [250 7500], clim2, 256, 256, 236);
% set(axes_all(1), 'YDir', 'normal'); colormap(axes_all(1), cmap1);
axis(axes_all(1), 'off');

% imagesc(axes_all(2), permute(spec2, [2, 1, 3]));
% [ax2, ~, ~,  ~, ~,  t2] = showAudioSpectrogramZZ_flexible_v1_flip(sound2, fs, axes_all(2), [250 7500], clim2, 256, 256, 236);
imagesc(axes_all(2), spec1', [0.15 max(spec1(:))]);
% set(axes_all(2), 'YDir', 'normal'); 
colormap(axes_all(2), cmap1);
axis(axes_all(2), 'off');
  
% computate distance between every two time points of the syllable
% d1 = comp.dvae{idx_rd(1)};
% d2 = comp.dvae{idx_rd(2)};
% distm = pdist2(d1, d2, 'euclidean');
distm = pdist2(d1, d2, 'cosine');
distm = distm'; 
% decompose 
% [B, C, u, v] = mmd_matrix_decompose(distm);
% imagesc(axes_all(3), t1, t2, distm');
imagesc(axes_all(3), distm, [0 0.5]);
colormap(axes_all(3), cmap2);
axis(axes_all(3), 'off');
end

