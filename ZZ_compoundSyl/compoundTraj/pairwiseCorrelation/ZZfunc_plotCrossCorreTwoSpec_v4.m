function [fig, axes_all, distm] = ZZfunc_plotCrossCorreTwoSpec_v4(d1, d2, spec1, sound2, clim2, cmap1, cmap2, flim, clim_distm)
% differ from v3: set a frequency limits to remove blank region

% read wav then plot spectrogram
fs = 20000;

% close all;
fig_size = [10 10 900 900];
left_margin = 0.05;       % Space from left edge of figure
bottom_margin = 0.05;     % Space from bottom edge of figure
left_col_width = 0.16;    % Width of left column (for B)
right_col_width = 0.16;   % Width of right column (for A and C; A and C will have same width)
top_row_height = 0.16;    % Height of top row (for A)
% bottom_row_height = 0.5; 
bottom_row_height = right_col_width * (size(d2,1) / size(d1,1));

[fig, axes_all] = ZZ_emptyCrossCorrePlot_v1(fig_size, left_margin, bottom_margin, left_col_width, right_col_width, top_row_height, bottom_row_height);

% imagesc(axes_all(2), permute(spec2, [2, 1, 3]));
[ax2, ~, ~,  ~, ~,  t2] = showAudioSpectrogramZZ_flexible_v1_flip(sound2, fs, axes_all(2), flim, clim2, 256, 256, 236);
axis(axes_all(2), 'off');

% freq_i = find((f>=flim(1)) & (f<=flim(2)));
imagesc(axes_all(1), spec1, [0.15 max(spec1(:))]); 
set(axes_all(1), 'YDir', 'normal'); colormap(axes_all(1), cmap1);
axis(axes_all(1), 'off');


  
% computate distance between every two time points of the syllable
% d1 = comp.dvae{idx_rd(1)};
% d2 = comp.dvae{idx_rd(2)};
% distm = pdist2(d1, d2, 'euclidean');
distm = pdist2(d1, d2, 'cosine');
distm = distm'; 
% decompose 
% [B, C, u, v] = mmd_matrix_decompose(distm);
% imagesc(axes_all(3), t1, t2, distm');
imagesc(axes_all(3), distm, clim_distm);
colormap(axes_all(3), cmap2);
axis(axes_all(3), 'off');
end

