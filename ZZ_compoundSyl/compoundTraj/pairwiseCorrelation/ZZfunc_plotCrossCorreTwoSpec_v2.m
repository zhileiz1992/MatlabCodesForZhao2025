function [fig, axes_all, distm, t1, t2] = ZZfunc_plotCrossCorreTwoSpec_v2(comp, idx_rd)
% differ from v1: when extract sound, pad half window
% idx_rd = randsample(1:size(comp,1), 2);

% read wav then plot spectrogram
fs = 20000;
pad_pre = 0.016; pad_pre_pt = floor(fs*pad_pre);  % pad equal amount as in VAE sliding
pad_post = 0.016; pad_post_pt = floor(fs*pad_post);
sounds = cell(2, 1);
for ii=1:length(idx_rd)
  idx = idx_rd(ii);
  [signal, fs] = audioread(comp.fn_wav{idx});
  i_start = max([1 comp.istart(idx)-pad_pre_pt]);
  i_end = min([size(signal,1) comp.iend(idx)+pad_post_pt]);
  sounds{ii} = signal(i_start:i_end);
end
  
% close all;
fig_size = [10 10 800 800];
left_margin = 0.05;       % Space from left edge of figure
bottom_margin = 0.05;     % Space from bottom edge of figure
left_col_width = 0.16;    % Width of left column (for B)
% right_col_width = 0.4;   % Width of right column (for A and C; A and C will have same width)
top_row_height = 0.16;    % Height of top row (for A)
% determine bottom height proportionally
[power1, ~, ~, ~, ~, ~] = getAudioSpectrogramZZ_flexible_v1(sounds{1}, fs, 256, 256, 236, [250 7500], [12 23]);
[power2, ~, ~, ~, ~, ~] = getAudioSpectrogramZZ_flexible_v1(sounds{2}, fs, 256, 256, 236, [250 7500], [12 23]);
% bottom_row_height = 0.5; 
if size(power1,2)>size(power2,2)
  right_col_width = 0.54;
  bottom_row_height = right_col_width * size(power2,2) / size(power1,2);
else
  bottom_row_height = 0.54; 
  right_col_width = bottom_row_height * size(power1,2) / size(power2,2);
end

[fig, axes_all] = ZZ_emptyCrossCorrePlot_v1(fig_size, left_margin, bottom_margin, left_col_width, right_col_width, top_row_height, bottom_row_height);
[ax1, ~, ~,  ~, ~,  t1] = showAudioSpectrogramZZ_flexible_v1(sounds{1}, fs, axes_all(1), [250 7500], [12 23], 256, 256, 236);
axis(axes_all(1), 'off');
[ax2, ~, ~,  ~, ~,  t2] = showAudioSpectrogramZZ_flexible_v1_flip(sounds{2}, fs, axes_all(2), [250 7500], [12 23], 256, 256, 236);
axis(axes_all(2), 'off');
  
% computate distance between every two time points of the syllable
d1 = comp.dvae{idx_rd(1)};
d2 = comp.dvae{idx_rd(2)};
% distm = pdist2(d1, d2, 'euclidean');
distm = pdist2(d1, d2, 'cosine');
distm = distm'; 
% decompose 
% [B, C, u, v] = mmd_matrix_decompose(distm);
% imagesc(axes_all(3), t1, t2, distm');
imagesc(axes_all(3), t1, t2, distm, [0 0.5]);
colormap(axes_all(3), gray);
axis(axes_all(3), 'off');
end

