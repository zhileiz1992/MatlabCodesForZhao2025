function [fig2, ax_all2, pass_i, newMask] = ZZfunc_identifySimilarityStrand_v1(distm, med_filter, thre, min_dur)
% given a similarity matrix, identify similarity strands
distm2 = medfilt2(distm, med_filter);
% thre = 0.35; 
binary_mask = distm2 <= thre;
cc = bwconncomp(binary_mask);
% close all;
[fig2, ax_all2] = generatePanelGrid_v2(1, 4, [0.7], [], [0.05;0.05], [0.05;0.05], 0.05, [0], [10 10 2000 600]);
% plot the original
imagesc(ax_all2(1), distm, [0 0.5]);
colormap(ax_all2(1), gray); title(ax_all2(1), 'Original matrix'); 
% plot the smooth
imagesc(ax_all2(2), distm2, [0 0.5]);
colormap(ax_all2(2), gray); title(ax_all2(2), 'Median filtered'); 
% plot the hard mask
imagesc(ax_all2(3), binary_mask);
title(ax_all2(3), sprintf('Passed threshold %.2f', thre)); 
% filter out based on min dur
% min_dur = 25;  % at least 25 frames on each side, i.e. 25 ms
c_all = cc.PixelIdxList;
pass_i = {};
count = 0;
for ci=1:size(c_all, 2)
  linear_i = c_all{ci};
  % convert to xy index
  [xx, yy] = ind2sub(size(distm2), linear_i);
  x_range = max(xx) - min(xx); 
  y_range = max(yy) - min(yy); 
  if x_range>=min_dur &&  y_range>=min_dur
    count = count + 1;
    pass_i{count} = linear_i;
  end
end
% plot those that pass the size limit
newMask = zeros(size(distm2));
for ci=1:size(pass_i, 2)
  newMask(pass_i{ci}) = 1; 
end
imagesc(ax_all2(4), newMask);
title(ax_all2(4), sprintf('Passed duration limit %d ms', min_dur)); 
end

