function [dist] = ZZfunc_pairwiseDist_v2_twoData(d_this1, d_this2, num_sample, metric)
% Calculate pairwise distance between rows in d_this
% differ from v2: calcualte the distance between two data arrays
ii = randsample(1:size(d_this1,1), num_sample, true);
jj = randsample(1:size(d_this2,1), num_sample, true);
dist_all = [];
for i=1:length(ii)
  d1 = d_this1(ii(i),:);
  d2 = d_this2(jj(i),:);
  temp = pdist2(d1, d2, metric);
  dist_all = [dist_all temp];
end
dist = nanmean(dist_all);
end

