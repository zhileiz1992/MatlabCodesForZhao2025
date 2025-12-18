function [dist] = ZZfunc_pairwiseDist_v1(d_this, num_sample, metric)
% Calculate pairwise distance between rows in d_this
ii = randsample(1:size(d_this,1), num_sample, true);
jj = randsample(1:size(d_this,1), num_sample, true);
idx = find(ii~=jj);
ii = ii(idx);
jj = jj(idx);
dist_all = [];
for i=1:length(ii)
  d1 = d_this(ii(i),:);
  d2 = d_this(jj(i),:);
  temp = pdist2(d1, d2, metric);
  dist_all = [dist_all temp];
end
dist = mean(dist_all);
end

