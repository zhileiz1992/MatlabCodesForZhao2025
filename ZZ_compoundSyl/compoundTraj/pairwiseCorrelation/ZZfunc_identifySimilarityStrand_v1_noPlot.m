function [pass_i] = ZZfunc_identifySimilarityStrand_v1_noPlot(distm, med_filter, thre, min_dur)
% given a similarity matrix, identify similarity strands
distm2 = medfilt2(distm, med_filter);
% thre = 0.35; 
binary_mask = distm2 <= thre;
cc = bwconncomp(binary_mask);
c_all = cc.PixelIdxList;

% is_pass = ZZfunc_calcRange_v1(c_all{1}, size(distm), min_dur)
is_pass = cellfun(@(x) ZZfunc_calcRange_v1(x, size(distm), min_dur), c_all);
pass_i = c_all(is_pass);
end

