function [p] = ZZfunc_runSignRank_v1(x)
%ZZFUNC_RUNSIGNRANK_V1 Summary of this function goes here
%   Detailed explanation goes here
x2 = x(~isnan(x));
if length(x2)<2
  p = nan;
else
  p = signrank(x2, 0, 'tail', 'right');
end

end

