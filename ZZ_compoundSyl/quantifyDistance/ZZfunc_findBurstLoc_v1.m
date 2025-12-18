function [burst_loc] = ZZfunc_findBurstLoc_v1(spike_loc, max_gap, method)
%ZZFUNC_FINDBURSTLOC_V1 Summary of this function goes here

% if only one spike
if length(spike_loc)==1
  burst_loc = spike_loc;
else
  bursts = {};
  % Calculate inter-spike intervals
  intervals = diff(spike_loc);

  % Find burst boundaries (where intervals >= max_interval)
  burst_end_idx = find(intervals >= max_gap);
  burst_start_idx = [1; burst_end_idx + 1];
  burst_end_idx = [burst_end_idx; length(spike_loc)];

  % Group spikes into bursts
  for i = 1:length(burst_start_idx)
      bursts{i} = spike_loc(burst_start_idx(i):burst_end_idx(i));
  end
  
  % define the location of the burst
  if strcmp(method, 'first')
    burst_loc = cellfun(@(x) x(1), bursts);
  elseif strcmp(method, 'center')
    burst_loc = cellfun(@(x) floor(mean(x)), bursts);
  end
  burst_loc = burst_loc';
end


end

