function [onsets, offsets] = ZZfunc_GetProbOnsetOffset_v2(prob_amp, dt, param)
% Identify the onsets and offset of peaks in the flatness curve

% get peaks in flatness
[~,idxs,~,~] = findpeaks(prob_amp,'MinPeakDistance',param.gapSize,'MinPeakHeight',param.thresholdFlatness, 'MinPeakProminence', param.minProminence);

% extend from the peak location find the onsets and offset
% either reaches 20% of the peak height
% per = param.extendPercent; 
extThre = param.extendThreshold;
onsets = [];
offsets = [];
if ~isempty(idxs)
  % first pass, extend until reach criteria
  for ii=1:length(idxs)
    pidx = idxs(ii);
    % extend to the left
    left_i = pidx-1;
    while (left_i>0) && (prob_amp(left_i)>extThre)
      left_i = left_i-1;
    end
    % extend to the right
    right_i = pidx+1;
    while (right_i<length(prob_amp)) && (prob_amp(right_i)>extThre)
      right_i = right_i+1;
    end
    % add to list if duration met criteria
    dur = dt*(right_i-left_i);
    if (dur>=param.minDuration) && (dur<=param.maxDuration)
      onsets = [onsets left_i];
      offsets = [offsets right_i];
    end
  end
  
  % second pass: merge peaks that are too close
  minInterval = param.minInterval / dt;
  [mergedOnsets, mergedOffsets] = merge_adjacent_peaks(onsets, offsets, minInterval);
  
  % third pass: filter out peaks that are too short or long
  dur = (mergedOffsets - mergedOnsets) * dt; 
  passDur = (dur>=param.minDuration) & (dur<=param.maxDuration);
  onsets = mergedOnsets(passDur);
  offsets = mergedOffsets(passDur);

end


% 
% figure; plot(i1:i2, -flatness(i1:i2))
% hold on; plot(idxs, param.thresholdFlatness, 'ro');
% hold on; plot(onsets, param.extendFlatness, 'yo');
% hold on; plot(offsets, param.extendFlatness, 'yo');

% figure; 
% ax1=subplot(3,1,1);
% plot(prob_amp); hold on;
% scatter(idxs, prob_amp(idxs), 'ro');
% scatter(onsets, prob_amp(onsets), 'go');
% scatter(offsets, prob_amp(offsets), 'bo');
% ax2=subplot(3,1,2);
% plot(deriv); hold on; 
% yline(0);
% % ax3=subplot(3,1,3);
% plot(deriv2); hold on; 
% linkaxes([ax1 ax2 ax3], 'x');





