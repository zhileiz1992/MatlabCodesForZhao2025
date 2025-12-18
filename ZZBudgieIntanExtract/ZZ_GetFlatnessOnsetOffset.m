function [onsets, offsets] = ZZ_GetFlatnessOnsetOffset(flatness, dt, param)
% Identify the onsets and offset of sound segments that meet the
% specified requirements
% using a threshold on tonal flatness and post-hoc merging 

% param.minDuration = 0.02;  %unit is sec, squawks in warble is quite short
% param.maxDuration = 10; 
% param.minInterval = 0;  % minimal interval between two syllables
% param.thresholdFlatness = -0.6; 
% param.extendFlatness = -0.8; % after thresholding extend to lower threshold 
% param.maskFrequency = 1000;  %ignore lower freq when calculate amplitude
% param.ampIgnore = -7; %ignore where amplitude is very small
% param.gapSize = 5;  

% get peaks in flatness
[~,idxs,~,~] = findpeaks(-flatness,'MinPeakDistance',param.gapSize,'MinPeakHeight',param.thresholdFlatness);
% find the onsets and offset
onsets = [];
offsets = [];
if ~isempty(idxs)
    [onsets,offsets]=get_onset_offset_ZZ(-flatness,idxs,param.extendFlatness,param.gapSize);
    if ~isempty(onsets) 
        ind2 = [];
        % check if identified peaks meet requirements
        if length(onsets)>2
            % check distance between two syllabels larger than minInterval
            % useful to identify calls, but not warbles
            if  dt*(onsets(2)-offsets(1))>param.minInterval
                ind2 = [ind2 1];
            end
            for j=2:length(onsets)-1
                if dt*(onsets(j)-offsets(j-1))>param.minInterval && dt*(onsets(j+1)-offsets(j))>param.minInterval
                    ind2=[ind2 ,j];
                end
            end
            if dt*(onsets(end)-offsets(end-1))>param.minInterval
                ind2 = [ind2 length(onsets)]; 
            end        
        elseif length(onsets) == 1
            ind2 = 1;
        else
            if dt*(onsets(2)-offsets(1))>param.minInterval
                ind2 = [1 2];
            end
        end
    onsets = onsets(ind2);
    offsets = offsets(ind2);
    % check if the duration of syllables within range
    dur = (offsets-onsets)*dt>param.minDuration & (offsets-onsets)*dt<param.maxDuration;
    onsets= onsets(dur); offsets= offsets(dur);
    end
end
end


% 
% figure; plot(i1:i2, -flatness(i1:i2))
% hold on; plot(idxs, param.thresholdFlatness, 'ro');
% hold on; plot(onsets, param.extendFlatness, 'yo');
% hold on; plot(offsets, param.extendFlatness, 'yo');

