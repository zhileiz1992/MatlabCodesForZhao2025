function [mergedOnsets, mergedOffsets] = merge_adjacent_peaks(onsets, offsets, minInterval)
%merge_adjacent_peaks Merges adjacent peaks based on a minimum interval threshold.
%
%   [mergedOnsets, mergedOffsets] = merge_adjacent_peaks(onsets, offsets, minInterval)
%   iteratively merges adjacent or overlapping peaks. The merging process
%   continues until no more adjacent peaks meet the merge criteria.
%
%   Inputs:
%       onsets      - A 1D array of peak start indices.
%       offsets     - A 1D array of peak end indices.
%       minInterval - A scalar value. If the distance between the offset of
%                     one peak and the onset of the next is less than this
%                     value, the peaks will be merged. This also handles
%                     overlapping peaks where the distance is negative.
%
%   Outputs:
%       mergedOnsets  - The new 1D array of onsets after merging.
%       mergedOffsets - The new 1D array of offsets after merging.
%
%   Example:
%       onsets = [10, 35, 50, 80];
%       offsets = [20, 45, 65, 90];
%       minInterval = 10;
%       [newOnsets, newOffsets] = merge_adjacent_peaks(onsets, offsets, minInterval);
%       % newOnsets will be [10, 35, 80]
%       % newOffsets will be [20, 65, 90] (peak 2 and 3 merged)

% --- Input Validation ---
if isempty(onsets)
    mergedOnsets = [];
    mergedOffsets = [];
    return;
end

if numel(onsets) ~= numel(offsets)
    error('Onsets and offsets arrays must have the same number of elements.');
end

% Ensure inputs are column vectors for easier indexing
mergedOnsets = onsets(:);
mergedOffsets = offsets(:);

if numel(mergedOnsets) < 2
    return; % Nothing to merge if there's only one or zero peaks
end

% --- Iterative Merging Logic ---
wasMergedInPass = true; % Flag to control the iterative process
while wasMergedInPass
    wasMergedInPass = false; % Assume no merges will happen in this pass
    
    i = 1;
    while i < numel(mergedOnsets)
        % Calculate distance to the next peak. A negative distance indicates an overlap.
        distance = mergedOnsets(i+1) - mergedOffsets(i);
        
        % Check if the peaks are close enough to merge or if they overlap.
        if distance < minInterval
            
            % A merge is necessary. Set the flag to true.
            wasMergedInPass = true;
            
            % The new merged peak will span from the start of the first peak
            % to the end of the second peak. We use max() to handle cases
            % where the first peak completely contains the second one.
            mergedOffsets(i) = max(mergedOffsets(i), mergedOffsets(i+1));
            
            % Remove the second peak, as it's now merged into the first.
            mergedOnsets(i+1) = [];
            mergedOffsets(i+1) = [];
            
            % After a merge, we break the inner loop and restart the whole
            % process from the beginning. This ensures that new adjacent
            % peaks (like a chain of 3 or more) are correctly merged.
            break;
        else
            % No merge, move to the next peak
            i = i + 1;
        end
    end
end

end

