function [overlapRatio, intersectionLength] = calculateOverlapRatio(peak1, peak2)
% calculateOverlapRatio Calculates the overlap ratio of two peaks.
% The ratio is defined as the length of the intersection divided by the
% length of the union (Jaccard Index).
%
% @param peak1  A 1x2 array [onset, offset] for the first peak.
% @param peak2  A 1x2 array [onset, offset] for the second peak.
% @return overlapRatio The overlap ratio, a value between 0 and 1.
%
% Example 1 (Partial Overlap):
%   p1 = [10, 20]; % A peak from 10 to 20 (length 10)
%   p2 = [15, 25]; % A peak from 15 to 25 (length 10)
%   % Intersection is [15, 20], length 5
%   % Union is [10, 25], length 15
%   % Ratio should be 5 / 15 = 0.3333
%   ratio = calculateOverlapRatio(p1, p2)
%
% Example 2 (No Overlap):
%   p3 = [1, 5];
%   p4 = [6, 10];
%   % Intersection is 0
%   ratio_no_overlap = calculateOverlapRatio(p3, p4) % Expected: 0
%
% Example 3 (Complete Overlap / Identical):
%   p5 = [50, 100];
%   p6 = [50, 100];
%   ratio_identical = calculateOverlapRatio(p5, p6) % Expected: 1

% --- Input Validation ---
if ~isnumeric(peak1) || ~isnumeric(peak2) || numel(peak1) ~= 2 || numel(peak2) ~= 2
    error('Inputs must be 1x2 numeric arrays of the form [onset, offset].');
end
if peak1(1) > peak1(2) || peak2(1) > peak2(2)
    error('Onset cannot be greater than offset for a peak.');
end

% --- Calculate Intersection (The overlapping portion) ---
% The start of the overlap is the latest of the two start times.
intersectionStart = max(peak1(1), peak2(1));

% The end of the overlap is the earliest of the two end times.
intersectionEnd = min(peak1(2), peak2(2));

% Calculate the length of the intersection.
% If the intersection end is before the start, there is no overlap, so the
% length is 0.
intersectionLength = max(0, intersectionEnd - intersectionStart);

% --- Calculate Union (The total range covered by both peaks) ---
% The length of the union can be calculated as:
% Length(Peak1) + Length(Peak2) - Length(Intersection)
length1 = peak1(2) - peak1(1);
length2 = peak2(2) - peak2(1);
unionLength = length1 + length2 - intersectionLength;


% --- Calculate Overlap Ratio ---
% This is the final ratio of Intersection over Union.
% We must handle the case where the union length is zero to avoid
% dividing by zero. This occurs if both peaks are identical points in time
% (e.g., [5, 5] and [5, 5]).
if unionLength == 0
    if intersectionLength > 0
        overlapRatio = 1; % Peaks are identical zero-length points
    else
        overlapRatio = 0; % Non-identical zero-length points
    end
else
    overlapRatio = intersectionLength / unionLength;
end

end

