function out = detect_two_peaks(prob_amp, varargin)
%DETECT_TWO_PEAKS  Find two adjacent major peaks and their onsets/offsets.
%
% out = detect_two_peaks(prob_amp, 'X', x, 'SmoothW', w, 'MinProm', p, 'MinDist', d)
%  - R2021a-compatible: uses at most 4 outputs from findpeaks.
%
% Onsets/offsets are taken as the nearest local minima bracketing each peak.
% If a minimum is missing, falls back to the data edge or to the argmin
% between the two peaks for the shared boundary.

% ---------- parse ----------
y = prob_amp(:);
N = numel(y);
ip = inputParser;
addParameter(ip,'X',(1:N)',@(v)isvector(v)&&numel(v)==N);
defW = max(5, round(N/40)*2+1);          % odd window
addParameter(ip,'SmoothW',defW,@(v)isnumeric(v)&&isscalar(v)&&v>=3);
addParameter(ip,'MinProm',[],@(v)isnumeric(v)&&isscalar(v)&&v>=0);
addParameter(ip,'MinDist',max(3, round(N/20)),@(v)isnumeric(v)&&isscalar(v)&&v>=1);
parse(ip,varargin{:});
x        = ip.Results.X(:);
w        = ip.Results.SmoothW;
minDist  = ip.Results.MinDist;

% ---------- robust smoothing ----------
y_med = movmedian(y, w, 'Endpoints','shrink');
y_sm  = movmean(y_med, w, 'Endpoints','shrink');

% data-driven prominence threshold
rngy     = range(y_sm);
noise    = 1.4826*mad(diff(y_sm),1);     % robust noise scale
autoProm = max(3*noise, 0.02*rngy);
minProm  = ip.Results.MinProm; if isempty(minProm), minProm = autoProm; end

% ---------- peaks ----------
[pkVal, pkIdx, pkW, pkProm] = findpeaks(y_sm, ...
    'MinPeakProminence', minProm, ...
    'MinPeakDistance',   minDist, ...
    'WidthReference','halfheight'); %#ok<ASGLU>

if numel(pkIdx) < 2
    % slightly relax if needed
    [pkVal, pkIdx, pkW, pkProm] = findpeaks(y_sm, ...
        'MinPeakProminence', 0.5*minProm, ...
        'MinPeakDistance',   max(1,round(minDist*0.6)), ...
        'WidthReference','halfheight'); %#ok<ASGLU>
end
if isempty(pkIdx)
    warning('detect_two_peaks:NoPeaks','No peaks found.');
end

% keep the two most prominent peaks; sort left→right
[~, ord] = sort(pkProm, 'descend');
ord = ord(1:min(2,numel(ord)));
pkIdx = sort(pkIdx(ord));
pkVal = y(pkIdx);

% ---------- minima & boundaries ----------
isMin  = islocalmin(y_sm, 'FlatSelection','center');
minIdx = find(isMin);
% include edges as permissible minima
if isempty(minIdx) || minIdx(1) ~= 1, minIdx = [1; minIdx]; end
if minIdx(end) ~= N, minIdx = [minIdx; N]; end

valley_idx = NaN;
onIdx = NaN(2,1); offIdx = NaN(2,1);

if numel(pkIdx) >= 1
    % per-peak bracketing minima
    for k = 1:min(2,numel(pkIdx))
        left  = minIdx(minIdx < pkIdx(k));
        right = minIdx(minIdx > pkIdx(k));
        onIdx(k)  = ternary(~isempty(left),  left(end),  1);
        offIdx(k) = ternary(~isempty(right), right(1),   N);
    end
end

% shared valley between two adjacent peaks (if present)
if numel(pkIdx) >= 2
    mask = minIdx > pkIdx(1) & minIdx < pkIdx(2);
    if any(mask)
        mids = minIdx(mask);
        [~, ii] = min(y_sm(mids));
        valley_idx = mids(ii);
    else
        % fallback: strict argmin in-between
        [~, ii] = min(y_sm(pkIdx(1):pkIdx(2)));
        valley_idx = pkIdx(1) + ii - 1;
    end
    % enforce adjacency using the valley as the shared boundary
    onIdx(2)  = max(onIdx(2),  valley_idx);
    offIdx(1) = min(offIdx(1), valley_idx);
end

% ---------- package output ----------
K = min(2, numel(pkIdx));
peakTbl = table((1:K)', pkIdx(1:K)', x(pkIdx(1:K)), y(pkIdx(1:K)), ...
    'VariableNames', {'peak_no','idx','x','y'});

boundsTbl = table((1:K)', onIdx(1:K), offIdx(1:K), ...
    x(onIdx(1:K)), x(offIdx(1:K)), ...
    'VariableNames', {'peak_no','onset_idx','offset_idx','onset_x','offset_x'});

out = struct('x',x,'y',y,'y_smooth',y_sm, ...
             'peaks',peakTbl, ...
             'valley_idx',valley_idx, ...
             'bounds',boundsTbl);

end

% --- small helper ---
function z = ternary(cond,a,b)
if cond, z = a; else, z = b; end
end
