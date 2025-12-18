function [onsets, offsets] = ZZ_seg_flatness_NC_v1(fns)
% given a list of audio files (NC format)
% segment syllables out based on spectral flatness
% loop through each file, calculate flatness, then find onsets and offsets

% parameters
param.maskFrequency = 1000;  %ignore lower freq when calculate amplitude
param.ampIgnore = -7; %ignore where amplitude is very small
% threshold to identify peaks
param.thresholdFlatness = -0.6;
param.extendFlatness = -0.7;
param.gapSize = 5;  % merge peaks if gap is small
param.minDuration = 0.04;  %unit is sec, squawks in warble is quite short
param.maxDuration = 10;   %in case there is long element
param.minInterval = 0;  % minimal interval between two syllables

% loop through files, identify onsets and offsets
onsets = {};
offsets = {};

% read in audio in nc format
parfor fi=1:length(fns)
  fn = fns{fi};
  signal = ncread(fn, 'data');
  signal = double(signal/10);  % divide by range of Intan AI, convert to double format
  dt = ncread(fn, 'dt');
  fs = floor(1/dt);
  onsets{fi} = [];
  offsets{fi} = [];
  if size(signal,1)>256
    % calculate flatness
    [flatness, dt] = ZZ_CalculateFlatness(signal, fs, param.ampIgnore, param.maskFrequency);
    % calculate onsets and offsets
    [onsets_this, offsets_this] = ZZ_GetFlatnessOnsetOffset(flatness, dt, param);
    % change the unit from frames to data points 
    onsets{fi} = floor(onsets_this*dt*fs); 
    offsets{fi} = floor(offsets_this*dt*fs);
  end
end


end

