function [onsets, offsets, labels] = ZZ_seg_flatness_NC_v2(fns, param)
% given a list of audio files (NC format)
% segment syllables out based on spectral flatness
% loop through each file, calculate flatness, then find onsets and offsets

% loop through files, identify onsets and offsets
onsets = {};
offsets = {};
labels = {};

% read in audio in nc format
parfor fi=1:length(fns)
  fn = fns{fi};
  signal = ncread(fn, 'data');
  signal = double(signal/10);  % divide by range of Intan AI, convert to double format
  dt = ncread(fn, 'dt');
  fs = floor(1/dt);
  onsets{fi} = [];
  offsets{fi} = [];
  labels{fi} = {};
  if size(signal,1)>5000  % ignore very short files
    % calculate flatness
    [flatness, dt] = ZZ_CalculateFlatness(signal, fs, param.ampIgnore, param.maskFrequency);
    % calculate onsets and offsets
    [onsets_this, offsets_this] = ZZ_GetFlatnessOnsetOffset(flatness, dt, param);
    % change the unit from frames to data points
    onsets{fi} = floor(onsets_this*dt*fs);
    offsets{fi} = floor(offsets_this*dt*fs);
    if ~isempty(onsets_this)
      labels_temp = {};
      for si=1:length(onsets_this)
        labels_temp{si} = 'a';
      end
      labels{fi} = labels_temp;
    end
  end
end


end

