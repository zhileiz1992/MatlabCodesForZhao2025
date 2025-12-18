function [warbleMetaAll, otherMetaAll] = ZZ_IdentifyWarbleFromIntan_v7(fd_base_intan, fd_base_save, data_date, meta)
% MATLAB function to extract warble from a list of wav files
% Zhilei Zhao, 02/02/2024
% Modified from the ZZ_IdentifyWarbleFromWav_v8_forPython_Ephys.m script
% Find warble from Intan recording, rather than from PyVAQ data
% Return two struct arrays that contain information about the extracted
% warble episodes and other sounds that cross the threshold
% Save as wav seperately if saveWavPath is specificed
% 04/02/2024 modify, when calculating flatness, also save a copy of wav files, so later no need to read all large intan files
% again; also screen 'other' category to save calls only


% clear; close all;
% add folder to path
% addpath(genpath("/mnt/z4/zz367/LabSoftware/ZhileiMatlabAudioScripts"));
addpath(genpath("/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZBudgieIntanExtract"));
% also add the lab electro_gui path
addpath(genpath('/mnt/z4/zz367/LabSoftware/ElectroGui-main'));
addpath(genpath('/home/zz367/ProjectsU/WarbleAnalysis/Jupyter/MatlabCodes/ZZ_extractWarbleFromWav'));


% if ~isempty(saveWavPath)
%   % where to save the extracted audios
%   if exist(saveWavPath, 'dir')
%     rmdir(saveWavPath, 's');
%   end
%   saveWavPathWarble = fullfile(saveWavPath, 'warble');
%   saveWavPathOther = fullfile(saveWavPath, 'other');
%   mkdir(saveWavPathWarble);
%   mkdir(saveWavPathOther);
% end

%% grab files in the intan folder
fclose('all');
folderData = fullfile(fd_base_intan, data_date);
disp(folderData);

S = dir(fullfile(folderData,'*.rhd'));
S = S(~[S.isdir]);
% sort by name, since Intan file names are already sorted
[~,idx] = sort({S.name});
files = S(idx);
% sort audio files by date
% [~,idx] = sort([S.datenum]);
% ignore files during light off period, pad 15 mins
startTime = datenum(meta.light_cycle{1}, 'HH:MM:SS');
endTime = datenum(meta.light_cycle{2}, 'HH:MM:SS');
keep_idx = [];
for i = 1:length(files)
  filename = files(i).name;
  timeStr = filename(end-9:end-4); % Extracts '103938' from 'XXX_240331_103938.rhd'
  fileTime = datenum(timeStr, 'HHMMSS');
  % Check if the file's time is within the desired range
  if fileTime >= startTime && fileTime <= endTime
    keep_idx = [keep_idx i];
  end
end
fileStructAll = files(keep_idx);


%% Calculate flatness of all audio
% num audio files to process per batch, decrease this if have small RAM
% for intan each file has 1-min recording
batchSize = 10;

% split the files into batches
fileStructSplit = {};
for idx=1:batchSize:length(fileStructAll)
  if (idx+batchSize) > length(fileStructAll)
    endIdx = length(fileStructAll);
  else
    endIdx = idx + batchSize - 1;
  end
  thisBatch = fileStructAll(idx:endIdx);
  fileStructSplit{floor(idx/batchSize)+1} = thisBatch;
end

%% Loop through each batch, calculate flatness separately for two chanels of each pair
flatnessAll = [];
% parameters related to flatness calculation
param.maskFrequency = 1000;  %ignore lower freq when calculate amplitude
param.ampIgnore = -7; %ignore where amplitude is very small
% check all channels
channels = reshape(meta.audio_ch_intan.', 1, []);
parfor batchIdx=1:length(fileStructSplit)
  fileStruct = fileStructSplit{batchIdx};
  % read all files in at once to reduce the possibility of boundaries
  signal = [];
  timeFile = [];  % duration for each audio file
  numPoints = [];  % number of data points in each intan audio file
  for i = 1:length(fileStruct)
    %     [signalThis, fs] = audioread(fullfile(fileStruct(i).folder,fileStruct(i).name), 'double');
    d_intan = read_Intan_RHD2000_file_to_struct_2_noNotch(fileStruct(i).folder, fileStruct(i).name, 0);
    % save all channels including speaking copy into a separate wav file
    fs = d_intan.frequency_parameters.board_adc_sample_rate;
    signalAll = d_intan.board_adc_data;
    signalAll = signalAll' / 10; % devide by the range of AI
    for pi=1:length(meta.pairID)
      fd_save_wav_this = fullfile(fd_base_save, meta.pairID{pi}, data_date, 'AllWav');
      if ~exist(fd_save_wav_this, 'dir')
        mkdir(fd_save_wav_this);
      end
      fn_this = fullfile(fd_save_wav_this, strrep(fileStruct(i).name, 'rhd', 'wav'));
      ch_this = meta.audio_ch_intan(pi,:);
      audiowrite(fn_this, signalAll(:, ch_this), fs);
    end
    % also save the speaker playback and pyvaq trigger in the folder of 1st pair
    fd_save_speaker = fullfile(fd_base_save, meta.pairID{1}, data_date, 'SpeakerTrigger');
    if ~exist(fd_save_speaker, 'dir')
      mkdir(fd_save_speaker);
    end
    d_trigger = d_intan.board_dig_in_data;
    d_comb = [signalAll(:, meta.playback_ch_intan) d_trigger'];
    fn_this = fullfile(fd_save_speaker, strrep(fileStruct(i).name, 'rhd', 'wav'));
    audiowrite(fn_this, d_comb, fs);
    % then concatenate data points for flatness calculation focus on selected channels
    signalThis = signalAll(:, channels);
    signal = cat(1, signal, signalThis);
    timeFile = cat(1, timeFile, length(signalThis)/fs);
    numPoints = cat(1, numPoints, length(signalThis));
  end
  
  % go through the sound signal, calculate flatness
  flatness = [];
  for jj=1:size(signal,2)
    [flatness_this, dt] = ZZ_CalculateFlatness(signal(:,jj), fs, param.ampIgnore, param.maskFrequency);
    flatness = [flatness flatness_this];
  end
  %   [flatness2, dt] = ZZ_CalculateFlatness(signal(:,2), fs, param.ampIgnore, param.maskFrequency);
  flatnessAll(batchIdx).batchIdx = batchIdx;
  flatnessAll(batchIdx).fileStruct = fileStruct;
  flatnessAll(batchIdx).flatness = flatness;
  flatnessAll(batchIdx).timeFile = timeFile;
  flatnessAll(batchIdx).numFrames = size(flatness, 1);
  flatnessAll(batchIdx).numPoints = numPoints;
  flatnessAll(batchIdx).fs = fs;
  flatnessAll(batchIdx).dt = dt;
end
% save the flatness struct in case needed later
fn_mat = fullfile(fd_base_save, meta.pairID{1}, data_date, 'SpeakerTrigger', sprintf('%s_%s.flatness.mat', meta.expID, data_date));
save(fn_mat, 'flatnessAll');

% remove electro_gui path to avoid conflict on the findpeaks function
rmpath(genpath('/mnt/z4/zz367/LabSoftware/ElectroGui-main'));


%% Loop through pairs to identify warble
warbleMetaAll = {};
otherMetaAll = {};
for pi=1:length(meta.pairID)
  pairName = meta.pairID{pi};
  birdNames = meta.birdID{pi};
  if isempty(birdNames{1}) && isempty(birdNames{2})  % ignore empty pairs
    continue
  end
  % where to save the results
  saveWavPath = fullfile(fd_base_save, pairName, data_date, 'WarbleWav');
  if exist(saveWavPath)
    rmdir(saveWavPath, 's');
  end
  saveWavPathWarble = fullfile(saveWavPath, 'warble');
  mkdir(saveWavPathWarble);
  saveWavPathOther = fullfile(saveWavPath, 'other');
  mkdir(saveWavPathOther);
  %% identify syllable onset/offset based on flatness
  % parameters that define what count as syllables
  param.minDuration = 0.02;  %unit is sec, squawks in warble is quite short
  param.maxDuration = 10;   %in case there is long element
  param.minInterval = 0;  % minimal interval between two syllables
  % threshold to identify peaks: this is setup-specific and depends on the
  % quality of recording, use the 'troubleshoot.m' function to plot flatness
  % and check if the threshold makes sense
  param.thresholdFlatness = -0.75;
  param.extendFlatness = -0.85;
  param.gapSize = 5;  % two syllables with gaps smaller than this will be merged
  fs = flatnessAll(1).fs;
  dt = flatnessAll(1).dt;
  temp = vertcat(flatnessAll(:).flatness);
  flatnessFlat1 = temp(:,2*pi-1);
  flatnessFlat2 = temp(:,2*pi);
  % take the min of the channels
  flatnessFlat = min([flatnessFlat1 flatnessFlat2], [], 2);
  timeFileAll = vertcat(flatnessAll(:).timeFile);
  numPointsAll = vertcat(flatnessAll(:).numPoints);
  [onsets, offsets] = ZZ_GetFlatnessOnsetOffset(flatnessFlat, dt, param);
  
  %% window-scan all onsets and offsets to check if meet warble requirements
  warbleParam.maxGap = 10;  % gap between syllables can't be larger, unit is sec
  warbleParam.minDuration = 5; % total duration can't be smaller, unit is sec
  warbleParam.minNumSyllable = 5; % total number of syllables can't be smaller
  warbleParam.minPeakDensity = 3; % calculate the number of syllables within a 5sec sliding window, peak density can't be lower
  [warbleStart,warbleEnd] = ZZ_OnsetOffsetToWarbleEpisode_v1(onsets, offsets, dt, warbleParam);
  
  
  % pad a few seconds before and after the warble onset/offset
  padTimeBefore = 3;
  padTimeAfter = 3;
  padFramesBefore = floor(padTimeBefore/dt);
  padFramesAfter = floor(padTimeAfter/dt);
  
  %% loop through identified warble episodes, save files
  % calculate what file corresponds to what frame idx
  fileFrameIdx = [];
  prev = 0;
  % need to consider the boundary effect
  % since we split all files into batches when calculating flatness
      % find file name based on real time range
  timeFileCumSum = cumsum(timeFileAll);
  numFramesCumSum = cumsum(vertcat(flatnessAll(:).numFrames));
  numPointsCumSum = cumsum(numPointsAll);
  for bi=1:length(flatnessAll)
    % how many frame per data point
    framePerPoint = flatnessAll(bi).numFrames / sum(flatnessAll(bi).numPoints);
    for fi=1:length(flatnessAll(bi).fileStruct)
      % how many frames for this file depends on its nubmer of data points
      framesThisFile =  flatnessAll(bi).numPoints(fi) * framePerPoint;
      fileFrameIdx = [fileFrameIdx prev+framesThisFile];
      prev = prev+framesThisFile;
    end
    prev = numFramesCumSum(bi);
  end
  % need to round the last element 
  fileFrameIdx(end) = length(flatnessFlat);
  if ~isempty(warbleStart)
    %% locate the corresponding wav files
    % save information for warble episode in a struct array
    warbleMeta = [];
    for jj=1:length(warbleStart)
      warbleMeta(jj).fs = fs;
      warbleMeta(jj).dt = dt;
    end
    % determine the onsets and offsets of warble in the rhd files
    fs = flatnessAll(1).fs;
    parfor jj=1:length(warbleStart)
      % identify what files the onset and offset belong to
      % add padding, caution about boundary condition
      frameStart = max([1, onsets(warbleStart(jj))-padFramesBefore]);
      frameEnd = min([size(flatnessFlat,1), offsets(warbleEnd(jj))+padFramesAfter]);
%       fileIdxStart = 1; 
%       fileIdxEnd = 1; 
      for fii=1:length(fileFrameIdx)
        if fii==1
          prev = 1;
        else
          prev = fileFrameIdx(fii-1);
        end
        if (frameStart>=prev) && (frameStart<=fileFrameIdx(fii))
          fileIdxStart = fii;
        end
        if (frameEnd>=prev) && (frameEnd<=fileFrameIdx(fii))
          fileIdxEnd = fii;
          break;
        end
      end
      fileInRange = [fileStructAll(fileIdxStart:fileIdxEnd)];
      % calculate when the warble starts and ends in the signal
      if fileIdxStart==1
        prev = 0;
        prevAbsIdx = 0;
      else
        prev = fileFrameIdx(fileIdxStart-1);   % starting frame index
        prevAbsIdx = numPointsCumSum(fileIdxStart-1);  % starting point index
      end
      % relative time
      timeStartRel = dt*(frameStart - prev);
      timeEndRel = dt*(frameEnd - prev);
      % relative index in the signal array, deal with boundaries
      segStart = max([1 floor(timeStartRel*fs)]);
      segEnd = min([sum(numPointsAll(fileIdxStart:fileIdxEnd)) floor(timeEndRel*fs)]);
      % absolute index in all signals
      segStartAbs = prevAbsIdx + segStart;
      segEndAbs = prevAbsIdx + segEnd;
      % save the meta information
      warbleMeta(jj).count = jj;
      warbleMeta(jj).fs = fs;
      warbleMeta(jj).filenames = fileInRange;
      warbleMeta(jj).timeStartRel = timeStartRel;
      warbleMeta(jj).timeEndRel = timeEndRel;
      warbleMeta(jj).duration = timeEndRel-timeStartRel;
      warbleMeta(jj).segStartRel = segStart;
      warbleMeta(jj).segEndRel = segEnd;
      warbleMeta(jj).segStartAbs = segStartAbs;
      warbleMeta(jj).segEndAbs = segEndAbs;
      warbleMeta(jj).syllableOnsetsTime = (onsets(warbleStart(jj):warbleEnd(jj))-frameStart)*dt;
      warbleMeta(jj).syllableOffsetsTime = (offsets(warbleStart(jj):warbleEnd(jj))-frameStart)*dt;
      % when does the warble start, this is just an estimate
      temp = strsplit(fileInRange(1).name, '_');
      fileTime = sprintf('%s_%s', temp{end-1}, strrep(temp{end}, '.rhd', ''));
      warbleMeta(jj).absTimeEst = datetime(fileTime,"InputFormat",'yyMMdd_HHmmss') + seconds(timeStartRel);
      % save audio as wav if requested
      if ~isempty(saveWavPath)
        signalAll = [];
        for ii=fileIdxStart:fileIdxEnd
          %         d_intan = read_Intan_RHD2000_file_to_struct_2_noNotch(fileStructAll(ii).folder, fileStructAll(ii).name, 0);
          fn_wav = strrep(fileStructAll(ii).name, 'rhd', 'wav');
          [signalThis, fs2] = audioread(fullfile(fd_base_save, pairName, data_date, 'AllWav', fn_wav));
          % focus on selected channels
          %         signalThis = d_intan.board_adc_data(channels,:);
          %         fs = d_intan.frequency_parameters.board_adc_sample_rate;
          %         signalThis = signalThis';
          signalAll = [signalAll; signalThis];
        end
        signalWarble = signalAll(segStart:segEnd, :);
        dstr = datestr(datetime(warbleMeta(jj).absTimeEst), 'yyyymmddHHMMSS');
        fn_save = fullfile(saveWavPathWarble, sprintf('%s_warble_%05d_%s.wav', pairName, jj, dstr));
        %       audiowrite(fn_save, signalWarble/10, fs);
        audiowrite(fn_save, signalWarble, fs2);
      end
    end
  end
  
  %% screen remaining and save calls in a different struct array
  % identify calls by freq ratio
  criteria.minDur = 0.1;
  criteria.maxDur = 0.35;
  criteria.minFreq = 2000;
  criteria.maxFreq = 4000;
  criteria.minFreqRatio = 2;
  otherMeta = [];
  % what onset/offset index already included in warble
  warbleFrameIdx = [];
  for ii=1:length(warbleStart)
    warbleFrameIdx = [warbleFrameIdx warbleStart(ii):warbleEnd(ii)];
  end
  otherIdx = setdiff(1:length(onsets), warbleFrameIdx);
  % init the struct
  for j_idx=1:length(otherIdx)
    %   otherMeta(j_idx) = struct();
    otherMeta(j_idx).count = j_idx;
  end
  parfor j_idx=1:length(otherIdx)
%   for j_idx=1:10
    jj = otherIdx(j_idx);
    % identify what files the onset and offset belong to
    %     frameStart = onsets(jj);
    %     frameEnd = offsets(jj);
    % consider padding
    frameStart = max([1, onsets(jj)-padFramesBefore]);
    frameEnd = min([size(flatnessFlat,1)-1, offsets(jj)+padFramesAfter]);
    for fii=1:length(fileFrameIdx)
      if fii==1
        prev = 1;
      else
        prev = fileFrameIdx(fii-1);
      end
      if (frameStart>=prev) && (frameStart<=fileFrameIdx(fii))
        fileIdxStart2 = fii;
      end
      if (frameEnd>=prev) && (frameEnd<=fileFrameIdx(fii))
        fileIdxEnd2 = fii;
        break;
      end
    end
    fileInRange = [fileStructAll(fileIdxStart2:fileIdxEnd2)];
    % calculate when the warble starts and ends in the signal
    if fileIdxStart2==1
      prev = 0;
      prevAbsIdx = 0;
    else
      prev = fileFrameIdx(fileIdxStart2-1);   % starting frame index
      prevAbsIdx = numPointsCumSum(fileIdxStart2-1);  % starting point index
    end
    % relative time
    timeStartRel = dt*(frameStart - prev);
    timeEndRel = dt*(frameEnd - prev);
    % relative index in the signal array, deal with boundaries
    segStart = max([1 floor(timeStartRel*fs)]);
    segEnd = min([sum(numPointsAll(fileIdxStart2:fileIdxEnd2)) floor(timeEndRel*fs)]);
    % absolute index in all signals
    segStartAbs = prevAbsIdx + segStart;
    segEndAbs = prevAbsIdx + segEnd;
    % check if it meets the criteria for calls
    signalAll = [];
    for ii=fileIdxStart2:fileIdxEnd2
      fn_wav = strrep(fileStructAll(ii).name, 'rhd', 'wav');
      [signalThis, fs2] = audioread(fullfile(fd_base_save, pairName, data_date, 'AllWav', fn_wav));
      signalAll = [signalAll; signalThis];
    end
    % exclude the padding part
    realStart = floor(dt*(onsets(jj)-prev)*fs);
    realEnd = floor(dt*(offsets(jj)-prev)*fs);
    signalSeg = signalAll(realStart:realEnd, :);
%     disp(fn_wav);
    amp_sum = sum(abs(signalSeg));
    [mV, mi] = max(amp_sum);
    signalSegMax = signalSeg(:,mi);
    [power, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(signalSegMax, fs);
    % calculate the frequency ratio
    fidx_in = find((f>=criteria.minFreq) & (f<=criteria.maxFreq));
    fidx_out = find((f<criteria.minFreq) | (f>criteria.maxFreq));
    spec = powerGrey;
    mean_freq_in = mean(mean(spec(fidx_in, :)));
    mean_freq_out = mean(mean(spec(fidx_out, :)));
    freqRatio = mean_freq_in / mean_freq_out;
    dur = size(signalSeg,1) / fs;
    % save if meet criteria
    if (dur>=criteria.minDur) && (dur<=criteria.maxDur) && (freqRatio>=criteria.minFreqRatio) 
      % save the meta information
      %     otherCount = otherCount+1;
      %     otherMeta(j_idx).count = j_idx;
      otherMeta(j_idx).fs = fs;
      otherMeta(j_idx).filenames = fileInRange;
      otherMeta(j_idx).timeStartRel = timeStartRel;
      otherMeta(j_idx).timeEndRel = timeEndRel;
      otherMeta(j_idx).duration = timeEndRel-timeStartRel;
      otherMeta(j_idx).segStartRel = segStart;
      otherMeta(j_idx).segEndRel = segEnd;
      otherMeta(j_idx).segStartAbs = segStartAbs;
      otherMeta(j_idx).segEndAbs = segEndAbs;
      otherMeta(j_idx).syllableOnsetsTime = (onsets(jj)-frameStart)*dt;
      otherMeta(j_idx).syllableOffsetsTime = (offsets(jj)-frameStart)*dt;
      % when does the warble start, this is just an estmiate
      temp = strsplit(fileInRange(1).name, '_');
      fileTime = sprintf('%s_%s', temp{end-1}, strrep(temp{end}, '.rhd', ''));
      otherMeta(j_idx).absTimeEst = datetime(fileTime,"InputFormat",'yyMMdd_HHmmss') + seconds(timeStartRel);
      % save audio as wav if requested
      channels = meta.audio_ch_intan(pi,:);
      if ~isempty(saveWavPath)
        signalOther = signalAll(segStart:segEnd, :);
        dstr = datestr(datetime(otherMeta(j_idx).absTimeEst), 'yyyymmddHHMMSS');
        fn_save = fullfile(saveWavPathOther, sprintf('%s_other_%05d_%s.wav', pairName, j_idx, dstr));
        audiowrite(fn_save, signalOther, fs2);
      end
    end
  end
  % clean the other meta
  s = otherMeta;
  otherMeta = s(~arrayfun(@(x) isempty(x.fs), s));
  
  % save meta data
  warbleMetaAll{pi} = warbleMeta; 
  otherMetaAll{pi} = otherMeta; 
  fn_meta_warble = fullfile(fd_base_save,  pairName, data_date, sprintf('%s_%s_warble_meta.mat', pairName, data_date));
  save(fn_meta_warble, 'warbleMeta');
  fn_meta_other = fullfile(fd_base_save,  pairName, data_date, sprintf('%s_%s_other_meta.mat', pairName, data_date));
  save(fn_meta_other, 'otherMeta')
end
end
