function [warbleMeta] = ZZ_IdentifyWarbleFromIntan_v1(expID, folderData, channels, saveWavPath)
% MATLAB function to extract warble from a list of wav files
% Zhilei Zhao, 09/21/2023
% Modified from the ZZ_IdentifyWarbleFromWav_v8_forPython_Ephys.m script
% Find warble from Intan recording, rather than from PyVAQ data
% Return a dataframe that specifies the start and stop location of each
% warble episode. Save as wav seperately if requested


% clear; close all;
% add folder to path
% addpath(genpath("/mnt/z4/zz367/LabSoftware/ZhileiMatlabAudioScripts"));
addpath(genpath("/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZBudgieIntanExtract"));
% also add the lab electro_gui path
addpath(genpath('/mnt/z4/zz367/LabSoftware/ElectroGui-main'));


%% scan through the Intan folder
fclose('all');
disp(folderData);

S = dir(fullfile(folderData,'*.rhd'));
S = S(~[S.isdir]);
% sort by name
[~,idx] = sort({S.name});
% sort audio files by date
% [~,idx] = sort([S.datenum]);
fileStructAll = S(idx);

% where to save the extracted audios
if ~isempty(saveWavPath)
  if exist(saveWavPath, 'dir')
    rmdir(saveWavPath, 's');
  end
  mkdir(saveWavPath);
end

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

%% Loop through each batch, calculate flatness separately for two chanels
flatnessAll = [];
% parameters related to flatness calculation
param.maskFrequency = 1000;  %ignore lower freq when calculate amplitude
param.ampIgnore = -7; %ignore where amplitude is very small
parfor batchIdx=1:length(fileStructSplit)
  fileStruct = fileStructSplit{batchIdx};
  % read all files in at once to reduce the possibility of boundaries
  signal = [];
  timeFile = [];  % duration for each audio file
  numPoints = [];  % number of data points in each intan audio file
  for i = 1: length(fileStruct)
%     [signalThis, fs] = audioread(fullfile(fileStruct(i).folder,fileStruct(i).name), 'double');
    d_intan = read_Intan_RHD2000_file_to_struct_2_noNotch(fileStruct(i).folder, fileStruct(i).name, 0);
    % focus on selected channels
    signalThis = d_intan.board_adc_data(channels,:);
    fs = d_intan.frequency_parameters.board_adc_sample_rate;
    signalThis = signalThis';
    signal = cat(1, signal, signalThis);
    timeFile = cat(1, timeFile, length(signalThis)/fs);
    numPoints = cat(1, numPoints, length(signalThis));
    %temp save wav for check
%     fn_wav = fullfile(fd_save_warble, strrep(fileStruct(i).name, 'rhd', 'wav'));
%     audiowrite(fn_wav, signalThis/10, fs);   
  end

  % go through the sound signal, calculate flatness
  [flatness1, dt] = ZZ_CalculateFlatness(signal(:,1), fs, param.ampIgnore, param.maskFrequency);
  [flatness2, dt] = ZZ_CalculateFlatness(signal(:,2), fs, param.ampIgnore, param.maskFrequency);
  flatnessAll(batchIdx).batchIdx = batchIdx;
  flatnessAll(batchIdx).fileStruct = fileStruct;
  flatnessAll(batchIdx).flatness1 = flatness1;
  flatnessAll(batchIdx).flatness2 = flatness2;
  flatnessAll(batchIdx).timeFile = timeFile;
  flatnessAll(batchIdx).numFrames = length(flatness1);
  flatnessAll(batchIdx).numPoints = numPoints; 
  flatnessAll(batchIdx).fs = fs;
  flatnessAll(batchIdx).dt = dt;
end

% remove electro_gui path to avoid conflict on the findpeaks function
rmpath(genpath('/mnt/z4/zz367/LabSoftware/ElectroGui-main'));

%% identify syllable onset/offset based on flatness
% parameters that define what count as syllables
param.minDuration = 0.04;  %unit is sec, squawks in warble is quite short
param.maxDuration = 10;   %in case there is long element
param.minInterval = 0;  % minimal interval between two syllables
% param.thresholdFlatness = -0.8;  % loose criteria
% param.extendFlatness = -0.85; % after thresholding extend to lower threshold
% for Intan, the baseline of flatness change to 0.7, so use different
% threshold
param.thresholdFlatness = -0.4;
param.extendFlatness = -0.45;

param.gapSize = 5;
fs = flatnessAll(1).fs;
dt = flatnessAll(1).dt;
flatnessFlat1 = vertcat(flatnessAll(:).flatness1);
flatnessFlat2 = vertcat(flatnessAll(:).flatness2);
% take the min of the channels
flatnessFlat = min([flatnessFlat1 flatnessFlat2], [], 2);
timeFileAll = vertcat(flatnessAll(:).timeFile);
numPointsAll = vertcat(flatnessAll(:).numPoints);
[onsets, offsets] = ZZ_GetFlatnessOnsetOffset(flatnessFlat, dt, param);

%% window-scan all onsets and offsets to check if meet warble requirements
warble.maxGap = 10;  % gap between syllables can't be larger, unit is sec
warble.minDuration = 5; % total duration can't be smaller, unit is sec
warble.minProportion = 0.1; % syllable:total duration ratio can't be smaller
% go through the syllable onsets lists, identify warble episodes
warbleStart = [];
warbleEnd = [];
if length(onsets)>1
  idxStart = 1;
  idxPrevious = 1;
  idxEnd = 2;
  while idxEnd<=length(onsets)
    % check gap size
    gap = dt * (onsets(idxEnd)-offsets(idxPrevious));
    % gap small, could be a warble
    if gap <= warble.maxGap && idxEnd<length(onsets)
      idxPrevious = idxPrevious+1;
      idxEnd = idxEnd + 1;
      % gap large, terminate, check duration and proportion
    else
      dur = dt * (offsets(idxPrevious) - onsets(idxStart));
      syllableDur = 0;
      for jj=idxStart:(idxPrevious)
        syllableDur = syllableDur + dt*(offsets(jj)-onsets(jj));
      end
      if dur>=warble.minDuration && (syllableDur/dur>=warble.minProportion)
        warbleStart = [warbleStart idxStart];
        warbleEnd = [warbleEnd idxPrevious];
      end
      % update the index number
      idxStart = idxEnd;
      idxPrevious = idxEnd;
      idxEnd = idxEnd + 1;
    end
  end
end

%% loop through identified warble episodes, save files
if ~isempty(warbleStart)
  %% locate the corresponding wav files
  % pad a few seconds before and after
  padTime = 0;
  % find file name based on real time range
  timeFileCumSum = cumsum(timeFileAll);
  numFramesCumSum = cumsum(vertcat(flatnessAll(:).numFrames));
  numPointsCumSum = cumsum(numPointsAll);
  % save information for warble episode in a struct array
  warbleMeta = [];
  for jj=1:length(warbleStart)
    warbleMeta(jj).fs = fs;
    warbleMeta(jj).dt = dt;
  end
  % compare amplitude at different frequency range to filter our noise
  param.lowFreqRange = [500, 1000];
  param.otherFreqRange = [1000, 8000];
  param.ratioThreshod = 0.01;
  param.passRatio = 0.01; %low to be sensitive, high to be selective
  % calculate what file corresponds to what frame idx
  fileFrameIdx = [];
  prev = 0;
  % need to consider the boundary effect 
  % since we split all files into batches when calculating flatness
%   for jj=1:length(flatnessAll)
%     eachFileFrame = flatnessAll(jj).numFrames / length(flatnessAll(jj).fileStruct);
%     thisRange = (prev+eachFileFrame):eachFileFrame:numFramesCumSum(jj);
%     fileFrameIdx = [fileFrameIdx thisRange];
%     prev = numFramesCumSum(jj);
%   end
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
  
  % read in the warble episodes
  parfor jj=1:length(warbleStart)
    % identify what files the onset and offset belong to
    frameStart = onsets(warbleStart(jj));
    frameEnd = offsets(warbleEnd(jj));
    for fii=1:length(fileFrameIdx)
      if fii==1
        prev = 1;
      else
        prev = fileFrameIdx(fii-1);
      end
      if (frameStart>=prev) && (frameStart<=fileFrameIdx(fii))
        fileIdxStart = fii;
      end
      if (frameEnd>prev) && (frameEnd<=fileFrameIdx(fii))
        fileIdxEnd = fii;
        break;
      end
    end
    % read intan audio in
    signalAll = [];
    fileInRange = [];
    for ii=fileIdxStart:fileIdxEnd
%       [signalThis, fs]=audioread(fullfile(fileStructAll(ii).folder,fileStructAll(ii).name));
      d_intan = read_Intan_RHD2000_file_to_struct_2_noNotch(fileStructAll(ii).folder, fileStructAll(ii).name, 0);
      % focus on selected channels
      signalThis = d_intan.board_adc_data(channels,:);
      fs = d_intan.frequency_parameters.board_adc_sample_rate;
      signalThis = signalThis';
      signalAll = [signalAll; signalThis];
      fileInRange = [fileInRange; fileStructAll(ii)];
    end
    
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
    segEnd = min([size(signalAll,1) floor(timeEndRel*fs)]);
    signalWarble = signalAll(segStart:segEnd, 1:length(channels));
    % absolute index in all signals
    segStartAbs = prevAbsIdx + segStart;
    segEndAbs = prevAbsIdx + segEnd; 
    % check ratio between two frequency range, filter out cage noise
    signalWarbleToCheck = signalWarble(:,1);
    [diffAmp] = ZZ_CalculateFreqDiff(signalWarbleToCheck, fs, param.ampIgnore, param.lowFreqRange, param.otherFreqRange);
    % num of points that pass the ratio Threshold
    passThreRatio = sum((diffAmp>=param.ratioThreshod))/length(diffAmp);
    passThre = (passThreRatio>=param.passRatio);
    % save the meta information
    warbleMeta(jj).count = jj;
    warbleMeta(jj).fs = fs;
    warbleMeta(jj).filenames = fileInRange;
    warbleMeta(jj).timeStartRel = timeStartRel;
    warbleMeta(jj).timeEndRel = timeEndRel;
    warbleMeta(jj).duration = timeEndRel-timeStartRel;
    warbleMeta(jj).absTime = fileInRange(1).date;
    warbleMeta(jj).segStartRel = segStart;
    warbleMeta(jj).segEndRel = segEnd;
    warbleMeta(jj).segStartAbs = segStartAbs;
    warbleMeta(jj).segEndAbs = segEndAbs;
    warbleMeta(jj).passThreRatio = passThreRatio;
    warbleMeta(jj).passThre = passThre;
    % save audio 
    if ~isempty(saveWavPath)
      dstr = datestr(datetime(warbleMeta(jj).absTime), 'yyyymmddHHMMSS');
      fn_save = fullfile(saveWavPath, sprintf('%s_warble_%05d_%s.wav', expID, jj, dstr));
      audiowrite(fn_save, signalWarble/10, fs); 
    end
  end
  
  % save meta information
  if ~isempty(warbleMeta)
    [a1, analysisDate, a2] = fileparts(folderData);
    fn_mat = fullfile(saveWavPath, sprintf('%s_%s_warble_meta.mat', expID, analysisDate));
    save(fn_mat, 'warbleMeta');
    % total duration of warble
    disp(strcat('Total warble duration (sec): ', num2str(sum([warbleMeta.duration]))));
  end
end
end