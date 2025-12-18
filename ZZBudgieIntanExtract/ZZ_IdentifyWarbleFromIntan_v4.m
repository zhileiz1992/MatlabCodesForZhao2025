function [warbleMeta, otherMeta] = ZZ_IdentifyWarbleFromIntan_v4(expID, folderData, channels, saveWavPath)
% MATLAB function to extract warble from a list of wav files
% Zhilei Zhao, 02/02/2024
% Modified from the ZZ_IdentifyWarbleFromWav_v8_forPython_Ephys.m script
% Find warble from Intan recording, rather than from PyVAQ data
% Return two struct arrays that contain information about the extracted
% warble episodes and other sounds that cross the threshold 
% Save as wav seperately if saveWavPath is specificed


% clear; close all;
% add folder to path
% addpath(genpath("/mnt/z4/zz367/LabSoftware/ZhileiMatlabAudioScripts"));
addpath(genpath("/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZBudgieIntanExtract"));
% also add the lab electro_gui path
addpath(genpath('/mnt/z4/zz367/LabSoftware/ElectroGui-main'));
addpath(genpath('/home/zz367/ProjectsU/WarbleAnalysis/Jupyter/MatlabCodes/ZZ_extractWarbleFromWav'));


if ~isempty(saveWavPath)
  % where to save the extracted audios
  if exist(saveWavPath, 'dir')
    rmdir(saveWavPath, 's');
  end
  saveWavPathWarble = fullfile(saveWavPath, 'warble');
  saveWavPathOther = fullfile(saveWavPath, 'other');
  mkdir(saveWavPathWarble);
  mkdir(saveWavPathOther);
end

%% scan through the Intan folder
fclose('all');
disp(folderData);

S = dir(fullfile(folderData,'*.rhd'));
S = S(~[S.isdir]);
% sort by name, since Intan file names are already sorted
[~,idx] = sort({S.name});
% sort audio files by date
% [~,idx] = sort([S.datenum]);
fileStructAll = S(idx);


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
  for i = 1:length(fileStruct)
    %     [signalThis, fs] = audioread(fullfile(fileStruct(i).folder,fileStruct(i).name), 'double');
    d_intan = read_Intan_RHD2000_file_to_struct_2_noNotch(fileStruct(i).folder, fileStruct(i).name, 0);
    % focus on selected channels
    signalThis = d_intan.board_adc_data(channels,:);
    fs = d_intan.frequency_parameters.board_adc_sample_rate;
%     signalThis = signalThis';
    signalThis = signalThis' / 10;  % devide by the range of AI
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
param.minDuration = 0.02;  %unit is sec, squawks in warble is quite short
param.maxDuration = 10;   %in case there is long element
param.minInterval = 0;  % minimal interval between two syllables
% threshold to identify peaks: this is setup-specific and depends on the
% quality of recording, use the 'troubleshoot.m' function to plot flatness
% and check if the threshold makes sense
param.thresholdFlatness = -0.6;
param.extendFlatness = -0.7;
param.gapSize = 5;  % two syllables with gaps smaller than this will be merged
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
if ~isempty(warbleStart)    
  %% locate the corresponding wav files
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
  % calculate what file corresponds to what frame idx
  fileFrameIdx = [];
  prev = 0;
  % need to consider the boundary effect
  % since we split all files into batches when calculating flatness
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
  
  % determine the onsets and offsets of warble in the rhd files
  parfor jj=1:length(warbleStart)
%   parfor jj=1:24
    fs = flatnessAll(1).fs;
    % identify what files the onset and offset belong to
    % add padding, caution about boundary condition 
    frameStart = max([1, onsets(warbleStart(jj))-padFramesBefore]);
    frameEnd = min([size(flatnessFlat,1), offsets(warbleEnd(jj))+padFramesAfter]);
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
        d_intan = read_Intan_RHD2000_file_to_struct_2_noNotch(fileStructAll(ii).folder, fileStructAll(ii).name, 0);
        % focus on selected channels
        signalThis = d_intan.board_adc_data(channels,:);
        fs = d_intan.frequency_parameters.board_adc_sample_rate;
        signalThis = signalThis';
        signalAll = [signalAll; signalThis];
      end
      signalWarble = signalAll(segStart:segEnd, :);
      dstr = datestr(datetime(warbleMeta(jj).absTimeEst), 'yyyymmddHHMMSS');
      fn_save = fullfile(saveWavPathWarble, sprintf('%s_warble_%05d_%s.wav', expID, jj, dstr));
      audiowrite(fn_save, signalWarble/10, fs);
    end
  end
end

%% save remaining in a different struct array
otherMeta = [];
otherCount = 0; 
% what onset/offset index already included in warble
warbleFrameIdx = [];
for ii=1:length(warbleStart)
  warbleFrameIdx = [warbleFrameIdx warbleStart(ii):warbleEnd(ii)];
end
for jj=1:length(onsets)
  if ~ismember(jj, warbleFrameIdx)
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
        fileIdxStart = fii;
      end
      if (frameEnd>prev) && (frameEnd<=fileFrameIdx(fii))
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
    otherCount = otherCount+1; 
    otherMeta(otherCount).count = otherCount;
    otherMeta(otherCount).fs = fs;
    otherMeta(otherCount).filenames = fileInRange;
    otherMeta(otherCount).timeStartRel = timeStartRel;
    otherMeta(otherCount).timeEndRel = timeEndRel;
    otherMeta(otherCount).duration = timeEndRel-timeStartRel;
    otherMeta(otherCount).segStartRel = segStart;
    otherMeta(otherCount).segEndRel = segEnd;
    otherMeta(otherCount).segStartAbs = segStartAbs;
    otherMeta(otherCount).segEndAbs = segEndAbs;
    otherMeta(otherCount).syllableOnsetsTime = (onsets(jj)-frameStart)*dt;
    otherMeta(otherCount).syllableOffsetsTime = (offsets(jj)-frameStart)*dt;
    % when does the warble start, this is just an estmiate
    temp = strsplit(fileInRange(1).name, '_');
    fileTime = sprintf('%s_%s', temp{end-1}, strrep(temp{end}, '.rhd', ''));
    otherMeta(otherCount).absTimeEst = datetime(fileTime,"InputFormat",'yyMMdd_HHmmss') + seconds(timeStartRel);
    % save audio as wav if requested
    if ~isempty(saveWavPath)
      signalAll = [];
      for ii=fileIdxStart:fileIdxEnd
        d_intan = read_Intan_RHD2000_file_to_struct_2_noNotch(fileStructAll(ii).folder, fileStructAll(ii).name, 0);
        % focus on selected channels
        signalThis = d_intan.board_adc_data(channels,:);
        fs = d_intan.frequency_parameters.board_adc_sample_rate;
        signalThis = signalThis';
        signalAll = [signalAll; signalThis];
      end
      signalOther = signalAll(segStart:segEnd, :);
      dstr = datestr(datetime(otherMeta(otherCount).absTimeEst), 'yyyymmddHHMMSS');
      fn_save = fullfile(saveWavPathOther, sprintf('%s_other_%05d_%s.wav', expID, otherCount, dstr));
      audiowrite(fn_save, signalOther/10, fs);
    end
  end
end
end