% given a dbase with syllable segmentations (manual, flatness, or WhisperSeg)
% audotmatically select contact call-like syllables, so we can train 
% a discrete VAE to cdifferentiate different call types
% 2024/02/15

clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes'));

%% inputs
fd_z4 = '/mnt/z4/';
pairID = 'pair1'; 
birdID = 'pair1CU21RigB';
data_dates = {'2023-11-04', '2023-11-05', '2023-11-08', '2023-11-12'};
data_date = data_dates{4};
fd_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles', pairID, data_date, birdID, 'warble');
% where to save results
fd_save = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'Figures', pairID, data_date, birdID, 'warble');
if ~exist(fd_save)
  mkdir(fd_save);
end
fn_dbase = fullfile(fd_dbase, sprintf('%s.%s.warble.focal.dbase.mat', birdID, strrep(data_date,'-','')));
load(fn_dbase);

%% extract sound data for the segments
% where original NC files are saved, need to change the path format, since
% sorting is done in Win 10, but analysis is in Ubuntu system
nc_path = ZZ_winPathToLinux_v1(dbase.PathName, 'Y:', '/mnt/z4');
% loop through each file
segments = [];
seg_count = 0;
% fi = 3;
for fi=1:length(dbase.SoundFiles)
  fn_audio = fullfile(nc_path, dbase.SoundFiles(fi).name);
  audio = ncread(fn_audio, 'data');
  audio = double(audio/10);  % divide by range of Intan AI, convert to double format
  fs = dbase.Fs;
  % loop through each segment
  syl_time = dbase.SegmentTimes{1, fi};
  syl_list = dbase.SegmentTitles{fi};
  
  for si=1:size(syl_time)
    seg_count = seg_count+1;
    segments(seg_count).title = syl_list{si};
    segments(seg_count).title_str = sprintf('%d-%d-%d-%s', seg_count, fi, si, syl_list{si});
    segments(seg_count).fn_audio = fn_audio;
    segments(seg_count).seg_start = syl_time(si,1);
    segments(seg_count).seg_end = syl_time(si,2);
    segments(seg_count).fs = fs;
    segments(seg_count).signalRaw = audio(segments(seg_count).seg_start : segments(seg_count).seg_end);
    segments(seg_count).duration = length(segments(seg_count).signalRaw) / fs;
    % band pass the signal
    segments(seg_count).signalFiltered = ZZ_butterFilterFunc_v1(segments(seg_count).signalRaw, fs, 500, 10000, 2);
    % normalize amplitude
    segments(seg_count).signalNorm = audioNormalization_YW(segments(seg_count).signalFiltered, 0.25);
    % calculate spectrogram
    [power, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(segments(seg_count).signalFiltered, fs);
    segments(seg_count).spec = powerGrey;
    segments(seg_count).spec_colored = power;
    segments(seg_count).spec_raw = powerRaw;
    segments(seg_count).f = f;
    segments(seg_count).t = t; 
  end
end


%% select contact calls based on duration and frequency range
criteria.minDur = 0.1;
criteria.maxDur = 0.35;
criteria.minFreq = 2000;
criteria.maxFreq = 4000;
criteria.minFreqRatio = 2;
% calculate the freq ratio betwen the 2-4k and other ranges
for si=1:length(segments)
  % si=10;
  % find the frequency index within range
  fidx_in = find((segments(si).f>=criteria.minFreq) & (segments(si).f)<=criteria.maxFreq);
  fidx_out = find((segments(si).f<criteria.minFreq) | (segments(si).f)>criteria.maxFreq);
  spec = segments(si).spec;
  mean_freq_in = mean(mean(spec(fidx_in, :)));
  mean_freq_out = mean(mean(spec(fidx_out, :)));
  segments(si).freqRatio = mean_freq_in / mean_freq_out;
end
% find all syllables that meet the criteria 
idx_pass_all = find(([segments.duration]>=criteria.minDur) & ([segments.duration]<=criteria.maxDur) & ([segments.freqRatio]>=criteria.minFreqRatio));
% plot and check 
save_folder = fullfile(fd_save, 'auto_call'); 
maxDur = 0.3; 
plotSpectrogramListZZ_flexible_v1(segments(idx_pass_all), 'signalFiltered', 1, maxDur, 3, 10, save_folder, 'auto_call', 'title_str', [500, 7500], [8, 20]);
% save the extracted data
fn_save = fullfile(save_folder, sprintf('%s_%s_auto_call.mat', birdID, data_date));
auto_calls = segments(idx_pass_all);
save(fn_save, 'auto_calls', '-v7.3');


%% export as h5 file to train the VAE
for di=3:4
  data_date = data_dates{di}
  fd_save = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'Figures', pairID, data_date, birdID, 'warble');
  fn_auto_call = fullfile(fd_save, 'auto_call', sprintf('%s_%s_auto_call.mat', birdID, data_date));
  load(fn_auto_call);
  fd_h5 = fullfile(fd_save, 'auto_call_h5');
  spec_field = 'spec'; 
  spec_size = [128, 128]; 
  prefix = sprintf('%s_%s_autocalls', birdID, data_date);
  fd_h5 = ZZ_exportAsH5Func_v1(auto_calls, spec_field, fd_h5, spec_size, prefix);
end



