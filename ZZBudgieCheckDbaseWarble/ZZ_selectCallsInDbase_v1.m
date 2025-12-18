% given a dbase with syllable segmentations (manual, flatness, or WhisperSeg)
% select contact call-like syllables, so we can train a discrete VAE to
% differentiate different call types
% 2024/02/15

clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes'));

%% inputs
fd_z4 = '/mnt/z4/';
birdID = 'pair1CU21RigB';
data_date = '2023-11-05';
fd_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles', 'pair1', data_date, birdID, 'warble');
% where to save results
fd_save = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'Figures', 'pair1', data_date, birdID, 'warble');
if ~exist(fd_save)
  mkdir(fd_save);
end
fn_dbase = fullfile(fd_dbase, 'pair1CU21RigB.20231104.warble.labelZZ.dbase.mat');
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


%% plot selected segments to check labeling
% find all 'w' segments: contact calls
close all;
target_syl = 'w'; 
s_idx = find(strcmp(target_syl, {segments.title}));
seg_selected = segments(s_idx);
maxDur = 0.3; 
save_folder = fullfile(fd_save, target_syl);
prefix = sprintf('%s_%s', birdID, target_syl);
plotSpectrogramListZZ_flexible_v1(seg_selected, 'signalFiltered', 1, maxDur, 3, 10, save_folder, prefix, 'title_str', [500, 7500], [8, 20], ); 
% plotSpectrogramListZZ_flexible_v1(segments, 'signalFiltered', 1, maxDur, 3, 5, save_folder, suffix, 'title_str',0.1, 0.15); 
% save the extracted data
fn_save = fullfile(save_folder, sprintf('%s_target_syl_%s.mat', birdID, target_syl));
save(fn_save, 'seg_selected');


%% find ephys trace and sortd spikes
d_anno = dbase; 
ech = 8; % what ephys 
% is it a top spike (1) or bottom spike (2)
spike_shape = 2; 
fn_dbase_ephys = fullfile(fd_dbase, 'pair1CU21RigB.20231104.warble.chan8.dbase.mat');
load(fn_dbase_ephys); 
d_ephys = dbase; 
spike_all = d_ephys.EventTimes{1, 1};
spike_selected =  d_ephys.EventIsSelected{1, 1}; % note that spikes may be manually excluded
% go through selected syllables, find ephys and spiking data
ephys_selected = seg_selected;
for si=1:length(ephys_selected)
%   si = 5; 
  fn_ephys = strrep(ephys_selected(si).fn_audio, 'chan0', sprintf('chan%d', ech));
  ephys = ncread(fn_ephys, 'data');
  ephys_selected(si).e_trace = ephys(ephys_selected(si).seg_start:ephys_selected(si).seg_end); 
  % find the spikes within the range
  [a,b,c] = fileparts(ephys_selected(si).fn_audio);
  f_idx = find(strcmp([b c], {d_ephys.SoundFiles.name})); 
  % spikes need to pass threshold and not manually excluded
  spike_pass = spike_all{spike_shape,f_idx};
  spike_sel = spike_selected{spike_shape,f_idx};
  spike_t = spike_pass(spike_sel==1);
  spike_iv = zeros(size(ephys));  % indicator variable for spikes
  spike_iv(spike_t) = 1; 
  ephys_selected(si).spike_iv = spike_iv(ephys_selected(si).seg_start:ephys_selected(si).seg_end);
end


%% overlay spikes on spectrograms
save_folder = fullfile(fd_save, [target_syl '_ephys']);
plotSpectrogramListSpikes_v1(ephys_selected, 'signalFiltered', 1, maxDur, 3, 10, save_folder, prefix, 'title_str', [500, 7500], [8, 20]); 
% save the extracted data
fn_save = fullfile(save_folder, sprintf('%s_target_ephys_%s.mat', birdID, target_syl));
save(fn_save, 'ephys_selected');


%% select contact calls based on duration and frequency range
criteria.minDur = 0.1;
criteria.maxDur = 0.35;
criteria.minFreq = 2000;
criteria.maxFreq = 4000;
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
% determine what cutoff values to use
% plot distribution of freq ratio: calls vs non-calls in the labeled
% not all files are labeled
labeled_files = {d_anno.SoundFiles(1:100).name}; 
idx_call = [];
idx_other = [];
for si=1:length(segments)
  [a,b,c] = fileparts(segments(si).fn_audio);
  if ismember([b c], labeled_files)
    if strcmp(segments(si).title, 'w') || strcmp(segments(si).title, 'v')
      idx_call = [idx_call si];
    else
      idx_other = [idx_other si];
    end
  else
    continue
  end
end
% compare the distribution
figure; 
x = [segments(idx_call).freqRatio]; 
y = [segments(idx_other).freqRatio]; 
xlimits = [min([x y]), max([x y])];
subplot(2,1,1); % 1 row, 2 columns, first subplot
histogram(x);
title('Distribution of call-like');
subplot(2,1,2); % 1 row, 2 columns, first subplot
histogram(y);
title('Distribution of other vocalizations');
xlabel('Freq. ratio: 2-4k vs. other', 'FontSize', 16);

% seems that freqRatio 2-3 is a good choice?
criteria.minFreqRatio = 2; 
% calculate the true positive and false positive rate
seg_labeled = segments(sort([idx_call idx_other]));
idx_pass = find(([seg_labeled.duration]>=criteria.minDur) & ([seg_labeled.duration]<=criteria.maxDur) & ([seg_labeled.freqRatio]>=criteria.minFreqRatio));
% true positive
true_pos = length(intersect(idx_call, idx_pass))/length(idx_pass);
false_pos = 1 - true_pos; 
% false negative
false_neg = 1-length(intersect(idx_call, idx_pass))/length(idx_call);
disp([criteria.minFreqRatio, false_pos, false_neg]);
% which manually labeled calls are not included
a = setdiff(idx_call, idx_pass);


%% select all call-like vocalizations in the dbase
idx_pass_all = find(([segments.duration]>=criteria.minDur) & ([segments.duration]<=criteria.maxDur) & ([segments.freqRatio]>=criteria.minFreqRatio));
% plot and check 
save_folder = fullfile(fd_save, 'auto_call'); 
plotSpectrogramListZZ_flexible_v1(segments(idx_pass_all), 'signalFiltered', 1, maxDur, 3, 10, save_folder, 'auto_call', 'title_str', [500, 7500], [8, 20]);
% save the extracted data
fn_save = fullfile(save_folder, sprintf('%s_%s_auto_call.mat', birdID, '2023-11-04'));
auto_calls = segments(idx_pass_all);
save(fn_save, 'auto_calls', '-v7.3');


%% export as h5 file to train the VAE
fd_h5 = fullfile(fd_save, 'auto_call_h5');
spec_field = 'spec'; 
spec_size = [128, 128]; 
prefix = sprintf('%s_%s_autocalls', birdID, '2023-11-04');
fd_h5 = ZZ_exportAsH5Func_v1(auto_calls, spec_field, fd_h5, spec_size, prefix);



