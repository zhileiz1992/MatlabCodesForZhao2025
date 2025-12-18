% for neurons that are already sorted, given the dbase
% 1. automatically identify call-like vocalizations 
% 2. plot the spetrograms
% 3. manually select the calls belong to the same type
% 4. for those calls, align the spikes to the spectrograms
% Zhilei Zhao, 02/28/2024

clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes'));

%% inputs
fd_z4 = '/mnt/z4/';
birdID = 'pair1CU21RigB';
data_date = '2023-11-04'; ch = 8; sorted_range = [1:200];
fn = 'pair1CU21RigB.20231104.warble.partnerChan8.1-200.dbase.mat'; % dbase file name
fd_data_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles', 'pair1', data_date, birdID, 'warble');
fd_save_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'Figures', 'pair1', data_date, birdID, 'warble');
% is it a top spike (1) or bottom spike (2)
spike_shape = 2; 
% load the dbase with sorted spikes
fn_dbase = fullfile(fd_data_base, fn);
load(fn_dbase); 
% where to save
fd_save = fullfile(fd_save_base, 'partner_call');
if ~exist(fd_save)identifyVocalBird)(info)
  mkdir(fd_save)
end


%% extract sound data for the segments
% need to change the dbase path format, since
% sorting is done in Win 10, but analysis is in Ubuntu system
nc_path = ZZ_winPathToLinux_v1(dbase.PathName, 'Y:', '/mnt/z4');
% loop through each file
segments = [];
seg_count = 0;
% fi = 3;
for fi=sorted_range
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
idx_pass_all = find(([segments.duration]>=criteria.minDur) & ([segments.duration]<=criteria.maxDur) & ([segments.freqRatio]>=criteria.minFreqRatio));
seg_selected = segments(idx_pass_all);
% fix the title str to reflect syllable order in the new selected dbase
for si=1:length(seg_selected)
  temp = strsplit(seg_selected(si).title_str, '-');
  % replace with the new order id
  temp{1} = num2str(si); 
  seg_selected(si).title_str_new = strjoin(temp, '-');
end
% plot and check 
save_folder = fullfile(fd_save, 'auto_call'); 
maxDur = 0.4; 
plotSpectrogramListZZ_flexible_v1(seg_selected, 'signalFiltered', 1, maxDur, 3, 10, save_folder, 'auto_call', 'title_str_new', [500, 7500], [8, 20]);
% save the extracted data
fn_save = fullfile(save_folder, sprintf('%s_%s_auto_call.mat', birdID, data_date));
auto_calls = seg_selected;
save(fn_save, 'auto_calls', '-v7.3');


%% find the ephys trace and sorted spikes
d_ephys = dbase; 
spike_all = d_ephys.EventTimes{1, 1};
spike_selected =  d_ephys.EventIsSelected{1, 1}; % note that spikes may be manually excluded
% go through selected syllables, find ephys and spiking data
ephys_selected = seg_selected;
for si=1:length(ephys_selected)
%   si = 5; 
  [a,b,c] = fileparts(ephys_selected(si).fn_audio);
  temp = strsplit(b, '_');
  temp{end} = sprintf('chan%d', ch);
  fn_ephys = fullfile(a, [strjoin(temp, '_') c]);
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


%% Manually pick some syllables, plot ephys voltage trace and sorted by number of spikes
% sort by number of spikes
p1 = [8,20,30,39,44,49,55,65,99,105,129,137,140,172,187,196,216,227,251,256,266,272,275,282,298,326,331,343,403,409,411,425,447,451,462,467,472,476,482,485,487,503,512,513,522,524,528,539,550,555];
p = p1; 
ephys_this = ephys_selected(p); 
num_spikes = cellfun(@sum, {ephys_this(:).spike_iv});
[a, sort_i] = sort(num_spikes, 'descend');
ephys_this = ephys_this(sort_i); 
fd_save_this = fullfile(fd_save, 'selected');
prefix = sprintf('%s_%s', birdID, data_date);
plotSpectrogramListSpikes_v1(ephys_this, 'signalFiltered', 1, 0.3, 3, 10, fd_save_this, prefix, 'title_str_new', [500, 7500], [8, 20]); 
% plot raw traces as well
% FIR pass the raw trace
for si=1:length(ephys_this)
  ephys_this(si).e_trace_FIR = ZZ_FIRbandpass(double(ephys_this(si).e_trace), ephys_this(si).fs,  400, 9000, 80);
end
fd_save_this = fullfile(fd_save, 'selected_ephys'); 
plotSpectrogramListTraceSpikes_v2(ephys_this, 'signalFiltered', 1, 0.24, 2, 10, fd_save_this, prefix, '', [500, 7500], [8, 20]); 

