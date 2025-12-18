% for neurons that are already sorted, manually pick a bunch similar
% syllables, then align the spikes 
% Zhilei Zhao, 04/29/2024

clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes'));

%% inputs
fd_z4 = '/mnt/z4/';
birdID = 'pair2RigBCU25';
pairID = 'pair2CU20CU25';
data_date = '2024-01-19'; ch = 16; sorted_range = [1:120];
% data_date = '2023-11-08'; ch = 15; sorted_range = [1:36];
fd_data_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles', pairID, data_date, birdID, 'warble');
fd_save_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'Figures', pairID, data_date, birdID, 'warble');
% is it a top spike (1) or bottom spike (2)
spike_shape = 2; 
% load the dbase with sorted spikes
fn_dbase = fullfile(fd_data_base, sprintf('%s.%s.warble.good.ch%d.dbase.mat', birdID, strrep(data_date,'-',''), ch));
load(fn_dbase); 
% where to save
fd_save = fullfile(fd_save_base, 'manual_picked');
if ~exist(fd_save)
  mkdir(fd_save)
end


%% extract sound data for the selected syllables
syl_label = 'z';
pad = 0.075;  % how many seconds to pad before and after the syllables
% where original NC files are saved, need to change the path format, since
% sorting is done in Win 10, but analysis is in Ubuntu system
nc_path = ZZ_winPathToLinux_v1(dbase.PathName, 'Y:', '/mnt/z4');
% get all the spiking time and selection data
spike_all = dbase.EventTimes{1, 1};
spike_selected =  dbase.EventIsSelected{1, 1}; % note that spikes may be manually excluded
% loop through each file, extract data for the segments
segments = [];
seg_count = 0;
for fi=sorted_range
  % find where the selected syllables are located
  % loop through each segment
  syl_time = dbase.SegmentTimes{1, fi};
  syl_list = dbase.SegmentTitles{fi};
  s_idx = [];
  for si=1:length(syl_list)
    if ismember(syl_label, syl_list{si})
      s_idx = [s_idx si];    
    end
  end
  if ~isempty(s_idx)
    % get the audio data
    fn_audio = fullfile(nc_path, dbase.SoundFiles(fi).name);
    audio = ncread(fn_audio, 'data');
    audio = double(audio/10);  % divide by range of Intan AI, convert to double format
    fs = dbase.Fs;
    % get the ephys data
    fn_ephys = strrep(fn_audio, 'chan0', sprintf('chan%d', ch));
    ephys = ncread(fn_ephys, 'data');
    for si=s_idx
      seg_count = seg_count+1;
      segments(seg_count).title = syl_list{si};
      segments(seg_count).title_str = sprintf('%d-%d-%d-%s', seg_count, fi, si, syl_list{si});
      segments(seg_count).fn_audio = fn_audio;
      segments(seg_count).fs = fs;
      padPoints = floor(fs*pad);
      % calculate the segment start and end
      segments(seg_count).seg_start_ori = syl_time(si,1);
      segments(seg_count).seg_end_ori = syl_time(si,2);
      segments(seg_count).seg_start = max([1 syl_time(si,1)-padPoints]);
      segments(seg_count).seg_end = min([size(ephys,1) syl_time(si,2)+padPoints]);
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
      % get the ephys trace and spikes
      segments(seg_count).e_trace = ephys(segments(seg_count).seg_start:segments(seg_count).seg_end); 
      % spikes need to pass threshold and not manually excluded
      spike_pass = spike_all{spike_shape,fi};
      spike_sel1 = spike_selected{1,fi};
      spike_sel2 = spike_selected{2,fi};
      spike_t = spike_pass(spike_sel1==1 & spike_sel2==1);
      spike_iv = zeros(size(ephys));  % indicator variable for spikes
      spike_iv(spike_t) = 1; 
      segments(seg_count).spike_iv = spike_iv(segments(seg_count).seg_start:segments(seg_count).seg_end);
    end
  end
end


%% Manually pick some syllables, plot ephys voltage trace and sorted by number of spikes
% sort by number of spikes
p = 1:length(segments);
ephys_this = segments(p); 
% num_spikes = cellfun(@sum, {ephys_this(:).spike_iv});
% [a, sort_i] = sort(num_spikes, 'descend');
% ephys_this = ephys_this(sort_i); 
fd_save_this = fullfile(fd_save, syl_label); 
prefix = syl_label;
plotSpectrogramListSpikes_v1(ephys_this, 'signalFiltered', 1, 0.35, 3, 6, fd_save_this, prefix, 'title_str', [500, 7500], [8, 20], 0.13); 

% plot raw traces as well
% FIR pass the raw trace
for si=1:length(ephys_this)
  ephys_this(si).e_trace_FIR = ZZ_FIRbandpass(double(ephys_this(si).e_trace), ephys_this(si).fs,  400, 9000, 80);
end
fd_save_this = fullfile(fd_save, 'p1_sep'); 
plotSpectrogramListTraceSpikes_v2(ephys_this, 'signalFiltered', 1, 0.4, 2, 6, fd_save_this, prefix, '', [500, 7500], [8, 20], 0.13); 

