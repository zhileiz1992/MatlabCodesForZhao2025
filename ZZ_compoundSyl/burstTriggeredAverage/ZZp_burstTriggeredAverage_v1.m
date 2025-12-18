% Goal: calculate burst-triggered average spectrogram for compound syllables

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl/quantifyDistance'));


%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
clims = {[10 21]; [10 21]; [10.5 21.5]; [12 23]};
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';

% a custom colormap for spectrogram
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


%% Bird-specific setting
bi = 2;
birdID = birdIDs{bi};
pairID = pairIDs{bi};

% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% information regarding spikes and sliding windows (zero lags)
fd_field = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_firingField']);
fn_spike = fullfile(fd_field, sprintf('%s.spike_win.mat', birdID));
load(fn_spike);
% subset to compound syllables
comp = spike_win((ismember(spike_win.category, {'b','x'})) & (spike_win.dur>=0.3),:);

% where to save results
fd_save = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_burstTriggerAveSpec']);
if ~exist(fd_save)
  mkdir(fd_save);
end

% loop through neurons
% ni = 24;
for ni=1:size(info_neu, 1)
neuronID = info_neu.neuronID{ni};
spike_this = comp(strcmp(comp.neuronID, neuronID),:);
% only look at syllables that have spikes
idx = cellfun(@(x) ~isempty(x), spike_this.win_loc);
spike_has = spike_this(idx,:);

% Find the location of the burst
% try first spike of the burst or center of mass
% burst_method = 'first';
burst_method = 'center';
% convert gap to data points
fs = 20000; 
max_gap = 10 / 1000 * fs;
spike_has.burst_loc = cellfun(@(x) ZZfunc_findBurstLoc_v1(x, max_gap, burst_method), spike_has.spike_in, 'UniformOutput', false);


% loop through syllables, read sound file, align to burst center, save audio data
pad_pre = 0.1; 
pad_post = 0.2; 
seg_d_pre = nan(floor(pad_pre*fs), 1);
seg_d_post = nan(floor(pad_post*fs), 1);
seg_all = [];
for si=1:size(spike_has, 1)
  % read audio
  fn_a = spike_has.fn_audio{si};
  signal = double(ncread(fn_a, 'data'));
  seg_start = spike_has.seg_start(si);  % rel. origin of syllable segment
  burst_loc = spike_has.burst_loc{si};
  % loop through burst locs
  for bi=1:size(burst_loc,1)
    seg_d = [seg_d_pre; seg_d_post];
    b_center = burst_loc(bi) + seg_start;
    % calculate the start and end
    i_start = max([1 b_center-size(seg_d_pre,1)]);
    i_end = min([size(signal,1) b_center+size(seg_d_post,1)-1]);
    % calculate an offset
    offset = size(seg_d_pre,1) - (b_center - i_start);
%     offset = max([1 offset]);
    seg_d((offset+1):(i_end-i_start+offset+1), 1) = signal(i_start:i_end,1);
    seg_all = [seg_all; seg_d'];
  end
end


% calculate spectrograms then average
% ignore rows that have NaN (at file boundaries)
seg_all_clean = seg_all(~any(isnan(seg_all), 2), :);
spec_aligned = cell(size(seg_all_clean,1),1);
for si=1:size(seg_all_clean, 1)
  % calculate spectrogram
  sound = seg_all_clean(si,:);
  [power, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(sound, fs, 256, 256, 236, [250 7500], clims{bi});
  spec_aligned{si} = powerGrey/1024;  % normalize
end

% limit to certain frequency range
freq_range = [500 5500];
freq_i = find((f>=freq_range(1)) & (f<=freq_range(2)));


% plot some exmaples to check
% close all;
% [fig, axs] = generatePanelGrid_v2(1, 10, [0.75], [], [0.08;0.02], [0.1;0.05], 0.01, [0], [50 50 1200 300]);
% rel_t = t - t(1) - pad_pre;
% for si=1:10
%   ax = axs(1,si);
%   spec = spec_aligned{si};
%   spec = spec(freq_i, :);
%   imagesc(ax, rel_t, f(freq_i), spec, [0, 1]);
%   set(ax, 'YDir', 'Normal');
%   colormap(ax, custom_map);
% end

% calculate and plot averaged spectrogram
close all;
stacked_array = cat(3, spec_aligned{:});
stacked_array = stacked_array(freq_i, :, :);
spec_mean = mean(stacked_array, 3);
fig2 = ZZfunc_newFigurePDFsize_v1([50 50 700 350]);
ax = gca;
imagesc(ax, rel_t, f(freq_i), spec_mean, [0.1 max(spec_mean(:))]);
set(ax, 'YDir', 'Normal');
set(ax, 'YTick', [2000 4000]);
set(ax, 'YTickLabel', {'2k', '4k'}, 'FontSize', 10);
colormap(ax, 'jet');
% colormap(ax, custom_map);
title(ax, sprintf('%s %d syls, %d bursts', neuronID, size(spike_this,1), size(seg_all_clean, 1)));
fn_pdf = fullfile(fd_save, sprintf('%s.%s.%s.burstTriggeredAveSpec.pdf', birdID, neuronID, burst_method));
print(fig2, fn_pdf, '-dpdf', '-painters');

end








