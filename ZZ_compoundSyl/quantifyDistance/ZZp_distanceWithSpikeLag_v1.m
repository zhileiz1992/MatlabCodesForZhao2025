% quantity acoustic distance on VAE latents with different time lags to the spike onset
% Zhilei, 10/03/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/overlaySpikes'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl/compoundTraj'));


%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what's the window size
win_frame = 32;
ms_per_frame = 1;
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
% what's the range to include spikes: default only consider spikes within 1-window of syllable boundaries
epad_frame = win_frame;


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird specific folder setting
% where is VAE data located
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', sprintf('latents.%s', vae_run));
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));
  
% information regarding spikes and sliding windows (zero lags)
fd_field = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_firingField']);
fn_spike = fullfile(fd_field, sprintf('%s.spike_win.mat', birdID));
load(fn_spike);

% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% where to save results
fd_save = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_distance']);
if ~exist(fd_save, 'dir'); mkdir(fd_save); end



%% Loop through neurons, calculate acoustic distance with different lags to spike
lags = -100:5:100;

ni=44;
neuronID = info_neu.neuronID{ni};
spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
% only look at compound syllables 
spike_this = spike_this((ismember(spike_this.category, {'b','x'})) & (spike_this.dur>=0.3),:);
% only look at syllables that have spikes
idx = cellfun(@(x) ~isempty(x), spike_this.win_loc);
spike_has = spike_this(idx,:);


%% 1. Grab the VAE latents with different lags to spikes
% loop through syllables, then loop through spikes, grab the VAE data with different lags
d_vae = cell(size(spike_has,1),1);
parfor ri=1:size(spike_has,1)
    % load the vae latent 
    sID = spike_has.syl_ID{ri};
    vae = h5read(fn_vae, ['/' sID]);
    vae = vae';
    % for each spike, determine what window index of vae 
    % assign to the first or last window if out of range
%     win_loc = spike_has.win_loc{ri};
    win_loc = spike_has.win_loc{ri} + win_frame;  % convert to the index in the vae matrix
    [res, i_all] = ZZfunc_getVaeWithLag_v2(vae, win_loc, lags);
    d_vae{ri} = res;
end
win_count = cellfun(@(x) size(x,1), d_vae);

% stack into a large 3d matrix
d_vae = cat(1, d_vae{:});
% get the start/end index of each syllable in this 3d matrix
cum_count = cumsum(win_count);
spike_has.vstart = [1; cum_count(1:(end-1))+1];
spike_has.vend = cum_count;


%% 2. Calculate then plot acoustic distance
num_round = 20;  % how many round of sampling to take
dist_all = nan(length(lags), num_round);
rng(1992);
for ti=1:length(lags)
  d_this = squeeze(d_vae(:,ti,:));
  for round_i=1:num_round
    dist = ZZfunc_pairwiseDist_v2(d_this, 1000, 'cosine');
    dist_all(ti, round_i) = dist;
  end
end
% dist_mean = mean(dist_all, 2);
% plot results
close all; 
fig = ZZfunc_newFigurePDFsize_v1([10 10 500 500]);
ax = gca; hold(ax, 'on');
plot_mean_with_ci(ax, lags, dist_all', [0 0 0], 2, '-', [0.5 0.5 0.5], 0.25);
xlabel('Time lag from spike (ms)', 'FontSize', 14);
ylabel('Acoustic distance (cosine)', 'FontSize', 14);

fn_fig = fullfile(fd_save, sprintf('%s.%s.comp_dist.pdf', birdID, neuronID));
print(fig, fn_fig, '-dpdf', '-painters');








