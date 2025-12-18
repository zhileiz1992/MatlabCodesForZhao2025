% Ensemble model to predict neuronal firing rate from VAE and latency to syllable onset/offset
% Step 1: prepare data to train XGBoost model in python
% Zhilei, 10/05/2025


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
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';


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
fd_save_base = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel');



%% 1. Loop through neuron, prepare predictor and firing rate data
% predictor: 32 VAE latents + latency to onset + latency to offset + rel latency to onset + rel latency to offset
% group id is the syl_ID
fs = 20000;
% window size and hop length, unit is frame (ms)
model_win = 20;  
model_hop = 5; 
lag = 20;  % the lag between the center of the firing rate window and the center of the spectrogram sliding window
% parameters for calculate IFR
ifr_sigma = 0.003;  % unit is sec
ifr_binsz = 0.001; 
% format of the data table
varNames = {'syl_ID', 'category', 'label', 'i', 't_onset', 't_offset', 'rel_onset', 'rel_offset', 't_vae'};
temp = compose('vae%d', 1:32);
varNames = [varNames temp 'ifr' 'spike_count'];
varTypes = {'string', 'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double'};
temp = repmat({'double'}, 1, 32);
varTypes = [varTypes temp 'double' 'double'];

% loop through neurons, assemble predictor variables
ni = 44; 

neuronID = info_neu.neuronID{ni};
spike_this = spike_win(strcmp(spike_win.neuronID, neuronID), :);
% what syllable types to analyze
calls = spike_this(strcmp(spike_this.category, 'v'),:);
comp = spike_this((ismember(spike_this.category, {'b','x'})) & (spike_this.dur>=0.3),:);
spike = [calls; comp];

d_all = cell(size(spike,1), 1);
parfor ri=1:size(spike, 1)
   % contruct the spike data
   total_len = spike.seg_end(ri) - spike.seg_start(ri) +1;  % total len of the padded segment
   rel_ori = spike.seg_start_ori(ri) - spike.seg_start(ri) +1;  % index of syllable onset
   syl_len = spike.seg_end_ori(ri) - spike.seg_start_ori(ri) +1;  % len of the syllable
   dur = syl_len / fs;
   spike_iv = zeros(total_len, 1);
   if ~isempty(spike.spike_i{ri})
     spike_iv(spike.spike_i{ri}) = 1;
   end
   rt_spike = ((1:length(spike_iv))-rel_ori) / fs * 1000;
   
   % calculate IFR
   [ifr, kt] = ZZfunc_IFRwithGaussian_v2(spike_iv, fs, ifr_sigma, ifr_binsz);
   % anything smaller than 0.1Hz is set to 0
   ifr(ifr<0.01) = 0;
   % convert time relative to syllable onset, unit is ms
   rt_kt = (kt - (rel_ori-1) / fs)* 1000;
   
   % grab the VAE latents
   sID = spike.syl_ID{ri};
   vae = h5read(fn_vae, ['/' sID]);
   vae = vae';
   % calcualte relative time of each VAE rows, consider the center of the sliding window
   rt_vae = (1:size(vae,1)) - 1 - win_frame/2;
   
   % loop through time points in the IFR, determine predictor x variables
   % define the time grids of the models: start with half of the model win
   t_grid = (rt_kt(1) + model_win/2) : model_hop: (rt_kt(end) - model_win/2);
   d_this = table('Size', [length(t_grid) size(varNames,2)],'VariableTypes', varTypes, 'VariableNames', varNames);
   for kt_i=1:length(t_grid)
     t = t_grid(kt_i);
     d_this.syl_ID{kt_i} = sID; 
     d_this.category{kt_i} = spike.category{ri};
     d_this.label{kt_i} = spike.label{ri};
     d_this.i(kt_i) = kt_i; 
     % timing predictors
     d_this.t_onset(kt_i) = t; % abs time to syllable onset
     d_this.t_offset(kt_i) = dur*1000 - t;  % abs time to syllable offset
     d_this.rel_onset(kt_i) = t / (dur*1000); % relative time to syllable onset
     d_this.rel_offset(kt_i) = d_this.t_offset(kt_i) / (dur*1000);
     % VAE latent predictors
     t_v = t + lag;
     d_this.t_vae(kt_i) = t_v; 
%      i_v = find(rt_vae==t_v);
%      [~, i_v] = min(abs(rt_vae-t_v));
     if (t_v<rt_vae(1)) || (t_v>rt_vae(end))
       v_this = nan(1, 32);
     else
       [~, i_v] = min(abs(rt_vae-t_v));
       v_this = vae(i_v,:);
     end
     d_this(kt_i, 10:41) = num2cell(v_this);
     % what to be predicted: IFR, take the average of the window range
     i_e = find((rt_kt>=(t-model_win/2)) & (rt_kt<=(t+model_win/2)));
     d_this.ifr(kt_i) = nanmean(ifr(i_e)); 
     % also count the number of spikes in the window
     i_s = find((rt_spike>=(t-model_win/2)) & (rt_spike<=(t+model_win/2)));
     d_this.spike_count(kt_i) = sum(spike_iv(i_s));
     
   end
   
   % add to the master table
   d_all{ri} = d_this;
end

% save into disk for future use
fd_save_this = fullfile(fd_save_base, 'CuratedData');
if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end
fn_mat = fullfile(fd_save_this, sprintf('%s.%s.call_comp.modelData.mat', birdID, neuronID));
save(fn_mat, 'd_all');



% plot to check examples (optional)
%   ri = 1;
%  close all; figure; 
%  ax1 = subplot(2,1,1);
%  plot(ax1, spike_iv); hold on; 
%  xline(ax1, rel_ori, 'LineStyle', '--', 'Color', 'green');
%  xline(ax1, spike.seg_end_ori(ri) - spike.seg_start(ri) +1, 'LineStyle', '--', 'Color', 'red');
%  xlim(ax1, [1 total_len]);
%  ax2 = subplot(2,1,2);
%  plot(ax2, rt_kt, ifr); 
%  xlim(ax2, [rt_kt(1) rt_kt(end)]);








