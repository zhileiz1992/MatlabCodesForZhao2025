% try using GLM to predict the firing rate
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



%% Loop through neurons, build GLM model
ni = 44;
neuronID = info_neu.neuronID{ni};

% load the curated data
fd_save_this = fullfile(fd_save_base, 'CuratedData');
fn_d = fullfile(fd_save_this, sprintf('%s.%s.call_comp.modelData.mat', birdID, neuronID));
a = load(fn_d);
d_all = a.d_all;

% concatenate into a large table
d = vertcat(d_all{:});
% remove those that have nan in the VAE latents
d_m = d(~isnan(d.vae1),:);



%% Build the GLM model: try on different inputs
% on both calls and compound syllables
d_comp = d_m(ismember(d_m.category, {'b', 'x'}),:);
d_call = d_m(strcmp(d_m.category, 'v'),:);
d_v1 = d_m(strcmp(d_m.label, 'v1'),:);

d_list = {d_m, d_comp, d_call, d_v1};
suffix = {'Calls+Compounds', 'Compounds', 'Calls', 'V1'};

r2_all = cell(size(d_list,2), 1);
for type_i = 1:size(d_list,2)
  d_plot = d_list{type_i};
  
  % how many times to repeat
  n = 20;
  rng(1992);
  r2 = nan(n, 3);
  pred_time = {'t_onset', 't_offset', 'rel_onset'};
  pred_vae = compose('vae%d', 1:32);
  
  for ii=1:n
    
    % full model
    pred_names = [pred_time pred_vae];
    % output_name = 'ifr';
    output_name = 'spike_count';
    group_name = 'syl_ID';
    % dist = 'gamma';
    dist = 'poisson';
    link = 'log';
    train_test_ratio = 0.8;
    [mdl, metrics, coefTbl, testPredTbl] = glm_grouped_pipeline2(d_plot, pred_names, output_name, group_name, dist, link, train_test_ratio, '');
    
    % Acoustic only
    [mdl2, metrics2, coefTbl2, testPredTbl2] = glm_grouped_pipeline2(d_plot, pred_vae, output_name, group_name, dist, link, train_test_ratio, '');
    
    % Time only
    pred_names3 = {'t_onset', 't_offset', 'rel_onset'};
    [mdl3, metrics3, coefTbl3, testPredTbl3] = glm_grouped_pipeline2(d_plot, pred_time, output_name, group_name, dist, link, train_test_ratio, '');
    
    % r2(ii,:) = [metrics.R2 metrics2.R2 metrics3.R2];
    r2(ii,:) = [metrics.R2_dev_test metrics2.R2_dev_test metrics3.R2_dev_test];
    
  end
  
  r2_all{type_i} = r2;
end


%% plot results
close all;
[fig, axs] = generatePanelGrid_v2(1, 4, [0.75], [], [0.075;0.05], [0.075;0.05], 0.075, [0], [10 10 1600 500]);
for type_i=1:size(d_list,2)
  ax = axs(1, type_i);
  r2 = r2_all{type_i};
%   axes(ax);
%   violinplot(r2, 'Parent', ax);
  boxplot(ax, r2);
  xticks(ax, [1 2 3]);
  xticklabels(ax, {'Full', 'Acoustic only', 'Time only'});
  ylabel(ax, 'Deviance explained', 'FontSize', 14);
  title(ax, sprintf('GLM: %s', suffix{type_i}), 'FontSize', 14);
end










