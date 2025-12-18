% try using LSBoost to predict the firing rate
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
% d_list = {d_comp, d_call, d_v1};
% suffix = {'Compounds', 'Calls', 'V1'};

r2_all = cell(size(d_list,2), 1);
for type_i = 1:size(d_list,2)
  d_plot = d_list{type_i};
  
  % how many times to repeat
  n = 10;
  rng(1992);
  r2 = nan(n, 3);
  pred_time = {'t_onset', 't_offset', 'rel_onset'};
  pred_vae = compose('vae%d', 1:32);
  train_test_ratio = 0.8;
  
  for ii=1:n
    fprintf('%d %d\n', type_i, ii);
    
    % full model
    pred_names = [pred_time pred_vae];
    output_name = 'ifr';
%     output_name = 'spike_count';
    group_name = 'syl_ID';
    % dist = 'gamma';
    [mdl, metrics, importTbl, testPredTbl] = lsboost_grouped_pipeline(d_plot, pred_names, output_name, group_name, train_test_ratio);
    
    % Acoustic only
    [mdl2, metrics2, importTbl2, testPredTbl2] = lsboost_grouped_pipeline(d_plot, pred_vae, output_name, group_name, train_test_ratio);
    
    % Time only
    [mdl3, metrics3, importTbl3, testPredTbl3] = lsboost_grouped_pipeline(d_plot, pred_time, output_name, group_name, train_test_ratio);
    
    % r2(ii,:) = [metrics.R2 metrics2.R2 metrics3.R2];
    r2(ii,:) = [metrics.R2_SSE metrics2.R2_SSE metrics3.R2_SSE];
    
  end
  
  r2_all{type_i} = r2;
end


%% plot results
close all;
[fig, axs] = generatePanelGrid_v2(1, 4, [0.75], [], [0.075;0.05], [0.075;0.05], 0.075, [0], [10 10 1500 500]);
for type_i=1:size(d_list,2)
  ax = axs(1, type_i);
  r2 = r2_all{type_i};
%   axes(ax);
%   violinplot(r2, 'Parent', ax);
  boxplot(ax, r2);
  xticks(ax, [1 2 3]);
  xticklabels(ax, {'Full', 'Acoustic only', 'Time only'});
  ylabel(ax, 'Deviance explained', 'FontSize', 14);
  title(ax, sprintf('LSBoost: %s', suffix{type_i}), 'FontSize', 14);
end

r2 = r2_all{3};
[p12,~,~] = signrank(r2(:,1), r2(:,2));
[p13,~,~] = signrank(r2(:,1), r2(:,3));


 
close all;
[fig, axs] = generatePanelGrid_v2(1, 3, [0.75], [], [0.075;0.05], [0.075;0.05], 0.075, [0], [10 10 1500 500]);
test_all = {testPredTbl, testPredTbl2, testPredTbl3};
t_str = {'Full', 'Acoustic only', 'Time only'};
for type_i=1:size(test_all,2)
  ax = axs(type_i);
  t_this = test_all{type_i};
  scatter(ax, t_this.actual, t_this.predicted, 10, 'filled');
  xlim(ax, [0 300]);
  ylim(ax, [0 300]);
%   axis(ax,'equal'); 
  xlabel(ax, 'Actual IFR');
  ylabel(ax, 'Predicted IFR');
  title(ax, t_str{type_i});
end







