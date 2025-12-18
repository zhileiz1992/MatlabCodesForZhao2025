% check the generalization using GLM modeling
% acoustics only predictors

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

% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% where to save results
fd_save_base = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel');
fd_save_glm = fullfile(fd_save_base, 'glm');


% loop over neurons
% ni = 44;
r2_list = cell(size(info_neu,1), 1);
for ni=1:size(info_neu, 1)
  neuronID = info_neu.neuronID{ni};
  disp(neuronID);
  
  % load the curated data
  fd_data = fullfile(fd_save_base, 'CuratedData2');
  fn_d = fullfile(fd_data, sprintf('%s.%s.call_comp.modelData.mat', birdID, neuronID));
  a = load(fn_d);
  d_all = a.d_all;
  
  % concatenate into a large table
  d = vertcat(d_all{:});
  % remove those that have nan in the VAE latents
  d_m = d(~isnan(d.vae1),:);
  
  
  %% Build the GLM model: try on different inputs
  % try on different inputs
  d_comp = d_m(ismember(d_m.category, {'b', 'x'}),:);
  d_call = d_m(strcmp(d_m.category, 'v'),:);
  
  pred_vae = compose('vae%d', 1:32);
  output_name = 'spike_count';
  group_name = 'syl_ID';
  % dist = 'gamma';
  dist = 'poisson';
  link = 'log';
  train_test_ratio = 0.8;
  n = 10;  % how many times to repeat
  rng(1992);
  r2_all = nan(n, 6);  % comp_on_comp; comp_on_call;shuffleCall_on_call;  call_on_call; call_on_comp; huffleComp_on_comp;
  
  for ii=1:n
    % build a model on compounds
    [mdl1, metrics1, ~, ~] = glm_grouped_pipeline2(d_comp, pred_vae, output_name, group_name, dist, link, train_test_ratio, '');
    r2_all(ii, 1) = metrics1.R2_dev_test;
    % evaluate on calls
    X_te = d_call(:, pred_vae);
    y_te = d_call(:, output_name);
    r2_all(ii, 2) = ZZfunc_getGLMdeviance_v1(mdl1, X_te, y_te);
    
    % build a model on shuffled compounds
    [mdl1s, metrics1s, ~, ~] = glm_grouped_pipeline2_shuffle(d_comp, pred_vae, output_name, group_name, dist, link, train_test_ratio, '');
    r2_all(ii, 6) = metrics1s.R2_dev_test;
    
    
    % build a model on calls
    [mdl2, metrics2, ~, ~] = glm_grouped_pipeline2(d_call, pred_vae, output_name, group_name, dist, link, train_test_ratio, '');
    r2_all(ii, 4) = metrics2.R2_dev_test;
    % evaluate on calls
    X_te = d_comp(:, pred_vae);
    y_te = d_comp(:, output_name);
    r2_all(ii, 5) = ZZfunc_getGLMdeviance_v1(mdl2, X_te, y_te);
    
    % % build a model on shuffled calls
    [mdl2s, metrics2s, ~, ~] = glm_grouped_pipeline2_shuffle(d_call, pred_vae, output_name, group_name, dist, link, train_test_ratio, '');
    r2_all(ii, 3) = metrics2s.R2_dev_test;
  end
  % save results
  r2_list{ni} = r2_all;
  
  res.ni = ni;
  res.neuronID = neuronID;
  res.r2_all = r2_all;
  res.models = {mdl1, mdl1s, mdl2, mdl2s};
  fn_res = fullfile(fd_save_glm, sprintf('%s.%s.glm_res.mat', birdID, neuronID));
  save(fn_res, 'res');
  
  % plot results as violin plot
  close all;
  fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
  violinplot(r2_all);
  ylabel('Deviance from null model');
  xticks([1 2 3 4 5 6]);
  xticklabels({'comp_on_comp', 'comp_on_call', 'shuffle_call', 'call_on_call', 'call_on_comp', 'shuffle_comp'});
  ax = gca;  % get current axes
  ax.TickLabelInterpreter = 'none';
  title(neuronID);
  fn_fig = fullfile(fd_save_glm, sprintf('%s.%s.glm_res.pdf', birdID, neuronID));
  print(fig, fn_fig, '-dpdf', '-painters');
%   close(fig);
  
end


%% Plot results
r2_a = cat(3, r2_list{:});
figure; 
violinplot(squeeze(r2_a(:, 2, :)));
ylim([-0.5 0.5]);
yline(0);
ylabel('Deviance from null model', 'FontSize', 14);
title('Compound syl. model on calls', 'FontSize', 14);


figure; 
violinplot(squeeze(r2_a(:, 5, :)));
ylim([-0.5 0.5]);
yline(0);
ylabel('Deviance from null model', 'FontSize', 14);
title('Call model on compound syl.', 'FontSize', 14);










