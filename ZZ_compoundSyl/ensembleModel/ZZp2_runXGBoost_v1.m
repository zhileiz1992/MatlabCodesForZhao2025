% call python to run XGBoost

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



%% Build the XGBoost model: try on different inputs
% try on different inputs
d_comp = d_m(ismember(d_m.category, {'b', 'x'}),:);
d_call = d_m(strcmp(d_m.category, 'v'),:);
d_v1 = d_m(strcmp(d_m.label, 'v1'),:);

% what predictors to use
pred_time = {'t_onset', 't_offset', 'rel_onset', 'rel_offset'};
pred_vae = compose('vae%d', 1:32);

% the same output and group for all models
% output_name = 'ifr';
output_name = 'spike_count';
group_name = 'syl_ID';
n_fold_split = '10'; 
fn_py = '/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl/ensembleModel/ZZ_runXGBoost_v1.py';

% loop through different inputs
input_names = {'Compound', 'Call', 'v1'};
d_list = {d_comp, d_call, d_v1};
% then different types of models
type_names = {'Full', 'Acoustic', 'Time'};
pred_list = {[pred_vae pred_time], pred_vae, pred_time};


metrics_all = struct();
for di=1:size(d_list, 2)
  for mi=1:size(type_names, 2)
    run_name = sprintf('%s_%s', input_names{di}, type_names{mi});
    fprintf('Building model for %s\n', run_name);
    d_input = d_list{di};
    pred_name = pred_list{mi};

    fd_save_res = fullfile(fd_save_base, 'XGB_res', neuronID, run_name);
    % fit the model
    metrics = ZZfunc_runXGBoostPython_v1(d_input, pred_name, output_name, group_name, n_fold_split, fd_save_res, fn_py, run_name);
    metrics_all.(run_name) = metrics;
  end
end
% save the master metrics
fn_m = fullfile(fd_save_base, 'XGB_res', neuronID, sprintf('%s.%s.metrics_all.mat', birdID, neuronID));
save(fn_m, 'metrics_all');


%% Plot results
close all; 
[fig, axs] = generatePanelGrid_v2(1, 3, [0.75], [], [0.075;0.05], [0.075;0.05], 0.075, [0], [10 10 1500 500]);
m_plot = 'improvement_over_null';
% m_plot = 'r2s_model';
col_types = {'#66c2a5','#fc8d62','#8da0cb'};
for di=1:size(d_list, 2)
  ax = axs(di); cla(ax);
  v_plot = [];
  grp_plot = [];
  v_stat = [];
  for mi=1:size(type_names, 2)
    run_name = sprintf('%s_%s', input_names{di}, type_names{mi});
    v_this = metrics_all.(run_name).(m_plot);
    v_plot = [v_plot; v_this];
    v_stat = [v_stat v_this];
    grp_this = repmat(type_names(mi), size(v_this,1), 1);
    grp_plot = [grp_plot; grp_this];
  end
%   boxplot(ax, v_plot);
  violinplot(v_plot, grp_plot, 'Parent', ax, 'GroupOrder', type_names, 'MedianMarkerSize', 50, 'MarkerSize', 30);
  xlim(ax, [0.5 3.5]);
  ylabel(ax, '% improvement over null');
  title(ax, input_names{di});
  ax.FontSize =14;
  
  % perform statistic tests
%   [p12,~,~] = signrank(v_stat(:,1), v_stat(:,2));
%   [p13, ~,~] = signrank(v_stat(:,1), v_stat(:,3));
%   % add a line to show the p-value
%   y_loc = max(v_plot) * 1.1;
%   plot(ax, [1 2], [y_loc y_loc], 'LineStyle', '-', 'LineWidth', 1', 'Color', '#737373');
%   text(ax, 1.5, y_loc*1.1, sprintf('%.4f', p12), 'FontSize', 12);
%   plot(ax, [1 3], [y_loc*1.2 y_loc*1.2], 'LineStyle', '-', 'LineWidth', 1', 'Color', '#737373');
%   text(ax, 2, y_loc*1.25, sprintf('%.4f', p13), 'FontSize', 12);
  
end
fd_save_plot = fullfile(fd_save_base, 'XGB_res', neuronID);
fn_fig = fullfile(fd_save_plot, sprintf('%s.%s.%s.pdf', birdID, neuronID, m_plot));
print(fig, fn_fig, '-dpdf', '-painters');













