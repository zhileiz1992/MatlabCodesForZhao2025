% pull and plot the results from XGBoost models

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


bi = 4;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird specific folder setting
% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% where to save results
fd_save_base = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel');
fd_save_plot = fullfile(fd_save_base, 'plots');
if ~exist(fd_save_plot, 'dir'); mkdir(fd_save_plot); end


%% 1. Load the modeling results
% what models to check 
model_names = {'Compound_Full', 'Compound_Acoustic', 'Compound_Time', ...
           'Call_Full', 'Call_Acoustic', 'Call_Time', ...
           'v1_Full', 'v1_Acoustic', 'v1_Time'};
input_names = {'Compound', 'Call', 'v1'};
type_names = {'Full', 'Acoustic', 'Time'};
res = cell(size(info_neu, 1), size(model_names, 2));
num_split_fold = 10;  % how many folds were used for modeling
nan_a = nan(num_split_fold, 5);  % a default metrics if no results
for ni=1:size(info_neu, 1)
  neuronID = info_neu.neuronID{ni};
  fd_res = fullfile(fd_save_base, 'XGB_res', neuronID);
  fn_metric = fullfile(fd_res, sprintf('%s.%s.metrics_all.mat', birdID, neuronID));
  a = load(fn_metric);
  m = a.metrics_all;
  for mi=1:size(model_names, 2)
    mn = model_names{mi};
    if ismember(mn, fieldnames(m))
      res{ni, mi} = table2array(m.(mn));
    else
      res{ni, mi} = nan_a;
    end
  end
end
  

%% 2. Plot the distribution of model performance
field_names = m.Compound_Full.Properties.VariableNames;
field_check = 'improvement_over_null';
field_i = find(strcmp(field_names, field_check));

% grab the metrics values
m_v = cellfun(@(x) x(:,field_i), res, 'UniformOutput', false);
% calculate the mean metrics over folds
m_mean = cellfun(@(x) mean(x), m_v);

% 2.1 plot the overall trend for all neurons
close all; 
[fig, axs] = generatePanelGrid_v2(1, 3, [0.75], [], [0.075;0.05], [0.075;0.05], 0.075, [0], [10 10 1200 500]);
hold(axs, 'on');
for ni=1:size(m_mean,1)
  % first on compound syllable
  plot(axs(1), [1 2 3], m_mean(ni, 1:3), 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 5, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
  % then on calls
  plot(axs(2), [1 2 3], m_mean(ni, 4:6), 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 5, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
  % then only v1
  plot(axs(3), [1 2 3], m_mean(ni, 7:9), 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 5, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
end
ylim(axs(1), [-10 50]); ylim(axs(2), [-10 75]); ylim(axs(3), [-10 85]);
for ai=1:3
  xlim(axs(ai), [0.8 3.2]);
  title(axs(ai), input_names{ai});
  xticks(axs(ai), [1 2 3]);
  xticklabels(axs(ai), type_names);
  axs(ai).FontSize = 14;
  ylabel(axs(ai), field_check);
end


% 2.2 Only those that have significant improvement over null for the full model
p_v = cellfun(@(x) ZZfunc_runSignRank_v1(x), m_v);
thre = 0.05;
close all; 
[fig, axs] = generatePanelGrid_v2(1, 3, [0.75], [], [0.075;0.05], [0.075;0.05], 0.075, [0], [10 10 1200 500]);
hold(axs, 'on');
pass_count = zeros(3,1);
for ni=1:size(m_mean,1)
  % first on compound syllable
  if ~isnan(p_v(ni, 1)) && p_v(ni,1)<=thre
    plot(axs(1), [1 2 3], m_mean(ni, 1:3), 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 5, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(1) = pass_count(1)+1;
  end
  % then on calls
  if ~isnan(p_v(ni, 4)) && p_v(ni,4)<=thre
    plot(axs(2), [1 2 3], m_mean(ni, 4:6), 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 5, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(2) = pass_count(2)+1;
  end
  % then only v1
  if ~isnan(p_v(ni, 7)) && p_v(ni,7)<=thre
    plot(axs(3), [1 2 3], m_mean(ni, 7:9), 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 5, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(3) = pass_count(3)+1;
  end
end
ylim(axs(1), [-5 50]); ylim(axs(2), [-10 70]); ylim(axs(3), [-5 85]);
for ai=1:3
  xlim(axs(ai), [0.8 3.2]);
  title(axs(ai), sprintf('%s n=%d', input_names{ai}, pass_count(ai)));
  xticks(axs(ai), [1 2 3]);
  xticklabels(axs(ai), type_names);
  axs(ai).FontSize = 14;
  ylabel(axs(ai), field_check, 'Interpreter', 'none');
end
fn_fig = fullfile(fd_save_plot, sprintf('%s.improvement_over_null.pdf', birdID));
print(fig, fn_fig, '-dpdf', '-painters');


% 2.3 Plot the contribution of Acoustics and Timing predictors
% defined as the drop in performance after removing the predictor
close all; 
[fig, axs] = generatePanelGrid_v2(1, 3, [0.75], [], [0.075;0.05], [0.1;0.05], 0.075, [0], [10 10 800 500]);
hold(axs, 'on');
pass_count = zeros(3,1);
contr1 = []; contr2 = []; contr3 = [];
for ni=1:size(m_mean,1)
  % first on compound syllable
  if ~isnan(p_v(ni, 1)) && p_v(ni,1)<=thre
    y = [m_mean(ni,1)-m_mean(ni,3), m_mean(ni,1)-m_mean(ni,2)];
    plot(axs(1), [1, 2], y, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(1) = pass_count(1)+1;
    contr1 = [contr1; y];
  end
  % then on calls
  if ~isnan(p_v(ni, 4)) && p_v(ni,4)<=thre
    y = [m_mean(ni,4)-m_mean(ni,6), m_mean(ni,4)-m_mean(ni,5)];
    plot(axs(2), [1, 2], y, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(2) = pass_count(2)+1;
    contr2 = [contr2; y];
  end
  % then only v1
  if ~isnan(p_v(ni, 7)) && p_v(ni,7)<=thre
    y = [m_mean(ni,7)-m_mean(ni,9), m_mean(ni,7)-m_mean(ni,8)];
    plot(axs(3), [1, 2], y, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(3) = pass_count(3)+1;
    contr3 = [contr3; y];
  end
end
for ai=1:3
  xlim(axs(ai), [0.8 2.2]);
  ylim(axs(ai), [-5 55]);
  title(axs(ai), sprintf('%s n=%d', input_names{ai}, pass_count(ai)));
  xticks(axs(ai), [1 2]);
  xticklabels(axs(ai), {'Rm. Acoustics', 'Rm. Time'});
  axs(ai).FontSize = 12;
  ylabel(axs(ai), 'Drop in model performance');
end
fn_fig = fullfile(fd_save_plot, sprintf('%s.contribution.pdf', birdID));
print(fig, fn_fig, '-dpdf', '-painters');

% perform simple stats
p1 = signrank(contr1(:,1), contr1(:,2));
p2 = signrank(contr2(:,1), contr2(:,2));
p3 = signrank(contr3(:,1), contr3(:,2));



%% 3. Replot model performance for selected example neurons
ni = 44; 
neuronID = info_neu.neuronID{ni};
fd_res = fullfile(fd_save_base, 'XGB_res', neuronID);
fn_metric = fullfile(fd_res, sprintf('%s.%s.metrics_all.mat', birdID, neuronID));
a = load(fn_metric);
metrics_all = a.metrics_all;

% plot the model performance as violin plots
m_plot = 'improvement_over_null';
close all; 
[fig, axs] = generatePanelGrid_v2(1, 3, [0.75], [], [0.075;0.05], [0.075;0.05], 0.075, [0], [10 10 800 500]);
% y_lim = [0 42]; 
for di=1:size(input_names,2)
  ax = axs(di); cla(ax);
  v_plot = [];
  grp_plot = [];
  v_stat = [];
  for mi=1:size(type_names, 2)
    run_name = sprintf('%s_%s', input_names{di}, type_names{mi});
    if ismember(run_name, fieldnames(metrics_all))
      v_this = metrics_all.(run_name).(m_plot);
      v_plot = [v_plot; v_this];
      v_stat = [v_stat v_this];
      grp_this = repmat(type_names(mi), size(v_this,1), 1);
      grp_plot = [grp_plot; grp_this];
    end
  end
  if isempty(v_plot); continue; end
%   boxplot(ax, v_plot);
  violinplot(v_plot, grp_plot, 'Parent', ax, 'GroupOrder', type_names, 'MedianMarkerSize', 50, 'MarkerSize', 30);
  xlim(ax, [0.5 3.5]);
%   ylim(ax, [0 max(v_plot)+2]);
  ylim(ax, [0 42]);
  ylabel(ax, '% improvement over null');
  title(ax, input_names{di});
  ax.FontSize =14;
end
fn_fig = fullfile(fd_save_plot, sprintf('%s.exampleSameY.%s.pdf', birdID, neuronID));
print(fig, fn_fig, '-dpdf', '-painters');

% perform stats
p112 = signrank(metrics_all.Compound_Full.(m_plot), metrics_all.Compound_Acoustic.(m_plot));
p113 = signrank(metrics_all.Compound_Full.(m_plot), metrics_all.Compound_Time.(m_plot));

p212 = signrank(metrics_all.Call_Full.(m_plot), metrics_all.Call_Acoustic.(m_plot));
p213 = signrank(metrics_all.Call_Full.(m_plot), metrics_all.Call_Time.(m_plot));









