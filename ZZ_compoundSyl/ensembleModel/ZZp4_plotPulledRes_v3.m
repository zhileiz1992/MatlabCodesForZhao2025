% pull and plot the results from XGBoost models
% differ from v2: use median, instead of mean as the averaged measure

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

% save results to bird M1
fd_save = fullfile(fd_base, 'Figures', pairIDs{1}, 'CompoundSyl', birdIDs{1}, 'ensembleModel');
fd_save_plot = fullfile(fd_save, 'plotsCombined2');
if ~exist(fd_save_plot, 'dir'); mkdir(fd_save_plot); end



%% Load results from all birds
model_names = {'Compound_Full', 'Compound_Acoustic', 'Compound_Time', ...
           'Call_Full', 'Call_Acoustic', 'Call_Time', ...
           'v1_Full', 'v1_Acoustic', 'v1_Time'};
input_names = {'Compound', 'Call', 'v1'};
type_names = {'Full', 'Acoustic', 'Time'};
num_split_fold = 10;  % how many folds were used for modeling
nan_a = nan(num_split_fold, 5);  % a default metrics if no results

res_all = cell(size(birdIDs, 2), 1);
% loop through birds
for bi=1:size(birdIDs,2)
  birdID = birdIDs{bi};
  pairID = pairIDs{bi};
  
  % load information regarding sparse neurons
  fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
  a = load(fn_info_neu); info_neu = a.info;
  
  % Load the modeling results
  res = cell(size(info_neu, 1), size(model_names, 2));
  for ni=1:size(info_neu, 1)
    neuronID = info_neu.neuronID{ni};
    fd_res = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel', 'XGB_res', neuronID);
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
  
  res_all{bi} = res;
end
% concatenate into a master data
res = vertcat(res_all{:});



%% 2. Plot the distribution of model performance
field_names = m.Compound_Full.Properties.VariableNames;
field_check = 'improvement_over_null';
field_i = find(strcmp(field_names, field_check));

% grab the metrics values
m_v = cellfun(@(x) x(:,field_i), res, 'UniformOutput', false);
% calculate the median metrics over folds
% m_mean = cellfun(@(x) mean(x), m_v);
m_mean = cellfun(@(x) median(x), m_v);


%% 2.1 plot the overall trend for all neurons
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


%% 2.2 Only those that have significant improvement over null for the full model
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



%% 2.3 Plot the contribution of Acoustics and Timing predictors
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



%% 2.4 Plot the ratio instead
close all; 
[fig, axs] = generatePanelGrid_v2(1, 3, [0.75], [], [0.075;0.05], [0.1;0.05], 0.075, [0], [10 10 800 500]);
hold(axs, 'on');
pass_count = zeros(3,1);
contr1 = []; contr2 = []; contr3 = [];
% set threshold on full model, only check those that pass thresholds
p_thre = 0.05;  % improvements significantly larger than 0
full_thre = 15;  % mean improvements larger than 15%
temp2 = [];
for ni=1:size(m_mean,1)
  % first on compound syllable
  if ~isnan(p_v(ni, 1)) && p_v(ni,1)<=p_thre && m_mean(ni,1)>=full_thre
    % if improvement lower than zero, set to zero
    y_base = m_mean(ni,1); 
    y1 = max([0 m_mean(ni,2)]);
    y2 = max([0 m_mean(ni,3)]);
    y = [y1/y_base y2/y_base];
    plot(axs(1), [1, 2], y, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(1) = pass_count(1)+1;
    contr1 = [contr1; y];
  end
  % then on calls
  if ~isnan(p_v(ni, 4)) && p_v(ni,4)<=p_thre && m_mean(ni,4)>=full_thre    
    y_base = m_mean(ni,4); 
    y1 = max([0 m_mean(ni,5)]);
    y2 = max([0 m_mean(ni,6)]);
    y = [y1/y_base y2/y_base];
    plot(axs(2), [1, 2], y, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(2) = pass_count(2)+1;
    contr2 = [contr2; y];
    temp2 = [temp2; ni];
  end
  % then only v1
  if ~isnan(p_v(ni, 7)) && p_v(ni,7)<=p_thre && m_mean(ni,7)>=full_thre
    y_base = m_mean(ni,7); 
    y1 = max([0 m_mean(ni,8)]);
    y2 = max([0 m_mean(ni,9)]);
    y = [y1/y_base y2/y_base];
    plot(axs(3), [1, 2], y, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373')
    pass_count(3) = pass_count(3)+1;
    contr3 = [contr3; y];
  end
end
for ai=1:3
  xlim(axs(ai), [0.8 2.2]);
  ylim(axs(ai), [-0.05 1.2]);
  title(axs(ai), sprintf('%s n=%d', input_names{ai}, pass_count(ai)));
  xticks(axs(ai), [1 2]);
  xticklabels(axs(ai), {'Rm. Acoustics', 'Rm. Time'});
  axs(ai).FontSize = 12;
  ylabel(axs(ai), 'Drop in model performance');
end
fn_fig = fullfile(fd_save_plot, sprintf('%s.contributionRatio.pdf', birdID));
print(fig, fn_fig, '-dpdf', '-painters');

% perform simple stats
p1 = signrank(contr1(:,1), contr1(:,2));
p2 = signrank(contr2(:,1), contr2(:,2));
p3 = signrank(contr3(:,1), contr3(:,2));



%% 2.5 Plot the performance of the full model as violin plot
% only plot the significant ones
p_thre = 0.05;  % improvements significantly larger than 0
full_thre = 15;  
% for calls
i_p_call = find((p_v(:,4)<=p_thre) & (m_mean(:,4)>=full_thre));
% then for compounds
i_p_comp = find((p_v(:,1)<=p_thre) & (m_mean(:,1)>=full_thre));
% concatenate
v_conc = [m_mean(i_p_call,4); m_mean(i_p_comp,1)];
grp_conc = [repmat({'Calls'},length(i_p_call),1); repmat({'Compound'},length(i_p_comp),1)];

% plot
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 250 500]);
violinplot(v_conc, grp_conc, 'ViolinColor', [0 0.45 0.74], 'MedianMarkerSize', 400, 'MarkerSize', 30);
ylim([0 90]);
xlim([0.5 2.5]);
ylabel('% improvement over null', 'FontSize', 14);
fn_fig = fullfile(fd_save_plot, sprintf('%s.fullModelPerformance.pdf', birdID));
print(fig, fn_fig, '-dpdf', '-painters');

p = ranksum(m_mean(i_p_call,4), m_mean(i_p_comp,1));









