% find example neurons that have good model performance

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
fd_save_plot = fullfile(fd_save, 'plotsExample');
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
neuID_all = {};
birID_all = {};
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
    neuID_all = [neuID_all; neuronID];
    birID_all = [birID_all; birdID];
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



%% 2. Calculate median of metrics over all folds
field_names = m.Compound_Full.Properties.VariableNames;
field_check = 'improvement_over_null';
field_i = find(strcmp(field_names, field_check));

% grab the metrics values
m_v = cellfun(@(x) x(:,field_i), res, 'UniformOutput', false);
% calculate the median metrics over folds
% m_mean = cellfun(@(x) mean(x), m_v);
m_mean = cellfun(@(x) median(x), m_v);



%% 3. Search example neurons for calls




%% 4. Plot predicted vs real
% ni=65;
ni=1;
% ni = 27;
birdID = birID_all{ni};
pairID = pairIDs{strcmp(birdIDs, birdID)};
% neuronID = '20240917-ch13';
neuronID = neuID_all{ni};
fd_res = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel', 'XGB_res', neuronID);
model_type = 'Call_Full';
m_this = readtable(fullfile(fd_res, model_type, sprintf('%s.metrics.csv', model_type)));
% what fold to plot
% fold_num = 3;
[ipv_score, fold_num] = max(m_this{:,4});
% load predicted vs actual
d_pred = readtable(fullfile(fd_res, model_type, sprintf('%s.pred_actual.fold%d.csv', model_type, fold_num-1)));
d_pred = table2array(d_pred);

% plot
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 500 500]);
col = '#0072BD';
scatter(d_pred(1,:), d_pred(2,:), 30, 'Marker', 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
xlim([-0.5 max(d_pred(1,:))+0.5]);
ylim([-0.5 max(d_pred(1,:))+0.5]);
xlabel('Observed no. of spikes', 'FontSize', 14);
ylabel('Predicted no. of spikes', 'FontSize', 14);
title(sprintf('%s ipv=%.3f', neuronID, ipv_score), 'FontSize', 14);
fn_fig = fullfile(fd_res, model_type, sprintf('pred_vs_actual.fold%d.pdf', fold_num));
print(fig, fn_fig, '-dpdf', '-painters');


%% test significance
m = load(fullfile(fd_res, sprintf('%s.%s.metrics_all.mat', birdID, neuronID)));
m = m.metrics_all;
m_full = m.Call_Full.improvement_over_null;
m_aco = m.Call_Acoustic.improvement_over_null;
m_time = m.Call_Time.improvement_over_null;
% perform tests
p1 = signrank(m_full, m_aco);
p2 = signrank(m_full, m_time);





