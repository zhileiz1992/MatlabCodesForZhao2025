% Check if models trained on calls/compounds can predict firing in compound syllables / calls
% only on acoustic models
% differ from v1: loop over all birds


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
fd_save_plot = fullfile(fd_save, 'plotGeneralization2');
if ~exist(fd_save_plot, 'dir'); mkdir(fd_save_plot); end


% bi = 1;
for bi=2:size(birdIDs,2)
  
birdID = birdIDs{bi};
pairID = pairIDs{bi};


% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% loop through neurons
% ni = 44;
metrics_all = cell(size(info_neu,1), 2);
ref_list = {'Compound_Acoustic', 'Call_Acoustic'};
test_list = {'Call_Acoustic', 'Compound_Acoustic'};
fn_py = '/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl/ensembleModel/ZZ_checkAccuracyOnTest_v1.py';

for ni=1:size(info_neu,1)
  neuronID = info_neu.neuronID{ni};
  fprintf('%d %d %s\n', bi, ni, neuronID);
  fd_res = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel', 'XGB_res', neuronID);
  
  for di=1:size(ref_list, 2)
    ref_name = ref_list{di};
    test_name = test_list{di};
    % grab the test dataset
    fn_test_predictor = fullfile(fd_res, test_name, sprintf('%s.predictor.csv', test_name));
    fn_test_output = fullfile(fd_res, test_name, sprintf('%s.output.csv', test_name));
    % what models have been trained, should be in .json files
    fns_model = dir(fullfile(fd_res, ref_name, '*.json'));
    if isempty(fns_model); continue; end
    
    fd_save_this = fullfile(fd_res, sprintf('Generalize.%s.%s', ref_name, test_name));
    if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end
    metrics_this = [];
    for mi=1:size(fns_model, 1)
      model_path = fullfile(fns_model(mi).folder, fns_model(mi).name);
      fn_save = fullfile(fd_save_this, sprintf('model.fold%d.csv', mi));
      fn_temp_py = fullfile(fd_save_this, 'temp2.py');
      metrics = ZZfunc_checkXGBoostOnTest_v1(model_path, fn_test_predictor, fn_test_output, fn_save, fn_py, fn_temp_py);
      if ~isempty(metrics)
        metrics = table2array(metrics);
        metrics_this = [metrics_this; metrics'];
      end
    end
    metrics_all{ni, di} = metrics_this;
  end
end

% concatenate into one array
% metrics_a = cat(3, metrics_all{:});
% fn_m = fullfile(fd_save_plot, 'generalization_call_compound.metrics_all.mat');
fn_m = fullfile(fd_save_plot, sprintf('%s.generalization.metrics_all.mat', birdID));
save(fn_m, 'metrics_all');

end



%% Combine results into one data table
% read the cross-modeling metrics and within-modeling metrics
m_full = []; 
for bi=1:size(birdIDs,2)
  % load information regarding sparse neurons
  fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairIDs{bi}, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdIDs{bi}));
  a = load(fn_info_neu); info_neu = a.info;
  % load cross-modeling metrics
  fn_m = fullfile(fd_save_plot, sprintf('%s.generalization.metrics_all.mat', birdIDs{bi}));
  b = load(fn_m); m_cross = b.metrics_all;
  % then load the within-modeling results
  % only look at the improvements over null
  m = info_neu; 
  fd_within = fullfile(fd_base, 'Figures', pairIDs{bi}, 'CompoundSyl', birdIDs{bi}, 'ensembleModel', 'XGB_res');
  for ni=1:size(info_neu,1)
    if ~isempty(m_cross{ni,1}); m.comp_on_call{ni} = m_cross{ni,1}(:, 5); end
    if ~isempty(m_cross{ni,2}); m.call_on_comp{ni} = m_cross{ni,2}(:, 5); end
    fn_m = fullfile(fd_within, info_neu.neuronID{ni}, sprintf('%s.%s.metrics_all.mat', birdIDs{bi}, info_neu.neuronID{ni}));
    c = load(fn_m); m_within = c.metrics_all;
    if ismember('Compound_Acoustic', fieldnames(m_within)); m.comp_on_comp{ni} = m_within.Compound_Acoustic.improvement_over_null; end
    if ismember('Call_Acoustic', fieldnames(m_within)); m.call_on_call{ni} = m_within.Call_Acoustic.improvement_over_null; end
  end
  m.birdID = repmat(birdIDs(bi), size(m,1), 1);
  m_full = [m_full; m];
end

% add columns to indicate if models perform sigificantly better than null
field_n = {'comp_on_call', 'call_on_comp',  'comp_on_comp', 'call_on_call'};
for fi=1:size(field_n,2)
  temp = cellfun(@(x) ZZfunc_runSignRank_v1(x), m_full.(field_n{fi}));
  m_full.(['p_' field_n{fi}]) = temp;
  temp2 = cellfun(@(x) max([0 median(x)]), m_full.(field_n{fi}));  % set to 0 is smaller than 0
  m_full.(['median_' field_n{fi}]) = temp2;
end

% save the data
fn_full = fullfile(fd_save_plot, 'combined.m_full.mat');
save(fn_full, 'm_full');


%% plot results
p_thre = 0.05; 
v_thre = 0; 
close all; 
[fig, axs] = generatePanelGrid_v2(1, 2, [0.8], [], [0.075;0.05], [0.075;0.05], 0.075, [0], [10 10 800 400]);
hold(axs, 'on');
% first plot comp_on_comp vs comp_on_call, only plot neurons with comp_on_comp been significantly better
ax1 = axs(1);
i_sig = find((m_full.p_comp_on_comp<=p_thre) & (m_full.median_comp_on_comp>=v_thre));
x = m_full.median_comp_on_comp(i_sig); 
y = m_full.median_comp_on_call(i_sig);
p = m_full.p_comp_on_call(i_sig);
mask = (p<=p_thre);
scatter(ax1, x(~mask), y(~mask), 30, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none');
scatter(ax1, x(mask), y(mask), 30, 'Marker', 'o', 'MarkerFaceColor', '#e41a1c', 'MarkerEdgeColor', 'none');
ylim(ax1, [0 max(y)+2]);
xlim(ax1, [0 max(x)+2]);
xlabel(ax1, '% improvement on compound syls.', 'FontSize', 12);
ylabel(ax1, '% improvement on calls', 'FontSize', 12);
title(ax1, sprintf('Compound model -> calls (%d/%d)', sum(mask), length(i_sig)), 'FontSize', 12);

% then plot call_on_call vs call_on_comp, only plot neurons with call_on_call been significantly better
ax2 = axs(2);
i_sig = find((m_full.p_call_on_call<=p_thre) & (m_full.median_call_on_call>=v_thre));
x = m_full.median_call_on_call(i_sig); 
y = m_full.median_call_on_comp(i_sig);
p = m_full.p_call_on_comp(i_sig);
mask = (p<=p_thre);
scatter(ax2, x(~mask), y(~mask), 30, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none');
scatter(ax2, x(mask), y(mask), 30, 'Marker', 'o', 'MarkerFaceColor', '#e41a1c', 'MarkerEdgeColor', 'none');
ylim(ax2, [0 max(y)+2]);
xlim(ax2, [0 max(x)+2]);
xlabel(ax2, '% improvement on calls', 'FontSize', 12);
ylabel(ax2, '% improvement on compound syls.', 'FontSize', 12);
title(ax2, sprintf('Call model -> compound syls. (%d/%d)', sum(mask), length(i_sig)), 'FontSize', 12);

fn_fig = fullfile(fd_save_plot, 'combined.cross_model.pdf');
print(fig, fn_fig, '-dpdf', '-painters');










