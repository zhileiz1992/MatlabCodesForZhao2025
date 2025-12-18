% Check if models trained on calls can predict firing in compound syllables
% only on acoustic models


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
fd_save_plot = fullfile(fd_save, 'plotGeneralization');
if ~exist(fd_save_plot, 'dir'); mkdir(fd_save_plot); end


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% loop through neurons
% ni = 44;
metrics_all = cell(size(info_neu,1), 1);
for ni=61:size(info_neu,1)
  neuronID = info_neu.neuronID{ni};
  disp(neuronID);
  fd_res = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel', 'XGB_res', neuronID);
  
  %% first assess: call -> compound model
%   ref_name = 'Call_Acoustic';
%   test_name = 'Compound_Acoustic';
  ref_name = 'Compound_Acoustic';
  test_name = 'Call_Acoustic';
  % grab the test dataset
  fn_test_predictor = fullfile(fd_res, test_name, sprintf('%s.predictor.csv', test_name));
  fn_test_output = fullfile(fd_res, test_name, sprintf('%s.output.csv', test_name));
  % what models have been trained, should be in .json files
  fns_model = dir(fullfile(fd_res, ref_name, '*.json'));
  if isempty(fns_model); continue; end
  % where to save results
%   fd_save_this = fullfile(fd_res, 'General_call_compound');
  fd_save_this = fullfile(fd_res, 'General_compound_call');
  if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end
  fn_py = '/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl/ensembleModel/ZZ_checkAccuracyOnTest_v1.py';
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
  metrics_all{ni} = metrics_this;
end

% concatenate into one array
metrics_a = cat(3, metrics_all{:});
% fn_m = fullfile(fd_save_plot, 'generalization_call_compound.metrics_all.mat'); 
 fn_m = fullfile(fd_save_plot, 'generalization_compound_call.metrics_all.mat'); 
save(fn_m, 'metrics_a');


%% plot results
% check the improvement over null
m_improv = squeeze(metrics_a(:, 5, :));
m_median = median(m_improv, 1);
close all; figure; 
histogram(m_median, 'BinEdges', -50:5:50);

% plot as violin plots
figure; 
violinplot(m_improv);
ylim([-50 50]);
xlim([0.5 size(m_improv,2)+0.5]);
ylabel('% improvements over null', 'FontSize', 14);
% title('Calls model -> Compound syls', 'FontSize', 14);
title('Compound syls -> Calls', 'FontSize', 14);

% only those that are signficantly better than 0
idx_sig = [];
for ni=1:size(metrics_all, 1)
  m = metrics_all{ni};
  if ~isempty(m)
    p = ZZfunc_runSignRank_v1(m(:,5));
    if p<=0.05
      idx_sig = [idx_sig ni];
    end
  end
end
m_sig = metrics_all(idx_sig);
m_sig_a = cat(3, m_sig{:});
m_sig_improv = squeeze(m_sig_a(:, 5, :));
figure; 
violinplot(m_sig_improv);
ylabel('% improvements over null', 'FontSize', 14);
% title('Calls model -> Compound syls', 'FontSize', 14);
title('Compound syls -> Calls', 'FontSize', 14);

% what are these neurons
neuron_sig = info_neu.neuronID(idx_sig);











