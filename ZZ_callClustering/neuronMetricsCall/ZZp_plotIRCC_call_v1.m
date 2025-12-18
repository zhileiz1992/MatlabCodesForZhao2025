% script to plot IRCC values that are already calculated
% Zhilei, 06/30/2025
% a criteria metric struct array has been calculated when running the ZZp_plotMOSeq_call_v3 script
% unit of averaging is each neuron-syllable pair


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_callClustering'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};


%% 1. Load the metrics
d_all = cell(size(birdIDs,2), 1);
for bi=1:size(birdIDs, 2)
  bd = birdIDs{bi};
  fns_d = dir(fullfile(fd_home, 'Figures', pairIDs{bi}, 'HahnloserNew',  bd, 'v*', '*.criteria.mat'));
  d_this = [];
  % loop through call types
  for ci=1:size(fns_d,1)
    fn_d = fullfile(fns_d(ci).folder, fns_d(ci).name);
    load(fn_d);
    % only look at neurons that pass the criteria
    p = criteria([criteria.isPass]);
    d_this = [d_this p.ircc];
    fprintf('%s c%d %.3f\n', bd, ci, mean([p.ircc]));
  end
  d_all{bi} = d_this;
end


%% 2. Plot distribution
% GPT
% Convert d_all into a single numeric vector and group labels
all_data = [];
group_labels = [];

for k = 1:numel(d_all)
    this_data = d_all{k}(:); % ensure column vector
    fprintf('%s mean IRCC: %.3f\n', pretty_ids{k}, mean(this_data));
    all_data = [all_data; this_data];
    group_labels = [group_labels; repmat(k, numel(this_data), 1)];
end
fprintf('Overall mean IRCC: %.3f\n', mean(all_data));

% Create plot
close all; figure;
boxplot(all_data, group_labels);

xlabel('Bird ID');
ylabel('Value');
title('Box Plot of Distributions in d\_all');

% save figure
fn_fig = fullfile(fd_home, 'Figures', 'CombinedAnalysis', 'IRCC.boxplot.fig');
savefig(fn_fig);




