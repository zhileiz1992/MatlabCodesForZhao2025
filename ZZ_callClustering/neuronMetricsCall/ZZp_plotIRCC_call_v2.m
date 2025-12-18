% script to plot IRCC values that are already calculated
% Zhilei, 06/30/2025
% a criteria metric struct array has been calculated when running the ZZp_plotMOSeq_call_v3 script
% differ from v1: unit of averaging is each neuron

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
  d_this = struct([]);
  % loop through call types
  for ci=1:size(fns_d,1)
    fn_d = fullfile(fns_d(ci).folder, fns_d(ci).name);
    load(fn_d);
    syl = strsplit(fns_d(ci).name, '-');
    % only look at neurons that pass the criteria
    p = criteria([criteria.isPass]);
    [p.syl] = deal(syl{2});
    d_this = [d_this p];
  end
  
  % calculate the mean ircc of each neuron
  T = struct2table(d_this);
  % Compute mean ircc per neuronID
  result = groupsummary(T, 'neuronID', 'mean', 'ircc');
  d_all{bi} = result.mean_ircc;
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
close all; 
fig = ZZfunc_newFigurePDFsize_v1([50 50 300 400]); 
% boxplot(all_data, group_labels);
violinplot(all_data, group_labels, 'ViolinColor', [0.21 0.49 0.72], 'EdgeColor', [1 1 1], 'MarkerSize', 7, 'MedianMarkerSize', 100);
% add a sample size number to each violin 
for bi=1:size(d_all,1)
  text(bi, 0, sprintf('%d', length(d_all{bi})), 'FontSize', 8, 'HorizontalAlignment', 'center');
end
ylim([-1 1]);
xlim([0.5 4.5]);
xticks([1 2 3 4]);
xticklabels(pretty_ids);
xlabel('Bird ID');
ylabel('Inter-rendition correlation coeffients');
title('IRCC');

% save figure
fn_fig = fullfile(fd_home, 'Figures', 'CombinedAnalysis', 'IRCC.boxplot.fig');
savefig(fn_fig);
fn_pdf = strrep(fn_fig, '.fig', '.pdf');
print(fig, fn_pdf, '-dpdf', '-painters');


%% 3. Check why some neurons are not included
bi = 2;
fn_meta = fullfile(fd_home, 'DbaseFiles', pairIDs{bi}, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdIDs{bi}));
load(fn_meta);
bd = birdIDs{bi};
fns_d = dir(fullfile(fd_home, 'Figures', pairIDs{bi}, 'HahnloserNew',  bd, 'v*', '*.criteria.mat'));
d_this = struct([]);
% loop through call types
for ci=1:size(fns_d,1)
  fn_d = fullfile(fns_d(ci).folder, fns_d(ci).name);
  load(fn_d);
  syl = strsplit(fns_d(ci).name, '-');
  % only look at neurons that pass the criteria
  p = criteria([criteria.isPass]);
  [p.syl] = deal(syl{2});
  d_this = [d_this p];
end
% get unique neuron ID
neu_in = unique({d_this.neuronID}, 'sorted');
% check what neurons are not included
neu_all = info.neuronID;
neu_out = setdiff(neu_all, neu_in);
















