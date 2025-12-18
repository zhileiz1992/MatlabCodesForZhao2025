% analyze the burst converage for shared and syllable-specific MO neurons for call subtype pairs
% Zhilei, 07/12/2025
% differ from v3: do it two pass, first strict, then relaxed for checking specific


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what call subtype pairs to analyze
to_analyze = {{{'v4', 'v5'}, {'v1', 'v7'}}, {}, {}, {}};

% loop through birds
bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
fd_save_base = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew', birdID);

% loop through analyze pairs
pi = 2;
p = to_analyze{bi}{pi};
% save plots in subfolder
fd_save_this = fullfile(fd_save_base, 'Coverage', [p{1} '-' p{2}]);
if ~exist(fd_save_this)
  mkdir(fd_save_this);
end
% grab the criteria data struct for the call pairs
fn1 = fullfile(fd_save_base, p{1}, sprintf('Hahnloser-%s-chan0.criteria5.mat', p{1}));
c1 = load(fn1); c1 = c1.criteria;
fn2 = fullfile(fd_save_base, p{2}, sprintf('Hahnloser-%s-chan0.criteria5.mat', p{2}));
c2 = load(fn2); c2 = c2.criteria;
% also grab the original spike data
fn1 = fullfile(fd_save_base, 'extractedPull', sprintf('%s.%s.segments_all.pull.mat', birdID, p{1}));
d1 = load(fn1); d1 = d1.seg_selected;
fn2 = fullfile(fd_save_base, 'extractedPull', sprintf('%s.%s.segments_all.pull.mat', birdID, p{2}));
d2 = load(fn2); d2 = d2.seg_selected;


%% determine if a neuron is shared or specific
% if missing data or too few trials, add NaN and exclude from comparison
all_neuron = union({c1.neuronID}, {c2.neuronID});
t_all = struct();
% loop through each neuron
for ni=1:size(all_neuron, 2)
  neuronID = all_neuron{ni};
  % get the already calculated criteria info, note that not all neurons are in both criteria data struct
  cn1 = c1(strcmp({c1.neuronID}, neuronID));
  cn2 = c2(strcmp({c2.neuronID}, neuronID));
  
  % save basic info to the data struct
  t_all(ni).neuronID = neuronID;
  t_all(ni).syl1 = p{1};
  t_all(ni).syl2 = p{2};
  t_all(ni).criteria1 = cn1;
  t_all(ni).criteria2 = cn2;
  
  % if both syllables have enough trials
  min_trial = 10;
  seg_this1 = d1(strcmp({d1.neuronID}, neuronID));
  seg_this2 = d2(strcmp({d2.neuronID}, neuronID));
  % default to not pass
  t_all(ni).isPass1 = 0;
  t_all(ni).isPass2 = 0;
  t_all(ni).total_rend1 = size(seg_this1, 2); 
  t_all(ni).total_rend2 = size(seg_this2, 2);
  if (size(seg_this1, 2)>=min_trial) && (size(seg_this2, 2)>=min_trial)
    % calculate sparseness
    psth_bin_size = 0.01;  % unit is seconds
    pad = 0.05;  % how much to look before and after the syllable, unit is seconds
    [spa1, aligned_spike1] = ZZfunc_calcSparseness_v1(seg_this1, pad, psth_bin_size);
    t_all(ni).sparse1 = spa1;
    [spa2, aligned_spike2] = ZZfunc_calcSparseness_v1(seg_this2, pad, psth_bin_size);
    t_all(ni).sparse2 = spa2;
    t_all(ni).t1 = cn1.psth_max_smooth_t;
    t_all(ni).t2 = cn2.psth_max_smooth_t;
    % check if the neuron passed using a strict criteria
    prop_thre = 0.25;  % include neuron if more than this proportion of renditions have spikes around the peak PSTH, set 0.25
    psth_thre = 10; 
    ircc_thre = 0.1; 
    min_rend = 10;  % minimal number of renditions that have spikes around the peak PSTH
    sparse_thre = 0.2; 
    p1 = (cn1.prop_psth_fire>=prop_thre) & (cn1.peak_psth>=psth_thre) & (cn1.ircc>=ircc_thre) ...
      & (cn1.num_rend_psth_fire>=min_rend) & (spa1>=sparse_thre);
    p2 = (cn2.prop_psth_fire>=prop_thre) & (cn2.peak_psth>=psth_thre) & (cn2.ircc>=ircc_thre) ...
      & (cn2.num_rend_psth_fire>=min_rend) & (spa2>=sparse_thre);
    
    % if both pass strict criteria
    if p1 && p2
      t_all(ni).isPass1 = 1;
      t_all(ni).isPass2 = 1;
    end
    
    % if only one pass, use a relaxed criteria to check for the other
    prop_thre = 0;  % include neuron if more than this proportion of renditions have spikes around the peak PSTH, set 0.25
    psth_thre = 5; 
    ircc_thre = 0.1; 
    min_rend = 4;  % minimal number of renditions that have spikes around the peak PSTH
    sparse_thre = 0.1; 
    if p1 && (~p2)
      t_all(ni).isPass1 = 1;
      t_all(ni).isPass2 = (cn2.prop_psth_fire>=prop_thre) & (cn2.peak_psth>=psth_thre) & (cn2.ircc>=ircc_thre) ...
      & (cn2.num_rend_psth_fire>=min_rend) & (spa2>=sparse_thre);
    end
    if (~p1) && p2
      t_all(ni).isPass1 = (cn1.prop_psth_fire>=prop_thre) & (cn1.peak_psth>=psth_thre) & (cn1.ircc>=ircc_thre) ...
      & (cn1.num_rend_psth_fire>=min_rend) & (spa1>=sparse_thre);
      t_all(ni).isPass2 = 1;
    end
  end
end


%% find the shared and specific neurons
i_share = find([t_all.isPass1] & [t_all.isPass2]);
i_specific1 = find([t_all.isPass1] & (~[t_all.isPass2]));
i_specific2 = find((~[t_all.isPass1]) & [t_all.isPass2]);
{t_all(i_specific1).neuronID}
{t_all(i_specific2).neuronID}

% plot where these bursts are
t_shared1 = [t_all(i_share).t1];
t_shared2 = [t_all(i_share).t2];
t_specific1 = [t_all(i_specific1).t1];
t_specific2 = [t_all(i_specific2).t2]; 
close all; 
fig = ZZfunc_newFigurePDFsize_v1([50 50 600 200]);
% plot shared as blue, distinguish birds
jit = 0.05;
plot(t_shared1, randn(size(t_shared1))*jit+1, 'o',  'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
hold on;
plot(t_shared2, randn(size(t_shared2))*jit+2, 'o',  'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
hold on;
% for specific, plot red dots
plot(t_specific1, randn(size(t_specific1))*jit+1.25, 'o',  'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
hold on;
plot(t_specific2, randn(size(t_specific2))*jit+2.25, 'o',  'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
hold on;

ylim([0.5 2.5]);
yticks([1 2]);
yticklabels(p);
ax = gca;
ax.YAxis.FontSize = 12;
temp = [t_shared1 t_shared2 t_specific1 t_specific2];
xlim([min(temp)-0.01 max(temp)+0.01]);
xlabel('Relative time (sec)', 'FontSize', 12);
title(sprintf('%s %s-%s', birdID, p{1}, p{2}), 'FontSize', 12);
% save figure
fn_pdf = fullfile(fd_save_this, sprintf('%s.%s-%s.tDist.pdf', birdID, p{1}, p{2}));
print(fig, fn_pdf, '-dpdf', '-painters');
% also save the t_all 
fn_t = fullfile(fd_save_this, sprintf('%s.%s-%s.t_all.mat', birdID, p{1}, p{2}));
save(fn_t, 't_all');


%% calculate a simple coverage metrics
% Compute histogram counts
t_shared = ([t_all(i_share).t1]+[t_all(i_share).t2])/2;  % take the average for shared
t_specific = [t_all(i_specific1).t1 t_all(i_specific2).t2]; 
% divide the time to 10-ms bins
tmin = min([t_shared t_specific]);
tmax = max([t_shared t_specific]);
bin_width = 0.01;
% Define bin edges
edges = tmin:bin_width:tmax;
% Optionally, get bin centers for plotting
bin_centers = edges(1:end-1) + bin_width/2;
% Compute histogram counts
count_shared = histcounts(t_shared, edges);
count_specific = histcounts(t_specific, edges);
% convert to covered or not
x = count_shared>0;
y = count_specific>0;
% smooth data
x_smooth = smoothdata(x, 'gaussian', 3);
y_smooth = smoothdata(y, 'gaussian', 3);

% plot coverage
figure;
plot(bin_centers, x_smooth);
hold on; 
plot(bin_centers, y_smooth);
  
% test for correlation
[rho, pval] = corr(x_smooth', y_smooth', 'Type', 'Pearson');
fprintf('Pearson r = %.4f\n', rho);
fprintf('p-value = %.4f\n', pval);

  
%% calculate separately by these two call subtypes, contanate the data
% Compute histogram counts
temp = max([t_shared1 t_specific1]);
t_shared = [t_shared1 t_shared2+temp+0.05];
t_specific = [t_specific1 t_specific2+temp+0.05];

% divide the time to 10-ms bins
tmin = min([t_shared t_specific]);
tmax = max([t_shared t_specific]);
bin_width = 0.01;
% Define bin edges
edges = tmin:bin_width:tmax;
% Optionally, get bin centers for plotting
bin_centers = edges(1:end-1) + bin_width/2;
% Compute histogram counts
count_shared = histcounts(t_shared, edges);
count_specific = histcounts(t_specific, edges);
% convert to covered or not
x = count_shared>0;
y = count_specific>0;
% smooth data
x_smooth = smoothdata(x, 'gaussian', 3);
y_smooth = smoothdata(y, 'gaussian', 3);

% plot coverage
figure;
plot(bin_centers, x_smooth);
hold on; 
plot(bin_centers, y_smooth);
  
% test for correlation
[rho, pval] = corr(x_smooth', y_smooth', 'Type', 'Pearson');
fprintf('Pearson r = %.4f\n', rho);
fprintf('p-value = %.4f\n', pval);
  
  
  
  
  
  
