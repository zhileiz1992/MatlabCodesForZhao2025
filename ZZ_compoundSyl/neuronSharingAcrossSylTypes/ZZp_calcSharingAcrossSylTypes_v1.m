% calculate a coarse metrics of how many syllable types a neuron can be tuned to

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));



%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};

% where to save results
fd_save = fullfile(fd_base, 'Figures', 'CombinedAnalysis', 'sharedFiring');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end

info_all = [];
syls = {'v', 'b', 'h', 'e', 'x', 'comp'}; 
syls_pretty = {'Call', 'Chortle', 'Click', 'Squawk', 'Other', 'Comp'};
for bi=1:size(birdIDs,2)
% bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
fprintf('Calculating for %s...\n', birdID);

% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;
% load data on spikes
vae_run = 'traj_chop_32_1_32';
fd_field = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_firingField']);
fn_spike = fullfile(fd_field, sprintf('%s.spike_win.mat', birdID));
a = load(fn_spike); spike_win = a.spike_win;



%% Loop through neurons, count spikes for each syllable category
% set thresholds on fraction of renditions that have spikes
p_thre = 0.25; 
% also thresholds on mean number of spikes
m_thre = 0.5;

info_count = info_neu;
fs = 20000;
pad_pre = 0.05; pad_pre_pt = floor(pad_pre*fs); % how much to consider before syllable onset
pad_post = 0;  pad_post_pt = floor(pad_post*fs);% how much after syllable offset
% ni = 44; 
for ni=1:size(info_neu,1)
  neuronID = info_neu.neuronID{ni};
  % subset the spike_win table
  spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
  
  % loop through syllable types
  s_count = [];
  for si=1:size(syls,2)
    v = syls{si};
    if strcmp(v, 'comp')
      spike = spike_this(ismember(spike_this.category, {'b', 'x'}),:);
      spike = spike(spike.dur>=0.3, :);
    else
      spike = spike_this(strcmp(spike_this.category, v),:);
    end
    s_count(si).syl = v;
    s_count(si).num_syl = size(spike,1);
    % how many syllable renditions have at least one spike
    s_c = zeros(size(spike,1), 1);
    for ri=1:size(spike,1)
      ori = spike.seg_start_ori(ri) - spike.seg_start(ri);
      rel_end = spike.seg_end_ori(ri) - spike.seg_start(ri);
      total_len = spike.seg_end(ri) - spike.seg_start(ri);
      % considering pad
      w_left = max([1 ori-pad_pre_pt]);
      w_right = min([total_len rel_end+pad_post_pt]);
      % check how many spikes in this window?
      spike_i = spike.spike_i{ri};
      in_idx = find((spike_i>=w_left) & (spike_i<=w_right));
      s_c(ri) = length(in_idx);
    end
    s_count(si).spike_counts = s_c;
    s_count(si).has_spike = sum(s_c>0);
    s_count(si).frac_has_spike = s_count(si).has_spike / s_count(si).num_syl; 
    s_count(si).mean_num_spike = mean(s_c);
    
    % check if pass criteria
    s_count(si).is_pass = (s_count(si).frac_has_spike>=p_thre) & (s_count(si).mean_num_spike>=m_thre);
  end
  
  % add info to the table
  info_count.s_count{ni} = s_count; 
  info_count.is_pass{ni} = [s_count.is_pass];
  info_count.birdID{ni} = birdID;
end

info_all = [info_all; info_count];

end


fn_save = fullfile(fd_save, 'info_all.mat');
save(fn_save, 'info_all');



%% Plot the co-occurrence of firing across syllable types
X = vertcat(info_all.is_pass{:});
check_i = [1 3 4 6];
% check_i = [1 3 4];
X = X(:, check_i);
% labels = syls(check_i);
labels = syls_pretty(check_i);

close all;
[fig, axs] = ZZfunc_UpSetPlot_v1(X, labels);
title(axs(1), sprintf('Total n=%d', size(X,1)));
% save plots
fn_fig = fullfile(fd_save, 'co-occurence.category.pdf');
print(fig, fn_fig, '-dpdf', '-painters');


%% calculate mean number of syllable types
sum_type = sum(X, 2); 
mean_type = mean(sum_type);
sem_type = std(sum_type) / sqrt(length(sum_type));
fprintf('Mean number of types: %.2f +- %.2f\n', mean_type, sem_type);


%% how many neurons fire to both calls and click/squawk
both_idx = find(X(:,1) & (X(:,2) |  X(:,3)));
n = length(both_idx);
fprintf('Both call and click/squak: n=%d, percent=%.2f\n', n, n/size(X,1));


% examine example neurons that fire both to calls and squawks
both_idx2 = find(X(:,1) & X(:,3));
a = info_all(both_idx2, :);






