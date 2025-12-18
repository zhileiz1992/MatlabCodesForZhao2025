% quantify the variation of spike latency in compound syllables

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));



%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};

% where to save results
fd_save = fullfile(fd_base, 'Figures', 'CombinedAnalysis', 'latencyVariation');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end

cv_all = [];
%% Loop through birds, load data
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

% only look at compound syllables
comp = spike_win(ismember(spike_win.category, {'b', 'x'}),:);
comp = comp(comp.dur>=0.3,:);

% calculate latency to the first spike
fs = 20000;
for ri=1:size(comp, 1)
  if ~isempty(comp.spike_i{ri})
    istart = comp.r_start{ri};
    ispike = comp.spike_i{ri};
    lat = (ispike - istart) / fs;
    % pick the first spike after syllable onset
    ipick = find(lat>=0 & lat<=comp.dur(ri));
    if ~isempty(ipick)
      comp.latency(ri) = lat(ipick(1));
    else
      comp.latency(ri) = nan;
    end
  else
    comp.latency(ri) = nan;
  end
end


% loop through each neuron, calculate CV
cv_list = info_neu; 
for ni=1:size(info_neu,1)
  neuID = info_neu.neuronID{ni};
  comp_this = comp(strcmp(comp.neuronID, neuID),:);
  % calculate cv
  lat_list = comp_this.latency;
  % filter out the nan
  lat_good = lat_list(~isnan(lat_list));
  cv = nan;
  if length(lat_good)>=10
    cv = std(lat_good) / mean(lat_good);
  end
  cv_list.cv(ni) = cv;
  cv_list.num_rends(ni) = length(lat_good);
  cv_list.birdID{ni} = birdID;
end
cv_all = [cv_all; cv_list];
end
  
% calculate mean CV and SEM
cv_array = cv_all.cv;
cv_array = cv_array(~isnan(cv_array));
mean_cv = mean(cv_array);
sem = std(cv_array) / sqrt(length(cv_array));
fprintf('CV mean += SEM: %.3f += %.3f\n', mean_cv, sem);

























