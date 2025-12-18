% calculate and plot the firing field of MO neurons in the acoustic space
% step 3: Run UMAP using selected inputs
% differ from v3: for very large datasets (e.g. bird M2), run UMAP on half of the datasets


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what's the sliding window size
win_frame = 32;
% what VAE run to use
vae_run = 'traj_chop_32_1_32';


% bi = 1;
for bi = 2:4
  
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird specific folder setting
fd_vae = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', vae_run);

% where to save UMAP results
fd_save =  fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', sprintf('%s_umap', vae_run));
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save UMAP results to %s\n', fd_save);


%% 1. Select inputs to the UMAP
% only include those that have ephys data

% 1.1. Use all syllables
% umap_run = 'all_syl';
umap_run = 'all_syl2';
syls = {'v', 'b', 'x', 'h', 'e'};

% 1.2 Exclude h/e
% umap_run = 'no_eh';
% syls = {'v', 'b', 'x'};

% load information about sorted neuron
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% what neurons to analyze
neu_analyze = info_neu(82:121,:);
date_unique = unique(neu_analyze.date);


info = [];
d_input = [];
for si=1:size(syls,2)
  fprintf('Loading for %s...\n', syls{si});
  fn_vae = fullfile(fd_vae, sprintf('%s.%s.dvae.mat', birdID, syls{si}));
  load(fn_vae);
  dvae = dvae([dvae.has_ephys]==1);
  % only keep those that in the date range
  date_info = cellfun(@(x) strsplit(x,'_'), {dvae.syl_ID}, 'UniformOutput', false);
  date_info = cellfun(@(x) x{end-5}, date_info, 'UniformOutput', false);
  date_info = cellfun(@(x) ['20' x], date_info, 'UniformOutput', false);
  idx = find(ismember(date_info, date_unique));
  dvae = dvae(idx);
  lens = cellfun(@(x) size(x,1), {dvae.vae}, 'UniformOutput', false);
  info_this = struct('syl_ID', {dvae.syl_ID}, 'lens', lens, 'category', syls{si});
  d = cat(1, dvae.vae);
  d_input = [d_input; d];
  info = [info info_this];
  clear dvae d info_this
end

% save info for later use
fn_info = fullfile(fd_save, sprintf('%s.%s.info.mat', birdID, umap_run));
save(fn_info, 'info', '-v7.3');
% save data as a matrix for UMAP input
fn_d = fullfile(fd_save, sprintf('%s.%s.d_input.csv', birdID, umap_run));
writematrix(d_input, fn_d);



%% 2. Run UMAP
param.n_components = 2;
param.n_neighbors = 25; 
param.min_dist = 0;
% param.metric = 'cosine';
param.metric = 'euclidean';
param.random_state = 1118;
fn_save = fullfile(fd_save, sprintf('%s.%s.umap_res.csv', birdID, umap_run));
% don't save the umap model 
fn_p = '';
umap_res = ZZfunc_runUMAP_v1(fn_d, fn_save, fn_p, param);



end












