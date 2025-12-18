% calculate and plot the firing field of MO neurons in the acoustic space
% step 3: Run UMAP using selected inputs
% differ from v1: use the unified VAE data struct, instead of separate struct for each syllable


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


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird specific folder setting
fd_vae = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', vae_run);
fn_vae = fullfile(fd_vae, sprintf('%s.dvae.mat', birdID));

% where to save UMAP results
fd_save =  fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', sprintf('%s_umap', vae_run));
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save UMAP results to %s\n', fd_save);


%% 1. Select inputs to the UMAP
% load the VAE data struct
load(fn_vae);
% only uses those have ephys data
d = dvae([dvae.has_ephys]==1);

% 1.1. Use all syllables
umap_run = 'all_syl';
lens = cellfun(@(x) size(x,1), {d.vae});
d_input = cat(1, d.vae);



% save as a matrix
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

% also save the ID info
d_info = struct('ri', {d.ri}, 'syl_ID', {d.syl_ID});
fn_d_info = fullfile(fd_save, sprintf('%s.%s.d_info.csv', birdID, umap_run));
writetable(struct2table(d_info), fn_d_info);















