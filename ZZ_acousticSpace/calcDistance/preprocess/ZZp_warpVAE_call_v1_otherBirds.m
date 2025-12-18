% For calls, time-warp the VAE latents 
clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/plotTraj'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
to_removes = {{'v6'}; {'v5', 'v6'}; {'v5','v7'}; {'v3','v4','v5'}};  % what call subtypes to remove (oldID)
% rename relationship
syls_all = {{'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}};
             

for bi=3:size(birdIDs,2)

birdID = birdIDs{bi};
pairID = pairIDs{bi};
% what syllable subtypes this bird has
syls = syls_all{bi};
% what's the window size
win_frame = 32;
% where is VAE/UMAP embedding located
vae_run = 'traj_chop_32_1_32';
% where is retrieved VAE data stored
fd_home = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed');
% fd_vaeRes = fullfile(fd_home, 'sylOnly', 'batch1');
% fd_vaeRes = fullfile(fd_home, 'sylOnly2', 'batch1');
fd_vaeRes = fullfile(fd_home, 'sylOnly3', 'batch1');
fprintf('Save MMD results to %s\n', fd_vaeRes);
% what suffix to use when saving
suffix = 'dvae_warp';


%% 1. Retrieve VAE data, time-warp the VAE latents (32-dim trajectories)
% do it separately for each syllable type
for vi=1:size(syls,2)
  fn = fullfile(fd_vaeRes, sprintf('%s.%s.%s.dvae_raw.mat', birdID, vae_run, syls{vi}));
  a = load(fn);
  d_syl = a.dvae; 
  % warp to the mean len
  lens = cellfun(@(x) size(x,1), {d_syl.vae});
  mean_len = floor(mean(lens));
  fprintf('Intepolate %s to %d frames (ms)\n', syls{vi}, mean_len);
  parfor ri=1:size(d_syl,2)
    d = d_syl(ri).vae;
    interp_d = ZZfunc_interp_window(d, mean_len);
    d_syl(ri).vae_interp = interp_d;
  end
  % save as separate .mat for later reuse
  dvae = d_syl;
  fn_save = fullfile(fd_vaeRes, sprintf('%s.%s.%s.%s.mat', birdID, vae_run, syls{vi}, suffix));
  save(fn_save, 'dvae', '-v7.3');
end


end

