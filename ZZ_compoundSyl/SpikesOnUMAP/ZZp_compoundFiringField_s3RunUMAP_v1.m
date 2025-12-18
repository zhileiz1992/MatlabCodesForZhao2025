% calculate and plot the firing field of MO neurons in the acoustic space
% step 3: Run UMAP using selected inputs


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
% where is the ephys data struct for call syllables
fd_ephys_call = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedReplaced2');
fns_ephys_call = dir(fullfile(fd_ephys_call, sprintf('%s.*.segments_all.replaced2.mat', birdID)));
% where is the ephys data struct for non-v syllables
fd_ephys_other = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'extractedReplaced2');
fns_ephys_other = dir(fullfile(fd_ephys_other, sprintf('%s.*.segments_all.replaced2.mat', birdID)));
fns_ephys = [fns_ephys_call; fns_ephys_other];

% where is VAE latents located
fd_vae_call = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
fns_vae_call = dir(fullfile(fd_vae_call, sprintf('%s.*.vae_latents.replaced2.mat', birdID)));
fd_vae_other = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', vae_run);
fns_vae_other = dir(fullfile(fd_vae_other, sprintf('%s.*.vae_latents.replaced2.mat', birdID)));
fns_vae = [fns_vae_call; fns_vae_other];

% where to save results
fd_save =  fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', sprintf('%s_umap', vae_run));
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save UMAP results to %s\n', fd_save);



%% 1. Select inputs to the UMAP
% note that the same syllable may occur many times in the ephys struct and VAE struct
% need to de-duplicate: use the syl_ID as an identifier 
ID_vae = [];
for si=1:size(fns_vae, 1)
  fn = fullfile(fns_vae(si).folder, fns_vae(si).name);
  v = strsplit(fns_vae(si).name, '.');
  v = v{2};
  fprintf('Processing for %s\n', v);
  a = load(fn); latent = a.latent; 
  syl_IDs = unique({latent.syl_ID}, 'stable');
  syl_ID_ori = {latent.syl_ID};
  ID_vae_this = [];
  parfor ri=1:size(syl_IDs,2)
    syl_ID = syl_IDs{ri};
    idx = find(strcmp(syl_ID_ori, syl_ID));
    idx = idx(1);  % de-duplicate
    ID_vae_this(ri).syl_ID = syl_ID;
    ID_vae_this(ri).label = v; 
    ID_vae_this(ri).category = v(1);
    ID_vae_this(ri).vae = latent(idx).vae; 
  end
  ID_vae = [ID_vae ID_vae_this];
end  
% save for direct loading in the future




















