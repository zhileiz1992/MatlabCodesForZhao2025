% retrive the VAE data for selecte call subtypes in the first sorting batch
% save as .mat for later reuse
% differ from v1: use the replaced names where low-quality clusters are filtered out

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
names_olds = {{'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}};
             
names_news = {{'v3', 'v6', 'v5', 'v1', 'v2', 'v7', 'v4'}, ...
             {'v1', 'v4', 'v3', 'v2', 'v5', 'v6'}, ...
             {'v4', 'v1', 'v2', 'v3', 'v7', 'v5', 'v6'}, ...
             {'v1', 'v4', 'v5', 'v6', 'v7', 'v3', 'v2'}};


for bi=3:size(birdIDs,2)

birdID = birdIDs{bi};
pairID = pairIDs{bi};
% what's the window size
win_frame = 32;
% where is VAE/UMAP embedding located
vae_run = 'traj_chop_32_1_32';
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', sprintf('latents.%s', vae_run));
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));
fn_info = fullfile(fd_vae, sprintf('%s.latents.%s.info.csv', birdID, vae_run));
info = readtable(fn_info, 'Delimiter', ',');
% use data from batch1, since batch 2 may have calls from partner bird
batch = 'batch1';
% where to save
fd_save =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'sylOnly3', batch);
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save MMD results to %s\n', fd_save);


%% 1. Replace the old call subtype names in the table
names_old = names_olds{bi};            
names_new = names_news{bi};
           
info2 = info;
for ri=1:size(info,1)
  no = info.call_subtype{ri};
  ni = find(strcmp(names_old, no));
  if ~isempty(ni)
    nn = names_new{ni};
    info2.call_subtype{ri} = nn;
  end
end


%% 2. Retrieve VAE raw data for all syllable type, save as .mat
% what syllable subtypes this bird has
temp = info(strcmp(info.batch, batch),:);
syls = unique(temp.call_subtype, 'sorted');
disp(syls);
for si=2:size(syls,1)
  ss = syls{si};
  fprintf('Retrieve VAE for: %s\n', ss);
  d_vae = ZZfunc_retrieveVAE_v1(fn_vae, info2, {batch}, {ss});
  dvae = d_vae{1};
  % save as .mat
  fn_save = fullfile(fd_save, sprintf('%s.%s.%s.dvae_raw.mat', birdID, vae_run, ss));
  save(fn_save, 'dvae', '-v7.3');
end
  
  
end
  
  
  
  
  
  
  