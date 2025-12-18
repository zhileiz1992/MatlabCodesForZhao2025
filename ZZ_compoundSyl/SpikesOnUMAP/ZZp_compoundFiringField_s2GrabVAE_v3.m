% calculate and plot the firing field of MO neurons in the acoustic space
% step 2: Grab the VAE latents for all syllables
% differ from v2: use syl_ID as unique ID; save separately for each syllable type; save for v0 as well


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


bi = 4;
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

% where is VAE embedding located
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', sprintf('latents.%s', vae_run));

% where to save results
fd_save =  fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', vae_run);
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save VAE latents to %s\n', fd_save);


%% 1. Grab syllable categories
% load the meta info about syllable sliding
fn_info = fullfile(fd_vae, sprintf('%s.latents.%s.info.csv', birdID, vae_run));
opts = detectImportOptions(fn_info);
opts = setvartype(opts, 'call_subtype', 'string');
opts.Delimiter = ',';
info = readtable(fn_info, opts);
% what vae results to use
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));

% add a label for the categories
temp = cellfun(@(x) x(1), info.call_subtype, 'UniformOutput', false);
info.category = temp;

% grab syllable categories
cate = unique(info.category, 'sorted');


%% 2. Go through the ephy struct, mark which syl_ID have ephys data
% this helps remove the non-focal syllables in the supplementary sort batch
syl_ID_ephys = {};
% loop through ephys struct
for si=1:size(fns_ephys, 1)
  fn_e = fullfile(fns_ephys(si).folder, fns_ephys(si).name);
  v = strsplit(fns_ephys(si).name, '.');
  v = v{2};
  fprintf('Loading data for %s\n', v);
  a = load(fn_e);
  seg_selected = a.seg_selected;
  
  id_all = {seg_selected.syl_ID};
  syl_ID_ephys = [syl_ID_ephys id_all];
end
clear seg_selected;



%% 3. Loop through each category, save vae data
for si=1:size(cate,1)
  syl = cate{si};
  info_syl = info(strcmp(info.category, syl),:);
  % check if there is any duplicate
  temp = unique(info_syl.syl_ID);
  if size(temp,1)~=size(info_syl,1)
    fprintf('Found duplicates for %s\n', syl);
  end
  
  % grab vae data
  dvae = [];
  parfor ri=1:size(info_syl, 1)
    syl_ID = info_syl.syl_ID{ri};
    this_vae = h5read(fn_vae, ['/' syl_ID]);
    this_vae = this_vae';
    dvae(ri).ri = ri;
    dvae(ri).syl_ID = syl_ID;
    dvae(ri).vae = this_vae;
    this_meta = info_syl(strcmp(info_syl.syl_ID, syl_ID),:);
    dvae(ri).vae_meta = this_meta;
    dvae(ri).has_ephys = ismember(syl_ID, syl_ID_ephys);
  end
  
  fn_mat = fullfile(fd_save, sprintf('%s.%s.dvae.mat', birdID, syl));
  save(fn_mat, 'dvae', '-v7.3');
end









