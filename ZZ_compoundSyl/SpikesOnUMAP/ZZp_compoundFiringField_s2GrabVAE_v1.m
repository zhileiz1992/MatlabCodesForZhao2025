% calculate and plot the firing field of MO neurons in the acoustic space
% step 2: Grab the VAE latents for each entry in the ephys struct
% this step has been done for all call syllables


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what syllables to extract
labels_extract = {'b', 'x', 'h', 'e'};
% what's the sliding window size
win_frame = 32;
% what VAE run to use
vae_run = 'traj_chop_32_1_32';


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird specific folder setting
% where is the ephys data struct
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'extractedReplaced2');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.*.segments_all.replaced2.mat', birdID)));
% where is VAE/UMAP embedding located
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', sprintf('latents.%s', vae_run));

% where to save results
fd_save =  fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', vae_run);
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save VAE latents to %s\n', fd_save);


%% 1. Loop through ephys struct, map VAE results
% save as a separate data struct, with matched row index
% load the meta info about syllable sliding 
fn_info = fullfile(fd_vae, sprintf('%s.latents.%s.info.csv', birdID, vae_run));
opts = detectImportOptions(fn_info);
opts = setvartype(opts, 'call_subtype', 'string');
opts.Delimiter = ',';
info = readtable(fn_info, opts);
% what vae results to use
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));

% loop through ephys struct
for si=2:size(fns_ephys, 1)
  fn_e = fullfile(fns_ephys(si).folder, fns_ephys(si).name);
  v = strsplit(fns_ephys(si).name, '.');
  v = v{2};
  fprintf('Loading data for %s\n', v);
  a = load(fn_e);
  seg_selected = a.seg_selected;

  latent = [];
  parfor ri=1:size(seg_selected,2)
%   for ri=1:100
    % locate the syllable by syl_ID
    seg = seg_selected(ri);
    fn_wav = strsplit(seg.fn_audio, '/');
    fn_wav = strrep(fn_wav{end}, '.nc', '.wav');
    i_loc = strsplit(seg.title_str, '-');
    i_loc = str2num(i_loc{2}) - 1;  % matlab is 1-based
    syl_ID = sprintf('%s_%d_%d_%d', fn_wav, i_loc, seg.seg_start_ori, seg.seg_end_ori);
    % read in the embeding data
    this_vae = h5read(fn_vae, ['/' syl_ID]);
    % locate the meta info on sliding as well
    this_meta = info(strcmp(info.syl_ID, syl_ID),:);
    latent(ri).index = ri;
    latent(ri).syl_ID = syl_ID; 
    latent(ri).vae = this_vae';
    latent(ri).vae_meta = this_meta;
  end
  
  % check if there is any entries where data is not found
  size_vae = cellfun(@(x) size(x, 2), {latent.vae});
  i_miss = find(size_vae==0);
  fprintf('Data not found for: %d\n', i_miss);
  
  % save results
  fn_save = fullfile(fd_save, sprintf('%s.%s.vae_latents.replaced2.mat', birdID, v));
  save(fn_save, 'latent', '-v7.3');
  
  clear latent seg_selected;
end
  




