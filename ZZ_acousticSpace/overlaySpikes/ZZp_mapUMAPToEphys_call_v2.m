% Map the results from UMAP to Ephys struct for calls
% Zhilei, 09/08/2025
% differ from v1: use the replaced call subtype identities

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
% where is the ephys data struct
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedReplaced2');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.replaced2.mat', birdID)));
% what's the window size
win_frame = 32;
% where is VAE/UMAP embedding located
vae_run = 'traj_chop_32_1_32';
% fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', sprintf('latents.%s', vae_run));
umap_run = 'umapAll.v4v5';
fd_umap = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', umap_run);
% what syllables to map, need to match what UMAP is trained on
syls = {'v4', 'v5'};
syls_new = {'v1', 'v2'};

% where to save results
fd_save =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save UMAP embedding to %s\n', fd_save);


%% 1. Loop through ephys struct, map UMAP results
% save as a separate data struct, with matched row index
% load the meta info about syllable sliding 
fn_info = fullfile(fd_umap, sprintf('%s.%s.info.csv', birdID, umap_run));
opts = detectImportOptions(fn_info);
opts = setvartype(opts, 'call_subtype', 'string');
opts.Delimiter = ',';
info = readtable(fn_info, opts);
% what umap results to use
fn_umap = fullfile(fd_umap, sprintf('%s.%s.embedding.csv', birdID, umap_run));
d_umap = readmatrix(fn_umap);

% loop through ephys struct
for si=2:size(syls_new, 2)
  v = syls_new{si};
  fn_e = fullfile(fd_ephys, sprintf('%s.%s.segments_all.replaced2.mat', birdID, v));
  fprintf('Loading data for %s\n', v);
  a = load(fn_e);
  seg_selected = a.seg_selected;

  umap = [];
  parfor ri=1:size(seg_selected,2)
%   for ri=1:100
    % locate the syllable by syl_ID
    seg = seg_selected(ri);
    fn_wav = strsplit(seg.fn_audio, '/');
    fn_wav = strrep(fn_wav{end}, '.nc', '.wav');
    i_loc = strsplit(seg.title_str, '-');
    i_loc = str2num(i_loc{2}) - 1;  % matlab is 1-based
    syl_ID = sprintf('%s_%d_%d_%d', fn_wav, i_loc, seg.seg_start_ori, seg.seg_end_ori);
    % read in the embeding data, determine start and end index
    this_meta = info(strcmp(info.syl_ID, syl_ID),:);
    ii_start = this_meta.count_start + 1;  % note python is 0-based
    ii_end = this_meta.count_end;
    this_umap = d_umap(ii_start:ii_end,:);
    % locate the meta info on sliding as well
    umap(ri).index = ri;
    umap(ri).syl_ID = syl_ID; 
    umap(ri).umap = this_umap;
    umap(ri).umap_meta = this_meta;
  end
  
  % check if there is any entries where data is not found
  size_umap = cellfun(@(x) size(x, 2), {umap.umap});
  i_miss = find(size_umap==0);
  fprintf('Data not found for: %d\n', i_miss);
  
  % save results
  fn_save = fullfile(fd_save, sprintf('%s.%s.%s.replaced2.mat', birdID, v, umap_run));
  save(fn_save, 'umap');
end
  




