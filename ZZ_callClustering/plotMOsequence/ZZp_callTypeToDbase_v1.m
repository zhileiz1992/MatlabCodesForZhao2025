% convert the results from VAE/UMAP/HDBSCAN call clustering into dbase of each date
% Zhilei, 06/17/2025
% Call subtype ID starts with V1
% Calls that don't belong to any cluster will be labeled as V0
% Calls that are excluded from clustering analysis (short) will be still labeled as V 

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% 1. Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
fd_save_base = fullfile(fd_home, 'vaeWav');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what spectrogram input dataset
input_rn = 'spec_goffinet_cutoff_256_176';
syl = 'v';  % what syllable to focus on
% what subfolder has the VAE/UMAP results
fd_vae = 'UMAPonVAE6'; 
vae_run = input_rn;
% what suffix to use for saving dbase
suffix = 'Wsp1Call';
% what dbase to use as inputs
suffix_in = 'Wsp1'; 


%% 2. Loop through birds, convert clustering result to dbase
bi = 4;
bd = birdIDs{bi};

% read the clustering results
fd_d = fullfile(fd_save_base, bd, fd_vae, syl, vae_run);
fn_embed = fullfile(fd_d, sprintf('%s.%s.embedding.csv', bd, vae_run)); 
embed = readtable(fn_embed);
% add a date column and filename column
splitParts = cellfun(@(s) strsplit(s, '/'), embed.fn_wav, 'UniformOutput', false);
% Step 2: Extract the second-to-last element from each split array.
embed.date = cellfun(@(parts) parts{end-1}, splitParts, 'UniformOutput', false);
embed.fn = cellfun(@(parts) parts{end}, splitParts, 'UniformOutput', false);

% loop through dates
dates_unique = unique(embed.date, 'stable');
for di=1:size(dates_unique,1)
  dd = dates_unique{di};
  disp(dd);
  if strcmp(dd, 'unsorted')
    continue
  end
  dd_short = strrep(dd, '-', '');
  % load the input dbase
  fd_dbase = fullfile(fd_home, 'DbaseFiles', pairIDs{bi}, dd, bd, 'warble');
  fn_dbase = sprintf('%s.%s.warble.good.%s.dbase.mat', bd, dd_short, suffix_in);
  a = load(fullfile(fd_dbase, fn_dbase)); 
  dbase_in = a.dbase;
  
  % replace the syllable labels with call clustering results
  dbase_out = dbase_in; 
  embed_s = embed(strcmp(embed.date, dd), :);
  % loop through each file
  fns_unique = unique(embed_s.fn, 'stable');
  fns_nc = {dbase_in.SoundFiles.name};
  for fi=1:size(fns_unique, 1)
    % replace .wav with .nc
    fn_this = strrep(fns_unique{fi}, '.wav', '.nc');
    f_loc = find(strcmp(fns_nc, fn_this));
    % replace the syllable labels
    segs = dbase_in.SegmentTitles{f_loc};
    rows = embed_s(strcmp(embed_s.fn, fns_unique{fi}),:);
    % note that the syllable index in embedding table starts with 0
    for ri=1:size(rows,1)
      s_loc = rows.s_idx(ri)+1;  
      segs{s_loc} = sprintf('v%d', rows.hdbscan_cluster(ri));
    end
    dbase_out.SegmentTitles{f_loc} = segs;
  end
  
  % save new dbase
  dbase = dbase_out;
  fn_save = sprintf('%s.%s.warble.good.%s.dbase.mat', bd, dd_short, suffix);
  save(fullfile(fd_dbase, fn_save), 'dbase');
end
  
      
    
  























