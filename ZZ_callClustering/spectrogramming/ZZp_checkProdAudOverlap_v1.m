% check if I accidentally mislabeled the same nc file as both 'bProd' and 'bAud'
% Zhilei, 06/24/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% where call embedding results are stored
fd_embed_base = fullfile(fd_home, 'vaeWav');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what dbase type to check
suffix = 'Wsp1';


% loop through birds
bi = 4;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load information about sorted neurons
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
load(fn_info);
% get unique dates
date_unique = unique(info.date_long, 'sorted');

%% 2. loop through dates
for di=1:size(date_unique, 1)
  dd = date_unique{di};
  d_short = strrep(dd, '-', '');
  fd_dbase = fullfile(fd_home, 'DbaseFiles', pairID, dd, birdID, 'warble');
  fn_dbase = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, d_short, suffix);
  load(fullfile(fd_dbase, fn_dbase));
  % which field is bProd or bAud
  p_prod = find(strcmp(dbase.PropertyNames, 'bProd'));
  p_aud = find(strcmp(dbase.PropertyNames, 'bAud'));
  % find files that are both bProd and bAud
  p = dbase.Properties;
  i_both = find(p(:,p_prod) & p(:,p_aud));
  if ~isempty(i_both)
    fprintf('%s %s found overlap\n', birdID, dd);
    fprintf('%d ', i_both);
    disp('\n');
  else
    fprintf('%s %s no overlap!\n', birdID, dd);
  end
end













