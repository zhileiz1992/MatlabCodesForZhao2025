% rename the ephy struct to the new call subtype names
% read in the old segments_all.pull data struct, then change to the new naming
% 09/02/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% 1. Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% rename relationship
names_old = {{'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}};
             
names_new = {{'v3', 'v7', 'v5', 'v1', 'v2', 'v6', 'v4'}, ...
             {'v1', 'v6', 'v3', 'v2', 'v4', 'v5'}, ...
             {'v7', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6'}, ...
             {'v1', 'v4', 'v5', 'v6', 'v7', 'v3', 'v2'}};
% what suffix to use for saving dbase
suffix = 'Wsp3Call';



%% 2. Loop through birds, convert clustering result to dbase
% bi = 3;
for bi=3:4
bd = birdIDs{bi};
no = names_old{bi};
nn = names_new{bi};
fd_save = fullfile(fd_home, 'Figures', pairIDs{bi}, 'HahnloserNew', bd, 'extractedReplaced');
if ~exist(fd_save); mkdir(fd_save); end

% loop through ephys struct
for syl_i=1:size(no,2)
  v_old = no{syl_i};
  v_new = nn{syl_i};
  % read in old ephys struct
  fd_old = fullfile(fd_home, 'Figures', pairIDs{bi}, 'HahnloserNew', bd, 'extractedPull');
  fn_old = fullfile(fd_old, sprintf('%s.%s.segments_all.pull.mat', bd, v_old));
  load(fn_old);
  for ri=1:size(seg_selected,2)
    seg_selected(ri).title = v_new;
    seg_selected(ri).title_str = strrep(seg_selected(ri).title_str, v_old, v_new);
  end
  % save
  fn_new = fullfile(fd_save, sprintf('%s.%s.segments_all.replaced.mat', bd, v_new));
  save(fn_new, 'seg_selected', '-v7.3');
end  
  
end
  
  
  
  
  
  
  
  
