% rename the ephy struct to the new call subtype names
% read in the old segments_all.pull data struct, then change to the new naming
% 09/25/2025
% differ from v2: add a syl_ID field

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% 1. Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
to_removes = {{'v6'}; {'v5', 'v6'}; {'v5','v7'}; {'v3','v4','v5'}};  % what call subtypes to remove (oldID)
% rename relationship
names_old = {{'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}, ...
  {'v1', 'v2', 'v3', 'v4', 'v5', 'v6'}, ...
  {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}, ...
  {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'}};

names_new = {{'v3', 'v6', 'v5', 'v1', 'v2', 'v7', 'v4'}, ...
  {'v1', 'v4', 'v3', 'v2', 'v5', 'v6'}, ...
  {'v4', 'v1', 'v2', 'v3', 'v7', 'v5', 'v6'}, ...
  {'v1', 'v4', 'v5', 'v6', 'v7', 'v3', 'v2'}};



%% 2. Loop through birds, convert clustering result to dbase
% bi = 3;
for bi=2:4
  bd = birdIDs{bi};
  no = names_old{bi};
  nn = names_new{bi};
  to_rm = to_removes{bi};
  fd_save = fullfile(fd_home, 'Figures', pairIDs{bi}, 'HahnloserNew', bd, 'extractedReplaced2');
  if ~exist(fd_save); mkdir(fd_save); end
  
  % loop through ephys struct
  for syl_i=1:size(no,2)
    v_old = no{syl_i};
    v_new = nn{syl_i};
    %   if ismember(v_old, to_rm); continue; end
    % read in old ephys struct
    fd_old = fullfile(fd_home, 'Figures', pairIDs{bi}, 'HahnloserNew', bd, 'extractedPull');
    fn_old = fullfile(fd_old, sprintf('%s.%s.segments_all.pull.mat', bd, v_old));
    load(fn_old);
    for ri=1:size(seg_selected,2)
      seg_selected(ri).title = v_new;
      seg_selected(ri).title_str = strrep(seg_selected(ri).title_str, v_old, v_new);
      % add a syl_ID field
      fn_wav = strsplit(seg_selected(ri).fn_audio, '/');
      fn_wav = strrep(fn_wav{end}, '.nc', '.wav');
      i_loc = strsplit(seg_selected(ri).title_str, '-');
      i_loc = str2num(i_loc{2}) - 1;  % matlab is 1-based
      syl_ID = sprintf('%s_%d_%d_%d', fn_wav, i_loc, seg_selected(ri).seg_start_ori, seg_selected(ri).seg_end_ori);
      seg_selected(ri).syl_ID = syl_ID;
    end
    % save
    fn_new = fullfile(fd_save, sprintf('%s.%s.segments_all.replaced2.mat', bd, v_new));
    save(fn_new, 'seg_selected', '-v7.3');
  end
  
end








