% a script to pull extracted ephys data for calls
% first sort batch + supplemented + manual
% also filter out duplicated neurons
% Zhilei, 07/07/2025


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
to_remove = {{}, {'20240902-ch15'}, {}, {}};  % what neurons to remove due to duplication
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% save in subfolder for each bird-call subtype
save_dir = 'extractedPull';


% loop through birds
% bi = 1;
% for bi=1:size(birdIDs, 2)
for bi=3:3
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load information about extracted ephys data
fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
% grab what calls have been extracted for the first sort batch
fds = dir(fullfile(fd_save_master, birdID, 'extracted', '*segments_all.mat'));
% where to save pulled results
fd_save_this = fullfile(fd_save_master, birdID, 'extractedPull');
if ~exist(fd_save_this, 'dir')
  mkdir(fd_save_this)
end


%% 2. Loop through call subtype, pull data and de-duplicate
for li=1:size(fds,1)
  % load previous saved data struct
  temp = strsplit(fds(li).name, '.');
  syl_label = temp{2};
  fprintf('Pulling for %s %s...\n', birdID, syl_label);
  load(fullfile(fds(li).folder, fds(li).name));
  % for production choose chan0, for auditory choose chan17
  pt = 'chan0';
  seg_selected = segments_all(strcmp({segments_all.aud_ch}, pt));
  % add supplemented sorting if exist
  fn_add = fullfile(fd_save_master, birdID, 'extractedAdded', sprintf('%s.%s.segments_all.added.mat', birdID, syl_label));
  if exist(fn_add, 'file')
    a = load(fn_add);
    seg_add = a.segments_all;
    seg_selected = [seg_selected seg_add];
  end
  % add manual sorting if exist
  fn_manual = fullfile(fd_save_master, birdID, 'extractedAdded', sprintf('%s.%s.segments_all.manual.mat', birdID, syl_label));
  if exist(fn_manual, 'file')
    a = load(fn_manual);
    seg_manual = a.segments_all;
    seg_selected = [seg_selected seg_manual];
  end
    
  % remove duplication
  if ~isempty(to_remove{bi})
    idx_pass = find(~ismember({seg_selected.neuronID}, to_remove{bi}));
    seg_selected = seg_selected(idx_pass);
  end

  % save
  fn_save = fullfile(fd_save_this, sprintf('%s.%s.segments_all.pull.mat', birdID, syl_label));
  save(fn_save, 'seg_selected', '-v7.3');
end


end













