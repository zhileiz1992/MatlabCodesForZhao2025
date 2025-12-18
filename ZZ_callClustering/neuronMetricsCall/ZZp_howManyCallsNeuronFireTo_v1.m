% count how many calls a neuron can fire to
clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what call clusters are good clusters
good_call = {{'v1', 'v2', 'v3', 'v4', 'v5', 'v6'}, ...
             {'v1', 'v2', 'v3', 'v4'}, ...
             {'v1', 'v2', 'v3'}, ...
             {'v1', 'v2', 'v3', 'v4'}};


% loop through birds
% bi = 2;
info_all = [];
for bi=1:size(birdIDs,2)
  birdID = birdIDs{bi};
  pairID = pairIDs{bi};
  g_c = good_call{bi};


  %% 1. Load information about sorted neurons and call subtypes
  fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
  fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
  load(fn_info);
  % add a column to count number of calls a neuron fire to
  info.count = zeros(size(info,1), 1);
  
  % load information on what neurons fire for what 
  fd_res = fullfile(fd_save_master, birdID, 'popRaster2');
  
  for vi=1:size(g_c,2)
    v = g_c{vi};
    fn_neu = fullfile(fd_res, v, sprintf('Hahnloser-%s-chan0.neuron_orderedPlotted7.mat', v));
    load(fn_neu);
    for ni=1:size(neuron_ordered,2)
      idx = find(strcmp(info.neuronID, neuron_ordered{ni}));
      info.count(idx) = info.count(idx) + 1;
    end
  end
  info_all = [info_all; info];
end
    
  
%% Calculate mean number
a = info_all.count;
fprintf('Mean number of calls (including zero: %.2f\n', mean(a));
b = a(a~=0);
m = mean(b);                    % mean
sem = std(b) / sqrt(length(b));
fprintf('Mean number of calls (excluding zero: %.2f +- %.2f\n', m, sem);

  
  
  
  
  
  
  
