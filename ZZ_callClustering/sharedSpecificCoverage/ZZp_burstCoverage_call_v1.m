% analyze the burst converage for shared and syllable-specific MO neurons for call subtype pairs
% Zhilei, 07/11/2025


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what call subtype pairs to analyze
to_analyze = {{{'v4', 'v5'}, {'v1', 'v7'}}, {}, {}, {}};

% loop through birds
bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
fd_save_base = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew', birdID);

% loop through analyze pairs
pi = 1; 
p = to_analyze{bi}{pi};
% grab the criteria data struct for the call pairs
fn1 = fullfile(fd_save_base, p{1}, sprintf('Hahnloser-%s-chan0.criteria5.mat', p{1})); 
c1 = load(fn1); c1 = c1.criteria; 
fn2 = fullfile(fd_save_base, p{2}, sprintf('Hahnloser-%s-chan0.criteria5.mat', p{2})); 
c2 = load(fn2); c2 = c2.criteria; 

% find the shared and specific neurons
idx1 = find([c1.isPass]);
n1 = {c1(idx1).neuronID};
idx2 = find([c2.isPass]);
n2 = {c2(idx2).neuronID};
shared = intersect(n1, n2);
specific1 = setdiff(n1, n2);
specific2 = setdiff(n2, n1);

% find the burst location of these neurons
ns = {shared; specific1; specific2};
t_all = struct();
count = 0;
for ii=1:size(ns, 1)
  n = ns{ii};
  for ni=1:size(n, 2)
    % if shared take average
    i1 = find(strcmp({c1.neuronID}, n{ni}));
    i2 = find(strcmp({c2.neuronID}, n{ni}));
    t1 = nan; t2 = nan;
    if ~isempty(i1)
      t1 = c1(i1).max_psth_t;
    end
    if ~isempty(i2)
      t2 = c2(i2).max_psth_t;
    end
    count = count+1;
    t_all(count).neuronID = n{ni};
    t_all(count).isPass1 = c1(i1).isPass;
    t_all(count).isPass2 = c2(i2).isPass;
    t_all(count).t1 = t1;
    t_all(count).t2 = t2;
  end
end






















