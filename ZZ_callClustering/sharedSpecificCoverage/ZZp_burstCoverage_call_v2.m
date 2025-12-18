% analyze the burst converage for shared and syllable-specific MO neurons for call subtype pairs
% Zhilei, 07/11/2025
% differ from v2: re-evaluate if a neuron pass or not by adding criteria on IRCC and sparseness


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
% also grab the original spike data
fn1 = fullfile(fd_save_base, 'extractedPull', sprintf('%s.%s.segments_all.pull.mat', birdID, p{1})); 
d1 = load(fn1); d1 = d1.seg_selected; 
fn2 = fullfile(fd_save_base, 'extractedPull', sprintf('%s.%s.segments_all.pull.mat', birdID, p{2})); 
d2 = load(fn2); d2 = d2.seg_selected; 


%% get data for all neurons
all_neuron = union({c1.neuronID}, {c2.neuronID});
t_all = struct();
% loop through each neuron
for ni=1:size(all_neuron, 2)
  neuronID = all_neuron{ni};
  % get the already calculated criteria info
  cn1 = c1(strcmp({c1.neuronID}, neuronID));
  cn2 = c2(strcmp({c2.neuronID}, neuronID)); 
  
  % save basic info to the data struct
  t_all(ni).neuronID = neuronID;
  t_all(ni).syl1 = p{1};
  t_all(ni).syl2 = p{2};
  t_all(ni).criteria1 = cn1;
  t_all(ni).criteria2 = cn2;
  if ~isempty(cn1) t_all(ni).t1 = cn1.max_psth_t; else t_all(ni).t1=nan; end
  if ~isempty(cn2) t_all(ni).t2 = cn2.max_psth_t; else t_all(ni).t2=nan; end
  
  % calculate sparseness
  psth_bin_size = 0.01;  % unit is seconds
  pad = 0.05;  % how much to look before and after the syllable, unit is seconds
  seg_this1 = d1(strcmp({d1.neuronID}, neuronID));
  if isempty(seg_this1)
    t_all(ni).sparse1 = nan;
  else
    [spa1, aligned_spike1] = ZZfunc_calcSparseness_v1(seg_this1, pad, psth_bin_size);
    t_all(ni).sparse1 = spa1;
  end
  seg_this2 = d2(strcmp({d2.neuronID}, neuronID));
  if isempty(seg_this2)
    t_all(ni).sparse2 = nan;
  else
    [spa2, aligned_spike2] = ZZfunc_calcSparseness_v1(seg_this2, pad, psth_bin_size);
    t_all(ni).sparse2 = spa2;
  end
  
  % determine if the neuron pass criteria for the syllable pair
  psth_thre = 10;
  min_ren = 4; 
  ircc_thre = 0.1; 
  sparse_thre = 0.1; 
  if isempty(cn1)
    t_all(ni).isPass1 = 0;
  else
    t_all(ni).isPass1 = (cn1.peak_psth>=psth_thre) &  (cn1.num_rend_psth_fire>=min_ren) & (cn1.ircc>=ircc_thre) & (spa1>=sparse_thre);
  end
  if isempty(cn2)
    t_all(ni).isPass2 = 0;
  else
    t_all(ni).isPass2 = (cn2.peak_psth>=psth_thre) &  (cn2.num_rend_psth_fire>=min_ren) & (cn2.ircc>=ircc_thre) & (spa2>=sparse_thre);
  end
end

%% determine if neuron is shared or specific
i_share = find([t_all.isPass1] & [t_all.isPass2]);
i_specific1 = find([t_all.isPass1] & (~[t_all.isPass2]));
i_specific2 = find((~[t_all.isPass1]) & [t_all.isPass2]);

{t_all(i_share).neuronID}
{t_all(i_specific1).neuronID}
{t_all(i_specific2).neuronID}











