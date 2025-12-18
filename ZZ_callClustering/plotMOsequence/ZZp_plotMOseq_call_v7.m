% given dbases with sorted neurons and call subtype information
% plot the Hahnloser-style plot to show MO sparse sequences
% 09/10/2025
% save results to the Figures/pairID/HahnloserNew/birdID/popRaster folder
% differ from v6: use the replaced names where low-quality clusters were identified


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



% loop through birds
% bi = 2;
for bi=2:size(birdIDs,2)
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load information about sorted neurons and call subtypes
fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
load(fn_info);
% get unique dates
date_unique = unique(info.date_long);
% where to save intermediate results and plots
fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
% grab what calls have been extracted
fds = dir(fullfile(fd_save_master, birdID, 'extractedReplaced2', '*segments_all.replaced2.mat'));
% where to save results
fd_save = fullfile(fd_save_master, birdID, 'popRaster2');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end


%% 2. Plot Hahnloser style plot
% for a given syllable, first determine what neurons can be classified as sparse neurons 
% then align spikes from these neurons, plot population raster
% two types of plots: tDiff and tSame for using the same time scale or not
for li=1:size(fds,1)
  % load previous saved data struct
  temp = strsplit(fds(li).name, '.');
  syl_label = temp{2};
  fprintf('Plotting for %s %s...\n', birdID, syl_label);
  a = load(fullfile(fds(li).folder, fds(li).name));
  seg_selected = a.seg_selected;
  pt = 'chan0';  % focus on production
  seg_selected = seg_selected(strcmp({seg_selected.aud_ch}, pt)); 
  clear a; 
  
  % save results in a subfolder
  fd_save_this = fullfile(fd_save, syl_label);
  if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end
  
  
  % calculating PSTH and firing metrics for each neuron, same parameters were used when comparing to broad neurons
  param.align_pad = 0.1;  % how much to look before and after the syllable when align renditions, unit is seconds
  param.fs = 20000;
  param.calc_win = [0.05 0.05];  % how much to include before syllable onset and after syllable offset when calculating metrics
  % for calculating PSTH
  param.psth_bin_size = 0.005;  % bin_size for calculating psth, unit is seconds
  % for IRCC
  param.ircc_sigma = 0.01; % kernel width when calculating IRCC, unit is sec
  param.ircc_bin_size = 0.01;  % temporal resolution for IRCC, unit is sec
  % for sparseness
  param.sparse_bin_size = 0.01;
  % how much to look when count number of firing renditions near PSTH peak
  param.count_win = 0.015; 
  
  criteria = [];
  for ni=1:size(info,1)
    seg_this = seg_selected(strcmp({seg_selected.neuronID}, info.neuronID{ni}));
    metrics = struct();
    if ~isempty(seg_this)
      % where to save PSTH
      fn_pdf = fullfile(fd_save_this, sprintf('PSTH.%s.%s.%s.pdf', birdID, syl_label, info.neuronID{ni}));
%       fn_pdf = '';
      metrics = ZZfunc_calcFiringMetrics_v1(seg_this, param, fn_pdf);
      metrics.birdID = birdID;
      criteria = [criteria; metrics];
    end
  end
    
  
  % select sparse neurons to plot according to criteria
  min_fire_rend = 7;  %minimal number of renditions that have spikes near the psth peak (+-15ms)
  prop_thre = 0.2;  % minimal fraction of all renditions that have spikes near the psth peak (+-15ms)
  psth_thre = 10;  % minimal psth peak height
  ircc_thre = 0.2;  % minimal IRCC values
  sparse_thre = 0.15; % minimal sparseness
  for ni=1:size(criteria, 1)
    criteria(ni).isPass = ((criteria(ni).num_rend_psth_fire>=min_fire_rend) & ...
                           (criteria(ni).prop_psth_fire>=prop_thre) & ...
                           (criteria(ni).peak_psth>=psth_thre) & ...
                           (criteria(ni).ircc>=ircc_thre) & ...
                           (criteria(ni).sparseness>=sparse_thre));
  end
  % then select neurons based on criteria
  passed = criteria([criteria.isPass]);
  % sort by smoothed psth max time
  [~, sort_idx] = sort([passed.psth_max_smooth_t], 'descend');
  passed = passed(sort_idx);
  neuron_ordered = {passed.neuronID};
  
  % save the criteria and metrics info
  fn_criteria = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.criteria7.mat', syl_label, pt));
  save(fn_criteria, 'criteria');
  fn_neu_order = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.neuron_ordered7.mat', syl_label, pt));
  save(fn_neu_order, 'neuron_ordered');
  
%   % compare to previous codes
%   fn_old = fullfile(fd_save_master, birdID, 'v1', sprintf('Hahnloser-%s-%s.neuron_orderedPlotted5.mat', 'v1', pt));
%   a = load(fn_old); nsold = a.neuron_ordered;
%   disp('Not in the old: ');
%   disp(setdiff(neuron_ordered, nsold));
%   disp('\n\n');
%   disp('In old but missing in current: ');
%   disp(setdiff(nsold, neuron_ordered));
  
  % subset the ephys struct to only include filtered neurons
  seg_selected = seg_selected(ismember({seg_selected.neuronID}, neuron_ordered));
  
  
  % plot Hahnloser stype plot: different abosolute time scale
  close all;
  % color neurons differently, either unique color for each neuron, or loop through a defined lists
  A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#737373'};
  N = length(neuron_ordered);
  neuron_color = A(mod(0:N-1, numel(A))+1);
  fn_plot = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.tDiff7.fig', syl_label, pt));
  to_sample = 20;  % max number of renditions to sample, plot all renditions if -1
  pad_grey = false;  % if a neuron doesn't have to_sample number of renditions, pad empty rows
  tick_width = 8; % size of spike ticks, default to 8 for old plots
  fig_size = [10 10 400 900];  % in unit of pixels, [10 10 500 900] for old plots
  r_seed = 1992; % randome seed to ensure reproductivity 
  sample_method = 'max';   % how to choose what renditions to plot: random, max; 
  sampled_loc = true;  % whether to determine neuron location based on sampled renditions
  [aligned_spike, aligned_sound, neuron_ordered, fig, sampled_rends_ordered] =  ZZfunc_alignSpikeCall_v10_mean(seg_selected, fn_plot, neuron_ordered, neuron_color, to_sample, pad_grey, tick_width, fig_size, r_seed, sample_method, sampled_loc);
  % also export as pdf
  fn_pdf = strrep(fn_plot, '.fig', '.pdf');
  print(gcf, fn_pdf, '-dpdf', '-painters');  
  
  % also save the new neuron_order and sampled rends
  fn_order_new = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.neuron_orderedPlotted7.mat', syl_label, pt));
  save(fn_order_new, 'neuron_ordered');
  fn_sampled_rends = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.sampled_rends_ordered7.mat', syl_label, pt));
  save(fn_sampled_rends, 'sampled_rends_ordered');
  
  % plot with the same absolute time scale
  close all;
  fn_plot = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.tSame7.fig', syl_label, pt));
  r_seed = 1992;
  tlims = []; % use the natural time range of this syllable type
%   tlims = [-0.1; 0.5];  % or use the range of other syllable types
  [aligned_spike, aligned_sound, neuron_ordered, fig, sampled_rends_ordered, xlims] =  ZZfunc_alignSpikeCall_v10_tSame(seg_selected, fn_plot, neuron_ordered, neuron_color, to_sample, pad_grey, tick_width, r_seed, sample_method, sampled_loc, tlims);
  % also export as pdf
  fn_pdf = strrep(fn_plot, '.fig', '.pdf');
  print(gcf, fn_pdf, '-dpdf', '-painters');  
  % save the time ranges as well
  fn_t = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.tSame.xlims7.mat', syl_label, pt));
  save(fn_t, 'xlims');
  
end

end
  
  
  
  
  
  
  
  
  
  
  
  