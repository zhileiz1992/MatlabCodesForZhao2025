% given dbases with sorted neurons and call subtype information
% plot the Hahnloser-style plot to show MO sparse sequences
% 07/05/2025
% save results to the Figures/pairID/HahnloserNew folder
% differ from v3: remove duplicated neurons; remove neurons with fewer than 5 renditions; 
% also add option to plot using same x-axis time scale


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
to_remove = {{}, {'20240902-ch15'}, {}, {}};  % what neurons to remove due to duplication
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what call subtype dbase to use
suffix = 'Wsp2Call';


% loop through birds
bi = 1;
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
fds = dir(fullfile(fd_save_master, birdID, 'extracted', '*segments_all.mat'));


%% 2. Plot Hahnloser style plot
% for a given syllable, align spikes from different neurons, plot raster
% two types of plots: tDiff and tSame for using the same time scale or not
for li=1:size(fds,1)
  % load previous saved data struct
  temp = strsplit(fds(li).name, '.');
  syl_label = temp{2};
  fprintf('Plotting for %s %s...\n', birdID, syl_label);
  load(fullfile(fds(li).folder, fds(li).name));
  % for production choose chan0, for auditory choose chan17
  pt = 'chan0';
  % choose specific dates OR choose all dates
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
 
  % determine what neurons to plot: filter out neurons that doesn't have burst
  pad = 0.05;  % pad a short dur before syllable onset and offset to be inclusive, unit is sec
  prop_thre = 0.25;  % include neuron if more than this proportion of renditions have spikes
  psth_bin_size = 0.005;  % bin size when calculating psth, determines time resolution, 1-10ms recommended for bursty neuron
  psth_thre = 10; 
  ircc_sigma = 0.01; % kernel width when calculating IRCC, unit is sec
  ircc_bin_size = 0.01;  % temporal resolution for IRCC, unit is sec
  ircc_thre = 0; 
  min_rend = 5;  % minimal number of renditions
%   prop_near_peak = 0.20; % at least 20% of renditions have at least one spike near the peak PSTH (-10ms to 10ms)
%   [neuron_ordered, criteria] = ZZfunc_identifyFiringNeuron_v5(seg_selected, pad, prop_thre, psth_bin_size, psth_thre, ircc_sigma, ircc_bin_size, ircc_thre);
  [neuron_ordered, criteria] = ZZfunc_identifyFiringNeuron_v6(seg_selected, pad, prop_thre, psth_bin_size, psth_thre, ircc_sigma, ircc_bin_size, ircc_thre, min_rend);
%   [neuron_ordered, criteria] = ZZfunc_identifyFiringNeuron_v7(seg_selected, pad, prop_thre, psth_bin_size, psth_thre, ircc_sigma, ircc_bin_size, ircc_thre, min_rend, prop_near_peak);
  
  % save the neuron_order and criteria metrics
  fd_save_this = fullfile(fd_save_master, birdID, syl_label);
  if ~exist(fd_save_this, 'dir')
    mkdir(fd_save_this);
  end
  fn_neu_order = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.neuron_ordered.mat', syl_label, pt));
  save(fn_neu_order, 'neuron_ordered');
  fn_criteria = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.criteria.mat', syl_label, pt));
  save(fn_criteria, 'criteria');
  % subset the ephys struct to only include filtered neurons
  seg_selected = seg_selected(ismember({seg_selected.neuronID}, neuron_ordered));
  
  % plot Hahnloser stype plot: different abosolute time scale
  close all;
  % color neurons differently, either unique color for each neuron, or loop through a defined lists
  A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#737373'};
  N = length(neuron_ordered);
  neuron_color = A(mod(0:N-1, numel(A))+1);
  fn_plot = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.tDiff.fig', syl_label, pt));
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
  fn_order_new = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.neuron_orderedPlotted.mat', syl_label, pt));
  save(fn_order_new, 'neuron_ordered');
  fn_sampled_rends = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.sampled_rends_ordered.mat', syl_label, pt));
  save(fn_sampled_rends, 'sampled_rends_ordered');
  
  % plot with the same absolute time scale
  close all;
  fn_plot = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.tSame.fig', syl_label, pt));
  r_seed = 1992;
  tlims = []; % use the natural time range of this syllable type
%   tlims = [-0.1; 0.5];  % or use the range of other syllable types
  [aligned_spike, aligned_sound, neuron_ordered, fig, sampled_rends_ordered, xlims] =  ZZfunc_alignSpikeCall_v10_tSame(seg_selected, fn_plot, neuron_ordered, neuron_color, to_sample, pad_grey, tick_width, r_seed, sample_method, sampled_loc, tlims);
  % also export as pdf
  fn_pdf = strrep(fn_plot, '.fig', '.pdf');
  print(gcf, fn_pdf, '-dpdf', '-painters');  
  % save the time ranges as well
  fn_t = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.tSame.xlims.mat', syl_label, pt));
  save(fn_t, 'xlims');
  
  
end
  
  
  
  
  
  
  
  
  
  
  
  