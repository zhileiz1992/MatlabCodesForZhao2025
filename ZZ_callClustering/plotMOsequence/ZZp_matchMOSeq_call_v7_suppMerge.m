% match MO sequences between call subtypes to a reference subtype
% 09/07/2025
% save results to the Figures/pairID/HahnloserNew/birdID/callMatch2 folder
% differ from v6: make a giant supplementary figure for all pairwise comparison


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


% loop through birds
% bi = 1;
for bi=1:size(birdIDs,2)
birdID = birdIDs{bi};
pairID = pairIDs{bi};

% grab what calls have been extracted for ephys
fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
fds = dir(fullfile(fd_save_master, birdID, 'extractedReplaced2', '*.v*.segments_all.replaced2.mat'));
% get the syllable labels
syls = cellfun(@(s) strsplit(s, '.'), {fds.name}, 'UniformOutput', false);
syls = cellfun(@(x) x{2}, syls, 'UniformOutput', false);
syls = sort(syls);
% save to a new subfolder
fd_save = fullfile(fd_save_master, birdID, 'callMatch7', 'merged');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end

% read the ephys struct into RAM so no need to re-load each time
struct_all = struct();
for sii=1:size(syls,2)
  ss = syls{sii};
  fprintf('Load data for %s\n', ss);
  fn_e = fullfile(fds(sii).folder, fds(sii).name);
  a = load(fn_e);
  struct_all.(ss) = a.seg_selected;
end



%% 1. Loop through syllables
for ref_i=1:size(syls,2)
  ref_v = syls{ref_i};
  % load the neuron order and renditions plotted
  fn_order = fullfile(fd_save_master, birdID, 'popRaster2', ref_v, sprintf('Hahnloser-%s-chan0.neuron_orderedPlotted7.mat', ref_v));
  load(fn_order);
  fn_ren = fullfile(fd_save_master, birdID, 'popRaster2', ref_v, sprintf('Hahnloser-%s-chan0.sampled_rends_ordered7.mat', ref_v));
  load(fn_ren);
  % set color
  A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#737373'};
  N = length(neuron_ordered);
  neuron_color = A(mod(0:N-1, numel(A))+1);
  close all; 
  fig_pos = [10 10 1800 700];
  [fig, axes] = generatePanelGrid_v2(2, size(syls,2), [0.15;0.75], [0.02], [0.03;0.02], [0.05;0.05], 0.02, [0;0], fig_pos);
  for tar_i=1:size(syls,2)
    tar_v = syls{tar_i};
    seg_selected = struct_all.(tar_v);
    % plot 
    tick_width = 8; % size of spike ticks, default to 8 for old plots
    fig_size = [10 10 400 900];  % in unit of pixels, [10 10 500 900] for old plots
    r_seed = 1992; % randome seed to ensure reproductivity 
    sample_method = 'max'; 
    tlims = []; % use the natural time range of this syllable type
    ax_list = axes(:, tar_i);
    label_x = floor(tar_i/size(syls,2));
    [aligned_spike, aligned_sound, neuron_ordered, ax_list, sampled_rends_ordered2, xlims] = ZZfunc_alignSpikeCall_v10_match_tSameMerge(seg_selected, ax_list, neuron_ordered, neuron_color, tick_width, r_seed, sample_method, sampled_rends_ordered, tlims, label_x);
    % change the title
    title(axes(1,tar_i), sprintf('%s sorted by %s seq.', tar_v, ref_v), 'FontSize', 8);
  end
  fn_pdf = fullfile(fd_save, sprintf('%s.merged.%s.pdf', birdID, ref_v));
  print(fig, fn_pdf, '-dpdf', '-painters');
%      fn_tiff = fullfile(fd_save, sprintf('%s.merged.%s.tiff', birdID, ref_v));
%     exportgraphics(fig, fn_tiff, 'Resolution', 600);
end




end


















  
  
  
  
  
  
  
  