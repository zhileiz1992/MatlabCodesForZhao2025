% match MO sequences between call subtypes to a reference subtype
% 09/10/2025
% save results to the Figures/pairID/HahnloserNew/birdID/callMatch2 folder
% differ from v6: use the renaming where low-quality clusteres were identified


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
bi = 4;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load extracted data for one call subtype, then loop over different different reference order
% grab what calls have been extracted for ephys
fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
fds = dir(fullfile(fd_save_master, birdID, 'extractedReplaced2', '*.v*.segments_all.replaced2.mat'));
% get the syllable labels
syls = cellfun(@(s) strsplit(s, '.'), {fds.name}, 'UniformOutput', false);
syls = cellfun(@(x) x{2}, syls, 'UniformOutput', false);
syls = sort(syls);

% loop through each call subtype
for vi=1:size(fds,1)
  % load the data
  fn_d = fullfile(fds(vi).folder, fds(vi).name);
  a = load(fn_d);
  syl_this = strsplit(fds(vi).name, '.');
  syl_this = syl_this{2};
  
  seg_selected = a.seg_selected;
  % loop through different reference order
  for ri=1:size(syls,2)
    ref_v = syls{ri};
    % where to save results
    fd_save = fullfile(fd_save_master, birdID, 'callMatch7', sprintf('ref_%s', ref_v));
    if ~exist(fd_save, 'dir')
      mkdir(fd_save);
    end
    % load the neuron order and renditions plotted
    fn_order = fullfile(fd_save_master, birdID, 'popRaster2', ref_v, sprintf('Hahnloser-%s-chan0.neuron_orderedPlotted7.mat', ref_v));
    load(fn_order);
    fn_ren = fullfile(fd_save_master, birdID, 'popRaster2', ref_v, sprintf('Hahnloser-%s-chan0.sampled_rends_ordered7.mat', ref_v));
    load(fn_ren);
    
    % plot the population raster with the defined neuron order and rendition numbers
    close all;
    % color neurons differently, either unique color for each neuron, or loop through a defined lists
    A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#737373'};
    N = length(neuron_ordered);
    neuron_color = A(mod(0:N-1, numel(A))+1);
    fn_plot = fullfile(fd_save, sprintf('Hahnloser-ref_%s-%s.tDiff.fig', ref_v, syl_this));
%     to_sample = 20;  % max number of renditions to sample, plot all renditions if -1
%     pad_grey = false;  % if a neuron doesn't have to_sample number of renditions, pad empty rows
    tick_width = 8; % size of spike ticks, default to 8 for old plots
    fig_size = [10 10 400 900];  % in unit of pixels, [10 10 500 900] for old plots
    r_seed = 1992; % randome seed to ensure reproductivity 
    sample_method = 'max';   % how to choose what renditions to plot: random, max; 
%     sampled_loc = true;  % whether to determine neuron location based on sampled renditions
    [aligned_spike, aligned_sound, neuron_ordered, fig, sampled_rends] = ZZfunc_alignSpikeCall_v10_match(seg_selected, fn_plot, neuron_ordered, neuron_color, tick_width, fig_size, r_seed, sample_method, sampled_rends_ordered);
    % change the title
    axes = findobj(fig, 'Type', 'axes');
    title(axes(2), sprintf('%s sorted by %s seq.', syl_this, ref_v), 'FontSize', 12);
    % also export as pdf
    fn_pdf = strrep(fn_plot, '.fig', '.pdf');
    print(gcf, fn_pdf, '-dpdf', '-painters');
    
    % plot with the same absolute time scale
    close all;
    fn_plot = fullfile(fd_save, sprintf('Hahnloser-ref_%s-%s.tSame.fig', ref_v, syl_this));
    r_seed = 1992;
    tlims = []; % use the natural time range of this syllable type
  %   tlims = [-0.1; 0.5];  % or use the range of other syllable types
    [aligned_spike, aligned_sound, neuron_ordered, fig, sampled_rends_ordered, xlims] = ZZfunc_alignSpikeCall_v10_match_tSame(seg_selected, fn_plot, neuron_ordered, neuron_color, tick_width, r_seed, sample_method, sampled_rends_ordered, tlims);
    % change the title
    axes = findobj(fig, 'Type', 'axes');
    title(axes(2), sprintf('%s sorted by %s seq.', syl_this, ref_v), 'FontSize', 12);
    % also export as pdf
    fn_pdf = strrep(fn_plot, '.fig', '.pdf');
    print(gcf, fn_pdf, '-dpdf', '-painters');  
    % save the time ranges as well
    fn_t = fullfile(fd_save, sprintf('Hahnloser-%s-%s.tSame.xlims.mat', ref_v, syl_this));
    save(fn_t, 'xlims');
    
    
  end
end



















  
  
  
  
  
  
  
  