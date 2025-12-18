% given dbases with sorted neurons and call subtype information
% plot the Hahnloser-style plot to show MO sparse sequences
% 06/18/2025
% save results to the Figures/pairID/HahnloserNew folder
% differ from v1: use the more coarse clustering

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
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what call subtype dbase to use
suffix = 'Wsp2Call';
% suffix = 'Wsp.disVAEcall';


% loop through birds
bi = 4;
birdID = birdIDs{bi};
pairID = pairIDs{bi};

%% 1. Load information about sorted neurons and call subtypes
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
load(fn_info);
% get unique dates
date_unique = unique(info.date_long);
% where to save intermediate results and plots
fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
if ~exist(fd_save_master, 'dir')
  mkdir(fd_save_master);
end
disp(fd_save_master);

% determine what call subtypes have been identified
labels_all = {};
parfor di=1:size(date_unique,1)
  data_date = date_unique{di};
  date_short = strrep(data_date,'-','');
  fd_master = fullfile(fd_home, 'DbaseFiles', pairID, data_date, birdID, 'warble');
  % what segmentation dbase to use
  fn_seg = fullfile(fd_master, sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, date_short, suffix));
  a = load(fn_seg);
  labels = a.dbase.SegmentTitles;
  labels_all{di} = unique([labels{:}]);
end
label_unique = unique([labels_all{:}]);
% grab call subtypes
is_match = ~cellfun('isempty', regexp(label_unique, '^v(0*[1-9]\d*)$'));
% Extract matching elements
labels_extract = label_unique(is_match);


%% 2. Extract ephys data for annotated syllable (voltage trace & sorted spike)
% skip if already extracted
% this may take a few minutes, save results to disk for later use
% how many seconds to extract before and after the annotated syllable
pad_pre = 0.1;
pad_post = 0.1;
% loop through syllable types
for li=1:size(labels_extract,2)
  syl_label = labels_extract{li};
  fprintf('Extracting for %s %s...\n', birdID, syl_label);
  % loop through all sorted neurons, save results in a master data struct
  segments_all = {};
  parfor di=1:size(date_unique,1)
    data_date = date_unique{di};
    %   fprintf('%d: %s', di, data_date);
    fd_master = fullfile(fd_home, 'DbaseFiles', pairID, data_date, birdID, 'warble');
    % what segmentation dbase to use
    fn_seg = fullfile(fd_master, sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, strrep(data_date,'-',''), suffix));
    a = load(fn_seg);
    dbase_seg = a.dbase;
    % find all sorted neurons in this date
    sorted_idx = find(strcmp(info.date_long, data_date));
    segments_this = struct([]);
    for ci=1:length(sorted_idx)
      ch = info.channel{sorted_idx(ci)};
      spike_shape = str2num(info.spike_shape{sorted_idx(ci)});
      % load the dbase with sorted spikes
      fn_spike = fullfile(fd_master, sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, strrep(data_date,'-',''), ch));
      temp = load(fn_spike);
      dbase_spike = temp.dbase;
      % only examine file entries with spikes sorted
      ch_num = regexp(ch, '\d+', 'match');
      ch_num = str2num(ch_num{1});
      segments = ZZ_extractSyllableSpikeFromDbaseFunc_v4(dbase_spike, dbase_seg, {syl_label}, ch_num, spike_shape, pad_pre, pad_post);
      if ~isempty(segments)
        [segments.data_date] = deal(data_date);  % add extra columns of date and neuron ID
        [segments.neuronID] = deal(sprintf('%s-%s', info.date{sorted_idx(ci)}, ch));
        [segments.fn_dbase] = deal(fn_spike);
        segments_this = [segments_this segments];
      end
    end
    segments_all{di} = segments_this;
  end
  % Remove empty elements
  nonempty_segments = segments_all(~cellfun(@isempty, segments_all));
  % Concatenate all struct arrays
  segments_all = [nonempty_segments{:}];
  
  % save the results into disk for late use
  fd_save_this = fullfile(fd_save_master, birdID, 'extracted');
  if ~exist(fd_save_this, 'dir')
    mkdir(fd_save_this);
  end
  fn_save_seg = fullfile(fd_save_this, sprintf('%s.%s.segments_all.mat', birdID, syl_label));
  save(fn_save_seg, 'segments_all', '-v7.3');
end


%% 3. Plot Hahnloser style plot
% for a given syllable, align spikes from different neurons, plot raster
for li=1:size(labels_extract,2)
  % load previous saved data struct
  syl_label = labels_extract{li};
  fprintf('Plotting for %s %s...\n', birdID, syl_label);
  fd_extracted = fullfile(fd_save_master, birdID, 'extracted');
  fn_save_seg = fullfile(fd_extracted, sprintf('%s.%s.segments_all.mat', birdID, syl_label));
  load(fn_save_seg);
  % for production choose chan0, for auditory choose chan17
  pt = 'chan0';
  % choose specific dates OR choose all dates
  seg_selected = segments_all(strcmp({segments_all.aud_ch}, pt));
  
  % determine what neurons to plot: filter out neurons that doesn't have burst
  pad = 0.05;
  psth_thre = 0.05;
  neuron_ordered = ZZfunc_identifyFiringNeuron_v3(seg_selected, pad, psth_thre);  % peak of PSTH larger than 0.05
  % save the neuron_order
  fd_save_this = fullfile(fd_save_master, birdID, syl_label);
  if ~exist(fd_save_this, 'dir')
    mkdir(fd_save_this);
  end
  fn_neu_order = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.neuron_ordered.mat', syl_label, pt));
  save(fn_neu_order, 'neuron_ordered');
  % subset the ephys struct to only include filtered neurons
  seg_selected = seg_selected(ismember({seg_selected.neuronID}, neuron_ordered));
  
  % plot Hahnloser stype plot
  % color neurons differently, either unique color for each neuron, or loop
  A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#737373'};
  N = length(neuron_ordered);
  neuron_color = A(mod(0:N-1, numel(A))+1);
  fn_plot = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.fig', syl_label, pt));
  to_sample = 20;  % max number of renditions to sample, plot all renditions if -1
  pad_grey = true;  % if a neuron doesn't have to_sample number of renditions, pad empty rows
  tick_width = 8; 
  fig_size = [50 80 500 1000];  % in unit of pixels
  r_seed = 1992; 
  reorder_method = '';
  sample_method = 'max';   % how to choose what renditions to plot: random; max; 
  [aligned_spike, aligned_sound, neuron_ordered, fig] = ZZfunc_alignSpikeCall_v8_mean(seg_selected, fn_plot, neuron_ordered, neuron_color, to_sample, reorder_method, pad_grey, tick_width, fig_size, r_seed, sample_method);
  
  % also export as pdf
  fn_pdf = strrep(fn_plot, '.fig', '.pdf');
  print(gcf, fn_pdf, '-dpdf', '-painters');
  
end
  
  
  
  
  
  
  
  
  
  
  
  