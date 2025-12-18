% A pipeline script to examine the properties of sorted neurons
% Zhilei Zhao, 2025-02-12
% Input: list of dbases with one channel sorted in each dbase file
% dbase is constructed by the ZZ_modern4 of electro_gui
% already run through the WhisperSeg and VAE/UMAP/HDBSCAN
% 1. Extract sorted spikes for annotated call types
% 2. Alignment to a given call
% 3. Plot aligned spikes
% Differ from v1: apply to a list of sorted dates


close all; clear;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase')); 
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils')); 

%% 0. Inputs
fd_z4 = '/mnt/z4'; 
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed'); 
pairID = 'pair5CU29CU55'; 
birdID = 'pair5RigCCU29'; 
partnerID = 'pair5Rig0CU25'; 
num_v_cluster = 12;  % how many v clusters for this bird
% what dates and neurons to analyze
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', 'pair5RigCCU29_sparseInfo.mat');
load(fn_info);
% get unique dates
date_unique = unique(info.date_long);
% where to save intermediate results and plots
fd_save_master = fullfile(fd_home, 'Figures', pairID);


%% 2. Flush the Discrete-VAE results into a new dbase to get call subtypes
% only need to run once for each sorted date
% already run the jupyter notebook: ProjectsU/EphysMONAO/Jupyter/EphysVAE/ZZp_EphysDiscreteVAEapply_v1_pair5CU29CU55_loop.ipynb
% then run codes below to flush the new annotations into a new dbase named **.Wsp.disVAEcall.dbase
% di = 1;
for di=1:size(date_unique, 1)
  data_date = date_unique{di};
  fd_master = fullfile(fd_home, 'DbaseFiles', pairID, data_date, birdID, 'warble');
  fn_seg_old = fullfile(fd_master, sprintf('%s.%s.warble.good.Wsp.dbase.mat', birdID, strrep(data_date,'-','')));
  fn_seg_new = fullfile(fd_master, sprintf('%s.%s.warble.good.Wsp.disVAEcall.dbase.mat', birdID, strrep(data_date,'-','')));
  % what files contain the sub-segment info
  fns_VAE = dir(fullfile(fd_home, 'Figures', pairID, 'VAE', 'ApplyRes', data_date, sprintf('VAE*%s*.cluster.csv', birdID)));
  load(fn_seg_old);
  dbase_old = dbase; dbase_new = dbase;
  sound_files = {dbase_old.SoundFiles(:).name};
  seg_title_new = dbase_new.SegmentTitles;
  % loop through each VAE sub-segment file, flush production and auditory into the same dbase
  for fi=1:length(fns_VAE)
    res_pd = readtable(fullfile(fns_VAE(fi).folder, fns_VAE(fi).name), 'Delimiter', ',');
    % loop through each entry in the table
    for ei=1:size(res_pd,1)
      [a,b,c] = fileparts(res_pd.fn_h5{ei});
      info_temp = strsplit(b, '_');
      fn_sound = strjoin({info_temp{2:(end-2)}}, '_');
      seg_idx = str2num(info_temp{end-1})+1;  %matlab is 1-based
      % locate the file
      sound_idx = find(strcmp(fn_sound, sound_files));
      % replace the label
      seg_this = seg_title_new{sound_idx};
      if length(seg_this)>1
        seg_this{seg_idx} = res_pd.label{ei};
      else
        seg_this = res_pd.label{ei};
      end
      seg_title_new{sound_idx} = seg_this;
    end
  end
  dbase_new.SegmentTitles = seg_title_new;
  dbase = dbase_new;
  save(fn_seg_new, 'dbase');
end


%% 3. Extract ephys data for annotated syllable (voltage trace & sorted spike)
% this may take a few minutes, save results to disk for later use
% how many seconds to extract before and after the annotated syllables
pad_pre = 0.1;  
pad_post = 0.1;
% what syllable to extract
% syl_label = 'v5';
syl_label = 'b';
% loop through all sorted neurons, save results in a master data struct
segments_all = struct([]);
% di = 1;
for di=1:size(date_unique,1)
  data_date = date_unique{di};
  fprintf('%d: %s', di, data_date);
  fd_master = fullfile(fd_home, 'DbaseFiles', pairID, data_date, birdID, 'warble');
  % what segmentation dbase to use
  fn_seg = fullfile(fd_master, sprintf('%s.%s.warble.good.Wsp.disVAEcall.dbase.mat', birdID, strrep(data_date,'-','')));
  load(fn_seg);
  dbase_seg = dbase;
  % find all sorted neurons in this date
  sorted_idx = find(strcmp(info.date_long, data_date));
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
      segments_all = [segments_all segments];
    end
  end
end
% save the results into disk for late use
fd_save_this = fullfile(fd_save_master, 'Hahnloser', syl_label);
if ~exist(fd_save_this)
  mkdir(fd_save_this);
end
fn_save_seg = fullfile(fd_save_this, sprintf('%s.%s.segments_all.mat', birdID, syl_label));
save(fn_save_seg, 'segments_all', '-v7.3');
  

%% 4. Plot Hahnloser style plot
% for a given syllable, align spikes from different neurons, plot raster
% load previous saved data struct
syl_label = 'v11';
fd_save_this = fullfile(fd_save_master, 'Hahnloser', syl_label);
fn_save_seg = fullfile(fd_save_this, sprintf('%s.%s.segments_all.mat', birdID, syl_label));
load(fn_save_seg);
% for production choose chan0, for auditory choose chan17
pt = 'chan0';
% pt = 'chan17';
% choose specific dates OR choose all dates
% sel_dates = {'2024-09-12', '2024-09-13', '2024-09-14', '2024-09-15'};
% seg_selected = segments_all(strcmp({segments_all.aud_ch}, pt) & ismember({segments_all.data_date}, sel_dates));
seg_selected = segments_all(strcmp({segments_all.aud_ch}, pt));

% the list of neurons to plot: all neurons in the struct, or filter out neurons that doesn't have burst, or load from .mat
neuron_ordered = ZZfunc_identifyFiringNeuron_v2(seg_selected, 0.05, 0.05);  % peak of PSTH larger than 0.05
% save the neuron_order
fn_neu_order = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.neuron_ordered.mat', syl_label, pt));
save(fn_neu_order, 'neuron_ordered');
% subset the ephys struct to only include filtered neurons
seg_selected = seg_selected(ismember({seg_selected.neuronID}, neuron_ordered));
% neuron_ordered = unique({seg_selected.neuronID});  % this is order by neuron ID name, i.e. date of recording

% color neurons differently, either unique color for each neuron, or loop 
% neuron_color = distinguishable_colors(length(neuron_ordered));
% A = {'#6ED6E2','#E64AA4','#8FCA46','#FF3017','#2857AE','#FFC44B','#69C248','#544DA8'};  %Okubo color codes
A = {'#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#f781bf','#a65628','#737373'};  % ColorBrewer colors
% A = {'#4477AA', '#AA3377', '#66CCEE', '#EE6677', '#228833', '#CCBB44'}; % old color codes
N = length(neuron_ordered);  
neuron_color = A(mod(0:N-1, numel(A))+1);
fn_plot = fullfile(fd_save_this, sprintf('Hahnloser-%s-%s.fig', syl_label, pt));
% seg_selected = seg_selected([seg_selected.ch]==1);
% plot by reordering the neurons according to burst time
to_sample = 50;  % max number of renditions to sample, plot all renditions if -1
reorder_method = 'peak_fr';  % reorder neurons by burst time, methods: peak_fr; center_of_mass
pad_grey = true;  % if a neuron doesn't have to_sample number of renditions, pad empty rows
tick_width = 8; 
[aligned_spike, aligned_sound, neuron_ordered] = ZZfunc_alignSpikeCall_v6_mean(seg_selected, fn_plot, neuron_ordered, neuron_color, to_sample, reorder_method, pad_grey, tick_width);
% also export as pdf
fn_pdf = strrep(fn_plot, '.fig', '.pdf');
print(gcf, fn_pdf, '-dpdf', '-painters');










  