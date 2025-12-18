% calculate and plot the firing field of MO neurons in the acoustic space
% step 1: Extract ephys struct for non-V syllables
% this step has been done for all call syllables

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what WhisperSeg dbase to use
suffix = 'Wsp2Call';
% what syllables to extract
labels_extract = {'b', 'x', 'h', 'e'};


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load data
% load neuron information 
fn_info = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
load(fn_info);
% get unique dates
date_unique = unique(info.date_long);
% where to save extracted ephys struct
fd_save = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'extractedReplaced2');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end
disp(fd_save);



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
%   for di=1:size(date_unique,1)
    data_date = date_unique{di};
%     disp(data_date);
    %   fprintf('%d: %s', di, data_date);
    fd_master = fullfile(fd_base, 'DbaseFiles', pairID, data_date, birdID, 'warble');
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
      segments = ZZ_extractSyllableSpikeFromDbaseFunc_v4_simple(dbase_spike, dbase_seg, {syl_label}, ch_num, spike_shape, pad_pre, pad_post);
      % only focus on production for now
      segments = segments(strcmp({segments.aud_ch}, 'chan0'));
      if ~isempty(segments)
        [segments.data_date] = deal(data_date);  % add extra columns of date and neuron ID
        [segments.neuronID] = deal(sprintf('%s-%s', info.date{sorted_idx(ci)}, ch));
        [segments.fn_dbase] = deal(fn_spike);
        segments_this = [segments_this segments];
      end
    end
    segments_all{di} = segments_this;
%     clear segments_this;
%     clear segments;
  end
  % Remove empty elements
  nonempty_segments = segments_all(~cellfun(@isempty, segments_all));
  % Concatenate all struct arrays
%   segments_all = [nonempty_segments{:}];
  seg_selected = [nonempty_segments{:}];
  
  % save the results into disk for late use
  fn_save_seg = fullfile(fd_save, sprintf('%s.%s.segments_all.replaced2.mat', birdID, syl_label));
  save(fn_save_seg, 'seg_selected', '-v7.3');
  
  clear nonempty_segments;
  clear seg_selected;
end









