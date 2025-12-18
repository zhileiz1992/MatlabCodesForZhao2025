% when matching sequences between call subtypes, some key neurons don't have enough renditions
% may need to sort a few more files for that particular neuron-call subtype pairs
% 1) load the original Wsp dbase, then segment/annotate all audios files that were not annotated yet -> Wsp2Med2 dbase
% 2) duplicate oginal Ephys dbase, flush the new annotation in, no overwritting of already annotated 
% 3) mark sound files that have the desired the call subtype -> addMissNosort_v5 dbase


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
% what clims are used for each bird when making spectrograms for VAE training
clims = {[1.5,7], [1.5,7], [1.5,8.5], [1.5,8]};
% what max duration is used for each bird when embedding the call spectrogram in 128*128 matrix
% check the ZZp1_makeSpectrogram_Goffinet_MO_v2_getMaxDur.ipynb for calculation
max_durs = {'0.2800', '0.3250', '0.3230', '0.3300'};
% what call subtype pairs to cross match, each bird can have several pairs
% crossM = {{{'v4', 'v5'}, {'v1', 'v7'}}, ...
%           {{'v1', 'v4'}}};
crossM = {{{'v2', 'v4'}, {'v2', 'v5'}, {'v2', 'v1'}}};
% what WhisperSeg model to use
fd_ckpt = '/mnt/z4/zz367/WarbleAnalysis/Results/ModelBackup/20250514_MOsortedAll/final_checkpoint_ct2';
% what WhisperSeg dbase to use
suffix_in = 'Wsp2Call';
% an intermediate dbase to store all segmentation
suffix_med1 = 'Wsp2nonsort';
suffix_med2 = 'Wsp2Med2';
% what output suffix
suffix_out = 'addMissNosort';


% loop through birds
bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
% what VAE/UMAP/HDBSCAN models to use
fd_models = fullfile(fd_home, 'vaeWav', birdID, 'UMAPonVAE7', 'v', 'spec_goffinet_nn_256_176');
fn_vae = fullfile(fd_models, sprintf('%s_checkpoint_final.tar', birdID));
fn_umap = fullfile(fd_models, sprintf('UMAPmodel_%s.p', birdID));
fn_hdbscan = fullfile(fd_models, sprintf('HDBSCANmodel_%s.p', birdID));


%% 1. Automatically generate an info file specifying what call subtype-neuron pairs need more data
% what pairs for this bird
cM = crossM{bi};
T = table();
for ci=1:size(cM, 2)
  % calculate both ways
  for flip_i=1:2
    ref_v = cM{ci}{flip_i};
    syl_v = cM{ci}{3-flip_i};
    fprintf('ref %s <- target %s\n', ref_v, syl_v);
    % read in neuron order and number of renditions from reference syllable
    fd_hahn = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew', birdID, ref_v);
    fn_order = fullfile(fd_hahn, sprintf('Hahnloser-%s-chan0.neuron_orderedPlotted.mat', ref_v));
    load(fn_order);
    fn_ren = fullfile(fd_hahn, sprintf('Hahnloser-%s-chan0.sampled_rends_ordered.mat', ref_v));
    load(fn_ren);
    % load extracted data struct for target syllable
    fd_d = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew', birdID, 'extracted');
    fn_d = sprintf('%s.%s.segments_all.mat', birdID, syl_v);
    load(fullfile(fd_d, fn_d));
    seg_selected = segments_all(strcmp({segments_all.aud_ch}, 'chan0'));
    % go through each neuron, check if target syllable has enough renditions
    for ni=1:size(neuron_ordered, 2)
      nidx = find(strcmp({seg_selected.neuronID}, neuron_ordered{ni}));
      if length(nidx) < length(sampled_rends_ordered{ni})
        T = [T; {neuron_ordered{ni}, syl_v}];
      end
    end
  end
end
% write results to a txt file
fd_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo');
fn_info = fullfile(fd_info, sprintf('%s_addMissingCall2.txt', birdID));
writetable(T, fn_info, 'Delimiter', '\t', 'WriteVariableNames', false);
  

% load the info regarding missing neuron-call subtype pairs
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_addMissingCall2.txt', birdID));
% or load manually found neuron-call subtype pair
% fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_addManualCall2.txt', birdID));
info = readtable(fn_info, 'Delimiter', '\t', 'ReadVariableNames', false);
% get unique dates
temp = cellfun(@(x) strsplit(x, '-'), info.Var1, 'UniformOutput', false);
date_unique = cellfun(@(x) x{1}, temp, 'UniformOutput', false);
date_unique = unique(date_unique, 'stable');


%% 2. Generate complete annotation at the date level
% loop through date, generate a new segmentation dbase with all audio file labeled 
% no overwriting of the annotations that alreay exist
% this will take quite a while for running WhisperSeg over all sound files 
for di=6:size(date_unique, 1)
% for di=[2,8]
  dd = date_unique{di};
  date_long = [dd(1:4) '-' dd(5:6) '-' dd(7:8)];
  disp(date_long);
  % load the original Wsp dbase
  fd_dbase = fullfile(fd_home, 'DbaseFiles', pairID, date_long, birdID, 'warble');
  fn_seg = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, dd, suffix_in);
  load(fullfile(fd_dbase, fn_seg)); 
  dseg = dbase; dseg.Fs = 20000;
  
  % what files don't have the annotation yet and are labeled as good episodes
  % or just all non-annotated
%   p_idx = find(strcmp(dseg.PropertyNames, 'bGood'));
%   p = dseg.Properties;
%   idx_good = find(p(:, p_idx));
  idx_non = find(cellfun(@isempty, dseg.SegmentTitles));
  idx_export = idx_non;
%   idx_export = intersect(idx_good, idx_non); 
  % export audio as wav files in temp folder: cares about the default chan0 for now
  fd_temp = fullfile(fd_home, 'tempWav', birdID, date_long); 
  if exist(fd_temp, 'dir')
    rmdir(fd_temp, 's');
  end
  mkdir(fd_temp);
  fns_sound = dseg.SoundFiles;
  parfor fii=1:length(idx_non)
    si = idx_non(fii);
    fn_nc = fullfile(fns_sound(si).folder, fns_sound(si).name);
    nc_info = ncinfo(fn_nc);
    if ~isempty(nc_info.Variables)
      signal = ncread(fn_nc, 'data');
      fn_save = fullfile(fd_temp, strrep(fns_sound(si).name, '.nc', '.wav'));
      audiowrite(fn_save, signal, dseg.Fs);
    end
  end
  
  % run the latest WhisperSeg model
  seg_res = ZZfunc_runWhisperSeg_v1(fd_temp, fd_ckpt); 
  
  % flush results into an intermediate dbase
  dbase_seg = dseg;
  replace_str = '_chan0.wav';
  [dbase, dbaseOld] = ZZ_flushDbaseSegmentationFunc_v3(dbase_seg, seg_res, replace_str); 
  fn_temp = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, dd, suffix_med1); 
  fn_seg = fullfile(fd_dbase, fn_temp);
  save(fn_seg, 'dbase');
  
  % apply the trained VAE/UMAP/HSBSCAN model to get call subtype labels
  fn_dbase = fn_seg; 
  clim = clims{bi};   
  clim = sprintf('[%.3f, %.3f]', clim(1), clim(2));
  max_dur = max_durs{bi};
  syl = 'v';
  [cmd, fn_embed, embed] = ZZfunc_applyClusterModel_v1(fn_dbase, fn_vae, fn_umap, fn_hdbscan, fd_temp, clim, max_dur, syl);
  
  % flush the all subtype labels to a new dbase
  a = load(fn_dbase); 
  dbase = a.dbase;
  dbase_out = ZZfunc_callTypeToDbase_v1(dbase, embed); 
  
   % add the original annotation back for already annotated files
  idx_already = find(~cellfun(@isempty, dseg.SegmentTitles));
  dbase_out.SegmentTitles(idx_already) = dseg.SegmentTitles(idx_already);
  dbase_out.SegmentTimes(idx_already) = dseg.SegmentTimes(idx_already);
  dbase_out.SegmentIsSelected(idx_already) = dseg.SegmentIsSelected(idx_already);
  
  % save the dbase
  fn_med2 =  fullfile(fd_dbase, sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, dd, suffix_med2)); 
  dbase = dbase_out;
  save(fn_med2, 'dbase');
end


%% 3. Mark desired call subtype at the neuron level 
% each WspMed2 -> addMissing_vx dbase
not_applied = {};
for ri=7:size(info, 1)
  temp = strsplit(info.Var1{ri}, '-');
  dd = temp{1};
  date_long = [dd(1:4) '-' dd(5:6) '-' dd(7:8)];
  % what call subtype
  call_type = info.Var2{ri};
  % get channel name and number
  ch = temp{2};
  ch_num = regexp(ch, '\d+', 'match');
  ch_num = str2num(ch_num{1});
  % load the complete segmentation dbase
  fd_dbase = fullfile(fd_home, 'DbaseFiles', pairID, date_long, birdID, 'warble');
  fn_seg = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, dd, suffix_med2);
  if ~exist(fullfile(fd_dbase,fn_seg), 'file')
    not_applied = [not_applied dd];
    continue;
  end
  load(fullfile(fd_dbase, fn_seg)); dseg = dbase; 
  % load the ephys dbase
  fn_ephys = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, dd, ch);
  load(fullfile(fd_dbase, fn_ephys)); 
  
  % find sound files that don't have annotation yet (should be all of them)
%   idx_non = find(cellfun(@isempty, dbase.SegmentTitles));
  idx_non = 1:size(dbase.SegmentTitles,2);
  % write the annotations in
  dbase.SegmentTitles(idx_non) = dseg.SegmentTitles(idx_non);
  dbase.SegmentTimes(idx_non) = dseg.SegmentTimes(idx_non);
  dbase.SegmentIsSelected(idx_non) = dseg.SegmentIsSelected(idx_non);
  
  % add a new property field to mark the existence of desired call subtype
  p_name = sprintf('bM%s', call_type);
  idx_has = find(cellfun(@(x) ismember(call_type, x), dbase.SegmentTitles));
  p = dbase.Properties;
  p(:, end+1) = deal(false);
  p(idx_has, end) = deal(true);
  dbase.Properties = p;
  dbase.PropertyNames = [dbase.PropertyNames p_name];
  
  % save as new dbase
  fn_final = sprintf('%s.%s.warble.good.%s.%s_%s.dbase.mat', birdID, dd, ch, suffix_out, call_type);
  save(fullfile(fd_dbase, fn_final), 'dbase');
  
end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
