% given a ephsy data struct, go through each spike, assign a sliding window index to it
% Zhilei, 09/08/2025
% differ from v1: use the replaced call subtype names

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedReplaced2');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.replaced2.mat', birdID)));
% what's the window size
win_frame = 32;
ms_per_frame = 1; 
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
fd_save =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
% what syllable types to analyze
syls = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
% what's the range to include spikes: default only consider spikes within 1-window of syllable boundaries
epad_frame = win_frame;


% loop through each syllables
% syl_i = 1;
for syl_i=2:size(syls,2)
  disp(syl_i);


%% 1. Load data
ss = syls{syl_i};
fn_e = fullfile(fd_ephys, sprintf('%s.%s.segments_all.replaced2.mat', birdID, ss));
load(fn_e);
fn_vae = fullfile(fd_save, sprintf('%s.%s.vae_latents.replaced2.mat', birdID, ss));
load(fn_vae);


%% 2. Loop through each rendition, then loop spike
% how to much to include outside syllable boundaries to consider
fs = 20000; 
epad_pt = floor(epad_frame*ms_per_frame/1000*fs); 
% how much extra premotor lag to add
lag_frame = 0;
% save results to a new struct
spike_embed = [];
for ri=1:size(seg_selected,2)
  seg = seg_selected(ri);
  spike_embed(ri).syl_ID = latent(ri).syl_ID;
  spike_embed(ri).neuronID = seg.neuronID;
  % find where the spikes are
  spike_i = find(seg.spike_iv);
  if ~isempty(spike_i)
    % determine the relative start/ends
    spike_embed(ri).r_start = seg.seg_start_ori - seg.seg_start+1;
    spike_embed(ri).r_end = seg.seg_end_ori - seg.seg_start+1;
    spike_embed(ri).r_stop = seg.seg_end - seg.seg_start+1;
     % boundaries where spikes will be included
    spike_embed(ri).e_start = max([1 spike_embed(ri).r_start-epad_pt]);
    spike_embed(ri).e_end = min([spike_embed(ri).r_stop spike_embed(ri).r_end+epad_pt]);
    % what spikes are inside
    spike_in = spike_i((spike_i>=spike_embed(ri).e_start) & (spike_i<=spike_embed(ri).e_end));
    if ~isempty(spike_in)
      % locate the vae/umap data
      latent_s = latent(ri);
      % calculate the relative start of each sliding window, origin is syllable onset
      slide_start = -spec_frame + latent_s.vae_meta.i_s_first; 
      slide_end = -spec_frame + latent_s.vae_meta.i_s_last; 
      slide_pos = slide_start:(latent_s.vae_meta.hop_frame):slide_end;
      
      % loop through each spike,find the sliding window
      spike_embed(ri).spike_in = spike_in;
      spike_embed(ri).sp_frame = zeros(size(spike_in));
      spike_embed(ri).win_loc = zeros(size(spike_in));
      spike_embed(ri).mat_loc = zeros(size(spike_in));
      for sp_i=1:size(spike_in,1)
        sp = spike_in(sp_i);
        % convert position to spectrogram frames, relative to the syllable onset
        sp_frame = (sp - spike_embed(ri).r_start) / fs * 1000 / ms_per_frame;
        % choose the window that just passed the spike
        win_loc = ceil(sp_frame + lag_frame); 
        spike_embed(ri).sp_frame(sp_i) = sp_frame;
        spike_embed(ri).win_loc(sp_i) = win_loc;
        % what's the index of this window in the vae/umap matrix
        matrix_loc = find(slide_pos==win_loc);
        if isempty(matrix_loc)  % outside the range of sliding windows, assign to first or last window
          if win_loc<slide_pos(1); matrix_loc=1; end
          if win_loc>slide_pos(end); matrix_loc=length(slide_pos); end
        end
        spike_embed(ri).mat_loc(sp_i) = matrix_loc;
      end
    end
  end
end

% save results
fn_res = fullfile(fd_save, sprintf('%s.%s.sliding_loc.replaced2.mat', birdID, ss));
save(fn_res, 'spike_embed');
    
end
    
%     close all; figure; 
%     seg = seg_selected(ri);
%     subplot(2, 1, 1);
%     plot(seg.signalFiltered); hold on;
%     xline(spike_embed(ri).r_start, '--', 'Color', 'green');
%     xline(spike_embed(ri).r_end, '--', 'Color', 'red');
%     xline(spike_embed(ri).e_start, '-', 'Color', 'green');
%     xline(spike_embed(ri).e_end, '-', 'Color', 'red');
%     xlim([0 spike_embed(ri).r_stop]);
%     subplot(2, 1, 2);
%     plot(seg.spike_iv); hold on;
%     xlim([0 spike_embed(ri).r_stop]);












