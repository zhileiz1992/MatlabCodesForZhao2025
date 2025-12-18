% from the high-similarity bands, identify acoustic primitives for each compound syllable

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/calcDistance/MMD'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZBudgieIntanExtract'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';
% size of sliding window
win_size = 32;  % in unit of frames
sec_per_frame = 0.001; 
% coloring for syllables
syl_full = {'v', 'h', 'e', 'b', 'x'};
col_full = {'#e41a1c', '#984ea3', '#4daf4a', '#377eb8', '#737373'};
col_dict = struct;
for si=1:size(syl_full,2) col_dict.(syl_full{si}) = col_full{si}; end

% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load data
% where to save results
fd_save = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre');
% if ~exist(fd_save, 'dir'); mkdir(fd_save); end
fs = 20000;
% load the info and VAE data for compound syllable
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
load(fn_comp);

% where the cross-distance matrix is located
fd_mat = fullfile(fd_save, 'ref_tar_loop');
% where to save result
fd_save_res = fullfile(fd_save, 'ref_tar_res');
fd_save_plot = fullfile(fd_save_res, 'plots');
if ~exist(fd_save_plot, 'dir'); mkdir(fd_save_plot); end


%% 2. Calculate density of similarity bands, identify primitives
% parameters for identifying peaks in the prob*amp curve
param.gapSize = 10;  
param.thresholdFlatness = 0.2;  % first threshold to identify peak
param.extendPercent = 0.4; % after thresholding extend to lower threshold to certain percentage of the peak height
param.minInterval = 0.01;  % min distance between two adjacent peaks, smaller will be merged
param.minDuration = 0.025;  % minimal duration of the identified primitives
% param.maxDuration = 0.2; 
param.maxDuration = 0.075; 
param.minProminence = 0.1;  % use in the internal findpeaks function

% pad when extracting sound
pad_pre = win_size*sec_per_frame/2; pad_pre_pt = floor(fs*pad_pre);  % pad equal amount as in VAE sliding
pad_post = win_size*sec_per_frame/2; pad_post_pt = floor(fs*pad_post);

tic;
prim = []; 
% is_plot = true;
is_plot = false;

% for batch_i=1:10
for batch_i=1:250
% batch_i = 1;

fn_mat = fullfile(fd_mat, sprintf('batch%d.pass_list.mat', batch_i));
a = load(fn_mat); pass_list=a.pass_list;
% loop through the reference syllable
for ri=1:size(pass_list,2)
  ref_i = pass_list(ri).ref_i;
  pass_all = pass_list(ri).pass_all;
  n = pass_all(1).size_d;
  m = cellfun(@(x) size(x,2), {pass_all.pass_i});
  % locate the strands on top of the reference syllable
  xall = zeros(sum(m), n(2));
  count = 0;
  for si=1:size(pass_all,2)
    p = pass_all(si);
    pass_i = p.pass_i;
    for pi=1:size(pass_i,2)
      [yy, xx] = ind2sub(p.size_d, pass_i{pi});
      x_range = [min(xx) max(xx)];
      % add the half-width of sliding windows back, so considering the center of sliding window
      x_range = x_range + win_size/2;
      count = count + 1;
      iend = min([n(2) x_range(2)]);  % ignore those that pass the boundaries
      xall(count, x_range(1):iend) = 1;
    end   
  end
  
  % calculate density
  % x_prob = sum(x_cat, 1) / size(x_cat, 1);  % divided by the total number of strands
  x_prob = sum(xall, 1) / size(pass_all,2);
  
  % calculate amplitude
  [signal, fs] = audioread(comp.fn_wav{ref_i});
  i_start = max([1 comp.istart(ref_i)-pad_pre_pt]);
  i_end = min([size(signal,1) comp.iend(ref_i)+pad_post_pt]);
  sound = signal(i_start:i_end);
  win_ms = 10; 
  cut_env = 60;
  [amp, ~] = ZZfunc_soundAmplitude_v1(sound, fs, [250 5000], win_ms, cut_env);
  
  % use binarized amplitude to seg the probablity curve
  ref_db = -30;
  amp2 = zeros(size(amp));
  amp2(amp>=ref_db) = 1;
  
  % or use the orignal amplitude curve: clipped and scaled 
%   ref_db = 40;
%   amp2 = amp + ref_db;
%   amp2(amp2<0) = 0;
%   amp2 = amp2 / ref_db; 
  
  % interpolate to the same length
  x_original = 1:length(amp2);
  x_new = linspace(1, length(amp2), length(x_prob));
  amp_interp = interp1(x_original, amp2, x_new);
  % multiply probablity with amplitude
  prob_amp = x_prob .* amp_interp;
  prob_amp = smoothdata(prob_amp, 'gaussian', 7);
  
  % find peaks in the prob*amp curve
  % [onsets, offsets] = ZZ_GetFlatnessOnsetOffset(flatness, dt, param);
  [onsets, offsets] = ZZfunc_GetProbOnsetOffset_v1(prob_amp, sec_per_frame, param);
  
  % save result
  prim(ref_i).ref_i = ref_i;
  prim(ref_i).xall = xall;
  prim(ref_i).x_prob = x_prob;
  prim(ref_i).sound = sound;
  prim(ref_i).amp = amp;
  prim(ref_i).amp_interp = amp_interp;
  prim(ref_i).x_prob = x_prob;
  prim(ref_i).prob_amp = prob_amp;
  prim(ref_i).onsets = onsets;
  prim(ref_i).offsets = offsets;
  
  
  if is_plot
    % overlay the location of the similarity strands on reference syllable
%     close all;
    [fig3, axes3] = generatePanelGrid_v3(5, 1, [0.12;0.3;0.12;0.12;0.12], [0.02;0.02;0.02;0.02], [0.05;0.05], [0.12;0.05], 0.05, [0;0;0;0;0], [10 10 700 900], false);
    ax1 = axes3(1);
    [ax1, ~, ~,  ~, ~,  t1] = showAudioSpectrogramZZ_flexible_v1(sound, fs, ax1, [250 7500], [12 23], 256, 256, 236);
    title(ax1, sprintf('Compound syl. %d', ref_i), 'FontSize', 14);
    xticks(ax1, []);

    % then plot where the strands are
    ax2 = axes3(2); cla(ax2); hold(ax2, 'on'); 
    for row_i=1:size(xall,1)
      x_strand = find(xall(row_i,:)==1);
      plot(ax2, [x_strand(1) x_strand(end)], [row_i row_i], 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black');
    end
    xlim(ax2, [0 size(xall,2)]);
    ylim(ax2, [0 size(xall,1)]);
    ylabel(ax2, 'Strand ID', 'FontSize', 12); 
    % xlabel(ax2, 'Rel. time (ms)', 'FontSize', 12);
    xticks(ax2, []);

    % then calculate and plot probability
    ax3 = axes3(3); 
    plot(ax3, x_prob); 
    xlim(ax3, [0 length(x_prob)]);
    ylabel(ax3, 'Probability', 'FontSize', 12);     
    % xlabel(ax3, 'Rel. time (ms)', 'FontSize', 12);
    xticks(ax3, []);

    % plot the sound amplitude
    ax4 = axes3(4);
    plot(ax4, amp);
    xlim(ax4, [1 length(amp)]);
    ylabel(ax4, 'Amplitude', 'FontSize', 12);
    xticks(ax4, []);

    % use binarized amplitude to seg the probablity curve
    ax5 = axes3(5); cla(ax5); hold(ax5, 'on');
    plot(ax5, prob_amp);
    xlim(ax5, [1 length(prob_amp)]);
    ylabel(ax5, 'Prob * Amp', 'FontSize', 12);

    % indicate where the peaks are
    if ~isempty(onsets)
      for ii=1:length(onsets)
          xline(ax5, onsets(ii), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 2);
          xline(ax5, offsets(ii), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
          % also plot in the spectrogram
          xline(ax1, onsets(ii)/length(prob_amp), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 2);
          xline(ax1, offsets(ii)/length(prob_amp), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
      end
    end

    % save figure
    fn_fig = fullfile(fd_save_plot, sprintf('ProbHigh.ref%d.pdf', ref_i));
    print(fig3, fn_fig, '-dpdf', '-painters');
    close(fig3);
  
  end
end

end
% toc;

% save results
fn_prim = fullfile(fd_save_res, sprintf('%s.prim.mat', birdID));
save(fn_prim, 'prim', '-v7.3');










