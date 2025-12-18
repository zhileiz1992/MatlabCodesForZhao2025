% plot primitive example to use in figures
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
fd_save = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
% if ~exist(fd_save, 'dir'); mkdir(fd_save); end
fs = 20000;
% load the info and VAE data for compound syllable
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
load(fn_comp);

% where the cross-distance matrix is located
fd_mat = fullfile(fd_save, 'ref_tar_loop');
% where to save result
fd_save_res = fullfile(fd_save, 'ref_tar_res');
fd_save_plot = fullfile(fd_save_res, 'example');
if ~exist(fd_save_plot, 'dir'); mkdir(fd_save_plot); end

% data for identified peaks
fn_prim = fullfile(fd_save_res, sprintf('%s.prim.mat', birdID));
load(fn_prim);


%% 2. Replot the stacked plot
% re-adjust the relative size of each panel
ref_i = 103;

fs = 20000; 
sound = prim(ref_i).sound;
xall = prim(ref_i).xall;
x_prob = prim(ref_i).x_prob;
amp = prim(ref_i).amp;
prob_amp = prim(ref_i).prob_amp;
onsets = prim(ref_i).onsets;
offsets = prim(ref_i).offsets;

close all;
[fig3, axes3] = generatePanelGrid_v2(5, 1, [0.2;0.2;0.1;0.1;0.1], [0.01;0.01;0.01;0.01], [0.05;0.05], [0.12;0.05], 0.05, [0;0;0;0;0], [10 10 700 900]);
ax1 = axes3(1);
[ax1, ~, ~,  ~, ~,  t1] = showAudioSpectrogramZZ_flexible_v1(sound, fs, ax1, [250 7500], [12 23], 256, 256, 236);
title(ax1, sprintf('Compound syl. %d', ref_i), 'FontSize', 14);
xticks(ax1, []);
yticks(ax1, [1000 3000 5000 7000]);
yticklabels(ax1, {'1k', '3k', '5k', '7k'});

% then plot where the strands are
ax2 = axes3(2); cla(ax2); hold(ax2, 'on'); 
for row_i=1:size(xall,1)
  x_strand = find(xall(row_i,:)==1);
  plot(ax2, [x_strand(1) x_strand(end)], [row_i row_i], 'LineStyle', '-', 'LineWidth', 1, 'Color', 'black');
end
xlim(ax2, [1 size(xall,2)]);
ylim(ax2, [1 size(xall,1)]);
ylabel(ax2, 'Strand ID', 'FontSize', 12); 
% xlabel(ax2, 'Rel. time (ms)', 'FontSize', 12);
xticks(ax2, []);

% then calculate and plot probability
ax3 = axes3(3); 
plot(ax3, x_prob); 
xlim(ax3, [1 length(x_prob)]);
ylabel(ax3, 'Density', 'FontSize', 12);     
% xlabel(ax3, 'Rel. time (ms)', 'FontSize', 12);
xticks(ax3, []);

% plot the sound amplitude
ax4 = axes3(4); cla(ax4); hold(ax4, 'on');
plot(ax4, amp);
yline(ax4, -30, 'LineStyle', '--', 'LineWidth', 1, 'Color', '#737373');
xlim(ax4, [1 length(amp)]);
ylabel(ax4, 'Amplitude', 'FontSize', 12);
xticks(ax4, []);
ylim(ax4, [-35 1]);

% use binarized amplitude to seg the probablity curve
ax5 = axes3(5); cla(ax5); hold(ax5, 'on');
plot(ax5, prob_amp);
xlim(ax5, [1 length(prob_amp)]);
ylabel(ax5, 'Prob * Amp', 'FontSize', 12);
xlabel(ax5, 'Rel. time (ms)', 'FontSize', 12);

% indicate where the peaks are
if ~isempty(onsets)
  for ii=1:length(onsets)
      xline(ax5, onsets(ii), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 2);
      xline(ax5, offsets(ii), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
      % also plot in the spectrogram
      xline(ax1, onsets(ii)/length(prob_amp), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1);
      xline(ax1, offsets(ii)/length(prob_amp), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
  end
end

for ai=1:size(axes3,1)
  axes3(ai).Box = 'off';
end

% save figure
fn_fig = fullfile(fd_save_plot, sprintf('ProbHigh.ref%d.pdf', ref_i));
print(fig3, fn_fig, '-dpdf', '-painters');



%% 3. Plot example cross-distance matrix
ref_i = 103;
fd_save_distm = fullfile(fd_save_plot, 'distm_plot2');
if ~exist(fd_save_distm, 'dir'); mkdir(fd_save_distm); end

to_loop = 1000;
i_all = 1:size(comp,1);
rng(ref_i);
rd_list = randsample(i_all(i_all~=ref_i), to_loop);
pass_all = [];
d1 = comp.dvae{ref_i};
d2s = comp.dvae(rd_list);
% for rdi=1:length(rd_list)
for rdi=1:200
%   rdi = 11;
  idx_rd = [ref_i rd_list(rdi)];
  close all;
  % calculate then plot the cross-correlation
  [fig, axes_all, distm, t1, t2] = ZZfunc_plotCrossCorreTwoSpec_v2(comp, idx_rd);
  fn_fig = fullfile(fd_save_distm, sprintf('ref%d.rdi%d.tar%d.pdf', ref_i, rdi, idx_rd(2)));
  print(fig, fn_fig, '-dpdf', '-painters');
 
  % identify the similarity strands
%   med_filter = [11 11];
%   thre = 0.35;
%   min_dur = 25;
%   [fig2, ax_all2, pass_i, newMask] = ZZfunc_identifySimilarityStrand_v1(distm, med_filter, thre, min_dur);
%   fn_fig = fullfile(fd_save_distm, sprintf('ref%d.rdi%d.tar%d.strand.pdf', ref_i, rdi, idx_rd(2)));
%   print(fig2, fn_fig, '-dpdf', '-painters');

end


% plot for a specific target syllable with colorbar
rdi = 11;
idx_rd = [ref_i rd_list(rdi)];
% plot one with color bar
close all;
[fig, axes_all, distm, t1, t2] = ZZfunc_plotCrossCorreTwoSpec_v2(comp, idx_rd);
colorbar(axes_all(3));
fn_fig = fullfile(fd_save_distm, sprintf('ref%d.rdi%d.tar%d.colorbar.pdf', ref_i, rdi, idx_rd(2)));
print(fig, fn_fig, '-dpdf', '-painters');



%% 4. Plot the location of primitives in target compound syllables
% x axis plot the latency of the similarity strand on the reference syllable
% y axis plot the latency of the similarity strand on the target syllable
ref_i = 103;

fs = 20000; 
sound = prim(ref_i).sound;
xall = prim(ref_i).xall;
x_prob = prim(ref_i).x_prob;
amp = prim(ref_i).amp;
prob_amp = prim(ref_i).prob_amp;
onsets = prim(ref_i).onsets;
offsets = prim(ref_i).offsets;

% also load the raw location of the similarity strands
batch_size = 20;
batch_i = ceil(ref_i / batch_size);
fn_pass = fullfile(fd_save, 'ref_tar_loop', sprintf('batch%d.pass_list.mat', batch_i));
a = load(fn_pass); pass_list = a.pass_list;
i_loc = find([pass_list.ref_i]==ref_i);
pass_all = pass_list(i_loc).pass_all;

% loop over identified primitives
% for each primitive, loop through the high similarity strands, determine if it belongs to the primitive or not
% if so, record the center location of x(reference) and y(target)
% get the xy location of all similarity strands for the given reference syllable
strands = [];
count = 0;
for si=1:size(pass_all, 2)
    p = pass_all(si);
    pass_i = p.pass_i;
    for strand_i=1:size(pass_i,2)
      [yy, xx] = ind2sub(p.size_d, pass_i{strand_i});
      x_range = [min(xx) max(xx)];
      % add the half-width of sliding windows back, so considering the center of sliding window
      x_range = x_range + win_size/2;
      y_range = [min(yy) max(yy)];
      % add the half-width of sliding windows back, so considering the center of sliding window
      y_range = y_range + win_size/2;
      count = count + 1;
      strands(count).si = si;
      strands(count).tar_i = p.tar_i;
      strands(count).strand_i = strand_i;
      strands(count).x_range = x_range;
      strands(count).y_range = y_range;
      strands(count).x_center = mean(x_range);
      strands(count).y_center = mean(y_range);
    end   
end

% loop throgh strands, determine if it belongs to any identified primitives
% require: center within the boundary; overlap larger than a ratio
thre_overlap = 0;
for strand_i=1:size(strands,2)
  % loop through identified primitives
  prim_belong = nan;
  for prim_i=1:length(onsets)
    % calculate overlap
    a = [onsets(prim_i) offsets(prim_i)];
    b = strands(strand_i).x_range;
    x_center = strands(strand_i).x_center;
    [overlapRatio, intersectionLength] = calculateOverlapRatio(a, b);
    if (intersectionLength>=thre_overlap) && (x_center>=a(1)) && (x_center<=a(2))
%     if (x_center>=a(1)) && (x_center<=a(2))
      prim_belong = prim_i;
      break
    end
  end
  strands(strand_i).prim_belong = prim_belong;
end

col_prim = lines(length(onsets));
close all; 
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax=gca; hold(ax, 'on');
for prim_i=1:length(onsets)
  % prim_i = 2;
  strands_plot = strands([strands.prim_belong]==prim_i);
  scatter(ax, [strands_plot.x_center], [strands_plot.y_center], 20, 'Marker', 'o', 'MarkerFaceColor', col_prim(prim_i,:), 'MarkerEdgeColor', 'none');
end
xlim(ax, [1 p.size_d(2)]);
ylim(ax, [1 1000]);
xlabel('Latency of primitives in refererence syl. (ms)', 'FontSize', 12);
ylabel('Latency of primitives in other syl. (ms)', 'FontSize', 12);

fn_fig = fullfile(fd_save_plot, sprintf('ref%d.target_loc.pdf', ref_i));
print(fig, fn_fig, '-dpdf', '-painters');
  















