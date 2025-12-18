% generate cross-correlation plot for compound syllables comparing to self
% Zhilei, 09/22/2025


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

% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load compound syllable data
fs = 20000;
fd_save = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
load(fn_comp);



%% 5. Loop over compound syllables, calculate cross-correlation again self
fd_save_ref = fullfile(fd_save, 'self_comparison');
if exist(fd_save_ref, 'dir'); rmdir(fd_save_ref, 's'); end
mkdir(fd_save_ref);
pass_list = [];
% for ref_i=1:size(comp,1)
for ref_i=1:200
  %     tic;
  %   d1 = comp.dvae{ref_i};
  %   d2 = d1;
  % calculate then plot the cross-correlation
  idx_rd = [ref_i ref_i];
  close all;
  [fig, axes_all, distm, t1, t2] = ZZfunc_plotCrossCorreTwoSpec_v1(comp, idx_rd);
  fn_fig = fullfile(fd_save_ref, sprintf('ref%d.tar%d.pdf', ref_i, idx_rd(2)));
  print(fig, fn_fig, '-dpdf', '-painters');
  % no plotting
  %     distm = pdist2(d1, d2s{rdi}, 'cosine');
  %     distm = distm';
  
  % identify the similarity strands
  med_filter = [11 11];
  thre = 0.35;
  %     thre = 0.2;
  min_dur = 25;
  [fig2, ax_all2, pass_i, ~] = ZZfunc_identifySimilarityStrand_v1(distm, med_filter, thre, min_dur);
  fn_fig = fullfile(fd_save_ref, sprintf('ref%d.tar%d.strand.pdf', ref_i, idx_rd(2)));
  print(fig2, fn_fig, '-dpdf', '-painters');
  % no plotting
  %     pass_i = ZZfunc_identifySimilarityStrand_v1_noPlot(distm, med_filter, thre, min_dur);
  
  %       pass_all(rdi).ref_i = ref_i;
  %     toc;
  pass_list(ref_i).ref_i = ref_i;
  pass_list(ref_i).tar_i = ref_i;
  pass_list(ref_i).pass_i = pass_i;
end







