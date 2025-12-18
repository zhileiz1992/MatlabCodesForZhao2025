% calculate the acoutic distance before and after the identified primitive
% differ from v1: color pre-prim and post-prim differently

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
clims = {[10 21]; [10 21]; [10.5 21.5]; [12 23]};
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
clim = clims{bi};



%% 1. Load compound syllable data
fs = 20000;
fd_save = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
load(fn_comp);
% load information about calls
fn_call = fullfile(fd_save, sprintf('%s.call.mat', birdID));
load(fn_call); 
fn_call_mean = fullfile(fd_save, sprintf('%s.call_mean.mat', birdID));
load(fn_call_mean);
% load information about high-similarity strands
vi = 1;
fd_save_ref = fullfile(fd_save, 'compound_call2');
fn_strand = fullfile(fd_save_ref, sprintf('call_v%d.strands.mat', vi));
load(fn_strand);
% what primitive to look at
pi = 1;
strands2 = strands([strands.prim_belong]==pi);
onset = strands2(1).onsets(pi); 
offset = strands2(1).offsets(pi); 
mean_dur = strands2(1).size_d(2);

% where to save results
fd_dist = fullfile(fd_save_ref, 'distance', sprintf('v%d', vi));
if ~exist(fd_dist, 'dir'); mkdir(fd_dist); end


%% 2. Compare acoustic distance before, during and after the primitive
% walk through different time points from the similarity band onset on the compound syllable
for ri=1:size(strands2, 2)
  tar_i = strands2(ri).tar_i;
  strands2(ri).dvae = comp.dvae{tar_i};
end

% how much to calculate before and after the primitive region
pad = 50; 
t_range = -50:(offset-onset+pad);
dist_all = nan(1, length(t_range));
for ti=1:length(t_range)
  t = t_range(ti);
  d_this = [];
  for ri=1:size(strands2,2)
    t_this = strands2(ri).y_range(1) + t;
    d = strands2(ri).dvae;
    if (t_this>=1) & (t_this<=size(d, 1))
      d_this = [d_this; d(t_this,:)];
    end
  end
  % calcualte pairwise distance
 dist = ZZfunc_pairwiseDist_v1(d_this, 1000, 'cosine');
 dist_all(ti) = dist;
end
% save for later use
dist_struct.t_range = t_range; 
dist_struct.dist_all = dist_all;
fn_struct = fullfile(fd_dist, sprintf('%s.distance.v%d.p%d.mat', birdID, vi, pi));
save(fn_struct, 'dist_struct');


%%  plot results
dist_all = dist_struct.dist_all;
t_range = dist_struct.t_range;
close all; 
col_pre = '#377eb8';
col_prim = '#e7298a';
col_post = '#e6ab02';
fig = ZZfunc_newFigurePDFsize_v1([10 10 350 600]);
ax = gca; hold(ax, 'on');
plot(ax, t_range(t_range<0), dist_all(t_range<0), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_pre);
plot(ax, t_range(t_range>=0 & t_range<=(offset-onset)), dist_all(t_range>=0 & t_range<=(offset-onset)), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_prim);
plot(ax, t_range(t_range>(offset-onset)), dist_all(t_range>(offset-onset)), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_post);
xline(0, 'LineStyle', '--', 'LineWidth', 1);
xline(offset-onset, 'LineStyle', '--', 'LineWidth', 1);
xlim(ax, [t_range(1) t_range(end)]);
xlabel('Time from prim. onset (ms)', 'FontSize', 12);
ylabel('Mean cosine distance', 'FontSize', 12);
fn_fig = fullfile(fd_dist, sprintf('%s.distance2.v%d.p%d.pdf', birdID, vi, pi));
print(fig, fn_fig, '-dpdf', '-painters');


  
  
  
  
  
  
  
  
  
  





