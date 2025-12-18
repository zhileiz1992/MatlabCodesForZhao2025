% generate cross-correlation plot between pairs of compound syllables
% differ from v2: use a big loop to calculate the distance matrix and save to disk for later use

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


%% 1. Load VAE and UMAP data
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'ApplySylAll', sprintf('latents.%s', vae_run));
fn_info = fullfile(fd_vae, sprintf('%s.latents.%s.info.csv', birdID, vae_run));
info_vae = readtable(fn_info, 'Delimiter', ',');
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));
% UMAP
fd_umap_base = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'UMAPcomp');
umap_run = 'v.b.h.e.x.n-1';
syls = {'v', 'b', 'h', 'e', 'x'};
fd_umap = fullfile(fd_umap_base, umap_run);
fn_info_umap = fullfile(fd_umap, sprintf('%s.%s.info.csv', birdID, umap_run));
info_umap = readtable(fn_info_umap, 'Delimiter', ',');
fn_embed = fullfile(fd_umap, sprintf('%s.%s.embedding.csv', birdID, umap_run));
embed = readmatrix(fn_embed);

% where to save results
fd_save = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
if ~exist(fd_save, 'dir'); mkdir(fd_save); end


%% 2. Focus on compound syllables
fs = 20000;
info_vae.dur = (info_vae.iend - info_vae.istart) / fs;
dur_thre = 0.3;  % unit is sec
comp = info_vae(info_vae.dur>=dur_thre & ismember(info_vae.call_subtype, {'b', 'x'}) & ismember(info_vae.batch, {'batch1'}), :);

% grab the VAE latents
for ii=1:size(comp, 1)
  syl_ID = comp.syl_ID{ii};
  dvae = h5read(fn_vae, ['/' syl_ID]);
  comp.dvae{ii} = dvae';
end
% save for later use
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
save(fn_comp, 'comp', '-v7.3');


%% 5. Select one syllable as the reference, loop over other syllables
% then for a given position in the reference syllable, find similarity strands to all other target syllables
% overlay the location of these similarity strands onto the reference syllable
% ref_i = 8004;
% ref_i = 50;
% ref_all = [8004 10 50 3922 1000 89 350 800 555 3232];
% how many target syllables to choose
to_loop = 1000;
% do it in batches, save every batch to the disk
batch_size = 20;
num_batch = ceil(size(comp,1) / batch_size);

fd_save_ref = fullfile(fd_save, 'ref_tar_loop2');
if exist(fd_save_ref, 'dir'); rmdir(fd_save_ref, 's'); end
mkdir(fd_save_ref);

for batch_i=1:num_batch
  iis = (batch_i-1)*batch_size+1;
  iie = min([size(comp,1) batch_i*batch_size]);
  ref_list = iis:iie;
  pass_list = [];
  fprintf('Calculating for batch %d...\n', batch_i);
  for ref_ii=1:length(ref_list)
    ref_i = ref_list(ref_ii); 
    
%     tic;
    i_all = 1:size(comp,1);
    rng(ref_i);
    rd_list = randsample(i_all(i_all~=ref_i), to_loop);
    pass_all = [];
    d1 = comp.dvae{ref_i};
    d2s = comp.dvae(rd_list);
    parfor rdi=1:length(rd_list)
      %   idx_rd = [ref_i rd_list(rdi)];
      %   idx_rd = [ref_i 3922];
      %   close all;
      % calculate then plot the cross-correlation
      %   [fig, axes_all, distm, t1, t2] = ZZfunc_plotCrossCorreTwoSpec_v1(comp, idx_rd);
      %   fn_fig = fullfile(fd_save_ref, sprintf('ref%d.tar%d.pdf', ref_i, idx_rd(2)));
      %   print(fig, fn_fig, '-dpdf', '-painters');
      % no plotting
      %   d1 = comp.dvae{idx_rd(1)};
      %   d2 = comp.dvae{idx_rd(2)};
      % distm = pdist2(d1, d2, 'euclidean');
      distm = pdist2(d1, d2s{rdi}, 'cosine');
      distm = distm';
      
      % identify the similarity strands
      med_filter = [11 11];
%       thre = 0.35;
      thre = 0.2; 
      min_dur = 25;
      %   [fig2, ax_all2, pass_i, newMask] = ZZfunc_identifySimilarityStrand_v1(distm, med_filter, thre, min_dur);
      %   fn_fig = fullfile(fd_save_ref, sprintf('ref%d.tar%d.strand.pdf', ref_i, idx_rd(2)));
      %   print(fig2, fn_fig, '-dpdf', '-painters');
      % no plotting
      pass_i = ZZfunc_identifySimilarityStrand_v1_noPlot(distm, med_filter, thre, min_dur);
      
%       pass_all(rdi).ref_i = ref_i;
      pass_all(rdi).tar_i = rd_list(rdi);
%       pass_all(rdi).distm = distm;
      % don't save distm, since it's too large, only save its size
      pass_all(rdi).size_d = size(distm);
      pass_all(rdi).pass_i = pass_i;
    end
%     toc;
    pass_list(ref_ii).ref_i = ref_i;
    pass_list(ref_ii).pass_all = pass_all;
  end
  % save to disk
  fn_struct = fullfile(fd_save_ref, sprintf('batch%d.pass_list.mat', batch_i));
  save(fn_struct, 'pass_list');
  
  clear pass_list;
  clear distm; 
  clear pass_all;
end
    
    
    
    
    
    
