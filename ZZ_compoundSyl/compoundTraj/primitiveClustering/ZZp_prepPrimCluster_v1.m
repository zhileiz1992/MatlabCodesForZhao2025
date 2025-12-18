% preprocessing to prepare for VAE training and UMAP/HDBSCAN clustering
% Zhilei, 09/21/2025
% from the onsets/offsts of identified primitives, export to .wav, segmentation .label.txt and .time.txt files

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


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Load data on identified primitives
fd_save_base = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre');
fn_prim = fullfile(fd_save_base, 'ref_tar_res', sprintf('%s.prim.mat', birdID));
load(fn_prim);
% load information regarding compound syllables
fn_comp = fullfile(fd_save_base, sprintf('%s.comp.mat', birdID));
load(fn_comp);

% where to save results
fd_vae = fullfile(fd_save_base, 'prim_VAE');
fd_vae_wav = fullfile(fd_vae, 'audio');
if exist(fd_vae_wav, 'dir')
  rmdir(fd_vae_wav, 's');
end
mkdir(fd_vae_wav);


%% Loop through identified primitives, export wav and segmentation
fs = 20000;
for si=1:size(prim,2)
% for si=1:10
  if ~isempty(prim(si).onsets)
    fn_wav = fullfile(fd_vae_wav, sprintf('%s.ref%d.wav', birdID, si));
    sound = prim(si).sound;
    audiowrite(fn_wav, sound, fs);
    % export segmentation info
    seg_t = [prim(si).onsets prim(si).offsets];  
    % convert to data points; originally in ms
    seg_t = floor(seg_t / 1000 * fs);
    fn_t = strrep(fn_wav, '.wav', '.time.txt');
    writematrix(seg_t, fn_t, 'Delimiter', ',');
    
    fn_n = strrep(fn_wav, '.wav', '.label.txt');
    fid = fopen(fn_n, 'w');
    for i = 1:numel(prim(si).onsets)
      fprintf(fid, 'c\n');
    end
    fclose(fid);
    
%     close all;
%     [fig, axes] = generatePanelGrid_v2(1, 10, [0.5], [], [0.05;0.05], [0.05;0.05], 0.01, [0], [10 10 1000 400]);
%     for ii=1:size(seg_t, 1)
%       ax = axes(ii);
%       [ax, ~, ~,  ~, ~,  ~] = showAudioSpectrogramZZ_flexible_v1(sound(seg_t(ii,1):seg_t(ii,2)), fs, ax, [250 7500], [12 23], 256, 256, 236); 
%     end
%     title(axes(1), sprintf('ref_i %d', si));
  end
end











