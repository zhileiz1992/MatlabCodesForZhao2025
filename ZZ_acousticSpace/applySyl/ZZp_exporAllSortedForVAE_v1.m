% for a given bird, export the wav and segmentation info for applying the trained VAE model
% including the first sort batch, and later added supplementary sorting
% Zhilei, 07/31/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% 1. Folder inputs
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
% what kind of dbase to go through
% first sort batch
wsp_suffix = 'Wsp1';
% supplemenary sorting
wsp_supp = 'addMiss_*'; 

% where to save
fd_save = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', 'AudioAll1');


%% 2. export all first bach 
fns1 = dir(fullfile(fd_base, 'DbaseFiles', pairID, '*', birdID, 'warble', sprintf('%s.*.warble.good.%s.dbase.mat', birdID, wsp_suffix)));
% loop through each dbase, export production files (bProd==1)
for fi=1:size(fns1, 1)
  temp = strsplit(fns1(fi).folder, '/');
  data_date = temp{end-2};
  fd_save_this = fullfile(fd_save, 'batch1', data_date);
  if ~exist(fd_save_this, 'dir')
    mkdir(fd_save_this);
  end
  fn_d = fullfile(fns1(fi).folder, fns1(fi).name);
  load(fn_d);
  dbase.Fs = 20000;
  % which field is the bProd, should be the 3rd, but double check
  p_idx = find(strcmp(dbase.PropertyNames, 'bProd'));
  p = dbase.Properties;
  inc_idx = find(p(:,p_idx));
  fns = dbase.SoundFiles;
  % loop through sound files, export wav and segment info
  parfor sii=1:length(inc_idx)
    si = inc_idx(sii);
    % read nc, export as wav
    fn_nc = fullfile(fns(si).folder, fns(si).name);
    signal = ncread(fn_nc, 'data');
    fn_save = fullfile(fd_save_this, strrep(fns(si).name, '.nc', '.wav'));
    audiowrite(fn_save, signal, dbase.Fs);
    % export segment information
    seg_t = dbase.SegmentTimes{si};
    fn_t = strrep(fn_save, '.wav', '.time.txt');
    writematrix(seg_t, fn_t, 'Delimiter', ',');
    seg_n = dbase.SegmentTitles{si};
    fn_n = strrep(fn_save, '.wav', '.label.txt');
    fid = fopen(fn_n, 'w');
    for i = 1:numel(seg_n)
      fprintf(fid, '%s\n', seg_n{i});
    end
    fclose(fid);
  end
end


%% 3. Export supplementary sorting
fns2 = dir(fullfile(fd_base, 'DbaseFiles', pairID, '*', birdID, 'warble', sprintf('%s.*.warble.good.ch*.%s.dbase.mat', birdID, wsp_supp)));
% loop through each dbase, export production files (bProd==1)
for fi=1:size(fns2, 1)
  temp = strsplit(fns2(fi).folder, '/');
  data_date = temp{end-2};
  fd_save_this = fullfile(fd_save, 'batch2', data_date);
  if ~exist(fd_save_this, 'dir')
    mkdir(fd_save_this);
  end
  fn_d = fullfile(fns2(fi).folder, fns2(fi).name);
  load(fn_d);
  dbase.Fs = 20000;
  % select what sound files to use: choose bSorted==1 && bProd==0 && bM_vX==1
  p = dbase.Properties;
  if size(p, 2)<7
    fprintf('Not sorted: %s\n', fn_d);
  end
  inc_idx = find(p(:,1) & (~p(:,3)) & p(:,7));
  fns = dbase.SoundFiles;
  % loop through sound files, export wav and segment info
  parfor sii=1:length(inc_idx)
    si = inc_idx(sii);
    % read nc, export as wav
    fn_nc = fullfile(fns(si).folder, fns(si).name);
    signal = ncread(fn_nc, 'data');
    fn_save = fullfile(fd_save_this, strrep(fns(si).name, '.nc', '.wav'));
    audiowrite(fn_save, signal, dbase.Fs);
    % export segment information
    seg_t = dbase.SegmentTimes{si};
    fn_t = strrep(fn_save, '.wav', '.time.txt');
    writematrix(seg_t, fn_t, 'Delimiter', ',');
    seg_n = dbase.SegmentTitles{si};
    fn_n = strrep(fn_save, '.wav', '.label.txt');
    fid = fopen(fn_n, 'w');
    for i = 1:numel(seg_n)
      fprintf(fid, '%s\n', seg_n{i});
    end
    fclose(fid);
  end
end





















