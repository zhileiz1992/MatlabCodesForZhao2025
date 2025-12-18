% a script to merge adjacent segments that are too close to each other
% Zhilei Zhao, 05/13/2025

close all; clear;


%% 1. Input
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'WarbleAnalysis');
% where wav data is stored
fd_base = fullfile(fd_home, 'DataNew', '20250512_MOsorted');
birdIDs = {'pair5RigCCU29', 'pair5Rig0CU55', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
% what's the threshold to merge in unit of ms
thre = 0.02;  
fs = 20000; 
% what symbol to give for the merged segment 
sym = 'x'; 


%% 2. Loop through birds and merge
thre_pt = floor(thre*fs);
% bi = 1; 
for bi=1:size(birdIDs,2)
  fn_dbase = fullfile(fd_base, sprintf('%s.applyJulie.dbase.mat', birdIDs{bi}));
  load(fn_dbase);
  % loop through each sound file
  % si = 1;
  for si=1:size(dbase.SoundFiles,1)
    t = dbase.SegmentTimes{si};
    n = dbase.SegmentTitles{si};
    [t_merged, n_merged] = GPT_merge_segments(t, n, thre_pt);
    dbase.SegmentTimes{si} = t_merged;
    dbase.SegmentTitles{si} = n_merged';
  end
  fn_dbase = fullfile(fd_base, sprintf('%s.applyJulieMerged.dbase.mat', birdIDs{bi}));
  save(fn_dbase, 'dbase');
end