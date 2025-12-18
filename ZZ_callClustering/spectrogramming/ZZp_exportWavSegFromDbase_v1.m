% export wav files and segmentations to a temporary folder
% for convenient network training and access

clear; close all;


%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
fd_save_base = fullfile(fd_home, 'vaeWav');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
% what WhisperSeg results to use
suffix = 'Wsp1';


%% Loop through bird and dates
% bi = 1;
for bi=2:size(birdIDs,2)
  birdID = birdIDs{bi};
  pairID = pairIDs{bi};
  % load meta info
  fn_meta = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
  load(fn_meta);
  % what date have sorted neurons
  date_unique = unique(info.date, 'sorted');
  
  % Loop through each date, export wav files and segments
  for date_i=1:size(date_unique, 1)
    date_short = date_unique{date_i};
    data_date = strjoin({date_short(1:4), date_short(5:6), date_short(7:8)}, '-'); % add a dash in between
    % save in separate folder
    fd_save_this = fullfile(fd_save_base, birdID, data_date);
    if exist(fd_save_this, 'dir')
      rmdir(fd_save_this, 's');
    end
    mkdir(fd_save_this);
    % load dbase
    fn_base = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, date_short, suffix);
    fn_d = fullfile(fd_home, 'DbaseFiles', pairID, data_date, birdID, 'warble', fn_base);
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
end









