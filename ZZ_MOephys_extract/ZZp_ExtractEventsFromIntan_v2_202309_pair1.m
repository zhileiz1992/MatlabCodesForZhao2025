% A script to extract interesting events from budgie Intan recordings
% Then convert data to NC files so spikes can be sorted in electro_gui
% Zhilei Zhao, 2024/03/28
% modified to be a pipeline; deal with several pairs at the same time to avoid readings files twice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intan recording setup for a pair:
% 1. Two budgies in a dyadic chamber
% 2. Two directional mics record audio
% 3. Videos (60fps) and audios (50k) are recorded by PyVAQ
% 4. Another copy of audio (20k) are also recorded by Intan board adc
% 5. One bird or both birds may have ephys recordings, each has 16 channels
% 6. Each ephys recording has 3-channel accelerometer data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis goals:
% 1. Focus on the audios in the Intan data, no need to deal with sync issue
% 2. Start by finding all warble episodes, by either bird
% 3. Then also find standalone contact calls outside warble
% 4. Finally find gestures without vocalizations in the accelerometer data
% (Wait for Han on this)
% 5. Save these events into different subfolders, need to avoid redundancy

close all; clear;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/'));
addpath(genpath('/home/zz367/ProjectsU/WarbleAnalysis/Jupyter/MatlabCodes/ZZ_extractWarbleFromWav'));

%% Input settings
% Folder setting
fd_z4 = '/mnt/z4';
fd_base_intan = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'EphysData');
fd_base_pyvaq = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'PyVAQData');
fd_base_save = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'ExtractedIntan');
fd_base_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles');
% Dataset setting
meta.expID = 'CCU21';
% bird ID in the order of ephys channels in the intan file, if only one bird is recorded, use '' for the other bird
meta.birdID = {{'', ''}, {'pair1RigCCU21', ''}}; 
meta.pairID = {'pairX', 'pair1CU21CUX'};
meta.cam = {'CamD', 'CamE'};  % what camera for each pair
meta.audio_ch_intan = [[2, 3]; [4, 5]];  % what audio channels for each bird in intan file, 1st AI channel of Intan is speaker copy
meta.playback_ch_intan = 1; % what channel for the contact call playback copy in intan
meta.trigger_ch_intan = 1;  % what channel for the trigger signal of pyVAQ
meta.ephys_ch_list = {{[], []}, {1:16, []}}; % what ephys channel for these two birds in intan, [] if only one bird is recorded 
meta.acc_ch_list = {{[], []}, {1:3, []}}; % what accelerometer channel for these two birds in intan
meta.light_cycle = {'10:45:00', '23:15:00'}; % what time range to check for, ignore other time
meta.fs_pyvaq = 50000;  % default pyvaq sampling rate
meta.fs_intan = 20000;  % default intan sampling rate
% what date to analyze 
data_dates = {'2023-09-11', '2023-09-12', '2023-09-13', '2023-09-14', '2023-09-15', '2023-09-16', '2023-09-17', '2023-09-18', '2023-09-19', '2023-09-20'};


%% 1. Extract warble and standalone contact calls from Intan recordings
% di = 1;
for di=1:length(data_dates)
  data_date = data_dates{di};
  % extract vocalizations first, differentiate warble and standalone calls, save as wav files
  [warbleMetaAll, otherMetaAll] = ZZ_IdentifyWarbleFromIntan_v7(fd_base_intan, fd_base_save, data_date, meta);
  % then export the ephys & accelerometer channels to nc files for each pair, for the warble first
  for pi=1:length(warbleMetaAll)
    [fds_save_nc] = ZZ_exportIntanNC_v2_pairs(warbleMetaAll, meta, pi, data_date, fd_base_save, 'warble');
  end
  % then for the calls
  for pi=1:length(otherMetaAll)
    [fds_save_nc] = ZZ_exportIntanNC_v2_pairs(otherMetaAll, meta, pi, data_date, fd_base_save, 'other');
  end
end


%% 2. Create the dbase using electro_gui as the starting point of analysis
% use the ZZnew default setting
% save as birdID.date.warble.start.dbase.mat
di = 1;
data_date = data_dates{di};
for pi=1:length(meta.pairID)
  for bi=1:length(meta.birdID{pi})
    if ~isempty(meta.birdID{pi}{bi})
      fd_dbase = fullfile(fd_base_dbase, meta.pairID{pi}, data_date, meta.birdID{pi}{bi}, 'warble');
      disp(fd_dbase);
      if ~exist(fd_dbase)
        mkdir(fd_dbase);
        % create an empty text file
        fn_txt = fullfile(fd_dbase, sprintf('%s.%s.warble.good.txt', meta.birdID{pi}{bi}, strrep(data_date, '-','')));
        fileID = fopen(fn_txt, 'w'); % Open a new text file for writing
        fclose(fileID);
        fn_mat = fullfile(fd_dbase, sprintf('%s.%s.warble.start.dbase.mat', meta.birdID{pi}{bi}, strrep(data_date, '-','')));
        dbase = struct();
        save(fn_mat, 'dbase'); 
      end
    end
  end
end


%% 3. Select good episodes produced by the focal bird for initial sorting
% go through the WarbleWav folder, open in Audacity, choose good warble episodes for the focal bird
% save id in a txt file in the dbase folder: birdID.date.warble.good.txt 
% subset the dbase
di = 1;
data_date = data_dates{di};
for pi=1:length(meta.pairID)
  for bi=1:length(meta.birdID{pi})
    if ~isempty(meta.birdID{pi}{bi})
      fd_base = fullfile(fd_base_dbase, meta.pairID{pi}, data_date, meta.birdID{pi}{bi}, 'warble');
      fn_dbase = fullfile(fd_base, sprintf('%s.%s.warble.start.dbase.mat', meta.birdID{pi}{bi}, strrep(data_date, '-','')));
      fn_select = fullfile(fd_base, sprintf('%s.%s.warble.good.txt', meta.birdID{pi}{bi}, strrep(data_date, '-','')));
      % subset
      if exist(fn_dbase) && exist(fn_select)
        [dbase_subset, dbase_old] = ZZ_selectDbaseFocalWarbleFunc_v2(fn_dbase, fn_select);
        fn_save = strrep(fn_dbase, 'start', 'good');
        dbase = dbase_subset;
        save(fn_save, 'dbase');
      end
    end
  end
end


%% 4. Select good episodes produced by the partner bird for initial sorting
% go through the WarbleWav folder, open in Audacity, choose good warble episodes for the focal bird
% save id in a txt file in the dbase folder: birdID.date.warble.good.txt 
% subset the dbase
di = 1;
data_date = data_dates{di};
for pi=1:length(meta.pairID)
  for bi=1:length(meta.birdID{pi})
    if ~isempty(meta.birdID{pi}{bi})
      fd_base = fullfile(fd_base_dbase, meta.pairID{pi}, data_date, meta.birdID{pi}{bi}, 'warble');
      fn_dbase = fullfile(fd_base, sprintf('%s.%s.warble.start.dbase.mat', meta.birdID{pi}{bi}, strrep(data_date, '-','')));
      fn_select = fullfile(fd_base, sprintf('%s.%s.warble.partner.txt', meta.birdID{pi}{bi}, strrep(data_date, '-','')));
      % subset
      if exist(fn_dbase) && exist(fn_select)
        [dbase_subset, dbase_old] = ZZ_selectDbaseFocalWarbleFunc_v2(fn_dbase, fn_select);
        fn_save = strrep(fn_dbase, 'start', 'partner');
        dbase = dbase_subset;
        save(fn_save, 'dbase');
      end
    end
  end
end
    




