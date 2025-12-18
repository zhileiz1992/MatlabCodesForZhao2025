% A script to extract interesting events from budgie Intan recordings
% Then convert data to NC files so spikes can be sorted in electro_gui
% Zhilei Zhao, 2024/02/01
% for MO ephys budgie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intan recording setup:
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
% Dataset setting
expID = 'pair1';
birdID = {'pair1CU21RigB', 'pair1CU35'}; % in the order of audio channels
data_dates = {'2023-11-04', '2023-11-05', '2023-11-08', '2023-11-12'};
cam = 'CamE';  % what camera for this pair of birds
audio_ch_pyvaq = [3, 4];  % what audio channels for this pair of birds in pyvaq
audio_ch_intan = [4, 5];  % what audio channels for this pair of birds in intan, 1st AI channel of Intan is speaker copy
% bird_pos = 1;  % which bird if the focal bird, 1 means the 1st in the audio channel list
% partner_pos = 2;
playback_ch_intan = 1; % what channel for the contact call playback copy
ephys_ch_list = {1:16, []}; % what ephys channel for these two birds in intan, [] if only one bird is recorded 
acc_ch_list = {1:3, []}; % what accelerometer channel for these two birds in intan
light_cycle = {'11:00:00', '23:00:00'}; % light cycle in the room
fs_pyvaq = 50000;  % default pyvaq sampling rate
fs_intan = 20000;  % default intan sampling rate


%% Extract warble and standalone contact calls, save as NC files
% loop through each date
% for di=1:length(data_dates)
for di=[1]
  data_date = data_dates{di};
  % scan for warble and other vocalization data
  % only do once for both birds
  % save in a subfolder
  fd_save_warble = fullfile(fd_base_save, expID, data_date);
  if exist(fd_save_warble)
    rmdir(fd_save_warble, 's');
  end
  mkdir(fd_save_warble);
  % locate the warble episodes in the Intan data from either bird
  % save the rest vocalization events in 'other'
  fd_data_intan = fullfile(fd_base_intan, data_date);
  saveWavPath = fullfile(fd_save_warble, 'WarbleWav');
  [warbleMeta, otherMeta] = ZZ_IdentifyWarbleFromIntan_v4(expID, fd_data_intan, audio_ch_intan, saveWavPath);
  % save the meta information
  fn_meta_warble = fullfile(fd_save_warble, sprintf('%s_%s_warble_meta.mat', expID, data_date));
  save(fn_meta_warble, 'warbleMeta');
  fn_meta_other = fullfile(fd_save_warble, sprintf('%s_%s_other_meta.mat', expID, data_date));
  save(fn_meta_other, 'otherMeta')
  
  % loop through each bird
  for bi=1:length(birdID)
    if isempty(ephys_ch_list{bi}) % skip if not recorded for ephys
      continue
    end
    % what audio channel is current focal bird
    bird_adc_ch = audio_ch_intan(bi);
    partner_adc_ch = audio_ch_intan(audio_ch_intan~=bird_adc_ch);
    ephys_ch = ephys_ch_list{bi};
    acc_ch = acc_ch_list{bi};
    
    % save warble episodes
    prefix = 'warble';
    fd_save_nc = fullfile(fd_save_warble, birdID{bi}, 'NCFiles', prefix);
    metaInfo = warbleMeta;
    [fd_save_nc] = ZZ_exportIntanNC_v2(metaInfo, fd_save_nc, bird_adc_ch, partner_adc_ch, ephys_ch, acc_ch, prefix);
    % save other vocalizations
    prefix = 'other';
    fd_save_nc = fullfile(fd_save_warble, birdID{bi}, 'NCFiles', prefix);
    metaInfo = otherMeta;
    [fd_save_nc] = ZZ_exportIntanNC_v2(metaInfo, fd_save_nc, bird_adc_ch, partner_adc_ch, ephys_ch, acc_ch, prefix);
    
    % create dbase folder for this bird and date
    fd_base_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles');
    fd_dbase_this = fullfile(fd_base_dbase, expID, data_date, birdID{bi}, 'warble');
    if ~exist(fd_dbase_this)
      mkdir(fd_dbase_this);
    end
    fd_dbase_this = fullfile(fd_base_dbase, expID, data_date, birdID{bi}, 'other');
    if ~exist(fd_dbase_this)
      mkdir(fd_dbase_this);
    end
  end
end

%% Create the dbase using electro_gui as the starting point of analysis
% use the ZZnew default setting
% save as birdID.date.warble.start.dbase.mat








