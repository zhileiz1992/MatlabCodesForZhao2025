% A script to extract interesting events in budgie Intan recordings
% Then convert data to NC files so spikes can be sorted in electro_gui
% Zhilei Zhao, 2023/09/21
% for MO ephys bird 2: CU25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intan recording setup:
% 1. Two budgies in a dyadic chamber
% 2. Two directional mics record audio
% 3. Videos (60fps) and audios (50k) are recorded by PyVAQ
% 4. Another copy of audio (20k) are also recorded by Intan board adc
% 5. One bird or both birds may have ephys recordings, each has 16 channels
% 6. Each ephys recording has a 3-channel accelerometer data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis goals:
% 1. Focus on the audios in the Intan data, no need to deal with sync issue
% 2. Start by finding all warble episodes, by either bird
% 3. Then find contact calls outside warble
% 4. Finally find gestures without vocalizations in the accelerometer data
% 5. Save these events into different subfolders, need to avoid redundency

close all; clear;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/'));

%% Input settings
% Folder setting
fd_z4 = '/mnt/z4';
fd_base_intan = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'EphysDataNew', '2024Ephys1');
fd_base_pyvaq = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'PyVAQDataNew', '2024PyVAQ1');
fd_base_save = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'ExtractedNC');
% Dataset setting
expID = 'pair2CU25RigB';
% data_date = '2023-09-18_test';
% data_date = '2023-10-01';
data_dates = {'2024-01-19', '2024-01-20', '2024-01-22', '2024-01-26'};
cam = 'CamD';  % what camera for this pair of birds
audio_ch_pyvaq = [1, 2];  % what audio channels for this pair of birds in pyvaq
audio_ch_intan = [2, 3];  % what audio channels for this pair of birds in intan
bird_pos = 2;  % which bird, 1 means the 1st in the audio channel list
partner_pos = 1;
playback_ch_intan = 1; % what channel for the contact call playback copy
ephys_ch = {1:16, []}; % what ephys channel for these two birds in intan
acc_ch = {1:3, []}; % what accelerometer channel for these two birds in intan
light_cycle = {'11:00:00', '23:00:00'}; % light cycle in the room
fs_pyvaq = 50000;  % default pyvaq sampling rate
fs_intan = 20000;  % default intan sampling rate


%% Extract warble and contact calls out, save as NC files
% loop through each date
% for di=1:length(data_dates)
for di=[4]
  data_date = data_dates{di};
  % save in a subfolder
  fd_save_warble = fullfile(fd_base_save, data_date);
  if exist(fd_save_warble)
    rmdir(fd_save_warble, 's');
  end
  mkdir(fd_save_warble);
  % locate the warble episodes in the Intan data
  fd_data_intan = fullfile(fd_base_intan, data_date);
  saveWavPath = fullfile(fd_save_warble, 'WarbleWav');
  [warbleMeta, otherMeta] = ZZ_IdentifyWarbleFromIntan_v3(expID, fd_data_intan, audio_ch_intan, saveWavPath);
  % save the meta information
  fn_meta_warble = fullfile(fd_save_warble, sprintf('%s_%s_warble_meta.mat', expID, data_date));
  save(fn_meta_warble, 'warbleMeta');
  fn_meta_other = fullfile(fd_save_warble, sprintf('%s_%s_other_meta.mat', expID, data_date));
  save(fn_meta_other, 'otherMeta');
  
  %% Export corresponding ephys data to nc files
  bird_adc_ch = audio_ch_intan(bird_pos);
  partner_adc_ch = audio_ch_intan(partner_pos);
  ephys_ch = 1:16;
  acc_ch = 1:3;
  
  % where to save
  prefix = 'warble';
  fd_save_nc = fullfile(fd_save_warble, 'NCFiles_warble', prefix);
  metaInfo = warbleMeta;
  [fd_save_nc] = ZZ_exportIntanNC_v1(metaInfo, fd_save_nc, bird_adc_ch, partner_adc_ch, ephys_ch, acc_ch, prefix);
  % where to save
  prefix = 'other';
  fd_save_nc = fullfile(fd_save_warble, 'NCFiles_warble', prefix);
  metaInfo = otherMeta;
  [fd_save_nc] = ZZ_exportIntanNC_v1(metaInfo, fd_save_nc, bird_adc_ch, partner_adc_ch, ephys_ch, acc_ch, prefix);
end


%% Manual: create the dbase for NC files using electro_gui
% Save as XXX.warble.start.dbase.mat or XXX.other.start.dbase.mat file
fd_base_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles');
for di=[2,3,4,5]
  fd_dbase = fullfile(fd_base_dbase, data_dates{di});
  if ~exist(fd_dbase)
    mkdir(fd_dbase);
  end
end


%% Segment the syllables, add annotations to the dbase file
% different segmentation algorithms can be used here
% the simplest will be based on spectral flatness
% but can use TweetyNet as well
% write a general function that flush segmentation into dbase

fd_base_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles');
di = 3; 
prefix = 'warble';
% load the starting dbase
fn_dbase_temp = dir(fullfile(fd_base_dbase, data_dates{di}, sprintf('*.%s.start.dbase.mat', prefix)));
fn_dbase = fullfile(fn_dbase_temp.folder, fn_dbase_temp.name);
load(fn_dbase);

% segmentation based on spectral flatness
% loop through each file, find onsets and offsets of syllables
dbase_path = ZZ_winPathToLinux_v1(dbase.PathName, 'Y:', '/mnt/z4');
fns = fullfile(dbase_path, {dbase.SoundFiles(:).name});
[onsets, offsets] = ZZ_seg_flatness_NC_v1(fns); 
% flush the segmentation into dbase
[dbase, dbaseOld] = ZZ_flushDbaseSegmentation_v1(dbase, onsets, offsets); 
fn_save_dbase = strrep(fn_dbase, 'start', 'seg');
save(fn_save_dbase, 'dbase');














