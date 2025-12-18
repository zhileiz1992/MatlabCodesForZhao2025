% a temp script to extract wav files from Intan recordings
clear; close all; 
addpath(genpath("/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZBudgieIntanExtract"));
% also add the lab electro_gui path
addpath(genpath('/mnt/z4/zz367/LabSoftware/ElectroGui-main'));

fd_base_intan = '/mnt/z4/zz367/EphysMONAO/EphysDataNew/2024Ephys1';
data_date = '2024-08-15';
fd_save = '/mnt/z4/zz367/EphysMONAO/Analyzed/tempWav';
fd_save_this = fullfile(fd_save, data_date);
if exist(fd_save_this)
  rmdir(fd_save_this, 's');
end
mkdir(fd_save_this);

%% loop through all intan files
% what audio channel to extract
channels = [1,2,3,4,5,6,7];
% 
fns_intan = dir(fullfile(fd_base_intan, data_date, '*.rhd'));
% for fi=1:3
parfor fi=1:length(fns_intan)
  fn_save = fullfile(fd_save_this, strrep(fns_intan(fi).name, 'rhd', 'wav'));
  d_intan = read_Intan_RHD2000_file_to_struct_2_noNotch(fns_intan(fi).folder, fns_intan(fi).name, 0);
  % focus on selected channels
  signalThis = d_intan.board_adc_data(channels,:);
  fs = d_intan.frequency_parameters.board_adc_sample_rate;
  signalThis = signalThis';
  audiowrite(fn_save, signalThis/10, fs);
end

