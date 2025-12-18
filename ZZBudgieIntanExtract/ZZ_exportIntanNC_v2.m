% a function to export intan data as NC files to visualize in electro_gui
% when there are > 1 birds, the channel names shouldn't be > 21
% use the modulus operation 
function [fd_save_nc] = ZZ_exportIntanNC_v2(metaInfo, fd_save_nc, bird_adc_ch, partner_adc_ch, ephys_ch, acc_ch, prefix)
%% settings
% save the audio of focal bird in chan 0, audio of other bird in chan 17
song_channel = 0;
partner_channel = 17;
acc_channel = [18 19 20];
% use Brian's writer
writer = @writeIntanNcFile;
% clear the folder if already exists
if exist(fd_save_nc, 'dir')
  rmdir(fd_save_nc, 's');
end
mkdir(fd_save_nc);

%% go through each row of the metaInfo struct
% ri = 37;
parfor ri=1:length(metaInfo)
  seg_start = metaInfo(ri).segStartRel;
  seg_end = metaInfo(ri).segEndRel;
  % read data in
  dataAll = [];
  fns_rhd = metaInfo(ri).filenames;
  for fi=1:length(fns_rhd)
    data = read_Intan_RHD2000_file_to_struct_2_ZZ(fns_rhd(fi).folder, fns_rhd(fi).name, 0);
    dataAll = [dataAll; data];
  end
  % get start time of this episode (approximate)
  fn0_rhd = fullfile(fns_rhd(1).folder, fns_rhd(1).name);
  t_stamp0 = get_rhd_filename_timestamp(fn0_rhd);
  t_stamp = t_stamp0 + metaInfo(ri).timeStartRel;
  [year0, month0, day0, hour0, minute0, second0] = get_rhd_time_components(t_stamp);
  
  % save audio data for focal bird
  % devide files into chunk if it's too long
  chunk_dur = 15;  % duration of each chunk, unit is 30 sec
  audio_all = vertcat([dataAll(:).board_adc_data]);
  ephys_all = vertcat([dataAll(:).amplifier_data]);
  aux_all = vertcat([dataAll(:).aux_input_data]);
  fs = data.frequency_parameters.board_adc_sample_rate;
  delta_t = 1/fs;
  fs_aux = data.frequency_parameters.aux_input_sample_rate;
  delta_t_aux = 1/fs_aux;
  
  for chunk_num = 1:ceil(metaInfo(ri).duration/chunk_dur)
    % chunk_num = 1;
    % audio data for focal bird
    chunk_start = seg_start + chunk_dur*fs*(chunk_num-1);
    chunk_end = min([seg_end, chunk_start+chunk_dur*fs-1]);
    song_data = struct();
    song_data.timeVector = abs([year0, month0, day0, hour0, minute0, second0]);
    song_data.metaData = sprintf('%s\t%d\t%d\t%d', fns_rhd(1).name, chunk_num, chunk_start, chunk_end);
    song_data.deltaT = delta_t;
    song_data.data = audio_all(bird_adc_ch, chunk_start:chunk_end);
    % audio data for partner bird
    partner_data = song_data;
    partner_data.data = audio_all(partner_adc_ch, chunk_start:chunk_end);
    % calculate amplitude difference
    plot_bool = false;
    [amp_diff, delta_t_amp, pos_idx, neg_idx] = ZZ_audioAmpDiff_v1(song_data.data, partner_data.data, fs, plot_bool);
    if length(pos_idx)>=length(neg_idx)
      sing_bird = 1;
    else
      sing_bird = 2;
    end
    % save to nc files
    fn_prefix = sprintf('%s_%05d_%05d_%d_%s', prefix, ri, chunk_num, sing_bird, strrep(fns_rhd(1).name,'.rhd',''));
    fn_nc = fullfile(fd_save_nc, sprintf('%s_chan%d', fn_prefix, song_channel));
    writer(fn_nc, song_data.timeVector, song_data.deltaT, song_channel, song_data.metaData, song_data.data, true);
    % save audio data for its partner
    fn_nc = fullfile(fd_save_nc, sprintf('%s_chan%d', fn_prefix, partner_channel));
    writer(fn_nc, partner_data.timeVector, partner_data.deltaT, partner_channel, partner_data.metaData, partner_data.data, true);
    
    % save ephys data
    for ei=ephys_ch
      ephys_data = song_data;
      ephys_data.data = ephys_all(ei, chunk_start:chunk_end);
      % Notch filter data for 60 Hz noise
      ephys_data.data = notch_filter(ephys_data.data, fs, data.notch_filter_frequency, 10);
      ei_rel = mod(ei, 16);
      if ei_rel==0
        ei_rel = 16;
      end
      fn_nc = fullfile(fd_save_nc, sprintf('%s_chan%d', fn_prefix, ei_rel));
      writer(fn_nc, ephys_data.timeVector, ephys_data.deltaT, ei_rel, ephys_data.metaData, ephys_data.data, true);
    end
    
    % save accelerometer data, note the sampling rate is different
    chunk_start_aux = round(fs_aux * (chunk_start - 1)/ fs) + 1;
    chunk_end_aux = round(fs_aux * (chunk_end - 1) / fs) + 1;
    for aii=1:length(acc_ch)
      ai = acc_ch(aii);
      acc_data = song_data;
      acc_data.deltaT = delta_t_aux;
      acc_data.data = aux_all(ai, chunk_start_aux:chunk_end_aux);
      fn_nc = fullfile(fd_save_nc, sprintf('%s_chan%d', fn_prefix, acc_channel(aii)));
      writer(fn_nc, acc_data.timeVector, acc_data.deltaT, acc_channel(aii), acc_data.metaData, acc_data.data, true);
    end
    
    % add a chan to calculate amplitude difference between two mics
    ampDiff_chan = 21;
    ampDiff_data = song_data;
    ampDiff_data.deltaT = delta_t_amp;
    ampDiff_data.data = amp_diff;
    fn_nc = fullfile(fd_save_nc, sprintf('%s_chan%d', fn_prefix, ampDiff_chan));
    writer(fn_nc, ampDiff_data.timeVector, ampDiff_data.deltaT, ampDiff_chan, ampDiff_data.metaData, ampDiff_data.data, true);
    
  end
end
end




