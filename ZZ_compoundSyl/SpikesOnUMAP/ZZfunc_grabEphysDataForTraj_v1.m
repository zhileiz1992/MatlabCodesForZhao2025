function [sound_spike] = ZZfunc_grabEphysDataForTraj_v1(d_plot, birdID, pairID, spike_shape, pad_sound, fs, fd_base)
% given a table that contains dbase info and segmentation info for syllable renditions
% grab the sound and spike data

% read dbase

fns_dbase = unique(d_plot.fn_dbase, 'sorted');
dbase_list = cell(size(fns_dbase, 1), 1);
%     fd_dbase = fullfile(fd_base, 'DbaseFiles', pairID, info_neu.date_long{ni}, birdID, 'warble');
for fi=1:size(fns_dbase, 1)
  fn = strsplit(fns_dbase{fi}, '/');
  fn_temp = strsplit(fn{end}, '.');
  date_str = fn_temp{2};
  date_long = [date_str(1:4) '-' date_str(5:6) '-' date_str(7:8)];
  fd_dbase = fullfile(fd_base, 'DbaseFiles', pairID, date_long, birdID, 'warble');
  a = load(fullfile(fd_dbase, fn{end}));
  dbase_list{fi} = a.dbase;
end

% read sound and spike data
sound_spike = [];
for idx=1:size(d_plot,1)
  fn = d_plot.fn_dbase{idx};
  fi = find(strcmp(fns_dbase, fn));
  dbase_spike = dbase_list{fi};
  % locate sound files
  sf = dbase_spike.SoundFiles;
  fn_sound = d_plot.fn_audio{idx};
  temp = strsplit(fn_sound, '/');
  sf_i = find(strcmp({sf.name}, temp{end}));
  istart = d_plot.seg_start_ori(idx);
  iend = d_plot.seg_end_ori(idx);
  nID = d_plot.neuronID{1};
  temp_str = strsplit(nID, '-');  
  ch_num = regexp(temp_str{2}, '\d+', 'match');
  ch_num = str2double(ch_num{1}); % Convert to numeric
%   spike_shape = str2double(info_neu.spike_shape{ni});
  % retrieve sound and spike data
  [sound, e_trace, spike_iv] = ZZfunc_getDataFromDbase_v1(dbase_spike, sf_i, istart, iend, pad_sound, fs, ch_num, spike_shape);
  sound_spike(idx).syl_ID = d_plot.syl_ID{idx};
  sound_spike(idx).sound = sound;
  sound_spike(idx).e_trace = e_trace;
  sound_spike(idx).spike_iv = spike_iv;
end



end

