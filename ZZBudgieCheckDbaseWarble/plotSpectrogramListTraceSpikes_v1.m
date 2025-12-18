function [fig_count]=plotSpectrogramListTraceSpikes_v1(ephys_struct, sound_field, fig_start, maxDur, plot_row, plot_col, save_folder, prefix, title_field, flim, clim, width, height, xgap, ygap, left_margin, bottom_margin, fig_pos, NFFT, windowSize, windowOverlap)
% the sound_field in ephys_struct will be the audio input
% ephys_struct has to have the following field: fs,
% use Brian's spectrogram function to plot a list of syllables arranged in
% several multi-panel plots
% width and height are size of each panel relative to the figure
% Overlay sorted spikes to the spectrograms
if ~exist('width', 'var')
    width = 0.075; 
end
if ~exist('height', 'var')
    height = 0.21; 
end
if ~exist('xgap', 'var')
    xgap = 0.004; 
end
if ~exist('ygap', 'var')
    ygap = 0.06; 
end
if ~exist('left_margin', 'var')
    left_margin = 0.1; 
end
if ~exist('bottom_margin', 'var')
    bottom_margin = 0.1; 
end
if ~exist('fig_pos', 'var')
    fig_pos = [10 10 1200 800]; 
end
if ~exist('flim', 'var') || isempty(flim)
    flim = [500, 7500];
end
if ~exist('clim', 'var') || isempty(clim)
    clim = [10, 20];
end
if ~exist('NFFT', 'var')
    NFFT = 512; 
end
if ~exist('windowSize', 'var')
    windowSize = 512; 
end
if ~exist('windowOverlap', 'var')
    windowOverlap = windowSize-floor(0.001*ephys_struct(1).fs); % each spec window slides by 1 ms
end

if ~isempty(save_folder) && ~exist(save_folder)
  mkdir(save_folder);
end

%% pad zero to the short segments
fs = ephys_struct.fs;
numPoint = floor(maxDur*fs);
ephys_structPad = ephys_struct;
for ii=1:length(ephys_structPad)
  s = ephys_structPad(ii).(sound_field);
  sPad = zeros(numPoint, 1);
  idxStart = max(1, floor(numPoint/2 - length(s)/2));
  idxEnd = min(numPoint, idxStart+length(s)-1);
  sPad(idxStart:idxEnd) = s(1:(idxEnd-idxStart+1));
  ephys_structPad(ii).(sound_field) = sPad;
  ephys_structPad(ii).onsetTimeOfFile = 0;
  ephys_structPad(ii).offsetTimeOfFile = maxDur;
  % pad the spikes as well
  e = ephys_structPad(ii).spike_iv;
  ePad = zeros(numPoint, 1);
  ePad(idxStart:idxEnd) = e(1:(idxEnd-idxStart+1));
  ephys_structPad(ii).spike_iv = ePad;
  % pad the ephys trace as well
  et = ephys_structPad(ii).e_trace_FIR;
  etPad = zeros(numPoint, 1);
  etPad(idxStart:idxEnd) = et(1:(idxEnd-idxStart+1));
  ephys_structPad(ii).e_trace_FIR = etPad;
end

%% loop through all syllables
plot_count = 0;
fig_count = fig_start;
total_panels = plot_row * plot_col * 3;  % each syllable will have 3 plots: spectrogram, ephys trace, and spikes
% calculate the y axis limit for ephys trace
y_max = max(cellfun(@max, {ephys_structPad(:).e_trace_FIR}))+10;
y_min = min(cellfun(@min, {ephys_structPad(:).e_trace_FIR}))-10;
for cdx=1:length(ephys_structPad)
  % determine where this syllable belong to what panel group 
  group_count = mod(cdx, plot_row*plot_col); 
  if group_count==1
    if cdx>1
      if ~isempty(save_folder)
        fn_fig = fullfile(save_folder, strcat(prefix, '_', num2str(fig_count), '.fig'));
        savefig(fn_fig);
      end
      fig_count = fig_count + 1;
    end
    figure(fig_count);
    set(gcf,'Position', fig_pos);
  end
  % where to place the spectrograms
  ax1_count = group_count + 2*plot_col*floor((group_count-0.1)/plot_col);
  if group_count==0
    ax1_count = plot_row*plot_col*3 - plot_col - plot_col; 
  end
  % where to place the ephys trace 
  ax2_count = ax1_count + plot_col; 
  % where to place the spike
  ax3_count = ax2_count + plot_col; 
  temp = ephys_structPad(cdx);
  ax1 = subplot(plot_row*3, plot_col, ax1_count);
%   ax = gca;
  [ax, power, powerGrey, S, f, t]=showAudioSpectrogramZZ_flexible_v1(temp.(sound_field), temp.fs, ax1, flim, clim, NFFT, windowSize, windowOverlap);
  if ~isempty(title_field)
    title(temp.(title_field));
  end
  %     title(num2str(cdx));
  x_tick = [0, 1];
  xticks(x_tick);
  xticklabels(round(maxDur*x_tick, 2));
  yticks(1000:2000:7000);
  yticklabels([1:2:7]);
  if mod(plot_count, plot_col)~=1
    set(ax1, 'YTickLabel',[]);
  end
  set(ax1,'TickDir','out');
%   % change subplot location
%   plot_i = ceil(group_count/plot_col);
%   plot_j = group_count - plot_col*(plot_i-1);
%   plot_x = left_margin + width*(plot_j-1) + xgap*(plot_j-1);
%   plot_y = 1 - bottom_margin - height*plot_i - ygap*plot_i;
%   set(ax1, 'Xtick', []);
%   set(ax1, 'Position', [plot_x plot_y width heights(1)]);
  % plot the ephys trace
  ax2 = subplot(plot_row*3, plot_col, ax2_count);
  plot(ax2, temp.e_trace_FIR, 'black'); 
  xlim(ax2, [0 length(temp.e_trace_FIR)]);
  ylim(ax2, [y_min y_max]);
  axis off; 
  % overlay spikes
  % need to convert between data points and spectrogram frames
  ax3 = subplot(plot_row*3, plot_col, ax3_count);
  spike_idx = find(temp.spike_iv==1);
  if ~isempty(spike_idx)
    spike_frame = spike_idx / length(temp.spike_iv);
    xline(ax3, spike_frame, 'black');
    % add a small text to show how many spikes are in the bursts
%     text(ax3, spike_frame(end), 7100, num2str(length(spike_frame)), 'Color', 'b', 'FontSize', 12);
  end
  xlim(ax3, [0 1]);
  axis off; 

end
% save the last figure
if ~isempty(save_folder)
  fn_fig = fullfile(save_folder, strcat(prefix, '_', num2str(fig_count), '.fig'));
  savefig(fn_fig);
end
end

