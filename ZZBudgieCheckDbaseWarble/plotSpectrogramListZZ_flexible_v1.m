function [fig_count]=plotSpectrogramListZZ_flexible_v1(sound_struct, sound_field, fig_start, maxDur, plot_row, plot_col, save_folder, prefix, title_field, flim, clim, width, height, xgap, ygap, left_margin, bottom_margin, fig_pos, NFFT, windowSize, windowOverlap)
% the sound_field in sound_struct will be the audio input
% sound_struct has to have the following field: fs,
% use Brian's spectrogram function to plot a list of syllables arranged in
% several multi-panel plots
% width and height are size of each panel relative to the figure
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
    windowOverlap = windowSize-floor(0.001*sound_struct(1).fs); % each spec window slides by 1 ms
end

if ~isempty(save_folder) && ~exist(save_folder)
  mkdir(save_folder);
end

%% pad zero to the short segments
fs = sound_struct.fs;
numPoint = floor(maxDur*fs);
sound_structPad = sound_struct;
for ii=1:length(sound_structPad)
  s = sound_structPad(ii).(sound_field);
  sPad = zeros(numPoint, 1);
  idxStart = max(1, floor(numPoint/2 - length(s)/2));
  idxEnd = min(numPoint, idxStart+length(s)-1);
  sPad(idxStart:idxEnd) = s(1:(idxEnd-idxStart+1));
  sound_structPad(ii).(sound_field) = sPad;
  sound_structPad(ii).onsetTimeOfFile = 0;
  sound_structPad(ii).offsetTimeOfFile = maxDur;
end

%% loop through all syllables
plot_count = 0;
fig_count = fig_start;
for cdx=1:length(sound_structPad)
  if mod(cdx, plot_col*plot_row)==1
    if cdx>1
      if ~isempty(save_folder)
        fn_fig = fullfile(save_folder, strcat(prefix, '_', num2str(fig_count), '.fig'));
        savefig(fn_fig);
      end
      fig_count = fig_count + 1;
    end
    figure(fig_count);
    set(gcf,'Position', fig_pos);
    plot_count = 0;
  end
  plot_count = plot_count + 1;
  temp = sound_structPad(cdx);
  ax = subplot(plot_row, plot_col, plot_count);
%   ax = gca;
  showAudioSpectrogramZZ_flexible_v1(temp.(sound_field), temp.fs, gca(), flim, clim, NFFT, windowSize, windowOverlap);
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
    set(ax, 'YTickLabel',[]);
  end
  set(ax,'TickDir','out');
  % change subplot location
  plot_i = ceil(plot_count/plot_col);
  plot_j = plot_count - plot_col*(plot_i-1);
  plot_x = left_margin + width*(plot_j-1) + xgap*(plot_j-1);
  plot_y = 1 - bottom_margin - height*plot_i - ygap*plot_i;
  set(ax, 'Position', [plot_x plot_y width height]);
%   if plot_i~=plot_row
%     set(ax, 'XTickLabel',[]);
%   end
end
% save the last figure
if ~isempty(save_folder)
  fn_fig = fullfile(save_folder, strcat(prefix, '_', num2str(fig_count), '.fig'));
  savefig(fn_fig);
end
end

