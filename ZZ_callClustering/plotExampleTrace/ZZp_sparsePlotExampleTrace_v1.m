% plot example ephys traces to use in figures
% Zhilei, 09/14/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% where to save
fd_save = fullfile(fd_home, 'Figures', 'CombinedAnalysis', 'SparseExampleTrace');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end


%% 1. Inputs
% how much time to plot, unit is second
dur = 10;  
fs = 20000; 
% choose a few broad neurons as examples
bis = [1];
neuronIDs = {'20240903-ch6'};
fis = [104];
ts = [1.2];


%% 2. Loop through neurons and plot
for ni=1:size(neuronIDs, 2)
  n = neuronIDs{ni};
  birdID = birdIDs{bis(ni)}; 
  pairID = pairIDs{bis(ni)}; 
  info = strsplit(n, '-');
  date_short = info{1};
  date_long = [date_short(1:4) '-' date_short(5:6) '-' date_short(7:8)];
  ch = info{2};
  ch_num = str2num(strrep(ch, 'ch', ''));
  % load dbase
  fd_e = fullfile(fd_home, 'DbaseFiles', pairID, date_long, birdID, 'warble');
  fn_e = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, date_short, ch);
  load(fullfile(fd_e, fn_e));
 
  % read the sound and ephys data
  i = fis(ni);
  fns = dbase.SoundFiles;
  fn_aud = fullfile(fns(i).folder, fns(i).name);
  aud = double(ncread(fn_aud, 'data'));
  fn_e = strrep(fn_aud, 'chan0', sprintf('chan%d', ch_num));
  ephys = double(ncread(fn_e, 'data'));
  
  % filtering to remove noise
  % band pass audio signal
  aud_f = ZZ_butterFilterFunc_v1(aud, fs, 250, 7500, 2);
  % FIR pass ephys signal
  ephys_f = ZZ_FIRbandpass(ephys, fs,  400, 9000, 80);
  
  % find the plot range
  istart = floor(ts(ni)*fs);
  iend = floor((ts(ni)+dur)*fs)-1;
  aud_s = aud_f(istart:iend);
  ephys_s = ephys_f(istart:iend);
  
  % plot
  close all; 
  [fig, axes] = generatePanelGrid_v2(2, 1, [0.4;0.4], [0.05], [0.1;0.05], [0.05;0.05], 0.05, [0 0], [10 10 1200 250]);
  set(fig, 'Renderer', 'painters');
  % first spectrogram
  ax1 = axes(1);
  [power1, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(aud_s, fs, 512, 512, 492, [250 6500], [12 25]);
%   [power1, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(aud_s, fs, 256, 256, 236, [250 7500], [12 25]);
  imagesc(ax1, t, f, power1);
  colormap jet;
  set(ax1, 'Xlim', [t(1) t(end)]);
%   set(ax1, 'Xlim', [2 4]);
  set(ax1, 'YDir', 'normal');
  set(ax1, 'YTick', [2000 4000 6000]);
  set(ax1, 'YTickLabel', {'2k', '4k', '6k'});
  set(ax1, 'XTick', []);
  title(ax1, n);
  ax1.FontSize = 10;
  % then ephys trace
  ax2 = axes(2);
  te = (1:length(ephys_s))/fs;
  plot(ax2, te, ephys_s, 'LineStyle', '-', 'Color', '#000000', 'LineWidth', 1);
  set(ax2, 'Xlim', [te(1) te(end)]);
  set(ax2, 'Ylim', [min(ephys_s)-10 max(ephys_s)+10]);
  set(ax2, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'Box', 'off', 'Visible', 'off');
  
  % save fig
  fd_save_this = fullfile(fd_save, birdID);
  if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end
  fn_fig = fullfile(fd_save_this, sprintf('%s.%s.t%.3f.pdf', birdID, n, ts(ni)));
  print(fig, fn_fig, '-dpdf', '-painters', '-r600');
  
  % save as matlab fig
  fn_fig = fullfile(fd_save_this, sprintf('%s.%s.t%.3f.fig', birdID, n, ts(ni)));
  savefig(fn_fig);
end
  
  
  
  

  
  
  
  
  
  
  
