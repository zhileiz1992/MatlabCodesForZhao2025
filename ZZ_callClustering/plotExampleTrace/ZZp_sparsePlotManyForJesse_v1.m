% plot example ephys traces to use in figures
% Zhilei, 09/21/2025
% plot many example neuron traces for Jesse

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
fd_save = fullfile(fd_home, 'Figures', 'CombinedAnalysis', 'SparseExampleForJesse');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end



%% 1. Inputs
% how much time to plot, unit is second
dur = 10;
fs = 20000;
% choose a few broad neurons as examples
% bis = [1, 1, 2, 2, 4, 3];
bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
clim = [12 22.5];
% Load information about sorted neurons and call subtypes
fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
load(fn_info);



%% 2. Loop through neurons and plot
for ni=4:size(info, 1)
  n = info.neuronID{ni};
  disp(ni);
  infotemp = strsplit(n, '-');
  date_short = infotemp{1};
  date_long = [date_short(1:4) '-' date_short(5:6) '-' date_short(7:8)];
  ch = infotemp{2};
  %   ch_num = str2num(strrep(ch, 'ch', ''));
  ch_num = regexp(ch, '\d+', 'match');
  ch_num = str2num(ch_num{1});
  % load dbase
  fd_e = fullfile(fd_home, 'DbaseFiles', pairID, date_long, birdID, 'warble');
  fn_e = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, date_short, ch);
  load(fullfile(fd_e, fn_e));
  
  % find files that are sorted
  p = dbase.Properties;
  pn = dbase.PropertyNames;
  pi = find(strcmp(pn, 'bProd'));
  fidx = find(p(:, pi)==1);
  
  % read the sound and ephys data
  %   i = fis(ni);
  % select 10 files to plot
  fns = dbase.SoundFiles;
  to_sample = min([20 length(fidx)]);
  %   i_rd = randsample(fidx, to_sample);
  % plot
  close all;
  [fig, axes] = generatePanelGrid_v2(40, 1, zeros(40,1)+0.02, zeros(39,1)+0.002, [0.05;0.05], [0.05;0.05], 0.05, zeros(40,1), [10 10 1200 200*to_sample]);
  set(fig, 'Renderer', 'painters');
  for ii=1:to_sample
    i = fidx(ii);
    
    fn_aud = fullfile(fns(i).folder, fns(i).name);
    aud = double(ncread(fn_aud, 'data'));
    fn_e = strrep(fn_aud, 'chan0', sprintf('chan%d', ch_num));
    ephys = double(ncread(fn_e, 'data'));
    
    % filtering to remove noise
    % band pass audio signal
    aud_f = ZZ_butterFilterFunc_v1(aud, fs, 250, 7500, 2);
    % FIR pass ephys signal
    ephys_f = ZZ_FIRbandpass(ephys, fs,  400, 9000, 80);
    
    
    % first spectrogram
    ax1 = axes(2*ii-1);
    [power1, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(aud_f, fs, 512, 512, 492, [250 6500], clim);
    %   [power1, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(aud_s, fs, 256, 256, 236, [250 7500], [12 25]);
    imagesc(ax1, t, f, power1);
    colormap jet;
    set(ax1, 'Xlim', [t(1) t(end)]);
    %   set(ax1, 'Xlim', [2 4]);
    set(ax1, 'YDir', 'normal');
    set(ax1, 'YTick', [2000 4000 6000]);
    set(ax1, 'YTickLabel', {'2k', '4k', '6k'});
    set(ax1, 'XTick', []);
    %     title(ax1, sprintf('%s %s', birdID, n));
    ax1.FontSize = 10;
    
    % then ephys trace
    ax2 = axes(2*ii);
    te = (1:length(ephys_f))/fs;
    plot(ax2, te, ephys_f, 'LineStyle', '-', 'Color', '#000000', 'LineWidth', 1);
    set(ax2, 'Xlim', [te(1) te(end)]);
    set(ax2, 'Ylim', [min(ephys_f)-10 max(ephys_f)+10]);
    set(ax2, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'Box', 'off', 'Visible', 'off');
  end
  % save as pdf
  fd_save_this = fullfile(fd_save, birdID);
  if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end
  fn_fig = fullfile(fd_save_this, sprintf('%s.n%d.%s.pdf', birdID, ni, n));
  print(fig, fn_fig, '-dpdf', '-painters', '-r600');
end












