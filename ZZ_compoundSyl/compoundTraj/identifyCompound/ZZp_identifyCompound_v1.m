% identify and visualize compound syllables

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load the information about syllables
fd_data = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'ApplySylAll', sprintf('latents.%s', vae_run));
fn_info = fullfile(fd_data, sprintf('%s.latents.%s.info.csv', birdID, vae_run));
info = readtable(fn_info, 'Delimiter', ',');
% for acoustic analysis, focus on 'batch 1', since 'batch 2' may contains sound from partner bird
info = info(strcmp(info.batch, 'batch1'), :);


%% 2. Set a duration threshold for category 'b' and 'x'
% add a duration column to the table
fs = 20000;
info.dur = (info.iend - info.istart) / fs;
% vary threshold, check how many 'b' and 'x' pass
thre_list = [0.2, 0.3, 0.4, 0.5, 0.6];
pass_all = zeros(length(thre_list), 2);

for ti=1:length(thre_list)
  dur_thre = thre_list(ti);  % unit is sec
  % find 'b' and 'x' that pass the threshold
  comp = info(info.dur>=dur_thre & ismember(info.call_subtype, {'b', 'x'}), :);
  % what's the percentage of syllables that pass the threshold
  pass_b = sum(strcmp(comp.call_subtype, 'b')) / sum(strcmp(info.call_subtype, 'b'));
  pass_x = sum(strcmp(comp.call_subtype, 'x')) / sum(strcmp(info.call_subtype, 'x'));
  pass_all(ti, :) = [pass_b pass_x];
end

% plot res
close all;
figure;
plot(thre_list, pass_all(:,1), 'ro', 'MarkerFaceColor', 'red', 'DisplayName', 'b'); hold on;
plot(thre_list, pass_all(:,2), 'bo', 'MarkerFaceColor', 'blue', 'DisplayName', 'x'); hold on;
legend;
xlabel('Duration threshold (sec)');
ylabel('Proportion that pass threshold');


%% 3. Plot spectrograms to get intuition
dur_thre = 0.3;  % unit is sec
comp = info(info.dur>=dur_thre & ismember(info.call_subtype, {'b', 'x'}), :);
fd_save = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ExampleSpec');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end

syls = {'b', 'x'};
rng(1992);
for syl_i=1:size(syls,2)
  idx = find(strcmp(comp.call_subtype, syls{syl_i}));
  % randomly sample 300
  idx_rd = randsample(idx, 600);
  % load wav data
  d = cell(length(idx_rd), 1);
  pad = 0.015; pad_pt = floor(fs*pad);
  parfor ii=1:length(idx_rd)
    idx = idx_rd(ii);
    [signal, fs] = audioread(comp.fn_wav{idx});
    i_start = max([1 comp.istart(idx)-pad_pt]);
    i_end = min([size(signal,1) comp.iend(idx)+pad_pt]);
    d{ii} = signal(i_start:i_end);
  end
  % plot spectrograms
  for ii=1:20
    close all;
    [fig, axes] = generatePanelGrid_v2(3, 10, [0.25;0.25;0.25], [0.05;0.05], [0.05;0.05], [0.05;0.05], 0.01, [0;0;0], [10 10 2000 800]);
    for ri=1:30
      plot_i = floor((ri-1)/10) + 1;
      plot_j = mod(ri-1, 10) + 1;
      ax = axes(plot_i, plot_j);
      audio = d{30*(ii-1)+ri};
      [ax, ~, ~,  ~, ~,  ~] = showAudioSpectrogramZZ_flexible_v1(audio, fs, ax, [250 7500], [12 23], 256, 256, 236);
      axis(ax, 'off');
    end
    fn_fig = fullfile(fd_save, sprintf('ExampleSpec.%s.%s.%d.pdf', birdID, syls{syl_i}, ii));
    print(fig, fn_fig, '-dpdf', '-painters');
  end
end

















