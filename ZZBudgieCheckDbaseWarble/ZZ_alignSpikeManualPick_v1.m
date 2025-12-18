% for neurons that are already sorted, manually pick a bunch similar
% syllables, then align the spikes 
% Zhilei Zhao, 02/19/2024

clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes'));

%% inputs
fd_z4 = '/mnt/z4/';
birdID = 'pair1CU21RigB';
data_date = '2023-11-04';
fd_save_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'Figures', 'pair1', data_date, birdID, 'warble');
% load the dbase with segmented syllables and sorted spikes
fn_dbase = fullfile(fd_save_base, 'w_ephys', 'pair1CU21RigB_target_ephys_w.mat');
load(fn_dbase); 
% where to save
fd_save = fullfile(fd_save_base, 'w_ephys_picked');
if ~exist(fd_save)
  mkdir(fd_save)
end


%% fix the title str to reflect syllable order in the new selected dbase
for si=1:length(ephys_selected)
  temp = strsplit(ephys_selected(si).title_str, '-');
  % replace with the new order id
  temp{1} = num2str(si); 
  ephys_selected(si).title_str_new = strjoin(temp, '-');
end
% plot the spectrograms for manual check
prefix = sprintf('%s_w', birdID); 
plotSpectrogramListSpikes_v1(ephys_selected, 'signalFiltered', 1, 0.3, 3, 10, fd_save, prefix, 'title_str_new', [500, 7500], [8, 20]); 


%% manually pick similar syllables
% go through the plotted spectrograms, identify similar ones
p1 = [30, 31, 80, 86, 100, 102, 177, 192, 237, 240, 265, 279, 282, 285, 292, 296, 309, 318, 320, 341, 394, 406, 441, 450, 452];
p2 = [39, 43, 46, 56, 62, 63, 65, 75, 76, 77, 82, 90, 98, 116, 117, 131, 133, 134, 135, 136, 141, 142, 146, 190, 198, 202, 316, 317, 330, 333, 339, 352, 359, 376, 393]; 
p3 = []; 
fd_save_this = fullfile(fd_save, 'p1'); 
plotSpectrogramListSpikes_v1(ephys_selected(p1), 'signalFiltered', 1, 0.3, 3, 10, fd_save_this, prefix, 'title_str_new', [500, 7500], [8, 20]); 


%% Plot ephys voltage trace and sorted by number of spikes
% sort by number of spikes
p = p1; 
ephys_this = ephys_selected(p); 
num_spikes = cellfun(@sum, {ephys_this(:).spike_iv});
[a, sort_i] = sort(num_spikes, 'descend');
ephys_this = ephys_this(sort_i);
plotSpectrogramListSpikes_v1(ephys_this, 'signalFiltered', 1, 0.3, 3, 10, '', prefix, 'title_str_new', [500, 7500], [8, 20]); 
% plot raw traces as well
% FIR pass the raw trace
for si=1:length(ephys_this)
  ephys_this(si).e_trace_FIR = ZZ_FIRbandpass(double(ephys_this(si).e_trace), ephys_this(si).fs,  400, 9000, 80);
end
fd_save_this = fullfile(fd_save, 'p1_sep'); 
plotSpectrogramListTraceSpikes_v2(ephys_this, 'signalFiltered', 1, 0.22, 2, 10, fd_save_this, prefix, '', [500, 7500], [8, 20]); 

 


