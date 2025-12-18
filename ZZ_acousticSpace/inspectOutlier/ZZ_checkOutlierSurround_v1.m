% check the surrounding of outlier/jumping spectrogram windows
% figure out what's causing the jump
clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


% pathname = '/mnt/z4/zz367/EphysMONAO/Analyzed/vaeWav/pair5RigCCU29/Traj/applySyl1/callUMAPone';
% filename = 'pair5RigCCU29.v4.rd300.embedding.csv';
pathname = '/mnt/z4/zz367/EphysMONAO/Analyzed/vaeWav/pair5RigCCU29/Traj/applySyl1/paramSearch1/v4v5.nn25_md0_euclidean_vs0_ke0/';
filename = 'pair5RigCCU29.v4v5.nn25_md0_euclidean_vs0_ke0.embedding.csv';

fn_embed = fullfile(pathname, filename);
embed_data = readtable(fn_embed, 'Delimiter', ',');

fn_h5 = '/mnt/z4/zz367/EphysMONAO/Analyzed/vaeWav/pair5RigCCU29/Traj/Spectrogram1/pair5RigCCU29.spec_goffinet_traj_256_236.h5';
fn_info = strrep(fn_h5, '.h5', '.info.csv');
spec_info = readtable(fn_info, 'Delimiter', ',');
% add a syl_id to the table
spec_info.syl_id = spec_info.fn_wav + "_" + string(spec_info.s_idx);


%% Check windows surrounding the jump
% jump = 2296;
jump = 3808;
vaeCols = find(startsWith(embed_data.Properties.VariableNames, 'vae'));
h5_idx_list = [3806 3807 3808 3809 3810 3867];
vt = embed_data(h5_idx_list, vaeCols);
vt = table2array(vt);
figure(11); set(gcf, 'Position', [50 50 300 900]);
imagesc(vt', [-1 1]); 
colormap('turbo');
xlim([0.5 length(h5_idx_list)+0.5]);
ylim([0.5 32.5]);
xticks(1:length(h5_idx_list));
xticklabels(h5_idx_list);

% calculate distance matrix
D_euc = squareform(pdist(vt, 'euclidean'));
D_cos = squareform(pdist(vt, 'cosine'));
D_cor = squareform(pdist(vt, 'correlation'));
figure(12); set(gcf, 'Position', [50 50 1200 300]);
D_all = {D_euc, D_cos, D_cor};
str_all = {'Euclidean', 'Cosine', 'Correlation'};
for ii=1:size(D_all,2)
    subplot(1,3,ii);
    imagesc(D_all{ii}); 
    colormap('turbo');
    xticks(1:size(D_euc,1));
    xticklabels(h5_idx_list);
    yticks(1:size(D_euc,1));
    yticklabels(h5_idx_list);
    title(str_all{ii});
    colorbar;
end


% how are the latents distributed? Do some latents dominate?
close all;
xmin = min(min(table2array(embed_data(:, vaeCols))));
xmax = max(max(table2array(embed_data(:, vaeCols))));
[fig, axes] = generatePanelGrid_v2(4, 8, [0.15;0.15;0.15;0.15], [0.05;0.05;0.05], [0.05;0.05], [0.05;0.05], 0.02, [0;0;0;0], [10 10 1800 800]);
for ii=0:31
  plot_i = idivide(int32(ii), 8)+1;
  plot_j = mod(ii, 8) + 1;
  ax = axes(plot_i, plot_j);
  vae_d = sprintf('vae%d', ii);
  histogram(ax, embed_data.(vae_d));
  title(ax, vae_d);
  xlim(ax, [xmin xmax]);
end







  
  
  
  
  
  
  
  

