% plot summary metric from call clustering analysis
% 2025/06/16

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% 1. Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
fd_save_base = fullfile(fd_home, 'vaeWav');
% where to save the combined results
fd_save_comb = fullfile(fd_save_base, 'combined', 'CallClusteringMetric');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what spectrogram input dataset
% input_rn = 'spec_goffinet_cutoff_256_176';
input_rn = 'spec_goffinet_nn_256_176';
% default color list
% col_list = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#737373'};
syl = 'v';  % what syllable to focus on
% what subfolder has the VAE/UMAP results
% UMAP on VAE
% fd_vae = 'UMAPonVAE6'; vae_run = input_rn;
fd_vae = 'UMAPonVAE7'; vae_run = input_rn;


%% 2. Read data
metric_real = {};
metric_random = {};
for bi=1:size(birdIDs,2)
  bd = birdIDs{bi};
  fd_d = fullfile(fd_save_base, bd, fd_vae, syl, vae_run);
  fn_real = fullfile(fd_d, sprintf('%s.UMAPonVAE.metrics.csv', bd));
  metric_real{bi} = readtable(fn_real, 'delimiter', ',');
  fn_random = fullfile(fd_d, sprintf('%s.random.UMAPonVAE.metrics.csv', bd));
  metric_random{bi} = readtable(fn_random, 'delimiter', ',');
end

%% 3. Plot data
% from GPT
close all;
% Metric names and configuration
metric_names = {'silhouette_umap', 'calinski_umap', 'calinski_vae'};
pretty_metrics = {'Silhouette score (UMAP)', 'Calinski-Harabasz index (UMAP)', 'Calinski-Harabasz index (VAE)'}; 
num_metrics = numel(metric_names);
num_birds = numel(metric_real);
x_gap = 3;
x_offset = 0.4;
violin_width = 0.8;

% fig_pos = [100 100 1000 500];
fig_pos = [100 100 800 200];
[fig, axes] = generatePanelGrid_v2(1, 3, [0.75], [], [0.15;0.02], [0.05;0.05], 0.1, [0], fig_pos);% loop through each cluster

for m = 1:num_metrics
%     nexttile;
    ax = axes(m);
    hold(ax, 'on');

    metric = metric_names{m};
    is_log = contains(metric, 'calinski');
    x_labels = {};
    x_ticks = [];

    real_ys = {};  % to store real y-values
    sim_ys = {};   % to store all simulated y-values

    for b = 1:num_birds
        y_sim_all = metric_random{b}.(metric);
        y_real = metric_real{b}.(metric);

        if is_log
            y_sim = y_sim_all(y_sim_all > 0);
            if isempty(y_sim)
                warning('Bird%d: All simulated %s values are non-positive.', b, metric);
                continue;
            end
            support_min = max(min(y_sim), 1e-6);
        else
            y_sim = y_sim_all;
            support_min = min(y_sim);
        end

        support_max = max(y_sim);
        y_grid = linspace(support_min, support_max, 200);
        [f, y_vals] = ksdensity(y_sim, y_grid, 'Function', 'pdf');

        if is_log
            valid = (y_vals > 0);
            y_vals = y_vals(valid);
            f = f(valid);
            y_grid = y_grid(valid);
        end

        if isempty(f) || all(f == 0)
            warning('Bird%d: Failed to generate valid density for %s.', b, metric);
            continue;
        end

        f = f / max(f) * violin_width / 2;

        x0 = (b - 1) * x_gap + 1;
        x_ticks(end+1) = x0;
%         x_labels{end+1} = sprintf('Bird%d', b);
        x_labels{end+1} = pretty_ids{b};
        fill(ax, [x0 - f, fliplr(x0 + f)], [y_grid, fliplr(y_grid)], ...
             [0.22 0.49 0.72], 'EdgeColor', 'none', 'FaceAlpha', 0.75);

        if ~is_log || y_real > 0
            scatter(ax, x0 + x_offset, y_real, 120, 'r', 'filled', 'MarkerEdgeColor', 'none');
            real_ys{end+1} = y_real;
        end

        sim_ys{end+1} = y_sim;
    end

    xticks(ax, x_ticks + x_offset / 2);
    xticklabels(ax, x_labels);
    title(ax, pretty_metrics{m});
    xlim(ax, [min(x_ticks) - 1, max(x_ticks) + 2]);
    if is_log
        set(ax, 'YScale', 'log');
    end
    ylabel(ax, pretty_metrics{m});

    % Add 10% headroom above max value
    sim_vals_all = vertcat(sim_ys{:});
    real_vals_all = vertcat(real_ys{:});
    all_vals = [sim_vals_all; real_vals_all];

    if ~isempty(all_vals)
        ymax = max(all_vals);
        ymin = min(all_vals);
        if is_log
            ylim(ax, [max(ymin, 1e-6)*0.9, ymax * 1.5]);
        else
            ylim(ax, [ymin-0.02, ymax + 0.1 * abs(ymax)]);
        end
    end


    box off;
end

% save fig
fn_pdf = fullfile(fd_save_comb, 'callClusteringMetric.pdf');
print(fig, fn_pdf, '-dpdf', '-painters');
  


  
  
  
  
  
  
  
  
  
  
  
  
  