% Calculte in the MMD matrix, how many high-similarity strands exist
% Zhilei, 10/21/2025


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl/compoundTraj/pairwiseCorrelation'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};

% rename relationship
syls_all_list = {{'v1', 'v2', 'v3', 'v4', 'v5', 'v6'}, ...
             {'v1', 'v2', 'v3', 'v4'}, ...
             {'v1', 'v2', 'v3'}, ...
             {'v1', 'v2', 'v3', 'v4'}};

% raw color limits for neural and acoustic plots
clims_n_all = {[0.5 0.85]; [0.5 0.85]; [0.5 0.85]; [0.5 0.85]};
clims_a_all = {[-0.2 -0.03]; [-0.25 -0.03]; [-0.15 -0.015]; [-0.2 -0.03]};
% after convert to similarity index
clims_n2_all = {[0.35 0.95]; [0.3 0.95]; [0.4 0.95]; [0.35 0.95]};
clims_a2_all = {[0 1]; [0 1]; [0 1]; [0 1]};

           
bi = 4;
  
birdID = birdIDs{bi};
pairID = pairIDs{bi};
disp(birdID);
% what good subtypes to analyze
syls = syls_all_list{bi};

col_list = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#ccad25'};
for si=1:size(syls,2); col_dict.(syls{si})=col_list{si}; end
% where is calculated decoding matrix located
vae_run = 'traj_chop_32_1_32';
fd_data = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix5'], 'intermediate');
% what's the window size and hop size when calculating SVM decoding matrix, unit is seconds
fs = 20000;
pad = 0.08;
bin = 0.01;
hop = 0.002;
% where is the MMD result
fd_mmd = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'MMD', 'MatrixCall3', 'all', 'intermediate');
sigma = 0.9;
win_frame = 32;
hop_frame = 1;
sec_per_frame = 0.001;
% where to save the plots
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix5'], 'shareFracMMD');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end
fprintf('Save SVM plots to %s\n', fd_save);


%% 1.2 Read acoustic MMD matrix in
acoRaw = cell(size(syls,2), size(syls,2));
aco = cell(size(syls,2), size(syls,2));  % mean accuracy matrix of the real dataset
% also trim to desired time range
rel_t_aco = (0:1000)*hop_frame*sec_per_frame - win_frame*sec_per_frame;  % a long time cooridiante axis for the sliding window start
% rel_t_aco = rel_t_aco + win_frame*sec_per_frame/2;  % convert to window center
% tpad_a = [0.032 0];  % how much time before syllable onset and after syllable offset
tpad_a = [0 -0.032];
aco2 = cell(size(syls,2), size(syls,2));  % time-trimmed matrix
aco_shuff = cell(size(syls,2), size(syls,2)); % shuffled data
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    run_name = sprintf('%s_%s', syls{ii}, syls{jj});
    fn_res = fullfile(fd_mmd, sprintf('%s.%s.%.2f.mmd.mat', birdID, run_name, sigma));
    a = load(fn_res); a_this = a.mmd.C;
    aco{ii, jj} = a_this;
    acoRaw{ii,jj} = a.mmd.mmd;
    % trim to desired time range
    ty = rel_t_aco(1:size(a_this,1));
    tx = rel_t_aco(1:size(a_this,2));
    tendy = ty(end) + tpad_a(2);  % note that for MMD, no more extra slide at the end, sliding window stops at the syllable offset
    iy = find((ty>=-tpad_a(1)) & (ty<=tendy));
    tendx = tx(end) + tpad_a(2);
    ix = find((tx>=-tpad_a(1)) & (tx<=tendx));
    aco2{ii,jj} = a_this(iy, ix);
    aco_shuff{ii,jj} = a.mmd.mmd11; 
  end
end
% tile into a large 2d matrix
lens_a = cellfun(@(x) size(x,2), aco2(1,:));
aco3 = ZZfunc_tileArrayCell_v1(aco2);
% save for later use
rel_t_a_chop = rel_t_aco(rel_t_aco>=-tpad_a(1));
d_aco.aco=aco; d_aco.aco2=aco2; d_aco.aco3=aco3; d_aco.acoRaw=acoRaw; d_aco.lens_a=lens_a; d_aco.aco_shuff=aco_shuff;
d_aco.rel_t_aco=rel_t_aco; d_aco.rel_t_a_chop=rel_t_a_chop;
fn_aco = fullfile(fd_save, sprintf('%s.d_aco.mat', birdID));
save(fn_aco, 'd_aco');


%% 2.1 Plot the full matrix of raw values
% acoustic MMD
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
v_loc = cumsum(lens_a)-lens_a/2;
line_loc = cumsum(lens_a);
clim_a = clims_a_all{bi};
ZZfunc_showMMDmatrix_v5(ax, aco3, clim_a, v_loc, syls, 'Acoustic dMMD', true, line_loc, '#1f78b4', '--', 1, gray, false);
fn_pdf = fullfile(fd_save, sprintf('%s.acousticMMDmatrix.t%.3f_%.3f.c%.2f_%.2f.pdf', birdID, tpad_a(1), tpad_a(2), clim_a(1), clim_a(2)));
print(fig, fn_pdf, '-dpdf', '-painters');


%% 2.2. Change to a uniform similarity index
% acoustic MMD
vmin = clims_a_all{bi}(1); 
vmax = clims_a_all{bi}(2);
aco4 = (vmax - aco3) / (vmax - vmin);
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
v_loc = cumsum(lens_a)-lens_a/2;
line_loc = cumsum(lens_a);
ZZfunc_showMMDmatrix_v5(ax, aco4, [0 1], v_loc, syls, 'Acoustic MMD', true, line_loc, '#1f78b4', '--', 1, flipud(gray), false);
fn_pdf = fullfile(fd_save, sprintf('%s.acousticMMDmatrix.t%.3f_%.3f.c%.2f_%.2f.converted.pdf', birdID, tpad_a(1), tpad_a(2), vmin, vmax));
print(fig, fn_pdf, '-dpdf', '-painters');



%% 3. Identify the high-similarity strands
v_thre = 0.05;  % acoustic similarity higher
t_thre = 20;  % duration longer, unit is ms
vmin = clims_a_all{bi}(1); 
vmax = clims_a_all{bi}(2);
pass_list = [];
count = 0;
n = size(aco2,1); 
close all;
[fig, axs] = generatePanelGrid_v2(n, n, zeros(n,1)+0.12, zeros(n-1,1)+0.02, [0.05;0.05], [0.05;0.05], 0.02, zeros(n,1), [10 10 900 900]);
for ii=1:size(aco2,1)
  for jj=(ii+1):size(aco2,2)
    count = count + 1;
    % convert MMD distance to similarity
    distm = aco2{ii, jj};
%     distm = (vmax - distm) / (vmax - vmin);
%     % identify strands
%     pass_i = ZZfunc_identifySimilarityStrand_v1_noPlot_larger(distm, [3 3], v_thre, t_thre);
    pass_i = ZZfunc_identifySimilarityStrand_v1_noPlot(distm, [3 3], vmax*1.25, t_thre);
    mask = zeros(size(distm));
    if ~isempty(pass_i)
      for pi=1:size(pass_i,2)
        mask(pass_i{pi}) = 1; 
      end
      ax = axs(ii, jj);
      imagesc(ax, mask);
    end
    pass_list(count).birdID = birdID; 
    pass_list(count).ii = ii;
    pass_list(count).jj = jj;
    pass_list(count).pass_i = pass_i;
    pass_list(count).num_strand = size(pass_i, 2);
    % get duration as well
    pass_list(count).sylDur_ii = size(distm,1);
    pass_list(count).sylDur_jj = size(distm,2);
    % duration of the high-similarity strands
    pass_list(count).strandDur_ii = 0;
    pass_list(count).strandDur_jj = 0;
    for pi=1:size(pass_i,2)
      [yy, xx] = ind2sub(size(distm), pass_i{pi});
      pass_list(count).strandDur_ii = pass_list(count).strandDur_ii + max(xx) - min(xx);
      pass_list(count).strandDur_jj = pass_list(count).strandDur_jj + max(yy) - min(yy);
    end
    % calculate fraction
    pass_list(count).durFrac_ii = pass_list(count).strandDur_ii / pass_list(count).sylDur_ii;
    pass_list(count).durFrac_jj = pass_list(count).strandDur_jj / pass_list(count).sylDur_jj;
    
  end
end
fn_fig = fullfile(fd_save, 'identifiedStrands_mmd.pdf');
print(fig, fn_fig, '-dpdf', '-painters');

% save data
fn_data = fullfile(fd_save, sprintf('%s.pass_list.mat', birdID));
save(fn_data, 'pass_list');




%% 4. Pull results from all birds
pass_comb = [];
for bi=1:size(birdIDs, 2)
  fd_save = fullfile(fd_base, 'Figures', pairIDs{bi}, 'AcousticSpace', birdIDs{bi}, 'embed', [vae_run '_decodePopMatrix5'], 'shareFracMMD');
  fn =  fullfile(fd_save, sprintf('%s.pass_list.mat', birdIDs{bi}));
  a = load(fn); 
  pass_comb = [pass_comb a.pass_list];
end

% total number of pairwise comparison
fprintf('Total number of pairwise comparison: %d\n', size(pass_comb, 2));

% what fraction of comparison have strands
has_strand = ([pass_comb.num_strand]>0);
fprintf('Percent that has high-similarity strands: %.2f\n', sum(has_strand)/size(pass_comb,2) * 100);

% how much in terms of call duration?
frac_dur = mean([[pass_comb.durFrac_ii] [pass_comb.durFrac_jj]]);
fprintf('Occupy percent of all call duration: %.2f\n', frac_dur * 100);
















