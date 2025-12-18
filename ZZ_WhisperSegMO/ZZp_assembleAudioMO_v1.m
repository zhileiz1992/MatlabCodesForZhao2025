% A script to assemble audio dataset from the MO ephys project
% to use for training WhisperSeg model and VAE/UMAP model
% Zhilei, 05/12/2025

close all; clear;


%% 1. Input
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'WarbleAnalysis');
% where the ephys data is stored
fd_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles');
fd_nc = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'ExtractedIntan');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
partnerIDs = {'pair5Rig0CU55', '', '', ''};
% whether to extract partner bird as well, only for birds with auditory files sorted
extract_partner = {1, 0, 0, 0};
% where to store the extracted audio
fd_save_base = fullfile(fd_home, 'DataNew', '20250512_MOsorted');


%% 2. Collect audio file names
fns_unique = {};
fns_p_unique = {};
% for focal bird, bProd=1; for partner bird, bAud=1
% bi = 1;
for bi=1:size(birdIDs,2)
  bd = birdIDs{bi};
  % load the neuron information
  fn_meta = fullfile(fd_dbase, pairIDs{bi}, 'MetaInfo', sprintf('%s_sparseInfo.mat', bd));
  load(fn_meta);
  % loop through files, collect audio file names
  fns = {};
  fns_p = {};
  parfor ni=1:size(info,1)
    fn_dbase = fullfile(fd_dbase, pairIDs{bi}, info.date_long{ni}, bd, 'warble', sprintf('%s.%s.warble.good.%s.dbase.mat', bd, info.date{ni}, info.channel{ni}));
    a = load(fn_dbase); dbase = a.dbase;
    % find the index of desired properties
    prod_idx = find(strcmp(dbase.PropertyNames, 'bProd'));
    aud_idx = find(strcmp(dbase.PropertyNames, 'bAud'));
    % what sound files have the desired property
    p = dbase.Properties;
    prod = find(p(:,prod_idx));
    aud = find(p(:,aud_idx));
    f = dbase.SoundFiles;
    % get production files
    fns_this = {};
    if ~isempty(prod)
      for ii=1:length(prod)
        fns_this = [fns_this fullfile(f(prod(ii)).folder, f(prod(ii)).name)];
      end
    end
    fns{ni} = fns_this;
    % get auditory files
    if extract_partner{bi}==1 && ~isempty(aud)
      fns_this = {};
      if ~isempty(aud)
        for ii=1:length(aud)
          temp = fullfile(f(aud(ii)).folder, f(aud(ii)).name);
          fn = strrep(temp, 'chan0.nc', 'chan17.nc');  % use the other audio channel
          fns_this = [fns_this fn];
        end
      end
      fns_p{ni} = fns_this;
    end
  end
  
  % remove redundent files
  fns_unique{bi} = sort(unique([fns{:}]));
  fns_p_unique{bi} = sort(unique([fns_p{:}]));
end


%% 3. Randomly sample files and save
num_sample = 200;
fs = 20000;
% bi = 1;
for bi=1:size(birdIDs,2)
  rng(1992);
  f = fns_unique{bi};
  act_sample = min([num_sample size(f,2)]);
  fn_rd = randsample(f, act_sample);
  % save into a new folder
  fd_save = fullfile(fd_save_base, birdIDs{bi});
  if exist(fd_save)
    rmdir(fd_save, 's');
  end
  mkdir(fd_save);
  % loop through files
  parfor fi=1:size(fn_rd,2)
    % fi = 1;
    % read data
    d = ncread(fn_rd{fi}, 'data');
%     d = d/10.0;
    % d = d;  % normalize to similar range
    [a,b,c] = fileparts(fn_rd{fi});
    fn_save = fullfile(fd_save, [b '.wav']);
    audiowrite(fn_save, d, fs);
  end
  
  % save partner bird if requested
  if extract_partner{bi}==1
    f = fns_p_unique{bi};
    act_sample = min([num_sample size(f,2)]);
    fn_rd = randsample(f, act_sample);
    fd_save = fullfile(fd_save_base, partnerIDs{bi});
    if exist(fd_save)
      rmdir(fd_save, 's');
    end
    mkdir(fd_save);
    % loop through files
    parfor fi=1:size(fn_rd,2)
      % fi = 1;
      % read data
      d = ncread(fn_rd{fi}, 'data');
      % d = d;  % normalize to similar range
      [a,b,c] = fileparts(fn_rd{fi});
      fn_save = fullfile(fd_save, [b '.wav']);
      audiowrite(fn_save, d, fs);
    end
  end
end

















