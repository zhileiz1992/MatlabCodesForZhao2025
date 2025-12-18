%% This script takes in budgie contact call spectrogram and outputs hdf5 files for neural network training
%% also contain simple features representing the calls
% Zhilei modified, 02/14/2023

clear all; close all;
addpath(genpath('Y:\zz367\LabSoftware\ZhileiMatlabAudioScripts'));


%% Folder setting
folderBase = 'Y:\zz367\BudgieContactCalls\Cage52ColonyNoiseCalls';
folderSave = fullfile(folderBase, 'VAE_h5');
if exist(folderSave)
    rmdir(folderSave, 's');
end
mkdir(folderSave);
% inputs
fnCSV1 = fullfile(folderBase, 'CombinedColonyNoiseAprilCage52_dataInfo.csv');
fnCSV2 = fullfile(folderBase, 'CombinedColonyNoiseJulyCage52_dataInfo.csv');
dataInfo1 = readtable(fnCSV1, 'Delimiter', ',');
dataInfo1 = table2struct(dataInfo1);
dataInfo2 = readtable(fnCSV2, 'Delimiter', ',');
dataInfo2 = table2struct(dataInfo2);
dataInfo = [dataInfo1; dataInfo2];

%% Loop through each bird and dataset
for bdx=1:size(dataInfo,1)
    [a,expID,b] = fileparts(dataInfo(bdx).WavFolder);
    fprintf('analyze bird: %d', bdx);
    % find the postproof mat from replicates
    fnsResp = dir(fullfile(folderBase, sprintf('%s_response1', expID), '*postproof.mat'));
    disp(fnsResp.name);
    load(fullfile(fnsResp.folder, fnsResp.name));
    resp = callStructProof;
    % select calls that are in time limit, save as .h5
    folderSaveThis = fullfile(folderSave, expID);
    if exist(folderSaveThis)
        rmdir(folderSaveThis,'s');
    end
    mkdir(folderSaveThis);
%     callCount = 0;
    parfor ii=1:size(resp,1)
        callType = resp(ii).title;
        % only look at calls
        if (contains(callType, 'c')) || (contains(callType, 'm'))
%             callCount = callCount + 1;
            this_call = resp(ii);
            filename = strcat(expID, '_call_',num2str(ii, '%05.f'),'.h5');
            fn = fullfile(folderSaveThis, filename);
            h5create(fn, '/spec', [128 128]);
            % exp ID
            h5writeatt(fn,'/spec','budgieID', expID);
            % re-calculate the spectrogram using amplitude normalized signal
            [spec, dt, f ,T]=get_spec(this_call.signalNorm, this_call.fs);
            % duration
            h5writeatt(fn,'/spec','duration',this_call.duration);
            % fundamental freq
            [val_s,ind_f]=max(spec(f>1000&f<5000,:),[],1);
            call_f  = f(ind_f+find(f>1000,1));
            % Mean freq. and std freq.
            mean_f = nanmean(call_f(call_f>500));
            sd_f  = nanstd(call_f(call_f>500));
            % Max freq. and min freq.
            max_f = nanmax(call_f);
            min_f = nanmin(call_f(call_f>500));
            h5writeatt(fn,'/spec','pitch',call_f);
            h5writeatt(fn,'/spec','mean_f',mean_f);
            h5writeatt(fn,'/spec','sd_f',sd_f);
            h5writeatt(fn,'/spec','max_f',max_f);
            h5writeatt(fn,'/spec','min_f',min_f);
            % syllable Entropy
            spectral_entropy = spectralEntropy(this_call.signalNorm, this_call.fs);
            % geometric mean/ arithmetic mean , measure of flatness
            wiener_entropy = geomean(exp(spec),'all','omitnan')/nanmean(exp(spec),'all');
            h5writeatt(fn,'/spec','spectral_entropy',spectral_entropy);
            h5writeatt(fn,'/spec','wiener_entropy',wiener_entropy);
            % syllable spectrogram resize to 128x128
            spec_intp = interpolate_spec(f,dt,[0, size(spec,2)],spec);
            h5write(fn,'/spec',spec_intp);
        end
    end
end
