function [warbleMeta, callMeta, fd_nc_warble, fd_nc_call] = ZZ_extractWarbleCallIntan_func_v11(fd_base_intan, fd_base_save, meta, data_date)
% extract warble and calls from intan recordings
% the algorithm to identify warble is from v11
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/'));
addpath(genpath('/home/zz367/ProjectsU/WarbleAnalysis/Jupyter/MatlabCodes/ZZ_extractWarbleFromWav'));

% fd_save_warble = fullfile(fd_base_save, expID, data_date);
% if exist(fd_save_warble)
%   rmdir(fd_save_warble, 's');
% end
% mkdir(fd_save_warble);
% locate the warble episodes in the Intan data from either bird
% save the rest vocalization events in 'other'
folderData = fullfile(fd_base_intan, data_date);
% saveWavPath = fullfile(fd_save_warble, 'WarbleWav');
[warbleMetaAll, otherMetaAll] = ZZ_IdentifyWarbleFromIntan_v6(fd_base_intan, fd_base_save, data_date, meta);


end

